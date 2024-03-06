from pathlib import Path
import numpy as np
import pandas as pd
from dry_pipe import DryPipe
pd.options.mode.chained_assignment = None


def parse_sam_file(samfile, threshold=30) -> pd.DataFrame:
    # import the data
    df = pd.read_csv(
        samfile,
        dtype={
            "chr_name": str,
            "pos": "Int64",
            "flag": "Int64",
            "name": str,
            "score": "Int64",
            "cigar": str
        },
        sep="\t",
        header=None,
        usecols=range(6),
        on_bad_lines='warn',
        names=["name", "flag", "chr_name", "pos", "score", "cigar"],
        engine="python"
    )
    # reads first in pair are foward ou reverse? if foward keep those reads, if else, drop them:
    strand= []
    first_in_pair = []
    for i in df["flag"]:
        binairy_number = bin(i)[2:]
        while len(binairy_number) < 7:
            binairy_number = "0"+ binairy_number
        binairy_number = binairy_number[::-1]
        #if first in pair
        if binairy_number[6] == "1":
            first_in_pair.append(1)
            if binairy_number[4] == "0":
                strand.append("+")
            else:
                strand.append("-")
        else:
            first_in_pair.append(0)
            if binairy_number[4] == "0":
                strand.append("+")
            else:
                strand.append("-")

    df["strand"] = strand
    df["first_in_pair"] = first_in_pair

    df = df.loc[df["first_in_pair"] == 1]

    # use regex to extract each int/marker and build a secondary df holding these info.
    #
    # typic input:
    # ------------
    #
    # 20      42M
    # 25    72M3S
    # Name: cigar
    #
    # typic output:
    # -------------
    #
    #           size marker
    #    match
    # 20 0        42      M
    # 25 0        72      M
    #    1         3      S
    cigar_df = df.cigar.str.extractall(r"(?P<size>\d+)(?P<marker>[A-Z]{1})")
    cigar_df["size"] = cigar_df["size"].astype(int)

    # because needed to test validity, extract first and last cigar marker.
    df["cigar_first_marker"] = cigar_df.marker.groupby(level=0).first()[df.index]
    df["cigar_last_marker"] = cigar_df.marker.groupby(level=0).last()[df.index]

    # the size can be computed with the sum of the size associated with the
    # markers 'M', 'D', 'N', '=' and 'X'

    df["aligned_read_length"] = (
        cigar_df.loc[cigar_df.marker.str.contains("[MDN=X]")].groupby(level=0).sum()
    )
    # two case of validity: positive strand means we need first 'M' (match) marker,
    # negative strand means we need last 'M' (match) marker
    is_valid_positive = (df.strand == "+") & (df.cigar_first_marker == "M")
    is_valid_negative = (df.strand == "-") & (df.cigar_last_marker == "M")
    # we accept reads that are either a valid positive, negative positive and which
    # length are more than a threshold.
    df["valid"] = (is_valid_positive | is_valid_negative) & (
        df.aligned_read_length > threshold
    )

    pos_cor = (3,-6) # ajout d'un paramètre pour la correction de la position dépendante de la taille de la duplication (initialement 5) SJ 18/02/2022
    valid_df = df[df.valid]
    valid_df["start"] = np.where(
        valid_df.strand == "+",
        valid_df.pos + pos_cor[0],
        valid_df.pos + valid_df.aligned_read_length + pos_cor[1],
    )

    # if second in pair: WARNING THE ORIENTATION OF resistance gene IS OPPOSED TO THE ADAPTER.May not be the case in every design!
    # the inversion is appropriate with pDM023 = Meso, but not with pFG051 = E. coli
    # valid_df = valid_df.replace(['+','-'],['-','+'])

    return valid_df


def sam_to_site(sam_df: pd.DataFrame) -> pd.DataFrame:
    # we aggregate the scores by start
    score_pos = (
        sam_df[sam_df.strand == "+"].groupby("start").score.count().rename("score_pos")
    )
    score_neg = (
        sam_df[sam_df.strand == "-"].groupby("start").score.count().rename("score_neg")
    )
    score_total = sam_df.groupby("start").score.count().rename("score_total")
    # we take the first chr_name. Should be the same for each occurence
    chr_name = sam_df.groupby("start").chr_name.first()

    # Concatenate the aggregates
    site_df = pd.concat(
        [chr_name, score_pos, score_neg, score_total], axis=1
    ).reset_index()
    # fill the last values
    site_df["end"] = site_df.start + 1
    site_df["name"] = "i"
    # set the empry score to 0
    site_df = site_df.fillna({"score_pos": 0, "score_neg": 0, "score_total": 0})
    # reorder the columns and sort by start values
    site_df = site_df[
        ["chr_name", "start", "end", "name", "score_pos", "score_neg"]
    ].sort_values(by=["chr_name","start"])
    return site_df


def site_to_interval(
    site_df: pd.DataFrame, *, normalization_value: int, score_threshold=1
) -> pd.DataFrame:
    total_reads = site_df["score_pos"].sum() + site_df["score_neg"].sum()

    pos_interval_df = (
        site_df[site_df["score_pos"] > score_threshold]
        .drop("score_neg", axis=1)
        .rename(columns={"score_pos": "score"})
    )
    pos_interval_df["strand"] = "+"

    neg_interval_df = (
        site_df[site_df["score_neg"] > score_threshold]
        .drop("score_pos", axis=1)
        .rename(columns={"score_neg": "score"})
    )
    neg_interval_df["strand"] = "-"

    interval_df = pd.concat([pos_interval_df, neg_interval_df])
    interval_df["normalized_score"] = (
        (interval_df.score / interval_df.score.sum() * normalization_value)
        .astype(int)
        .clip(1)
    )
    interval_df["fake_score"] = 999
    interval_df = interval_df.sort_values(by=["chr_name","start"])
    return interval_df


def build_bed_file(interval_df: pd.DataFrame, filename):
    """
    bed_file
    --------

    "chr_name", "start", "end", "name", "fake_score", "strand"

    pVCR94deltaX3deltaacr2FRT	26	27	i	999	+
    pVCR94deltaX3deltaacr2FRT	47	48	i	999	+
    pVCR94deltaX3deltaacr2FRT	68	69	i	999	+
    pVCR94deltaX3deltaacr2FRT	97	98	i	999	-
    pVCR94deltaX3deltaacr2FRT	245	246	i	999	-
    """
    interval_df[["chr_name", "start", "end", "name", "fake_score", "strand"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_stranded_unnormalized_bedgraph_file(
    interval_df: pd.DataFrame, filename):
    """
    stranded_unnormalized_bedgraph_file
    --------

    "chr_name", "start", "end", "score"

    pVCR94deltaX3deltaacr2FRT	26	27	2400
    pVCR94deltaX3deltaacr2FRT	47	48	60
    pVCR94deltaX3deltaacr2FRT	68	69	60
    pVCR94deltaX3deltaacr2FRT	97	98	1380
    pVCR94deltaX3deltaacr2FRT	245	246	3780
    """
    interval_df[["chr_name", "start", "end", "score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_stranded_bedgraph_file(interval_df: pd.DataFrame, filename):
    """
    stranded_bedgraph_file
    --------

    "chr_name", "start", "end", "normalized_score"

    pVCR94deltaX3deltaacr2FRT	26	27	2189
    pVCR94deltaX3deltaacr2FRT	47	48	54
    pVCR94deltaX3deltaacr2FRT	47	48	54
    pVCR94deltaX3deltaacr2FRT	68	69	54
    pVCR94deltaX3deltaacr2FRT	97	98	1258
    pVCR94deltaX3deltaacr2FRT	245	246	3448
    """
    interval_df[["chr_name", "start", "end", "normalized_score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )

def build_unstranded_bedgraph_file(interval_df: pd.DataFrame, filename):
    """
    unstranded_bedgraph_file
    --------

    "chr_name", "start", "end", "normalized_score"

    pVCR94deltaX3deltaacr2FRT	26	27	2189
    pVCR94deltaX3deltaacr2FRT	47	48	108<-(score for insertions on the same position, but different strands are merged)
    pVCR94deltaX3deltaacr2FRT	68	69	54
    pVCR94deltaX3deltaacr2FRT	97	98	1258
    pVCR94deltaX3deltaacr2FRT	245	246	3448
    """

    interval_df = interval_df.groupby(["chr_name", "start", "end"],sort=False).agg({'normalized_score':'sum'}).reset_index()
    interval_df[["chr_name", "start", "end", "normalized_score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_scored_bed_file(interval_df: pd.DataFrame, filename):
    """
    scored_bed_file
    --------

    "chr_name", "start", "end", "name", "score", "strand"

    pVCR94deltaX3deltaacr2FRT	26	27	i	2400	+
    pVCR94deltaX3deltaacr2FRT	47	48	i	60	+
    pVCR94deltaX3deltaacr2FRT	68	69	i	60	+
    pVCR94deltaX3deltaacr2FRT	97	98	i	1380	-
    pVCR94deltaX3deltaacr2FRT	245	246	i	3780	-
    """
    interval_df[["chr_name", "start", "end", "name", "normalized_score", "strand"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )

create_scored_bedfile = True
create_stranded_bedgraph = True
create_unstranded_bedgraph = True
create_bedfile = True

@DryPipe.python_call()
def sam2sites(
    input_filename,
    normalization_value,
    read_len_threshold,
    score_threshold,
    __task_output_dir
):

    create_stranded_unnormalized_bedgraph = True
    
    basename = Path(input_filename.name).basename().stripext()

    sam_df = parse_sam_file(input_filename, read_len_threshold)
    site_df = sam_to_site(sam_df)
    interval_df = site_to_interval(
        site_df,
        normalization_value=normalization_value,
        score_threshold=score_threshold,
    )

    if create_bedfile:
        build_bed_file(interval_df, f"{__task_output_dir}/{basename}.bed")
    if create_scored_bedfile:
        build_scored_bed_file(interval_df, f"{__task_output_dir}/{basename}_scored.bed")
    if create_stranded_bedgraph:
        build_stranded_bedgraph_file(
            interval_df, f"{__task_output_dir}/{basename}_stranded.bg")
    if create_unstranded_bedgraph:
        build_unstranded_bedgraph_file(
            interval_df, f"{__task_output_dir}/{basename}_unstranded.bg")
    if create_stranded_unnormalized_bedgraph:
        build_stranded_unnormalized_bedgraph_file(
            interval_df, f"{__task_output_dir}/{basename}_stranded_unnormalized.bg"
        )
