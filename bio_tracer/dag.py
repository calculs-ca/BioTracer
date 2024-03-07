import os.path
from itertools import groupby

from bio_tracer.sam2sites import sam2sites
from bio_tracer.sites2genes import sites2genes
from dry_pipe import TaskConf

adapter_seq = 'ACTGTCTCTTATACACATCT'

genome_fasta = '/acl/jacquesp-rodrigu1/Snakemake2Drypipe/2_References/BW25113_KAN_pFG051.fa'
genome_features = '/acl/jacquesp-rodrigu1/Snakemake2Drypipe/2_References/BW25113_features.bed'
mappable_features = '/acl/jacquesp-rodrigu1/Snakemake2Drypipe/2_References/BW25113_mappable_features.bed'

samples = [
    'RK1723_20231106_R00_R1.fastq.gz', 'RK1723_20231106_R00_R2.fastq.gz'
]


task_conf = TaskConf(
    executer_type="process",
    extra_env={
        "MUGQIC_INSTALL_HOME": "/cvmfs/soft.mugqic/CentOS6"
    }
    # implicit:
    #python_bin="/home/maxl/miniconda3/envs/Genomics/bin/python"
)



def name_and_read_number_from_filename(f):
    try:
        # strip .fastaq.gz from  'RK1723_20231106_R00_R1.fastq.gz'
        prefix, _, _ = f.split('.', 2)

        # prefix -> 'RK1723_20231106_R00_R1'

        name, read_number = prefix.rsplit('_', 1)

        # name -> 'RK1723_20231106_R00'
        # read_number -> 'R1'

        return name, read_number
    except ValueError as e:
        raise Exception(f"failed to parse file {f} {e}")

def pair_samples(files):

    for _, pair in groupby(sorted(files), key=lambda f: name_and_read_number_from_filename(f)[0]):

        read1, read2 = sorted(list(pair), key=lambda f: name_and_read_number_from_filename(f)[1])
        yield read1, read2



def dag_gen(dsl):

    genome_fasta_basename = os.path.basename(genome_fasta)

    index_genome_task = dsl.task(
        key="index_genome"
    ).inputs(
        genome_fasta=dsl.file(genome_fasta)
    ).outputs(
        local_genome_fasta=dsl.file(genome_fasta_basename),
        i1=dsl.file(f'{genome_fasta_basename}.amb'),
        i2=dsl.file(f'{genome_fasta_basename}.bwt'),
        i3=dsl.file(f'{genome_fasta_basename}.sa'),
        i4=dsl.file(f'{genome_fasta_basename}.ann'),
        i5=dsl.file(f'{genome_fasta_basename}.pac')
    ).calls("""
        #!/usr/bin/bash
        set -xe
        
        cp $genome_fasta $local_genome_fasta

        export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6        
        module use $MUGQIC_INSTALL_HOME/modulefiles            
        module add mugqic/bwa/0.7.17
        
        bwa index $local_genome_fasta
    """)()

    yield index_genome_task

    for read_fastq_1, read_fastq_2 in pair_samples(samples):

        sample_pair_name = name_and_read_number_from_filename(read_fastq_1)[0]

        fastp_task = dsl.task(
            key=f"fastp.{sample_pair_name}"
        ).inputs(
            read1=dsl.file(read_fastq_1),
            read2=dsl.file(read_fastq_2),
            adapter_seq=adapter_seq
        ).outputs(
            fast_qc=dsl.file(f"{sample_pair_name}.fast_qc.zip"),
            cutadapt_r1=dsl.file(f"{sample_pair_name}_R1_cutadapt.fq"),
            cutadapt_r2=dsl.file(f"{sample_pair_name}_R2_cutadapt.fq"),
            fastp_out1=dsl.file(f"{sample_pair_name}_R1_qtrim.fq"),
            fastp_out2=dsl.file(f"{sample_pair_name}_R2_qtrim.fq"),
            fastp_json_out=dsl.file(f"{sample_pair_name}.fastp.json"),
            fastp_html_out=dsl.file(f"{sample_pair_name}.fastp.html"),
            fail_out=dsl.file(f"{sample_pair_name}.failed"),
            unpaired_1_out=dsl.file(f"{sample_pair_name}_R1.unpaired"),
            unpaired_2_out=dsl.file(f"{sample_pair_name}_R2.unpaired")
        ).calls("""
            #!/usr/bin/bash
            set -xe                 
            module add fastqc/0.11.9            
            fastqc $read1 $read2 -t 24 -o $__task_output_dir
        """).calls("""
            #!/usr/bin/bash
            set -xe                 

            export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
            module use $MUGQIC_INSTALL_HOME/modulefiles
                        
            cutadapt --adapter $adapter_seq -A $adapter_seq --minimum-length 5 --cores 4 \\
              --output $cutadapt_r1 \\
              --paired-output $cutadapt_r2 \\
              $read1 $read2

        """).calls("""
            #!/usr/bin/bash
            set -xe
            
            export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
            module use $MUGQIC_INSTALL_HOME/modulefiles
            module add mugqic/fastp/0.23.2
            
            fastp --disable_adapter_trimming --cut_right \\
                --cut_right_window_size 4 \\
                --cut_mean_quality 20 \\
                --length_required 10 \\
                -w 24 \\
                -i $read1 -o $fastp_out1 \\
                -I $read2 -O $fastp_out2 \\
                -j $fastp_json_out -h $fastp_html_out \\ 
                --failed_out $fail_out \\
                --unpaired1 $unpaired_1_out --unpaired2 $unpaired_2_out
        """)()

        yield fastp_task

        map_reads_task = dsl.task(
            key=f"map_reads.{sample_pair_name}"
        ).inputs(
            index_genome_task.outputs.i1,
            index_genome_task.outputs.i2,
            index_genome_task.outputs.i3,
            index_genome_task.outputs.i4,
            index_genome_task.outputs.i5,
            genome_fasta=index_genome_task.outputs.local_genome_fasta,
            r1=fastp_task.outputs.fastp_out1,
            r2=fastp_task.outputs.fastp_out2,
            norm_value=10000000,
            read_len_threshold=20,
            score_threshold=1,
            genome_features=dsl.file(genome_features),
            mappable_features=dsl.file(mappable_features)
        ).outputs(
            mapped_sam=dsl.file(f"mapped_{sample_pair_name}.sam"),
            flagstat=dsl.file(f"{sample_pair_name}_flagstat.txt"),
            q30_sam=dsl.file(f"{sample_pair_name}.q30.sam"),
            q30_stranded_bg=dsl.file(f"{sample_pair_name}.q30_stranded.bg"),
            intersect_bed=dsl.file(f"{sample_pair_name}_allfeat_intersect.bed"),
            mappable_intersect_bed=dsl.file(f"{sample_pair_name}_mappable_intersect.bed")
        ).calls("""
            #!/usr/bin/bash

            set -xe
            
            export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
            module use $MUGQIC_INSTALL_HOME/modulefiles            
            module add mugqic/bwa/0.7.17        
            
            bwa mem -k 10 -t 24 $genome_fasta $r1 $r2 -o $mapped_sam        
            
            module add mugqic/samtools/1.14

            # flagstat            
            samtools flagstat -@ 24 $mapped_sam > $flagstat
            
            # filter
            samtools view --threads 24 -S -q 30 $mapped_sam -o $q30_sam
        """).calls(
            sam2sites
        ).calls(
            """
            #!/usr/bin/bash
            set -xe
            
            export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
            module use $MUGQIC_INSTALL_HOME/modulefiles                                            
            module add mugqic/bedtools/2.30
            
            # intersect
            bedtools intersect -a $genome_features -b $q30_stranded_bg -wo > $intersect_bed
            
            # mappability
            bedtools intersect -a $mappable_features -b $q30_stranded_bg -wo > $mappable_intersect_bed"
            """
        ).calls(
            sites2genes
        )()

        yield map_reads_task

    for _ in dsl.query_all_or_nothing("map_reads.*"):
        yield dsl.task(
            key="multiqc"
        ).calls("""
            #!/usr/bin/bash
            set -xe
            
            export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
            module use $MUGQIC_INSTALL_HOME/modulefiles                                
            module add mugqic/MultiQC/1.14
            
            multiqc \\
                $__pipeline_instance_dir/output/fastqc.*/*.zip \\
                $__pipeline_instance_dir/output/fastp.*/* \\                    
                $__pipeline_instance_dir/flagstat.*/* \\
                -o multiQC                        
        """)()



if __name__ == '__main__':

    for read_fastq_1, read_fastq_2 in pair_samples(samples):
        sample_pair_name = name_and_read_number_from_filename(read_fastq_1)[0]
        print(f"sample_pair_name: {sample_pair_name}")