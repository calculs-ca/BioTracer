#ajout d'un score_E correspondant au degrés d'essentialité estimé sous forme de
#variable continue de 0 à l'inf où 0 == E, >1==F, >2==N

#!/usr/bin/env Rscript

# PODNAME: tradis_essentiality.R
# ABSTRACT: tradis_essentiality.R

library("MASS")
options(warn=-1)
options(width=80)

args <- commandArgs(trailingOnly = TRUE)
input = args[1]

if( is.na(input) ){
	cat(paste("Usage: tradis_essentiality.R data.tab\n\n"))
	cat(strwrap("Produces calls of gene essentiality using an adaptation of the
	method described in Langridge et al. Genome Research 2009 and Barquist et al.
	NAR 2013. A loess curve is fit to the distribution of insertion indices, and
	sed to identify the minima between the 'essential' and 'non-essential'
	distributions. These distributions are then used to fit gamma distributions,
	which are then used to calculate log-odds ratios, which are used to determine
	an insertion-index threshold for gene essentiality. Note that this analysis
	requires a saturated mutant library, and is not suitable for the analysis of
	data sets with low insertion density. The script produces a number of
	diagnostic plots which can be used to verify that this condition has been
	met.\n"), fill=80)
	q(status=1)
}

STM_baseline <- read.table(input, sep="\t",header=TRUE,stringsAsFactors=F, quote="\"")

ii <- STM_baseline$ins_index

#identify second maxima
h <- hist(ii, breaks=200,plot=FALSE)
maxindex <- which.max(h$density[10:length(h$density)])
maxval <- h$mids[maxindex+3]

# print pdf of loess curve and later on, histogram
pdf(paste(input, "QC_and_changepoint_plots", "pdf", sep = "."))

#find inter-mode minima with loess
nG <- length(STM_baseline$read_count)
r <- floor(maxval *2000)
I = ii < r / 2000
h1 = hist(ii[I],breaks=(0:r/2000))
lo <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density
plot(h1$density, main="Density")
lines(predict(lo),col='red',lwd=2)
m = h1$mids[which.min(predict(lo))]
I1 = ((ii < m)&(ii >= 0))

h = hist(ii, breaks="FD",plot=FALSE)
I2 = ((ii >= m)&(ii < h$mids[max(which(h$counts>5))]))
f1 = (sum(I1) + sum(ii == 0))/nG
f2 = (sum(I2))/nG

d1 = fitdistr(ii[I1], "exponential")
d2 = fitdistr(ii[I2], "gamma") #fit curves

# print pdf of histogram
#pdf("Loess_and_changepoint_estimation.pdf")

#plots
hist(ii,breaks="FD", xlim=c(0,max(ii)), freq=FALSE,xlab="Insertion index", main="Gamma fits")
lines(0:200/500, f1*dgamma(0:200/500, 1, d1$estimate[1])) # was [2]
lines(0:200/500, f2*dgamma(0:200/500, d2$estimate[1], d2$estimate[2]))
# print changepoint

#calculate log-odds ratios to choose thresholds
lower <- max(which(log((pgamma(1:1000/10000, d2$e[1],d2$e[2])*
	(1-pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)))/
	(pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)*
	(1-pgamma(1:1000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
upper <- min(which(log((pgamma(1:1000/10000, d2$e[1],d2$e[2])*
	(1-pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)))/
	(pgamma(1:1000/10000, 1,d1$e[1], lower.tail=FALSE)*
	(1-pgamma(1:1000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))

essen <- lower/10000
ambig <- upper/10000

lines(c(lower/10000, lower/10000), c(0,20), col="red")
lines(c(upper/10000, upper/10000), c(0,20), col="red")

mtext(paste(essen, ":", "Essential changepoint"), side=3, adj=1, padj=2)
mtext(paste(ambig, ":", "Ambiguous changepoint"), side=3, adj=1, padj=3.75)
dev.off()

# read_index computation
STM_baseline$read_index <- STM_baseline$read_count / STM_baseline$gene_length

# score_E computation
l=c()
for (iindex in STM_baseline$ins_index) {
if (iindex<=essen) {l=c(l,iindex/essen)}
else if (iindex>essen & iindex<=ambig) {l=c(l,(iindex-essen)/(essen)+1)}
else if (iindex>ambig) {l=c(l,((iindex-ambig)/ambig)+2)}}
STM_baseline$score_E=l
STM_baseline$F_threshold=essen
STM_baseline$N_threshold=ambig

# Output writing
write.csv(STM_baseline, file=paste(input, "all", "csv", sep="."), row.names = FALSE, col.names= TRUE, quote=FALSE)
write.csv(STM_baseline[STM_baseline$ins_index < essen,], file=paste(input, "essen", "csv", sep="."), row.names = FALSE, col.names= TRUE, quote=FALSE)
write.csv(STM_baseline[STM_baseline$ins_index >= ambig,], file=paste(input, "nonessen", "csv", sep="."), row.names = FALSE, col.names= TRUE, quote=FALSE)
write.csv(STM_baseline[STM_baseline$ins_index >= essen & STM_baseline$ins_index < ambig,], file=paste(input, "ambig", "csv", sep="."), row.names = FALSE, col.names= TRUE, quote=FALSE)
