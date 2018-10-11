# @MPRA
Analysis Tools for Massively Parallel Reporter Assay data

This R package is under development. Only the Beta version is available for testing purposes. 

To run this package, please install using the source file:

install.packages("atMPRA_1.0.tar.gz", repos=NULL, type="SOURCE")

library(atMPRA)

## An example: power calculation using simulated data and specified tests
nsim = 10 

ntag = 10

slope=c(rep(1, ntag*nsim), rep(1.5, ntag*nsim))

data(datMean)

results = getPower(nsim = nsim, ntag = ntag, nrep = 1, slope , scenario=c("fixTotalDepth"), method=c("MW", "Matching", "Adaptive", "Fisher", "QuASAR"), sigma2_DNA_a0=0.001, sigma2_DNA_a1=0.23, sigma2_RNA_a0=0.18, sigma2_RNA_a1=35,  fixTotalD=200000, datMean=datMean, cutoff=-1, cutoffo=-1)

## Another example: simulating and analyzing MPRA data using specified tests
nrep=5

simData = sim_fixTotalD(datMean=datMean, totalDepth=200000, sigma2_DNA_a0=0.001, sigma2_DNA_a1=0.23, sigma2_RNA_a0=0.18, sigma2_RNA_a1=35, ntag=ntag, nsim=nsim, nrep=nrep, slope=slope)

results2 = analyzeMPRA(simData, nrep, rnaCol=2+nrep+1, nrep, nsim, ntag, method=c("MW","mpralm", "edgeR", "DESeq2"), cutoff=0, cutoffo=0)

calEnrichment(results2)
