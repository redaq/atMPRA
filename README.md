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

results = getPower(nsim = nsim, ntag = ntag, nrep = 1, slope=slope , scenario=c("fixInputDist"), method="MW", sigma2_DNA_a0=0.001, 
sigma2_DNA_a1=0.23, sigma2_RNA_a0=0.18, sigma2_RNA_a1=35, fixInput  = c(20, 120, 20, 120))

## Another example: simulating and analyzing MPRA data using specified tests
nrep=5

data(datMean)

simData = sim_fixTotalD(datMean=datMean, totalDepth=200000, sigma2_DNA_a0=0.001, sigma2_DNA_a1=0.23, sigma2_RNA_a0=0.18, sigma2_RNA_a1=35, ntag=ntag, nsim=nsim, nrep=nrep, slope=slope)

results2 = analyzeMPRA(simData, nrep, rnaCol=2+nrep+1, nrep, nsim, ntag, method=c("MW","mpralm", "edgeR", "DESeq2"), cutoff=0, cutoffo=0)


