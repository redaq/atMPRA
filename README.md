# @MPRA
Analysis Tools for Massively Parallel Reporter Assay data

A newer version of this package has been uploaded, though we are still actively working on the functions.

This R package is under development. Only the Beta version is available for testing purposes. 


# Installation

Please download the atMPRA_1.2.tar.gz file.

In R, type: 
install.packages("atMPRA_1.2.tar.gz", repos=NULL, type="source")

library(atMPRA)


# Examples
data(GSE70531_params) 

inputDispFunc=GSE70531_params[[1]]

outputDispFunc=GSE70531_params[[2]]

totalDepth = 200000

ntag= 10

nsim= 10

nrepIn=5

nrepOut = 5

inputProp = GSE70531_params[[3]](runif(ntag*nsim*2))

slopel = GSE70531_params[[4]](runif(nsim*2))

slope = rep(slopel, each=ntag)

### simulate data that resembles the GSE70531
datt=sim_fixDepth(inputProp, ntag, nsim, nrepIn,  nrepOut, slope, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc, sampleDepth=totalDepth) 

rnaCol=8

### estimate model parameters for the simulated data
result=estimateMPRA(datt, nrepIn, rnaCol, nrepOut, nsim, ntag)

### test allele-specific effects using specified methods
results2 = analyzeMPRA(datt, nrepIn, rnaCol, nrepOut, nsim, ntag, method=c("MW", "Matching", "Adaptive", "Fisher", "QuASAR", "T-test", "mpralm", "edgeR", "DESeq2"), cutoff=0, cutoffo=0)

### simulate data with fixed mean tag counts for the two alleles
inputDist= GSE70531_params[[3]](runif(nsim*ntag*2))

datt=sim_fixInputMean(mean_A=10, mean_B=100,  ntag=ntag, nsim=nsim, nrepIn=nrepIn, nrepOut=nrepOut, slope=slope, inputDist=inputDist, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc)

### Power calculation
nrepIn=2

nrepOut = 2

slopel = GSE70531_params[[4]](runif(nsim))

slopel = c(slopel, slopel+1)

slope = rep(slopel, each=ntag)

result3 = getPower(nsim, ntag, nrepIn, nrepOut, slope, scenario="fixInputDist", method=c("MW","T-test", "mpralm", "edgeR", "DESeq2"), fixInput  = c(20, 100), inputDist=inputDist, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc,  cutoff=-1, cutoffo=-1)

result4 = getPower(nsim, ntag, nrepIn, nrepOut, slope=1, scenario="fixTotalDepth", method=c("MW", "Matching", "Adaptive", "Fisher", "QuASAR", "T-test", "mpralm", "edgeR", "DESeq2"), fixTotalD= 200000, inputDist=inputDist,inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc,  cutoff=-1, cutoffo=-1)







