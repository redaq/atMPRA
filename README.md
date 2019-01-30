# @MPRA
Analysis Toolset for Massively Parallel Reporter Assay data

A newer version of this package has been uploaded, though we are still actively working on the functions.

This R package is under development. Only the Beta version is available for testing purposes. 


# Installation

This package is dependent on the R packages: 

coin, MatchIt (Available on CRAN)

DESeq2,   mpra, edgeR (Avaiable on Bioconductor)

QuASAR (Available here: https://github.com/piquelab/QuASAR)


After the packages above are successfully installed, please download the atMPRA_1.3.tar.gz file.

In R, type: 
```r
install.packages("atMPRA_1.3.tar.gz", repos=NULL, type="source")

library(atMPRA)
```

More tutorials can be found in the [vignette](https://github.com/redaq/atMPRA/blob/master/atMPRA_vignettes.pdf). 

# Quick example
```r
nsim=10

ntag=10

result = getPower(nsim=nsim, ntag=ntag, nrepIn=3, nrepOut=3, slope = c(rep(1, ntag*nsim), rep(2,ntag*nsim)), method=c("MW", "mpralm"), scenario="fixTotalDepth")
```

### Simulate data that resembles the GSE70531
```r
data(GSE70531_params) 

totalDepth = 200000

ntag= 10

nsim= 10

nrepIn=5

nrepOut = 5

inputProp = GSE70531_params[[3]](runif(ntag*nsim*2))

slopel = GSE70531_params[[4]](runif(nsim*2))

inputDispFunc=GSE70531_params[[1]]

outputDispFunc=GSE70531_params[[2]]

slope = rep(slopel, each=ntag)

datt=sim_fixDepth(inputProp, ntag, nsim, nrepIn,  nrepOut, slope, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc, sampleDepth=totalDepth) 
```

### Estimate model parameters for the simulated data
```r
new_params=estimateMPRA(datt, nrepIn, rnaCol=8, nrepOut, nsim, ntag)

datt=sim_fixDepth(inputProp=new_params[[3]](runif(ntag*nsim*2)), ntag, nsim, nrepIn,  nrepOut, slope, inputDispFunc=new_params[[1]], outputDispFunc=new_params[[2]], sampleDepth=totalDepth) 

```

### Test allele-specific effects using specified methods
```r
results = analyzeMPRA(datt, nrepIn, rnaCol, nrepOut, nsim, ntag, method=c("MW", "Adaptive", "QuASAR", "T-test", "mpralm", "DESeq2"), cutoff=0, cutoffo=0)
```

### Simulate data with fixed mean tag counts for the two alleles
```r
inputDist= GSE70531_params[[3]](runif(nsim*ntag*2))

datt=sim_fixInputMean(mean_A=10, mean_B=100,  ntag=ntag, nsim=nsim, nrepIn=nrepIn, nrepOut=nrepOut, slope=slope, inputDist=inputDist, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc)
```

### Power calculation
```r
nrepIn=2

nrepOut = 2

slopel = GSE70531_params[[4]](runif(nsim))

slopel = c(slopel, slopel+1)

slope = rep(slopel, each=ntag)

result2 = getPower(nsim, ntag, nrepIn, nrepOut, slope, scenario="fixInputDist", method=c("MW","T-test", "mpralm", "edgeR", "DESeq2"), fixInput  = c(20, 100), inputDist=inputDist, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc,  cutoff=-1, cutoffo=-1)

result3 = getPower(nsim, ntag, nrepIn, nrepOut, slope=1, scenario="fixTotalDepth", method=c("MW", "Matching", "Adaptive", "Fisher", "QuASAR", "T-test", "mpralm", "edgeR", "DESeq2"), fixTotalD= 200000, inputDist=inputDist,inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc,  cutoff=-1, cutoffo=-1, p.adjust.method="none")

```





