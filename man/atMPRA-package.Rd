\name{Analysis Toolset for MPRA Data}
\alias{atMPRA-package}
\alias{atMPRA}
\docType{package}
\title{
\packageTitle{atMPRA}
}
\description{
\packageDescription{atMPRA}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{atMPRA}
\packageIndices{atMPRA}
This package implements the adaptive test and a number of other methods for analyzing MPRA data, and also functions of power calculations for the design of MPRA experiments.
The main functions are:
\code{getPower}: function for power calculation using simulated data and specified test
\code{sim_fixInputDist}: function for simulating MPRA data with specified mean and standard deviation for the input distributions 
\code{sim_fixTotalD}: function for simulating MPRA data with specified input distribution and total depth
\code{sim_fixMeanD}: function for simulating MPRA data with specified input distribution and mean depth
\code{analyzeMPRA}: function for analyzing MPRA data (either real or simulated) using specified test

}
\author{
\packageAuthor{atMPRA}

Maintainer: \packageMaintainer{atMPRA}
}
\references{
Qiao, D., Zigler, C., Cho, M.H., Silverman, E.K., Zhou, X., et al. (2018). \emph{Statistical considerations for the analysis of massively parallel reporter assays data.}
}
\keyword{MPRA}
\seealso{
\code{\link[atMPRA:atMPRA-package]{atMPRA}}
}
\examples{
nsim = 10
ntag = 10
slope=c(rep(1, ntag*nsim), rep(1.5, ntag*nsim)) 
results = getPower(nsim = nsim, ntag = ntag, nrep = 1, slope=slope , scenario=c("fixInputDist"), method="MW", sigma2_DNA_a0=0.001, sigma2_DNA_a1=0.23, sigma2_RNA_a0=0.18, sigma2_RNA_a1=35, fixInput  = c(20, 120, 20, 120))

nrep=5
data(datMean)
simData = sim_fixTotalD(datMean=datMean, totalDepth=200000, sigma2_DNA_a0=0.001, sigma2_DNA_a1=0.23, sigma2_RNA_a0=0.18, sigma2_RNA_a1=35, ntag=ntag, nsim=nsim, nrep=nrep, slope=slope)

results2 = analyzeMPRA(simData, nrep, rnaCol=2+nrep+1, nrep, nsim, ntag, method=c("MW","mpralm", "edgeR", "DESeq2"), cutoff=0, cutoffo=0)

}
