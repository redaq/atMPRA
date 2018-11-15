\name{getPower}
\alias{getPower}
\title{
Power calculation for designing MPRA experiments
}
\description{
This function computes the power (or type I error if under the null) of the specified tests using simulated MPRA data with user-specified parameters for the MPRA experiment.
}
\usage{


getPower(nsim = 100, ntag = 10, nrepIn = 2, nrepOut=2, slope = 1, scenario=c("fixInputDist", "fixTotalDepth", "fixMeanDepth"), method=c("MW", "Matching", "Adaptive", "Fisher", "QuASAR", "T-test", "mpralm", "edgeR", "DESeq2"), fixInput  = c(20, 100), fixMeanD = 70 , fixTotalD= 20000000, std_A=mean_A, std_B=mean_B, inputDist=NULL, inputDispFunc=NULL, outputDispFunc=NULL, inputDispParam=NULL, outputDispParam=NULL, cutoff=-1, cutoffo=-1)
}
\arguments{
  \item{nsim}{
An integer indicating the number of simulations or number of SNPs included in the dataset
}
  \item{ntag}{
An integer indicating the number of tags/barcodes for each oligonucleotide (oligos) or each allele
}
  \item{nrepIn}{
An integer indicating the number of DNA replicates 
}
  \item{nrepOut}{
An integer indicating the number of RNA replicates 
}
  \item{slope}{
A numeric vector of length 1 or 2*\code{nsim}*\code{ntag} specifying the transfection efficiencies across the oligos in the simulated data. The first half is for the reference allele, and the second half is for the mutant allele. To compute the power, the first half  would be different from the values in the second half of the vector (under the alternative). This is referred as 'b' in the paper. See examples below.
}
  \item{scenario}{
Three different simulation scenarios are included. This includes \\
\code{fixInputDist} - fix the mean and sd of the true input counts across the tags for the two alleles. The same input distribution is used for all simulations/SNPs;
\code{fixTotalDepth} - Sample the true input counts from the default/given distribution, and scale to the specified total depth.
\code{fixMeanDepth} - sample the true input counts from the default/given distrubtion, and scale it so the mean depth across all tags equals the specified mean depth.
}
  \item{method}{
Accepts a vector of characters specifying the tests to be used. The possible options are: MW, Matching, Adaptive, Fisher, QuASAR, T-test, mpralm, edgeR and DESeq2.
}
  \item{fixInput}{
A vector of two numeric values specifying the mean of input count for the reference and alternative alleles for all SNPs.
}
 
  \item{fixMeanD}{
If scenario is \code{fixTotalDepth}, a numeric number is expected, which is the total reads across all tags.
}
  \item{fixTotalD}{
If scenario is \code{fixMeanDepth}, a numeric value is expected, which is the mean reads across all tags. 
}

  \item{std_A}{
An optional parameter specifying the standard devaition of the mean input counts across tags for the reference allele.
}
  \item{std_B}{
An optional parameter specifying the standard devaition of the mean input counts across tags for the mutant allele.

}

  \item{inputDist}{
A vector of numerical proportions that can be used to generate the DNA proportions across the tags.
}
  \item{inputDispFunc}{
Optional parameter that provides a dispersion function estimated for the input replicates using DESeq2.
}
  \item{outputDispFunc}{
Optional parameter that provides a dispersion function estimated for the output replicates using DESeq2.
}
  \item{inputDispParam}{
This parameter is required if inputDispFunc is not provided. It should give three parameter estimates for the dispersion function of the DNA input replicates. The three parameters correspond to a0, a1, and d2, which specify that the dispersion parameter is a lognormal distribution with mean log(a0+a1/mu) and sd d2, where mu is the mean of DNA count across the replicates.
}
  \item{outputDispParam}{
This parameter is required if outputDispFunc is not provided. It should give three parameter estimates for the dispersion function of the RNA output replicates. The three parameters correspond to a0, a1, and d2, which specify that the dispersion parameter is a lognormal distribution with mean log(a0+a1/mu) and sd d2, where mu is the mean of RNA count across the replicates.
}

  \item{cutoff}{
A numeric or integer value. Tags with DNA count less than or equal to \code{cutoff} in any of the DNA replicates will be removed.
}
  \item{cutoffo}{
A numeric or integer value. Tags with RNA count less than or equal to \code{cutoffo} in any of the RNA replicates will be removed.
}
}
\details{
This function simulates MPRA data according to user-specified parameters, and computes the power/type I error of the tests requested using the simulated data. Before analyzing the simulated data, normalization and filtering is performed.
}
\value{
\item{Power}{The power/type I error of the specified tests.}
\item{simData}{The simulated data frame.}
\item{results}{The actual p-values of all the SNPs for the specified tests in the simulated data.}
}
\references{
Qiao, D., Zigler, C., Cho, M.H., Silverman, E.K., Zhou, X., et al. (2018). \emph{Statistical considerations for the analysis of massively parallel reporter assays data.}

}
\author{
Dandi Qiao
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[atMPRA:analyzeMPRA]{atMPRA}}
}
\examples{
data(GSE70531_params) 
ntag= 10
nsim= 10
nrepIn=2
nrepOut = 2
inputDist= GSE70531_params[[3]](runif(nsim*ntag*2))
inputDispFunc=GSE70531_params[[1]]
outputDispFunc=GSE70531_params[[2]]
result = getPower(nsim, ntag, nrepIn, nrepOut, slope = 1, scenario="fixInputDist", method=c("MW","T-test", "mpralm", "edgeR", "DESeq2"), fixInput  = c(20, 100), inputDist=inputDist, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc,  cutoff=-1, cutoffo=-1)

}
\keyword{power}
\keyword{MPRA}
