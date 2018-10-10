\name{getPower}
\alias{getPower}
\title{
Power calculation for designing MPRA experiments
}
\description{
This function computes the power (or type I error if under the null) of the specified tests using simulated MPRA data with user-specified parameters for the MPRA experiment.
}
\usage{
getPower(nsim = 100, ntag = 10, nrep = 1, slope = NA, scenario = c("fixInputDist", "fixTotalDepth", "fixMeanDepth"), method = c("MW", "Matching", "Adaptive", "Fisher", "QuASAR", "T-test", "mpralm", "edgeR", "DESeq2"), sigma2_DNA_a0, sigma2_DNA_a1, sigma2_RNA_a0, sigma2_RNA_a1, fixInput = NA, fixMeanD = 70, fixTotalD = 2e+07, datMean = NA, cutoff = -1, cutoffo = -1)
}
\arguments{
  \item{nsim}{
An integer indicating the number of simulations or number of SNPs included in the dataset
}
  \item{ntag}{
An integer indicating the number of tags/barcodes for each oligonucleotide (oligos) or each allele
}
  \item{nrep}{
An integer indicating the number of replicates 
}
  \item{slope}{
A numeric vector of length 2*\code{nsim}*\code{ntag} specifying the transfection efficiencies across the oligos in the simulated data. The first half is for the reference allele, and the second half is for the mutant allele. To compute the power, the first half  would be different from the values in the second half of the vector (under the alternative). This is referred as 'b' in the paper. See examples below.
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
  \item{sigma2_DNA_a0}{
This parameter and the \code{sigma2_DNA_a1} parameter specify the relationship between the dispersion parameter and the mean count across replicates for the DNA counts in the simulation. Given the true input count for the tag, the observed input counts across the replicates were generated using a negative binomial distribution with mean equal to the true input count, and dispersion parameter = a0 + a1/mean.
}
  \item{sigma2_DNA_a1}{
This parameter and the \code{sigma2_DNA_a0} parameter specify the relationship between the dispersion parameter and the mean count across replicates for the DNA counts in the simulation. Given the true input count for the tag, the observed input counts for this tag across the replicates were generated using a negative binomial distribution with mean equal to the true input count, and the dispersion parameter equals to a0 + a1/mean.

}
  \item{sigma2_RNA_a0}{
This parameter and the \code{sigma2_RNA_a1} parameter specify the relationship between the dispersion parameter and the mean count across replicates for the RNA counts in the simulation. Given the true input count for the tag, the observed RNA counts across the replicates were generated using a negative binomial distribution with mean equal to the transfection efficiency times the true input count, and dispersion parameter = a0 + a1/mean.
}
  \item{sigma2_RNA_a1}{
This parameter and the \code{sigma2_RNA_a0} parameter specify the relationship between the dispersion parameter and the mean count across replicates for the RNA counts in the simulation. Given the true input count for the tag, the observed RNA counts across the replicates were generated using a negative binomial distribution with mean equal to the transfection efficiency (as specified using slope) times the true input count, and dispersion parameter = a0 + a1/mean.
}
  \item{fixInput}{
If scenario is \code{fixInputDist}, a numeric vector of length 4 is expected. The four values correspond to the mean and sd of the true input counts for the reference allele, and the mean and sd of the true input counts for the mutant allele.
}
  \item{fixMeanD}{
If scenario is \code{fixTotalDepth}, a numeric number is expected, which is the total reads across all tags.
}
  \item{fixTotalD}{
If scenario is \code{fixMeanDepth}, a numeric value is expected, which is the mean reads across all tags. 
}
  \item{datMean}{
If scenario is specified to be \code{fixTotalDepth} or \code{fixMeanDepth}, a vector of numeric values is expected, otherwise it can be left as NA. There is no requirement in length. The true input counts for each tag will be sampled from this vector (before scaling).
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
nsim = 10
ntag = 10
slope=c(rep(1, ntag*nsim), rep(1.5, ntag*nsim)) 
results = getPower(nsim = nsim, ntag = ntag, nrep = 1, slope=slope , scenario=c("fixInputDist"), method="MW", sigma2_DNA_a0=0.001, sigma2_DNA_a1=0.23, sigma2_RNA_a0=0.18, sigma2_RNA_a1=35, fixInput  = c(20, 120, 20, 120))

}
\keyword{power}
\keyword{MPRA}
