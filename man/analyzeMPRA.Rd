\name{analyzeMPRA}
\alias{analyzeMPRA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyze MPRA data using specified tests
}
\description{
This function pre-processes the given MPRA data, and analyzes it using specified tests.
}
\usage{
analyzeMPRA(datt, nrepIn, rnaCol, nrepOut, nsim, ntag, method = c("MW", "Matching", "Adaptive", "Fisher", "QuASAR", "T-test", "mpralm", "edgeR", "DESeq2"), cutoff = -1, cutoffo = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datt}{
A data frame containing the MPRA dataset. It should have nsim*ntag*2 rows and 2+nrepIn+nrepOut columns. The first column should be named 'allele', and the second column should be named 'simN'. The 'allele' columns should contain only two possible values 'Ref' and 'Mut' to refer to the two versions of alleles for each SNP.
}
  \item{nrepIn}{
An integer indicating the number of DNA replicates.
}
  \item{rnaCol}{
An integer indicating the starting column of the RNA replicates in \code{datt}.
}
  \item{nrepOut}{
An integer indicating the number of RNA replicates.
}
  \item{nsim}{
An integer indicating the number of SNPs/comparisons in the MPRA data. A comparison refer to the unit with two alleles for testing allele-specific expression.
}
  \item{ntag}{
An integer indicating the number of tags/barcodes for each allele.
}
  \item{method}{
A vector of characters specifying the tests to be used. The possible options are: MW, Matching, Adaptive, Fisher, QuASAR, T-test, mpralm, edgeR and DESeq2.

}
  \item{cutoff}{
A numeric or integer value. Tags with DNA count less than or equal to \code{cutoff} in any of the DNA replicates will be removed.
}
  \item{cutoffo}{
A numeric or integer value. Tags with RNA count less than or equal to \code{cutoffo} in any of the RNA replicates will be removed.
}
}
\details{
This function first normalizes the replicates using the maximum depth across all replicates, and filters the tags according to the specified cutoffs. Then it analyzes the processed MPRA data using the specified tests.
}
\value{
\item{results}{The actual p-values of all the SNPs for the specified tests using the given dataset.}
}
\references{
Qiao, D., Zigler, C., Cho, M.H., Silverman, E.K., Zhou, X., et al. (2018). \emph{Statistical considerations for the analysis of massively parallel reporter assays data.}
}
\author{
Dandi Qiao
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[atMPRA:getPower]{atMPRA}}
\code{\link[atMPRA:sim_fixTotalD]{atMPRA}}
\code{\link[atMPRA:sim_fixMeanD]{atMPRA}}
\code{\link[atMPRA:sim_fixInputDist]{atMPRA}}
}
\examples{
nsim = 10
ntag = 10
slope=c(rep(1, ntag*nsim), rep(1.5, ntag*nsim)) 
nrep=5
data(datMean)
simData = sim_fixTotalD(datMean=datMean, totalDepth=200000, sigma2_DNA_a0=0.001, sigma2_DNA_a1=0.23, sigma2_RNA_a0=0.18, sigma2_RNA_a1=35, ntag=ntag, nsim=nsim, nrep=nrep, slope=slope)
results = analyzeMPRA(simData, nrep, rnaCol=2+nrep+1, nrep, nsim, ntag, method=c("MW","mpralm", "edgeR", "DESeq2"), cutoff=0, cutoffo=0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MPRA}

