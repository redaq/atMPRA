\name{calEnrichment}
\alias{calEnrichment}
\title{
Compute bias towards SNPs with allelic DNA imbalance}
\description{
This function uses hypergeometric tests to compute enrichment in SNPs with significant allelic DNA imbalance among the SNPs with significant allele-specific effect defined by each method.
}
\usage{
calEnrichment(resultAll, threshold = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultAll}{
A data frame containing the p-values for all SNPs using a number of analyses methods. This can be generated using /code{analyzeMPRA}. The first column is named 'simN' which is the SNP name, the second column is the p-value indicating whether there is allelic DNA-imbalance for the SNP. The other columns are the p-values for testing allele-specific effects using different methods.
}
  \item{threshold}{
The significance level used. 0.05 is the default.
}
}
\details{
This method to used to check for bias in the specified tests caused by allelic DNA imbalance, i.e. difference in the DNA count distributions between the two alleles.
}
\value{
\item{enrichmen}{A data frame containing the enrichment p-values for each method. Column 'q' refers to the number of significant SNPs with allelic DNA imbalance. Column 'm' refers to the number of significant tests with allele-specific effect defined by each method in each row. 'enrichP' is the p-value from the hypergeometric test. Only SNPs with valid test p-values were considered here.}
}
\references{
Qiao, D., Zigler, C., Cho, M.H., Silverman, E.K., Zhou, X., et al. (2018). \emph{Statistical considerations for the analysis of massively parallel reporter assays data.}

}
\author{
Dandi Qiao
}

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
calEnrichment(results)

}
