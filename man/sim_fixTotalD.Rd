\name{sim_fixTotalD}
\alias{sim_fixTotalD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating MPRA data
}
\description{
This function simulates an MPRA dataset with specified input distribution and total depth across tags/barcodes.

}
\usage{
sim_fixTotalD(datMean, totalDepth = total, sigma2_DNA_a0, sigma2_DNA_a1, sigma2_RNA_a0, sigma2_RNA_a1, ntag = ntagsV, nsim = nelements, nrep = 5, measurementERROR = TRUE, slope)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datMean}{
A vector of numeric values indicating the input DNA count distribution, which is a pool of relative abundance of reads across tags. There is no requirement in length. The true input counts for each tag will be sampled from this vector (before scaling).
}
  \item{totalDepth}{
A numeric number indicating the total read depth over all tags.
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
  \item{ntag}{
An integer indicating the number of tags/barcodes for each oligonucleotide (oligos) or each allele.
}
  \item{nsim}{
An integer indicating the number of simulations or number of SNPs included in the dataset.
}
  \item{nrep}{
An integer indicating the number of replicates. 
}
  \item{measurementERROR}{
A logical value indicating whether there is randomness in the observed input counts for each tag/barcode. If FALSE, the true input count is used for all DNA replicates. Otherwise, TRUE (default). 
}
  \item{slope}{
A numeric vector of length 2*\code{nsim}*\code{ntag} specifying the transfection efficiencies across the oligos in the simulated data. The first half is for the reference allele, and the second half is for the mutant allele. To compute the power, the first half  would be different from the values in the second half of the vector (under the alternative). This is referred as 'b' in the paper. See examples below.
}

}
\value{
\item{datt}{A simulated data frame with ntag*nsim*2 number of rows and 2+nrep*2 number of columns. The first two columns are the allele and SNP name for each tag. The other columns are the generated DNA or RNA counts for the nrep replicates.}
}
\references{
Qiao, D., Zigler, C., Cho, M.H., Silverman, E.K., Zhou, X., et al. (2018). \emph{Statistical considerations for the analysis of massively parallel reporter assays data.}
}
\author{
Dandi Qiao
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
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
