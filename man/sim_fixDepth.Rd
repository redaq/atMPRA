\name{sim_fixDepth}
\alias{sim_fixDepth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating MPRA data
}
\description{
This function simulates an MPRA dataset with specified input distribution and total depth across tags/barcodes.

}
\usage{
sim_fixDepth(inputProp, ntag, nsim, nrepIn,  nrepOut, slope, inputDispFunc=NULL, outputDispFunc=NULL, sampleDepth=NULL, inputDispParam=NULL, outputDispParam=NULL,  meanDepth=NULL) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputProp}{
A vector of numeric values indicating the input DNA count distribution, which is a pool of relative abundance of reads across tags. It should include 2*ntag*nsim proportions to indicate the proportion of reads correspond to each tag.
}
  \item{ntag}{
An integer indicating the number of tags/barcodes for each oligonucleotide (oligos) or each allele.
}
  \item{nsim}{
An integer indicating the number of simulations or number of SNPs included in the dataset.
}
  \item{nrepIn}{
An integer indicating the number of replicates for the DNA input.
}

  \item{nrepOut}{
An integer indicating the number of replicates for the RNA output.
}

  \item{slope}{
A numeric vector of length 2*\code{nsim}*\code{ntag} specifying the transfection efficiencies across the oligos in the simulated data. The first half is for the reference allele, and the second half is for the mutant allele. To compute the power, the first half  would be different from the values in the second half of the vector (under the alternative). This is referred as 'b' in the paper. See examples below.
}
  \item{inputDispFunc}{
Optional parameter that provides a dispersion function estimated for the input replicates using DESeq2. 
}
  \item{outputDispFunc}{
Optional parameter that provides a dispersion function estimated for the output replicates using DESeq2. 
}
  \item{sampleDepth}{
An integer vector specifying the total read depth over all tags. It could be of length 1 or length nrepIn+nrepOut. If it is of length 1, the same total depth is used for all DNA and RNA replicates. If this is specified, values for \code{meanDepth} is ignored.
}
  \item{inputDispParam}{
This parameter is required if inputDispFunc is not provided. It should give three parameter estimates for the dispersion function of the DNA input replicates. The three parameters correspond to a0, a1, and d2, which specify that the dispersion parameter is a lognormal distribution with mean log(a0+a1/mu) and sd d2, where mu is the mean of DNA count across the replicates.
}
  \item{outputDispParam}{
This parameter is required if outputDispFunc is not provided. It should give three parameter estimates for the dispersion function of the RNA output replicates. The three parameters correspond to a0, a1, and d2, which specify that the dispersion parameter is a lognormal distribution with mean log(a0+a1/mu) and sd d2, where mu is the mean of RNA count across the replicates.
}

  \item{meanDepth}{
An integer vector specifying the mean read depth over all tags. It could be of length 1 or length nrepIn+nrepOut. If it is of length 1, the same mean depth is used for all DNA and RNA replicates.
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
data(GSE70531_params) 
inputDispFunc=getParam[[1]]
outputDispFunc=getParam[[2]]
totalDepth = 200000
ntag= 10
nsim= 10
nrepIn=5
nrepOut = 5
inputProp = getParam[[3]](runif(ntag*nsim*2))
slopel = getParam[[4]](runif(nsim*2))
slope = rep(slopel, each=ntag)

datt=sim_fixDepth(inputProp, ntag, nsim, nrepIn,  nrepOut, slope, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc, sampleDepth=totalDepth) 
}
