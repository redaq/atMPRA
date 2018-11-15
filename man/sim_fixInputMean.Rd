\name{sim_fixInputMean}
\alias{sim_fixInputMean}
\title{
Simulating MPRA data
}
\description{
This function simulates an MPRA dataset with specified input distribution and mean depth across tags/barcodes.
}
\usage{
sim_fixInputMean(mean_A, mean_B,  ntag, nsim, nrepIn, nrepOut, slope, inputDist=NULL, std_A=mean_A, std_B=mean_B, inputDispFunc=NULL, outputDispFunc=NULL, inputDispParam=NULL, outputDispParam=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{mean_A}{
The mean of the true input counts across tags for the reference allele.
}
  \item{mean_B}{
The mean of the true input counts across tags for the mutant allele.

}
  \item{ntag}{
An integer indicating the number of tags/barcodes for each oligonucleotide (oligos) or each allele.
}
  \item{nsim}{
An integer indicating the number of simulations or number of SNPs included in the dataset.
}
  \item{nrepIn}{
An integer indicating the number of DNA replicates.
}
  \item{nrepOut}{
An integer indicating the number of RNA output replicates.
}


  \item{inputDist}{
This parameter is required if std_A and std_B were not specified. This should be a vector of proportions that can be used to sample the input proportions across tags.
}
  \item{std_A}{
An optional parameter specifying the standard devaition of the mean input counts across tags for the reference allele.
}
  \item{std_B}{
An optional parameter specifying the standard devaition of the mean input counts across tags for the mutant allele.

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


}
\value{
\item{datt}{A simulated data frame with ntag*nsim*2 number of rows and 2+nrep*2 number of columns. The first two columns are the allele and SNP name for each tag. The other columns are the generated DNA or RNA counts for the nrep replicates.
}}
\references{
Qiao, D., Zigler, C., Cho, M.H., Silverman, E.K., Zhou, X., et al. (2018). \emph{Statistical considerations for the analysis of massively parallel reporter assays data.}
}
\author{
Dandi Qiao
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[atMPRA:sim_fixTotalD]{atMPRA}}
\code{\link[atMPRA:sim_fixInputDist]{atMPRA}}
}
\examples{
data(GSE70531_params) 
ntag= 10
nsim= 10
nrepIn=5
nrepOut = 5
slopel = getParam[[4]](runif(nsim*2))
slope = rep(slopel, each=ntag)
inputDist= getParam[[3]](runif(nsim*ntag*2))
inputDispFunc=getParam[[1]]
outputDispFunc=getParam[[2]]

datt=sim_fixInputMean(mean_A=10, mean_B=100,  ntag=ntag, nsim=nsim, nrepIn=nrepIn, nrepOut=nrepOut, slope=slope, inputDist=inputDist, inputDispFunc=inputDispFunc, outputDispFunc=outputDispFunc)

}
