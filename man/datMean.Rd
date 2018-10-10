\name{datMean}
\alias{datMean}
\docType{data}
\title{
An example of input DNA count distribution across tags.
}
\description{
An example of input DNA count distribution across tags. Generated from the first DNA replicate in the GSE70531 MPRA data described in Ulirsche et al, 2016. 
}
\usage{data("datMean")}
\format{
  The format is:
 num [1:231476] 1.86e-06 3.28e-06 3.22e-06 6.44e-06 2.03e-06 ...
}
\details{
The first DNA replicate in the GSE70531 data (column name 'K562_minP_DNA1' in the raw data file). The counts were normalized using the sum of reads in this replicate. Counts 0 were set to be 1e-10 to avoid NA counts to be generated in the simulation.
}
\source{
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70531}
}
\references{
Ulirsch JC, Nandakumar SK, Wang L, Giani FC et al. \emph{Systematic Functional Dissection of Common Genetic Variation Affecting Red Blood Cell Traits.} Cell 2016 Jun 2;165(6):1530-1545. PMID: 27259154
}
\examples{
data(datMean)
}
\keyword{datasets}
