\name{RCA-package}
\alias{RCA-package}
\docType{package}
\title{
Relational Class Analysis
}
\description{
Relational Class Analysis (RCA) is a method for
detecting heterogeneity in attitudinal data. 
(as described in Goldberg A., 2011, Am. J. Soc, 116(5)) .
src also contains the source files of igraph-C 
which was written by Gabor.
}
\details{
\tabular{ll}{
Package: \tab RCA\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2013-05-15\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
RCA C Library: Amir Goldberg

igraph : Gabor Csardi

R Wrapper: Jinjian Zhai
}
\references{
Goldberg A., 2011, Am. J. Soc. 116(5) 
}
\keyword{ statistics, R }
\seealso{
igraph
}
\examples{
#attach("../data/SampleInputMatrix.rda")
#rcaout=RCA(SampleInputMatrix, 100, 0.075)
#rcaout$mod
#table(rcaout$member)
#rcaout$merge
#rcaout<-RCA(SampleInputMatrix,1000,0.25)
#rcaout$mod
#table(rcaout$member)
}
