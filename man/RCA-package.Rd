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
src/igraph-R contains the header files of igraph-C 
which was written by Gabor.
}
\details{
\tabular{ll}{
Package: \tab RCA\cr
Type: \tab Package\cr
Version: \tab 1.3.3\cr
Date: \tab 2013-03-28\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Amir Goldberg
igraph header files: Gabor Csardi
Maintainer: Jinjian Zhai
}
\references{
The src contains some header files of igraph-C package. 
}
\keyword{ statistics, R }
\seealso{
igraph
}
\examples{
## attach("RCA/data/SampleInputMatrix.rda")
## Result=RCA(SampleInputMatrix, 100, 0.075)
}
