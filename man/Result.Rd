\name{Result}
\alias{Result}
\docType{data}
\title{
Relational Class Analysis for R
}
\description{
Relational Class Analysis (RCA) is a method for detecting
 heterogeneity in attitudinal data. (as described in 
Goldberg A., 2011, Am. J. Soc, 116(5)). src/igraph-R 
contains the header files of igraph-C which was written 
by Gabor.
}
\usage{data(Result)}
\format{
  The format is:
List of 3
 $ member: int [1:1551] 0 0 0 1 0 0 2 0 1 2 ...
 $ mod   : num 0.299
 $ merge : num [1:2, 1:2] 0 3 2 1
}
\details{
Results include:

a. [your_returned_data]$member: Assignment vector (size 1xN)
b. [your_returned_data]$mod: Modularity (double)
c. [your_returned_data]$merge: Merge tree (array of varying size)
d. [your_returned_data]$stats: Stats (2 doubles)
}
\examples{
## data(Result)
}
\keyword{datasets}
