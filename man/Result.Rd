\name{Result}
\alias{Result}
\docType{data}
\title{
Relational Class Analysis for R
}
\description{
Relational Class Analysis (RCA) is a method for detecting
 heterogeneity in attitudinal data. (as described in 
Goldberg A., 2011, Am. J. Soc, 116(5)). src also contains the files of igraph-C which was written 
by Gabor.
}
\format{
  The format is:

List of 3

 $ member: int [1:1551] 0 0 0 1 0 0 2 0 1 2 ...

 $ mod   : num 0.299

 $ merge : num [1:2, 1:2] 0 3 2 1
}
\details{
Results include:

(a).[your_returned_data]$member: Assignment vector (size 1xN)

Each cell in the vector includes the group number the corresponding 
observation was assigned to. Group numbers range from 0 to G-1.

(b)[your_returned_data]$mod: Modularity (double)

A number ranging from 0 to 1 indicating the strength of the partition 
into groups.

(c)[your_returned_data]$merge:

Cluster tree, reported as a merge matrix (array of varying size). 
See igraph/leading.eigenvector.community for more information.

}
\keyword{datasets}
