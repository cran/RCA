\name{RCA}
\alias{RCA}
%- Also NEED an '\alias' for EACH other topic documented here.
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
\usage{
RCA(inputMatrix, bootstrap = 1000, p_value = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputMatrix}{
    Input matrix. A dataset of size N x K. N: observation. K: variables.

}
  \item{bootstrap}{
    Bootstrap. Default is 1000.
}
  \item{p_value}{
    p_value. Default is 0.05.
}
}
\details{
%This is a R function using the C implementation of Relational Class Analysis (RCA).

Let X be a dataset of size N x K. X must not include missing data.
RCA finds an optimal division of X into G groups, such that each group
of observations follows a distinctive pattern of relationships between the
K variables. Each group comprises an exclusive subset of X, such that the K variables within each group have a unique covariance structure. 

%Group assignment is reported as a vector M. Modulairty, Q, measures the strength of the division into groups.

RCA divides X into a variable number of G groups, assigning each observation in X to one group. It reports a measure called modularity, which ranges from 0 to 1 and which measures the strength of the division into groups. It also reports a cluster tree which traces the process by which the data were partitioned into G groups. 

%RCA-R is a C-based package that implements the RCA algorithm and generates a simple-to-use command line executable program. You will need to follow the steps below in order to compile the RCA-R code and generate a local copy of the executable.

}
\value{
Results include:

\item{}{[your_returned_data]$member: Assignment vector (size 1xN)

Each cell in the vector includes the group number the corresponding observation 

was assigned to. Group numbers range from 0 to G-1. 

}

\item{}{[your_returned_data]$mod: Modularity (double)

A number ranging from 0 to 1 indicating the strength of the partition into groups. 

}

\item{}{[your_returned_data]$merge: 

Cluster tree, reported as a merge matrix (array of varying size). 

See igraph/leading.eigenvector.community for more information. 

}

}
\references{
Goldberg, A., AJS 116(5): 1397-1436 (http://http://www.jstor.org/stable/10.1086/657976)}
\author{
Amir Goldberg, Gabor Csardi, Jinjian Zhai
}
\note{
Depended on igraph library.
}


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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RCA , statistics, R}
