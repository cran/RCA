\name{SampleInputMatrix}
\alias{SampleInputMatrix}
\docType{data}
\title{
Relational Class Analysis
%%   ~~ data name/kind ... ~~
}
\description{
Sample input matrix.
}
\usage{data(SampleInputMatrix)}
\format{
  The format is:
 int [1:1551, 1:17] 4 2 2 3 2 1 1 1 3 2 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:17] "V1" "V2" "V3" "V4" ...
}
\details{
Let X be a dataset of size N x K. X must not include missing data.
RCA finds an optimal division of X into G groups, such that each group
of observations follows a distinctive pattern of raltionships between the
K variables. Group assignment is reported as a vector M. Modulairty, Q,
measures the strength of the division into groups.
}
\examples{
## data(SampleInputMatrix)
## maybe str(SampleInputMatrix) ; plot(SampleInputMatrix) ...
}
\keyword{datasets}
