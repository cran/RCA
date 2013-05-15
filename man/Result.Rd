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

Each cell in the vector includes the group number the corresponding observation was assigned to. Group numbers range from 0 to G-1.

(b)[your_returned_data]$mod: Modularity (double)

A number ranging from 0 to 1 indicating the strength of the partition into groups.

(c)[your_returned_data]$merge:

Cluster tree, reported as a merge matrix (array of varying size). See igraph/leading.eigenvector.community for more information.

(d)[your_returned_data]$stats: Stats (2 doubles)

%Amir Goldberg April 12, 2013 10:46 PM: Please add a new line after each of the numbered items. Change the text as follows:
%a. [your_returned_data]$member: assignment vector (size 1xN). 
%b. [your_returned_data]$mod: Modularity (double). 
%c. [your_returned_data]$merge: 

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\keyword{datasets}
