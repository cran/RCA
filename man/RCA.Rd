\name{RCA}
\alias{RCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relational Class Analysis
%%  ~~function to do ... ~~
}
\description{
Relational Class Analysis (RCA) is a method for
detecting heterogeneity in attitudinal data. RCA-R is a
C-based package that implements the RCA algorithm and
generates a simple-to-use command line executable program.
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
This is a R function using the C implementation of Relational Class Analysis (RCA).

Let X be a dataset of size N x K. X must not include missing data.
RCA finds an optimal division of X into G groups, such that each group
of observations follows a distinctive pattern of raltionships between the
K variables. Group assignment is reported as a vector M. Modulairty, Q,
measures the strength of the division into groups.

RCA-R is a C-based package that implements the RCA algorithm and generates
a simple-to-use command line executable program.
You will need to follow the steps below in order to compile the RCA-R code
and generate a local copy of the executable.

}
\value{
Results include:

a. [your_returned_data]$member: Assignment vector (size 1xN)
b. [your_returned_data]$mod: Modularity (double)
c. [your_returned_data]$merge: Merge tree (array of varying size)
d. [your_returned_data]$stats: Stats (2 doubles)
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Goldberg, A., AJS 116(5): 1397-1436 (http://http://www.jstor.org/stable/10.1086/657976)}
\author{
Amir Goldberg, Jinjian Zhai, Gabor Csardi
}
\note{
Depended on igraph library.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
igraph
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(inputMatrix, bootstrap=1000, p_value=0.05) { 
    # If the library hasn't been loaded yet, load it 
    if (!is.loaded('RCA_R_wrapper')) { 
        dyn.load('RCA.so') 
    } 
    # Call the C function. A list of parameters values after the function is called 
    # is returned, assigned to the same names as are given before the 
    # = signs in the arguments. 
    # Change 2D matrix data to 1D data
	Dim=dim(inputMatrix) # Dim=[row,col]
	inputData=double(Dim[1]*Dim[2])
	for(i in 1:Dim[1]){
        for(j in 1:Dim[2]){
            inputData[j+(i-1)*Dim[2]]=inputMatrix[i,j]
        }
    }
    # Calling RCA-R function
    returned_data = .C('RCA_R_wrapper', inputData=inputData, numObs=as.integer(Dim[1]), numVars=as.integer(Dim[2]), bootstrap=as.integer(bootstrap), z_score=as.double(qnorm(1-p_value/2)), member=integer(Dim[1]), mergeDim=integer(2), mod=as.double(0), stats=double(2), merge=matrix(0,1,1) )
    # get membership tag one by one because it is too big an array for R to handle at the same time. 
    for(i in 0:(Dim[1]-1)){
        each=.C('get', i=as.integer(i), element=as.integer(1))
        returned_data$member[i+1]=each$element
    }
    
    # C only accepts 1D input. We have to transform merge result here.
    row=returned_data$mergeDim[1]
    col=returned_data$mergeDim[2]
    merges2D=matrix(0,row,col)
    mergeResult=.C('getMerges', merges=integer(row*col))
    for(i in 1:row){
        for(j in 1:col){
            merges2D[i,j]=mergeResult$merges[(j-1)*row+i]
        }
    }
    result=list(member=returned_data$member, mod=returned_data$mod, merge=merges2D)
    # Return the value of the result parameter 
    return(result) 
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RCA , statistics, R}
