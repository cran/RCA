\name{RCA-package}
\alias{RCA}
\alias{print.RCA}
\alias{plot.RCA}
\docType{package}

\title{
Relational Class Analysis
}

\description{
Relational Class Analysis (RCA) is a method for detecting heterogeneity in attitudinal data (as described in Goldberg A., 2011, Am. J. Soc, 116(5)).

The \code{RCA} function produces an object of class \code{RCA}. This object includes: (1) a vector of class membership ------ \code{$membership}; (2) a list of modules ------ \code{$modules}; and (3) the relationality matrix ------ \code{$R}.
}

\usage{
RCA(matrix, max = NULL, min = NULL, num = 1000, 
  alpha = 0.05)

\method{print}{RCA}(x, ...)

\method{plot}{RCA}(x, module = NULL, colorblind = FALSE, 
  heatmap = TRUE, heat_labels = FALSE, 
  drop_neg_ties = TRUE, layout = layout.kamada.kawai, 
  edge_color = "gray", vertex_color = "white", 
  vertex_frame_color = "black", vertex_size = 20, 
  vertex_label_color = "black", vertex_label_cex = 0.8, 
  margin = 0, ...)
}

\arguments{
\item{matrix}{ a matrix of size \code{m} x \code{n}. Rows correspond to observations and columns correspond to variables.}

\item{max}{ either a single integer or a vector of length \code{n}. If an integer is specified, its value is used as the maximum for all variables in the matrix. If a vector is specified, then each element is used as the maximum value for each variable. The default value is \code{NULL}, in which case \code{max} is assumed to be the maximum of each variable in the input matrix.}

\item{min}{ either a single integer or a vector of length \code{n}. If an integer is specified, its value is used as the minimum for all variables in the matrix. If a vector is specified, then each element is used as the minimum value for each variable. The default value is \code{NULL}, in which case \code{min} is assumed to be the minimum of each variable in the input matrix.}

\item{num}{ the number of bootstrap samples to be used in testing for significant relationalities. The default is \code{1000}.}

\item{alpha}{ a value between \code{0} and \code{1} specifying the significance level to be used in testing for significant relationalities. The default is \code{0.05}.}

\item{x}{ an object of class \code{RCA}.}

\item{module}{ the module number to graph.}

\item{colorblind}{ change the graph color schemes to be colorblind friendly. Default is \code{FALSE}.}

\item{heatmap}{ choose which type of graph to produce. Default is \code{TRUE} for heatmap. \code{FALSE} produces an igraph.}

\item{heat_labels}{ add correlation labels to the heatmap. Default is \code{FALSE}.}

\item{drop_neg_ties}{ drop negative graph edges. Default is \code{TRUE}.}

\item{layout}{ choose a layout for plotting graphs. Default is \code{layout.kamada.kawai}. Additional layouts can be chosen from \code{layout\{igraph\}}.}

\item{edge_color}{ color of the graph edges. Default is \code{gray}.}

\item{vertex_color}{ color of the vertices. Default is \code{white}.}

\item{vertex_frame_color}{ color of the vertex frames. Default is \code{black}.}

\item{vertex_size}{ vertex frame size. Default is \code{20}.}

\item{vertex_label_color}{ color of the vertex labels. Default is \code{black}.}

\item{vertex_label_cex}{ size of the vertex labels. Default is \code{0.8}.}

\item{margin}{ adjust the margins of the graph. Default is \code{0}. Negative values reduce the margins.}

\item{...}{ additional arguments}
}

\details{
RCA computes the relationality between all pairs of observations and partitions the sample into subgroups of schematically overlapping respondents such that those who subscribe to the same logic are clustered together. It accomplishes this by performing the following sequence: 

1.) RCA computes a measure of relationality for each pair of observations in the data set, resulting in a complete, undirected, and weighted graph. 

2.) Graph edges whose weights are statistically insignificant are removed. Remaining edges are transformed by their absolute value. 

3.) RCA partitions the graph into subgroups of schematically similar observations using a graph-partitioning algorithm. 
}

\author{
Amir Goldberg, Sarah K. Stein

Maintainer: Amir Goldberg <amirgo@stanford.edu>
}

\references{
Goldberg, Amir. "Mapping shared understandings using relational class analysis: the case of the cultural omnivore reexamined." American Journal of Sociology 116.5 (2011): 1397-1436.
}

\keyword{ package }

\examples{
# Test matrix
set.seed <- 2
matrix <- round(matrix(runif(100, min = 0, max = 6), 
  ncol = 5))

rca <- RCA(matrix)
rca <- RCA(matrix, max = rep(6, 5), min = rep(0, 5))
rca <- RCA(matrix, max = 6, min = 0)
rca <- RCA(matrix, num = 10000, alpha = 0.1)

summary(rca)
print(rca)

rca$membership
rca$modules[[1]]$matrix
rca$R

# Heatmap
plot(rca, module = 1)
plot(rca, module = 1, colorblind = TRUE, 
  heat_labels = TRUE)

# IGraph
plot(rca, module = 1, heatmap = FALSE)
plot(rca, module = 1, heatmap = FALSE, 
  layout = layout.circle, edge_color = "red", 
  vertex_color = "gray", vertex_frame_color = "red", 
  vertex_size = 30, vertex_label_color= "red", 
  vertex_label_cex = 1, margin = 0.2)
}
