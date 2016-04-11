\name{matrix.graph}
\alias{matrix.graph}
\title{
Returning a matrix with information on connections between patches
}
\description{
Based on a landscape graph, this function allows the creation of a matrix of Euclidean distances (straight-line pairwise distance between the margins of all the patches), matrix of topological distances (minimum number of connections between any two patches) and adjacency matrix (this a matrix of 0 and 1, showing the adjacency between any two patches, where 0 means that the patches are not connected and 1 means that the patches are connected).
}
\usage{
matrix.graph(rl, mat)
}
\arguments{
  \item{rl}{
Object of class 'landscape'. 
}
  \item{mat}{
mat - one of the following:
\itemize{
 \item 'euc_distance' - euclidian distance between patches (edge-to-edge).
 \item 'centr_distance' - euclidian distance between patches (centroid-to-centroid). 
 \item 'adjacency' - adjacency matrix, with values d_ij, taking value 0 if patches i and j are not connected and value 1 if those patches are connected.
 \item 'top_matrix' - topological distance, with values d_ij, where the value d is the minimum number of connections between the patches i and j. Topological distance is defined as the minimum number of links between patches i and j.
 }
}
}
\value{
This function returns a matrix (each one of the specified matrices: Euclidean distance, topological distance and adjacency matrix).
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}
}
\examples{

data(rland)

#Computing matrix of topological distances:

matrix.graph(rl=rland, mat="top_matrix")

}
