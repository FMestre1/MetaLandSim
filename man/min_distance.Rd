\name{min_distance}
\alias{min_distance}
\title{
Computes topological distance
}
\description{
Function to compute topological distance between patches. Topological distance is defined as the minimum number of links between any two patches.
}
\usage{
min_distance(rl)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
}
\value{
Returns a matrix with the topological distance between the nodes.
}
\author{
Frederico Mestre and Fernando Canovas.
}
\seealso{
\code{\link{rland.graph}}
}
\examples{

data(rland)

min_distance(rl=rland)

}
