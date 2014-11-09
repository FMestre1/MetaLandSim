\name{edge.graph}
\alias{edge.graph}
\title{
Produce an edge (links) data frame
}
\description{
Returns a data frame with the information on the connections between patches (assuming binary connections). 
}
\usage{
edge.graph(rl)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
}
\value{
Produces a data frame with the information on the edges (links): the IDs of both patches, the area, the coordinates and the Euclidean distance.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}
}
\examples{

data(rland)

edge_df <- edge.graph(rl=rland)

}