\name{cluster.graph}
\alias{cluster.graph}
\title{
Delivers the number of patches per cluster
}
\description{
Returns a data frame with the number of nodes (habitat patches) in each component of the landscape graph (in this case a component is a group of patches connected by the species dispersal distance).
}
\usage{
cluster.graph(rl)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
}
\details{
The components are defined based on the species mean dispersal ability. This implies that the connectivity model between patches is binary (connected/not connected) as opposed to probabilistic.
}
\value{
This function returns a data frame with the number of patches of each component (group of patches).   
The returned data frame has two fields: cluster (Id of the component) and number of nodes (the number of nodes of the respective component).

}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}
}
\examples{

data(rland)

cluster.graph(rl=rland)

#Output:

#  cluster number of nodes
#1        1              11
#2        2               1
#3        3              13
#4        4               1
#5        5               1
#6        6              15
#7        7               2
#8        8               1
#9        9               3
#10      10               1
#11      11               1
#12      12               2
#13      13               4
#14      14               1
#15      15               1
#16      16               1
#17      17               1


}