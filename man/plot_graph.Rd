\name{plot_graph}
\alias{plot_graph}
\title{
Graphical display of the landscape
}
\description{
Plots the landscape graph, with or without the species occupation (respectively lists returned by \code{\link{species.graph}} or \code{\link{rland.graph}}) and with or without the links between patches.
}
\usage{
plot_graph(rl, species, links)
}
\arguments{
  \item{rl}{
Object of class 'landscape' (species=FALSE) or 'metapopulation' (species=TRUE).
}
  \item{species}{
TRUE/FALSE, TRUE if 'x' is of class 'metapopulation' or 'FALSE' if x is of class 'landscape'.
}
  \item{links}{
TRUE/FALSE, show links between connected patches.
}
}
\value{
Graphical display of the landscape.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}, \code{\link{species.graph}}
}
\examples{

data(rland)
data(occ.landscape)

#Without the species occupancy:
plot_graph(rl=rland, species=FALSE, links=FALSE)

#With the species occupancy:
plot_graph(rl=occ.landscape, species=TRUE, links=FALSE)

}
