\name{plotL.graph}
\alias{plotL.graph}
\title{
Plot one landscape of the list created by span.graph
}
\description{
Plots a given landscape of a landscape sequence from \code{\link{span.graph}}.
}
\usage{
plotL.graph(rl, rlist, nr, species, links, ...)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
  \item{rlist}{
List returned by \code{\link{span.graph}}.
}
  \item{nr}{
index of the landscape to display graphically.
}
  \item{species}{
TRUE/FALSE, TRUE if 'rl' is of class 'metapopulation' or 'FALSE' if rl is of class 'landscape'.
}
  \item{links}{
TRUE/FALSE, show links between connected patches.
}
  \item{\dots}{
Other arguments.
}
}
\value{
Graphical display of the landscape.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{plot_graph}}, \code{\link{span.graph}}, \code{\link{rland.graph}}
}
\examples{

data(rland)
data(landscape_change)

plotL.graph(rl=rland, rlist=landscape_change, nr=50, species=FALSE, links=FALSE)

}
