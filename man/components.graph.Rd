\name{components.graph}
\alias{components.graph}
\title{
Number of components of a landscape
}
\description{
Returns the number of components in the landscape graph (in this case a component is a group of patches connected by the species dispersal distance).
}
\usage{
components.graph(rl)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
}
\value{
Returns the number of components (groups of connected patches) of a landscape.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}
}
\examples{

data(rland)

components.graph(rl=rland)

#21

}
