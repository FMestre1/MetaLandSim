\name{remove.species}
\alias{remove.species}
\title{
Remove the species occupancy from the landscape
}
\description{
This function converts an object of class 'metapopulation' (with the species occupancy) in a object of class 'landscape' (without the species occupancy).
}
\usage{
remove.species(sp)
}
\arguments{
  \item{sp}{
Object of class 'metapopulation'.
}
}
\value{
Delivers an object of class 'landscape'.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}, \code{\link{species.graph}}
}
\examples{

data(occ.landscape)

rl1 <- remove.species(sp=occ.landscape)

}
