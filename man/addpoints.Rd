\name{addpoints}
\alias{addpoints}
\title{
Add a given number of patches to a landscape
}
\description{
Adds a given number of patches to the landscape.
}
\usage{addpoints(rl, nr)}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
  \item{nr}{
Number of patches to be added (see 'note').
}
}
\value{
Returns an object of class 'landscape'.
}
\author{
Frederico Mestre and Fernando Canovas
}
\note{
The number of patches to be added might be impaired by the minimum distance between points.
}
\seealso{
\code{\link{rland.graph}}, \code{\link{removepoints}} 
}
\examples{

data(rland)

#Checking the number of patches in the starting landscape:

rland$number.patches

#60

#Adding 10 patches to a landscape:

rl1 <- addpoints(rl=rland, nr=10)

#Checking the number of patches in the output landscape:

rl1$number.patches

#70

}
