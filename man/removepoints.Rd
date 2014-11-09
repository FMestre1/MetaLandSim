\name{removepoints}
\alias{removepoints}
\title{
Remove a given number of patches from the landscape
}
\description{
Randomly removes a given number of patches from the landscape.
}
\usage{
removepoints(rl, nr)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
  \item{nr}{
Number of patches to remove.
}
}
\value{
Returns an object of class 'landscape'.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}, \code{\link{addpoints}} 
}
\examples{

data(rland)

#Checking the number of patches in the starting landscape:

rland$number.patches

#60

#Removing 10 patches from the landscape:

rl1 <- removepoints(rl=rland, nr=10)

#Checking the number of patches in the output landscape:

rl1$number.patches

#50

}