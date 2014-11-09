\name{extract.graph}
\alias{extract.graph}
\title{
Extract landscape from span.graph generated list
}
\description{
Extracts a landscape from an object delivered by \code{\link{span.graph}}. The output is an object of class 'landscape'.
}
\usage{
extract.graph(rl, rlist, nr)
}
\arguments{
  \item{rl}{
Object of class 'landscape' used to generate the list, with \code{\link{span.graph}}.
}
  \item{rlist}{
Object delivered by \code{\link{span.graph}}.
}
  \item{nr}{
Position of the landscape in the list (rlist).
}
}
\value{
Delivers an object of class 'landscape'.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{span.graph}}, \code{\link{rland.graph}} 
}
\examples{

data(rland)
data(landscape_change)

#Extracting the landscape of the 50th time step:

rl1 <- extract.graph(rl=rland, rlist=landscape_change, nr=50)

}
