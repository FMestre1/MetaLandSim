\name{cluster.id}
\alias{cluster.id}
\title{
Classify patches in clusters
}
\description{
reclassify clusters of a landscape according to a given mean dispersal distance.
}
\usage{
cluster.id(rl)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
}
\details{
After changing the landscape some components (groups of connected patches) might suffer changes (e.g. the removal of patches might split components). This function re-attributes a code to each patch, identifying the groups of connected patches (components), after this type of disturbance to the habitat network.Mainly to be used internally. 
}
\value{
Returns the same landscape object, with the clusters reclassified.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}
}
\examples{

data(rland)

#After removing 30 (50\%) of the patches of a landscape:

rland2 <- removepoints(rl=rland, nr=35)

#A reclassification might be needed to identify components: 

rland2 <- cluster.id(rl=rland2)

#After removing 35 patches, there's a different number of components:

components.graph(rl=rland) 

#21

components.graph(rl=rland2) 

#16

}