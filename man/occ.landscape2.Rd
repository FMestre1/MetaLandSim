\name{occ.landscape2}
\alias{occ.landscape2}
\docType{data}
\title{
Sample landscape with 10 simulated occupancy snapshots 
}
\description{
Sample species occupancy in a network during 10 time steps. Simulated data.
}
\usage{data(occ.landscape2)}
\format{
A list with the following elements:
\itemize{
 \item mapsize - landscape mosaic side length, in meters.
 \item minimum.distance - minimum distance between patches centroids.
 \item mean.area - mean area, in hectares.
 \item SD.area - standard deviation of the area.
 \item number.patches - number of patches.
 \item dispersal - mean dispersal ability of the species.
 \item distance.to.neighbours - data frame with pairwise distance between patches.
 \item nodes.characteristics - data frame with the characteristics of each patch, (species 1 to 10 - occupancy snapshots).
}
}
\examples{
data(occ.landscape2)
}
\keyword{datasets}
