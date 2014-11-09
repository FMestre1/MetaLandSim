\name{occ.landscape}
\alias{occ.landscape}
\docType{data}
\title{
Sample landscape with one simulated occupancy snapshot 
}
\description{
Sample random landscape graph, with species occupancy data (occupancy rate - 50\%). Simulated data.
}
\usage{data(occ.landscape)}
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
 \item nodes.characteristics - data frame with the characteristics of each patch.
}
}
\examples{

data(occ.landscape)

}
\keyword{datasets}
