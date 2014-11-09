\name{landscape}
\alias{landscape}

\title{
Class 'landscape'
}
\description{

Class representing a landscape graph, as produced by \code{\link{rland.graph}}, \code{\link{convert.graph}} and \code{\link{import.shape}}.
}

\section{Slots}{
\itemize{
  \item mapsize - Side of the landscape in meters.
  \item minimum.distance - Minimum distance between patches centroids, in meters.
  \item mean.area - Mean patch area in hectares.
  \item SD.area - Standard deviation of patches area.
  \item number.patches - Total number of patches.
  \item dispersal - Species mean dispersal ability, in meters.
  \item nodes.characteristics - Data frame with patch (node) information (coordinates, area, radius, cluster, distance to nearest neighbor and ID). 
}
}

\author{
Frederico Mestre and Fernando Canovas
}
