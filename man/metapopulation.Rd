\name{metapopulation}
\alias{metapopulation}

\title{
Class 'metapopulation'
}
\description{

Class representing a landscape graph with species' patch occupancy data, as produced by \code{\link{species.graph}}, \code{\link{convert.graph}} and \code{\link{import.shape}}.
}

\section{Slots}{
\itemize{
  \item mapsize - Landscape mosaic side length, in meters.
  \item minimum.distance - Minimum distance between patches centroids, in meters.
  \item mean.area - Mean patch area in hectares.
  \item SD.area - Standard deviation of patches area.
  \item number.patches - Total number of patches.
  \item dispersal - Species mean dispersal ability, in meters.
  \item distance.to.neighbours - Data frame with pairwise distance between patches, in meters.
  \item nodes.characteristics - Data frame with patch (node) information (coordinates, area, radius, cluster, distance to nearest neighbor, ID and species). 
} 
}

\author{
Frederico Mestre and Fernando Canovas
}
