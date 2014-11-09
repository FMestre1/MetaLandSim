\name{cabrera}
\alias{cabrera}
\docType{data}
\title{
Modified patch occupancy data of Cabrera vole
}
\description{
One season patch occupancy dataset for \emph{Microtus cabrerae} in SW 
Portugal. This dataset is in the format produced by \link{species.graph}, 
\link{convert.graph} or \link{import.shape} (class 'metapopulation'), and it was created by converting a 
data frame using the function convert.graph. The data frame had the 
information of one snapshot of patch occupancy data of Cabrera vole 
(Microtus cabrera) in southwestern Portugal.
}
\usage{data(cabrera)}
\format{A list with the following elements:
\itemize{
 \item mapsize - 8200 (landscape mosaic side length, in meters).
 \item minimum.distance - 10.04 (minimum distance between patches centroids).
 \item mean.area - 0.46 (mean area, in hectares).
 \item SD.area - 1.05 (SD of the area).
 \item number.patches - 793 (number of patches).
 \item dispersal - 800 (mean dispersal ability of the species).
 \item distance.to.neighbours - data frame with pairwise distance between patches.
 \item nodes.characteristics - data frame with the characteristics of each patch.
}
}
\details{
To create this sample dataset the occupancy status of patches was 
scrambled, however the proportion of occupied patches was kept.
}
\source{
Original field data was obtained during project PERSIST (PTDC/BIA-BEC/105110/2008).
}
\examples{
data(cabrera)
}
\keyword{datasets}
