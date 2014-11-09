\name{summary_landscape}
\alias{summary_landscape}
\title{
Summarize 'landscape' class objects
}
\description{
This function summarizes a \code{\link{landscape}} class object. 
}
\usage{
summary_landscape(object)
}
\arguments{
  \item{object}{
Object of class \code{\link{landscape}}
}
}
\details{
This function can be used to retrieve basic information on the objects of class 'landscape'.
}
\value{
Returns a data frame with the following information on a \code{\link{landscape}} class object:

\item{landscape area (hectares) }{Landscape mosaic area, in hectares}
\item{number of patches }{Number of patches in the landscape}
\item{mean patch area (hectares) }{Mean patch area, in hectares}
\item{SD patch area }{SD of the patch area}
\item{mean distance amongst patches (meters) }{Mean inter-patch distance, in meters}
\item{minimum distance amongst patches (meters) }{Minimum inter-patch distance, in meters}
}
\author{
Frederico Mestre and Fernando Canovas
}
\note{
The minimum distance between patches is different from that given in the object of class 'landscape', in the slot 'minimum.distance'. This is because this output is computed from the landscape structure and the one in the 'landscape' object was the parameter used to built the landscape. The minimum inter-patch distance given as a parameter in the function \code{\link{rland.graph}} will consider distance between patch centroids. The minimum inter-patch distance returned here considers the edge-to-edge distance, so this might be smaller that the parameter of \code{\link{rland.graph}}. In order to see the difference between centroid-to-centroid and edge-to-edge inter-patch distance compute both using the \code{\link{matrix.graph}} function (methods are 'centr_distance' and 'euc_distance', respectively).
}
\seealso{
\code{\link{rland.graph}}, \code{\link{landscape}}, \code{\link{matrix.graph}}
}
\examples{

data(rland)

summary_landscape(object=rland)

#                                            Value
#landscape area (hectares)                 100.000
#number of patches                          60.000
#mean patch area (hectares)                  0.061
#SD patch area                               0.041
#mean distance amongst patches (meters)    528.345
#minimum distance amongst patches (meters)  51.780

}
