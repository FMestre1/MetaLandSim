\name{summary_metapopulation}
\alias{summary_metapopulation}
\title{
Summarize 'metapopulation' class objects
}
\description{
This function summarizes a \code{\link{metapopulation}} class object. 
}
\usage{
summary_metapopulation(object)
}
\arguments{
  \item{object}{
Object of class \code{\link{metapopulation}}
}
}
\details{
This function can be used to retrieve basic information on the objects of class 'metapopulation'.
}
\value{
Returns a data frame with the following information on a \code{\link{metapopulation}} class object:

\item{landscape area (hectares) }{Landscape mosaic area, in hectares}
\item{number of patches }{Number of patches in the landscape}
\item{mean patch area (hectares) }{Mean patch area, in hectares}
\item{SD patch area }{SD of the patch area}
\item{mean distance amongst patches (meters) }{Mean inter-patch distance, in meters}
\item{minimum distance amongst patches (meters) }{Minimum inter-patch distance, in meters}
\item{species occurrence - snapshot }{Occupation data of the focal species, numbered from 1 to the number of snapshots}
}
\author{
Frederico Mestre and Fernando Canovas
}
\note{
The minimum distance between patches is different from that given in the object of class 'landscape', in the slot 'minimum.distance'. This is because this output is computed from the landscape structure and the one in the 'landscape' object was the parameter used to built the landscape. The minimum inter-patch distance given as a parameter in the function \code{\link{rland.graph}} will consider distance between patch centroids. The minimum inter-patch distance returned here considers the edge-to-edge distance, so this might be smaller that the parameter of \code{\link{rland.graph}}. In order to see the difference between centroid-to-centroid and edge-to-edge inter-patch distance compute both using the \code{\link{matrix.graph}} function (methods are 'centr_distance' and 'euc_distance', respectively).
}
\seealso{
\code{\link{species.graph}}, \code{\link{metapopulation}}, \code{\link{matrix.graph}}
}
\examples{

data(occ.landscape)
data(occ.landscape2)

summary_metapopulation(object=occ.landscape)

#                                            Value
#landscape area (hectares)                 100.000
#number of patches                          60.000
#mean patch area (hectares)                  0.061
#SD patch area                               0.041
#mean distance amongst patches (meters)    528.345
#minimum distance amongst patches (meters)  51.780
#species occurrence - snapshot 1            50.000


summary_metapopulation(object=occ.landscape2)

#                                            Value
#landscape area (hectares)                 100.000
#number of patches                          60.000
#mean patch area (hectares)                  0.069
#SD patch area                               0.039
#mean distance amongst patches (meters)    521.717
#minimum distance amongst patches (meters)  45.905
#species occurrence - snapshot 1            50.000
#species occurrence - snapshot 2            58.333
#species occurrence - snapshot 3            61.667
#species occurrence - snapshot 4            61.667
#species occurrence - snapshot 5            58.333
#species occurrence - snapshot 6            60.000
#species occurrence - snapshot 7            70.000
#species occurrence - snapshot 8            68.333
#species occurrence - snapshot 9            68.333
#species occurrence - snapshot 10           56.667

}
