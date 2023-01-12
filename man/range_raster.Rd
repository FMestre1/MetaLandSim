\name{range_raster}
\alias{range_raster}
\title{
Probability of occupancy, dispersal model 
}
\description{
This function creates the raster map with the expansion simulation, estimating probability of occupancy, at a given time step, based on species dispersal and landscape configuration. \code{range_raster} uses the output from \code{\link{range_expansion}} and a raster map with the species current occupancy.
}
\usage{
range_raster(presences.map, re.out, mask.map=NULL, plot=TRUE)
}
\arguments{
  \item{presences.map}{
string of the raster file name with species occurrence.
}
  \item{re.out}{
object of class list \code{expansion}. Output from \code{\link{range_expansion}}.
}
  \item{mask.map}{
default NULL. String of the raster file name with the mask. Usually, 1 over the area where the analyses should be done.
}
  \item{plot}{
default \code{TRUE}. Whether It will (\code{TRUE}) or will not (\code{FALSE}) return a graphics for the expansion model functions and raster maps with expansion probabilities in all four cardinal points.
}
}
\details{
The function automatically reads the raster input files (\code{presences.map} and \code{mask.map}, if present). Usually, 0 for absence and 1 for presence in every square cell over a given resolution. Note that the projection for the raster layer should be one of those supporting metric units (i.e., linear scale is equal in all directions around any point such as Transverse Mercator; see \url{https://spatialreference.org/}).

Then, it computes and fits a sigmoidal function for the expansion probability.

The user might have to manually adjust the starting values of the function \code{fit.sigmoid}, (defined internally in this function) if it has difficulty adjusting to the output of \code{\link{range_expansion}}.
}
\value{
Produces the spatial realization of the dispersal model, which is saved in the working directory (named 'PROB').
}
\references{
Mestre, F., Risk, B., Mira, A., Beja, P., Pita, R. (2017) <doi:10.1016/j.ecolmodel.2017.06.013>
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{range_expansion}}
}
\examples{

\dontrun{

#Loading required packages
library(MetaLandSim)

#Loading the range expansion simulation output and required rasters
data(rg_exp)

presences <- system.file("examples/presences.asc", package="MetaLandSim")
mask <- system.file("examples/landmask.asc", package="MetaLandSim")

range.map <- range_raster(presences.map=presences, re.out=rg_exp, mask.map=mask, plot = FALSE)

#Ploting the results with the terra package 
plot(range.map)

#Ploting the results with the rasterVis package 
require(rasterVis)
levelplot(range.map, contour=TRUE)

}
}
