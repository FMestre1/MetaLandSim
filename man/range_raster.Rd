\name{range_raster}
\alias{range_raster}
\title{
Probability of occupancy, dispersal model 
}
\description{
This function intends to create a raster map, estimating probability of occupancy, at a given time step, based on species dispersal and landscape configuration. \code{range_raster} uses the output from \code{\link{range_expansion}} and a raster map with the species current occupancy.
}
\usage{
range_raster(presences.map, re.out, mask.map=NULL, plot.directions=TRUE)
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
  \item{plot.directions}{
default \code{TRUE}. Whether It will (\code{TRUE}) or will not (\code{FALSE}) return a graphics for the expansion model functions and raster maps with expansion probabilities in all four cardinal points.
}
}
\details{
The function automatically reads the raster input files (\code{presences.map} and \code{mask.map}, if present). Usually, 0 for absence and 1 for presence in every square cell over a given resolution. Supported file types are those that can be read via rgdal (see \href{http://www.gdal.org/formats_list.html}{gdal}. Note that the projection for the raster layer should be one of those supporting metric units (i.e., linear scale is equal in all directions around any point such as Transverse Mercator; see \url{http://spatialreference.org/}).

Then, it computes and fits single sigmoidal functions for every direction on the expansion movements, as previously calculated in \code{\link{range_expansion}} and provided here. Four different raster maps are generated, every one of them estimating the probability of expansion for north, south, east and west directions. The four maps are finally summarized into a single range expansion map, which is returned to the user as an object of class \code{RasterLayer}.\cr

Additionally, a raster file is computed, showing the time steps at which each distance is attained. The user should combine the dispersal model (the probability raster) with the time step raster in order to produce a map to any given time period.\cr

This function internally uses a connection to \href{http://grass.osgeo.org/}{GRASS GIS software} through the package \code{rgrass7-package}, in order to increase the performance for geographical calculations.\cr

Finally, the user might have to manually adjust the starting values of the function \code{fit.sigmoid}, defined internally in this function if it has difficulty to adjust to the output of \code{\link{range_expansion}}.
}
\value{
Produces the spatial realization of the dispersal model, composed by a stack of two objects of the class \code{RasterLayer} (see \code{\link{Raster-class}} package for further description), with the probability of occupancy and the time step a given distance is occupied considering only dispersal ability and landscape configuration. This version of MetaLandSim uses GRASS, version 7 through the package \href{http://cran.r-project.org/package=rgrass7}{rgrass7}. Additionally these rasters are saved in the working directory (files 'PROB' and 'TSTEP') defined by the user and can be directly imported to any GIS software.
}
\references{
The same as range_expansion.
}
\author{
Frederico Mestre and Fernando Canovas
}
\note{
This function depends on \href{http://cran.r-project.org/package=rgrass7}{rgrass7}.
}
\seealso{
\code{\link{range_expansion}}, \code{\link{Raster-class}}, \href{http://cran.r-project.org/package=rgrass7}{rgrass7}, \code{\link{initGRASS}}
}
\examples{

\dontrun{

data(rg_exp)

presences <- system.file("examples/presences.asc", package="MetaLandSim")
mask <- system.file("examples/landmask.asc", package="MetaLandSim")

require(rgrass7)

#Initializing a GRASS session in a temporal directory:
#### Under Linux systems:
initGRASS("/usr/bin/grass", home=tempdir())
#### Under Windows systems:
initGRASS("C:/GRASS", home=tempdir())

range_raster(presences.map=presences, re.out=rg_exp, mask.map=mask) -> range.map

plot(range.map)

require(rasterVis)

levelplot(range.map, contour=TRUE)
}
}
