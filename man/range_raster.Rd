\name{range_raster}
\alias{range_raster}
\title{
Probability of occupancy, dispersal model 
}
\description{
This function creates the raster map with the expansion simulation, estimating probability of occupancy, at a given time step, based on species dispersal and landscape configuration. \code{range_raster} uses the output from \code{\link{range_expansion}} and a raster map with the species current occupancy.
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

Then, it computes and fits single sigmoidal functions for every direction on the expansion movements (building four sub-models, one to each main cardinal direction), as previously computed by \code{\link{range_expansion}}. Four different raster maps are generated (the sub-models), each estimating the probability of expansion for north, south, east and west directions. The four maps are finally summarized into a single range expansion map, which is returned to the user as an object of class \code{RasterLayer} and saved in the working directory. These four maps do not directly express the probabilities in the output of \code{\link{range_expansion}}. Rather, the outputs are weighted by directionality, so that e.g. the north model favours the dispersal towards the north while penalizing dispersal in every other direction. As such, the resulting dispersal model is not a direct spatial transcription of the four data frames provided by \code{\link{range_expansion}} but an interpretation weighted by the spatial context given by directionality.\cr

Additionally, a raster file is computed, showing the time steps at which each distance is reached. This output depicts the adjustment of a linear model to the output of the four sub-models together, considering the four equally. This output should not be used as guideline to mask the model if running to several time periods. For example if projections to 2050 and 2080 are to be made than the dispersal model should be run twice, adjusting the time steps to the desired date.\cr

This function internally uses a connection to \href{http://grass.osgeo.org/}{GRASS GIS software} through the package \code{rgrass7-package}, in order to increase the performance for geographical calculations.\cr

Finally, the user might have to manually adjust the starting values of the function \code{fit.sigmoid}, (defined internally in this function) if it has difficulty adjusting to the output of \code{\link{range_expansion}}.
}
\value{
Produces the spatial realization of the dispersal model, composed by a stack of two objects of the class \code{RasterLayer} (see \code{\link{Raster-class}} package for further description), with the probability of occupancy and the time step a given distance is occupied. This version of MetaLandSim uses GRASS, version 7 through the package \href{https://CRAN.R-project.org/package=rgrass7}{rgrass7}. Additionally these rasters are saved in the working directory (files 'PROB' and 'TSTEP') defined by the user and can be directly imported to any GIS software.
}
\references{
The same as range_expansion.
}
\author{
Frederico Mestre and Fernando Canovas
}
\note{
This function depends on \href{https://CRAN.R-project.org/package=rgrass7}{rgrass7}.
}
\seealso{
\code{\link{range_expansion}}, \code{\link{Raster-class}}, \href{https://CRAN.R-project.org/package=rgrass7}{rgrass7}, \code{\link{initGRASS}}
}
\examples{

\dontrun{

#Installing the rgrass7 development version from GitHub (v. 0.2-1)
#(for now only in GitHub)
library(devtools)
install_github("rsbivand/rgrass7")

#Loading required packages
library(MetaLandSim)
library(rgrass7)

#Loading the range expansion simulation output and required rasters
data(rg_exp)
presences <- system.file("examples/presences.asc", package="MetaLandSim")
mask <- system.file("examples/landmask.asc", package="MetaLandSim")

if (packageVersion("rgrass7") >= "0.2.1") use_sp()

#Initializing a GRASS session in a temporal directory
#(The user should insert the correct path to the executable GRASS file)

#### Under Linux systems:
initGRASS("/usr/bin/grass", home=tempdir())

#### Under Windows systems:
initGRASS("C:/Program Files/GRASS GIS 7.6", home=tempdir())

range.map <- range_raster(presences.map=presences, re.out=rg_exp, mask.map=mask)

#Ploting the results with the raster package 
plot(range.map)

#Ploting the results with the rasterVis package 
require(rasterVis)
levelplot(range.map, contour=TRUE)

}
}
