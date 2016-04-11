\name{species.graph}
\alias{species.graph}
\title{
Simulate landscape occupation
}
\description{
Given a set of parameters, this function allows to simulate the occupation of an empty landscape, class "metapopulation".
}
\usage{
species.graph(rl, method = 'percentage', parm, nsew = 'none', plotG = TRUE)
}
\arguments{
  \item{rl}{
     Object of class "landscape".
}
  \item{method}{
     One of the following (default 'percentage'):
     \strong{click} - individually select the patches with occurrence of the species by clicking on the map. Use only for individual landscape simulations.
     \strong{percentage} - percentage of the patches to be occupied by the species.
     \strong{number} - number of patches to be occupied by the species.
}
  \item{parm}{
     Parameter to specify the species occurrence - either percentage of occupied patches or number of occupied patches, depending on the method chosen.
}
  \item{nsew}{
     'N', 'S', 'E', 'W' or none - point of entry of the species in the landscape. By default set to "none".
}
  \item{plotG}{
     TRUE/FALSE, to show graphic output.
}
}
\value{
Returns a list, with the following elements:
\itemize{
  \item mapsize - Landscape mosaic side length, in meters.
  \item minimum.distance - Minimum distance between patches centroids, in meters.
  \item mean.area - Mean patch area in hectares.
  \item SD.area - Standard deviation of patches area.
  \item number.patches - Total number of patches.
  \item dispersal - Species mean dispersal ability, in meters.
  \item distance.to.neighbours - Data frame with pairwise distance between patches, in meters.
  \item nodes.characteristics - Data frame with patch (node) information (coordinates, area, radius, cluster, distance to nearest neighbour, ID and species). 
        }
An additional field, colour, has only graphical purposes.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}, \code{\link{simulate_graph}}, \code{\link{remove.species}}
}
\examples{

data(rland)

##Creating a 50\% occupation in an empty landscape (using the "landscape" dataset):

sp1 <- species.graph(rl=rland, method="percentage", parm=50, nsew="none", plotG=TRUE)
}


