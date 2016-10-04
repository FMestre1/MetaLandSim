\name{range_expansion}
\alias{range_expansion}
\title{
Produce a range expansion model
}
\description{
This function returns the expansion probability, from a landscape with a given set of parameters, into the four cardinal directions. This can subsequently be converted in a dispersal model by the function \code{\link{range_raster}}. The dispersal model can be combined with an ecological niche model.
}
\usage{
range_expansion(rl, percI, param, b, tsteps, iter)
}
\arguments{
  \item{rl}{
Object of class 'landscape'. Starting landscape for the expansion procedure.
}
  \item{percI}{
Pecentage of patch occupancy in the starting landscape.
}
  \item{param}{
Parameter data frame delivered by \code{\link{parameter.estimate}}, including:
\itemize{
 \item alpha - Parameter relating extinction with distance.
 \item y - Parameter y in the colonization probability.  
 \item e - Parameter defining the extinction probability in a patch of unit area.
 \item x - Parameter scaling extinction risk with patch area.
 }
}
  \item{b}{
Parameter scaling emigration with patch area (if conn='op1' or 'op2') in \code{\link{spom}}. By default, equal to 1.
}
  \item{tsteps}{
Number of time steps to simulate (e.g. years).
}
  \item{iter}{
Number of iterations of the simulation procedure.
}
}
\details{
The expansion algorithm has been improved, since the paper Mestre et al. (2016) describing the package was published. Now, instead of the transition between adjacent landscape units being dictated by the occupation of a spurious node (representing the margin through which the expansion takes place) a somewhat more realistic approach is followed. If, during the metapopulational dynamics simulation, any patch located between the landscape unit (LU) margin and a parallel line placed at a distance equivalent to half of the mean dispersal ability of the species is occupied, than the algorithm assumes that the species will have the ability to go across to the next LU. In this new empty LU initial occupation is defined as follows: a new line is placed, with a spacing equivalent to half the dispersal ability of the species. In the area defined by the margins of the LU and this line the species will occupy in the same proportion as in the preceding LU.\cr
An example, with the expansion eastward (the process is repeated 4 times, one in each cardinal direction): One LU with 200 patches, occupation of 50\%, mapsize (length of the LU side) 1000 species mean dispersal ability of 200. Metapopulational dynamics are simulated until one patch is occupied in a area defined by the north, south, east margins of the LU and a vertical line placed at x=800 (1000-200). Then, if any of this patches is occupied, a new LU (a random realization of the same parameter set) is created and initial occupation is defined at an area defined between the west, north and south margins and a line placed at x=200. This occupation level has the same percentage as the previous landscape.
}
\value{
This function returns a list, of class 'expansion', of four data frames with the proportion of occupations at several distances from the closest occupied landscape mosaic. These four data frames correspond to the proportion of occupation to the north, south, east and west.
Each data frame has the following columns:
  \itemize{
  \item DISTANCE - Distance (mapsize x number of landscapes).
  \item OCCUPATION - How many times did the landscape at this distance got occupied by the species (from a total of 'iter' repetitions). 
  \item PROPORTION - Proportion of occupation for the landscape at this distance (OCCUPATION/iter).
  \item TIME STEP - The average time steps at which a given distance is occupied.
  }
}
\author{
Frederico Mestre and Fernando Canovas
}
\note{
Depending on computing power and number of iterations (parameter 'iter') this function can take from a few hours to several days to run.
}
\seealso{
\code{\link{range_raster}} 
}
\examples{

\dontrun{
#Produce a model of range expansion:
#Note: this function should be run with >100 iterations (parameter "iter").

data(rland)
data(param2)

rg_exp1 <- range_expansion(rl=rland, percI=50, param=param2, b=1, tsteps=100, iter=100)
}
}
