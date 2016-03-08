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
\value{
This function returns a list, of class 'expansion', of four data frames with the proportion of occupations at several distances from the closest occupied landscape mosaic. These four data frames correspond to the proportion of occupation to the north, south, east and west.
Each data frame has the following columns:
  \itemize{
  \item DISTANCE - Distance (mapsize x number of landscapes).
  \item OCCUPATION - How many times did the landscape at this distance got occupied by the species (from a total of 'iter' repetitions). 
  \item PROPORTION - Proportion of occupation for the landscape at this distance (OCCUPATION/iter).
  }
}
\author{
Frederico Mestre and Fernando Canovas
}
\note{
Depending on computing power and number of iterations (parameter iter) this function can take from a few hours to several days to run.
}
\seealso{
\code{\link{range_raster}} 
}
\examples{

\dontrun{
#Produce a model of range expansion:
#Note: this function should be run with >100 iterations (parameter 'iter').

data(rland)
data(param1)

rg_exp1 <- range_expansion(rl=rland, percI=50, param=param1, b=1, tsteps=100, iter=100)
}
}
