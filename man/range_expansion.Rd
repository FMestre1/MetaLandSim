\name{range_expansion}
\alias{range_expansion}
\title{
Computes a range expansion model
}
\description{
This function returns the expansion probability, from a landscape with a given set of parameters. This can subsequently be converted in a dispersal model by the function \code{\link{range_raster}}.
}
\usage{
range_expansion(rl, percI, param, b, tsteps, iter, plot)
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
  \item{plot}{
Plot results.
}
}
\details{
The expansion algorithm has been improved, since the paper Mestre et al. (2017) describing the package was published. Now, instead of the transition between adjacent landscape units being dictated by the occupation of a spurious node (representing the margin through which the expansion takes place) a somewhat more realistic approach is followed. If, during the metapopulational dynamics simulation, any patch located between the landscape unit (LU) margin and a parallel line placed at a distance equivalent to half of the mean dispersal ability of the species is occupied, than the algorithm assumes that the species will have the ability to go across to the next LU. In this new empty LU initial occupation is defined as follows: a new line is placed, with a spacing equivalent to half the dispersal ability of the species. In the area defined by the margins of the LU and this line the species will occupy in the same proportion as in the preceding LU.\cr
After version 2.0.0,  the output, rather than considering distinct dispersal probabilities in all four cardinal directions (as in previous versions), considers the same probability of dispersal from a current presence in all directions. This does not change the results
in any meaningfull way given that these kinds of simulations require many iterations in which the distinctions between the dispersal to all four directions was diluted.
}
\value{
This function returns a data frame with the proportion of occupations at several distances from the closest occupied landscape mosaic. After version 2.0.0 the package uses the same dispersal probability in all directions relative to current presences.
the data frame has the following columns:
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
\references{
Mestre, F., Risk, B., Mira, A., Beja, P., Pita, R. (2017) <doi:10.1016/j.ecolmodel.2017.06.013>
}
\note{
Depending on computing power and number of iterations (parameter 'iter') this function can take some time to run.
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

rg_exp1 <- range_expansion(rl=rland, percI=80, param=param2, b=1, tsteps=100, iter=100)

}
}