\name{manage_expansion_sim}
\alias{manage_expansion_sim}
\title{
Simulate range expansion simulation
}
\description{
This function produces dispersal scenarios, considering different habitat networks properties, evaluating the variation in
dispersal speed and dispersal maximum distance (of range expansion).
}
\usage{
manage_expansion_sim(mapsize, dist_m, areaM, areaSD, Npatch,percI, 
                    param, b=1, tsteps, iter, variable,var_min,var_max,by)
}
\arguments{
  \item{mapsize}{
Landscape mosaic side length, in meters. To be internally passed to \code{\link{rland.graph}}
}
  \item{dist_m}{
Minimum distance between patches (centroid).To be internally passed to \code{\link{rland.graph}}
}
  \item{areaM}{
Mean area (in hectares). To be internally passed to \code{\link{rland.graph}}
}
  \item{areaSD}{
SD of the area of patches, in order to give variability to the patches area. To be internally passed to \code{\link{rland.graph}}
}
  \item{Npatch}{
Number of patches. To be internally passed to \code{\link{rland.graph}}
}
  \item{percI}{
Percentage of patch occupancy in the starting landscape. To be internally passed to  \code{\link{range_expansion}}
}
  \item{param}{
Parameter data frame delivered by \code{\link{parameter.estimate}}. To be internally passed to  \code{\link{range_expansion}}  It includes:
\itemize{
 \item alpha - Parameter relating extinction with distance.
 \item y - Parameter y in the colonization probability.  
 \item e - Parameter defining the extinction probability in a patch of unit area.
 \item x - Parameter scaling extinction risk with patch area.
 }
}  
  \item{b}{
Parameter scaling emigration with patch area (if conn='op1' or 'op2') in \code{\link{spom}}. By default, equal to 1.To be internally passed to  \code{\link{range_expansion}}
}
  \item{tsteps}{
Number of time steps to simulate (e.g. years).
}
  \item{iter}{
Number of iterations of the simulation procedure.
}
  \item{variable}{
Landscape graph characteristic to be varied. One of the following:
\itemize{
 \item area - Mean patch area (in hectares).
 \item dist - Minimum distance between patches (centroid).  
 \item npatch - Number of patches.
 \item sizevar - SD of the area of patches.
 }
}
  \item{var_min}{
Minimum value the changing variable can assume (beware of units used in each case).
}
  \item{var_max}{
Maximum value the changing variable can assume (beware of units used in each case).
}
  \item{by}{
Rate of variation of the changing variable.  
}
}
\details{
For details regarding the arguments that are to be internally passed to other functions, see the respective functions. Any of the arguments dist_m, areaM, areaSD, Npatch would be unnecessary if the respective variable is the one to be evaluated (it depends on the parameter variable). 
}
\value{
Returns a list of eight data frames. For the first four data frames (NORTH, SOUTH, EAST and WEST) each data frame's first column is the name of the variable to be changed. The other two columns are:

\item{MEAN EXPANSION SPEED }{Expansion speed in each simulated scenario. Speed given in km/time step}
\item{MAXIMUM EXPANSION DISTANCE }{Maximum distance of the expanded range, from an occupied site. Given in km.}\cr
The other four data frames have detailed information on the simulations for each of the values of parameter "variable". The first column has the distance (in km), and each of the following columns has the time step at which each distance was colonized for each of the simulations. 

}
\section{Warning}{
This function might be time consuming, and the code is experimental and should be improved in future versions of MetaLandSim.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}, \code{\link{range_expansion}}, \code{\link{expansion}}
}
\examples{

\dontrun{

data(param1)


sim_range <- manage_expansion_sim(mapsize=1000, dist_m=0, areaM, areaSD=0.001, 
                                  patch=300,percI=50, param=param1, b=1, 
                                  tsteps=100, iter=100,variable="area",var_min=0.01,
                                  var_max=0.6,by=0.1)
}
}
