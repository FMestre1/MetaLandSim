\name{iterate.graph}
\alias{iterate.graph}
\title{
Simulate landscape series occupation
}
\description{
Repeats the process of simulation by \code{\link{simulate_graph}} as many times as required (argument 'iter'). 
}
\usage{
iterate.graph(iter, mapsize, dist_m, areaM, areaSD, Npatch, disp, 
		span, par1 = "none", par2 = NULL, par3 = NULL, par4 = NULL, 
		par5 = NULL, method = "percentage", parm, nsew = "none", 
		succ="none", param_df, kern, conn, colnz, ext, beta1, 
		b = 1, c1 = NULL, c2 = NULL, z = NULL, R = NULL, graph)
}
\arguments{
  \item{iter}{
Number of repetitions of the simulation.
}
  \item{mapsize}{
Landscape mosaic side length, in meters. To be internally passed to \code{\link{rland.graph}}. 
}
  \item{dist_m}{
Minimum distance between patches (centroid). To be internally passed to \code{\link{rland.graph}}.
}
  \item{areaM}{
Mean area (in hectares). To be internally passed to \code{\link{rland.graph}}.
}
  \item{areaSD}{
SD of the area of patches, in order to give variability to the patches area. To be internally passed to \code{\link{rland.graph}}.
}
  \item{Npatch}{
Number of patches (might be impaired by the dist_m, see the "Note" section). To be internally passed to \code{\link{rland.graph}}.
}
  \item{disp}{
Species mean dispersal ability, in meters. To be internally passed to \code{\link{rland.graph}}.
}
  \item{span}{
Number of time steps (e.g. years) to simulate. To be internally passed to \code{\link{span.graph}}.
}
  \item{par1}{
One of the following (default 'none'):
\itemize{
 \item 'hab' percentage of the number of patches to eliminate.
 \item 'dincr' minimal distance (between centroids of patches) increase over the simulation (in meters).
 \item 'darea' percentage of increase/decrease of the mean area of patches, without changing SD.
 \item 'stoc' simultaneous creation and destruction of patches.
 \item 'ncsd' simultaneous creation and destruction of patches to the north and south of the landscape.
 \item 'aggr' correlated habitat destruction.
 \item 'none' no change. 
}
To be internally passed to \code{\link{span.graph}}.
}
  \item{par2}{
Parameter specifying details for the options in par1: percentage of patches do delete (if par1 = 'hab'); distance, in meters (if par1 = 'dincr'); percentage of increase/decrease of the mean area of patches (if par1 = 'area'); percentage of new patches (if par1 = 'stoc'); 'northerndness' of created patches (if par1 = 'ncsd'); percentage of destroyed patches (if par1 = 'aggr').
To be internally passed to \code{\link{span.graph}}. Default NULL.
}
  \item{par3}{
Additional parameter specifying details for the options in par1: percentage of destroyed patches (if par1 = 'stoc'); 'southerndness' of destroyed patches (if par1 = 'ncsd'); aggregation of destruction (if par1 = 'aggr'). Minimum area for patch deletion, in hectares (if par1='darea').
To be internally passed to \code{\link{span.graph}}. Default NULL.
}
  \item{par4}{
Percentage of created patches (if par1 = 'ncsd'). To be internally passed to \code{\link{span.graph}}. Default NULL.
}
  \item{par5}{
Percentage of destroyed patches (if par1 = 'ncsd'). To be internally passed to \code{\link{span.graph}}. Default NULL.
}
  \item{method}{
One of the following (default 'percentage'):
     \strong{click} - individually select the patches with occurrence of the species by clicking on the map. Use only for individual landscape simulations. However, this option should not be used with iterate.graph.
     \strong{percentage} - percentage of the patches to by occupied by the species.
     \strong{number} - number of patches to be occupied by the species.
To be internally passed to \code{\link{species.graph}}.}
  \item{parm}{
parameter to specify the species occurrence - either percentage of occupied patches or number of occupied patches, depending on the method chosen. To be internally passed to \code{\link{species.graph}}.
}
  \item{nsew}{
'N', 'S', 'E', 'W' or none - point of entry of the species in the landscape. By default set to "none". To be internally passed to \code{\link{species.graph}}.
}
  \item{succ}{
Set the preference of the species for patch successional stage: 'none', 'early', 'mid' and 'late'.
}
  \item{param_df}{
Parameter data frame delivered by \code{\link{parameter.estimate}}, including:
\itemize{
 \item alpha - Parameter relating extinction with distance.
 \item y - Parameter y in the colonization probability.  
 \item e - Parameter defining the extinction probability in a patch of unit area.
 \item x - Parameter scaling extinction risk with patch area.
 }
To be internally passed to \code{\link{simulate_graph}}.
}
  \item{kern}{
'op1' or 'op2'. Dispersal kernel. See details in the \code{\link{spom}} function. To be internally passed to \code{\link{spom}}.
}
  \item{conn}{
'op1' or 'op2'. Connectivity function. See details in the \code{\link{spom}} function. To be internally passed to \code{\link{spom}}.
}
  \item{colnz}{
'op1', 'op2' or 'op3'. Colonization function. See details in the \code{\link{spom}} function. To be internally passed to \code{\link{spom}}.
}
  \item{ext}{
'op1', 'op2' or 'op3'. Extinction function. See details in the \code{\link{spom}} function. To be internally passed to \code{\link{spom}}.
}
  \item{beta1}{
Parameter affecting long distance dispersal probability (if the Kern='op2'). To be internally passed to \code{\link{spom}}.
}
  \item{b}{
Parameter scaling emigration with patch area (if conn='op1' or 'op2'). To be internally passed to \code{\link{spom}}. By default set to 1.
}
  \item{c1}{
Parameter scaling immigration with the focal patch area (if conn='op2'). To be internally passed to \code{\link{spom}}.
}
  \item{c2}{
Parameter c in the option 3 of the colonization probability (if colnz='op3'). To be internally passed to \code{\link{spom}}.
}
  \item{z}{
Parameter giving the strength of the Allee effect (if colnz='op3'). To be internally passed to \code{\link{spom}}.
}
  \item{R}{
Parameter giving the strength of the Rescue effect (if ext='op3'). To be internally passed to \code{\link{spom}}.
}
  \item{graph}{
TRUE/FALSE, to show graphic output.
}
}
\value{
Returns a list of five data frames with information regarding the values of mean area, mean inter-patch distance, number of patches occupancy and patch occupancy turnover in each of the iterations, as well as the mean values and SD.     
}
\references{
References in the \code{\link{spom}} function.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}, \code{\link{span.graph}}, \code{\link{species.graph}}, \code{\link{simulate_graph}}, \code{\link{spom}}
}
\examples{
\dontrun{
data(param1)

#Example with 2 iterations (ideally >100):

it1 <- iterate.graph(iter = 2, mapsize =10000, dist_m = 10, areaM = 0.05, 
			areaSD = 0.02, Npatch = 250, disp = 800, span = 100, 
			par1 = "hab", par2 = 2, par3 = NULL, par4 = NULL, 
			par5 = NULL, method = "percentage", parm = 50, 
			nsew = "none", succ="none", param_df = param1,kern = "op1", 
			conn = "op1", colnz = "op1", ext = "op1", 
			beta1 = NULL, b = 1, c1 = NULL, c2 = NULL, z = NULL, 
			R = NULL, graph =TRUE)
}
}
