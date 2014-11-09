\name{span.graph}
\alias{span.graph}
\title{
Simulate landscape dynamics over a number of time steps
}
\description{
This function gets an initial landscape graph and gradually applies changes. For a good review and classification of such changes see Bogaert et al. (2004) (not all described changes have been applied here). Future versions of the package should include other methods to change the landscape.
}
\usage{
span.graph(rl, span = 100, par1 = 'none', par2 = NULL, 
           par3 = NULL, par4 = NULL, par5 = NULL)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
  \item{span}{
Number of time steps (e.g. years) to simulate.
}
  \item{par1}{
One of the following (default 'none'):
\itemize{
 \item 'hab' percentage of the number of patches to eliminate.
 \item 'dincr' minimal distance (between centroids of patches) increase over the simulation (in meters).
 \item 'darea' percentage of increase/decrease of the mean area of patches, without changing SD. Patches with area <1 square meter are deleted.
 \item 'stoc' simultaneous creation and destruction of patches.
 \item 'ncsd' simultaneous creation and destruction of patches to the north and south of the landscape.
 \item 'aggr' correlated habitat destruction.
 \item 'none' no change. 
}
}
  \item{par2}{
Parameter specifying details for the options in par1: percentage of patches do delete (if par1='hab'); distance, in meters (if par1='dincr'); percentage of increase/decrease (increase with negative sign) of the mean area of patches (if par1='darea'); percentage of new patches (if par1='stoc'); 'northerndness' of created patches (if par1='ncsd'); percentage of destroyed patches (if par1='aggr'). 
}
  \item{par3}{
Additional parameter specifying details for the options in par1: percentage of destroyed patches (if par1='stoc'); 'southerndness' of destroyed patches (if par1='ncsd'); aggregation of destruction (if par1='aggr'). Minimum area for patch deletion, in hectares (if par1='darea').
}
  \item{par4}{
Percentage of created patches (if par1='ncsd').
}
  \item{par5}{
Percentage of destroyed patches (if par1='ncsd'). 
}
}
\value{
Returns a list of data frames with the nodes characteristics of a given number of landscapes that suffer a specified change. The fields of these data frames are the same as those from the nodes characteristics resulting from \code{\link{rland.graph}}.
}
\references{
Bogaert, J., Ceulemans, R., & Salvador-Van Eysenrode, D. (2004). Decision tree algorithm for detection of spatial processes in landscape transformation. Environmental Management, 33(1): 62-73.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}, \code{\link{simulate_graph}}, \code{\link{iterate.graph}}
}
\examples{

data(rland)

#Simulating a decrease of 5\% in the patches area through 100 time steps:

span1 <- span.graph(rl=rland, span=100, par1="darea", par2=5, par3=NULL, par4=NULL, par5=NULL)

}
