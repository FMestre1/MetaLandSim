\name{simulate_graph}
\alias{simulate_graph}
\title{
Simulate species occupancy in one dynamic landscape
}
\description{
Simulates the species' occupation on a landscape sequence, resorting to the \code{\link{spom}} function.
}
\usage{
simulate_graph(rl, rlist, simulate.start, method, parm, nsew="none", succ="none", 
param_df, kern, conn, colnz, ext, beta1, b, c1, c2, z, R)
}
\arguments{
  \item{rl}{
Object of class 'landscape' or 'metapopulation'.
}
  \item{rlist}{
List delivered by \code{\link{span.graph}}.
}
  \item{simulate.start}{
TRUE (rl is of class 'landscape') or FALSE (rl is of class 'metapopulation')
}
  \item{method}{
One of the following:
     \strong{click} - individually select the patches with occurrence of the species by clicking on the map. Use only for individual landscape simulations. However, this option should not be used with iterate.graph.
     \strong{percentage} - percentage of the patches to by occupied by the species.
     \strong{number} - number of patches to be occupied by the species.
To be internally passed to \code{\link{species.graph}}.}
  \item{parm}{
Parameter to specify the species occurrence - either percentage of occupied patches or number of occupied patches, depending on the method chosen. To be internally passed to \code{\link{species.graph}}.
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
Parameter afecting long distance dispersal probability (if the Kern='op2'). To be internally passed to \code{\link{spom}}.
}
  \item{b}{
Parameter scaling emigration with patch area (if conn='op1' or 'op2'). To be internally passed to \code{\link{spom}}.
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
}
\value{
Returns a list of occupied landscapes, representing the same occupied landscape at different time steps.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{spom}}, \code{\link{span.graph}}, \code{\link{rland.graph}}, \code{\link{iterate.graph}}
}
\examples{

data(rland)
data(landscape_change)
data(param1)

sim1 <- simulate_graph(rl=rland, 
			rlist=landscape_change, 
			simulate.start=TRUE, 
			method="percentage", 
			parm=50, 
			nsew="none",
			succ = "none", 
			param_df=param1, 
			kern="op1", 
			conn="op1", 
			colnz="op1", 
			ext="op1", 
			beta1=NULL, 
			b=1, 
			c1=NULL, 
			c2=NULL, 
			z=NULL, 
			R=NULL
			)


}
