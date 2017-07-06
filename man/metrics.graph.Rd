\name{metrics.graph}
\alias{metrics.graph}
\title{
Computes landscape connectivity metrics
}
\description{
Computes several landscape metrics, mostly derived from graph theory or assuming a graph representation of the landscape.
}
\usage{
metrics.graph(rl, metric, dispersal.dist = NULL)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
  \item{metric}{
one of the following connectivity metrics:
\itemize{
 \item 'NC' - Number of components.
 \item 'LNK' - Number of links connecting the patches.
 \item 'SLC' - Area (in hectares) of the largest group of patches.
 \item 'MSC'- Mean area (in hectares) of the components.
 \item 'HI' - Harary Index.
 \item 'NH' - Normalized Harary Index. 
 \item 'ORD' - Landscape (graph) order.
 \item 'GD' - Landscape (graph) diameter.
 \item 'CCP' - Class coincidence probability. 
 \item 'LCP' - Landscape coincidence probability.
 \item 'ECS' - Expected cluster size. 
 \item 'AWF' - Area-weighted flux.
 \item 'IIC' - Integral index of connectivity.
 \item 'PC' - Probability of connectivity. 
 \item 'ECA' - Equivalent connected area. 
 }
}
  \item{dispersal.dist}{
Maximum dispersal distance for the binary indexes (NC, LNK, SLC, MSC, HI, NH, ORD, GD, CCP, LCP, ECS, IIC) and mean dispersal distance for the probabilistic indexes (AWF, PC, ECA). When no value is provided the function will assume the dispersal value provided by the 'landscape' object.
}
}
\details{
These metrics assume different types of links between nodes (patches). Some assume probabilistic connections between nodes (e.g. PC) while others assume binary connections (e.g. NC, SLC, LNK, IIC). Also, these metrics have several degrees of complexity, from the simpler ones (such as NC and LNK) to the more complex (such as IIC and PC). Some are purely structural; the same landscape has the same index whatever the species, while others are measures of functional, where the connectivity of a given landscape is dependent on the species (dispersal ability). Precaution must be taken when looking at the outputs produced by some of these metrics (particularly the simpler, structural ones). Regardless of being simpler to compute, the outputs might be misleading. This metrics can however be used as exploratory tools.\cr
This function was improved by the collaboration of Dr. Santiago Saura (Universidad Politecnica de Madrid).

Detail about each of the metrics:

\itemize{
 \item 'NC' - Number of components,groups of connected patches, in the landscape graph (Urban and Keitt, 2001). Patches in the same component are accessible, while patches in different components are not connected. More connected landscapes have less components. Threshold dependent (dispersal distance).
 \item 'LNK' - Number of links connecting the patches (considering that the maximum distance is the species dispersal distance and that these graphs are are binary, which means that nodes are either connected or unconnected) (Pascual-Hortal and Saura, 2006). Higher LNK implies higher connectivity. Threshold dependent (dispersal distance).
 \item 'SLC' - Area (in hectares) of the largest group of patches, or component (Pascual-Hortal and Saura, 2006). Threshold dependent (dispersal distance).
 \item 'MSC'- Mean area (in hectares) of a group of patches, or component (Pascual-Hortal and Saura, 2006). Threshold dependent (dispersal distance).
 \item 'HI' - Harary Index. Originally developed to characterize molecular graphs by Plavsic et al. (1993) it was later transposed to the landscape context by Ricotta et al. (2000). This index was considered by Ricotta et al. (2000) to be more effective from a statistical and ecological perspective. Higher HI implies higher connectivity. Threshold dependent (dispersal distance).
 \item 'NH' - Normalization of the Harary Index, facilitates analysis because this normalization will set the values between 0 and 1 and allow direct comparison of different habitat networks(Ricotta et al. 2000). Threshold dependent (dispersal distance). 
 \item 'ORD' - Order. Index originated in the graph theory and later translated into the landscape context by Urban and Keitt (2001) provides a simple structural evaluation of the graph: it is the number of patches of the component (group of patches) with more patches. Threshold dependent (dispersal distance).
 \item 'GD' - Graph diameter. Another index directly derived from graph theory, providing a simple quantification of the graph structure. The graph diameter is the maximum of all the shortest paths between the patches of an habitat network. It is computed in meters (euclidean distance), instead of number of links (such as HI, NH and IIC)(Bunn et al. 2000, Urban and Keith, 2001). Shorter diameter implies faster movement in the habitat network (Minor and Urban, 2008). Threshold dependent (dispersal distance).
 \item 'CCP' - Class coincidence probability. It is defined as the probability that two randomly chosen points within the habitat belong to the same component. Ranges between 0 and 1 (Pascual-Hortal and Saura 2006). Higher CCP implies higher connectivity. Threshold dependent (dispersal distance).
 \item 'LCP' - Landscape coincidence probability. It is defined as the probability that two randomly chosen points in the landscape (whether in an habitat patch or not) belong to the same habitat component. Ranges between 0 and 1 (Pascual-Hortal and Saura 2006). Threshold dependent (dispersal distance).
 \item 'CPL' - Characteristic path length. Mean of all the shortest paths between the network nodes (patches) (Minor and Urban, 2008). The shorter the CPL value the more connected the patches are. Threshold dependent (dispersal distance).
 \item 'ECS' - Expected cluster size. Mean cluster size of the clusters weighed by area. (O' Brien et al.,2006 and Fall et al, 2007). This represents the size of the component in which a randomly located point in an habitat patch would reside. Although it is informative regarding the area of the component, it does not provide any ecologically meaningful information regarding the total area of habitat, as an example: ECS increases with less isolated small components or patches, although the total habitat decreases(Laita et al. 2011). Threshold dependent (dispersal distance).
 \item 'AWF' - Area-weighted Flux. Evaluates the flow, weighted by area, between all pairs of patches (Bunn et al. 2000 and Urban and Keitt 2001). The probability of dispersal between two patches (pij), required by the AWF formula, was computed using pij=exp(-k*dij), where k is a constant making pij=0.5 at half the dispersal distance defined by the user. Does not depend on any distance threshold (probabilistic).  
 \item 'IIC' - Integral index of connectivity. Index developed specifically for landscapes by Pascual-Hortal and Saura (2006). It is based on habitat availability and on a binary connection model (as opposed to a probabilistic). It ranges from 0 to 1 (higher values indicating more connectivity). Threshold dependent (dispersal distance).
 \item 'PC' - Probability of connectivity. Probability that two points randomly placed in the landscape are in habitat patches that are connected, given the number of habitat patches and the connection probabilities (pij). Similar to IIC, although assuming probabilistic connections between patches (Saura and Pascual-Hortal 2007). Probability of inter-patch dispersal is computed in the same way as for AWF. Does not depend on any distance threshold (probabilistic).
\item 'ECA' -  The Equivalent Connected Area is the square root of the numerator in PC, not accounting for the total landscape area (AL) (Saura 2011a, 2011b). It is defined as '...the size of a single habitat patch (maximally connected) that would provide the same value of the probability of connectivity than the actual habitat pattern in the landscape' (Saura 2011a).
}
}
\value{
Returns the numeric value(s), corresponding to the chosen connectivity metric(s) for a given landscape.
}
\references{
Bunn, A. G., Urban, D. L., and Keitt, T. H. (2000). Landscape connectivity: a conservation application of graph theory. Journal of Environmental Management, 59(4): 265-278.

Fall, A., Fortin, M. J., Manseau, M., and O' Brien, D. (2007). Spatial graphs: principles and applications for habitat connectivity. Ecosystems, 10(3): 448-461.
 
Ivanciuc, O., Balaban, T. S., and Balaban, A. T. (1993). Design of topological indices. Part 4. Reciprocal distance matrix, related local vertex invariants and topological indices. Journal of Mathematical Chemistry, 12(1): 309-318.

Laita, A., Kotiaho, J.S., Monkkonen, M. (2011). Graph-theoretic connectivity measures: what do they tell us about connectivity? Landscape Ecology, 26: 951-967.

Minor, E. S., and Urban, D. L. (2007). Graph theory as a proxy for spatially explicit population models in conservation planning. Ecological Applications, 17(6): 1771-1782.
 
Minor, E. S., and Urban, D. L. (2008). A Graph-Theory Framework for Evaluating Landscape Connectivity and Conservation Planning. Conservation Biology, 22(2): 297-307.

O'Brien, D., Manseau, M., Fall, A., and Fortin, M. J. (2006). Testing the importance of spatial configuration of winter habitat for woodland caribou: an application of graph theory. Biological Conservation, 130(1): 70-83.

Pascual-Hortal, L., and Saura, S. (2006). Comparison and development of new graph-based landscape connectivity indices: towards the priorization of habitat patches and corridors for conservation. Landscape Ecology, 21(7): 959-967. 

Plavsic, D., Nikolic, S., Trinajstic, N., and Mihalic, Z. (1993). On the Harary index for the characterization of chemical graphs. Journal of Mathematical Chemistry, 12(1): 235-250.

Ricotta, C., Stanisci, A., Avena, G. C., and Blasi, C. (2000). Quantifying the network connectivity of landscape mosaics: a graph-theoretical approach. Community Ecology, 1(1): 89-94.

Saura, S., and Pascual-Hortal, L. (2007). A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and Urban Planning, 83(2): 91-103. 

Saura, S., Estreguil, C., Mouton, C. & Rodriguez-Freire, M. (2011a). Network analysis to assess landscape connectivity trends: application to European forests (1990-2000). Ecological Indicators 11: 407-416.

Saura, S., Gonzalez-Avila, S. & Elena-Rossello, R. (2011b). Evaluacion de los cambios en la conectividad de los bosques: el indice del area conexa equivalente y su aplicacion a los bosques de Castilla y Leon. Montes, Revista de Ambito Forestal 106: 15-21

Urban, D., and Keitt, T. (2001). Landscape connectivity: a graph-theoretic perspective. Ecology, 82(5): 1205-1218.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{rland.graph}}
}
\examples{

data(rland)

#Compute the Integral index of connectivity of a landscape:

metrics.graph (rl=rland, metric="AWF")

}
