\name{metrics.graph}
\alias{metrics.graph}
\title{
Computes landscape connectivity metrics
}
\description{
Computes several landscape metrics, mostly derived from graph theory or assuming a graph representation of the landscape.
}
\usage{
metrics.graph(rl, metric)
}
\arguments{
  \item{rl}{
Object of class 'landscape'.
}
  \item{metric}{
one of the following connectivity metrics:
\itemize{
 \item 'NC' - Number of components,groups of connected patches, in the landscape graph (Urban and Keitt, 2001).
 \item 'LNK' - Number of links connecting the patches (considering that the maximum distance is the species dispersal distance and that these graphs are are binary, which means that nodes are either connected or unconnected) (Pascual-Hortal and Saura, 2006).
 \item 'SLC' - Size (in hectares) of the largest group of patches, or component (Pascual-Hortal and Saura, 2006). 
 \item 'MSC'- Mean size (in hectares) of a group of patches, or component (Pascual-Hortal and Saura, 2006).
 \item 'HI' - Harary Index. Originally developed to characterize molecular graphs by Plavsic et al. (1993) it was later transposed to the landscape context by Ricotta et al. (2000). This index was considered by Ricotta et al. (2000) to be more effective from a statistical and ecological perspective.
 \item 'NH' - Normalization of the Harary Index, facilitates analysis because this normalization will set the values between 0 and 1 (Ricotta et al. 2000). 
 \item 'ORD' - Order. Index originated in the graph theory and later translated into the landscape context by Urban and Keitt (2001) provides a simple structural evaluation of the graph: it is the number of patches of the component (group of patches) with more patches.
 \item 'GD' - Graph diameter. Another index directly derived from graph theory, providing a simple quantification of the graph structure. The graph diameter or 'longest shortest path' between the two most distant patches in a network (Bunn et al. 2000, Urban and Keith, 2001).
 \item 'CCP' - Class coincidence probability. It is defined as the probability that two randomly chosen points within the habitat belong to the same component. Ranges between 0 and 1 (Pascual-Hortal and Saura 2006).
 \item 'LCP' - Landscape coincidence probability. It is defined as the probability that two randomly chosen points in the landscape (whether in an habitat patch or not) belong to the same habitat component. Ranges between 0 and 1 (Pascual-Hortal and Saura 2006).
 \item 'CPL' - Characteristic path length. Mean of all the shortest paths between all reachable nodes (patches) (Minor and Urban, 2008).
 \item 'ECS' - Expected cluster size. Mean cluster size of the clusters weighed by area. (O' Brien et al.,2006 and Fall et al, 2007).
 \item 'AWF' - Area-weighted Flux. Evaluates the flow, weighted by area, between all pairs of patches (Bunn et al. 2000 and Urban and Keitt 2001).
 \item 'IIC' - Integral index of connectivity. Index developed specifically for landscapes by Pascual-Hortal and Saura 2006. It is based on habitat availability and on a binary connection model (as opposed to a probabilistic). This index is considered to be the most adequate because it allows the evaluation of spatial changes and patch prioritization. It ranges from 0 to 1 (higher values indicating more connectivity).
 \item 'PC' - Probability of connectivity. Similar to IIC, although assuming probabilistic connections between patches (Saura and Pascual-Hortal 2007).
 }
}
}
\details{
These metrics assume different types of links between nodes (patches). Some assume probabilistic connections between nodes (e.g. PC) while others assume binary connections (e.g. NC, SLC, LNK, IIC). Some of these metrics are very simple, while others are more complex.  From the simpler ones (such as NC and LNK) to the more complex (such as IIC and PC). Some of these measures of connectivity are purely structural; the same landscape has the same index whatever the species, while others are measures of functional, where the connectivity of a given landscape is dependent on the species. Precaution must be taken when looking at the outputs produced by some of these metrics (particularly the simpler, structural ones). Regardless of being simpler to compute, the outputs might be misleading. This metrics can however be used as exploratory tools.
}
\value{
Returns a numeric value, which corresponds to the value of the chosen connectivity metric for the given landscape.
}
\references{
Bunn, A. G., Urban, D. L., and Keitt, T. H. (2000). Landscape connectivity: a conservation application of graph theory. Journal of Environmental Management, 59(4): 265-278.

Fall, A., Fortin, M. J., Manseau, M., and O' Brien, D. (2007). Spatial graphs: principles and applications for habitat connectivity. Ecosystems, 10(3): 448-461.
 
Ivanciuc, O., Balaban, T. S., and Balaban, A. T. (1993). Design of topological indices. Part 4. Reciprocal distance matrix, related local vertex invariants and topological indices. Journal of Mathematical Chemistry, 12(1): 309-318.

Minor, E. S., and Urban, D. L. (2007). Graph theory as a proxy for spatially explicit population models in conservation planning. Ecological Applications, 17(6): 1771-1782.
 
Minor, E. S., and Urban, D. L. (2008). A Graph-Theory Framework for Evaluating Landscape Connectivity and Conservation Planning. Conservation Biology, 22(2): 297-307.

O'Brien, D., Manseau, M., Fall, A., and Fortin, M. J. (2006). Testing the importance of spatial configuration of winter habitat for woodland caribou: an application of graph theory. Biological Conservation, 130(1): 70-83.

Pascual-Hortal, L., and Saura, S. (2006). Comparison and development of new graph-based landscape connectivity indices: towards the priorization of habitat patches and corridors for conservation. Landscape Ecology, 21(7): 959-967. 

Plavsic, D., Nikolic, S., Trinajstic, N., and Mihalic, Z. (1993). On the Harary index for the characterization of chemical graphs. Journal of Mathematical Chemistry, 12(1): 235-250.

Ricotta, C., Stanisci, A., Avena, G. C., and Blasi, C. (2000). Quantifying the network connectivity of landscape mosaics: a graph-theoretical approach. Community Ecology, 1(1): 89-94.

Saura, S., and Pascual-Hortal, L. (2007). A new habitat availability index to integrate connectivity in landscape conservation planning: comparison with existing indices and application to a case study. Landscape and Urban Planning, 83(2): 91-103. 

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

#0.03703083

}
