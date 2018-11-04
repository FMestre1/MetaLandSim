\name{MetaLandSim-package}
\alias{MetaLandSim-package}
\alias{MetaLandSim}
\docType{package}
\title{
Landscape And Range Expansion Simulation
}
\description{
The package MetaLandSim is a simulation environment, allowing the generation of random landscapes, represented as graphs, the simulation of landscape dynamics, metapopulation dynamics and range expansion.\cr
The package was developed as part of the Ph.D. thesis of Frederico Mestre (SFRH/BD/73768/2010), funded by European Social Funds and the Portuguese Foundation for Science and Technology, and included in the project NETPERSIST (PTDC/AAG-MAA/3227/2012), funded by European Regional Development Fund (ERDF) through COMPETE programme and Portuguese national funds through the Portuguese Foundation for Science and Technology.\cr
MetaLandSim is intended to provide a virtual environment, enabling the experimentation and simulation of processes at two scales: landscape and range. The simulation approach, taken by MetaLandSim, presents several advantages, like allowing the test of several alternatives and the knowledge of the full system (Peck, 2004; Zurell et al. 2009). The role of simulation in landscape ecology is fundamental due to the spatial and temporal scale of the studied phenomena, which frequently hinders experimentation (Ims, 2005).\cr
Here, graph and metapopulation theories are combined, which is a broadly accepted strategy to provide a modelling framework for metapopulation dynamics (Cantwell & Forman, 1993; Bunn et al. 2000; Ricotta et al. 2000; Minor & Urban, 2008; Galpern et al. 2011). Also, several graph-based connectivity metrics can be computed from the landscape graphs. This set of metrics have been proven useful elsewhere (Urban & Keitt, 2001; Calabrese & Fagan, 2004). The graph representation of landscape has one major advantage: it effectively summarizes spatial relationships between elements and facilitates a multi-scale analysis integrating patch and landscape level analysis (Calabrese & Fagan, 2004).\cr
MLS operates at two scales, providing researchers with the possibility of:
\itemize{
 \item Landscape scale - Simulation of metapopulation occupation on a dynamic  landscape, computation of connectivity metrics.
 \item Range scale - Computes dispersal model and range expansion scenario simulation.
 }
The landscape unit, an object of class \code{\link{landscape}}, is the basic simulation unit at both these scales. At the landscape scale, the persistence of the metapopulation in a dynamic landscape is evaluated through the simulation of landscape dynamics using the function \code{\link{iterate.graph}} or \code{\link{manage_landscape_sim}}.
At the range scale the metapopulation is allowed to expand to other, empty, landscape units using \code{\link{range_expansion}}, producing an object of class \code{\link{expansion}}. The function \code{\link{range_raster}} allows the conversion of the dispersal model obtained with the previous function into a raster. Finally, also at the range scale, the user can analyse the outcome of several alternative landscapes in range expansion speed and maximum dispersal distance, using the function \code{\link{manage_expansion_sim}}.\cr
Since version 1.0 new IFM parameter estimation capabilities are available, which based upon Bayesian statistics, using the functions first developed for the paper Risk et al.(2011).\cr\cr
We thank Dr. Santiago Saura (Universidad Politecnica de Madrid) for the very useful inputs and for the R script which greatly improved the connectivity metrics capabilities of MetaLandSim.\cr \cr
Reference paper: \href{http://www.sciencedirect.com/science/article/pii/S1364815216300718}{Mestre, F.; Canovas, F.; Pita, R.; Mira, A.; Beja. P. (2016). An R package for simulating metapopulation dynamics and range expansion under environmental change. Environmental Moldelling and Software, 81: 40-44}.
}
\details{
\tabular{ll}{
Package: \tab MetaLandSim\cr
Type: \tab Package\cr
Version: \tab 1.0.5\cr
Date: \tab 2018-11-04\cr
License: GPL (>=2)\cr
}
}
\author{
Frederico Mestre, Fernando Canovas, Benjamin Risk, Ricardo Pita, Antonio Mira and Pedro Beja.

Maintainer: Frederico Mestre <mestre.frederico@gmail.com>

}
\references{
Bunn, A. G., Urban, D. L. and Keitt, T. H. (2000). Landscape connectivity: a conservation application of graph theory. Journal of Environmental Management, 59(4), 265-278.

Calabrese, J. M. and Fagan, W. F. (2004). A comparison-shopper's guide to connectivity metrics. Frontiers in Ecology and the Environment, 2(10), 529-536.

Cantwell, M. D. and Forman, R. T. (1993). Landscape graphs: ecological modelling with graph theory to detect configurations common to diverse landscapes. Landscape Ecology, 8(4), 239-255.

Galpern, P., Manseau, M. and Fall, A. (2011). Patch-based graphs of landscape connectivity: a guide to construction, analysis and application for conservation. Biological Conservation, 144(1), 44-55.

Ims, R.A. (2005). The role of experiments in landscape ecology. In: Wiens, J.A., and Moss, M.R. (eds.). Issues and Perspectives in Landscape Ecology. Cambridge University Press. pp. 70-78.

Mestre, F., Pita, R., Pauperio, J., Martins, F. M., Alves, P. C., Mira, A., & Beja, P. (2015). Combining distribution modelling and non-invasive genetics to improve range shift forecasting. Ecological Modelling, 297, 171-179.

Mestre, F., Risk, B. B., Mira, A., Beja, P., & Pita, R. (2017). A metapopulation approach to predict species range shifts under different climate change and landscape connectivity scenarios. Ecological Modelling, 359, 406-414.

Minor, E. S. and Urban, D. L. (2008). A Graph Theory Framework for Evaluating Landscape Connectivity and Conservation Planning. Conservation Biology, 22(2), 297-307.

Peck, S. L. (2004). Simulation as experiment: a philosophical reassessment for biological modelling. Trends in Ecology & Evolution, 19(10), 530-534.

Ricotta, C., Stanisci, A., Avena, G. C., and Blasi, C. (2000). Quantifying the network connectivity of landscape mosaics: a graph-theoretical approach. Community Ecology, 1(1), 89-94.

Risk, B. B., De Valpine, P., Beissinger, S. R. (2011). A robust design formulation of the incidence function model of metapopulation dynamics applied to two species of rails. Ecology, 92(2), 462-474.

Urban, D. and Keitt, T. (2001). Landscape connectivity: a graph-theoretic perspective. Ecology, 82(5), 1205-1218.

Zurell, D., Berger, U., Cabral, J.S., Jeltsch, F., Meynard, C.N., Munkemuller, T., Nehrbass, N., Pagel, J., Reineking, B., Schroder, B. and Grimm, V. (2009). The virtual ecologist approach: simulating data and observers. Oikos, 119(4), 622-635.
}
