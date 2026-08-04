# Landscape And Range Expansion Simulation

The package MetaLandSim is a simulation environment, allowing the
generation of random landscapes, represented as graphs, the simulation
of landscape dynamics, metapopulation dynamics and range expansion.  
The package was developed as part of the Ph.D. thesis of Frederico
Mestre (SFRH/BD/73768/2010), funded by European Social Funds and the
Portuguese Foundation for Science and Technology, and included in the
project NETPERSIST (PTDC/AAG-MAA/3227/2012), funded by European Regional
Development Fund (ERDF) through COMPETE programme and Portuguese
national funds through the Portuguese Foundation for Science and
Technology.  
It is intended to provide a virtual environment, enabling the
experimentation and simulation of processes at two scales: landscape and
range. The simulation approach, taken by MetaLandSim, presents several
advantages, like allowing the test of several alternatives and the
knowledge of the full system (Peck, 2004; Zurell et al. 2009). The role
of simulation in landscape ecology is fundamental due to the spatial and
temporal scale of the studied phenomena, which frequently hinders
experimentation (Ims, 2005).  
Here, graph and metapopulation theories are combined, which is a broadly
accepted strategy to provide a modelling framework for metapopulation
dynamics (Cantwell & Forman, 1993; Bunn et al. 2000; Ricotta et al.
2000; Minor & Urban, 2008; Galpern et al. 2011). Also, several
graph-based connectivity metrics can be computed from the landscape
graphs. This set of metrics have been proven useful elsewhere (Urban &
Keitt, 2001; Calabrese & Fagan, 2004). The graph representation of
landscape has one major advantage: it effectively summarizes spatial
relationships between elements and facilitates a multi-scale analysis
integrating patch and landscape level analysis (Calabrese & Fagan,
2004).  
MetaLandSim operates at two scales, providing researchers with the
possibility of:

- Landscape scale - Simulation of metapopulation occupation on a dynamic
  landscape, computation of connectivity metrics.

- Range scale - Computes dispersal model and range expansion scenario
  simulation.

The landscape unit, an object of class
[`landscape`](https://fmestre1.github.io/MetaLandSim/reference/landscape.md),
is the basic simulation unit at both these scales. At the landscape
scale, the persistence of the metapopulation in a dynamic landscape is
evaluated through the simulation of landscape dynamics using the
function
[`iterate.graph`](https://fmestre1.github.io/MetaLandSim/reference/iterate.graph.md)
or
[`manage_landscape_sim`](https://fmestre1.github.io/MetaLandSim/reference/manage_landscape_sim.md).
At the range scale the metapopulation is allowed to expand to other,
empty, landscape units using
[`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md),
producing an object of class
[`expansion`](https://fmestre1.github.io/MetaLandSim/reference/expansion.md).
The function
[`range_raster`](https://fmestre1.github.io/MetaLandSim/reference/range_raster.md)
allows the conversion of the dispersal model obtained with the previous
function into a raster. Finally, also at the range scale, the user can
analyse the outcome of several alternative landscapes in range expansion
speed and maximum dispersal distance, using the function
[`manage_expansion_sim`](https://fmestre1.github.io/MetaLandSim/reference/manage_expansion_sim.md).  
Since version 1.0 new IFM parameter estimation capabilities are
available, which based upon Bayesian statistics, using the functions
first developed for the paper Risk et al.(2011).  
  
We thank Dr. Santiago Saura (Universidad Politecnica de Madrid) for the
very useful inputs and for the R script which greatly improved the
connectivity metrics capabilities of MetaLandSim.  
  
After version 2.0.0 MetaLandSim had a few major changes: 1) There is no
Graphic User Interface, the user will have to resort solely to the usual
R user interface; 2) It does not use GRASS, resorting uniquely to R
packages to conduct the simulations (mainly terra); 3) It depends on
much less packages (after removing rgrass7, maptools, rgeos, raster,
tcltk and fgui); 4) There were some major changes to the functions
[`range_raster`](https://fmestre1.github.io/MetaLandSim/reference/range_raster.md)
and
[`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md).
In what concerns
[`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md)
the output, rather than considering distinct dispersal probabilities in
all four cardinal directions (as in previous versions), considers the
same probability of dispersal from a current presence in all directions.
This has implications in the
[`range_raster`](https://fmestre1.github.io/MetaLandSim/reference/range_raster.md)
function, that converts the dispersal probability to a raster. However,
this does not change the results in any meaningfull way given that these
kinds of simulations require many iterations in which the distinctions
between the dispersal to all four directions was diluted.

## Details

|                     |             |
|---------------------|-------------|
| Package:            | MetaLandSim |
| Type:               | Package     |
| Version:            | 2.0.0       |
| Date:               | 2022-01-12  |
| License: GPL (\>=2) |             |

## Author

Frederico Mestre, Fernando Canovas, Benjamin Risk, Ricardo Pita, Antonio
Mira and Pedro Beja.

Maintainer: Frederico Mestre \<mestre.frederico@gmail.com\>

## References

Bunn, A. G., Urban, D. L. and Keitt, T. H. (2000). Landscape
connectivity: a conservation application of graph theory. Journal of
Environmental Management, 59(4), 265-278.

Calabrese, J. M. and Fagan, W. F. (2004). A comparison-shopper's guide
to connectivity metrics. Frontiers in Ecology and the Environment,
2(10), 529-536.

Cantwell, M. D. and Forman, R. T. (1993). Landscape graphs: ecological
modelling with graph theory to detect configurations common to diverse
landscapes. Landscape Ecology, 8(4), 239-255.

Galpern, P., Manseau, M. and Fall, A. (2011). Patch-based graphs of
landscape connectivity: a guide to construction, analysis and
application for conservation. Biological Conservation, 144(1), 44-55.

Ims, R.A. (2005). The role of experiments in landscape ecology. In:
Wiens, J.A., and Moss, M.R. (eds.). Issues and Perspectives in Landscape
Ecology. Cambridge University Press. pp. 70-78.

Mestre, F., Pita, R., Pauperio, J., Martins, F. M., Alves, P. C., Mira,
A., & Beja, P. (2015). Combining distribution modelling and non-invasive
genetics to improve range shift forecasting. Ecological Modelling, 297,
171-179.

Mestre, F., Risk, B. B., Mira, A., Beja, P., & Pita, R. (2017). A
metapopulation approach to predict species range shifts under different
climate change and landscape connectivity scenarios. Ecological
Modelling, 359, 406-414.

Mestre, F., Pita, R., Mira, A., Beja, P. (2020). Species traits, patch
turnover and successional dynamics: When does intermediate disturbance
favour metapopulation occupancy?. BMC Ecology.

Minor, E. S. and Urban, D. L. (2008). A Graph Theory Framework for
Evaluating Landscape Connectivity and Conservation Planning.
Conservation Biology, 22(2), 297-307.

Peck, S. L. (2004). Simulation as experiment: a philosophical
reassessment for biological modelling. Trends in Ecology & Evolution,
19(10), 530-534.

Ricotta, C., Stanisci, A., Avena, G. C., and Blasi, C. (2000).
Quantifying the network connectivity of landscape mosaics: a
graph-theoretical approach. Community Ecology, 1(1), 89-94.

Risk, B. B., De Valpine, P., Beissinger, S. R. (2011). A robust design
formulation of the incidence function model of metapopulation dynamics
applied to two species of rails. Ecology, 92(2), 462-474.

Urban, D. and Keitt, T. (2001). Landscape connectivity: a
graph-theoretic perspective. Ecology, 82(5), 1205-1218.

Zurell, D., Berger, U., Cabral, J.S., Jeltsch, F., Meynard, C.N.,
Munkemuller, T., Nehrbass, N., Pagel, J., Reineking, B., Schroder, B.
and Grimm, V. (2009). The virtual ecologist approach: simulating data
and observers. Oikos, 119(4), 622-635.
