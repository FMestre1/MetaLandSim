# MetaLandSim

**Metapopulations Persistence and Range Expansion Simulation**

This package is described in detail in the following publication:

> Mestre, F., Cánovas, F., Pita, R., Mira, A., & Beja, P. (2016). MetaLandSim: An R package for simulating metapopulation dynamics and range expansion under landscape scenarios. *Environmental Modelling & Software*, 75, 402–406. <https://doi.org/10.1016/j.envsoft.2016.03.007>

MetaLandSim is an R package designed to simulate metapopulation dynamics and range expansion within dynamic landscapes. It provides a virtual environment that enables experimentation and simulation of ecological processes at two scales: landscape and range. By integrating concepts from metapopulation and graph theories, MetaLandSim facilitates the modeling of species persistence and dispersal across heterogeneous and changing habitats.

## Features

-   **Random Landscape Generation**: Create random landscape graphs representing habitat patches and their spatial relationships.
-   **Dynamic Landscape Simulation**: Model changes in landscape structure over time, including habitat loss or gain.
-   **Metapopulation Dynamics**: Simulate species occupancy and persistence within dynamic landscapes.
-   **Range Expansion Modeling**: Assess species range expansion scenarios, particularly in response to environmental changes like climate change.
-   **Connectivity Metrics**: Compute various graph-based connectivity metrics to evaluate landscape structure and its influence on species movement.

## Installation

MetaLandSim is available on [CRAN](https://cran.r-project.org/package=MetaLandSim) and can be installed using the following command in R:

``` r
install.packages("MetaLandSim")
```

To install the development version from GitHub:

``` r
# Install devtools if not already installed
install.packages("devtools")

# Install MetaLandSim from GitHub
devtools::install_github("FMestre1/MetaLandSim")
```

## Usage

Here's a basic example of how to use MetaLandSim:

``` r
# Load the package
library(MetaLandSim)

# Generate a random landscape with 100 patches
landscape <- rland.graph(npatch = 100, areaM = 50, distM = 1000)

# Simulate metapopulation dynamics over 20 time steps
simulation <- iterate.graph(rl = landscape, npatch = 100, nyears = 20)

# Plot the resulting landscape
plot_graph(landscape)

# Compute connectivity metrics
metrics <- metrics.graph(rl = landscape, metric = "NC")  # Number of components
```

For detailed tutorials and examples, refer to the package vignettes:

-   [Landscape Occupation Simulation in Dynamic Landscapes](https://cran.r-project.org/web/packages/MetaLandSim/vignettes/Landscape_Occupation_Simulation.html)
-   [Model Parametrization](https://cran.r-project.org/web/packages/MetaLandSim/vignettes/Model_Parametrization.html)
-   [Range Expansion Simulation](https://cran.r-project.org/web/packages/MetaLandSim/vignettes/Range_Expansion_Simulation.html)

## Documentation

Comprehensive documentation is available:

-   [CRAN Package Page](https://cran.r-project.org/package=MetaLandSim)
-   [Reference Manual (PDF)](https://cran.r-project.org/web/packages/MetaLandSim/MetaLandSim.pdf)

## Citation

If you use MetaLandSim in your research, please cite the following publication:

> Mestre, F., Cánovas, F., Pita, R., Mira, A., & Beja, P. (2016). An R package for simulating metapopulation dynamics and range expansion under landscape scenarios. *Environmental Modelling & Software*, 75, 402–406. <https://doi.org/10.1016/j.envsoft.2016.03.007>

## Authors

-   **Frederico Mestre** – [GitHub](https://github.com/FMestre1) \| [Website](https://sites.google.com/view/fmestre)
-   **Fernando Cánovas**
-   **Ricardo Pita**
-   **António Mira**
-   **Pedro Beja**

## License

MetaLandSim is licensed under the [GPL-2 or GPL-3](https://www.gnu.org/licenses/) license.
