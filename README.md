# MetaLandSim

**Metapopulations Persistence and Range Expansion Simulation**

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
#1. Installing the CRAN version

install.packages("MetaLandSim")

#1. Installing the development version

# Install devtools if not already installed
install.packages("devtools")

# Install MetaLandSim from GitHub
devtools::install_github("FMestre1/MetaLandSim")
```
