# Simulate range expansion simulation

This function produces dispersal scenarios, considering different
habitat networks properties, evaluating the variation in dispersal speed
and dispersal maximum distance (of range expansion).

## Usage

``` r
manage_expansion_sim(mapsize, dist_m, areaM, areaSD, Npatch,percI, 
                    param, b=1, tsteps, iter, variable,var_min,var_max,by)
```

## Arguments

- mapsize:

  Landscape mosaic side length, in meters. To be internally passed to
  [`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

- dist_m:

  Minimum distance between patches (centroid).To be internally passed to
  [`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

- areaM:

  Mean area (in hectares). To be internally passed to
  [`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

- areaSD:

  SD of the area of patches, in order to give variability to the patches
  area. To be internally passed to
  [`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

- Npatch:

  Number of patches. To be internally passed to
  [`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

- percI:

  Percentage of patch occupancy in the starting landscape. To be
  internally passed to
  [`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md)

- param:

  Parameter data frame delivered by
  [`parameter.estimate`](https://fmestre1.github.io/MetaLandSim/reference/parameter.estimate.md).
  To be internally passed to
  [`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md)
  It includes:

  - alpha - Parameter relating extinction with distance.

  - y - Parameter y in the colonization probability.

  - e - Parameter defining the extinction probability in a patch of unit
    area.

  - x - Parameter scaling extinction risk with patch area.

- b:

  Parameter scaling emigration with patch area (if conn='op1' or 'op2')
  in [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).
  By default, equal to 1.To be internally passed to
  [`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md)

- tsteps:

  Number of time steps to simulate (e.g. years).

- iter:

  Number of iterations of the simulation procedure.

- variable:

  Landscape graph characteristic to be varied. One of the following:

  - area - Mean patch area (in hectares).

  - dist - Minimum distance between patches (centroid).

  - npatch - Number of patches.

  - sizevar - SD of the area of patches.

- var_min:

  Minimum value the changing variable can assume (beware of units used
  in each case).

- var_max:

  Maximum value the changing variable can assume (beware of units used
  in each case).

- by:

  Rate of variation of the changing variable.

## Details

For details regarding the arguments that are to be internally passed to
other functions, see the respective functions. Any of the arguments
dist_m, areaM, areaSD, Npatch would be unnecessary if the respective
variable is the one to be evaluated (it depends on the parameter
variable).

## Value

Returns a list of eight data frames. For the first four data frames
(NORTH, SOUTH, EAST and WEST) each data frame's first column is the name
of the variable to be changed. The other two columns are:

- MEAN EXPANSION SPEED :

  Expansion speed in each simulated scenario. Speed given in km/time
  step

- MAXIMUM EXPANSION DISTANCE :

  Maximum distance of the expanded range, from an occupied site. Given
  in km.

  
The other four data frames have detailed information on the simulations
for each of the values of parameter "variable". The first column has the
distance (in km), and each of the following columns has the time step at
which each distance was colonized for each of the simulations.

## Warning

This function might be time consuming, and the code is experimental and
should be improved in future versions of MetaLandSim.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md),
[`expansion`](https://fmestre1.github.io/MetaLandSim/reference/expansion.md)

## Examples

``` r

if (FALSE) { # \dontrun{

data(param1)


sim_range <- manage_expansion_sim(mapsize=1000, dist_m=0, areaM, areaSD=0.001, 
                                  patch=300,percI=50, param=param1, b=1, 
                                  tsteps=100, iter=100,variable="area",var_min=0.01,
                                  var_max=0.6,by=0.1)
} # }
```
