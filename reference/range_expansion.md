# Computes a range expansion model

This function returns the expansion probability, from a landscape with a
given set of parameters. This can subsequently be converted in a
dispersal model by the function
[`range_raster`](https://fmestre1.github.io/MetaLandSim/reference/range_raster.md).

## Usage

``` r
range_expansion(rl, percI, param, b, tsteps, iter, plot)
```

## Arguments

- rl:

  Object of class 'landscape'. Starting landscape for the expansion
  procedure.

- percI:

  Pecentage of patch occupancy in the starting landscape.

- param:

  Parameter data frame delivered by
  [`parameter.estimate`](https://fmestre1.github.io/MetaLandSim/reference/parameter.estimate.md),
  including:

  - alpha - Parameter relating extinction with distance.

  - y - Parameter y in the colonization probability.

  - e - Parameter defining the extinction probability in a patch of unit
    area.

  - x - Parameter scaling extinction risk with patch area.

- b:

  Parameter scaling emigration with patch area (if conn='op1' or 'op2')
  in [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).
  By default, equal to 1.

- tsteps:

  Number of time steps to simulate (e.g. years).

- iter:

  Number of iterations of the simulation procedure.

- plot:

  Plot results.

## Details

The expansion algorithm has been improved, since the paper Mestre et al.
(2017) describing the package was published. Now, instead of the
transition between adjacent landscape units being dictated by the
occupation of a spurious node (representing the margin through which the
expansion takes place) a somewhat more realistic approach is followed.
If, during the metapopulational dynamics simulation, any patch located
between the landscape unit (LU) margin and a parallel line placed at a
distance equivalent to half of the mean dispersal ability of the species
is occupied, than the algorithm assumes that the species will have the
ability to go across to the next LU. In this new empty LU initial
occupation is defined as follows: a new line is placed, with a spacing
equivalent to half the dispersal ability of the species. In the area
defined by the margins of the LU and this line the species will occupy
in the same proportion as in the preceding LU.  
After version 2.0.0, the output, rather than considering distinct
dispersal probabilities in all four cardinal directions (as in previous
versions), considers the same probability of dispersal from a current
presence in all directions. This does not change the results in any
meaningfull way given that these kinds of simulations require many
iterations in which the distinctions between the dispersal to all four
directions was diluted.

## Value

This function returns a data frame with the proportion of occupations at
several distances from the closest occupied landscape mosaic. After
version 2.0.0 the package uses the same dispersal probability in all
directions relative to current presences. the data frame has the
following columns:

- DISTANCE - Distance (mapsize x number of landscapes).

- OCCUPATION - How many times did the landscape at this distance got
  occupied by the species (from a total of 'iter' repetitions).

- PROPORTION - Proportion of occupation for the landscape at this
  distance (OCCUPATION/iter).

- TIME STEP - The average time steps at which a given distance is
  occupied.

## Author

Frederico Mestre and Fernando Canovas

## References

Mestre, F., Risk, B., Mira, A., Beja, P., Pita, R. (2017)
\<doi:10.1016/j.ecolmodel.2017.06.013\>

## Note

Depending on computing power and number of iterations (parameter 'iter')
this function can take some time to run.

## See also

[`range_raster`](https://fmestre1.github.io/MetaLandSim/reference/range_raster.md)

## Examples

``` r

if (FALSE) { # \dontrun{
#Produce a model of range expansion:
#Note: this function should be run with >100 iterations (parameter "iter").

data(rland)
data(param2)

rg_exp1 <- range_expansion(rl=rland, percI=80, param=param2, b=1, tsteps=100, iter=100)

} # }
```
