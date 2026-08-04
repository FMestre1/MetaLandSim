# Simulate landscape dynamics over a number of time steps

This function gets an initial landscape graph and gradually applies
changes. For a good review and classification of such changes see
Bogaert et al. (2004) (not all described changes have been applied
here). Future versions of the package should include other methods to
change the landscape.

## Usage

``` r
span.graph(rl, span = 100, par1 = 'none', par2 = NULL, 
           par3 = NULL, par4 = NULL, par5 = NULL)
```

## Arguments

- rl:

  Object of class 'landscape'.

- span:

  Number of time steps (e.g. years) to simulate.

- par1:

  Parameter determining the dynamism type. One of the following (default
  'none'):

  - 'hab' percentage of the number of patches to eliminate.

  - 'dincr' minimal distance (between centroids of patches) increase
    over the simulation (in meters).

  - 'darea' percentage of increase/decrease of the mean area of patches,
    without changing SD. Patches with area \<1 square meter are deleted.

  - 'stoc' simultaneous creation and destruction of patches with
    variation in the number of created and destroyed patches.

  - 'stoc2' simultaneous creation and destruction of patches with the
    same percentage of created and destroyed patches derived from the
    number of patches of the landscape in the preceding time step.

  - 'ncsd' simultaneous creation and destruction of patches to the north
    and south of the landscape.

  - 'aggr' correlated habitat destruction.

  - 'none' no change.

  The percentage of patches to be generated or destroyed at each time
  step is not fixed (except for 'stoc2' in which case the percentage of
  created and destroyed patches is the same and directly computed from
  the number of patches in the preceeding time step, allowing to have
  landscape dynamism without change in the number of patches). For
  example if the landscape at the time step t-1 has 200 patches and the
  user wishes to set up a destruction rate of 5%, than the number of
  destroyed patches is given by a random number obtained from a Poisson
  distribution with mean 10 (5% of 200).

- par2:

  Parameter specifying details for the options in par1: percentage of
  patches do delete (if par1='hab'); distance, in meters (if
  par1='dincr'); percentage of increase/decrease (increase with negative
  sign) of the mean area of patches (if par1='darea'); percentage of
  created/destroyed patches (if par1='stoc'); percentage of created
  patches (if par1='stoc2'); 'northerndness' of created patches (if
  par1='ncsd'); percentage of destroyed patches (if par1='aggr').

- par3:

  Additional parameter specifying details for the options in par1:
  percentage of destroyed patches (if par1='stoc2'); 'southerndness' of
  destroyed patches (if par1='ncsd'); aggregation of destruction (if
  par1='aggr'). Minimum area for patch deletion, in hectares (if
  par1='darea').

- par4:

  Percentage of created patches (if par1='ncsd').

- par5:

  Percentage of destroyed patches (if par1='ncsd').

## Value

Returns a list of data frames with the nodes characteristics of a given
number of landscapes that suffer a specified change. The fields of these
data frames are the same as those from the nodes characteristics
resulting from
[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md).

## References

Bogaert, J., Ceulemans, R., & Salvador-Van Eysenrode, D. (2004).
Decision tree algorithm for detection of spatial processes in landscape
transformation. Environmental Management, 33(1): 62-73.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`simulate_graph`](https://fmestre1.github.io/MetaLandSim/reference/simulate_graph.md),
[`iterate.graph`](https://fmestre1.github.io/MetaLandSim/reference/iterate.graph.md)

## Examples

``` r

data(rland)

#Simulating a decrease of 5% in the number of patches through 100 time steps:

span1 <- span.graph(rl=rland, span=100, par1="hab", par2=5, par3=NULL, par4=NULL, par5=NULL)
#> After time step 33 all landscapes had the same number of patches.
```
