# Classify patches in clusters

reclassify clusters of a landscape according to a given mean dispersal
distance.

## Usage

``` r
cluster.id(rl)
```

## Arguments

- rl:

  Object of class 'landscape'.

## Details

After changing the landscape some components (groups of connected
patches) might suffer changes (e.g. the removal of patches might split
components). This function re-attributes a code to each patch,
identifying the groups of connected patches (components), after this
type of disturbance to the habitat network.Mainly to be used internally.

## Value

Returns the same landscape object, with the clusters reclassified.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

## Examples

``` r

data(rland)

#After removing 30 (50%) of the patches of a landscape:

rland2 <- removepoints(rl=rland, nr=35)

#A reclassification might be needed to identify components: 

rland2 <- cluster.id(rl=rland2)

#After removing 35 patches, there's a different number of components:

components.graph(rl=rland) 
#> [1] 19

#21

components.graph(rl=rland2) 
#> [1] 20

#16
```
