# Summarize 'metapopulation' class objects

This function summarizes a
[`metapopulation`](https://fmestre1.github.io/MetaLandSim/reference/metapopulation.md)
class object.

## Usage

``` r
summary_metapopulation(object)
```

## Arguments

- object:

  Object of class
  [`metapopulation`](https://fmestre1.github.io/MetaLandSim/reference/metapopulation.md)

## Details

This function can be used to retrieve basic information on the objects
of class 'metapopulation'.

## Value

Returns a data frame with the following information on a
[`metapopulation`](https://fmestre1.github.io/MetaLandSim/reference/metapopulation.md)
class object:

- landscape area (hectares) :

  Landscape mosaic area, in hectares

- number of patches :

  Number of patches in the landscape

- mean patch area (hectares) :

  Mean patch area, in hectares

- SD patch area :

  SD of the patch area

- mean distance amongst patches (meters) :

  Mean inter-patch distance, in meters

- minimum distance amongst patches (meters) :

  Minimum inter-patch distance, in meters

- species occurrence - snapshot :

  Occupation data of the focal species, numbered from 1 to the number of
  snapshots

## Author

Frederico Mestre and Fernando Canovas

## Note

The minimum distance between patches is different from that given in the
object of class 'landscape', in the slot 'minimum.distance'. This is
because this output is computed from the landscape structure and the one
in the 'landscape' object was the parameter used to built the landscape.
The minimum inter-patch distance given as a parameter in the function
[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)
will consider distance between patch centroids. The minimum inter-patch
distance returned here considers the edge-to-edge distance, so this
might be smaller that the parameter of
[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md).
In order to see the difference between centroid-to-centroid and
edge-to-edge inter-patch distance compute both using the
[`matrix.graph`](https://fmestre1.github.io/MetaLandSim/reference/matrix.graph.md)
function (methods are 'centr_distance' and 'euc_distance',
respectively).

## See also

[`species.graph`](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md),
[`metapopulation`](https://fmestre1.github.io/MetaLandSim/reference/metapopulation.md),
[`matrix.graph`](https://fmestre1.github.io/MetaLandSim/reference/matrix.graph.md)

## Examples

``` r

data(occ.landscape)
data(occ.landscape2)

summary_metapopulation(object=occ.landscape)
#>                                             Value
#> landscape area (hectares)                 100.000
#> number of patches                          60.000
#> mean patch area (hectares)                  0.061
#> SD patch area                               0.041
#> mean distance amongst patches (meters)    528.345
#> minimum distance amongst patches (meters)  51.780
#> species occurrence - snapshot 1            50.000

#                                            Value
#landscape area (hectares)                 100.000
#number of patches                          60.000
#mean patch area (hectares)                  0.061
#SD patch area                               0.041
#mean distance amongst patches (meters)    528.345
#minimum distance amongst patches (meters)  51.780
#species occurrence - snapshot 1            50.000


summary_metapopulation(object=occ.landscape2)
#>                                             Value
#> landscape area (hectares)                 100.000
#> number of patches                          60.000
#> mean patch area (hectares)                  0.069
#> SD patch area                               0.039
#> mean distance amongst patches (meters)    521.717
#> minimum distance amongst patches (meters)  45.905
#> species occurrence - snapshot 1            50.000
#> species occurrence - snapshot 2            58.333
#> species occurrence - snapshot 3            61.667
#> species occurrence - snapshot 4            61.667
#> species occurrence - snapshot 5            58.333
#> species occurrence - snapshot 6            60.000
#> species occurrence - snapshot 7            70.000
#> species occurrence - snapshot 8            68.333
#> species occurrence - snapshot 9            68.333
#> species occurrence - snapshot 10           56.667

#                                            Value
#landscape area (hectares)                 100.000
#number of patches                          60.000
#mean patch area (hectares)                  0.069
#SD patch area                               0.039
#mean distance amongst patches (meters)    521.717
#minimum distance amongst patches (meters)  45.905
#species occurrence - snapshot 1            50.000
#species occurrence - snapshot 2            58.333
#species occurrence - snapshot 3            61.667
#species occurrence - snapshot 4            61.667
#species occurrence - snapshot 5            58.333
#species occurrence - snapshot 6            60.000
#species occurrence - snapshot 7            70.000
#species occurrence - snapshot 8            68.333
#species occurrence - snapshot 9            68.333
#species occurrence - snapshot 10           56.667
```
