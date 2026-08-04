# Summarize 'landscape' class objects

This function summarizes a
[`landscape`](https://fmestre1.github.io/MetaLandSim/reference/landscape.md)
class object.

## Usage

``` r
summary_landscape(object)
```

## Arguments

- object:

  Object of class
  [`landscape`](https://fmestre1.github.io/MetaLandSim/reference/landscape.md)

## Details

This function can be used to retrieve basic information on the objects
of class 'landscape'.

## Value

Returns a data frame with the following information on a
[`landscape`](https://fmestre1.github.io/MetaLandSim/reference/landscape.md)
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

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`landscape`](https://fmestre1.github.io/MetaLandSim/reference/landscape.md),
[`matrix.graph`](https://fmestre1.github.io/MetaLandSim/reference/matrix.graph.md)

## Examples

``` r

data(rland)

summary_landscape(object=rland)
#>                                             Value
#> landscape area (hectares)                 100.000
#> number of patches                          60.000
#> mean patch area (hectares)                  0.223
#> SD patch area                               0.091
#> mean distance amongst patches (meters)    492.910
#> minimum distance amongst patches (meters)  20.872

#                                            Value
#landscape area (hectares)                 100.000
#number of patches                          60.000
#mean patch area (hectares)                  0.061
#SD patch area                               0.041
#mean distance amongst patches (meters)    528.345
#minimum distance amongst patches (meters)  51.780
```
