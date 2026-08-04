# Number of components of a landscape

Returns the number of components in the landscape graph (in this case a
component is a group of patches connected by the species dispersal
distance).

## Usage

``` r
components.graph(rl)
```

## Arguments

- rl:

  Object of class 'landscape'.

## Value

Returns the number of components (groups of connected patches) of a
landscape.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

## Examples

``` r

data(rland)

components.graph(rl=rland)
#> [1] 19

#21
```
