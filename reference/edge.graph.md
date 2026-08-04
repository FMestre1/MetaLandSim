# Produce an edge (links) data frame

Returns a data frame with the information on the connections between
patches (assuming binary connections).

## Usage

``` r
edge.graph(rl)
```

## Arguments

- rl:

  Object of class 'landscape'.

## Value

Produces a data frame with the information on the edges (links): the IDs
of both patches, the area, the coordinates and the Euclidean distance.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

## Examples

``` r

data(rland)

edge_df <- edge.graph(rl=rland)
```
