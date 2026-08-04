# Graphical display of the landscape

Plots the landscape graph, with or without the species occupation
(respectively lists returned by
[`species.graph`](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md)
or
[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md))
and with or without the links between patches.

## Usage

``` r
plot_graph(rl, species, links)
```

## Arguments

- rl:

  Object of class 'landscape' (species=FALSE) or 'metapopulation'
  (species=TRUE).

- species:

  TRUE/FALSE, TRUE if 'x' is of class 'metapopulation' or 'FALSE' if x
  is of class 'landscape'.

- links:

  TRUE/FALSE, show links between connected patches.

## Value

Graphical display of the landscape.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`species.graph`](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md)

## Examples

``` r

data(rland)
data(occ.landscape)

#Without the species occupancy:
plot_graph(rl=rland, species=FALSE, links=FALSE)


#With the species occupancy:
plot_graph(rl=occ.landscape, species=TRUE, links=FALSE)

```
