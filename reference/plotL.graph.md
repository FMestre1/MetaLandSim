# Plot one landscape of the list created by span.graph

Plots a given landscape of a landscape sequence from
[`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md).

## Usage

``` r
plotL.graph(rl, rlist, nr, species, links, ...)
```

## Arguments

- rl:

  Object of class 'landscape'.

- rlist:

  List returned by
  [`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md).

- nr:

  index of the landscape to display graphically.

- species:

  TRUE/FALSE, TRUE if 'rl' is of class 'metapopulation' or 'FALSE' if rl
  is of class 'landscape'.

- links:

  TRUE/FALSE, show links between connected patches.

- ...:

  Other arguments.

## Value

Graphical display of the landscape.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`plot_graph`](https://fmestre1.github.io/MetaLandSim/reference/plot_graph.md),
[`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md),
[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

## Examples

``` r

data(rland)
data(landscape_change)

plotL.graph(rl=rland, rlist=landscape_change, nr=50, species=FALSE, links=FALSE)

```
