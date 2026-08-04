# Extract landscape from span.graph generated list

Extracts a landscape from an object delivered by
[`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md).
The output is an object of class 'landscape'.

## Usage

``` r
extract.graph(rl, rlist, nr)
```

## Arguments

- rl:

  Object of class 'landscape' used to generate the list, with
  [`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md).

- rlist:

  Object delivered by
  [`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md).

- nr:

  Position of the landscape in the list (rlist).

## Value

Delivers an object of class 'landscape'.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md),
[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

## Examples

``` r

data(rland)
data(landscape_change)

#Extracting the landscape of the 50th time step:

rl1 <- extract.graph(rl=rland, rlist=landscape_change, nr=50)
```
