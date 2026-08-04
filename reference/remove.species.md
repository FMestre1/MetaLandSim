# Remove the species occupancy from the landscape

This function converts an object of class 'metapopulation' (with the
species occupancy) in a object of class 'landscape' (without the species
occupancy).

## Usage

``` r
remove.species(sp)
```

## Arguments

- sp:

  Object of class 'metapopulation'.

## Value

Delivers an object of class 'landscape'.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`species.graph`](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md)

## Examples

``` r

data(occ.landscape)

rl1 <- remove.species(sp=occ.landscape)
```
