# Add a given number of patches to a landscape

Adds a given number of patches to the landscape.

## Usage

``` r
addpoints(rl, nr)
```

## Arguments

- rl:

  Object of class 'landscape'.

- nr:

  Number of patches to be added (see 'note').

## Value

Returns an object of class 'landscape'.

## Author

Frederico Mestre and Fernando Canovas

## Note

The number of patches to be added might be impaired by the minimum
distance between points.

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`removepoints`](https://fmestre1.github.io/MetaLandSim/reference/removepoints.md)

## Examples

``` r

data(rland)

#Checking the number of patches in the starting landscape:

rland$number.patches
#> [1] 60

#60

#Adding 10 patches to a landscape:

rl1 <- addpoints(rl=rland, nr=10)

#Checking the number of patches in the output landscape:

rl1$number.patches
#> [1] 70

#70
```
