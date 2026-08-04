# Remove a given number of patches from the landscape

Randomly removes a given number of patches from the landscape.

## Usage

``` r
removepoints(rl, nr)
```

## Arguments

- rl:

  Object of class 'landscape'.

- nr:

  Number of patches to remove.

## Value

Returns an object of class 'landscape'.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`addpoints`](https://fmestre1.github.io/MetaLandSim/reference/addpoints.md)

## Examples

``` r

data(rland)

#Checking the number of patches in the starting landscape:

rland$number.patches
#> [1] 60

#60

#Removing 10 patches from the landscape:

rl1 <- removepoints(rl=rland, nr=10)

#Checking the number of patches in the output landscape:

rl1$number.patches
#> [1] 50

#50
```
