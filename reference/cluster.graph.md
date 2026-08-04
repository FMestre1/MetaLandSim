# Delivers the number of patches per cluster

Returns a data frame with the number of nodes (habitat patches) in each
component of the landscape graph (in this case a component is a group of
patches connected by the species dispersal distance).

## Usage

``` r
cluster.graph(rl)
```

## Arguments

- rl:

  Object of class 'landscape'.

## Details

The components are defined based on the species mean dispersal ability.
This implies that the connectivity model between patches is binary
(connected/not connected) as opposed to probabilistic.

## Value

This function returns a data frame with the number of patches of each
component (group of patches). The returned data frame has two fields:
cluster (Id of the component) and number of nodes (the number of nodes
of the respective component).

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

## Examples

``` r

data(rland)

cluster.graph(rl=rland)
#>    cluster number of nodes
#> 1        1               4
#> 2        2               1
#> 3        3               2
#> 4        4               9
#> 5        5               2
#> 6        6               3
#> 7        7              13
#> 8        8               2
#> 9        9               1
#> 10      10               2
#> 11      11               2
#> 12      12               3
#> 13      13               4
#> 14      14               1
#> 15      15               3
#> 16      16               4
#> 17      17               1
#> 18      18               2
#> 19      19               1

#Output:

#  cluster number of nodes
#1        1              11
#2        2               1
#3        3              13
#4        4               1
#5        5               1
#6        6              15
#7        7               2
#8        8               1
#9        9               3
#10      10               1
#11      11               1
#12      12               2
#13      13               4
#14      14               1
#15      15               1
#16      16               1
#17      17               1

```
