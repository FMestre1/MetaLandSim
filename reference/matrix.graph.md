# Returning a matrix with information on connections between patches

Based on a landscape graph, this function allows the creation of a
matrix of Euclidean distances (straight-line pairwise distance between
the margins of all the patches), matrix of topological distances
(minimum number of connections between any two patches) and adjacency
matrix (this a matrix of 0 and 1, showing the adjacency between any two
patches, where 0 means that the patches are not connected and 1 means
that the patches are connected).

## Usage

``` r
matrix.graph(rl, mat)
```

## Arguments

- rl:

  Object of class 'landscape'.

- mat:

  mat - one of the following:

  - 'euc_distance' - euclidian distance between patches (edge-to-edge).

  - 'centr_distance' - euclidian distance between patches
    (centroid-to-centroid).

  - 'adjacency' - adjacency matrix, with values d_ij, taking value 0 if
    patches i and j are not connected and value 1 if those patches are
    connected.

  - 'top_matrix' - topological distance, with values d_ij, where the
    value d is the minimum number of connections between the patches i
    and j. Topological distance is defined as the minimum number of
    links between patches i and j.

## Value

This function returns a matrix (each one of the specified matrices:
Euclidean distance, topological distance and adjacency matrix).

## Author

Frederico Mestre and Fernando Canovas

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md)

## Examples

``` r

data(rland)

#Computing matrix of topological distances:

matrix.graph(rl=rland, mat="top_matrix")
#>     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> 1   0  1  2  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 2  NA  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 3  NA NA  0  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 4  NA NA NA  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 5  NA NA NA NA  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 6  NA NA NA NA NA  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 7  NA NA NA NA NA NA  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 8  NA NA NA NA NA NA NA  0  2  2  1  1  3  2  2  3  0  0  0  0  0  0  0  0  0
#> 9  NA NA NA NA NA NA NA NA  0  4  3  1  1  4  2  5  0  0  0  0  0  0  0  0  0
#> 10 NA NA NA NA NA NA NA NA NA  0  1  3  5  2  4  3  0  0  0  0  0  0  0  0  0
#> 11 NA NA NA NA NA NA NA NA NA NA  0  2  4  1  3  2  0  0  0  0  0  0  0  0  0
#> 12 NA NA NA NA NA NA NA NA NA NA NA  0  2  3  1  4  0  0  0  0  0  0  0  0  0
#> 13 NA NA NA NA NA NA NA NA NA NA NA NA  0  5  3  6  0  0  0  0  0  0  0  0  0
#> 14 NA NA NA NA NA NA NA NA NA NA NA NA NA  0  4  1  0  0  0  0  0  0  0  0  0
#> 15 NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  5  0  0  0  0  0  0  0  0  0
#> 16 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  0  0  0  0  0  0  0  0  0
#> 17 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  1  0  0  0  0  0  0  0
#> 18 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  0  0  0  0  0  0  0
#> 19 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  1  1  0  0  0  0
#> 20 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  2  0  0  0  0
#> 21 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  0  0  0  0
#> 22 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  1  3  3
#> 23 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  4  4
#> 24 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  3
#> 25 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0
#> 26 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 27 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 28 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 29 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 30 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 31 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 32 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 33 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 34 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 35 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 36 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 37 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 38 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 39 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 40 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 41 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 42 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 43 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 44 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 45 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 46 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 47 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 48 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 49 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 50 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 51 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 52 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 53 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 54 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 55 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 56 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 57 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 58 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 59 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 60 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>    26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
#> 1   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 2   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 3   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 4   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 5   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 6   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 7   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 8   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 9   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 10  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 11  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 13  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 14  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 15  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 16  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 17  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 18  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 19  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 20  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 21  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 22  2  3  3  1  2  2  1  3  4  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 23  3  4  4  2  3  3  2  4  5  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 24  2  2  3  4  5  1  2  1  4  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 25  1  1  1  4  5  2  2  2  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 26  0  1  1  3  4  1  1  2  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 27 NA  0  2  4  5  2  2  1  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 28 NA NA  0  4  5  2  2  3  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 29 NA NA NA  0  1  3  2  4  5  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 30 NA NA NA NA  0  4  3  5  6  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 31 NA NA NA NA NA  0  1  1  3  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 32 NA NA NA NA NA NA  0  2  3  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 33 NA NA NA NA NA NA NA  0  3  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 34 NA NA NA NA NA NA NA NA  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 35 NA NA NA NA NA NA NA NA NA  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 36 NA NA NA NA NA NA NA NA NA NA  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 37 NA NA NA NA NA NA NA NA NA NA NA  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#> 38 NA NA NA NA NA NA NA NA NA NA NA NA  0  1  0  0  0  0  0  0  0  0  0  0  0
#> 39 NA NA NA NA NA NA NA NA NA NA NA NA NA  0  0  0  0  0  0  0  0  0  0  0  0
#> 40 NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  1  0  0  0  0  0  0  0  0  0
#> 41 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  0  0  0  0  0  0  0  0  0
#> 42 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  1  1  0  0  0  0  0  0
#> 43 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  1  0  0  0  0  0  0
#> 44 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  0  0  0  0  0  0
#> 45 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  1  3  2  0  0
#> 46 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  2  1  0  0
#> 47 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  1  0  0
#> 48 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  0  0
#> 49 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0  0
#> 50 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA  0
#> 51 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 52 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 53 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 54 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 55 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 56 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 57 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 58 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 59 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 60 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>    51 52 53 54 55 56 57 58 59 60
#> 1   0  0  0  0  0  0  0  0  0  0
#> 2   0  0  0  0  0  0  0  0  0  0
#> 3   0  0  0  0  0  0  0  0  0  0
#> 4   0  0  0  0  0  0  0  0  0  0
#> 5   0  0  0  0  0  0  0  0  0  0
#> 6   0  0  0  0  0  0  0  0  0  0
#> 7   0  0  0  0  0  0  0  0  0  0
#> 8   0  0  0  0  0  0  0  0  0  0
#> 9   0  0  0  0  0  0  0  0  0  0
#> 10  0  0  0  0  0  0  0  0  0  0
#> 11  0  0  0  0  0  0  0  0  0  0
#> 12  0  0  0  0  0  0  0  0  0  0
#> 13  0  0  0  0  0  0  0  0  0  0
#> 14  0  0  0  0  0  0  0  0  0  0
#> 15  0  0  0  0  0  0  0  0  0  0
#> 16  0  0  0  0  0  0  0  0  0  0
#> 17  0  0  0  0  0  0  0  0  0  0
#> 18  0  0  0  0  0  0  0  0  0  0
#> 19  0  0  0  0  0  0  0  0  0  0
#> 20  0  0  0  0  0  0  0  0  0  0
#> 21  0  0  0  0  0  0  0  0  0  0
#> 22  0  0  0  0  0  0  0  0  0  0
#> 23  0  0  0  0  0  0  0  0  0  0
#> 24  0  0  0  0  0  0  0  0  0  0
#> 25  0  0  0  0  0  0  0  0  0  0
#> 26  0  0  0  0  0  0  0  0  0  0
#> 27  0  0  0  0  0  0  0  0  0  0
#> 28  0  0  0  0  0  0  0  0  0  0
#> 29  0  0  0  0  0  0  0  0  0  0
#> 30  0  0  0  0  0  0  0  0  0  0
#> 31  0  0  0  0  0  0  0  0  0  0
#> 32  0  0  0  0  0  0  0  0  0  0
#> 33  0  0  0  0  0  0  0  0  0  0
#> 34  0  0  0  0  0  0  0  0  0  0
#> 35  0  0  0  0  0  0  0  0  0  0
#> 36  0  0  0  0  0  0  0  0  0  0
#> 37  0  0  0  0  0  0  0  0  0  0
#> 38  0  0  0  0  0  0  0  0  0  0
#> 39  0  0  0  0  0  0  0  0  0  0
#> 40  0  0  0  0  0  0  0  0  0  0
#> 41  0  0  0  0  0  0  0  0  0  0
#> 42  0  0  0  0  0  0  0  0  0  0
#> 43  0  0  0  0  0  0  0  0  0  0
#> 44  0  0  0  0  0  0  0  0  0  0
#> 45  0  0  0  0  0  0  0  0  0  0
#> 46  0  0  0  0  0  0  0  0  0  0
#> 47  0  0  0  0  0  0  0  0  0  0
#> 48  0  0  0  0  0  0  0  0  0  0
#> 49  0  0  0  0  0  0  0  0  0  0
#> 50  2  1  0  0  0  0  0  0  0  0
#> 51  0  1  0  0  0  0  0  0  0  0
#> 52 NA  0  0  0  0  0  0  0  0  0
#> 53 NA NA  0  2  1  1  0  0  0  0
#> 54 NA NA NA  0  1  3  0  0  0  0
#> 55 NA NA NA NA  0  2  0  0  0  0
#> 56 NA NA NA NA NA  0  0  0  0  0
#> 57 NA NA NA NA NA NA  0  0  0  0
#> 58 NA NA NA NA NA NA NA  0  1  0
#> 59 NA NA NA NA NA NA NA NA  0  0
#> 60 NA NA NA NA NA NA NA NA NA  0
```
