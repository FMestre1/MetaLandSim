# Sample parameter data frame number 1

Sample data frame, as produced by
[`parameter.estimate`](https://fmestre1.github.io/MetaLandSim/reference/parameter.estimate.md).
These parameters are to be passed to
[`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).
These are made up parameters, not related to any species.

## Usage

``` r
data(param1)
```

## Format

A data frame with 4 rows displaying the four parameters (alpha, x, y, e)
to be passed to
[`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md):

- alpha - Parameter relating extinction with distance.

- y - Parameter y in the colonization probability.

- e - Parameter defining the extinction probability in a patch of unit
  area.

- x - Parameter scaling extinction risk with patch area.

## Details

The four parameters are to be passed to
[`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

## Examples

``` r

data(param1)

param1
#>       par_output
#> alpha 0.00100000
#> x     0.50000000
#> y     2.00000000
#> e     0.04662827

#      par_output
#alpha 0.00100000
#x     0.50000000
#y     2.00000000
#e     0.04662827
```
