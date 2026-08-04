# Sample parameter data frame number 2

Sample data frame, as produced by
[`parameter.estimate`](https://fmestre1.github.io/MetaLandSim/reference/parameter.estimate.md).
These parameters are to be passed to
[`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).
These are based upon those computed for Cabrerae's vole in the paper
Mestre et al. (2017) (see references).

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

## References

Mestre, F., Risk, B., Mira, A., Beja, P., Pita, R. (2017)
\<doi:10.1016/j.ecolmodel.2017.06.013\>

## Details

The four parameters are to be passed to
[`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

## Examples

``` r

data(param2)

param2
#>       par_output
#> alpha    0.00047
#> x        0.44000
#> y       18.15000
#> e        0.00480

#      par_output
#alpha    0.00047
#x        0.44000
#y       18.15000
#e        0.00480

```
