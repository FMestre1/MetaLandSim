# Create parameter data frame

This function creates a parameter data frame, using parameter values
computed with the application available in the papers of Moilanen (1999)
and ter Braak and Etienne (2003).

## Usage

``` r
create.parameter.df(alpha, x, y, e)
```

## Arguments

- alpha:

  Alpha parameter

- x:

  x parameter

- y:

  y parameter

- e:

  e parameter

## Details

It is highly recommended that the user reads both papers, as well as the
help files.

## Value

Returns a data frame, with the same format as the one returned by
[`parameter.estimate`](https://fmestre1.github.io/MetaLandSim/reference/parameter.estimate.md)
for the methods 'Rsnap_1' and 'Rsnap_x'.

## References

Moilanen, A. (1999). Patch occupancy models of metapopulation dynamics:
efficient parameter estimation using implicit statistical inference.
Ecology, 80(3): 1031-1043.

ter Braak, C. J., & Etienne, R. S. (2003). Improved Bayesian analysis of
metapopulation data with an application to a tree frog metapopulation.
Ecology, 84(1): 231-241.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`parameter.estimate`](https://fmestre1.github.io/MetaLandSim/reference/parameter.estimate.md)

## Examples

``` r

param2 <- create.parameter.df(alpha=0.5, x=0.1, y=5, e=0.1)

param2
#>       par_output
#> alpha        0.5
#> x            0.1
#> y            5.0
#> e            0.1

#      par_output
#alpha        0.5
#x            0.1
#y            5.0
#e            0.1

```
