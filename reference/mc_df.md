# Modified patch occupancy data of Cabrera vole as a data frame

One season patch occupancy dataset for *Microtus cabrerae* in SW
Portugal (modified). This dataset is in a format directly used by
[convert.graph](https://fmestre1.github.io/MetaLandSim/reference/convert.graph.md)
and converted to an object class 'metapopulation'.

## Usage

``` r
data(mc_df)
```

## Format

A data frame with 685 observations on the following 5 variables.

- `ID`:

  Patch Id.

- `x`:

  X coordinate.

- `y`:

  Y coordinate.

- `area`:

  Patch area, in hectares.

- `mc`:

  Occupancy state (0/1).

## Details

To create this sample dataset the occupancy status of patches was
scrambled, however the proportion of occupied patches was kept.

## Source

Original field data was obtained during project PERSIST
(PTDC/BIA-BEC/105110/2008).

## Examples

``` r

##To be converted in a object of class "metapopulation":
#mc1 <- convert.graph(dframe=mc_df,mapsize=8200,dispersal=800)

data(mc_df)

#Check the columns:

head(mc_df)
#>   ID        x       y  area mc
#> 1  1 1248.254   0.000 0.079  0
#> 2  2 1420.857  46.725 0.781  1
#> 3  3 1278.912  52.629 1.053  0
#> 4  4 6370.625  62.637 0.788  0
#> 5  5 1151.337  97.140 0.079  0
#> 6  6 1295.796 104.839 0.137  1

#  ID        x       y  area mc
#1  1 1248.254   0.000 0.079  0
#2  2 1420.857  46.725 0.781  1
#3  3 1278.912  52.629 1.053  1
#4  4 6370.625  62.637 0.788  0
#5  5 1151.337  97.140 0.079  0
#6  6 1295.796 104.839 0.137  1

```
