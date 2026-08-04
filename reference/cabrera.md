# Modified patch occupancy data of Cabrera vole

One season patch occupancy dataset for *Microtus cabrerae* in SW
Portugal. This dataset is in the format produced by
[species.graph](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md),
[convert.graph](https://fmestre1.github.io/MetaLandSim/reference/convert.graph.md)
or
[import.shape](https://fmestre1.github.io/MetaLandSim/reference/import.shape.md)
(class 'metapopulation'), and it was created by converting a data frame
using the function convert.graph. The data frame had the information of
one snapshot of patch occupancy data of Cabrera vole (Microtus cabrera)
in southwestern Portugal.

## Usage

``` r
data(cabrera)
```

## Format

A list with the following elements:

- mapsize - 8200 (landscape mosaic side length, in meters).

- minimum.distance - 10.04 (minimum distance between patches centroids).

- mean.area - 0.46 (mean area, in hectares).

- SD.area - 1.05 (SD of the area).

- number.patches - 793 (number of patches).

- dispersal - 800 (mean dispersal ability of the species).

- distance.to.neighbours - data frame with pairwise distance between
  patches.

- nodes.characteristics - data frame with the characteristics of each
  patch.

## Details

To create this sample dataset the occupancy status of patches was
scrambled, however the proportion of occupied patches was kept.

## Source

Original field data was obtained during project PERSIST
(PTDC/BIA-BEC/105110/2008).

## Examples

``` r
data(cabrera)
```
