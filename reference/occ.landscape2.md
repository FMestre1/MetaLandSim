# Sample landscape with 10 simulated occupancy snapshots

Sample species occupancy in a network during 10 time steps. Simulated
data.

## Usage

``` r
data(occ.landscape2)
```

## Format

A list with the following elements:

- mapsize - landscape mosaic side length, in meters.

- minimum.distance - minimum distance between patches centroids.

- mean.area - mean area, in hectares.

- SD.area - standard deviation of the area.

- number.patches - number of patches.

- dispersal - mean dispersal ability of the species.

- distance.to.neighbours - data frame with pairwise distance between
  patches.

- nodes.characteristics - data frame with the characteristics of each
  patch, (species 1 to 10 - occupancy snapshots).

## Examples

``` r
data(occ.landscape2)
```
