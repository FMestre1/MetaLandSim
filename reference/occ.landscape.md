# Sample landscape with one simulated occupancy snapshot

Sample random landscape graph, with species occupancy data (occupancy
rate - 50%). Simulated data.

## Usage

``` r
data(occ.landscape)
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
  patch.

## Examples

``` r

data(occ.landscape)
```
