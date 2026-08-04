# Class 'metapopulation'

Class representing a landscape graph with species' patch occupancy data,
as produced by
[`species.graph`](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md),
[`convert.graph`](https://fmestre1.github.io/MetaLandSim/reference/convert.graph.md)
and
[`import.shape`](https://fmestre1.github.io/MetaLandSim/reference/import.shape.md).

## Slots

- mapsize - Landscape mosaic side length, in meters.

- minimum.distance - Minimum distance between patches centroids, in
  meters.

- mean.area - Mean patch area in hectares.

- SD.area - Standard deviation of patches area.

- number.patches - Total number of patches.

- dispersal - Species mean dispersal ability, in meters.

- distance.to.neighbours - Data frame with pairwise distance between
  patches, in meters.

- nodes.characteristics - Data frame with patch (node) information
  (coordinates, area, radius, cluster, distance to nearest neighbor, ID
  and species).

## Author

Frederico Mestre and Fernando Canovas
