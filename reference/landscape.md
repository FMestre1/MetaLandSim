# Class 'landscape'

Class representing a landscape graph, as produced by
[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`convert.graph`](https://fmestre1.github.io/MetaLandSim/reference/convert.graph.md)
and
[`import.shape`](https://fmestre1.github.io/MetaLandSim/reference/import.shape.md).

## Slots

- mapsize - Side of the landscape in meters.

- minimum.distance - Minimum distance between patches centroids, in
  meters.

- mean.area - Mean patch area in hectares.

- SD.area - Standard deviation of patches area.

- number.patches - Total number of patches.

- dispersal - Species mean dispersal ability, in meters.

- nodes.characteristics - Data frame with patch (node) information
  (coordinates, area, radius, cluster, distance to nearest neighbor and
  ID).

## Author

Frederico Mestre and Fernando Canovas
