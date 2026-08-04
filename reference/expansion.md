# Class 'expansion'

Class representing an expansion object, as produced by
[`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md).

## Slots

A list of four data frames with the proportion of occupation at several
distances from the closest occupied landscape mosaic. These four data
frames correspond to the proportion of occupation to the north, south,
east and west. Each data frame has the following columns:

- DISTANCE - Distance (mapsize x number of landscapes).

- OCCUPATION - How many times did the landscape at this distance got
  occupied by the species (from a total of 'iter' repetitions).

- PROPORTION - Proportion of occupation for the landscape at this
  distance (OCCUPATION/iter).

## Author

Frederico Mestre and Fernando Canovas
