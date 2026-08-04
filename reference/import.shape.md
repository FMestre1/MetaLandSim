# Import a shapefile

Imports a shapefile, converting it to an object of class
'metapopulation' or 'landscape'.

## Usage

``` r
import.shape(filename, path, species.col, ID.col, area.col, dispersal, 
class.landscape=FALSE)
```

## Arguments

- filename:

  Character vector with the shapefile name.

- path:

  Character vector with the path to the file.

- species.col:

  Character vector with the name of the column (in the shapefile) with
  the species occupancy data.

- ID.col:

  Character vector with the name of the column (in the shapefile) with
  the patch Id.

- area.col:

  Character vector with the name of the column (in the shapefile) with
  the patch area, in hectares.

- dispersal:

  Species mean dispersal ability, in meters.

- class.landscape:

  Should the output belong to the class 'metapopulation' or 'landscape'.

## Value

Delivers an object of class 'metapopulation' or 'landscape'.

## Author

Frederico Mestre and Fernando Canovas

## Note

The shapefile must be in project coordinates (units=meters and
hectares).

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`convert.graph`](https://fmestre1.github.io/MetaLandSim/reference/convert.graph.md)

## Examples

``` r

if (FALSE) { # \dontrun{

rl1 <- import.shape(filename = "yourshapefile.shp"
      ,path = "C:/yourpath..."
      ,species.col= "column with species"
      ,ID.col="column with patch Id"
      ,area.col="Column with area"
      ,dispersal=800#Mean dispersal ability of the species 
      #(used to generate patch clusters, or components)
      )

} # }
```
