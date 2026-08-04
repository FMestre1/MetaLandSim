# Batch landscape simulation

Runs a series of simulations, using
[`iterate.graph`](https://fmestre1.github.io/MetaLandSim/reference/iterate.graph.md),
allows changing the simulations parameters in several sequential
simulations.

## Usage

``` r
manage_landscape_sim(par_df, parameters_spom, full.output)
```

## Arguments

- par_df:

  Arguments data frame to be used by iterate.graph (each row of this
  data frame is a set of Arguments). The data frame has to have the
  following columns in this order (the name of the column is not
  relevant):

  - MDST - Minimum inter-patch distance (in meters).

  - NPATCH - Number of patches in the landscape.

  - AREA_M - Mean area of the patches (in hectares).

  - AREA_SD - SD of the patches' area.

  - MAPSIZE - Landscape mosaic side length (in meters).

  - SPAN - Number of time steps in the simulation.

  - ITER - Number of iterations of the simulation.

  - PAR1_SPAN - parm1 for the span.graph function.

  - PAR2_SPAN - parm2 for the span.graph function.

  - PAR3_SPAN - parm3 for the span.graph function.

  - PAR4_SPAN - parm4 for the span.graph function.

  - PAR5_SPAN - parm5 for the span.graph function.

  - NSEW_SPECIES - Argument nsew for the species.graph function.

  - PARM_SPECIES - Argument parm for the species.graph function.

  - METHOD_SPECIES - Argument method for the species.graph function.

  - KERN - Argument kern for the spom function.

  - CONN - Argument conn for the spom function.

  - COLNZ - Argument colnz for the spom function.

  - EXT - Argument ext for the spom function.

  - BETA1 - Argument beta1 for the spom function.

  - B - Argument b for the spom function.

  - C1 - Argument c1 for the spom function.

  - C2 - Argument c2 for the spom function.

  - Z - Argument z for the spom function.

  - R2 - Argument R for the spom function.

  - DISPERSAL - Species mean dispersal ability (in meters).

  - SUCCESSION - Species successional preference (early, mid or late).

- parameters_spom:

  Parameters data frame, as given by
  [`parameter.estimate`](https://fmestre1.github.io/MetaLandSim/reference/parameter.estimate.md).

- full.output:

  Creates a folder named 'output' to which it saves the full results of
  the simulations made with the parameters in each row of 'par_df'. It
  will generate as many objects as the number of rows in this data
  frame.

## Details

For details regarding the arguments see the respective functions.

## Value

Returns a data frame with the parameters used for the simulations and
the results (mean occupation, mean number of patches, mean turnover,
mean distance and mean area).

## Author

Frederico Mestre and Fernando Canovas

## Note

Depending on computing capacity, this function can take from several
hours to several days to run.

## See also

[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md),
[`species.graph`](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md),
[`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md)

## Examples

``` r

#Setup the parameters for each simulation:
PAR1_SPAN2 <- rep("ncsd",820)#parameter 1 for the span function 
PAR2_SPAN2 <- rep(seq(from=0,to=80,by=2), each=20)#parameter 2 for the span function 
PAR3_SPAN2 <- rep(seq(from=0,to=80,by=2),20)#parameter 3 for the span function 
PAR4_SPAN2 <- rep(2,820)#parameter 4 for the span function 
PAR5_SPAN2 <- rep(2,820)#parameter 5 for the span function 
NSEW_SPECIES2 <- rep("none",820)#where to start populating the landscape 
PARM_SPECIES2 <- rep(5,820)#parameter for the species function 
METHOD_SPECIES2 <- rep("percentage",820)#method for populating the landscape 
MAPSIZE2 <- rep(10000,820)#dimension of the landscape 
SPAN2 <- rep(100,820)#number of time steps of each simulation 
ITER2 <- rep(5,820)#number of iterations of each simulation 
NPATCH2 <- rep(800,820)#number of patches 
AREA_M2 <- rep(0.45,820)#mean area 
AREA_SD2 <- rep(0.2,820)#area sd 
MDST2 <- rep(0,820)#minimum distance between 
KERN <- rep("op1",820)#kernel
CONN <- rep("op1",820)#connectivity function 
COLNZ <- rep("op1",820)#colonization function 
EXT <- rep("op1",820)#extinction function 
BETA1 <- rep("NULL",820) 
B <- rep(1,820) 
C1 <- rep("NULL",820) 
C2 <- rep("NULL",820) 
Z <- rep("NULL",820) 
R2 <- rep("NULL",820) 
DISPERSAL2 <- rep(800,820)#mean dispersal ability of the species
SUCC <- rep("early",820)


#Build parameter data frame (keep the order of the parameters):
simulation <- data.frame(MDST2,NPATCH2,AREA_M2,AREA_SD2,
MAPSIZE2,SPAN2,ITER2,PAR1_SPAN2,PAR2_SPAN2,PAR3_SPAN2,PAR4_SPAN2,PAR5_SPAN2,
NSEW_SPECIES2,PARM_SPECIES2,METHOD_SPECIES2,KERN,CONN,COLNZ,EXT,BETA1,B,C1,C2,Z,R2,DISPERSAL2,SUCC)


#Delete vectors used for data frame creation:
rm('PAR1_SPAN2','PAR2_SPAN2','PAR3_SPAN2','PAR4_SPAN2','PAR5_SPAN2',
'NSEW_SPECIES2','PARM_SPECIES2','METHOD_SPECIES2','MAPSIZE2','SPAN2','ITER2',
'NPATCH2','AREA_M2','AREA_SD2','MDST2','KERN','CONN','COLNZ','EXT',
'BETA1','B','C1','C2','Z','R2','DISPERSAL2','SUCC')


if (FALSE) { # \dontrun{
data(param1)

ms2 <- manage_landscape_sim(par_df=simulation,parameters_spom=param1)
} # }
```
