### R code from vignette source 'range_expansion.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: range_expansion.Rnw:27-42
###################################################
library(MetaLandSim)

#Create starting landscape (the simulation will assume that 
#all subsequent landscapes are built with the same parameter combination).

rl1 <- rland.graph(mapsize=10000, dist_m=10, areaM=0.05, areaSD=0.02, 
Npatch=250, disp=800, plotG=TRUE)

#Create range expansion model. Here run only with two repetitions (iter=2). 
#Ideally it should be run with more repetitions to provide more robust results.

#Not run
#rg_exp1 <- range_expansion(rl=rl1, percI=50, amin=0, param=param1,
#b=1, tsteps=100, iter=2)
#End (Not run)


###################################################
### code chunk number 2: range_expansion.Rnw:58-73
###################################################
data(rg_exp)
presences <- paste(system.file(package="MetaLandSim"),
 "/examples/presences.asc", sep="")
landmask <- paste(system.file(package="MetaLandSim"), 
"/examples/landmask.asc", sep="")

#require(rgrass7)

#First, start GRASS from R: 
#initGRASS(gisBase = "grass folder", home = tempdir(), 
#gisDbase = "mapset location",override = TRUE)
#Create raster, using the sample dataset 
#rg_exp (generated with 100 repetitions)
#range_raster(presences.map = presences, re.out=rg_exp, 
#mask.map=landmask, plot.directions=FALSE)


