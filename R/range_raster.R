range_raster <-
function(presences.map, re.out, mask.map=NULL, plot.directions=TRUE)
  {
    if(class(re.out) != "expansion") stop(paste(re.out, " should be an object of class class 'expansion'", sep=""), call.=FALSE)
    execGRASS("r.mask", flags="r")
    if(!is.null(mask.map))
      {
        execGRASS("r.in.gdal", input=mask.map, output="mask", flags=c("overwrite", "o"))
        execGRASS("g.region", rast = "mask")
        execGRASS("r.null", map="mask", setnull="0")
        execGRASS("r.in.gdal", input=presences.map, output="presences", flags=c("overwrite", "o"))
      }
    else {
        execGRASS("r.in.gdal", input=presences.map, output="presences", flags=c("overwrite", "o", "e"))
        execGRASS("g.region", rast = "presences")
    }
    execGRASS("r.null", map="presences", setnull="0")
    wrapping <- function(x) strsplit(x, " ")
    extract.value <- function(parameters, cardinal="north")
      {
        agrep(parameters, pattern=cardinal, value=TRUE, max.distance=list(all=0)) -> temp
        wrapping(temp)[[1]] -> result
        as.numeric(result[length(result)]) -> result
        return(result)
      }
    create.vector <- function(in.name="north_horizon", out.name="north_horizon", type="line")
      {
        if(type == "line")
          {
            execGRASS("r.thin", input=in.name, output="tmp2", flags="overwrite")
            flag.vector <- c("overwrite")
            execGRASS("r.to.vect", input="tmp2", output=out.name, feature=type, flags=flag.vector)
            execGRASS("g.remove", rast="tmp2")
          } 
        if(type == "point")
          {
            flag.vector <- c("overwrite","z")
            execGRASS("r.to.vect", input=in.name, output=out.name, feature=type, flags=flag.vector)
          }
      }
    fit.sigmoid <- function(y, x, start.params = list(a = 1, b = 0.5, c = 50))
      {
        fitmodel <- nls(y ~ a / (1 + exp(b * (x - c))), start = start.params)
        return(coef(fitmodel))
      }
    predict.sigmoid <- function(params, x)
      {
        return(params[1] / (1 + exp(params[2] * (x - params[3]))))
      }
    execGRASS("g.region", flags="p",intern=T) -> region.parameters
    extract.value(region.parameters, cardinal="north") -> north
    extract.value(region.parameters, cardinal="south") -> south
    extract.value(region.parameters, cardinal="east") -> east
    extract.value(region.parameters, cardinal="west") -> west
    extract.value(region.parameters, cardinal="nsres") -> nsres
    extract.value(region.parameters, cardinal="ewres") -> ewres
    execGRASS("r.mapcalculator", outfile="north_horizon", formula=paste("y() > ",
              north-nsres, sep=""), flags="overwrite")
    execGRASS("r.null", map="north_horizon", setnull="0")
    create.vector(in.name="north_horizon", out.name="north_horizon")
    execGRASS("r.mapcalculator", outfile="south_horizon", formula=paste("y() < ",
              south+nsres, sep=""), flags="overwrite")
    execGRASS("r.null", map="south_horizon", setnull="0")
    create.vector(in.name="south_horizon", out.name="south_horizon")
    execGRASS("r.mapcalculator", outfile="west_horizon", formula=paste("x() < ",
              west+ewres, sep=""), flags="overwrite")
    execGRASS("r.null", map="west_horizon", setnull="0")
    create.vector(in.name="west_horizon", out.name="west_horizon")
    execGRASS("r.mapcalculator", outfile="east_horizon", formula=paste("x() > ",
              east-ewres, sep=""), flags="overwrite")
    execGRASS("r.null", map="east_horizon", setnull="0")
    create.vector(in.name="east_horizon", out.name="east_horizon")
    execGRASS("r.mapcalculator", outfile="wd", formula="1", flags="overwrite")
    create.vector(in.name="presences", out.name="presences", type="point")
    execGRASS("r.cost", input="wd", output="temptrend", start_rast="north_horizon", flags="overwrite")
    execGRASS("r.cost", input="wd", output="tempdir", start_rast="presences", stop_points="north_horizon", flags="overwrite")
    execGRASS("r.mapcalculator", amap="temptrend", bmap="tempdir", outfile="tempout",
              formula=paste("(A + B)", sep=""), flags="overwrite")
    execGRASS("r.univar", map="tempout", intern=T) -> map.stats
    extract.value(map.stats,"minimum") -> minimum.map
    extract.value(map.stats,"maximum") -> maximum.map
    execGRASS("r.mapcalculator", amap="tempout", outfile="Ndirectionality",
              formula=paste("100 - ((A - ", minimum.map, ") * 100 / ", maximum.map-minimum.map, ")", sep=""),
              flags="overwrite")
    execGRASS("r.cost", input="wd", output="temptrend", start_rast="south_horizon", flags="overwrite")
    execGRASS("r.cost", input="wd", output="tempdir", start_rast="presences", stop_points="south_horizon", flags="overwrite")
    execGRASS("r.mapcalculator", amap="temptrend", bmap="tempdir", outfile="tempout",
              formula=paste("(A + B)", sep=""), flags="overwrite")
    execGRASS("r.univar", map="tempout", intern=T) -> map.stats
    extract.value(map.stats,"minimum") -> minimum.map
    extract.value(map.stats,"maximum") -> maximum.map
    execGRASS("r.mapcalculator", amap="tempout", outfile="Sdirectionality",
              formula=paste("100 - ((A - ", minimum.map, ") * 100 / ", maximum.map-minimum.map, ")", sep=""),
              flags="overwrite")
    execGRASS("r.cost", input="wd", output="temptrend", start_rast="east_horizon", flags="overwrite")
    execGRASS("r.cost", input="wd", output="tempdir", start_rast="presences", stop_points="east_horizon", flags="overwrite")
    execGRASS("r.mapcalculator", amap="temptrend", bmap="tempdir", outfile="tempout",
              formula=paste("(A + B)", sep=""), flags="overwrite")
    execGRASS("r.univar", map="tempout", intern=T) -> map.stats
    extract.value(map.stats,"minimum") -> minimum.map
    extract.value(map.stats,"maximum") -> maximum.map
    execGRASS("r.mapcalculator", amap="tempout", outfile="Edirectionality",
              formula=paste("100 - ((A - ", minimum.map, ") * 100 / ", maximum.map-minimum.map, ")", sep=""),
              flags="overwrite")
    execGRASS("r.cost", input="wd", output="temptrend", start_rast="west_horizon", flags="overwrite")
    execGRASS("r.cost", input="wd", output="tempdir", start_rast="presences", stop_points="west_horizon", flags="overwrite")
    execGRASS("r.mapcalculator", amap="temptrend", bmap="tempdir", outfile="tempout",
              formula=paste("(A + B)", sep=""), flags="overwrite")
    execGRASS("r.univar", map="tempout", intern=T) -> map.stats
    extract.value(map.stats,"minimum") -> minimum.map
    extract.value(map.stats,"maximum") -> maximum.map
    execGRASS("r.mapcalculator", amap="tempout", outfile="Wdirectionality",
              formula=paste("100 - ((A - ", minimum.map, ") * 100 / ", maximum.map-minimum.map, ")", sep=""),
              flags="overwrite")
    execGRASS("g.remove", rast="temptrend,tempdir,tempout")
    params.n <- fit.sigmoid(re.out$NORTH$PROPORTION, re.out$NORTH$DISTANCE/nsres)
    print(params.n)
    params.s <- fit.sigmoid(re.out$SOUTH$PROPORTION, re.out$SOUTH$DISTANCE/nsres)
    print(params.s)
    params.e <- fit.sigmoid(re.out$EAST$PROPORTION, re.out$EAST$DISTANCE/ewres)
    print(params.e)
    params.w <- fit.sigmoid(re.out$WEST$PROPORTION, re.out$WEST$DISTANCE/ewres)
    print(params.w)
    if(!is.null(mask.map))
      {
        execGRASS("r.mask", input="mask")
      }
    execGRASS("r.cost", input="wd", output="distances", start_rast="presences", flags="overwrite")
    execGRASS("r.mapcalculator", amap="distances", bmap="Ndirectionality", outfile="Nprobability",
               formula=paste("B * (", params.n[1]," / (1 + exp(", params.n[2]," * (A - ", params.n[3],"))))", sep=""), flags="overwrite")
    execGRASS("r.mapcalculator", amap="distances", bmap="Sdirectionality", outfile="Sprobability",
               formula=paste("B * (", params.s[1]," / (1 + exp(", params.s[2]," * (A - ", params.s[3],"))))", sep=""), flags="overwrite")
    execGRASS("r.mapcalculator", amap="distances", bmap="Edirectionality", outfile="Eprobability",
               formula=paste("B * (", params.e[1]," / (1 + exp(", params.e[2]," * (A - ", params.e[3],"))))", sep=""), flags="overwrite")
    execGRASS("r.mapcalculator", amap="distances", bmap="Wdirectionality", outfile="Wprobability",
               formula=paste("B * (", params.w[1]," / (1 + exp(", params.w[2]," * (A - ", params.w[3],"))))", sep=""), flags="overwrite")
    execGRASS("r.mapcalculator", amap="Nprobability", bmap="Sprobability", cmap="Eprobability", dmap="Wprobability",
              outfile="range",
              formula="(A + B +C + D) / 4", flags="overwrite")
    output <- raster(readRAST6("range"))
    if(plot.directions == TRUE)
      {
        dev.new()
        par(mfrow=c(1,2))
        plot(re.out$NORTH$DISTANCE, predict.sigmoid(params.n, re.out$NORTH$DISTANCE/nsres),type="l",
             main="Northern direction", xlab="Distance (m)", ylab="Probability")
        plot(raster(readRAST6("Nprobability")),main="Northern probability")
        contour(raster(readRAST6("Nprobability")), add=T)
        dev.new()
        par(mfrow=c(1,2))
        plot(re.out$SOUTH$DISTANCE, predict.sigmoid(params.s, re.out$SOUTH$DISTANCE/nsres),type="l",
             main="Southern direction", xlab="Distance (m)", ylab="Probability")
        plot(raster(readRAST6("Sprobability")),main="Southern probability")
        contour(raster(readRAST6("Sprobability")), add=T)
        dev.new()
        par(mfrow=c(1,2))
        plot(re.out$EAST$DISTANCE, predict.sigmoid(params.e, re.out$EAST$DISTANCE/ewres),type="l",
             main="Eastern direction", xlab="Distance (m)", ylab="Probability")
        plot(raster(readRAST6("Eprobability")),main="Eastern probability")
        contour(raster(readRAST6("Eprobability")), add=T)
        dev.new()
        par(mfrow=c(1,2))
        plot(re.out$WEST$DISTANCE, predict.sigmoid(params.w, re.out$WEST$DISTANCE/ewres),type="l",
             main="Western direction", xlab="Distance (m)", ylab="Probability")
        plot(raster(readRAST6("Wprobability")),main="Western probability")
        contour(raster(readRAST6("Wprobability")), add=T)
      }
    vect.list <- execGRASS("g.mlist", type="vect", separator=",", intern=TRUE)
    rast.list <- execGRASS("g.mlist", type="rast", separator=",", intern=TRUE)
    if(!is.null(mask.map))
      {
        execGRASS("r.mask", flags="r")
      }
    execGRASS("g.remove", rast=rast.list, vect=vect.list)
    return(output)
    }
