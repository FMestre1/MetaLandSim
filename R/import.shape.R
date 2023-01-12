import.shape <-
function(filename, path=NULL, species.col, ID.col, area.col, dispersal, class.landscape=FALSE)
  {
    pathfile <- paste(path, filename,sep="")
    #sf <- readShapePoly(pathfile)
    sf <- terra::vect(pathfile)#NEW
    #df1 <- sf@data
    df1 <- as.data.frame(sf)#NEW
    #ctr <- gCentroid(sf, byid=TRUE)
    ctr <- terra::centroids(sf)
    ctr<-as.data.frame(terra::crds(ctr))#NEW
    ID <- df1[, ID.col]
    area <- df1[, area.col]
    if(is.character(species.col)) species <- df1[, species.col]
    if(is.character(species.col)) df3 <- cbind(ID, ctr, area, species)
    else df3 <- cbind(ID,ctr,area)
    mapsize <- max(c(max(ctr$x)-min(ctr$x),max(ctr$y)-min(ctr$y)))
    if (class.landscape==FALSE) return(convert.graph(dframe=df3,mapsize,dispersal))
    if (class.landscape==TRUE) return(remove.species(convert.graph(dframe=df3,mapsize,dispersal)))
  }
