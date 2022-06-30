extract.graph <-
function (rl,rlist,nr)
  {
  if(!inherits(rl, "landscape") & !inherits(rl, "metapopulation"))
  #if (class(rl)!="landscape" & class(rl)!="metapopulation") 
  {
    stop(paste(rl, " should be either, an object of class class 'landscape' or 'metapopulation'.", sep=""), call. = FALSE)
  }
  if(inherits(rl, "landscape")){
  #if(class(rl)=="landscape"){
        mapsize <- rl$mapsize
    dist_m <- rl$minimum.distance
    disp <- rl$dispersal
    dist_m <-rl$minimum.distance
    rl1 <- rlist[[nr]]
    rl2 <- list(mapsize=mapsize, minimum.distance=dist_m, mean.area=mean(rl1$areas),
                     SD.area=sd(rl1$areas), number.patches=nrow(rl1),
                     dispersal=disp, nodes.characteristics=rl1)
	class(rl2) <- "landscape"
    if(nrow(rl1)>1)rl3 <- cluster.id(rl2)
    if(nrow(rl1)==1)rl3 <- rl2
    return(rl3)
	}
  if(inherits(rl, "metapopulation")){
  #if(class(rl)=="metapopulation"){
	  mapsize <- rl$mapsize
	  dist_m <- rl$minimum.distance
	  disp <- rl$dispersal
	  dist_m <-rl$minimum.distance
	  rl1 <- rlist[[1]][[nr]]
	  
	  neigh <- pairdist(rl1[,1:2])
	  neigh <- as.data.frame(neigh)
	  
	  rownames(neigh) <- 	rl1[,8]
	  colnames(neigh) <- 	rl1[,8]
	  
	  species.out <- list(mapsize=mapsize, minimum.distance=dist_m, 
	                      mean.area=mean(rl1$areas), SD.area=sd(rl1$areas), number.patches=nrow(rl1),
	                      dispersal=disp, distance.to.neighbours=neigh,
	                      nodes.characteristics=rl1)
	  class(species.out) <- "metapopulation"
	  
	  #if(nrow(rl1)>1)rl3 <- cluster.id(rl2)
	  #if(nrow(rl1)==1)rl3 <- rl2
	  return(species.out)
	}
  }
