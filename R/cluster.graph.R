cluster.graph <-
function(rl)
  {
  if(!inherits(rl, "landscape"))
	#if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
    cluster <- as.vector (rl$nodes.characteristics)["cluster"]
    b <- as.data.frame(table(cluster))
    names(b)[names(b)=="Freq"] <- "number of nodes"
    return(b)
  }
