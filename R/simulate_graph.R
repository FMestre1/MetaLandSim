simulate_graph <-
function (rl, rlist, simulate.start, method, parm,
                            nsew="none", param_df,kern, conn, 
							colnz, ext,beta1,b, c1, c2, z, R)
  {
    span <- length(rlist)
    metpop.list <- as.list(rep("", span))
    turnover.list <- as.list(rep("", span))
    if(simulate.start==TRUE)
      {
        sp_0 <- species.graph(rl,method=method,parm=parm,nsew=nsew,plotG=FALSE)
      }
    if (simulate.start==FALSE)
      {
        sp_0 <- rl
      }
    sp_1 <- sp_0$nodes.characteristics
    metpop.list[[1]] <- sp_1
    turnover.list[[1]] <- 0
    for(i in 2:span)
      {
        prec.sp <- metpop.list[[i-1]]
        mapsize <- as.numeric(rl[[1]])
        minimum.distance <- as.numeric(rl[[2]])
        mean.area <- mean(prec.sp$areas)
        SD.area <- sd(prec.sp$areas)
        number.patches <- nrow(prec.sp)
        dispersal <- as.numeric(rl[[6]])
        neigh <- as.data.frame(pairdist(prec.sp[, 1:2]))
        names(neigh) <- prec.sp$ID
        rownames(neigh) <- prec.sp$ID
        prec.sp$nneighbour <- nndist(prec.sp[, 1:2])
        prec.sp_1 <- list(mapsize=mapsize, minimum.distance=minimum.distance,
                          mean.area=mean.area, SD.area=SD.area,
                          number.patches=number.patches,dispersal=dispersal,
                          distance.to.neighbours=neigh,nodes.characteristics=prec.sp)
		class(prec.sp_1) <- "metapopulation"
        out_0 <- spom(prec.sp_1, kern, conn, colnz, ext, param_df, beta1,
                      b, c1, c2, z, R)
        turnover.list[[i]] <- ((out_0$turnover*100)/nrow(out_0$nodes.characteristics))
        out_1 <- out_0$nodes.characteristics[, -c(9,11)]
        names(out_1)[names(out_1)=="species2"] <- "species"
        lands_i <- rlist[[i]]
        out_2 <- merge_order(lands_i, out_1, by.x = "ID", by.y = "ID",sort=FALSE,keep_order=TRUE,all.x=TRUE,all.y=TRUE)
        if(any(is.na(out_2[, 2:8]))==TRUE)
          {
            out_3 <- out_2[, c(1:8,16)]
            out_3 <- na.omit(out_3)
          }
		else out_3 <- out_2
		if(any(is.na(out_2[, 9:16]))==TRUE)
          {
		    out_4 <- out_2[is.na(out_2$species),]
            out_4 <- out_4[,-c(9:15)]
		    out_4[,9] <- 0 
		    out_3 <- rbind(out_3,out_4)
          }
        out_4 <- data.frame(out_3$x.x, out_3$y.x, out_3$areas.x, out_3$radius.x,
                             out_3$cluster.x, out_3$colour.x, out_3$nneighbour.x,
                             out_3$ID, out_3$species)
        names(out_4)[names(out_4)=="out_3.x.x"] <- "x"
        names(out_4)[names(out_4)=="out_3.y.x"] <- "y"
        names(out_4)[names(out_4)=="out_3.areas.x"] <- "areas"
        names(out_4)[names(out_4)=="out_3.radius.x"] <- "radius"
        names(out_4)[names(out_4)=="out_3.cluster.x"] <- "cluster"
        names(out_4)[names(out_4)=="out_3.colour.x"] <- "colour"
        names(out_4)[names(out_4)=="out_3.nneighbour.x"] <- "nneighbour"
        names(out_4)[names(out_4)=="out_3.ID"] <- "ID"
        names(out_4)[names(out_4)=="out_3.species"] <- "species"
    
        metpop.list[[i]] <- out_4
      }
    turnover <-as.numeric(turnover.list)
    output <- list(metpop.list,turnover=turnover)
    return(output)
  }
