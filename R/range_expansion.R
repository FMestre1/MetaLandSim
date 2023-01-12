range_expansion <-
  function(rl, percI, param, b=1, tsteps, iter, plot = FALSE)
  {
    if(!inherits(rl, "landscape"))  
    {
      stop(paste(rl, " should be an object of class class 'landscape'.", sep = ""), 
           call. = FALSE)
    }
    
    mapsize <- rl$mapsize
    dist_m <- rl$minimum.distance
    areaM <- rl$mean.area
    areaSD <- rl$SD.area
    Npatch <- rl$number.patches
    disp <- rl$dispersal
      
    node.expansion <- function(occ_landscape, param, b, tsteps) {
      output0 <- c()
      npatch <- occ_landscape$number.patches
      mapsize <- occ_landscape$mapsize
      ID_land <- max(occ_landscape$nodes.characteristics$ID)
      nrow_land <- nrow(occ_landscape$nodes.characteristics)
      areaM2 <- occ_landscape$mean.area
      areaSD2 <- occ_landscape$SD.area
      dispersal <- 1/param[1,1]
      ocupp <- "N"
      for(j in 1:tsteps) {
        
        ocupp <- "N"
        
        occ_landscape_new <- spom(sp = occ_landscape, kern = "op1", conn = "op1", colnz = "op1", ext = "op1", param_df = param, b, 
                                  c1 = NULL, c2 = NULL, z = NULL, R = NULL, succ="none")
        
        if(sum(occ_landscape_new$nodes.characteristics$species2)==0) {
          message(paste("Empty landscape. Next simulation... ","- time step ", j,sep=""))
          break		
        }
        
        perc_occup <- (sum(occ_landscape_new$nodes.characteristics$species2[1:nrow_land])*100)/nrow_land
        v0 <- occ_landscape_new$nodes.characteristics[ which(occ_landscape_new$nodes.characteristics$x < dispersal/2), ]$species2
        
        if(sum(v0)==0) ocupp <- "N"
        if(sum(v0)!=0) ocupp <- "Y"
        
        if(ocupp == "N") {
          message(paste("No transition between landscape units. ","- time step ", j,sep=""))
          occ_landscape <- occ_landscape_new
          occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[, -c(9, 11)]
          names(occ_landscape$nodes.characteristics)[names(occ_landscape$nodes.characteristics) == "species2"] <- "species"
          occ_landscape <- occ_landscape[-9]
          occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[-(nrow_land + 1), ]
          occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[-(nrow_land + 1), ]
          occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[, -(nrow_land + 1)]
          occ_landscape$number.patches <- occ_landscape$number.patches - 1
          class(occ_landscape) <- "metapopulation"
        }
        if(ocupp == "Y"){
          message(paste("Transition between landscapes. New landscape created. ","- time step ", j,sep=""))
          
          v2 <- length(occ_landscape$nodes.characteristics[ which(occ_landscape$nodes.characteristics$y < dispersal/2), ]$ID)
          
          number_patches <- v2*(nrow_land/100)
          
          rl1 <- rland.graph(mapsize = occ_landscape$mapsize, dist_m = occ_landscape$minimum.distance, 
                             areaM2, areaSD2, Npatch = npatch, disp = occ_landscape$dispersal, plotG = FALSE)

          occ_landscape <- suppressWarnings(species.graph(rl = rl1, method = "number", parm = number_patches, 
                            plotG = FALSE))
          
          message(paste("Saving in output file! ","- time step ", j,sep=""))
          output0 <- c(output0,j)
          
        }
        
      }
      
      if(sum(output0)==0) output0 <- 0
      output1 <- cbind(mapsize,output0)
      output1 <- cbind(cumsum(output1[,1]),output1[,2])
      colnames(output1) <- c("DISTANCE", "TIME STEP")
      output1 <- as.data.frame(output1)
      return(output1)
      
    }#END
    
    distance <- mapsize * 1:tsteps
    output_all <- as.data.frame(distance)
    
    for (i in 1:iter) {
      sp1 <- species.graph(rl = rl, method = "percentage", parm = percI, nsew = "none", plotG = FALSE)
      #
      nodeALL <- node.expansion(occ_landscape = sp1, param, b, tsteps)
      output_all <- suppressWarnings(cbind(output_all, c(nodeALL[, 2], rep(NA, nrow(output_all) - length(nodeALL[, 2])))))
    }
    
    output_all[is.na(output_all)] <- 0
    
    tstep_average_all <- rowMeans(output_all[,2:(iter+1)], na.rm=FALSE)
    
    for (x in 2:(iter + 1)) {
      i_all <- which(output_all[, x] != 0)
      output_all[i_all, x] <- 1
    }
    
    output_all <- cbind(output_all[, 1], rowSums(as.data.frame(output_all[, 2:(ncol(output_all))])))
    output_all <- cbind(output_all[, 1:2], output_all[, 2]/iter, tstep_average_all)
    output_all[is.na(output_all)] <- 0
    
    output_all <- as.data.frame(output_all)
    names(output_all) <- c("DISTANCE", "OCCUPATION", "PROPORTION", "TIME STEP")
    
    p1all <- gvisLineChart(output_all, 
                           xvar = "DISTANCE", 
                           yvar = "PROPORTION", 
                           options = list(title = "Dispersal Simulation", 
                           width = 500, 
                           height = 300, 
                           curveType = "function", 
                           legend = "none", 
                           titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
                           vAxis = "{title: 'proportion'}", 
                           hAxis = "{title: 'distance (m)'}", 
                           series = "[{color: '#006400'}]", 
                           backgroundColor = "#D1EEEE")
                           )

    p1all$html$caption <- paste("<div><span>Range expansion , considering that ", 
                                   percI, "% of the patches are occupied in the first landscape mosaic.</span><br />", 
                                   sep = "")
    
    p1all$html$footer <- paste("\n<!-- htmlFooter -->\n<span> \n  ",
                                  R.Version()$version.string,
                                  " &#8226; <a href=\"http://code.google.com/p/google-motion-charts-with-r/\">googleVis-", 
                                  packageVersion("googleVis"),
                                  "</a>\n  &#8226; MetaLandSim-",
                                  packageVersion("MetaLandSim"),
                                  "\n  &#8226; <a href=\"https://developers.google.com/terms/\">Google Terms of Use</a> &#8226; <a href=\"https://google-developers.appspot.com/chart/interactive/docs/gallery/linechart.html#Data_Policy\">Data Policy</a>\n</span></div>\n</body>\n</html>\n", 
                                  sep="")
    
    if(plot == TRUE) plot(p1all)
    #class(output_all) <- "expansion"
    return(output_all)
  }
