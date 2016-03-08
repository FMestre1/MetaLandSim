range_expansion <-
function(rl, percI, param, b=1, tsteps, iter)
  {
    if (class(rl) != "landscape")
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
    node.expansion <- function(occ_landscape, param, b, node, tsteps) {
        output0 <- as.data.frame(matrix(nrow = tsteps, ncol = 2))  
        output0[, 1] <- 1:nrow(output0)
        npatch <- occ_landscape$number.patches
        mapsize <- occ_landscape$mapsize
        ID_land <- max(occ_landscape$nodes.characteristics$ID)
        mapsize <- occ_landscape$mapsize  #recuperar o mapsize
        nrow_land <- nrow(occ_landscape$nodes.characteristics)
        clock <- 1
        if (node == "North") {
            xx1xx <- "S"
            xx2xx <- NA
            xx3xx <- mapsize
            xx4xx <- "North"
            xx7xx <- NA  
            xx8xx <- 0  
            xx9xx <- "South"  
        }
        if (node == "South") {
            xx1xx <- "N"
            xx2xx <- NA
            xx3xx <- 0
            xx4xx <- "South"
            xx7xx <- NA  
            xx8xx <- mapsize  
            xx9xx <- "North"  
        }
        if (node == "East") {
            xx1xx <- "W"
            xx2xx <- mapsize
            xx3xx <- NA
            xx4xx <- "East"
            xx7xx <- 0  
            xx8xx <- NA  
            xx9xx <- "West"  
        }
        if (node == "West") {
            xx1xx <- "E"
            xx2xx <- 0
            xx3xx <- NA
            xx4xx <- "West"
            xx7xx <- mapsize  
            xx8xx <- NA  
            xx9xx <- "East"  
        }
        repeat {
            if (node == "North") {
                xx5xx <- min(abs(mapsize - occ_landscape$nodes.characteristics[, 
                  2])[abs(mapsize - occ_landscape$nodes.characteristics[, 2]) != 
                  0])
                xx6xx <- abs(mapsize - occ_landscape$nodes.characteristics[, 2])
                xx10xx <- min(occ_landscape$nodes.characteristics[, 2][occ_landscape$nodes.characteristics[, 
                  2] != 0])  
                xx11xx <- occ_landscape$nodes.characteristics[, 2]  
            }
            if (node == "South") {
                xx5xx <- min(occ_landscape$nodes.characteristics[, 2][occ_landscape$nodes.characteristics[, 
                  2] != 0])
                xx6xx <- occ_landscape$nodes.characteristics[, 2]
                xx10xx <- min(abs(mapsize - occ_landscape$nodes.characteristics[, 
                  2])[abs(mapsize - occ_landscape$nodes.characteristics[, 2]) != 
                  0])  
                xx11xx <- abs(mapsize - occ_landscape$nodes.characteristics[, 2])  
            }
            if (node == "East") {
                xx5xx <- min(abs(mapsize - occ_landscape$nodes.characteristics[, 
                  1])[abs(mapsize - occ_landscape$nodes.characteristics[, 1]) != 
                  0])
                xx6xx <- abs(mapsize - occ_landscape$nodes.characteristics[, 1])
                xx10xx <- min(occ_landscape$nodes.characteristics[, 1][occ_landscape$nodes.characteristics[, 
                  1] != 0])  
                xx11xx <- occ_landscape$nodes.characteristics[, 1]  
            }
            if (node == "West") {
                xx5xx <- min(occ_landscape$nodes.characteristics[, 1][occ_landscape$nodes.characteristics[, 
                  1] != 0])
                xx6xx <- occ_landscape$nodes.characteristics[, 1]
                xx10xx <- min(abs(mapsize - occ_landscape$nodes.characteristics[, 
                  1])[abs(mapsize - occ_landscape$nodes.characteristics[, 1]) != 
                  0])  
                xx11xx <- abs(mapsize - occ_landscape$nodes.characteristics[, 1])  
            }
            ocupp <- "N"
            if (clock == 1 || (clock >= 2 & ocupp == "N")) {
                occ_landscape$nodes.characteristics[nrow_land + 1, ] <- NA
                occ_landscape$nodes.characteristics$ID[(nrow_land + 1)] <- (ID_land + 
                  1)
                area1 <- mean(occ_landscape$nodes.characteristics$areas[1:nrow_land])
                occ_landscape$nodes.characteristics$radius[(nrow_land + 1)] <- sqrt((area1 * 
                  10000)/pi)
                occ_landscape$nodes.characteristics$areas[(nrow_land + 1)] <- (occ_landscape$nodes.characteristics$radius[(nrow_land + 
                  1)]) * mapsize
                occ_landscape$nodes.characteristics$y[nrow_land + 1] <- xx3xx
                occ_landscape$nodes.characteristics$x[nrow_land + 1] <- xx2xx
                max_cluster <- max(occ_landscape$nodes.characteristics$cluster, na.rm = T)
                occ_landscape$nodes.characteristics$cluster[(nrow_land + 1)] <- c(max_cluster + 
                  1)
                occ_landscape$nodes.characteristics$colour <- as.character(occ_landscape$nodes.characteristics$colour)
                occ_landscape$nodes.characteristics$colour[(nrow_land + 1)] <- xx4xx
                occ_landscape$nodes.characteristics$species[(nrow_land + 1)] <- 0
                occ_landscape$nodes.characteristics$nneighbour[(nrow_land + 1)] <- xx5xx
                occ_landscape$number.patches <- occ_landscape$number.patches + 1
                dist_nodos <- occ_landscape$distance.to.neighbours
                dist_nodos[, (nrow_land + 1)] <- c(xx6xx)
                colnames(dist_nodos)[(nrow_land + 1)] <- xx4xx
                dist_nodos[(nrow_land + 1), ] <- c(xx6xx, 0)
                rownames(dist_nodos)[(nrow_land + 1)] <- xx4xx
                occ_landscape$distance.to.neighbours <- dist_nodos
            }
            if (clock >= 2 & ocupp == "Y") {
                rl1 <- rland.graph(mapsize = occ_landscape$mapsize, dist_m = occ_landscape$minimum.distance, 
                  areaM = occ_landscape$mean.area, areaSD = occ_landscape$SD.area, 
                  Npatch = npatch, disp = occ_landscape$dispersal, plotG = FALSE)
                occ_landscape <- suppressWarnings(species.graph(rl = rl1, method = "percentage", 
                  parm = 0, plotG = FALSE))
                occ_landscape$nodes.characteristics[nrow_land + 2, ] <- NA
                occ_landscape$nodes.characteristics$ID[(nrow_land + 1):(nrow_land + 
                  2)] <- (ID_land + 1):(ID_land + 2)
                area1 <- mean(occ_landscape$nodes.characteristics$areas[1:nrow_land])
                occ_landscape$nodes.characteristics$radius[(nrow_land + 1):(nrow_land + 
                  2)] <- sqrt((area1 * 10000)/pi)
                occ_landscape$nodes.characteristics$areas[(nrow_land + 1):(nrow_land + 
                  2)] <- (occ_landscape$nodes.characteristics$radius[(nrow_land + 
                  1)]) * mapsize
                occ_landscape$nodes.characteristics$y[nrow_land + 1] <- xx3xx
                occ_landscape$nodes.characteristics$x[nrow_land + 1] <- xx2xx
                occ_landscape$nodes.characteristics$y[nrow_land + 2] <- xx8xx
                occ_landscape$nodes.characteristics$x[nrow_land + 2] <- xx7xx
                max_cluster <- max(occ_landscape$nodes.characteristics$cluster, na.rm = T)
                occ_landscape$nodes.characteristics$cluster[(nrow_land + 1):(nrow_land + 
                  2)] <- c(max_cluster + 1)
                occ_landscape$nodes.characteristics$colour <- as.character(occ_landscape$nodes.characteristics$colour)
                occ_landscape$nodes.characteristics$colour[(nrow_land + 1)] <- xx4xx
                occ_landscape$nodes.characteristics$colour[(nrow_land + 2)] <- xx9xx
                occ_landscape$nodes.characteristics$species[(nrow_land + 1)] <- 0
                occ_landscape$nodes.characteristics$species[(nrow_land + 2)] <- 1
                occ_landscape$nodes.characteristics$nneighbour[(nrow_land + 1)] <- xx5xx
                occ_landscape$nodes.characteristics$nneighbour[(nrow_land + 2)] <- xx10xx
                occ_landscape$number.patches <- occ_landscape$number.patches + 2
                dist_nodos <- occ_landscape$distance.to.neighbours
                dist_nodos[, (nrow_land + 1)] <- NA
                dist_nodos[, (nrow_land + 2)] <- NA
                dist_nodos[(nrow_land + 1), ] <- NA
                dist_nodos[(nrow_land + 2), ] <- NA
                dist_nodos[, (nrow_land + 1)] <- c(xx6xx, 0, mapsize)
                dist_nodos[, (nrow_land + 2)] <- c(xx11xx, mapsize, 0)
                dist_nodos[(nrow_land + 1), ] <- c(xx6xx, 0, mapsize)
                dist_nodos[(nrow_land + 2), ] <- c(xx11xx, mapsize, 0)
                colnames(dist_nodos)[(nrow_land + 1)] <- xx4xx
                colnames(dist_nodos)[(nrow_land + 2)] <- xx9xx
                rownames(dist_nodos)[(nrow_land + 1)] <- xx4xx
                rownames(dist_nodos)[(nrow_land + 2)] <- xx9xx
                occ_landscape$distance.to.neighbours <- dist_nodos
            }
            class(occ_landscape) <- "metapopulation"
			occ_landscape_new <- spom(sp = occ_landscape, kern = "op1", 
                conn = "op1", colnz = "op1", ext = "op1", param_df = param, b = b, 
                c1 = NULL, c2 = NULL, z = NULL, R = NULL)
            if (occ_landscape_new$nodes.characteristics$species2[nrow_land + 1] == 
                0) 
                ocupp <- "N"
            if (occ_landscape_new$nodes.characteristics$species2[nrow_land + 1] == 
                1) 
                ocupp <- "Y"
            if (ocupp == "N") {
                occ_landscape <- occ_landscape_new
                occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[, 
                  -c(9, 11)]
                names(occ_landscape$nodes.characteristics)[names(occ_landscape$nodes.characteristics) == 
                  "species2"] <- "species"
                occ_landscape <- occ_landscape[-9]
            }
            if (ocupp == "Y") 
                value <- clock
            if (ocupp == "N") 
                value <- 0
            output0[clock, 2] <- value
            if (clock == tsteps) 
                break
            if (nrow(occ_landscape$nodes.characteristics) == (npatch + 1)) {
                occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[-(nrow_land + 
                  1), ]
                occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[-(nrow_land + 
                  1), ]
                occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[, 
                  -(nrow_land + 1)]
                occ_landscape$number.patches <- occ_landscape$number.patches - 1
            }
            if (nrow(occ_landscape$nodes.characteristics) == (npatch + 2)) {
                occ_landscape$nodes.characteristics <- occ_landscape$nodes.characteristics[-c((nrow_land + 
                  1):(nrow_land + 2)), ]
                occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[-c((nrow_land + 
                  1):(nrow_land + 2)), ]
                occ_landscape$distance.to.neighbours <- occ_landscape$distance.to.neighbours[, 
                  -c((nrow_land + 1):(nrow_land + 2))]
                occ_landscape$number.patches <- occ_landscape$number.patches - 2
            }
            clock <- clock + 1
        }
        if (sum(output0[, 2]) == 0) {
            out_vec <- 0  
            output1 <- cbind(0, out_vec)
            colnames(output1) <- c("DISTANCE", "TIME STEP")
            output1 <- as.data.frame(output1)
        }
        if (sum(output0[, 2]) != 0) {
            output1 <- output0[output0[, 2] != 0, ]  
            output1[, 1] <- output1[, 1] * mapsize
            colnames(output1) <- c("DISTANCE", "TIME STEP")
            output1 <- as.data.frame(output1)
            return(output1)
        }
    }
    distance <- mapsize * 1:tsteps
    outputN <- distance
    outputS <- distance
    outputE <- distance
    outputW <- distance
    for (i in 1:iter) {
        sp1 <- species.graph(rl = rl, method = "percentage", parm = percI, nsew = "none", 
            plotG = FALSE)
        nodeN <- node.expansion(occ_landscape = sp1, param, b, node = "North", 
            tsteps)
        outputN <- suppressWarnings(cbind(outputN, c(nodeN[, 2], rep(NA, length(outputN) - 
            length(nodeN[, 2])))))
        nodeS <- node.expansion(occ_landscape = sp1, param, b, node = "South", 
            tsteps)
        outputS <- suppressWarnings(cbind(outputS, c(nodeS[, 2], rep(NA, length(outputS) - 
            length(nodeS[, 2])))))
        nodeE <- node.expansion(occ_landscape = sp1, param, b, node = "East", 
            tsteps)
        outputE <- suppressWarnings(cbind(outputE, c(nodeE[, 2], rep(NA, length(outputE) - 
            length(nodeE[, 2])))))
        nodeW <- node.expansion(occ_landscape = sp1, param, b, node = "West", 
            tsteps)
        outputW <- suppressWarnings(cbind(outputW, c(nodeW[, 2], rep(NA, length(outputW) - 
            length(nodeW[, 2])))))
    }
    outputN[is.na(outputN)] <- 0
    outputS[is.na(outputS)] <- 0
    outputE[is.na(outputE)] <- 0
    outputW[is.na(outputW)] <- 0
    for (x in 2:(iter + 1)) {
        i_N <- which(outputN[, x] != 0)
        outputN[i_N, x] <- 1
        i_S <- which(outputS[, x] != 0)
        outputS[i_S, x] <- 1
        i_E <- which(outputE[, x] != 0)
        outputE[i_E, x] <- 1
        i_W <- which(outputW[, x] != 0)
        outputW[i_W, x] <- 1
    }
    outputN <- cbind(outputN[, 1], rowSums(outputN[, 2:(ncol(outputN))]))
    outputS <- cbind(outputS[, 1], rowSums(outputS[, 2:(ncol(outputS))]))
    outputE <- cbind(outputE[, 1], rowSums(outputE[, 2:(ncol(outputE))]))
    outputW <- cbind(outputW[, 1], rowSums(outputW[, 2:(ncol(outputW))]))
    outputN <- cbind(outputN[, 1:2], outputN[, 2]/iter)
    outputS <- cbind(outputS[, 1:2], outputS[, 2]/iter)
    outputE <- cbind(outputE[, 1:2], outputE[, 2]/iter)
    outputW <- cbind(outputW[, 1:2], outputW[, 2]/iter)
    outputN <- as.data.frame(outputN)
    outputS <- as.data.frame(outputS)
    outputE <- as.data.frame(outputE)
    outputW <- as.data.frame(outputW)
    names(outputN) <- c("DISTANCE", "OCCUPATION", "PROPORTION")
    names(outputS) <- c("DISTANCE", "OCCUPATION", "PROPORTION")
    names(outputE) <- c("DISTANCE", "OCCUPATION", "PROPORTION")
    names(outputW) <- c("DISTANCE", "OCCUPATION", "PROPORTION")
    p1 <- gvisLineChart(outputN, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the North", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#006400'}]", 
        backgroundColor = "#D1EEEE"))
    p2 <- gvisLineChart(outputS, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the South", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#0000FF'}]", 
        backgroundColor = "#D1EEEE"))
    p3 <- gvisLineChart(outputE, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the East", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#8B0000'}]", 
        backgroundColor = "#D1EEEE"))
    p4 <- gvisLineChart(outputW, xvar = "DISTANCE", yvar = "PROPORTION", options = list(title = "Dispersal to the West", 
        width = 500, height = 300, curveType = "function", legend = "none", titleTextStyle = "{colour:'black', fontName:'Courier', fontSize:16}", 
        vAxis = "{title: 'proportion'}", hAxis = "{title: 'distance(meters)'}", series = "[{color: '#FF4500'}]", 
        backgroundColor = "#D1EEEE"))
    ln1 <- gvisMerge(p1, p2, horizontal = TRUE)
    ln2 <- gvisMerge(p3, p4, horizontal = TRUE)
    ln.final <- gvisMerge(ln1, ln2, horizontal = FALSE)
    ln.final$html$caption <- paste("<div><span>Range expansion in the four cardinal directions, considering that ", 
        percI, "% of the patches are occupied in the first landscape mosaic.</span><br />", 
        sep = "")
    ln.final$html$footer <- paste("\n<!-- htmlFooter -->\n<span> \n  ",R.Version()$version.string,"&#8226; <a href=\"http://code.google.com/p/google-motion-charts-with-r/\">googleVis-", packageVersion("googleVis"),"</a>\n  &#8226; MetaLandSim-",packageVersion("MetaLandSim"),"\n  &#8226; <a href=\"https://developers.google.com/terms/\">Google Terms of Use</a> &#8226; <a href=\"https://google-developers.appspot.com/chart/interactive/docs/gallery/linechart.html#Data_Policy\">Data Policy</a>\n</span></div>\n</body>\n</html>\n", sep="")
    output <- list(NORTH = outputN, SOUTH = outputS, EAST = outputE, WEST = outputW)
    plot(ln.final)
    class(output) <- "expansion"
    return(output)
  }
