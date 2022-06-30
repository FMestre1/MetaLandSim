metrics.graph <-
function(rl, metric,dispersal.dist=NULL)
  {
  if(!inherits(rl, "landscape"))
  #if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
  
  if (!rl$dispersal > 0)stop("Some metrics require the dispersal of the focal species to be defined!")
  if (is.null(dispersal.dist))
  { 
  dispersal.dist <- rl$dispersal
  message("The function will assume the dispersal distance provided in the 'landscape' object.")
  }
  result <- c()
  
    if("NC" %in% metric)
      {
	    rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        result <- c(result,NC = components.graph(rl))
      }
    if("LNK" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        result <- c(result, LNK = nrow(edge.graph(rl)))
      }
    if("SLC" %in% metric)
      {
	    rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        df0 <- rl$nodes.characteristics
        NC <- components.graph(rl)
        area_sum <- rep(NA, NC)
        for(i in 1:NC)
          {
            component <- df0[df0$cluster==i, ]
            area_sum[i] <- sum(component$areas)
          }
        result <- c(result, SLC = max(area_sum))
      }
    if("MSC" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        df0 <- rl$nodes.characteristics
        NC <- components.graph(rl)
        area_sum <- rep(NA, NC)
        for(i in 1:NC)
          {
            component <- df0[df0$cluster==i, ]
            area_sum[i] <- sum(component$areas)
          }
        result <- c(result, MSC = mean(area_sum))
      }
    if("HI" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
	    topological_matrix <-matrix.graph(rl, mat="top_matrix")
		topological_matrix[topological_matrix==0] <- Inf
		reciprocal <- 1/topological_matrix
		result <- c(result, HI = sum(reciprocal,na.rm=TRUE))
     }
    if("NH" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
		topological_matrix <-matrix.graph(rl, mat="top_matrix")
		topological_matrix[topological_matrix==0] <- Inf
		diag(topological_matrix) <- NA
		reciprocal <- 1/topological_matrix
		H <- sum(reciprocal,na.rm=T)
		n1 <- nrow(rl$nodes.characteristics)
		H_planar <- (n1*(n1+5)/4)-3
		H_chain <- c()
		for(i in 1:n1){
		v2 <- (n1-i)/n1
		H_chain <- cbind(H_chain,v2)		
		}
		H_chain <- c(as.vector(H_chain),1/(n1-1))
		H_chain <- sum(H_chain)
		result <- c(result, NH = (H-H_chain)/(H_planar-H_chain))
	  }
    if("ORD" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        a <- cluster.graph(rl)
        result <- c(result, ORD = max(a[, 2]))
      }
    if("GD" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
		d0 <- rl$nodes.characteristics
        disp <- rl$dispersal
        m2 <- as.matrix(d0[, 1:2])
        m3 <- pairdist(m2)
        m3[m3>disp] <- 0
        m3[m3 == 0] <- NA
        m4 <- allShortestPaths(m3)
		result <- c(result, GD = max(m4$length, na.rm=TRUE))
	}
    if("CCP" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        df0 <- rl$nodes.characteristics
        NC <- components.graph(rl)
        Ac <- sum(df0$areas)
        r0 <- rep(NA, NC)
        for(i in 1:NC)
          {
            df1 <- df0[df0$cluster==i, ]
            ci <- sum(df1$areas)
            r0[i] <- (ci/Ac)^2
          }
        result <- c(result, CCP = sum (r0))
      }
    if("LCP" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        NC <- components.graph(rl)
        df0 <- rl$nodes.characteristics
        r0_vec <- rep(NA, NC)
        AL <- ((rl$mapsize)^2)/10000
        for(i in 1:NC)
          {
            df1 <- df0[df0$cluster==i, ]
            ci <- sum(df1$areas)
            r0_vec[i] <- (ci/AL)^2
          }
		result <- c(result, LCP = sum(r0_vec))
      }
    if("CPL" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        m0 <- min_distance(rl)
		diag(m0) <- NA
		result <- c(result, CPL = mean(m0,na.rm=TRUE))
      }
    if("ECS" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        NC <- components.graph(rl)
        df0 <- rl$nodes.characteristics
        ai <- rep(NA, NC)
        for(i in 1:NC)
          {
            df1 <- df0[df0$cluster==i, ]
            ai[i] <- sum((df1$areas)^2)
          }
        a_num <- sum(ai)
        a <- sum(df0$areas)
        result <- c(result, ECS = a_num/a)
      }
    if("AWF" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        disp <- rl$dispersal
        d0 <- rl$nodes.characteristics
        nnodes <-rl$number.patches
        distP <- pairdist (d0[,1:2])
        k <- -(log(0.5)/(disp/2))
		pij <- as.data.frame(exp(-k*(distP)))
        names(pij) <- d0$ID
        rownames(pij) <- d0$ID
        paths <- as.data.frame(t(combn(x=d0$ID, m=2)))
		paths2 <- as.data.frame(cbind(paths[,2],paths[,1]))
		paths <- rbind(paths, paths2)
		for(i in 1: nrow(paths)) paths[i, 3] <- d0$areas[d0$ID %in% paths[i, 1]]
        for(i in 1: nrow(paths)) paths[i, 4] <- d0$areas[d0$ID %in% paths[i, 2]]
        for(i in 1: nrow(paths))paths[i, 5] <- pij[as.character(paths[i,1]),as.character(paths[i,2])]
        result <- c(result, AWF = (sum(paths[, 3]*paths[, 4]*paths[, 5])))
      }
    if("IIC" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
        dist_tp <- min_distance(rl)
        df0 <- rl$nodes.characteristics
        paths <- as.data.frame(which(min_distance(rl)!=0, arr.ind = TRUE, useNames = FALSE))
		repl <- cbind(c(1:nrow(df0)),df0$ID)
		paths[,1] <- repl[paths[,1],2]
        paths[,2] <- repl[paths[,2],2]
		sameID <- cbind(as.numeric(df0$ID), as.numeric(df0$ID))
        paths <-rbind(paths,sameID)
        names(paths)[names(paths) == "V1"] <- "node_A"
        names(paths)[names(paths) == "V2"] <- "node_B"
		for(i in 1:nrow(paths)) paths[i, 3] <- dist_tp[as.character(paths[i,1]),as.character(paths[i,2])]
        names(paths)[names(paths) == "V3"] <- "top_distance"
		for(f in 1:nrow(paths))paths[f, 4] <- df0$areas[df0$ID %in% paths[f, 1]]
        names(paths)[names(paths) == "V4"] <- "area_A"
		for(x in 1:nrow(paths))paths[x, 5] <- df0$areas[df0$ID %in% paths[x, 2]]
        names(paths)[names(paths) == "V5"] <- "area_B"
		paths[, 6] <- (paths[, 4]*paths[, 5])/(1+(paths[, 3])) 
        Al2 <- (((rl$mapsize)^2)/10000)^2
		result <- c(result, IIC = sum(paths[, 6])/Al2)
	}
    if("PC" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
		nodes <- rl$nodes.characteristics
		nodes <- cbind(nodes$ID, nodes$areas)
		colnames(nodes) <- c("ID","area")
		edges <- edge.graph(rl)
		edges <- cbind(edges[,2], edges[,3], edges[,10])
		edges <- as.data.frame(edges)
		colnames(edges) <- c("node_A","node_B","distance")
		pijs <- edges  
		pijs <- as.data.frame(pijs)
		colnames(pijs) <- c('nodei','nodej','prob_ij')
		distmed <- rl$dispersal
		alfa <- distmed
		pijs$prob_ij <- 1-2*(Igamma(1/1,(edges$distance/alfa)^1)/(2*gamma(1/1)))
		mygraph <- graph_from_data_frame(pijs, directed=FALSE, vertices=nodes)
		E(mygraph)$weight <- -log(E(mygraph)$prob_ij)
		mxprobpath_PC <- shortest.paths(mygraph, weights = NULL)
		pijast <- exp(-mxprobpath_PC)
		ai_aj_pijast <- pijast
		i <- 1
		j <- 1
		while (i<=gorder(mygraph)) 
	{
			while (j<=gorder(mygraph)) 
				{
					ai_aj_pijast[i,j] <- ai_aj_pijast[i,j]*nodes[i,2]*nodes[j,2]
					j <- j+1
				}
		j <- 1
		i <- i+1
	}
		PCnum <- sum(ai_aj_pijast)
        Al2 <- (((rl$mapsize)^2)/10000)
        result <- c(result, PC = (PCnum/Al2^2))
     }
    if("ECA" %in% metric)
      {
	  	rl$dispersal  <- dispersal.dist
		rl <- cluster.id(rl)
		nodes <- rl$nodes.characteristics
		nodes <- cbind(nodes$ID, nodes$areas)
		colnames(nodes) <- c("ID","area")
		edges <- edge.graph(rl)
		edges <- cbind(edges[,2], edges[,3], edges[,10])
		edges <- as.data.frame(edges)
		colnames(edges) <- c("node_A","node_B","distance")
		pijs <- edges  
		pijs <- as.data.frame(pijs)
		colnames(pijs) <- c('nodei','nodej','prob_ij')
		distmed <- rl$dispersal
		alfa <- distmed
		pijs$prob_ij <- 1-2*(Igamma(1/1,(edges$distance/alfa)^1)/(2*gamma(1/1)))
		mygraph <- graph_from_data_frame(pijs, directed=FALSE, vertices=nodes)
		E(mygraph)$weight <- -log(E(mygraph)$prob_ij)
		mxprobpath_PC <- shortest.paths(mygraph, weights = NULL)
		pijast <- exp(-mxprobpath_PC)
		ai_aj_pijast <- pijast
		i <- 1
		j <- 1
		while (i<=gorder(mygraph)) 
	{
			while (j<=gorder(mygraph)) 
				{
					ai_aj_pijast[i,j] <- ai_aj_pijast[i,j]*nodes[i,2]*nodes[j,2]
					j <- j+1
				}
		j <- 1
		i <- i+1
	}
		PCnum <- sum(ai_aj_pijast)
		EC <- sqrt(sum(ai_aj_pijast))
        result <- c(result, ECA = EC)

	  
	  }    
	
	return(round(result,5))

}
