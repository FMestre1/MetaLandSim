metrics.graph <-
function(rl, metric)
  {
	if (class(rl)!="landscape") 
  {
  stop(paste(rl, " should be an object of class class 'landscape'.", sep=""), call. = FALSE)
  }
  
  if (!rl$dispersal>0)stop("Some metrics require the mean dispersal of the focal species to be defined!")
    
  result <- c()
  
    if("NC" %in% metric)
      {
        result <- c(result,NC = components.graph(rl))
      }
    if("LNK" %in% metric)
      {
        result <- c(result, LNK = nrow(edge.graph(rl)))
      }
    if("SLC" %in% metric)
      {
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
	    topological_matrix <-matrix.graph(rl, mat="top_matrix")
		topological_matrix[topological_matrix==0] <- Inf
		reciprocal <- 1/topological_matrix
		result <- c(result, HI = sum(reciprocal,na.rm=TRUE))
     }
    if("NH" %in% metric)
      {
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
        a <- cluster.graph(rl)
        result <- c(result, ORD = max(a[, 2]))
      }
    if("GD" %in% metric)
      {
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
        m0 <- min_distance(rl)
		diag(m0) <- NA
		result <- c(result, CPL = mean(m0,na.rm=TRUE))
      }
    if("ECS" %in% metric)
      {
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
        disp <- rl$dispersal
        d0 <- rl$nodes.characteristics
        nnodes <-rl$number.patches
        distP <- pairdist (d0[,1:2])
		ext <- 1/disp
        pij <- as.data.frame(exp(-ext*(distP)))
        names(pij) <- d0$ID
        rownames(pij) <- d0$ID
        paths <- as.data.frame(t(combn(x=d0$ID, m=2)))
		for(i in 1: nrow(paths)) paths[i, 3] <- d0$areas[d0$ID %in% paths[i, 1]]
        for(i in 1: nrow(paths)) paths[i, 4] <- d0$areas[d0$ID %in% paths[i, 2]]
        for(i in 1: nrow(paths))paths[i, 5] <- pij[as.character(paths[i,1]),as.character(paths[i,2])]
        result <- c(result, AWF = (sum(paths[, 3]*paths[, 4]*paths[, 5])))
				
      }
    if("IIC" %in% metric)
      {
        dist_tp <- min_distance(rl)
        df0 <- rl$nodes.characteristics
        paths <- as.data.frame(which(min_distance(rl)!=0, arr.ind = TRUE, useNames = FALSE))
		repl <- cbind(c(1:nrow(df0)),df0$ID)
		paths[,1] <- repl[paths[,1],2]
        paths[,2] <- repl[paths[,2],2]
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
        d0 <- rl$nodes.characteristics
        disp <- rl$dispersal
        m2 <- as.matrix(d0[, 1:2])
        m3 <- pairdist(m2)
        m3[m3>disp] <- 0
        m3[m3 == 0] <- NA
        dmin <- min_distance(rl)
        connect <- as.data.frame(which(dmin!=0, arr.ind = TRUE, useNames = FALSE))
    	repl <- cbind(c(1:nrow(d0)),d0$ID)
		connect[,1] <- repl[connect[,1],2]
        connect[,2] <- repl[connect[,2],2]
        paths <- connect
		nr_conn <- nrow(paths)
        dists <- as.data.frame(pairdist(m2))
        colnames(dists) <- d0$ID
        rownames(dists) <- d0$ID
		for(r in 1:nr_conn) paths[r, 3] <- dists[as.character(paths[r, 1]),as.character(paths[r, 2])]
        for(r in 1:nr_conn)paths[r, 4] <- dmin[as.character(paths[r, 1]),as.character(paths[r, 2])]
        colnames(paths)[1] <- "node A"
        colnames(paths)[2] <- "node B"
        colnames(paths)[3] <- "distance"
        colnames(paths)[4] <- "top_distance"
        m4 <- allShortestPaths(m3)
		m4$length <- as.data.frame(m4$length)
		rownames(m4$length)<-as.character(d0$ID)
		colnames(m4$length)<-as.character(d0$ID)
    	m4$middlePoints <- as.data.frame(m4$middlePoints)
        rownames(m4$middlePoints)<-as.character(d0$ID)
		colnames(m4$middlePoints)<-as.character(d0$ID)
		m5 <- rep(NA, nr_conn)
        for(i in 1:nr_conn)
          {
            if (paths$top_distance[i]>1)
              {
			    ID <- cbind(1:length(d0$ID),d0$ID) 
                a <- extractPath(m4, ID[ID[,2]==paths[i, 1],1], ID[ID[,2]==paths[i, 2],1])
				for(j in 1:length(a)) a[j] <- ID[a[j],2]
				b <- length (a)-1
                b2 <- length(a)
                c1 <- as.data.frame(matrix(NA, b, 4))
                a2 <- rep(NA, (2*(b2-2))+2)
                a2[1] <- a[1]
                a2[length(a2)] <- a[length(a)]
                a3 <- rep(a[2:(length(a)-1)] ,each=2)
                for(s in 1:length(a3)) a2[s+1] <- a3[s]
                newpairs <- as.data.frame(matrix(a2, ncol=2, byrow=TRUE))
                ext <- 1/disp
                for(g in 1: nrow(newpairs))newpairs[g, 3] <- dists[as.character(newpairs[g, 1]),as.character(newpairs[g, 2])]
                for(h in 1: nrow(newpairs))newpairs[h, 4] <- exp(-ext*(newpairs[h, 3]))
                colnames(newpairs)[1] <- "node_A"
                colnames(newpairs)[2] <- "node_B"
                colnames(newpairs)[3] <- "distance"
                colnames(newpairs)[4] <- "probability"
                d <- prod(newpairs[, 4])
                m5[i] <- d
              }
            if(paths$top_distance[i] == 1)
              {
                ext <- 1/disp
                m5[i] <- exp(-paths$distance[i]*ext)
              }
            if(paths[i, 1] == paths[i, 2])
              {
                m5[i] <- 1
              }
          }
        paths2 <- cbind(paths, m5)
        colnames(paths2)[5] <- "pij*"
        for(i in 1:nrow(paths2))paths2[i, 6] <- d0$areas[d0$ID %in% paths2[i, 1]]
        for(i in 1:nrow(paths2))paths2[i, 7] <- d0$areas[d0$ID %in% paths2[i, 2]]
        colnames(paths2)[6] <- "areaA"
        colnames(paths2)[7] <- "areaB"
        paths2[, 8] <- paths2[, 5]*paths2[, 6]*paths2[, 7]
        Al2 <- (((rl$mapsize)^2)/10000)
        result <- c(result, PC = (sum(paths2[, 8]))/Al2^2)
     }
    return(round(result,3))
}