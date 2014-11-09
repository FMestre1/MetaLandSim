list.stats <-
function (sim_list, stat, plotG)
  {
    sim_list1 <- sim_list[[1]]
    n <- length(sim_list1)
    stat.list <- as.list(rep("", n))
    if(stat == "mean_area")
      {
        for(i in 1:n)
          {
            df_nodes <- sim_list1[[i]]
            mean_area <- as.numeric(mean(df_nodes[, "areas"]))
            stat.list[[i]] <- mean_area
          }
      }
    if(stat == "sd_area")
      {
        for(i in 1:n)
          {
            df_nodes <- sim_list1[[i]]
            sd_area <- as.numeric(sd(df_nodes[, "areas"]))
            stat.list[[i]] <- sd_area
          }
      }
    if(stat == "mean_distance")
      {
        for(i in 1:n)
          {
            df_nodes <- sim_list1[[i]]
            xy <- df_nodes[,1:2]
            xy_dist <- pairdist(as.matrix(xy))
			radiuses <- df_nodes$radius
				for(j in 1:nrow(xy_dist))
					{
					xy_dist[,j] <- xy_dist[,j] - radiuses[j]
					xy_dist[j,] <- xy_dist[j,] - radiuses[j]
					}
            xy_dist <- replace(xy_dist, xy_dist<0, 0)
			xy_dist2 <- xy_dist[upper.tri(xy_dist)]
            mean_distance <- mean(xy_dist2)
            stat.list[[i]] <- mean_distance
          }
      }
    if(stat == "n_patches")
      {
        for(i in 1:n)
          {
            npatches <- nrow(sim_list1[[i]])
            stat.list[[i]] <- npatches
          }
      }
    if(stat == "occupation")
      {
        for(i in 1:n)
          {
            occupation <- ((sum(sim_list1[[i]]$species)*100)/nrow(sim_list1[[i]]))
            stat.list[[i]] <- occupation
          }
      }
    if(stat == "turnover")
      {
        stat.list <- as.vector(sim_list[[2]])
      }
    stat.list2 <- as.numeric(stat.list)
    if (plotG == TRUE)
      {
        plot(stat.list2, type="l", col="darkgreen", xlab="time", ylab=stat)
      }
    return(stat.list2)
  }
