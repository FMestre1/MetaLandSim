manage_landscape_sim <-
function(par_df,parameters_spom)
  {
    
	output <- data.frame(matrix(nrow=nrow(par_df),ncol=5))
    colnames(output) <- c("mean occupation","mean number of patches","mean turnover",
                          "mean distance","mean area")
	lines_df <- nrow(par_df)
    for(i in 1:nrow(par_df))
      {
        
       it <- iterate.graph(iter=par_df[i,8],
        mapsize=par_df[i,6],
        dist_m=par_df[i,2],
        areaM=par_df[i,4],
        areaSD=par_df[i,5],
        Npatch=par_df[i,3],
        disp=par_df[i,27],
        span=par_df[i,7],
        par1=as.vector(par_df[i,9]),
        par2=as.vector(par_df[i,10]),
        par3=as.vector(par_df[i,11]),
        par4=as.vector(par_df[i,12]),
        par5=as.vector(par_df[i,13]),
        method=as.vector(par_df[i,16]),
        parm=par_df[i,15],
        nsew=as.vector(par_df[i,14]),
        a_min=par_df[i,1],
	    param_df=parameters_spom,
	    kern=as.vector(par_df[i,17]),
	    conn=as.vector(par_df[i,18]),
	    colnz=as.vector(par_df[i,19]),
	    ext=as.vector(par_df[i,20]),
	    beta1=as.vector(par_df[i,21]),
	    b=par_df[i,22],
	    c1=as.vector(par_df[i,23]),
	    c2=as.vector(par_df[i,24]),
	    z=as.vector(par_df[i,25]),
	    R=as.vector(par_df[i,26]),
	    graph=FALSE)
        
		output[i,1] <- mean(it$occupancy[,"mean"])
        output[i,2] <- mean(it$number_patches[,"mean"])
        output[i,3] <- mean(it$turnover[,"mean"])
        output[i,4] <- mean(it$mean_distance[,"mean"])
        output[i,5] <- mean(it$mean_area[,"mean"])
		
		cat("Completed simulation",i," of ",lines_df,"\n")
      }
	  
    output2 <- cbind(par_df,output)
    return(output2)
  }
