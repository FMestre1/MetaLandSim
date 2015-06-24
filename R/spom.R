spom <-
function(sp,kern,conn,colnz,ext,param_df,beta1=NULL,b=1,c1=NULL,c2=NULL,z=NULL,R=NULL)
  {

  if (class(sp)!="metapopulation") 
  {
  stop(paste(sp, " should be an object of class class 'metapopulation'.", sep=""), call. = FALSE)
  }
    alpha <- param_df[1,]
    y <- param_df[3,]
    e <- param_df[4,]
    x <- param_df[2,]
    dfsp <- sp$nodes.characteristics
    A <- dfsp$areas
    p <- dfsp$species
    dist1 <- dist(dfsp[,1:2])
    if(kern == "op1")
      {
        kern_m <- as.matrix(exp(-alpha*dist1))
        diag(kern_m) <- 0
        kern_m <- as.data.frame(kern_m)
      }
    if(kern == "op2")
      {
        kern_m <- as.matrix(1/(1+alpha*(dist1^beta1)))
        diag(kern_m) <- 0
        kern_m <- as.data.frame(kern_m)
      }
    if(conn == "op1")
      {
        Si <- sweep(kern_m, 2, A^b, "*")
        S <- as.vector(rowSums(Si[, p > 0, drop=FALSE]))
      }
    if(conn == "op2")
      {
        Si <- sweep(kern_m, 2, A^b, "*")
        S <- as.vector(rowSums(Si[, p > 0, drop=FALSE]))
        S <- (A^c1)*S
      }
    if(colnz == "op1")
      {
        C <- S^2/((S^2)+y)
      }
    if(colnz == "op2")
      {
        C <- 1 - exp(-y*S)
      }
    if(colnz == "op3")
      {
        C <- S^z/(S^z+(1/c2))
      }
    if(ext == "op1")
      {
        E <- e/(A^x)
        E <- ifelse(E>1, 1, E)
      }
    if(ext == "op2")
      {
        E <- 1-((-e)/(A^x))
      }
    if(ext == "op3")
      {
        E <- (e/A^x)*(1-C)^R
        E <- ifelse(E>1, 1, E)
      }
    
	cond <- ifelse(p, (1-C)*E, C)
    species2 <- ifelse(runif(length(p)) < cond, !p, p)
    turn <- ifelse (dfsp$species!=species2, 1, 0)
    nr_turn <- sum(turn)
    p3 <- cbind(sp$nodes.characteristics, species2, turn)
    p4 <- list(mapsize=sp$mapsize, minimum.distance=sp$minimum.distance,
               mean.area=sp$mean.area, SD.area=sp$SD.area,
               number.patches=sp$number.patches, dispersal=sp$dispersal,
               distance.to.neighbours=sp$distance.to.neighbours,
               nodes.characteristics=p3,turnover=nr_turn)
    return(p4)
  }
