calcmode <- function(data,adjust=1) {
  # Function for mode estimation of a continuous variable
  # data: vector used to estimate density 
  # adjust: increase this value to make the density estimate smoother.
  
  x<-data
  hist(x,freq=FALSE)
  s<-density(x,adjust=adjust)
  ymode <- max(s$y)
  modex = s$x[which.max(s$y)]
  lines(s$x,s$y,col="red")
  points(modex,ymode,col="blue")
  modex
}