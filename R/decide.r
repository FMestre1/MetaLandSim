decide<-function(log.MH) {
  if(is.na(log.MH)){
    print(log.MH)
    log.MH=-Inf
  }
  if(runif(1,min=0,max=1) < exp(log.MH))	return(TRUE)
  else return(FALSE)
}