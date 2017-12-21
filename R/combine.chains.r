combine.chains<-function(x1,x2,nburnin,nthin=1,z.thin=TRUE) {
  x=x1
  chain.length=length(x$alpha.chain)
  nburnin.thin=round(nburnin/nthin)+1
  for(i in 1:length(x)){
    if(is.vector(x[[i]])) x[[i]]=c(x[[i]][nburnin.thin:chain.length],x2[[i]][nburnin.thin:chain.length])
    if(is.matrix(x[[i]]))  x[[i]]=cbind(x[[i]][,nburnin.thin:chain.length],x2[[i]][,nburnin.thin:chain.length])
    if(is.array(x[[i]]) & !is.matrix(x[[i]])) {
      nsite = dim(x[[i]])[1]
      nyear = dim(x[[i]])[2]
      if(z.thin==TRUE) {
        x[[i]]=c(x[[i]][,,seq(nburnin.thin,dim(x1[[i]])[3],by=5)],x2[[i]][,,seq(nburnin.thin,dim(x1[[i]])[3],by=5)])
      } else {
        x[[i]]=c(x[[i]][,,nburnin.thin:dim(x1[[i]])[3]],x2[[i]][,,nburnin.thin:dim(x1[[i]])[3]])
      }
      dim(x[[i]])=c(nsite,nyear,length(x$z.chain)/(nsite*nyear))
    }
  }
  x
}