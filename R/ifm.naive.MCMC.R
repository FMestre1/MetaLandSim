ifm.naive.MCMC <- function(niter=1000,init,z.data, site.distance, site.area, sd.prop.e=0.2, sd.prop.x=0.5,sd.prop.y=10, sd.prop.b=0.2, sd.prop.alpha=5,nthin=1,print.by=100) {

  iter.chain=seq(1,niter,by=nthin)
  iter.chain=sort(rep(iter.chain,nthin))
  iter.print=seq(1,niter,by=print.by)
  iter.print=sort(rep(iter.print,print.by))

  e.chain=NULL
  x.chain=NULL
  y.chain=NULL
  b.chain=NULL
  alpha.chain=NULL
  deviance.chain=NULL

  current.e=init$e
  current.x=init$x

  # code will use y^2 as the parameter rather than Hanski's parameterization.
  # The output will be corrected to correspond to Hanski's parameterization.
  current.y=init$y^2
  current.b=init$b
  current.alpha = init$alpha

  z.index.list=list()
  z.data.list=list()
  nyear=ncol(z.data)
  nsite = nrow(z.data)
  
  #         for(i in (nyear-1):1){
  #           temp1=c(1:nsite)[!is.na(z.data[,i]) & !is.na(z.data[,(i+1)])]
  #           assign(paste("z.",i,".",(i+1),".index",sep=""),temp1,z.index.list)
  #           temp2=z.data[temp1,c(i,(i+1))]
  #           assign(paste("z.",i,".",(i+1),sep=""),temp2,z.data.list)
  #         }
  #         z.index.list=as.list(z.index.list)
  #         z.data.list=as.list(z.data.list)
  #
  #         if(names(z.index.list)[1]!="z.1.2.index") stop("Z INDEX MIS-SORTED")
  #

  for(i in 1:(nyear-1)){
    temp1=c(1:nsite)[!is.na(z.data[,i]) & !is.na(z.data[,(i+1)])]
    z.index.list[[i]] = temp1
    names(z.index.list)[i] = paste("z.",i,".",(i+1),".index",sep="")

    temp2=z.data[temp1,c(i,(i+1))]
    z.data.list[[i]] = temp2
    names(z.data.list)[i] = paste("z.",i,".",(i+1),sep="")
  }



  if(names(z.index.list)[1]!="z.1.2.index") stop("Z INDEX MIS-SORTED")

  current.sim.alpha.distance.list=list(rep(NA,nyear))
  current.area.j.b=site.area^current.b
  current.area.occ.list=list(rep(NA,nyear))
  current.s.i.sq.list=list(rep(NA,nyear))

  z.data.1.nminus1.stack=NULL
  z.data.2.n.stack=NULL
  current.p.colon.stack=NULL
  current.p.extinct.stack=NULL

  current.site.area.x=site.area^current.x
  current.p.extinct=current.e/current.site.area.x
  current.p.extinct[current.p.extinct>1]=1
  current.sim.alpha.distance=exp(-current.alpha*site.distance)
  diag(current.sim.alpha.distance)=0

  for(i in 1:(nyear-1)) {
    z.data.1.nminus1.stack=c(z.data.1.nminus1.stack,z.data.list[[i]][,1])
    z.data.2.n.stack=c(z.data.2.n.stack,z.data.list[[i]][,2])
    current.sim.alpha.distance.list[[i]]=current.sim.alpha.distance[z.index.list[[i]],z.index.list[[i]]]
    current.area.occ.list[[i]]=z.data.list[[i]][,1]*current.area.j.b[z.index.list[[i]]]
    current.s.i.sq.list[[i]]=(current.sim.alpha.distance.list[[i]]%*%current.area.occ.list[[i]])^2
    current.p.colon.stack=c(current.p.colon.stack,(current.s.i.sq.list[[i]]/(current.s.i.sq.list[[i]]+current.y)))
    current.p.extinct.stack=c(current.p.extinct.stack,current.p.extinct[z.index.list[[i]]])
  }

  current.p.extinct.re.stack=current.p.extinct.stack*(1-current.p.colon.stack)
  current.prob.stack=z.data.1.nminus1.stack*(1-current.p.extinct.re.stack)+(1-z.data.1.nminus1.stack)*current.p.colon.stack
  log.denominator.2.n=sum(z.data.2.n.stack*log(current.prob.stack)+(1-z.data.2.n.stack)*log(1-current.prob.stack))

  proposed.sim.alpha.distance.list=current.sim.alpha.distance.list
  proposed.area.occ.list=current.area.occ.list
  proposed.s.i.sq.list=current.s.i.sq.list

  zero.stack=rep(0.00001,length(current.prob.stack))
  one.stack=rep(0.99999,length(current.prob.stack))


  #-----------------------
  # ITERATIONS LOOP
  #-----------------------
  for (n in 1:niter) {


    #-----------------
    #-----------e
    #-----------------

    if(sd.prop.e>0) {
      proposed.e=rnorm(1,current.e,sd=sd.prop.e)

      ## This is where you control the priors. Note the
      ## proposed.e is replaced by 2*(upper bound) - proposed.e
      ## when the proposal is outside the prior's support.
      if(proposed.e<0) {
        proposed.e=-proposed.e
      }
      if(proposed.e>1) {
        proposed.e=2-proposed.e
      }
      proposed.p.extinct=proposed.e/current.site.area.x
      proposed.p.extinct[proposed.p.extinct>1]=1
      proposed.p.extinct.stack=NULL
      for(i in 1:(nyear-1)) {
        proposed.p.extinct.stack=c(proposed.p.extinct.stack,proposed.p.extinct[z.index.list[[i]]])
      }

      proposed.p.extinct.re.stack=proposed.p.extinct.stack*(1-current.p.colon.stack)
      proposed.prob.stack=z.data.1.nminus1.stack*(1-proposed.p.extinct.re.stack)+(1-z.data.1.nminus1.stack)*current.p.colon.stack
      log.numerator.2.n=sum(z.data.2.n.stack*log(proposed.prob.stack)+(1-z.data.2.n.stack)*log(1-proposed.prob.stack))
      log.MH=log.numerator.2.n-log.denominator.2.n
      accept=decide(log.MH)
      if (accept) {
        current.e=proposed.e
        current.p.extinct.stack=proposed.p.extinct.stack
        log.denominator.2.n=log.numerator.2.n
      }
    }
    #-----------------
    #-----------x
    #-----------------


    if(sd.prop.x>0){
      proposed.x=rnorm(1,current.x,sd=sd.prop.x)
      if(proposed.x<0) {
        proposed.x=-proposed.x
      }
      if(proposed.x>5) {
        proposed.x=10-proposed.x
      }

      proposed.site.area.x=site.area^proposed.x
      proposed.p.extinct=current.e/proposed.site.area.x
      proposed.p.extinct[proposed.p.extinct>1]=1
      proposed.p.extinct.stack=NULL
      for(i in 1:(nyear-1)) {
        proposed.p.extinct.stack=c(proposed.p.extinct.stack,proposed.p.extinct[z.index.list[[i]]])
      }

      proposed.p.extinct.re.stack=proposed.p.extinct.stack*(1-current.p.colon.stack)
      proposed.prob.stack=z.data.1.nminus1.stack*(1-proposed.p.extinct.re.stack)+(1-z.data.1.nminus1.stack)*current.p.colon.stack
      log.numerator.2.n=sum(z.data.2.n.stack*log(proposed.prob.stack)+(1-z.data.2.n.stack)*log(1-proposed.prob.stack))
      log.MH=log.numerator.2.n-log.denominator.2.n
      accept=decide(log.MH)
      if (accept) {
        current.x=proposed.x
        current.site.area.x=proposed.site.area.x
        current.p.extinct.stack=proposed.p.extinct.stack
        log.denominator.2.n=log.numerator.2.n
      }
    }

    #-----------------
    #-----------y
    #-----------------

    ## HERE I DO NOT SQUARE y
    if(sd.prop.y>0) {
      proposed.y=rnorm(1,current.y,sd=sd.prop.y)
      if(proposed.y<0) {
        proposed.y=-proposed.y
      }
      if(proposed.y>400) {
        proposed.y=800-proposed.y
      }

      proposed.p.colon.stack=NULL
      for(i in 1:(nyear-1)) {
        proposed.p.colon.stack=c(proposed.p.colon.stack,(current.s.i.sq.list[[i]]/(current.s.i.sq.list[[i]]+proposed.y)))
      }

      proposed.p.extinct.re.stack=current.p.extinct.stack*(1-proposed.p.colon.stack)
      proposed.prob.stack=z.data.1.nminus1.stack*(1-proposed.p.extinct.re.stack)+(1-z.data.1.nminus1.stack)*proposed.p.colon.stack
      log.numerator.2.n=sum(z.data.2.n.stack*log(proposed.prob.stack)+(1-z.data.2.n.stack)*log(1-proposed.prob.stack))
      log.MH=log.numerator.2.n-log.denominator.2.n
      accept=decide(log.MH)
      if(accept) {
        current.y=proposed.y
        current.p.colon.stack=proposed.p.colon.stack
        log.denominator.2.n=log.numerator.2.n
      }
    }
    #-----------------
    #-----------b
    #-----------------
    if(sd.prop.b>0) {
      proposed.b=rnorm(1,current.b,sd=sd.prop.b)
      if(proposed.b<0) {
        proposed.b=-proposed.b
      }
      #INFORMATIVE PRIOR.
      if(proposed.b>5) {
        proposed.b=10-proposed.b
      }
      proposed.area.j.b=site.area^proposed.b
      proposed.p.colon.stack=NULL
      for(i in 1:(nyear-1)) {
        proposed.area.occ.list[[i]]=z.data.list[[i]][,1]*proposed.area.j.b[z.index.list[[i]]]
        proposed.s.i.sq.list[[i]]=(current.sim.alpha.distance.list[[i]]%*%proposed.area.occ.list[[i]])^2
        proposed.p.colon.stack=c(proposed.p.colon.stack,(proposed.s.i.sq.list[[i]]/(proposed.s.i.sq.list[[i]]+current.y)))
      }

      proposed.p.extinct.re.stack=current.p.extinct.stack*(1-proposed.p.colon.stack)
      proposed.prob.stack=z.data.1.nminus1.stack*(1-proposed.p.extinct.re.stack)+(1-z.data.1.nminus1.stack)*proposed.p.colon.stack
      log.numerator.2.n=sum(z.data.2.n.stack*log(proposed.prob.stack)+(1-z.data.2.n.stack)*log(1-proposed.prob.stack))
      log.MH=log.numerator.2.n-log.denominator.2.n
      accept=decide(log.MH)
      if(accept) {
        current.b=proposed.b
        current.area.occ.list=proposed.area.occ.list
        current.s.i.sq.list=proposed.s.i.sq.list
        current.p.colon.stack=proposed.p.colon.stack
        log.denominator.2.n=log.numerator.2.n
      }


    }
    #-----------------
    #-----------ALPHA
    #-----------------
    # 2 corresponds to 10000/2 or 5 km maximum dispersal distance.
    if(sd.prop.alpha>0){
      proposed.alpha=rnorm(1,current.alpha,sd=sd.prop.alpha)
      if(proposed.alpha<2) {
        proposed.alpha=4-proposed.alpha
      }
      if(proposed.alpha>50) {
        proposed.alpha=100-proposed.alpha
      }
      proposed.sim.alpha.distance=exp(-proposed.alpha*site.distance)
      diag(proposed.sim.alpha.distance)=0
      proposed.p.colon.stack=NULL
      for(i in 1:(nyear-1)) {
        proposed.sim.alpha.distance.list[[i]]=proposed.sim.alpha.distance[z.index.list[[i]],z.index.list[[i]]]
        proposed.s.i.sq.list[[i]]=(proposed.sim.alpha.distance.list[[i]]%*%current.area.occ.list[[i]])^2
        proposed.p.colon.stack=c(proposed.p.colon.stack,(proposed.s.i.sq.list[[i]]/(proposed.s.i.sq.list[[i]]+current.y)))
      }

      proposed.p.extinct.re.stack=current.p.extinct.stack*(1-proposed.p.colon.stack)
      proposed.prob.stack=z.data.1.nminus1.stack*(1-proposed.p.extinct.re.stack)+(1-z.data.1.nminus1.stack)*proposed.p.colon.stack
      if(sum(proposed.prob.stack<=0)>0){
        print("Zeros in Prob")
        proposed.prob.stack[proposed.prob.stack<=0]=zero.stack[proposed.prob.stack<=0]
      }
      log.numerator.2.n=sum(z.data.2.n.stack*log(proposed.prob.stack)+(1-z.data.2.n.stack)*log(1-proposed.prob.stack))
      log.MH=log.numerator.2.n-log.denominator.2.n
      accept=decide(log.MH)
      if(accept) {
        current.alpha=proposed.alpha
        current.sim.alpha.distance.list=proposed.sim.alpha.distance.list
        current.s.i.sq.list=proposed.s.i.sq.list
        current.p.colon.stack=proposed.p.colon.stack
        log.denominator.2.n=log.numerator.2.n
      }

    }

    if(n==iter.chain[n]) {
      deviance.chain = c(deviance.chain,-2*log.denominator.2.n)
      e.chain=c(e.chain,current.e)
      x.chain=c(x.chain,current.x)
      y.chain=c(y.chain,current.y)
      b.chain=c(b.chain,current.b)
      alpha.chain=c(alpha.chain,current.alpha)
    }
    if(n==iter.print[n]) {
      print(n)
    }
  }

  # Take square root of y to return Hanski's parameterization:
  return(list(e.chain=e.chain,x.chain=x.chain,y.chain=sqrt(y.chain),
              b.chain=b.chain,alpha.chain=alpha.chain,deviance.chain=deviance.chain))
}
