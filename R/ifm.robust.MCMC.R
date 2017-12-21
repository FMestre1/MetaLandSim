ifm.robust.MCMC <- function(niter=1000, init, det.data, site.distance, site.area, sd.prop.p=0.1,sd.prop.mupsi1=0.1, sd.prop.e=0.2, sd.prop.x=0.2,sd.prop.y=0.2, sd.prop.b=0.2, sd.prop.alpha=0.2,nthin=1,nsite.subset=5,print.by=100) {



	iter.chain=seq(1,niter,by=nthin)
	iter.chain=sort(rep(iter.chain,nthin))
        iter.print=seq(1,niter,by=print.by)
        iter.print=sort(rep(iter.print,print.by))

        nsite=dim(det.data)[1]
        nyear=dim(det.data)[2]

        z.chain=NULL
        p.chain=NULL
      	mupsi1.chain=NULL
      	e.chain=NULL
      	x.chain=NULL
      	y.chain=NULL
      	b.chain=NULL
      	alpha.chain=NULL
        latent.deviance.chain=NULL

        current.p=init$p
      	current.mupsi1=init$mupsi1
      	current.e=init$e
      	current.x=init$x

        # note: code will use y^2, but the output
        # will be transformed to use the parameterization in
        # Hanski.
      	current.y=init$y^2
      	current.b=init$b
      	current.alpha=init$alpha

        ntrial=rep(NA,nsite*nyear)
        nsucc=rep(NA,nsite*nyear)
        dim(ntrial)=c(nsite,nyear)
        dim(nsucc)=c(nsite,nyear)
        for(t in 1:nyear) {
          ntrial[,t]=apply(det.data[,t,],1,na.cnt)
          nsucc[,t]=apply(det.data[,t,],1,sum,na.rm=TRUE)
        }

	nobs=nsite*nyear
        z.data=init$z.data
	year.vector=as.vector(col(z.data))

	z.index=c(1:nobs)
        dim(z.index)=c(nsite,nyear)

        nbundles.yr=ceiling(nsite/nsite.subset)

# 	z.missing.index=z.index[is.na(z.data)]
#         z.missing.index.list=new.env()
#         for(i in nyear:1) {
#           temp=c(1:nobs)[is.na(z.data) & col(z.data)==i]
#           assign(paste("z.missing.index.",i,sep=""),temp,envir=z.missing.index.list)
# 	}
#         z.missing.index.list=as.list(z.missing.index.list)
#
	z.missing.index.list=list()
	for(i in 1:nyear) {
	  temp=c(1:nobs)[is.na(z.data) & col(z.data)==i]
    z.missing.index.list[[i]] = temp
	  names(z.missing.index.list)[i] = paste("z.missing.index.",i,sep="")
	}

        if(names(z.missing.index.list)[1]!="z.missing.index.1") stop("MISSING YEAR INDEXING MIS-SORTED")

#         z.index.list=new.env()
#         for(i in (nyear-1):1){
#           temp1=c(1:nsite)[!is.na(z.data[,i]) & !is.na(z.data[,(i+1)])]
#           assign(paste("z.",i,".",(i+1),".index",sep=""),temp1,z.index.list)
#         }
#         z.index.list=as.list(z.index.list)

      z.index.list=list()
      for(i in 1:(nyear-1)){
        temp1=c(1:nsite)[!is.na(z.data[,i]) & !is.na(z.data[,(i+1)])]
        z.index.list[[i]]= temp1
        names(z.index.list)[i] = paste("z.",i,".",(i+1),".index",sep="")
      }

       if(names(z.index.list)[1]!="z.1.2.index") stop("Z INDEX MIS-SORTED")
        n.pairs.obs=NULL
        for (i in 1:(nyear-1)){
        n.pairs.obs=c(n.pairs.obs,length(z.index.list[[i]]))
        }

	current.z=z.data
	current.z[is.na(z.data)]=init$z.missing
	current.z.1.nminus1=current.z[,1:(nyear-1)]
	current.z.2.n=current.z[,2:nyear]
	current.sim.alpha.distance=exp(-current.alpha*site.distance)
	diag(current.sim.alpha.distance)=0
	current.area.j.b=site.area^current.b
	current.area.occ.1.nminus1=current.z.1.nminus1*current.area.j.b
	current.s.i.sq.1.nminus1=(current.sim.alpha.distance%*%current.area.occ.1.nminus1)^2
	current.p.colon.2.n=current.s.i.sq.1.nminus1/(current.s.i.sq.1.nminus1+current.y)

	current.site.area.x=site.area^current.x
	current.p.extinct=current.e/current.site.area.x
	current.p.extinct[current.p.extinct>1]=1
	current.p.extinct.re.2.n=current.p.extinct*(1-current.p.colon.2.n)

	current.prob.mat.2.n=current.z.1.nminus1*(1-current.p.extinct.re.2.n)+(1-current.z.1.nminus1)*current.p.colon.2.n

	zero.m=rep(0.00001,nsite*(nyear-1))
	dim(zero.m)=c(nsite,(nyear-1))
	one.m=rep(0.99999,nsite*(nyear-1))
	dim(one.m)=c(nsite,(nyear-1))

	for (n in 1:niter) {
		#-------------------------
		#-----------LATENT STATES WITH MISSING DATA
		#-------------------------
          for(t in 1:nyear) {
		for(i in 1:nbundles.yr) {
				if(i<nbundles.yr) {
					temp.subset=z.index[(nsite.subset*(i-1)+1):(nsite.subset*i),t]
					}
				if(i==nbundles.yr) {
					temp.subset=z.index[(nsite.subset*(i-1)+1):nsite,t]
					}
				temp=prop.state.z(z.subset=temp.subset,year.vector=year.vector,current.p=current.p,
                                  current.mupsi1=current.mupsi1,current.prob.mat.2.n=current.prob.mat.2.n,
                                  current.z=current.z,current.area.j.b=current.area.j.b,
                                  current.sim.alpha.distance=current.sim.alpha.distance,current.y=current.y,
                                  current.p.extinct=current.p.extinct,nsucc=nsucc,ntrial=ntrial,nsite=nsite,nyear=nyear)
                                if (temp$accept) {
                                  current.z=temp$current.z
                                  current.prob.mat.2.n=temp$current.prob.mat.2.n
                                }
                                }
              }
        current.z.1.nminus1=current.z[,1:(nyear-1)]
        current.z.2.n=current.z[,2:nyear]
        current.area.occ.1.nminus1=current.z.1.nminus1*current.area.j.b
	current.s.i.sq.1.nminus1=(current.sim.alpha.distance%*%current.area.occ.1.nminus1)^2
	current.p.colon.2.n=current.s.i.sq.1.nminus1/(current.s.i.sq.1.nminus1+current.y)

                #---------
                # P
                #---------
          if(sd.prop.p!=0){
          for (t in 1:nyear) {
            log.denominator.p=log.lik.data.z.p(ntrial.temp=ntrial[,t],nsucc.temp=nsucc[,t],z.temp=current.z[,t],p.temp=current.p[t])
            proposed.p=rnorm(1,current.p[t],sd.prop.p)
         	if(proposed.p < 0.1) {
                  proposed.p=0.2-proposed.p
			}
		if(proposed.p >0.9999) {
			proposed.p=1.9998-proposed.p
		}
            #FOR THE UNLIKELY SCENARIO THAT IT REFLECTS TO A NEGATIVE VALUE
            if(proposed.p<0.1){
              proposed.p = 0.2 - proposed.p
            }

            log.numerator.p=log.lik.data.z.p(ntrial.temp=ntrial[,t],nsucc.temp=nsucc[,t],z.temp=current.z[,t],p.temp=proposed.p)
            log.MH=log.numerator.p-log.denominator.p
            accept=decide(log.MH)
            if(accept) current.p[t]=proposed.p
          }
        }
		#-----------------
		#-----------MUPSI1
		#-----------------
          if(sd.prop.mupsi1>0) {
                log.denominator.1=sum(current.z[,1]*log(current.mupsi1)+(1-current.z[,1])*log(1-current.mupsi1))
                proposed.mupsi1=rnorm(1,current.mupsi1,sd.prop.mupsi1)
		if(proposed.mupsi1 < 0.0001) proposed.mupsi1=0.0002-proposed.mupsi1
		if(proposed.mupsi1 >0.9999) proposed.mupsi1=1.9998-proposed.mupsi1
		log.numerator.1=sum(current.z[,1]*log(proposed.mupsi1)+(1-current.z[,1])*log(1-proposed.mupsi1))
		log.MH=log.numerator.1-log.denominator.1
		if(proposed.mupsi1<0) log.mh=-100
		accept=decide(log.MH)
		if(accept) current.mupsi1=proposed.mupsi1
          }


		#-----------------
		#-----------e
                #-----------------
          if(sd.prop.e>0) {
               	log.denominator.2.n=sum(current.z.2.n*log(current.prob.mat.2.n)+(1-current.z.2.n)*log(1-current.prob.mat.2.n))
		proposed.e=rnorm(1,current.e,sd.prop.e)
		if(proposed.e<0) proposed.e=-proposed.e
		if(proposed.e>1) proposed.e=2-proposed.e
		proposed.p.extinct=proposed.e/current.site.area.x
		proposed.p.extinct[proposed.p.extinct>1]=1
		proposed.p.extinct.re.2.n=proposed.p.extinct*(1-current.p.colon.2.n)
		proposed.prob.mat.2.n=current.z.1.nminus1*(1-proposed.p.extinct.re.2.n)+(1-current.z.1.nminus1)*current.p.colon.2.n

		log.numerator.2.n=sum(current.z.2.n*log(proposed.prob.mat.2.n)+(1-current.z.2.n)*log(1-proposed.prob.mat.2.n))
		log.MH=log.numerator.2.n-log.denominator.2.n
		accept=decide(log.MH)
		if (accept) {
			current.e=proposed.e
			current.p.extinct=proposed.p.extinct
			current.prob.mat.2.n=proposed.prob.mat.2.n
                        log.denominator.2.n=log.numerator.2.n
			}
              }
		#-----------------
		#-----------x
		#-----------------
          if(sd.prop.x>0) {
          proposed.x=rnorm(1,current.x,sd=sd.prop.x)
		if(proposed.x<0) proposed.x=-proposed.x
		if(proposed.x>5) proposed.x=10-proposed.x
		proposed.site.area.x=site.area^proposed.x
		proposed.p.extinct=current.e/proposed.site.area.x
		proposed.p.extinct[proposed.p.extinct>1]=1
		proposed.p.extinct.re.2.n=proposed.p.extinct*(1-current.p.colon.2.n)
		proposed.prob.mat.2.n=current.z.1.nminus1*(1-proposed.p.extinct.re.2.n)+(1-current.z.1.nminus1)*current.p.colon.2.n

		log.numerator.2.n=sum(current.z.2.n*log(proposed.prob.mat.2.n)+(1-current.z.2.n)*log(1-proposed.prob.mat.2.n))
		log.MH=log.numerator.2.n-log.denominator.2.n
		accept=decide(log.MH)
		if (accept) {
                  	current.x=proposed.x
			current.site.area.x=proposed.site.area.x
			current.p.extinct=proposed.p.extinct
			current.prob.mat.2.n=proposed.prob.mat.2.n
                        log.denominator.2.n=log.numerator.2.n
			}
        }
		#-----------------
		#-----------y
		#-----------------
          if(sd.prop.y>0) {
		proposed.y=rnorm(1,current.y,sd.prop.y)
		if(proposed.y<0) proposed.y=-proposed.y
		if(proposed.y>400) proposed.y=800-proposed.y
		proposed.p.colon.2.n=current.s.i.sq.1.nminus1/(current.s.i.sq.1.nminus1+proposed.y)
		proposed.p.extinct.re.2.n=current.p.extinct*(1-proposed.p.colon.2.n)
		proposed.prob.mat.2.n=current.z.1.nminus1*(1-proposed.p.extinct.re.2.n)+(1-current.z.1.nminus1)*proposed.p.colon.2.n
		log.numerator.2.n=sum(current.z.2.n*log(proposed.prob.mat.2.n)+(1-current.z.2.n)*log(1-proposed.prob.mat.2.n))
		log.MH=log.numerator.2.n-log.denominator.2.n
		accept=decide(log.MH)
		if(accept) {
			current.y=proposed.y
                        current.p.colon.2.n=proposed.p.colon.2.n
			current.prob.mat.2.n=proposed.prob.mat.2.n
                       	log.denominator.2.n=log.numerator.2.n
			}

              }

  #-----------------
		#-----------b
		#-----------------
          if(sd.prop.b>0) {
		proposed.b=rnorm(1,current.b,sd.prop.b)
		if(proposed.b<0) proposed.b=-proposed.b
                #ASSUMED THAT VALUES OF b GREATER THAN 1 DO NOT MAKE SENSE OVER THE
                #RANGE OF HA IN OUR STUDY.
		if(proposed.b>5) proposed.b=10-proposed.b
		proposed.area.j.b=site.area^proposed.b
		proposed.area.occ.1.nminus1=current.z.1.nminus1*proposed.area.j.b
		proposed.s.i.sq.1.nminus1=(current.sim.alpha.distance%*%proposed.area.occ.1.nminus1)^2
		proposed.p.colon.2.n=proposed.s.i.sq.1.nminus1/(proposed.s.i.sq.1.nminus1+current.y)
		proposed.p.extinct.re.2.n=current.p.extinct*(1-proposed.p.colon.2.n)
		proposed.prob.mat.2.n=current.z.1.nminus1*(1-proposed.p.extinct.re.2.n)+(1-current.z.1.nminus1)*proposed.p.colon.2.n
          if(sum(proposed.prob.mat.2.n>=1)>0){
			proposed.prob.mat.2.n[proposed.prob.mat.2.n>=1]=one.m[proposed.prob.mat.2.n>=1]
			}
		log.numerator.2.n=sum(current.z.2.n*log(proposed.prob.mat.2.n)+(1-current.z.2.n)*log(1-proposed.prob.mat.2.n))
		log.MH=log.numerator.2.n-log.denominator.2.n

		accept=decide(log.MH)
		if(accept) {
			current.b=proposed.b
			current.area.j.b=proposed.area.j.b
			current.area.occ.1.nminus1=proposed.area.occ.1.nminus1
                        current.s.i.sq.1.nminus1=proposed.s.i.sq.1.nminus1
                        current.p.colon.2.n=proposed.p.colon.2.n
                       	current.prob.mat.2.n=proposed.prob.mat.2.n
			log.denominator.2.n=log.numerator.2.n
			}
              }


		#-----------------
		#-----------ALPHA
		#-----------------
          if(sd.prop.alpha>0) {
		proposed.alpha=rnorm(1,current.alpha,sd.prop.alpha)
		if(proposed.alpha<1) proposed.alpha=2-proposed.alpha
		if(proposed.alpha>30) proposed.alpha=60-proposed.alpha
		proposed.sim.alpha.distance=exp(-proposed.alpha*site.distance)
		diag(proposed.sim.alpha.distance)=0
		proposed.s.i.sq.1.nminus1=(proposed.sim.alpha.distance%*%current.area.occ.1.nminus1)^2
		proposed.p.colon.2.n=proposed.s.i.sq.1.nminus1/(proposed.s.i.sq.1.nminus1+current.y)
		proposed.p.extinct.re.2.n=current.p.extinct*(1-proposed.p.colon.2.n)
		proposed.prob.mat.2.n=current.z.1.nminus1*(1-proposed.p.extinct.re.2.n)+(1-current.z.1.nminus1)*proposed.p.colon.2.n
				if(sum(proposed.prob.mat.2.n<=0)>0){
				print(proposed.prob.mat.2.n)
				proposed.prob.mat.2.n[proposed.prob.mat.2.n<=0]=zero.m[proposed.prob.mat.2.n<=0]
				}
		log.numerator.2.n=sum(current.z.2.n*log(proposed.prob.mat.2.n)+(1-current.z.2.n)*log(1-proposed.prob.mat.2.n))
		log.MH=log.numerator.2.n-log.denominator.2.n
		accept=decide(log.MH)
		if(accept) {
			current.alpha=proposed.alpha
			current.sim.alpha.distance=proposed.sim.alpha.distance
			current.s.i.sq.1.nminus1=proposed.s.i.sq.1.nminus1
			current.p.colon.2.n=proposed.p.colon.2.n
			current.prob.mat.2.n=proposed.prob.mat.2.n
                        log.denominator.2.n=log.numerator.2.n
			}
        }

	if(n==iter.chain[n]) {
          z.chain=c(z.chain,current.z)
                p.chain=matrix(c(p.chain,current.p),nrow=nyear)
                mupsi1.chain=c(mupsi1.chain,current.mupsi1)
		e.chain=c(e.chain,current.e)
		x.chain=c(x.chain,current.x)
		y.chain=c(y.chain,current.y)
		b.chain=c(b.chain,current.b)
		alpha.chain=c(alpha.chain,current.alpha)
		latent.deviance.chain = c(latent.deviance.chain,-2*log.denominator.2.n)
		}
                if(n==iter.print[n])  print(n)
		}

	#####################################################
        dim(z.chain)=c(nsite,nyear,length(alpha.chain))

        muz.chain=NULL
        muz.missing.chain=NULL
        n.extinct.chain=NULL
        n.colon.chain=NULL
        n.extinct.obs.chain=NULL
        n.colon.obs.chain=NULL
        n.occ.pairs.obs.chain=NULL

        for(t in 1:nyear) {
          muz.missing.chain=c(muz.missing.chain,apply(z.chain[(z.missing.index.list[[t]]-(t-1)*nsite),t,],2,mean))
          muz.chain=c(muz.chain,apply(z.chain[,t,],2,mean))
        }
        for(t in 1:(nyear-1)){
          n.extinct.mat.temp=(z.chain[,t,]==1 & z.chain[,(t+1),]==0)
          n.extinct.chain=c(n.extinct.chain,apply(n.extinct.mat.temp,2,sum))
          n.extinct.obs.chain=c(n.extinct.obs.chain,apply(n.extinct.mat.temp[z.index.list[[t]],],2,sum))

          n.colon.mat.temp=(z.chain[,t,]==0 & z.chain[,(t+1),]==1)
          n.colon.chain=c(n.colon.chain,apply(n.colon.mat.temp,2,sum))
          n.colon.obs.chain=c(n.colon.obs.chain,apply(n.colon.mat.temp[z.index.list[[t]],],2,sum))

          n.occ.pairs.obs.temp=apply(z.chain[z.index.list[[t]],t,],2,sum)
          n.occ.pairs.obs.chain=c(n.occ.pairs.obs.chain,n.occ.pairs.obs.temp)
        }
        muz.missing.chain=matrix(muz.missing.chain,nrow=nyear,byrow=TRUE)
        muz.chain=matrix(muz.chain,nrow=nyear,byrow=TRUE)

        n.extinct.chain=matrix(n.extinct.chain,nrow=(nyear-1),byrow=TRUE)
        prop.extinct.chain=n.extinct.chain/(nsite*muz.chain[1:(nyear-1),])

        n.colon.chain=matrix(n.colon.chain,nrow=(nyear-1),byrow=TRUE)
        prop.colon.chain=n.colon.chain/(nsite*(1-muz.chain[1:(nyear-1),]))

        n.extinct.obs.chain=matrix(n.extinct.obs.chain,nrow=(nyear-1),byrow=TRUE)
        n.colon.obs.chain=matrix(n.colon.obs.chain,nrow=(nyear-1),byrow=TRUE)
        n.occ.pairs.obs.chain=matrix(n.occ.pairs.obs.chain,nrow=(nyear-1),byrow=TRUE)

        prop.extinct.obs.chain=n.extinct.obs.chain/n.occ.pairs.obs.chain
        prop.colon.obs.chain=n.colon.obs.chain/(n.pairs.obs-n.occ.pairs.obs.chain)

        prop.extinct.missing.chain=(n.extinct.chain-n.extinct.obs.chain)/(nsite*muz.chain[1:(nyear-1),]-n.occ.pairs.obs.chain)
        prop.colon.missing.chain=(n.colon.chain-n.colon.obs.chain)/(nsite*(1-muz.chain[1:(nyear-1),])-(n.pairs.obs-n.occ.pairs.obs.chain))


## NOTE: code below changes y back to original scale:
		#return(list(z.chain=z.chain,muz.chain=muz.chain,muz.missing.chain=muz.missing.chain,prop.extinct.chain=prop.extinct.chain,prop.extinct.obs.chain=prop.extinct.obs.chain,prop.extinct.missing.chain=prop.extinct.missing.chain,prop.colon.chain=prop.colon.chain,prop.colon.obs.chain=prop.colon.obs.chain,prop.colon.missing.chain=prop.colon.missing.chain,p.chain=p.chain,mupsi1.chain=mupsi1.chain,e.chain=e.chain,x.chain=x.chain,y.chain=sqrt(y.chain),b.chain=b.chain,alpha.chain=alpha.chain,latent.deviance.chain=latent.deviance.chain))

        return(list(z.chain=z.chain,muz.chain=muz.chain,muz.missing.chain=muz.missing.chain,prop.extinct.chain=prop.extinct.chain,prop.colon.chain=prop.colon.chain,p.chain=p.chain,mupsi1.chain=mupsi1.chain,e.chain=e.chain,x.chain=x.chain,y.chain=sqrt(y.chain),b.chain=b.chain,alpha.chain=alpha.chain,latent.deviance.chain=latent.deviance.chain))
}
