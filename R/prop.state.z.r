prop.state.z<-function(z.subset,year.vector,current.p,current.mupsi1,current.prob.mat.2.n,current.z,current.area.j.b,current.sim.alpha.distance,current.y,current.p.extinct,nsucc,ntrial,nsite,nyear) {

  year.temp=max(year.vector[z.subset])
  if(year.temp==1) {
    proposed.z.vector=rbinom(length(z.subset),1,current.mupsi1)
  }
  else {
    proposed.z.vector=rbinom(length(z.subset),1,current.prob.mat.2.n[z.subset-nsite])
  }
  proposed.z.vector[nsucc[z.subset]>=1]=1
  log.lik.data.proposed=log.lik.data.z.p(ntrial.temp=ntrial[z.subset],nsucc.temp=nsucc[z.subset],z.temp=proposed.z.vector,
    p.temp=current.p[year.temp])
  log.lik.data.current=log.lik.data.z.p(ntrial.temp=ntrial[z.subset],nsucc.temp=nsucc[z.subset],z.temp=current.z[z.subset],
    p.temp=current.p[year.temp])
  if(year.temp==nyear) {
    log.MH=log.lik.data.proposed-log.lik.data.current
  }
  else {
  proposed.z=current.z[,year.temp]
  proposed.z[(z.subset-(year.temp-1)*nsite)]=proposed.z.vector
  proposed.area.occ=proposed.z*current.area.j.b
  proposed.s.i.sq=(current.sim.alpha.distance%*%proposed.area.occ)^2
  proposed.p.colon=proposed.s.i.sq/(proposed.s.i.sq+current.y)
  proposed.p.extinct.re=current.p.extinct*(1-proposed.p.colon)

  proposed.prob=proposed.z*(1-proposed.p.extinct.re)+(1-proposed.z)*proposed.p.colon
  log.prob.z.tplus1.proposed=sum(current.z[,(year.temp+1)]*log(proposed.prob)+(1-current.z[,(year.temp+1)])*log(1-proposed.prob))
  log.prob.z.tplus1.current=sum(current.z[,(year.temp+1)]*log(current.prob.mat.2.n[,year.temp])+
				(1-current.z[,(year.temp+1)])*log(1-current.prob.mat.2.n[,year.temp]))
  log.MH=log.lik.data.proposed+log.prob.z.tplus1.proposed-log.lik.data.current-log.prob.z.tplus1.current
  }
  accept=decide(log.MH)
  if (accept) {
    current.z[z.subset]=proposed.z.vector
    if(year.temp!=nyear) {
    current.prob.mat.2.n[,year.temp]=proposed.prob
  }
    return(list(accept=accept,current.z=current.z, current.prob.mat.2.n=current.prob.mat.2.n))
  }
  else{
    return(list(accept=accept))
  }
}
