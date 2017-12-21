prop.year.missing <- function(z.missing.index.subset,year.vector,current.mupsi1,current.prob.mat.2.n,
  current.z,current.area.j.b,current.sim.alpha.distance,current.y,current.p.extinct,nsite) {
  year.missing=year.vector[min(z.missing.index.subset)]
  if(year.missing==1) {
    proposed.z.missing.vector=rbinom(length(z.missing.index.subset),1,current.mupsi1)
  }
  else {
    proposed.z.missing.vector=rbinom(length(z.missing.index.subset),1,current.prob.mat.2.n[z.missing.index.subset-nsite])
    }
  proposed.z=current.z[,year.missing]
  proposed.z[(z.missing.index.subset-(year.missing-1)*nsite)]=proposed.z.missing.vector
  proposed.area.occ=proposed.z*current.area.j.b
  proposed.s.i.sq=(current.sim.alpha.distance%*%proposed.area.occ)^2
  proposed.p.colon=proposed.s.i.sq/(proposed.s.i.sq+current.y)
  proposed.p.extinct.re=current.p.extinct*(1-proposed.p.colon) 
  proposed.prob=proposed.z*(1-proposed.p.extinct.re)+(1-proposed.z)*proposed.p.colon	

  log.prob.data.proposed=sum(current.z[,(year.missing+1)]*log(proposed.prob)+(1-current.z[,(year.missing+1)])*log(1-proposed.prob))
  log.prob.data.current=sum(current.z[,(year.missing+1)]*log(current.prob.mat.2.n[,year.missing])+
    (1-current.z[,(year.missing+1)])*log(1-current.prob.mat.2.n[,year.missing]))
    log.MH=log.prob.data.proposed-log.prob.data.current
    accept=decide(log.MH)
  if(accept) {
    current.z[z.missing.index.subset]=proposed.z.missing.vector
    current.prob.mat.2.n[,year.missing]=proposed.prob
    return(list(accept=accept,current.z=current.z,current.prob.mat.2.n=current.prob.mat.2.n))
  }
  else return(list(accept=accept))
}
