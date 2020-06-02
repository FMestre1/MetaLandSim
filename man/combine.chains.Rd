\name{combine.chains}
\alias{combine.chains}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combines two chains into a single chain.
}
\description{
Combines two lists of chains from ifm.naive.MCMC, ifm.missing.MCMC, or ifm.robust.MCMC into one list where each element is the concatenated chains.
}
\usage{
combine.chains(x1, x2, nburnin, nthin = 1, z.thin = TRUE)
}
\arguments{
  \item{x1}{
  First list of chains.
}
  \item{x2}{
 Second list of chains.
 }
  \item{nburnin}{
Number of initial iterations to discard.
}
  \item{nthin}{
If nthin>1, subsets to every nthin^th sample
}
  \item{z.thin}{
logical; defaults to TRUE. Thinning for the posterior sample of the occupancy states. If true, uses thinning equal to 5. The posterior sample of occupancy states is a large nsite x nyear x niter array, and this option reduces memory usage. Ignored if the chain is from the ifm.naive.MCMC (where occupancy states are fixed).
}
}
\value{
Named list with the same names as the inputs x1 and x2
}
\author{
Benjamin Risk
}
\examples{
data(simulatedifm)

init1=list(alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1),x=runif(1,0,5))

inm1 <- ifm.naive.MCMC(niter=500,init=init1,z.data =
 z.sim,site.distance=sim.distance,site.area=sim.area,
  sd.prop.alpha=4,sd.prop.b=0.6,sd.prop.y=40,sd.prop.e=0.05,sd.prop.x=0.4,nthin=1,print.by=100)
  
init1=list(alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1),x=runif(1,0,5))

inm2 <- ifm.naive.MCMC(niter=500,init=init1,z.data =
 z.sim,site.distance=sim.distance,site.area=sim.area,
  sd.prop.alpha=4,sd.prop.b=0.6,sd.prop.y=40,sd.prop.e=0.05,sd.prop.x=0.4,nthin=1,print.by=100)

sim.inm=combine.chains(inm1,inm2,nburnin=0,nthin=1)
}