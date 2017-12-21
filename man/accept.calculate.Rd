\name{accept.calculate}
\alias{accept.calculate}
\title{
Calculate acceptance rates in MCMC chains
}
\description{
Calculate acceptance rates of parameters in the IFM.
}
\usage{
accept.calculate(x, model = c("naive", "missing", "robust"))
}
\arguments{
  \item{x}{
  A named list with the MCMC chains estimated by ifm.naive.MCMC, ifm.missing.MCMC, or ifm.robust.MCMC.
  }
    \item{model}{
  Either "naive", "missing", or "robust"
}
}
\value{
Named list containing MCMC chain acceptance rates. Names are built from the input list, e.g., for model=``naive":
\item{acc.b.chain}{Acceptance rates of parameter b}
\item{acc.e.chain}{Acceptance rates of parameter e}
\item{acc.y.chain}{Acceptance rates of parameter y}
\item{acc.alpha.chain}{Acceptance rates of parameter alpha}
\item{acc.x.chain}{Acceptance rates of parameter x}
}
\author{
Benjamin Risk
}
\examples{

data(simulatedifm)

# Here, we run a chain with random initial values:
init1=list(alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1),x=runif(1,0,5))

inm1 <- ifm.naive.MCMC(niter=1000,init=init1,z.data =
 z.sim,site.distance=sim.distance,site.area=sim.area,
  sd.prop.alpha=4,sd.prop.b=0.6,sd.prop.y=40,sd.prop.e=0.05,sd.prop.x=0.4,nthin=1,print.by=100)
accept.calculate(inm1,model='naive')

}
