\name{coda.create}
\alias{coda.create}
\title{
Create files for use with R-package coda.
}
\description{
Creates two text files, <filename.txt> and filename_Index.txt>, in the format used by OpenBUGS, which can then be read using read.coda() to create an mcmc object for subsequent use of coda diagnostic and plotting functions.
}
\usage{
coda.create(object, file.name, write.index = TRUE, par.list = list("mupsi1.chain", 
"e.chain", "x.chain", "b.chain", "y.chain", "alpha.chain"), niter = 101000, nthin = 10)
}
\arguments{
  \item{object}{
Output from ifm.naive.MCMC, ifm.missing.MCMC, or ifm.robust.MCMC, i.e., a named list with vectors e.chain. x.chain, y.chain, b.chain, alpha.chain, and other parameters depending on the model.
}
  \item{file.name}{
Name of the text file that will be output in the format used by BUGS and JAGS.
}
  \item{write.index}{
Logical indicating whether or not to also create file.name_Index.txt. Defaults to TRUE, which allows the subsequent use of read.coda(). 
}
  \item{par.list}{
List of parameters to include in the file. Defaults to mupsi1.chain, e.chain, x.chain, b.chain, y.chain, alpha.chain; modify to include additional parameters in ifm.missing.MCMC and ifm.robust.MCMC.}
  \item{niter}{
Number of iterations in the original chain (before thinning).
}
  \item{nthin}{
Thinning used in the original estimation function; this is just bookkeeping, as this function does not perform the thinning.
}
}
\details{
Writes text files to the current working directory or to the path specified with "file.name". 
}
\author{
Benjamin Risk
}
\examples{

\dontrun{

# quick run; actual estimation requires more iterations:

data(simulatedifm)

myniter=100
nsite=nrow(z.sim)
nyear=ncol(z.sim)
nthin=1
nburnin=0
init1=list(alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1),x=runif(1,0,5))

inm1 <- ifm.naive.MCMC(niter=myniter,init=init1,z.data =
 z.sim,site.distance=sim.distance,site.area=sim.area,
  sd.prop.alpha=4,sd.prop.b=0.6,sd.prop.y=40,sd.prop.e=0.05,sd.prop.x=0.4,nthin=1,print.by=1000)

# write files in OpenBUGS format to working directory:
coda.create(inm1,"sim_inm1",par.list=list("e.chain","x.chain","alpha.chain",
"b.chain","y.chain"),niter=myniter,nthin=nthin)
}
}
