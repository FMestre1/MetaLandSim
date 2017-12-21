\name{ifm.naive.MCMC}
\alias{ifm.naive.MCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the naive design incidence function model
}
\description{
Estimates the IFM assuming no false absences and omitting sites for particular years in which data were missing.
}
\usage{
ifm.naive.MCMC(niter=1000,init,z.data, site.distance, site.area, sd.prop.e=0.2,
 sd.prop.x=0.5,sd.prop.y=10, sd.prop.b=0.2, sd.prop.alpha=5,nthin=1,print.by=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{niter}{
  Number of iterations in the MCMC chain.
}
  \item{init}{
  %- this is the 'init' from the robust function is it the same here?
Named list with values to initialize the chain. E.g.:\cr \cr
 init1=list(alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),\cr \cr
 e=runif(1,0,1),x=runif(1,0,5)).\cr \cr
alpha: initial value for alpha in dispersal model; described as 1 / average dispersal distance \cr \cr
b: initial value for parameter b in colonization model \cr \cr
y: initial value for parameter y in colonization model \cr \cr
e: initial value for e in extinction model \cr \cr
x: initial value for x in extinction model
}
  \item{z.data}{
  nsite x nyears matrix. If contains NAs, the corresponding parts are omitted from the likelihood (the missing data are not estimated).
}
  \item{site.distance}{
nsite x nsite matrix of distances between sites. The tuning parameters in the example are set for distances less than one, with max distance approximately 0.5.
}
  \item{site.area}{
Vector of length nsite with areas. The tuning parameters in the example are set for average area approximately equal to 1.
}
  \item{sd.prop.e}{
Standard deviation of the proposal distribution for parameter e.
}
  \item{sd.prop.x}{
Standard deviation of the proposal distribution for parameter x.
}
  \item{sd.prop.y}{
Standard deviation of the proposal distribution for parameter y.
}
  \item{sd.prop.b}{
Standard deviation of the proposal distribution for parameter b.
}
  \item{sd.prop.alpha}{
Standard deviation of the proposal distribution for parameter alpha.
}
  \item{nthin}{
If specified, keeps only every nthin^th sample from the MCMC chain. Use to save memory or when the chain is moving slowly.
}
  \item{print.by}{
Specifies how often to print the number of the current iteration.}
}
\value{
\item{e.chain}{posterior sample of e}
\item{x.chain}{posterior sampmle of x}
\item{y.chain}{posterior sample of y}
\item{b.chain}{posterior sample of b}
\item{alpha.chain}{posterior sample of alpha}
\item{deviance.chain}{posterior sample of -2*loglik}
%% ...
}
\references{
Risk, B. B., De Valpine, P., Beissinger, S. R. (2011). A robust design formulation of the incidence function model of metapopulation dynamics applied to two species of rails. Ecology, 92(2), 462-474.
}
\author{
Benjamin Risk
}
\examples{
\dontrun{
data(simulatedifm)

library("coda")

myniter=5000
nsite=nrow(z.sim)
nyear=ncol(z.sim)
nthin=1
nburnin=1000
## NOTE! The notation used here corresponds to MetaLandSim and differs from Risk et al 2011
## Here
## e (in MetaLandSim) = mu
## x = chi
## y = gamma
## b = beta
## alpha = alpha
##
# Priors:
#         e: [0,1]
#         x: [0,5]
#         y^2: [0,400]
#         b: [0,5]
#         alpha: [1,30]

# NOTE: If posteriors are truncated at zero, then estimates are biased. Rescale
# distances (e.g., divide by 10,000) and/or areas so that parameters are larger.

# Here, we run two chains with random initial values:
init1=list(alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1),x=runif(1,0,5))

a = Sys.time()
inm1 <- ifm.naive.MCMC(niter=myniter,init=init1,z.data =
 z.sim,site.distance=sim.distance,site.area=sim.area,
  sd.prop.alpha=4,sd.prop.b=0.6,sd.prop.y=40,sd.prop.e=0.05,sd.prop.x=0.4,nthin=1,print.by=1000)
accept.calculate(inm1,model='naive')
Sys.time() - a

init2=list(alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1),x=runif(1,0,5))
inm2 <- ifm.naive.MCMC(niter=myniter,init=init2,z.data =
z.sim,site.distance=sim.distance,site.area=sim.area,
sd.prop.alpha=4,sd.prop.b=0.6,sd.prop.y=40,sd.prop.e=0.05,sd.prop.x=0.4,nthin=1,print.by=1000)
accept.calculate(inm2,model='naive')
Sys.time() - a

coda.create(inm1,"sim_inm1",par.list=list("e.chain","x.chain","alpha.chain",
"b.chain","y.chain"),niter=myniter,nthin=nthin)
coda.create(inm2,"sim_inm2",par.list=list("e.chain","x.chain","alpha.chain",
"b.chain","y.chain"),niter=myniter,nthin=nthin)
coda.sim.inm1=read.coda("sim_inm1.txt","sim_inm1_Index.txt")
coda.sim.inm2=read.coda("sim_inm2.txt","sim_inm2_Index.txt")
coda.sim.inm.list=mcmc.list(coda.sim.inm1,coda.sim.inm2)
sim.inm=combine.chains(inm1,inm2,nburnin=nburnin,nthin=1)
coda.create(sim.inm,"sim_inm",par.list=list("e.chain","x.chain","alpha.chain",
"b.chain","y.chain"),niter=(2*myniter-2*nburnin),nthin=nthin)
coda.sim.inm.long=read.coda("sim_inm.txt","sim_inm_Index.txt")

summary(coda.sim.inm.list)
summary(coda.sim.inm.long)

gelman.diag(coda.sim.inm.list)

plot(coda.sim.inm.list)
plot(coda.sim.inm.long)
cumuplot(coda.sim.inm.long)

# calculate maximum a posteriori estimates:
m1 <- as.matrix(sim.inm)
e <- calcmode(m1[,1][[1]])
x <- calcmode(m1[,1][[2]])
y <- calcmode(m1[,1][[3]])
b <- calcmode(m1[,1][[4]])
alpha <- calcmode(m1[,1][[5]])

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ifm}
\keyword{metapopulation}
