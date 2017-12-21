\name{ifm.robust.MCMC}
\alias{ifm.robust.MCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the robust design incidence function model
}
\description{
Estimates the IFM with imperfect detection and missing data.
}
\usage{
ifm.robust.MCMC(niter = 1000, init, det.data, site.distance, site.area, sd.prop.p = 0.1,
sd.prop.mupsi1 = 0.1, sd.prop.e = 0.2, sd.prop.x = 0.2, sd.prop.y = 0.2, sd.prop.b = 0.2,
sd.prop.alpha = 0.2, nthin = 1, nsite.subset = 5, print.by = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{niter}{
  Number of iterations in the MCMC chain.
}
  \item{init}{
Named list with values to initialize the chain. E.g.:\cr \cr
 init1=list(z.data=initocc,z.missing=runif(nmissing),p=runif(nyear,\cr \cr
0.1,1),mupsi1=runif(1),alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),\cr \cr
e=runif(1,0,1),x=runif(1,0,5)).\cr \cr
z.data: a matrix with nrows = number of sites and ncol = number of years. Contains NAs for missing values. Contains naive estimates of occupancy elsewhere.\cr \cr
z.missing: z.missing: a vector of initial occupancy states for the missing data with length equal to the number of NAs in z.data (i.e., vectorized across years). Can use runif(nmissing). \cr \cr
p: vector of length nyears with inital probability of detection in each year \cr \cr
mupsi1: probability of initial occupancy in year 1; runif(1) suffices \cr \cr
alpha: initial value for alpha in dispersal model; described as 1 / average dispersal distance \cr \cr
b: initial value for parameter b in colonization model \cr \cr
y: initial value for parameter y in colonization model \cr \cr
e: initial value for e in extinction model \cr \cr
x: initial value for x in extinction model
}
  \item{det.data}{
Detection data in an array with dimensions nsites x nyears x nvisits. For removal design, set all values after a detection equal to NA. For missing data in a given year, set all visits to NA.
}
  \item{site.distance}{
nsite x nsite matrix of distances between sites. The tuning parameters in the example are set for distances less than one, with max distance approximately 0.5. Input data should have a similar scaling.
}
  \item{site.area}{
Vector of length nsite with areas. The tuning parameters in the example are set for average area approximately equal to 1. Input data should have a similar scaling.
}
  \item{sd.prop.p}{
Scalar equal to the standard deviation of the proposal distribution for probability of detection, which is a normal distribution centered at current value in the mcmc chain. The same standard deviation is used for all years.
}
  \item{sd.prop.mupsi1}{
Standard deviation of the proposal distribution for occupancy in year 1.
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
  \item{nsite.subset}{
The number of sites to include in the block sampling, where nsite.subset is equal to the number of sites updated in the same step. Larger values decrease the probability of acceptance.
}
  \item{print.by}{
Specifies how often to print the number of the current iteration.}
}
\value{
\item{z.chain}{nsite x nyear x niter array sampled from the posterior distribution of occupancy in each year (if detection occurred at a given year and site, then the value is identically equal to one for all iterations).}
\item{muz.chain}{nyear x niter matrix posterior sample of the proportion of sites occupied in each year.}
\item{muz.missing.chain}{nyear x niter matrix posterior sample of the proportion of sites occupied for sites with missing data.}
\item{prop.extinct.chain}{Extinction rate for all sites.}
\item{prop.colon.chain}{Colonization rate.}
\item{p.chain}{nyear x niter sample of detection probabilities.}
\item{mupsi1.chain}{posterior sample of parameter for occupancy in year 1.}
\item{e.chain}{posterior sample of e}
\item{x.chain}{posterior sampmle of x}
\item{y.chain}{posterior sample of y}
\item{b.chain}{posterior sample of b}
\item{alpha.chain}{posterior sample of alpha}
\item{latent.deviance.chain}{posterior sample of -2*loglik}
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

# There are more parameters in this model
# and estimating the posterior requires more iterations:
niter=2000
nsite=nrow(z.sim)
nyear=ncol(z.sim)
nthin=1
nburnin=1000
## NOTE! The notation used here corresponds to MetaLandSim and differs from Risk et al 2011
## Here
## e (in MetaLandSim) = mu (in Risk et al 2011)
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

# Count number of times a site was never visited in a given year:
nmissing = sum(is.na(z.sim.20))

# Create a dataset with initial guess of true occupancy for sites with visits.
# This dataset should be number of sites by years
# one way of generating these initial values:
initocc <- suppressWarnings(apply(sim.det.20,c(1,2),max,na.rm=TRUE))
# produces warnings but that's okay
initocc[initocc==-Inf]=NA


init1=list(z.data=initocc,z.missing=runif(nmissing),p=runif(nyear,
0.1,1),mupsi1=runif(1),alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1)
,x=runif(1,0,5))

# for diagnosing acceptance rates:
# init1=list(z.data=initocc,z.missing=runif(nmissing),p=runif(nyear,0.1,1),
mupsi1=runif(1),alpha=20, b=0.5,y=7.5,e=0.25,x=0.25)

a = Sys.time()
ir1 <- ifm.robust.MCMC(niter=niter,init=init1, det.data = sim.det.20,
site.distance=sim.distance,site.area=sim.area, sd.prop.p=0.25,sd.prop.mupsi1=0.2, 
sd.prop.alpha=2, sd.prop.b=0.6, sd.prop.y=20, sd.prop.e=0.1, sd.prop.x=0.4, nthin=1)
accept.calculate(ir1,model='robust')
Sys.time() - a

init2=list(z.data=initocc,z.missing=runif(nmissing),p=runif(nyear,
0.1,1),mupsi1=runif(1),alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),
e=runif(1,0,1),x=runif(1,0,5))
ir2 <- ifm.robust.MCMC(niter=niter,init=init2, det.data = sim.det.20,
 site.distance=sim.distance,site.area=sim.area, sd.prop.mupsi1=0.2, sd.prop.alpha=2, sd.prop.b=0.6,
  sd.prop.y=20, sd.prop.e=0.1, sd.prop.x=0.4, sd.prop.p = 0.25, nthin=1)
accept.calculate(ir2,model='robust')

coda.create(ir1,"sim_ir1",par.list=list("mupsi1.chain","e.chain","x.chain",
"alpha.chain","b.chain","y.chain","p.chain"),niter=niter,nthin=nthin)
coda.create(ir2,"sim_ir2",par.list=list("mupsi1.chain","e.chain","x.chain",
"alpha.chain","b.chain","y.chain","p.chain"),niter=niter,nthin=nthin)
coda.sim.ir1=read.coda("sim_ir1.txt","sim_ir1_Index.txt")
coda.sim.ir2=read.coda("sim_ir2.txt","sim_ir2_Index.txt")
coda.sim.ir.list=mcmc.list(coda.sim.ir1,coda.sim.ir2)
sim.ir=combine.chains(ir1,ir2,nburnin=nburnin,nthin=1)
coda.create(sim.ir,"sim_ir",par.list=list("mupsi1.chain","e.chain",
"x.chain","alpha.chain","b.chain","y.chain","p.chain"),niter=(2*niter-2*nburnin),nthin=nthin)
coda.sim.ir.long=read.coda("sim_ir.txt","sim_ir_Index.txt")

summary(coda.sim.ir.list)
summary(coda.sim.ir.long)

gelman.diag(coda.sim.ir.list)

plot(coda.sim.ir.list)
plot(coda.sim.ir.long)
cumuplot(coda.sim.ir.long)

# calculate maximum a posteriori estimates:
m1 <- as.matrix(sim.ir)
e <- calcmode(m1[,1][[1]])
x <- calcmode(m1[,1][[2]])
y <- calcmode(m1[,1][[3]])
b <- calcmode(m1[,1][[4]])
alpha <- calcmode(m1[,1][[5]])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{robust}% use one of  RShowDoc("KEYWORDS")
\keyword{missing}% __ONLY ONE__ keyword per line
\keyword{occupancy}
\keyword{metapopulation}
\keyword{ifm}
