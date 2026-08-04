# Estimate the 'missing' design incidence function model

Estimates the IFM with no false absences but incorporating missing data.

## Usage

``` r
ifm.missing.MCMC(niter=1000,init,z.data, site.distance, site.area,
sd.prop.mupsi1=0.1, sd.prop.e=0.1, sd.prop.x=0.5,sd.prop.y=10, sd.prop.b=0.2,
sd.prop.alpha=5, nthin=1,nsite.subset=10,print.by=100)
```

## Arguments

- niter:

  Number of iterations in the MCMC chain.

- init:

  Named list with values to initialize the chain. E.g.:  
    
  init1=list(z.missing=runif(nmissing),mupsi1=runif(1),alpha=runif(1,1,30),  
  b=runif(1,0,5),y=runif(1,0,20),  
    
  e=runif(1,0,1),x=runif(1,0,5)).  
    
  z.missing: a vector of initial occupancy states for the missing data
  with length equal to the number of NAs in z.data (i.e., vectorized
  across years). Can use runif(nmissing).  
    
  mupsi1: probability of initial occupancy in year 1; runif(1)
  suffices  
    
  alpha: initial value for alpha in dispersal model; described as 1 /
  average dispersal distance  
    
  b: initial value for parameter b in colonization model  
    
  y: initial value for parameter y in colonization model  
    
  e: initial value for e in extinction model  
    
  x: initial value for x in extinction model

- z.data:

  nsite x nyears matrix containing NA for missing data. Occupancy at
  sites with missing data will be estimated.

- site.distance:

  nsite x nsite matrix of distances between sites. The tuning parameters
  in the example are set for distances less than one, with max distance
  approximately 0.5. Input data should have a similar scaling.

- site.area:

  Vector of length nsite with areas. The tuning parameters in the
  example are set for average area approximately equal to 1. Input data
  should have a similar scaling.

- sd.prop.mupsi1:

  Standard deviation of the proposal distribution for occupancy in year
  1.

- sd.prop.e:

  Standard deviation of the proposal distribution for parameter e.

- sd.prop.x:

  Standard deviation of the proposal distribution for parameter x.

- sd.prop.y:

  Standard deviation of the proposal distribution for parameter y.

- sd.prop.b:

  Standard deviation of the proposal distribution for parameter b.

- sd.prop.alpha:

  Standard deviation of the proposal distribution for parameter alpha.

- nthin:

  If specified, keeps only every nthin^th sample from the MCMC chain.
  Use to save memory or when the chain is moving slowly.

- nsite.subset:

  The number of sites to include in the block sampling, where
  nsite.subset is equal to the number of sites updated in the same step.
  Larger values decrease the probability of acceptance.

- print.by:

  Specifies how often to print the number of the current iteration.

## Value

- z.chain:

  nsite x nyear x niter array sampled from the posterior distribution of
  occupancy in each year (if detection occurred at a given year and
  site, then the value is identically equal to one for all iterations).

- muz.chain:

  nyear x niter matrix posterior sample of the proportion of sites
  occupied in each year.

- muz.missing.chain:

  nyear x niter matrix posterior sample of the proportion of sites
  occupied for sites with missing data.

- prop.extinct.chain:

  Extinction rate for all sites.

- prop.colon.chain:

  Colonization rate.

- mupsi1.chain:

  posterior sample of parameter for occupancy in year 1.

- e.chain:

  posterior sample of e

- x.chain:

  posterior sampmle of x

- y.chain:

  posterior sample of y

- b.chain:

  posterior sample of b

- alpha.chain:

  posterior sample of alpha

## References

Risk, B. B., De Valpine, P., Beissinger, S. R. (2011). A robust design
formulation of the incidence function model of metapopulation dynamics
applied to two species of rails. Ecology, 92(2), 462-474.

## Author

Benjamin Risk

## Examples

``` r
if (FALSE) { # \dontrun{

data(simulatedifm)
library("coda")

niter=2000
nsite=100
nyear=10
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

# NOTE: If posteriors are truncated at zero, then estimates may be biased. Rescale
# distances (e.g., divide by 10,000) and/or areas so that parameters are larger.

nmissing = sum(is.na(z.sim.20))

init1=list(z.missing=runif(nmissing),mupsi1=runif(1),alpha=runif(1,1,30),
 b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1),x=runif(1,0,5))

a = Sys.time()
im1 <- ifm.missing.MCMC(niter=niter,init=init1,z.data = z.sim.20,
 site.distance=sim.distance,site.area=sim.area, sd.prop.mupsi1=0.2, sd.prop.alpha=4, sd.prop.b=0.6,
 sd.prop.y=40, sd.prop.e=0.05, sd.prop.x=0.4, nthin=1, print.by=500)
accept.calculate(im1,model='missing')
Sys.time() - a


init2=list(z.missing = runif(nmissing), mupsi1 = runif(1), alpha=runif(1,1,30),
 b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1),x=runif(1,0,5))
im2 <- ifm.missing.MCMC(niter=niter,init=init2, z.data = z.sim.20, site.distance=sim.distance,
site.area=sim.area, sd.prop.mupsi1=0.2, sd.prop.alpha=4, sd.prop.b=0.6, sd.prop.y=40, 
sd.prop.e=0.05,sd.prop.x=0.4, nthin=1, print.by=1000)
accept.calculate(im2,model='missing')
Sys.time() - a

coda.create(im1,"sim_im1",par.list=list("mupsi1.chain","e.chain","x.chain","alpha.chain",
"b.chain","y.chain"),niter=niter,nthin=nthin)
coda.create(im2,"sim_im2",par.list=list("mupsi1.chain","e.chain","x.chain","alpha.chain",
"b.chain","y.chain"),niter=niter,nthin=nthin)
coda.sim.im1=read.coda("sim_im1.txt","sim_im1_Index.txt")
coda.sim.im2=read.coda("sim_im2.txt","sim_im2_Index.txt")
coda.sim.im.list=mcmc.list(coda.sim.im1,coda.sim.im2)
sim.im=combine.chains(im1,im2,nburnin=nburnin,nthin=1)
coda.create(sim.im,"sim_im",par.list=list("mupsi1.chain","e.chain","x.chain","alpha.chain",
"b.chain","y.chain"),niter=(2*niter-2*nburnin),nthin=nthin)
coda.sim.im.long=read.coda("sim_im.txt","sim_im_Index.txt")

summary(coda.sim.im.list)
summary(coda.sim.im.long)

gelman.diag(coda.sim.im.list)

plot(coda.sim.im.list)
plot(coda.sim.im.long)
cumuplot(coda.sim.im.long)

# calculate maximum a posteriori estimates:
m1 <- as.matrix(sim.im)
e <- calcmode(m1[,1][[1]])
x <- calcmode(m1[,1][[2]])
y <- calcmode(m1[,1][[3]])
b <- calcmode(m1[,1][[4]])
alpha <- calcmode(m1[,1][[5]])

} # }
```
