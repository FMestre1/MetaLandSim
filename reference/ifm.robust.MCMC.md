# Estimate the robust design incidence function model

Estimates the IFM with imperfect detection and missing data.

## Usage

``` r
ifm.robust.MCMC(niter = 1000, init, det.data, site.distance, site.area, sd.prop.p = 0.1,
sd.prop.mupsi1 = 0.1, sd.prop.e = 0.2, sd.prop.x = 0.2, sd.prop.y = 0.2, sd.prop.b = 0.2,
sd.prop.alpha = 0.2, nthin = 1, nsite.subset = 5, print.by = 100)
```

## Arguments

- niter:

  Number of iterations in the MCMC chain.

- init:

  Named list with values to initialize the chain. E.g.:  
    
  init1=list(z.data=initocc,z.missing=runif(nmissing),p=runif(nyear,  
    
  0.1,1),mupsi1=runif(1),alpha=runif(1,1,30),
  b=runif(1,0,5),y=runif(1,0,20),  
    
  e=runif(1,0,1),x=runif(1,0,5)).  
    
  z.data: a matrix with nrows = number of sites and ncol = number of
  years. Contains NAs for missing values. Contains naive estimates of
  occupancy elsewhere.  
    
  z.missing: z.missing: a vector of initial occupancy states for the
  missing data with length equal to the number of NAs in z.data (i.e.,
  vectorized across years). Can use runif(nmissing).  
    
  p: vector of length nyears with inital probability of detection in
  each year  
    
  mupsi1: probability of initial occupancy in year 1; runif(1)
  suffices  
    
  alpha: initial value for alpha in dispersal model; described as 1 /
  average dispersal distance  
    
  b: initial value for parameter b in colonization model  
    
  y: initial value for parameter y in colonization model  
    
  e: initial value for e in extinction model  
    
  x: initial value for x in extinction model

- det.data:

  Detection data in an array with dimensions nsites x nyears x nvisits.
  For removal design, set all values after a detection equal to NA. For
  missing data in a given year, set all visits to NA.

- site.distance:

  nsite x nsite matrix of distances between sites. The tuning parameters
  in the example are set for distances less than one, with max distance
  approximately 0.5. Input data should have a similar scaling.

- site.area:

  Vector of length nsite with areas. The tuning parameters in the
  example are set for average area approximately equal to 1. Input data
  should have a similar scaling.

- sd.prop.p:

  Scalar equal to the standard deviation of the proposal distribution
  for probability of detection, which is a normal distribution centered
  at current value in the mcmc chain. The same standard deviation is
  used for all years.

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

- p.chain:

  nyear x niter sample of detection probabilities.

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

- latent.deviance.chain:

  posterior sample of -2\*loglik

## References

Risk, B. B., De Valpine, P., Beissinger, S. R. (2011). A robust design
formulation of the incidence function model of metapopulation dynamics
applied to two species of rails. Ecology, 92(2), 462-474.

## Author

Benjamin Risk

## Examples
