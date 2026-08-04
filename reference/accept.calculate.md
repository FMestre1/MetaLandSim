# Calculate acceptance rates in MCMC chains

Calculate acceptance rates of parameters in the IFM.

## Usage

``` r
accept.calculate(x, model = c("naive", "missing", "robust"))
```

## Arguments

- x:

  A named list with the MCMC chains estimated by ifm.naive.MCMC,
  ifm.missing.MCMC, or ifm.robust.MCMC.

- model:

  Either "naive", "missing", or "robust"

## Value

Named list containing MCMC chain acceptance rates. Names are built from
the input list, e.g., for model=“naive":

- acc.b.chain:

  Acceptance rates of parameter b

- acc.e.chain:

  Acceptance rates of parameter e

- acc.y.chain:

  Acceptance rates of parameter y

- acc.alpha.chain:

  Acceptance rates of parameter alpha

- acc.x.chain:

  Acceptance rates of parameter x

## Author

Benjamin Risk

## Examples

``` r

data(simulatedifm)

# Here, we run a chain with random initial values:
init1=list(alpha=runif(1,1,30), b=runif(1,0,5),y=runif(1,0,20),e=runif(1,0,1),x=runif(1,0,5))

inm1 <- ifm.naive.MCMC(niter=1000,init=init1,z.data =
 z.sim,site.distance=sim.distance,site.area=sim.area,
  sd.prop.alpha=4,sd.prop.b=0.6,sd.prop.y=40,sd.prop.e=0.05,sd.prop.x=0.4,nthin=1,print.by=100)
#> [1] 1
#> [1] 101
#> [1] 201
#> [1] 301
#> [1] 401
#> [1] 501
#> [1] 601
#> [1] 701
#> [1] 801
#> [1] 901
accept.calculate(inm1,model='naive')
#> $e.chain
#> [1] 0.493
#> 
#> $x.chain
#> [1] 0.341
#> 
#> $y.chain
#> [1] 0.547
#> 
#> $b.chain
#> [1] 0.211
#> 
#> $alpha.chain
#> [1] 0.216
#> 
#> $deviance.chain
#> [1] 0.912
#> 
```
