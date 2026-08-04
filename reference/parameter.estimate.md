# Estimate parameters

Estimates the parameters of the Stochastic Patch Occupancy Model with
the following approaches: regression of snapshot data (Hanski, 1994);
Monte Carlo simulation (Moilanen, 1999) and Bayesian MCMC on the full
dataset (ter Braak and Etienne, 2003).

## Usage

``` r
parameter.estimate(sp, method, alpha = NULL, nsnap)
```

## Arguments

- sp:

  Object of class 'metapopulation' with real patch occupancy data of the
  focal species.

- method:

  Method to be used in parameter estimation. Available methods:

  - Rsnap_1 - Regression of snapshot data, using one snapshot (code
    based on Oksanen, 2004).

  - Rsnap_x - Regression of snapshot data, using more than one snapshot
    (code based on Oksanen, 2004).

  - MCsim - Monte Carlo simulation.

  - norescue - Bayesian MCMC, not considering Rescue effect.

  - rescue - Bayesian MCMC, considering Rescue effect.

- alpha:

  Bolean (TRUE/FALSE). Estimate the alpha parameter.

- nsnap:

  Number of snapshots considered.

## Details

Parameter alpha describes the effect of distance to dispersal (inverse
of the average dispersal distance). Parameter x describes de dependence
of the extinction risk on patch size, and consequently on population
dimension. Parameter y scales colonization with connectivity. Parameter
e is the intrinsic extinction rate of local populations, which is the
extinction rate not considering immigration. In the current version the
methods 'MCsim', 'rescue' and 'norescue' only create the files to be
used in the applications already available. Future versions should allow
the direct estimation of parameters without the need for the
applications of Moilanen (1999) and Ter Braak and Etienne (2003).  
Regarding the method 'MCsim' the settings file produced (.set) by
default has the method Nlr (non-linear regression) chosen. The user
should read the file readme.txt, available with the application, where a
three step estimation process is described. The objective is to produce
the priors for the Monte Carlo simulation to run.  
It is highly recommended that the user reads both papers that provide
the applications to compute the methods 'MCsim', 'rescue' and
'norescue'. Several editions to the settings and parameters files of
both applications might be needed in order to customize the estimation
process. This function only generates the input files with the basic
needed structure.  
Parameter estimation is not the main purpose of this package. As such,
the user can estimate the parameters using other available software
tools and then apply the estimated parameters in the simulations. The
function
[`create.parameter.df`](https://fmestre1.github.io/MetaLandSim/reference/create.parameter.df.md)
can be used to create the data frame of the basic spom parameters. Other
required parameters can be directly given as arguments to the
[`iterate.graph`](https://fmestre1.github.io/MetaLandSim/reference/iterate.graph.md),
[`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md) or
[`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md)
functions.  
The application of the Moilanen paper considers the kernel 'op1',
connectivity 'op1', colonization 'op1' and extinction 'op1'. This SPOM
(Stochastic Patch Occupancy Model) is known as Incidence Function Model
(Hanski,1994 and 1999). In the original version of the mode b=1.However
this might be an useful parameter as it scales emigration with patch
area. This parameter can be estimated with field data. Moilanen (1998)
obtained the value for this parameter by regressing the patch area with
known population size.

## Value

With the methods 'Rsnap_1' and 'Rsnap_x' eturns a data frame with 4 rows
displaying the four parameters (alpha, x, y, e) to be passed to
[`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md):

- alpha - Parameter relating extinction with distance.

- y - Parameter y in the colonization probability.

- e - Parameter defining the extinction probability in a patch of unit
  area.

- x - Parameter scaling extinction risk with patch area.

Regarding the methods 'MCsim', 'rescue' and 'norescue' it returns the
files to be used as input in the applications. The files will be saved
in the working directory. After running the applications, a data frame
can be created in R using the function
[`create.parameter.df`](https://fmestre1.github.io/MetaLandSim/reference/create.parameter.df.md).
This will return a data frame with the same structure as the first two
methods.

## References

Hanski, I. (1994). A practical model of metapopulation dynamics. Journal
of Animal Ecology, 63: 151-162.

Hanski, I. (1999). *Metapopulation Ecology*. Oxford University Press.
313 pp.

Hanski, I., Alho, J. and Moilanen, A. (2000) Estimating the parameters
of survival and migration of individuals in metapopulations. Ecology,
81, 239-251.

Moilanen, A. (1998). Long-term dynamics in a metapopulation of the
American Pika. The American Naturalist, 152(4), 530-542.

Moilanen, A. (1999). Patch occupancy models of metapopulation dynamics:
efficient parameter estimation using implicit statistical inference.
Ecology, 80(3): 1031-1043.

Oksanen, J. (2004). Incidence Function Model in R. url.:.
http://cc.oulu.fi/~jarioksa/opetus/openmeta/metafit.pdf.

ter Braak, C. J., & Etienne, R. S. (2003). Improved Bayesian analysis of
metapopulation data with an application to a tree frog metapopulation.
Ecology, 84(1): 231-241.

## Author

Frederico Mestre and Fernando Canovas

## Note

A vignette is available with detailed information about the computation
of the parameters using each method. The method 'MCsim' creates the
files (data and settings files) to be used with the application
available with the paper by Moilanen (1999). The methods 'rescue' and
'norescue' create the files (data, parameters and distance files)to be
used with the application available with the paper by ter Braak and
Etienne (2003).  
The application by Moilanen is available in
<http://www.esapubs.org/archive/ecol/E080/003/>. The application by ter
Braak and Etienne is available in
<http://www.esapubs.org/archive/ecol/E084/005/suppl-1.htm>.

## See also

[`create.parameter.df`](https://fmestre1.github.io/MetaLandSim/reference/create.parameter.df.md),
[`iterate.graph`](https://fmestre1.github.io/MetaLandSim/reference/iterate.graph.md),
[`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md)
and [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md)

## Examples

``` r

data(occ.landscape)

#Using the Regression of snapshot data:

param1 <- parameter.estimate (sp=occ.landscape, method="Rsnap_1")
```
