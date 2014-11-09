\name{parameter.estimate}
\alias{parameter.estimate}
\title{
Estimate parameters
}
\description{
Estimates the parameters of the Stochastic Patch Occupancy Model with the following approaches: regression of snapshot data (Hanski, 1994); Monte Carlo simulation (Moilanen, 1999) and Bayesian MCMC on the full dataset (ter Braak and Etienne, 2003). 
}
\usage{
parameter.estimate(sp, method, alpha = NULL, nsnap)
}
\arguments{
  \item{sp}{
Object of class 'metapopulation' with real patch occupancy data of the focal species.
}
  \item{method}{
Method to be used in parameter estimation. Available methods:
\itemize{
 \item Rsnap_1 - Regression of snapshot data, using one snapshot (code based on Oksanen, 2004).
 \item Rsnap_x - Regression of snapshot data, using more than one snapshot (code based on Oksanen, 2004).
 \item MCsim - Monte Carlo simulation.
 \item norescue - Bayesian MCMC, not considering Rescue effect.
 \item rescue - Bayesian MCMC, considering Rescue effect.
 }
}
  \item{alpha}{
Bolean (TRUE/FALSE). Estimate the alpha parameter.
}
  \item{nsnap}{
Number of snapshots considered.
}
}
\details{
Parameter alpha describes the effect of distance to dispersal (inverse of the average dispersal distance). Parameter x describes de dependence of the extinction risk on patch size, and consequently on population dimension. Parameter y scales colonization with connectivity. Parameter e is the intrinsic extinction rate of local populations, which is the extinction rate not considering immigration.
In the current version the methods 'MCsim', 'rescue' and 'norescue' only create the files to be used in the applications already available. Future versions should allow the direct estimation of parameters without the need for the applications of Moilanen (1999) and Ter Braak and Etienne (2003).\cr
Future versions should include the estimation of other parameters, using the virtual migration model (Hanski et al. 2000).\cr
Regarding the method 'MCsim' the settings file produced (.set) by default has the method Nlr (non-linear regression) chosen. The user should read the file readme.txt, available with the application, where a three step estimation process is described. The objective is to produce the priors for the Monte Carlo simulation to run.\cr
It is highly recommended that the user reads both papers that provide the applications to compute the methods 'MCsim', 'rescue' and 'norescue'. Several editions to the settings and parameters files of both applications might be needed in order to customize the estimation process. This function only generates the input files with the basic needed structure.\cr
Parameter estimation is not the main purpose of this package. As such, the user can estimate the parameters using other available software tools and then apply the estimated parameters in the simulations. The function \code{\link{create.parameter.df}} can be used to create the data frame of the basic spom parameters. Other required parameters can be directly given as arguments to the \code{\link{iterate.graph}}, \code{\link{spom}} or \code{\link{range_expansion}} functions.\cr
The application of the Moilanen paper considers the kernel 'op1', connectivity 'op1', colonization 'op1' and extinction 'op1'. This SPOM (Stochastic Patch Occupancy Model) is known as Incidence Function Model (Hanski,1994 and 1999). In the original version of the mode b=1.However this might be an useful parameter as it scales emigration with patch area. This parameter can be estimated with field data. Moilanen (1998) obtained the value for this parameter by regressing the patch area with known population size.
}
\value{
With the methods 'Rsnap_1' and 'Rsnap_x' eturns a data frame with 4 rows displaying the four parameters (alpha, x, y, e) to be passed to \code{\link{spom}}:
\itemize{
 \item alpha - Parameter relating extinction with distance.
 \item y - Parameter y in the colonization probability.  
 \item e - Parameter defining the extinction probability in a patch of unit area.
 \item x - Parameter scaling extinction risk with patch area.
 }
Regarding the methods 'MCsim', 'rescue' and 'norescue' it returns the files to be used as input in the applications. The files will be saved in the working directory. After running the applications, a data frame can be created in R using the function \code{\link{create.parameter.df}}. This will return a data frame with the same structure as the first two methods.
}
\references{

Hanski, I. (1994). A practical model of metapopulation dynamics. Journal of Animal Ecology, 63: 151-162.

Hanski, I. (1999). \emph{Metapopulation Ecology}. Oxford University Press. 313 pp.

Hanski, I., Alho, J. and Moilanen, A. (2000) Estimating the parameters of survival and migration of individuals in metapopulations. Ecology, 81, 239-251.

Moilanen, A. (1998). Long-term dynamics in a metapopulation of the American Pika. The American Naturalist, 152(4), 530-542.

Moilanen, A. (1999). Patch occupancy models of metapopulation dynamics: efficient parameter estimation using implicit statistical inference. Ecology, 80(3): 1031-1043.
 
Oksanen, J. (2004). Incidence Function Model in R. url.:. http://cc.oulu.fi/~jarioksa/opetus/openmeta/metafit.pdf.

ter Braak, C. J., & Etienne, R. S. (2003). Improved Bayesian analysis of metapopulation data with an application to a tree frog metapopulation. Ecology, 84(1): 231-241.

}
\author{
Frederico Mestre and Fernando Canovas
}
\note{
A vignette is available with detailed information about the computation of the parameters using each method. The method 'MCsim' creates the files (data and settings files) to be used with the application available with the paper by Moilanen (1999). The methods 'rescue' and 'norescue' create the files (data, parameters and distance files)to be used with the application available with the paper by ter Braak and Etienne (2003).\cr
The application by Moilanen is available in \url{http://www.esapubs.org/archive/ecol/E080/003/}.
The application by ter Braak and Etienne is available in \url{http://www.esapubs.org/archive/ecol/E084/005/suppl-1.htm}.

}
\seealso{
\code{\link{create.parameter.df}}, \code{\link{iterate.graph}}, \code{\link{range_expansion}} and \code{\link{spom}} 
}
\examples{

data(occ.landscape)

#Using the Regression of snapshot data:

param1 <- parameter.estimate (sp=occ.landscape, method="Rsnap_1")

}
