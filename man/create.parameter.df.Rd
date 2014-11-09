\name{create.parameter.df}
\alias{create.parameter.df}
\title{
Create parameter data frame
}
\description{
This function creates a parameter data frame, using parameter values computed with the application available in the papers of Moilanen (1999) and ter Braak and Etienne (2003).
}
\usage{
create.parameter.df(alpha, x, y, e)
}
\arguments{
  \item{alpha}{
Alpha parameter
}
  \item{x}{
x parameter
}
  \item{y}{
y parameter
}
  \item{e}{
e parameter
}
}
\details{
It is highly recommended that the user reads both papers, as well as the help files. 
}
\value{
Returns a data frame, with the same format as the one returned by \code{\link{parameter.estimate}} for the methods 'Rsnap_1' and 'Rsnap_x'.
}
\references{
Moilanen, A. (1999). Patch occupancy models of metapopulation dynamics: efficient parameter estimation using implicit statistical inference. Ecology, 80(3): 1031-1043.
 
ter Braak, C. J., & Etienne, R. S. (2003). Improved Bayesian analysis of metapopulation data with an application to a tree frog metapopulation. Ecology, 84(1): 231-241.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{parameter.estimate}}
}
\examples{

param2 <- create.parameter.df(alpha=0.5, x=0.1, y=5, e=0.1)

param2

#      par_output
#alpha        0.5
#x            0.1
#y            5.0
#e            0.1


}
