\name{calcmode}
\alias{calcmode}
\title{
Function for mode estimation of a continuous variable 
}
\description{
Derives the mode, estimating the value of a continuous variable.
}
\usage{calcmode(data,adjust=1)}
\arguments{
  \item{data}{
vector used to estimate the mode.
}
  \item{adjust}{
increase this value to make the density estimate smoother.
}
}
\value{
Returns the numeric value of the mode.
}
\author{
Adapted from \url{https://stat.ethz.ch/pipermail/r-help/2008-August/172323.html}.
}
\examples{
vect1 = rchisq(1000,df=3)
calcmode(vect1)
vect1
}
