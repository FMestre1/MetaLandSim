\name{param2}
\alias{param2}
\docType{data}
\title{
Sample parameter data frame number 2
}
\description{
Sample data frame, as produced by \code{\link{parameter.estimate}}. These parameters are to be passed to \code{\link{spom}}. These are made up parameters, not related to any species.
}
\usage{data(param1)}
\format{
A data frame with 4 rows displaying the four parameters (alpha, x, y, e) to be passed to \code{\link{spom}}:
\itemize{
 \item alpha - Parameter relating extinction with distance.
 \item y - Parameter y in the colonization probability.  
 \item e - Parameter defining the extinction probability in a patch of unit area.
 \item x - Parameter scaling extinction risk with patch area.
 }
}
\details{
The four parameters are to be passed to \code{\link{spom}}.
}
\examples{

data(param1)

param1

#      par_output
#alpha 0.00250000
#x     0.50000000
#y     2.00000000
#e     0.04662827

}
\keyword{datasets}
