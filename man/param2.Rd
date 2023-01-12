\name{param2}
\alias{param2}
\docType{data}
\title{
Sample parameter data frame number 2
}
\description{
Sample data frame, as produced by \code{\link{parameter.estimate}}. These parameters are to be passed to \code{\link{spom}}. 
These are based upon those computed for Cabrerae's vole in the paper Mestre et al. (2017) (see references).
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
\references{
Mestre, F., Risk, B., Mira, A., Beja, P., Pita, R. (2017) <doi:10.1016/j.ecolmodel.2017.06.013>
}
\details{
The four parameters are to be passed to \code{\link{spom}}.
}
\examples{

data(param2)

param2

#      par_output
#alpha    0.00047
#x        0.44000
#y       18.15000
#e        0.00480


}
\keyword{datasets}
