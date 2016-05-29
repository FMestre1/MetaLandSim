\name{rg_exp}
\alias{rg_exp}
\docType{data}
\title{
List with range.expansion output
}
\description{
Output of \code{\link{range_expansion}}. Object of class 'expansion'.
}
\usage{data(rg_exp)}
\format{
List of four data frames ('NORTH', 'SOUTH', 'EAST' and 'WEST') with the probability of occupations at several distances from the closest occupied landscape mosaic. These four data frames correspond to the probability of occupation to the north, south, east and west.
Each data frame has the following columns:
  \itemize{
  \item DISTANCE - Distance (mapsize x number of landscapes).
  \item OCCUPATION - How many times did the landscape at this distance got occupied by the species (from a total of 'iter' repetitions). 
  \item PROPORTION - Proportion of occupation for the landscape at this distance (OCCUPATION/iter).
  \item TIME STEP - The average time step during which a given distance is reached.
  }
}
\examples{

data(rg_exp)

}
\keyword{datasets}
