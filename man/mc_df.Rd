\name{mc_df}
\alias{mc_df}
\docType{data}
\title{
Modified patch occupancy data of Cabrera vole as a data frame
}
\description{
One season patch occupancy dataset for \emph{Microtus cabrerae} in SW 
Portugal (modified). This dataset is in a format directly used by {\link{convert.graph}}
and converted to an object class 'metapopulation'.
}
\usage{data(mc_df)}
\format{
  A data frame with 685 observations on the following 5 variables.
  \describe{
    \item{\code{ID}}{Patch Id.}
    \item{\code{x}}{X coordinate.}
    \item{\code{y}}{Y coordinate.}
    \item{\code{area}}{Patch area, in hectares.}
    \item{\code{mc}}{Occupancy state (0/1).}
  }
}
\details{
To create this sample dataset the occupancy status of patches was 
scrambled, however the proportion of occupied patches was kept.
}
\source{
Original field data was obtained during project PERSIST (PTDC/BIA-BEC/105110/2008).
}
\examples{

##To be converted in a object of class "metapopulation":
#mc1 <- convert.graph(dframe=mc_df,mapsize=8200,dispersal=800)

data(mc_df)

#Check the columns:

head(mc_df)

#  ID        x       y  area mc
#1  1 1248.254   0.000 0.079  0
#2  2 1420.857  46.725 0.781  1
#3  3 1278.912  52.629 1.053  1
#4  4 6370.625  62.637 0.788  0
#5  5 1151.337  97.140 0.079  0
#6  6 1295.796 104.839 0.137  1


}
\keyword{datasets}
