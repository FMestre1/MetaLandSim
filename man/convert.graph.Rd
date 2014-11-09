\name{convert.graph}
\alias{convert.graph}
\title{
Convert data frame to landscape
}
\description{
Converts a given data frame in a list which can be used in the following functions, an object of class 'metapopulation'.
}
\usage{
convert.graph(dframe, mapsize, dispersal)
}
\arguments{
  \item{dframe}{
data frame with the original data and the following columns, in this order:
\itemize{
 \item ID - patch Id.
 \item X - Coordinate.
 \item Y - Coordinate.
 \item Area - Patch area, in hectares.
 \item Occupation - Species presence status (0/1).
}
}
  \item{mapsize}{
Landscape mosaic side length, in meters.
}
  \item{dispersal}{
Species mean dispersal ability, in meters.
}
}
\value{
Delivers an object of class 'metapopulation'.
}
\author{
Frederico Mestre and Fernando Canovas
}
\seealso{
\code{\link{species.graph}}
}
\examples{

data(mc_df)

#Checking the columns of the data frame:

head(mc_df)

#  ID        x       y  area mc
#1  1 1248.254   0.000 0.079  0
#2  2 1420.857  46.725 0.781  1
#3  3 1278.912  52.629 1.053  1
#4  4 6370.625  62.637 0.788  0
#5  5 1151.337  97.140 0.079  0
#6  6 1295.796 104.839 0.137  1


#In order to import the data frame mc_df:

sp1 <- convert.graph(dframe=mc_df, mapsize=8300, dispersal=800)

#verify class

class(sp1)

# [1] "metapopulation"

}