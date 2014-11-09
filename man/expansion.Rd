\name{expansion}
\alias{expansion}
\title{
Class 'expansion'
}
\description{
Class representing an expansion object, as produced by \code{\link{range_expansion}}.
}
\section{Slots}{
A list of four data frames with the proportion of occupation at several distances from the closest occupied landscape mosaic. These four data frames correspond to the proportion of occupation to the north, south, east and west.
Each data frame has the following columns:
  \itemize{
  \item DISTANCE - Distance (mapsize x number of landscapes).
  \item OCCUPATION - How many times did the landscape at this distance got occupied by the species (from a total of 'iter' repetitions). 
  \item PROPORTION - Proportion of occupation for the landscape at this distance (OCCUPATION/iter).
}
}

\author{
Frederico Mestre and Fernando Canovas
}
