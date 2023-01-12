\name{plot_expansion}
\alias{plot_expansion}
\title{
Graphical display of the expansion simulations
}
\description{
Plots the expansion simulations resulting from \code{\link{range_expansion}}.
}
\usage{
plot_expansion(expansion_object)
}
\arguments{
  \item{expansion_object}{
Output of the function range_expansion.
}
}
\value{
Graphical display of the range_expansion simulation results.
}
\author{
Frederico Mestre and Fernando Canovas
}
\examples{

data(rg_exp)

plot_expansion(exp = rg_exp)

}
