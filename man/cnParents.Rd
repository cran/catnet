\name{cnParents-method}
\alias{cnParents}
\alias{cnParents,catNetwork,missing-method}
\alias{cnParents,catNetwork,vector-method}

\title{Network Parents}

\description{Returns the list of parents of a \code{catNetwork} object.}

\usage{
cnParents(object, which)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{which}{a vector of node indices}
}

\value{A list of parents for nodes with indices \code{which}, or all node parents if \code{which} is not specified.}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnParents(object=cnet)
}

\seealso{\code{\link{cnMatParents}}, \code{\link{cnEdges}} }

\keyword{methods}
\keyword{graphs}
