\name{cnNodeMarginalProb-method}
\alias{cnNodeMarginalProb}
\alias{cnNodeMarginalProb,catNetwork-method}

\title{Node Marginal Probability}

\description{For a given network and node, the function finds the marginal probability of that node.
}

\usage{
cnNodeMarginalProb(object, node)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{node}{an \code{integer}, a node index in \code{object}}
}

\value{a numeric \code{vector}}

\author{N. Balov}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnNodeMarginalProb(cnet, node=5)
}

\seealso{\code{\link{cnProb}} }

\keyword{methods}
\keyword{graphs}
