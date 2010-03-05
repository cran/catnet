\name{cnFind-method}
\alias{cnFind}
\alias{cnFind,list-method}
\alias{cnFind,catNetworkEvaluate-method}

\title{Find Network by Complexity}

\description{
	This is a model selection routine that finds a network in a set of networks for a given complexity.
}

\usage{
 cnFind(object, complexity)
}

\arguments{
\item{object}{\code{catNetworkEvaluate} or \code{list} of \code{catNetwork}s}
\item{complexity}{an integer, the desired complexity}
}

\details{The complexity must be at least the number of nodes of the networks. If no network with the requested complexity exists in the list, then the one with the closest complexity is returned. 
}

\value{A \code{catNetwork} object.}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxParents=2, numCategories=2)
  psamples <- cnSamples(object=cnet, numsamples=100)
  netlist <- cnEvaluate(object=cnet, data=psamples, 
		perturbations=matrix(rep(0,10*100), nrow=10))
  bnet <- cnFind(object=netlist, complexity=cnComplexity(cnet))
  bnet
}

\seealso{\code{\link{cnFindAIC}}, \code{\link{cnFindBIC}} }

\keyword{methods}
\keyword{graphs}
