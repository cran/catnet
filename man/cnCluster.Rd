\name{ccnCluster-method}
\alias{cnCluster}
\alias{cnCluster,catNetwork-method}

\title{Clusters in Network}

\description{
	The function retrieves the clusters, these are the connected sub-networks, in a given network.
}

\usage{
 cnCluster(object)
}

\arguments{
\item{object}{a \code{catNetwork}}
}

\details{The function constructs a list of subsets of nodes of the \code{object}, each representing a connected sub-network. Isolated nodes, these are nodes not connected to any other, are not reported. Thus, every element of the output list contains at least two nodes.
}

\value{A \code{list} of named nodes.}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=30, maxParents=2, numCategories=2)
  cnCluster(object=cnet)
}

\seealso{\code{\link{cnFindAIC}}, \code{\link{cnFindBIC}} }

\keyword{methods}
\keyword{graphs}
