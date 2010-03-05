\name{cnSetProb-method}
\alias{cnSetProb}
\alias{cnSetProb,catNetwork-method}
\alias{cnSetProb,catSampleNetwork-method}

\title{Set Probability from Data}

\description{The function sets the probability structure of a network from data according to the Maximum Likelihood criterion.}

\usage{
cnSetProb(object, data)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{data}{a data \code{matrix} or \code{data.frame}}
}

\details{The \code{data} can be a \code{matrix} of indices with rows representing the network nodes, 
or a \code{data.frame} in the alternative, node-column format. 
The indices range from 1 to the number of categories of the corresponding node.
The function changes the \code{probability} slot of the \code{object} and returns an updated \code{catNetwork}.
It is a procedure of fitting a network to data by keeping the graph structure unchanged. 
}

\value{The updated \code{object}}

\author{N. Balov}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=3)
  psamples <- matrix(1+rbinom(10*100, 2, 0.4), nrow=10)
  rownames(psamples) <- cnet@nodes
  newcnet <- cnSetProb(object=cnet, data=psamples)
}

\keyword{methods}
\keyword{graphs}
