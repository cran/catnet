\name{cnEvaluate-method}
\alias{cnEvaluate}
\alias{cnEvaluate,catNetwork-method}

\title{Network Evaluation from Data}

\description{For a \code{catNetwork} and some data, 
the function returns a list of networks, 'compatible' with the given network, 
that best fit the data and compare them.}

\usage{
cnEvaluate(object, data, perturbations=NULL, maxParentSet=0, maxComplexity=0)
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary \code{matrix} with the dimensions of \code{data}. A value of 1 marks the corresponding node in the \code{data} as perturbed. See \code{\link{cnSamplesPert}} for more details.}
\item{maxParentSet}{an integer, the maximal allowed parents per node for the search}
\item{maxComplexity}{an integer, the maximal network complexity for the search}
}

\value{A \code{catNetworkEvaluate} object.}

\details{
The \code{data} can be a matrix of \code{character} categories with rows specifying the node-variables and columns assumed to be independent samples, or a \code{data.frame} with columns specifying the nodes and rows being the samples.

The \code{data} can be a sample generated from the network \code{object}, or may not, but in both cases, it must have the same nodes as the given network. Therefore, the rows, if \code{data} is a \code{matrix}, or 
the columns, if \code{data} is a \code{data.frame}, should be named after the \code{object}'s nodes.
	
The function searches in the class of networks contingent with the topological order 
of the given network, \code{object}. 
The result is a \code{catNetworkEvaluate} object that contains the set of optimal networks 
within the specified complexity range, up to \code{maxComplexity}, and several distance measures between these networks and the \code{object}. See the description of \code{\link{cnCompare}} for more details about the used comparison criteria.

If \code{maxParentSet} is not specified, and so it is 0 by default, then the maximum number of parents of the \code{object} is used.

If \code{maxComplexity} is not specified, and so it is 0 by default, then the complexity of the \code{object} is used as a maximum value.

\code{cnEvaluate} function can be used for reducing a complex network to a less complex one. 
This can be achieved by generating a large sample from the original network, 
then evaluating it to obtain a list of networks with increasing complexity fitting this sample, 
and finally, selecting a smaller network from this list. 
}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=12, maxParents=3, numCategories=2)
  psamples <- cnSamples(object=cnet, numsamples=20, output="matrix")
  perturbations <- matrix(rep(0, length(psamples)), nrow = nrow(psamples))
  nets <- cnEvaluate(object=cnet, data=psamples, perturbations, 
	maxParentSet =2, maxComplexity=36)
  nets
  cnPlot(nets)
}

\seealso{,\code{\link{catNetworkEvaluate-class}}, \code{\link{cnSearchOrder}} }

\keyword{methods}
\keyword{graphs}
