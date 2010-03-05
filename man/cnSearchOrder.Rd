\name{cnSearchOrder}
\alias{cnSearchOrder}

\title{Network Search for Given Node Order}

\description{The function implements a MLE based algorithm to search for optimal networks complying 
with a given node order. It returns a list of  networks, with complexities up to some maximal value, that best fit the data. }

\usage{
cnSearchOrder(data, perturbations=NULL, 
	maxParentSet=2, maxComplexity=0, nodeOrder=NULL, 
	parentsPool=NULL, fixedParentsPool=NULL, echo=FALSE)
}

\arguments{
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary matrix with the dimensions of \code{data}. A value 1 marks that the node in the corresponding sample as perturbed}
\item{maxParentSet}{an integer, the maximal number of parents per node with default value 2}
\item{maxComplexity}{an integer, the maximal network complexity for the search}
\item{nodeOrder}{a vector of node order; the search is among the networks with this topological order}
\item{parentsPool}{a list of parent sets to choose from}
\item{fixedParentsPool}{a list of parent sets to choose from}
\item{echo}{a \code{logical} that sets on/off some functional progress and debug information}
}

\details{
The \code{data} can be a matrix of \code{character} categories with rows specifying the node-variables and columns assumed to be independent samples from an unknown network, or 
a \code{data.frame} with columns specifying the nodes and rows being the samples.

The number of node categories are obtained from the sample. The function returns a list of networks, one for each admissible complexity within the 
specified range. The networks in the list are the Maximum Likelihood estimates in the class of networks having the given topological order of the nodes and complexity.
When \code{maxComplexity} is not given, thus zero, its value is reset to the maximum possible complexity for the given parent set size. When \code{nodeOrder} is not given or \code{NULL}, the order of the nodes in the data is taken, \code{1,2,...}. 

Parameters \code{parentsPool} and \code{fixedParentsPool} allow the user to put some exclusion/inclusion 
constrains on the possible parenthood of the nodes. They should be given as lists of index vectors, 
one for each node in the data sample.
}

\value{A \code{catNetworkEvaluate} object}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=12, maxParents=3, numCategories=2)
  psamples <- cnSamples(object=cnet, numsamples=100)
  nodeOrder <- sample(1:12)
  nets <- cnSearchOrder(data=psamples, perturbations=NULL, 
		maxParentSet=2, maxComplexity=36, nodeOrder)
  ## next we find the network with complexity of the original one and plot it
  cc <- cnComplexity(object=cnet)
  cnFind(object=nets, complexity=cc)
}

\seealso{\code{\link{cnSearchSA}}, \code{\link{cnEvaluate}} }

\keyword{methods}
\keyword{graphs}
