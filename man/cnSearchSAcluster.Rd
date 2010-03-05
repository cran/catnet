\name{cnSearchSAcluster}
\alias{cnSearchSAcluster}

\title{Stochastic Network Search. Parallel Processing.}

\description{
This function provides a MLE based network search in the space of node orders by Simulated Annealing and also employing parallel processing.}

\usage{
cnSearchSAcluster(data, perturbations, 
	maxParentSet, maxComplexity=0, 
	parentsPool=NULL, fixedParentsPool=NULL, 
	selectMode = "BIC", 
	tempStart=1, tempCoolFact=0.9, tempCheckOrders=10, maxIter=200, 
	orderShuffles=1, stopDiff=0, 
	priorSearch = NULL, 
	clusterNodes=2, clusterHost="localhost", echo=FALSE)
}

\arguments{
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary matrix with the dimensions of \code{data}. A value \code{1} designates the node in the corresponding sample as perturbed}
\item{maxParentSet}{an integer, maximal number of parents per node}
\item{maxComplexity}{an integer, maximal network complexity for the search}
\item{parentsPool}{a list of parent sets to choose from}
\item{fixedParentsPool}{a list of parent sets to choose from}
\item{selectMode}{a character, optimization network selection criterion such as "AIC" and "BIC"}
\item{tempStart}{a numerical value, the initial temperature for the annealing}	
\item{tempCoolFact}{a numerical value, the temperature multiplicative decreasing factor}
\item{tempCheckOrders}{an integer, the number of iteration with constant temperature}
\item{maxIter}{an integer, the total number of iterations, thus orders, to be searched}
\item{orderShuffles}{an integer, the number of order shuffles per iteration with 0 indicating random order at each iteration}
\item{stopDiff}{a numerical value, stopping epsilon criterion}
\item{priorSearch}{a \code{catNetworkEvaluate} object from a previous search}
\item{clusterNodes}{an integer, the number of processing units to be used}
\item{clusterHost}{IP address of the cluster host; 'localhost' as default}
\item{echo}{a \code{logical} that sets on/off some functional progress and debug information}
}

\details{This is essentially the parallel version of \code{cnSearchSA}.}

\value{A \code{catNetworkEvaluate} object.}

\author{N. Balov, P. Salzman}
\examples{
  cnet <- cnRandomCatnet(numnodes=12, maxParents=3, numCategories=2)
  psamples  <- cnSamples(object=cnet, numsamples=100)
  nets <- cnSearchSAcluster(data=psamples, perturbations=NULL, 
		maxParentSet=2, maxComplexity=30)
  cnFind(object=nets, complexity=30)
}

\seealso{\code{\link{cnSearchOrder}}, \code{\link{cnSearchSA}} }

\keyword{methods}
\keyword{graphs}
