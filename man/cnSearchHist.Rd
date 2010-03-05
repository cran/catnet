\name{cnSearchHist}
\alias{cnSearchHist}

\title{Parent Histogram Matrix}

\description{Estimation of the parent matrix of nodes from data. 
The frequency of node edges is obtained by fitting networks consistent to randomly generated node orders.}

\usage{
cnSearchHist(data, perturbations, 
	maxParentSet, maxComplexity, 
	parentsPool=NULL, fixedParentsPool=NULL, 
	niter=32, echo=FALSE)
}
\arguments{
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary \code{matrix} with the dimensions of \code{data}. A value \code{1} designates the corresponding node in the sample as perturbed}
\item{maxParentSet}{an \code{integer}, the maximal number of parents per node}
\item{maxComplexity}{an \code{integer}, the maximal network complexity for the search}
\item{parentsPool}{a list of parent sets to choose from}
\item{fixedParentsPool}{a list of parent sets to choose from}
\item{niter}{an \code{integer}, the number of single order searches to be performed}	
\item{echo}{a boolean that sets on/off some functional progress and debug information}
}

\details{After performing \code{niter} calls of \code{cnSearchOrder} for randomly generated node orders (uniformly over the space of all possible node orders), 
the function returns the histogram matrix of the parent-child relations, as given by \code{cnMatParents} function, for the resulting optimal networks. Essentially, the result is sum of \code{cnMatParents} applied to all found networks, across the random orders and complexities up to the maximum specified. 
The entries in the output \code{matrix} shows how many times a given parent-child pair was found. 
\code{cnSearchHist} function can be useful for empirical estimation of the relationships in some multivariate categorical data.}

\value{A \code{matrix}}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=8, maxParents=3, numCategories=2)
  psamples  <- cnSamples(object=cnet, numsamples=100)
  mhisto <- cnSearchHist(data=psamples, perturbations=NULL, 
		maxParentSet=2, maxComplexity=20)
  mhisto
}

\seealso{\code{\link{cnMatParents}}, \code{\link{cnSearchOrder}}, \code{\link{cnSearchHistCluster}} }

\keyword{methods}
\keyword{graphs}
