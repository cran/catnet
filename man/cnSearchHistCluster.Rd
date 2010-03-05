\name{cnSearchHistCluster}
\alias{cnSearchHistCluster}

\title{Parent Histogram Matrix - Parallel Version}

\description{Estimation of the parent matrix of nodes from data by parallel processing.}

\usage{
cnSearchHistCluster(data, perturbations, 
	maxParentSet, maxComplexity, 
	parentsPool=NULL, fixedParentsPool=NULL, niter=32, 
	clusterNodes=2, clusterHost="localhost", echo=FALSE)
}

\arguments{
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary \code{matrix} with the dimensions of \code{data}. A value \code{1} designates the corresponding node in the sample as perturbed}
\item{maxParentSet}{an \code{integer}, the maximal number of parents per node}
\item{maxComplexity}{an \code{integer}, the maximal network complexity for the search}
\item{parentsPool}{a list of parent sets to choose from}
\item{fixedParentsPool}{a list of parent sets to choose from}
\item{niter}{an \code{integer}, the number of single order searches to be performed}	
\item{clusterNodes}{number of processing units to be used}
\item{clusterHost}{IP address of the cluster host; 'localhost' as default}
\item{echo}{a boolean that sets on/off some functional progress and debug information}
}

\details{
After performing \code{niter} calls of \code{cnSearchOrder} for randomly generated node orders (uniformly over the space of all possible node orders), 
the function returns the histogram matrix of the parent-child relations, as given by \code{cnMatParents} function, for the resulting optimal networks.
This is essentially the parallel version of \code{cnSearchHist}.}

\value{A \code{matrix}}

\author{N. Balov}

\examples{
  library(snow)
  cnet <- cnRandomCatnet(numnodes=8, maxParents=3, numCategories=2)
  psamples  <- cnSamples(object=cnet, numsamples=100)
  mhisto <- cnSearchHistCluster(data=psamples, perturbations=NULL, 
		maxParentSet=2, maxComplexity=20)
  mhisto
}

\seealso{\code{\link{cnMatParents}}, \code{\link{cnSearchOrder}}, \code{\link{cnSearchHist}} }

\keyword{methods}
\keyword{graphs}
