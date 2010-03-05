\name{cnSamplesPert-method}
\alias{cnSamplesPert}
\alias{cnSamplesPert,catNetwork-method}
\alias{cnSamplesPert,catNetworkFit-method}

\title{Perturbed Samples from Network}

\description{Generates perturbed samples from of a \code{catNetwork} object.}

\usage{
cnSamplesPert(object, perturbations, numsamples=1, output="frame")
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{perturbations}{a \code{matrix}}
\item{numsamples}{an \code{integer}, the number of samples to be generated}
\item{output}{a \code{character}, the output format. Can be a \code{data.frame} or \code{matrix}.}
}

\details{A perturbed sample is a sample having nodes with predefined, thus not random, values. 
Non-perturbed nodes, the nodes which values have to be set, are designated with zeros in the \code{perturbation} matrix and their values are generated conditional on the values of their parents.
The non-zero values in the perturbation matrix remain the same in the sample matrix
but the zero ones are filled according to the probability law of the network.
As in \code{cnSamples} function, the result can be in matrix or data frame format, 
depending on the \code{output} parameter.}

\value{A \code{matrix} or \code{data.frame} of node categories}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=3)
  nsamples <- 20
  perturbations <- matrix(rbinom(10*nsamples, 2, 0.4), nrow=10)
  ## generate a perturbed sample of size 100 from cnet
  psamples <- cnSamplesPert(object=cnet, perturbations, numsamples=nsamples)
}

\seealso{\code{\link{cnSamples}} }

\keyword{methods}
\keyword{graphs}
