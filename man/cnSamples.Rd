\name{cnSamples-method}
\alias{cnSamples}
\alias{cnSamples,catNetwork-method}
\alias{cnSamples,catNetworkFit-method}

\title{Samples from Network}

\description{Generates samples from of a \code{catNetwork} object.}

\usage{
cnSamples(object, numsamples = 1, output="frame")
}

\arguments{
\item{object}{a \code{catNetwork}}
\item{numsamples}{an \code{integer}, the number of samples to be generated}
\item{output}{a \code{character}, the output format. Can be a \code{data.frame} or \code{matrix}.}
}

\details{If the output format is "matrix" then the resulting sample matrix is in row-node format - the rows correspond to the object's nodes while 
the individual samples are represented by columns. 
If the output format is "frame", which is by default, 
the result is a data frame with columns representing the nodes 
and levels the set of categories of the respected nodes.}

\value{A \code{matrix} or \code{data.frame} of node categories}

\author{N. Balov}
\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=3)
  ## generate a sample of size 100 from cnet
  psamples <- cnSamples(object=cnet, numsamples=100, output="frame")
}

\seealso{\code{\link{cnSamplesPert}} }

\keyword{methods}
\keyword{graphs}
