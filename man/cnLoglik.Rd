\name{cnLoglik-method}
\alias{cnLoglik}
\alias{cnLoglik,catNetwork-method}

\title{Sample Log-likelihood}

\description{Calculate the log-likelihood of a sample with respect to a given \code{catNetwork} object}

\usage{
cnLoglik(object, data)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{data}{a data \code{matrix} given in the column-sample format, or a \code{data.frame} in the row-sample format}
}

\value{A numeric value}

\author{N. Balov}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents = 3, numCategories = 2)
  psamples <- cnSamples(object=cnet, numsamples=100)
  cnLoglik(object=cnet, data=psamples)
}

\seealso{\code{\link{cnNodeLoglik}} }

\keyword{methods}
\keyword{graphs}
\keyword{distribution}
