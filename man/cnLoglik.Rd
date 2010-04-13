\name{cnLoglik-method}
\alias{cnLoglik}
\alias{cnLoglik,catNetwork-method}

\title{Sample Log-likelihood}

\description{Calculate the log-likelihood of a sample with respect to a given \code{catNetwork} object}

\usage{
cnLoglik(object, data, bysample=FALSE)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{data}{a data \code{matrix} given in the column-sample format, or a \code{data.frame} in the row-sample format}
\item{bysample}{a \code{logical}}
}

\details{If \code{bysample} is set to TRUE, the function output is a vector of length the sample size with values the log-likelihood of each individual sample instance. Otherwise, the total sum of the log-likelihood of the sample instances is reported.
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
