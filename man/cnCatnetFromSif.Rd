\name{cnCatnetFromSif}
\alias{cnCatnetFromSif}
\alias{cnCatnetFromSif,character-method}

\title{catNetwork from Simple Interaction File (SIF)}

\description{Creates a \code{catNetwork} object from a SIF file.}

\usage{
cnCatnetFromSif(file, numCategories=2)
}

\arguments{
\item{file}{a file name}
\item{numCategories}{an \code{integer}, the number of node categories}
}

\details{
The function imports a graph structure from a SIF file by assigning equal number \code{numCategories} 
of categories for each of its nodes and a random probability model. 
Subsequently, the probability model can be changed by calling \code{cnSetProb} function.
}

\value{A \code{catNetwork} object}

\author{N. Balov, P. Salzman}

\seealso{\code{\link{cnNew}}, \code{\link{cnCatnetFromGraph}}, \code{\link{cnSetProb}} }

\keyword{graphs}
