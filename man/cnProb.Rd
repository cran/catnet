\name{cnProb-method}
\alias{cnProb}
\alias{cnProb,catNetwork,missing,missing-method}
\alias{cnProb,catNetwork,vector,character-method}

\title{Conditional Probability Structure of Network}

\description{Returns the list of node conditional probabilities of a \code{catNetwork} object.}

\usage{
cnProb(object, which, file)
}

\arguments{
\item{object}{a \code{catNetwork} object}
\item{which}{a \code{vector} of indices}
\item{file}{a \code{character}, an optional file for the output to be written to}
}

\value{A list of probabilities of nodes with indices \code{which} and 
all node probabilities if \code{which} is not specified. 
If \code{file} is given, the result is written to a file with that name. 
Node probabilities are reported in the following format. First, node name and its parents are given, 
then a list of probability values corresponding to 
all combination of parent categories (put in brackets) and node categories. 
For example, the conditional probability of a node with two parents, 
such that both the node and its parents have three categories, 
is given by 27 values, one for each of the 3*3*3 combination.
}

\author{N. Balov, P. Salzman}

\examples{
  library(catnet)
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnProb(object=cnet)
}

\keyword{methods}
\keyword{graphs}
\keyword{distribution}
