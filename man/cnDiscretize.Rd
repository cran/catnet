\name{cnDiscretize}
\alias{cnDiscretize}

\title{Data Categorization}

\description{Numerical data discretization using empirical quantiles.
}

\usage{
cnDiscretize(data, numCategories)
}

\arguments{
\item{data}{a numerical \code{matrix} or \code{data.frame}}
\item{numCategories}{an \code{integer}, the number of categories per node}
}

\details{
The numerical \code{data} is discretized into given number of categories, \code{numCategories}, using the empirical node quantiles. 
As in all functions of \code{catnet} package that accept data, if the \code{data} parameter is a \code{matrix} then it is organized in the row-node format. If it is a \code{data.frame}, the column-node format is assumed.

For a given set of numbers in some range interval, the quantiles break that range into sub-intervals, each containing equal count of the initial numbers.
A quantile-based discretization method is used as follows. 
For each node, the sample node distribution is constructed, 
which is then represented by a sum of non-intersecting classes separated by the quantile points of the sample distribution. Each node value is assigned the class index in which it falls into. 

Currently, the function assigns equal number of categories for each node of the data.
}

\value{A \code{matrix} or \code{data.frame} of indices.}

\seealso{\code{\link{cnSamples}}}

\author{N. Balov, P. Salzman}

