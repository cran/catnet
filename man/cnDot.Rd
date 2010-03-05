\name{cnDot-method}
\alias{cnDot}
\alias{cnDot,catNetwork-method}
\alias{cnDot,catNetwork,character-method}
\alias{cnDot,matrix-method}
\alias{cnDot,matrix,character-method}
\alias{cnDot,list-method}
\alias{cnDot,list,character-method}

\title{Network Description File}

\description{The function generates a dot-file, the native storage format for \code{Graphviz} 
software package, that describes the graph structure of a \code{catNetwork} object.}

\usage{
cnDot(object, file="")
}

\arguments{
\item{object}{a \code{catNetwork}, a list of \code{catNetwork}s or a parent matrix}
\item{file}{a \code{character}, an optional output file name}
}

\details{
The function generates a \code{dot}-text file. 
In order to view the corresponding graph, the user needs a \code{dot}-file converter and \code{pdf}/postscript viewer. Two environment variables have to be set also: \code{R_DOTVIEWER} and \code{R_PDFVIEWER} with 
the full path to the executable routines. 

If \code{Graphviz} is installed and environment variable \code{R_DOTVIEWER} is set with the full 
path to the \code{dot} executable file (the routine that converts a \code{dot}-text file to a postscript or \code{pdf}), a \code{pdf} or postscript file is created depending whether \code{object} is a \code{catNetwork} or a list of \code{catNetwork}.

Next, if a \code{pdf}-viewer is available, the created postscript or \code{pdf} file is shown.

If \code{file} is not specified, then the function just prints out the resulting string 
which otherwise would be written into a file.}

\value{A \code{character} or \code{dot}-file}

\author{N. Balov}

\examples{
  cnet <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
  cnDot(object=cnet, file="cnet")
}

\seealso{\code{\link{catnet-package}}, \code{\link{cnPlot}} }

\keyword{methods}
\keyword{graphs}
\keyword{aplot}
