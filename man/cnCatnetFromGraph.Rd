\name{cnCatnetFromGraph}
\alias{cnCatnetFromGraph}
\alias{cnCatnetFromGraph,graph,missing,missing-method}
\alias{cnCatnetFromGraph,graph,list,missing-method}
\alias{cnCatnetFromGraph,graph,list,list-method}

\title{catNetwork from Graph}

\description{Creates a \code{catNetwork} object with the edge structure inherited from a graph.}

\usage{
cnCatnetFromGraph(graph, numCategories = 2, cats=NULL, probs=NULL)
}

\arguments{
\item{graph}{a \code{graph} object}
\item{numCategories}{an \code{integer}, the number of categories per node}
\item{cats}{a \code{list} of categories}
\item{probs}{a \code{list} of probabilities}
}

\details{The user may specify the number of categories (default value is 2), 
the list of node categories and the list of conditional probabilities. 
If no categories or probabilities are given, random ones are generated.}

\value{A \code{catNetwork} object.}

\author{N. Balov}

\examples{
  library(graph)
  myNodes<-c("a","s","p","q","r","t","u")
  myEdges<-list(
		a=list(edges=NULL), 
		s=list(edges=c("p","q")), 
		p=list(edges=c("q")), 
		q=list(edges=("r")), 
		r=list(edges=c("u")), 
		t=list(edges=c("q")), 
		u=list(edges=NULL))
  gr<-new("graphNEL", nodes=myNodes, edgeL=myEdges, edgemode="directed")
  cnet <- cnCatnetFromGraph(graph=gr)
}

\seealso{\code{\link{cnNew}}, \code{\link{cnCatnetFromSif}} }

\keyword{graphs}
