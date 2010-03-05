\name{cnSearchSA}
\alias{cnSearchSA}

\title{Stochastic Network Search}

\description{This function provides a MLE based network search in the space of node orders by Simulated Annealing. For a given sample from an unknown categorical network, it returns a list of \code{catNetwork} objects, with complexity up to some maximal value, that best fit the data.}

\usage{
cnSearchSA(data, perturbations, 
	maxParentSet, maxComplexity=0, 
	parentsPool=NULL, fixedParentsPool=NULL, 
	selectMode = "BIC", 
	tempStart=1, tempCoolFact=0.9, tempCheckOrders=10, maxIter=100, 
	orderShuffles=1, stopDiff=0, priorSearch = NULL, echo=FALSE)
}

\arguments{
\item{data}{a \code{matrix} in row-nodes format or a \code{data.frame} in column-nodes format}
\item{perturbations}{a binary matrix with the dimensions of \code{data}. A value \code{1} designates the node in the corresponding sample as perturbed}
\item{maxParentSet}{an integer, maximal number of parents per node}
\item{maxComplexity}{an integer, maximal network complexity for the search}
\item{parentsPool}{a list of parent sets to choose from}
\item{fixedParentsPool}{a list of parent sets to choose from}
\item{selectMode}{a character, optimization network selection criterion such as "AIC" and "BIC"}
\item{tempStart}{a numerical value, the initial temperature for the annealing}	
\item{tempCoolFact}{a numerical value, the temperature multiplicative decreasing factor}
\item{tempCheckOrders}{an integer, the number of iteration, orders to be searched, with constant temperature}
\item{maxIter}{an integer, the total number of iterations, thus orders, to be processed}
\item{orderShuffles}{an integer, the number of order shuffles per iteration with 0 indicating random order at each iteration}
\item{stopDiff}{a numerical value, stopping epsilon criterion}
\item{priorSearch}{a \code{catNetworkEvaluate} object from a previous search}
\item{echo}{a \code{logical} that sets on/off some functional progress and debug information}
}

\details{
The \code{data} can be a matrix of \code{character} categories with rows specifying the node-variables and columns assumed to be independent samples from an unknown network, or 
a \code{data.frame} with columns specifying the nodes and rows being the samples.

The number of categories for each node is obtained from the data. It is the user responsibility to make sure the data can be categorized reasonably. If the data is numerical it will be forcibly coerced to integer one, which however may result to NA entries or too many node categories, and in both cases the function will fail. 

The function returns a list of networks, one for any possible complexity within the specified range. Stochastic optimization, based on the criterion of maximizing the likelihood, is carried on the network with complexity closest to, but not above, \code{maxComplexity}. 
If \code{maxComplexity} is not specified, thus the function is called with the default zero value, then \code{maxComplexity} is set to be the complexity of a network with all nodes having the maximum, \code{maxParentSet}, the number of parents.
The \code{selectMode} parameter sets the selection criterion for the network upon which the 
maximum likelihood optimization is carried on. "BIC" is the default choice, while any value 
different from "AIC" and "BIC" results in the maximum complexity criterion to be used, 
the one which selects the network with complexity given by \code{maxComplexity}. 

The parameters \code{tempStart}, \code{tempCoolFact} and \code{tempCheckOrders} control the Simulated Annealing schedule.

\code{tempStart} is the starting temperature of the annealing process.

\code{tempCoolFact} is the cooling factor from one temperature step to another. 
It is a number between 0 and 1, inclusively; 
For example, if \code{tempStart} is the temperature in the first step, 
\code{tempStart*tempCoolFact} will be temperature in the second.

\code{tempCheckOrders} is the number of proposals to be checked, 
or with other words, order selections from the current order's neighborhood, at each step before decreasing the temperature. 

\code{maxIter} is the maximum number of orders to be checked. 
If for example \code{maxIter} is 40 and \code{tempCheckOrders} is 4, 
then 10 temperature decreasing steps will be eventually performed.

\code{orderShuffles} is a number that controls the extend of the order neighborhoods.
Each new proposed order is obtained from the last accepted one by \code{orderShuffles} switches of two node indices. 
\code{stopDiff} is a stopping criterion. If at a current temperature, 
after \code{tempCheckOrders} orders being checked, no likelihood improvement of level at least \code{stopDiff} is found, then the SA stops and the function exists. 
Setting this parameter to zero guarantees exhausting all of the maximum allowed 
\code{maxIter} order searches.

\code{priorSearch} is a result from previous search. This parameters allows a new search to be initiated  
from the best order found so far. Thus a chain of searches can be constructed with varying parameters 
providing greater adaptability and user control.

See the vignettes for more details on the algorithm.}

\value{A \code{catNetworkEvaluate} object.}

\author{N. Balov, P. Salzman}

\examples{
  cnet <- cnRandomCatnet(numnodes=12, maxParents=3, numCategories=2)
  psamples  <- cnSamples(object=cnet, numsamples=100)
  nets <- cnSearchSA(data=psamples, perturbations=NULL, 
		maxParentSet=2, maxComplexity=36)
  cc <- cnComplexity(object=cnet)
  cnFind(object=nets, complexity=cc)
}

\seealso{\code{\link{cnSearchOrder}}, \code{\link{cnSearchSAcluster}} }

\keyword{methods}
\keyword{graphs}
