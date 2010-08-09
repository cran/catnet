########################################################################
# Categorical Network Class Methods
# Histograms

##setMethod("cnParHist", "list", 
##          function(objectlist) {
##            if(!is(objectlist, "list"))
##              stop("A list of catNetworks should be specified.")
##            if(length(objectlist)==0 || !is(objectlist[[1]], "catNetwork"))
##              stop("At least one catNetworks should be specified.")
##            return(parHisto(objectlist))
##          })

parHisto <- function(objectlist, norder = NULL) {
  if(is(objectlist, "catNetwork")) {
    n <- objectlist@numnodes
  if(is.null(norder))
    norder <- seq(1, n)
    return(matParents(objectlist, norder))
  }
  
  n <- objectlist[[1]]@numnodes
  if(is.null(norder))
    norder <- seq(1, n)
  
  i <- 1
  nnets <- length(objectlist)
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    if(object@numnodes != n)
      stop("Networks should have the same number of nodes.")
    mpar <- matParents(object, norder)
    if(i==1)
      mhisto <- mpar
    else
      mhisto <- mhisto + mpar

    i <- i + 1
  }
 
  return(mhisto)
}
 
cnSearchHist <- function(data, perturbations,  
                         maxParentSet, parentSizes = NULL,
                         maxComplexity=0,
                         parentsPool = NULL, fixedParentsPool = NULL,
                         selectMode = "BIC", 
                         maxIter = 32, numThreads = 2, echo=FALSE) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

  fast <- TRUE
  
  if(is.matrix(data)) {
    numnodes <- nrow(data)
    numsamples <- ncol(data)
    nodenames <- rownames(data)
  }
  else {
    numnodes <- ncol(data)
    numsamples <- nrow(data)
    nodenames <- colnames(data)
  }
  
  if(numnodes < 1 || numsamples < 1)
    stop("No valid sample is specified.")

  if(length(nodenames) < numnodes) {
    nodenames <- seq(1, numnodes)
  }

  maxParentSet <- as.integer(maxParentSet)
  if(maxParentSet < 1) {
    if(!is.null(parentSizes))
      maxParentSet <- as.integer(max(parentSizes))
    if(maxParentSet < 1) 
      maxParentSet <- 1
  }

  if(!is.null(parentSizes)) {
    parentSizes <- as.integer(parentSizes)
    parentSizes[parentSizes<0] <- 0
    parentSizes[parentSizes>maxParentSet] <- maxParentSet
  }
  
  r <- .categorizeSample(data, perturbations)
  data <- r$data
  perturbations <- r$perturbations
  categories <- r$categories
  maxCategories <- r$maxCategories

  if(maxComplexity <= 0)
    maxComplexity <- as.integer(numnodes * exp(log(maxCategories)*maxParentSet) * (maxCategories-1))
  minComplexity <- sum(sapply(categories, function(cat) (length(cat)-1)))
  if(maxComplexity < minComplexity)
    maxComplexity <- minComplexity
  
  if(is.null(perturbations)) {
    dims <- dim(data)
    perturbations <- matrix(rep(0,dims[1]*dims[2]), dims[1], dims[2])
  }

  numThreads <- as.integer(numThreads)
  if(numThreads < 1)
    numThreads <- 1

  maxIter <- as.integer(maxIter)
  if(maxIter < numThreads)
    maxIter <- numThreads
  
  if(fast) {
    ## call the C-function
    .Call("ccnReleaseCache", PACKAGE="catnet")
    vhisto <- .Call("ccnParHistogram", 
                    data, perturbations, 
                    as.integer(maxParentSet), as.integer(parentSizes),
                    as.integer(maxComplexity),
                    parentsPool, fixedParentsPool,
                    selectMode, as.integer(maxIter),
                    as.integer(numThreads), 
                    ## cache
                    TRUE, 
                    echo, 
                    PACKAGE="catnet")

    mhisto <- matrix(vhisto, numnodes, numnodes)
    rownames(mhisto)<-nodenames
    colnames(mhisto)<-nodenames
    return(mhisto)
  }

  ## R-implementation

  lfact <- 1/(numnodes*numsamples)
  lsum <- 0
  
  norders <- maxIter
  for(k in 1:norders) {

    order <- sample(1:numnodes)
    if(echo)
      cat("\n Order: ", order, "\n")

    t1 <- proc.time()
    
    res <- optimalNetsForOrder(data, perturbations, 
                               categories, as.integer(maxCategories),
                               as.integer(maxParentSet), NULL, 
                               as.integer(maxComplexity), order, 
                               parentsPool, fixedParentsPool, 
                               fast=TRUE, echo=FALSE, useCache=FALSE)
    
    t2 <- proc.time()
    mins <- floor((t2[1]-t1[1])/60)
    secs <- t2[1]-t1[1] - 60*mins
    if(echo)
      cat(k, "\\", norders, "search time: ", mins, "min, ", secs, "sec\n")

    cnet <- cnFind(res, maxComplexity)
    hh <- parHisto(cnet, order)
    laux <- exp(cnet@likelihood * lfact)
    lsum <- lsum + laux
    hh <- hh * laux
    if(k == 1)
      mhisto <- hh
    else
      mhisto <- mhisto + hh
  }

  ##if(lsum > 0)
  ##  mhisto <- mhisto / lsum
  
  rownames(mhisto)<-nodenames
  colnames(mhisto)<-nodenames
  
  return(mhisto)
}
