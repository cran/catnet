########################################################################
# Categorical Network Class Methods
# Searching

cnClass2SearchOrder <- function(
                                data1, perturbations1 = NULL,
                                data2, perturbations2 = NULL,
                                maxParentSet = 2, maxComplexity = 0,
                                nodeOrder = NULL,  
                                parentsPool = NULL, fixedParentsPool = NULL, 
                                echo = FALSE) {

  if(!is.matrix(data1) && !is.data.frame(data1))
    stop("data should be a matrix or data frame")

  if(!is.matrix(data2) && !is.data.frame(data2))
    stop("data should be a matrix or data frame")
  
  fast <- TRUE
  
  t1 <- proc.time()

  perturbations <- NULL
  if(is.matrix(data1)) {
    numnodes1 <- nrow(data1)
    numsamples1 <- ncol(data1)
    data <- data1
    if(!is.null(perturbations1))
       perturbations <- perturbations1
  }
  else {
    numnodes1 <- ncol(data1)
    numsamples1 <- nrow(data1)
    data <- as.matrix(t(data1))
    if(!is.null(perturbations1))
      perturbations <- as.matrix(t(perturbations1))
  }

  if(is.matrix(data2)) {
    numnodes2 <- nrow(data2)
    numsamples2 <- ncol(data2)
    data <- cbind(data, data2)
    if(!is.null(perturbations2))
      perturbations <- cbind(perturbations, perturbations2)
  }
  else {
    numnodes2 <- ncol(data2)
    numsamples2 <- nrow(data2)
    data <- cbind(data, as.matrix(t(data2)))
    if(!is.null(perturbations2))
      perturbations <- cbind(perturbations, as.matrix(t(perturbations2)))
  }
  
  if(maxParentSet < 1)
    maxParentSet <- 1
  if(numnodes1 != numnodes2)
    stop("Data shuld have equal number of variables")
  numnodes <- numnodes1
  if(numnodes < 1 || numsamples1 < 1 || numsamples2 < 1)
    stop("No valid sample is specified.")

  r <- .categorizeSample(data, perturbations)
  data1 <- r$data[,1:numsamples1]
  data2 <- r$data[,(numsamples1+1):(numsamples1+numsamples2)]
  perturbations1 <- r$perturbations[,1:numsamples1]
  perturbations2 <- r$perturbations[,(numsamples1+1):(numsamples1+numsamples2)]
  
  categories <- r$categories
  maxCategories <- r$maxCategories

  if(maxComplexity <= 0)
    maxComplexity <- as.integer(numnodes * exp(log(maxCategories)*maxParentSet) * (maxCategories-1))
  minComplexity <- sum(sapply(categories, function(cat) (length(cat)-1)))
  if(maxComplexity < minComplexity) {
    cat("set maxComplexity to ", minComplexity, "\n")
    maxComplexity <- minComplexity
  }
  
  if(is.null(perturbations1)) {
    dims <- dim(data1)
    perturbations11 <- matrix(rep(0,dims[1]*dims[2]), dims[1], dims[2])
  }

  if(is.null(perturbations2)) {
    dims <- dim(data2)
    perturbations12 <- matrix(rep(0,dims[1]*dims[2]), dims[1], dims[2])
  }
  
  if(is.null(nodeOrder))
    nodeOrder <- 1:numnodes

   ## call the C-function
  
  nodenames <- seq(1,numnodes)
  if(length(dim(data)) == 2)
    nodenames <- dimnames(data)[[1]]
  bestnets <- .Call("ccnOptimalClass2NetsForOrder", 
                    data1, perturbations1,
                    data2, perturbations2, 
                    maxParentSet, maxComplexity, nodeOrder,
                    parentsPool, fixedParentsPool, FALSE, echo, 
                    PACKAGE="catnet")
  if(length(nodenames) == numnodes && length(bestnets) > 0) {
    for(i in 1:length(bestnets)) {
      if(is.null(bestnets[[i]]))
        stop("Failed to find a network")
      bestnets[[i]]@nodes <- nodenames[nodeOrder]
    }
  }
  
  eval <- new("catNetworkEvaluate", numnodes, numsamples1+numsamples2, 1)

  eval@nets <- bestnets
  for(i in 1:length(bestnets)) {
    eval@complexity[i] <- bestnets[[i]]@complexity
    eval@loglik[i] <- bestnets[[i]]@likelihood
    eval@nets[[i]]@categories <- categories[nodeOrder]
  }

  t2 <- proc.time()
  eval@time <- eval@time + as.numeric(t2[3] - t1[3])
  
  return(eval)
}
