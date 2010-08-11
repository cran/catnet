########################################################################
# Categorical Network Class Methods
# Searching

cnGenOrder <- function(order, shuffles, bjump = FALSE) {
  if(shuffles > 0) {
    neworder <- order
    nn <- length(order)
    for(sh in 1:shuffles) {
      order <- neworder
      n1 <- floor(1+runif(1,0,1)*nn)
      if(bjump) {
      	n2 <- n1
        while(n2==n1)
          n2 <- floor(1+runif(1,0,1)*nn)
      }
      else {
        n2 <- n1 + 1
        if(n1 == nn)
          n2 <- 1
      }
      if(n1 < n2) {
        if(n1 > 1)
          neworder[1:(n1-1)] <- order[1:(n1-1)]
        neworder[n1:(n2-1)] <- order[(n1+1):n2]
        neworder[n2] <- order[n1]
        if(n2 < nn)
          neworder[(n2+1):nn] <- order[(n2+1):nn]
      }
      else {
        if(n2 > 1)
          neworder[1:(n2-1)] <- order[1:(n2-1)]
        neworder[(n2+1):n1] <- order[n2:(n1-1)]
        neworder[n2] <- order[n1]
        if(n1 < nn)
          neworder[(n1+1):nn] <- order[(n1+1):nn]
      }
    }
  }
  else{ 
    neworder <- sample(1:length(order))
  }
  neworder
}

cnDiscretize <- function(data, numCategories, mode = "uniform", qlevels = NULL) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

  if(is.matrix(data)) {
    samples <- data
    outdataframe <- FALSE
  }
  else {
    samples <- t(data)
    outdataframe <- TRUE
  }
  
  numnodes <- dim(samples)[1]
  numsamples <- dim(samples)[2]
  if(numnodes < 1 || numsamples < 1)
    stop("No valid sample is specified.")

  if(is.vector(numCategories)) {
    len <- length(numCategories)
    numcats <- floor(numCategories[len])
    numcats <- rep(numcats, numnodes)
    if(len <= numnodes) {
      for(j in 1:len)
        numcats[j] <- numCategories[j]
    }
    else {
      for(j in 1:numnodes)
        numcats[j] <- numCategories[j]
    }
    for(j in 1:numnodes)
        if(numcats[j] < 2)
          numcats[j] <- 2
  }
  else {
    numcats <- floor(numCategories)
    if(numcats < 2) {
      warning("set numCategories to 2")
      numcats <- 2
    }
    numcats <- rep(numcats, numnodes)
  }

  if(mode == "quantile") {
    qlevels <- vector("list", numnodes)
    for(j in 1:numnodes)
      qlevels[[j]] <- seq(1:(numcats[j]-1))/numcats[j]
  }
  else { 
    if(is.null(qlevels) || !is.list(qlevels) || length(qlevels) < numnodes) {
      qlevels <- vector("list", numnodes)
      for(j in 1:numnodes)
        qlevels[[j]] <- rep(1, numcats[j])
    }
    for(j in 1:numnodes) {
      sl <- sum(qlevels[[j]])
      if(length(qlevels[[j]]) < numcats[j] || sl <= 0) {
        qlevels[[j]] <- rep(1, numcats[j])
        sl <- sum(qlevels[[j]])
      }
      qlevels[[j]] <- qlevels[[j]]/sl
    }
  }

  data <- matrix(rep(0, length(samples)), nrow=dim(samples)[1])
  for(j in 1:numnodes) {
    col <- samples[j,]
    ##levs <- levels(col)
    col <- as.numeric(col)
    if(!is.numeric(col))
      stop("data should be numeric")
    if(mode == "quantile") {
      qq <- quantile(col, qlevels[[j]])
    }
    else { 
      minc <- min(col)
      maxc <- max(col)
      qq <- sapply(1:numcats[j], function(i) minc+(maxc-minc)*sum(qlevels[[j]][1:i]))
    }
    for(i in 1:length(col)) {
      id <- which(qq > col[i])
      if(length(id)>0)
        data[j,i] <- min(id)
      else
        data[j,i] <- numcats[j]
      data[j,] <- as.integer(data[j,])
    }
  }
  rownames(data) <- rownames(samples)
  colnames(data) <- colnames(samples)

  if(outdataframe) {
    data <- as.data.frame(t(data))
    ## R converts integer back to numeric, hate it
    for(j in 1:numnodes)
      data[,j] <- as.integer(data[,j])
  }
  
  return(data)
}


cnSearchOrder <- function(data, perturbations = NULL,
                          maxParentSet = 0, parentSizes = NULL, 
                          maxComplexity = 0,
                          nodeOrder = NULL,  
                          parentsPool = NULL, fixedParentsPool = NULL,
                          edgeProb = NULL, 
                          echo = FALSE) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

  fast <- TRUE
  
  t1 <- proc.time()
  
  if(is.matrix(data)) {
    numnodes <- nrow(data)
    numsamples <- ncol(data)
  }
  else {
    numnodes <- ncol(data)
    numsamples <- nrow(data)
  }

  if(numnodes < 1 || numsamples < 1)
    stop("Insufficient data")

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

  if(!is.null(edgeProb)) {
    if(!is.matrix(edgeProb) || nrow(edgeProb) != numnodes || ncol(edgeProb) != numnodes)
      stop("edgeProb should be square matrix of length the number of nodes")
    for(i in 1:numnodes) {
      edgeProb[i, edgeProb[i,] < 0] <- 0
      edgeProb[i, edgeProb[i,] > 1] <- 1
    }
  }
    
  r <- .categorizeSample(data, perturbations)
  data <- r$data
  perturbations <- r$perturbations
  categories <- r$categories
  maxCategories <- r$maxCategories

  if(maxComplexity <= 0)
    maxComplexity <- as.integer(numnodes * exp(log(maxCategories)*maxParentSet) * (maxCategories-1))
  minComplexity <- sum(sapply(categories, function(cat) (length(cat)-1)))
  if(maxComplexity < minComplexity) {
    cat("set maxComplexity to ", minComplexity, "\n")
    maxComplexity <- minComplexity
  }
  maxComplexity <- as.integer(maxComplexity)
  
  if(is.null(perturbations)) {
    dims <- dim(data)
    perturbations <- matrix(rep(0, dims[1]*dims[2]), dims[1], dims[2])
  }

  if(is.null(nodeOrder))
    nodeOrder <- 1:numnodes

  nodenames <- rownames(data)    
  if(length(nodenames) < numnodes)
    nodenames <- seq(1,numnodes)
  
  if(!fast) {
    bestnets <- optimalNetsForOrder(data, perturbations, 
                                    categories, as.integer(maxCategories),
                                    as.integer(maxParentSet), as.integer(parentSizes), 
                                    as.integer(maxComplexity),
                                    as.integer(nodeOrder),
                                    parentsPool, fixedParentsPool, 
                                    fast, echo, FALSE)
  }
  else {
    ## call the C-function
    .Call("ccnReleaseCache", PACKAGE="catnet")
    bestnets <- .Call("ccnOptimalNetsForOrder", 
                      data, perturbations, 
                      as.integer(maxParentSet), as.integer(parentSizes),
                      as.integer(maxComplexity),
                      as.integer(nodeOrder),
                      parentsPool, fixedParentsPool, edgeProb, 
                      ## no cache
                      FALSE, 
                      echo, 
                      PACKAGE="catnet")
    if(length(nodenames) == numnodes && length(bestnets) > 0) {
      for(i in 1:length(bestnets)) {
        if(is.null(bestnets[[i]])) {
          warning("No network")
          next
        }
        bestnets[[i]]@nodes <- nodenames[nodeOrder]
      }
    }
  }

  if(echo)
    cat("Collating ", length(bestnets), " networks...\n")
  
  eval <- new("catNetworkEvaluate", numnodes, numsamples, length(bestnets))

  for(i in 1:length(bestnets)) {
    if(is.null(bestnets[[i]]))
      next
    eval@complexity[i] <- bestnets[[i]]@complexity
    eval@loglik[i] <- bestnets[[i]]@likelihood
    ## bestnets are ordered according to `nodidx'
    ## must set their categorical values
    bestnets[[i]]@categories <- categories[nodeOrder]

    ## reorder eval@nets[[nn]]'s nodes to match data's nodes
    enetnodes <- bestnets[[i]]@nodes
    ##cat(enetnodes,"\n")
    if(length(nodenames) == numnodes) {
      ord <- sapply(nodenames, function(c) {
        id <- which(enetnodes==c)
        if(length(id)>0)
          return(id[1])
	stop("nodes do not match")
        })
      if(sum(ord != 1:numnodes) > 0) {
        bestnets[[i]] <- cnReorderNodes(bestnets[[i]], ord)
      }
    }
  }
  eval@nets <- bestnets

  t2 <- proc.time()
  eval@time <- eval@time + as.numeric(t2[3] - t1[3])
  
  return(eval)
}

.searchSA <- function(eval, 
                      numnodes, numsamples, 
                      data, perturbations, 
                      categories, maxCategories,
                      maxParentSet, parentSizes, 
                      maxComplexity,
                      startorder, 
                      parentsPool, fixedParentsPool, edgeProb, 
                      selectMode, 
                      tempStart, tempCoolFact, tempCheckOrders, 
                      maxIter, orderShuffles, stopDiff,
                      numThreads, echo, saved.seed) {

  ##set.seed(saved.seed)

  nodenames <- rownames(data)
  .Call("ccnReleaseCache", PACKAGE="catnet")
  optnets <- .Call("ccnOptimalNetsSA",
                   nodenames, 
                   data, perturbations, 
                   as.integer(maxParentSet), as.integer(parentSizes),
                   as.integer(maxComplexity),
                   parentsPool, fixedParentsPool, edgeProb, 
                   selectMode, as.integer(startorder),
                   tempStart, tempCoolFact, as.integer(tempCheckOrders), 
                   as.integer(maxIter), orderShuffles, stopDiff,
                   ## threads
                   as.integer(numThreads),
                   ## cache
                   TRUE, 
                   echo, 
                   PACKAGE="catnet")
  
  if(echo)
    cat("Collating ", length(optnets), " networks...\n")
  
  for(nn in 1:length(optnets)) {
    if(is.null(optnets[[nn]])) {
      warning("No network")
      next
    }
    enetnodes <- optnets[[nn]]@nodes
    ord <- sapply(nodenames, function(c) which(enetnodes==c))
    if(sum(ord != 1:numnodes) > 0)    
      optnets[[nn]] <- cnReorderNodes(optnets[[nn]], ord)
    ## now the nodes in the network are as that in the data
    optnets[[nn]]@categories <- categories
  }
  
  ## replace existing best networks    
  if(!is.null(eval) && length(eval@nets) > 0) {
    enets <- eval@nets
    for(nn in 1:length(enets)) {
      if(is.null(enets[[nn]]) || !is(enets[[nn]], "catNetwork"))
        next
      bnet <- cnFind(optnets, enets[[nn]]@complexity)
      if(!is.null(bnet) && bnet@complexity == enets[[nn]]@complexity) 
        if(bnet@likelihood > enets[[nn]]@likelihood)
          ##replace with better
          enets[[nn]] <- bnet
    }
    eval@nets <- enets
  }
  else {
    eval@nets <- optnets
  }
  
  return(eval)

  fast <- TRUE
  
  jumpShuffles <- FALSE
  if(orderShuffles < 0) {
    jumpShuffles <- TRUE
    orderShuffles <- -orderShuffles
  }
  minShuffles <- floor(orderShuffles)
  orderShuffles <- orderShuffles - minShuffles

  if(fast) {
    ## call the C-function
    .Call("ccnReleaseCache", PACKAGE="catnet")
  }

  ntrials <- length(eval@complexity)
  
  optprob <- -Inf
  optnets <- eval@nets
  if(selectMode == "AIC") {
    optnet <- cnFindAIC(optnets)
  }
  else {
    if(selectMode == "BIC") {
      optnet <- cnFindBIC(optnets, numsamples)
    }
    else {
      optnet <- cnFind(optnets, maxComplexity)
    }
  }
  if(!is.null(optnet)) {
    optprob <- optnet@likelihood
    ##cat("optprob = ", optprob, "\n")
  }
  else {
    ntrials <- 0
  }

  temp <- tempStart
  
  naccept <- 0
  nsteps <- 0

  optorder <- startorder
  if(echo)
    cat("\nStart Order: ", optorder, "\n")
  
  while(nsteps < maxIter) {

    curTempDelta <- 0
    for(i in 1:tempCheckOrders) {

      neworder <- cnGenOrder(optorder, minShuffles+rbinom(1, 1, orderShuffles), jumpShuffles)
      
      newnets <- optimalNetsForOrder(data, perturbations, 
                                     categories, as.integer(maxCategories),
                                     as.integer(maxParentSet), as.integer(parentSizes),
                                     as.integer(maxComplexity),
                                     as.integer(neworder), 
                                     parentsPool, fixedParentsPool, 
                                     fast=TRUE, echo=FALSE, useCache=TRUE)

      if(length(newnets) < 1)
        next
      if(selectMode == "AIC") {
        newnet <- cnFindAIC(newnets)
      }
      else {
        if(selectMode == "BIC")
          newnet <- cnFindBIC(newnets, numsamples)
        else {
          newnet <- cnFind(newnets, maxComplexity)
        }
      }
      if(is.null(newnet))
        stop("no network found")
        
      newprob <- newnet@likelihood

      if(is.null(optnet) && naccept < 1) {
        ## it's the first network set for the search
        naccept <- naccept + 1
        optnets <- newnets
        optnet <- newnet
        optprob <- optnet@likelihood
        
        eval@nets <- optnets
        eval@complexity[ntrials+naccept] <- optnet@complexity
        eval@loglik[ntrials+naccept] <- optnet@likelihood
        nsteps <- nsteps + 1
        next
      }

      if(optprob >= newprob) {
        acceptprob <- exp((newprob - optprob)/temp)
      }
      
      if(optprob < newprob || runif(1,0,1) < acceptprob) {
        ## cat("optprob = ", optprob, ", newprob = ", newprob, "\n")
        ## accept
        optorder <- neworder
        optnets <- newnets
        optnet <- newnet
        optprob <- newprob
        naccept <- naccept + 1

        eval@complexity[ntrials+naccept] <- optnet@complexity
        eval@loglik[ntrials+naccept] <- optprob
        
        if(naccept > 1 && curTempDelta < abs(eval@loglik[ntrials+naccept] - eval@loglik[ntrials+naccept-1]))
          curTempDelta <- abs(eval@loglik[ntrials+naccept] - eval@loglik[ntrials+naccept-1])
      }
      
      nsteps <- nsteps + 1
      
    } ## tempCheckOrders

    if(length(newnets) < 1)
        stop("no networks found; check function parameters")
    
    if(naccept > 1 && curTempDelta < stopDiff)
      break
    
    if(echo)
      cat("\nTemp = ", temp, "Iterations = ", nsteps, "Order: ", optorder, "\n")
    
    temp <- temp * tempCoolFact
  }

  if(echo)
    cat("Collating ", length(eval@nets), " networks...\n")
  
  nodenames <- rownames(data)
  for(nn in 1:length(eval@nets)) {
    if(is.null(eval@nets[[nn]]) || !is(eval@nets[[nn]], "catNetwork"))
      next
    bnet <- cnFind(optnets, eval@nets[[nn]]@complexity)
    if(!is.null(bnet) && bnet@complexity == eval@nets[[nn]]@complexity)
      if(bnet@likelihood > eval@nets[[nn]]@likelihood) {
        ##replace with better
        eval@nets[[nn]] <- bnet
      }
 
    ## reorder eval@nets[[nn]]'s nodes to match data's nodes
    enetnodes <- eval@nets[[nn]]@nodes
    ord <- sapply(nodenames, function(c) which(enetnodes==c))
    if(sum(ord != 1:numnodes) > 0)    
      eval@nets[[nn]] <- cnReorderNodes(eval@nets[[nn]], ord)
    eval@nets[[nn]]@categories <- categories
  }
  
  if(echo)
      cat("Accepted/Total = ", naccept/nsteps, "\n")

  return(eval)
}

cnSearchSA <- function(data, perturbations,
                       maxParentSet=0, parentSizes = NULL, 
                       maxComplexity = 0,
                       parentsPool = NULL, fixedParentsPool = NULL, edgeProb = NULL, 
                       selectMode = "BIC", 
                       tempStart = 1, tempCoolFact = 0.9, tempCheckOrders = 10, 
                       maxIter = 100, orderShuffles = 1, stopDiff = 0,
                       numThreads = 2, 
                       priorSearch = NULL,  ## catNetworkEvaluate
                       echo=FALSE) {

  t1 <- proc.time()
  
  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

  if(is.matrix(data)) {
    numnodes <- nrow(data)
    numsamples <- ncol(data)
  }
  else {
    numnodes <- ncol(data)
    numsamples <- nrow(data)
  }

  if(numnodes < 1 || numsamples < 1)
    stop("Insufficient data")

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

  if(!is.null(edgeProb)) {
    if(!is.matrix(edgeProb) || nrow(edgeProb) != numnodes || ncol(edgeProb) != numnodes)
      stop("edgeProb should be square matrix of length the number of nodes")
    for(i in 1:numnodes) {
      edgeProb[i, edgeProb[i,] < 0] <- 0
      edgeProb[i, edgeProb[i,] > 1] <- 1
    }
  }
  
  if(!is.null(priorSearch)) {
    if(!is(priorSearch, "catNetworkEvaluate")) 
      stop("'priorSearch' should be a valid catNetworkEvaluate object or NULL")
    if(priorSearch@numnodes != numnodes || priorSearch@numsamples != numsamples)
      stop("priorSearch's number of nodes and sample size are not compatible with those of the data")
  }
  
  if(tempStart <= 0) {
    warning("tempStart is set to 1")
    tempStart <- 1
  }
  if(tempCoolFact <= 0 || tempCoolFact > 1) {
    warning("tempCoolFact is set to 1")
    tempCoolFact <- 1
  }
  if(tempCheckOrders < 1) {
    warning("tempCheckOrders is set to 1")
    tempCheckOrders <- 1
  }
  if(maxIter < tempCheckOrders) {
    warning("maxIter is set to tempCheckOrders")
    maxIter <- tempCheckOrders
  }

  r <- .categorizeSample(data, perturbations)
  data <- r$data
  perturbations <- r$perturbations
  categories <- r$categories
  maxCategories <- r$maxCategories

  if(maxComplexity <= 0)
    maxComplexity <- as.integer(numnodes * exp(log(maxCategories)*maxParentSet) * (maxCategories-1))
  minComplexity <- sum(sapply(categories, function(cat) (length(cat)-1)))
  if(maxComplexity < minComplexity) {
    warning("set maxComplexity to ", minComplexity)
    maxComplexity <- minComplexity
  }

  nodenames <- rownames(data)
  if(length(nodenames) < numnodes)
    nodenames <- seq(1, numnodes)

  if(is.null(perturbations)) {
    dims <- dim(data)
    perturbations <- matrix(rep(0, dims[1]*dims[2]), dims[1], dims[2])
  }
  
  if(is.null(priorSearch) || length(priorSearch@nets) < 1) {
    optorder <- sample(1:numnodes)
    eval <- new("catNetworkEvaluate", numnodes, numsamples, 0)
    eval@time <- 0
  }
  else {
    eval <- priorSearch
    optnets <- priorSearch@nets
    if(length(optnets) < 1)
      stop("'priorSearch' has no networks")

    if(selectMode == "AIC") {
      optnet <- cnFindAIC(optnets)
    }
    else {
      if(selectMode == "BIC") {
        optnet <- cnFindBIC(optnets, numsamples)
      }
      else {
        optnet <- cnFind(optnets, maxComplexity)
      }
    }
    optorder <- cnOrder(optnet)

    ## make sure the optnet comes is data compatible
    if(prod(nodenames == optnet@nodes) == 0)
      stop("The data names should correspond to the prior search nodes.")
  }

  saved.seed <- .Random.seed
  eval <- .searchSA(eval,
                    as.integer(numnodes), as.integer(numsamples),
                    data, perturbations, 
                    categories, as.integer(maxCategories),
                    as.integer(maxParentSet), parentSizes, 
                    as.integer(maxComplexity),
                    as.integer(optorder),  
                    parentsPool, fixedParentsPool, edgeProb, 
                    selectMode, 
                    tempStart, tempCoolFact, as.integer(tempCheckOrders), 
                    as.integer(maxIter), orderShuffles, stopDiff,
                    as.integer(numThreads), echo, saved.seed)
  
  t2 <- proc.time()
  eval@time <- eval@time + as.numeric(t2[3] - t1[3])

  return(eval)
}
