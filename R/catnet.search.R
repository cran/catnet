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

cnDiscretize <- function(data, numCategories) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

  numcats <- floor(numCategories)
  if(numcats < 2) {
    warning("set numCategories to 2")
    numcats <- 2
  }

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
  
  qlevels <- seq(1:(numcats-1))/numcats
  
  data <- matrix(rep(0, length(samples)), nrow=dim(samples)[1])
  for(j in 1:numnodes) {

    col <- samples[j,]
    ##levs <- levels(col)
    col <- as.numeric(col)
    if(!is.numeric(col))
      stop("data should be numeric")
    
    qq <- quantile(col, qlevels)
    for(i in 1:length(col)) {
      id <- which(qq > col[i])
      if(length(id)>0)
        data[j,i] <- min(id)
      else
        data[j,i] <- numcats
    }
  }
  rownames(data) <- rownames(samples)
  colnames(data) <- colnames(samples)

  if(outdataframe) {
    data <- as.data.frame(t(data))
  }
  
  return(data)
}


cnSearchOrder <- function(data, perturbations = NULL,
                          maxParentSet = 2, maxComplexity = 0,
                          nodeOrder = NULL,  
                          parentsPool = NULL, fixedParentsPool = NULL, 
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
  
  if(maxParentSet < 1)
    maxParentSet <- 1
  if(numnodes < 1 || numsamples < 1)
    stop("No valid sample is specified.")

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
  
  if(is.null(perturbations)) {
    dims <- dim(data)
    perturbations <- matrix(rep(0,dims[1]*dims[2]), dims[1], dims[2])
  }

  if(is.null(nodeOrder))
    nodeOrder <- 1:numnodes

  if(!fast) {
    bestnets <- optimalNetsForOrder(data, perturbations, 
                                    categories, as.integer(maxCategories),
                                    as.integer(maxParentSet), as.integer(maxComplexity), nodeOrder,
                                    parentsPool, fixedParentsPool, 
                                    fast, echo, FALSE)
  }
  else {
    ## call the C-function
    .Call("ccnReleaseCache", PACKAGE="catnet")
    nodenames <- seq(1,numnodes)
    if(length(dim(data)) == 2)
      nodenames <- dimnames(data)[[1]]
    bestnets <- .Call("ccnOptimalNetsForOrder", 
                      data, perturbations, 
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
  }
  
  eval <- new("catNetworkEvaluate", numnodes, numsamples, 1)

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


.searchSA <- function(eval, 
                      numnodes, numsamples, 
                      data, perturbations, 
                      categories, maxCategories,
                      maxParentSet, maxComplexity, startorder, 
                      parentsPool, fixedParentsPool,
                      selectMode, 
                      tempStart, tempCoolFact, tempCheckOrders, 
                      maxIter, orderShuffles, stopDiff,
                      fast, echo, saved.seed) {

  set.seed(saved.seed)
  
  t1 <- proc.time()
  
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
      optnet <- optnets[[length(optnets)]]
    }
  }
  if(!is.null(optnet)) {
    optprob <- optnet@likelihood
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
  
  useCache <- TRUE

  while(nsteps < maxIter) {

    curTempDelta <- 0
    for(i in 1:tempCheckOrders) {

      neworder <- cnGenOrder(optorder, minShuffles+rbinom(1, 1, orderShuffles), jumpShuffles)
      
      newnets <- optimalNetsForOrder(data, perturbations, 
                                 categories, maxCategories,
                                 maxParentSet, maxComplexity, neworder, 
                                 parentsPool, fixedParentsPool, 
                                 fast, echo, useCache)

      if(length(newnets) < 1)
        next
      if(selectMode == "AIC") {
        newnet <- cnFindAIC(newnets)
      }
      else {
        if(selectMode == "BIC")
          newnet <- cnFindBIC(newnets, numsamples)
        else 
          newnet <- newnets[[length(newnets)]]
      }
        
      newprob <- newnet@likelihood

      if(naccept < 1) {
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

      if(optprob >= newprob)
        acceptprob <- exp((newprob - optprob)/temp)
      
      if(optprob < newprob || runif(1,0,1) < acceptprob) {
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

  nodenames <- rownames(data)
  for(nn in 1:length(eval@nets)) {
    enet <- eval@nets[[nn]]
    if(is.null(enet) || !is(enet, "catNetwork"))
      next
    bnet <- cnFind(optnets, enet@complexity)
    if(!is.null(bnet))
      if(bnet@likelihood > enet@likelihood) {
        ##eval@nets[[nn]] <- bnet
        enet <- bnet
      }
    
    ## reorder eval@nets[[nn]]'s nodes according to match data's nodes
    enetnodes <- enet@nodes
    ord <- sapply(nodenames, function(c) which(enetnodes==c))
    eval@nets[[nn]] <- cnReorderNodes(enet, ord)
  }
  
  t2 <- proc.time()
  eval@time <- eval@time + as.numeric(t2[3] - t1[3])
  
  if(echo)
      cat("Accepted/Total = ", naccept/nsteps, "\n")

  return(eval)
}

cnSearchSA <- function(data, perturbations, maxParentSet, maxComplexity = 0,
                       parentsPool = NULL, fixedParentsPool = NULL,
                       selectMode = "BIC", 
                       tempStart = 1, tempCoolFact = 0.9, tempCheckOrders = 10, 
                       maxIter = 100, orderShuffles = 1, stopDiff = 0,
                       priorSearch = NULL,  ## catNetworkEvaluate
                       echo=FALSE) {

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
  
  if(maxParentSet < 1)
    maxParentSet <- 1
  if(numnodes < 1 || numsamples < 1)
    stop("No valid sample is specified.")

  if(!is.null(priorSearch)) {
    if(!is(priorSearch, "catNetworkEvaluate")) 
      stop("'priorSearch' should be a valid catNetworkEvaluate object or NULL")
    if(priorSearch@numnodes != numnodes || priorSearch@numsamples != numsamples)
      stop("priorSearch's number of nodes and sample size are not compatible with those of the data")
  }  
    
  t1 <- proc.time()
  
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

  nodenames <- rownames(data)
  
  if(maxComplexity <= 0)
    maxComplexity <- as.integer(numnodes * exp(log(maxCategories)*maxParentSet) * (maxCategories-1))
  minComplexity <- sum(sapply(categories, function(cat) (length(cat)-1)))
  if(maxComplexity < minComplexity) {
    warning("set maxComplexity to ", minComplexity)
    maxComplexity <- minComplexity
  }

  if(is.null(perturbations)) {
    dims <- dim(data)
    perturbations <- matrix(rep(0, dims[1]*dims[2]), dims[1], dims[2])
  }
  
  if(is.null(priorSearch) || length(priorSearch@nets) < 1) {
    optorder <- sample(1:numnodes)
    eval <- new("catNetworkEvaluate", numnodes, numsamples, 1)
    eval@time <- 0
  }
  else {
    eval <- priorSearch
    optnets <- priorSearch@nets
    if(length(optnets) < 1)
      stop("No networks have been found; 'priorSearch' is not valid")

    if(selectMode == "AIC") {
      optnet <- cnFindAIC(optnets)
    }
    else {
      if(selectMode == "BIC") {
        optnet <- cnFindBIC(optnets, numsamples)
      }
      else {
        optnet <- optnets[[length(optnets)]]
      }
    }
    optorder <- cnOrder(optnet)
  }

  saved.seed <- .Random.seed

  eval <- .searchSA(eval,
                    numnodes, numsamples,
                    data, perturbations, 
                    categories, maxCategories,
                    maxParentSet, maxComplexity, optorder,  
                    parentsPool, fixedParentsPool,
                    selectMode, 
                    tempStart, tempCoolFact, tempCheckOrders, 
                    maxIter, orderShuffles, stopDiff,
                    fast=TRUE, echo,
                    saved.seed)

  for(i in 1:length(eval@nets)) {
    if(is.null(eval@nets[[i]]))
      next
    eval@nets[[i]]@categories <- categories
  }
  
  return(eval)
}

cnSearchSAcluster <- function(data, perturbations,  
                       maxParentSet, maxComplexity = 0,
                       parentsPool = NULL, fixedParentsPool = NULL,
                       selectMode = "BIC", 
                       tempStart = 1, tempCoolFact = 0.9, tempCheckOrders = 10, 
                       maxIter = 200, orderShuffles = 1, stopDiff = 0,
                       priorSearch = NULL,  ## catNetworkEvaluate
                       clusterNodes = 2,
                       clusterHost = "localhost",
                       echo = FALSE) {

  if(!require("snow")) {
    stop("Snow not found.")
  }

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
  
  if(maxParentSet < 1)
    maxParentSet <- 1
  if(numnodes < 1 || numsamples < 1)
    stop("No valid sample is specified.")

  if(!is.null(priorSearch)) {
    if(!is(priorSearch, "catNetworkEvaluate")) 
      stop("'priorSearch' should be a valid catNetworkEvaluate object or NULL")
    if(priorSearch@numnodes != numnodes || priorSearch@numsamples != numsamples)
      stop("priorSearch's number of nodes and sample size are not compatible with those of the data")
  }
  
  t1 <- proc.time()
  
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

  if(clusterNodes < 2)
    clusterNodes <- 2
  tempCheckOrders <- floor(0.5 + tempCheckOrders/clusterNodes)

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

  perturbations <- perturbations
  if(is.null(perturbations)) {
    dims <- dim(data)
    perturbations <- matrix(rep(0,dims[1]*dims[2]), dims[1], dims[2])
  }

  if(is.null(priorSearch) || length(priorSearch@nets) < 1) {
    eval <- new("catNetworkEvaluate", numnodes, numsamples, 1)
    ## start from a random order
    startorder <- lapply(1:clusterNodes, function(i) sample(1:numnodes))
  }
  else {
    eval <- priorSearch
    if(selectMode == "AIC") {
      optnet <- cnFindAIC(eval@nets)
    }
    else {
      if(selectMode == "BIC")
        optnet <- cnFindBIC(eval@nets, numsamples)
      else 
        optnet <- eval@nets[[length(eval@nets)]]
    }
    ## start from a prior order
    cnorder <- cnOrder(optnet)
    startorder <- lapply(1:clusterNodes, function(i) cnorder)
  }

  cl <- makeSOCKcluster(rep(clusterHost, clusterNodes))

  ## create new instances
  maxParentSetTemp <- as.integer(maxParentSet)
  maxComplexityTemp <- as.integer(maxComplexity)
  parentsPoolTemp <- parentsPool
  fixedParentsPoolTemp <- fixedParentsPool
  bFastTemp <- TRUE
  ## there is no way to catch the messages from the other processes
  bEchoTemp <- FALSE
  selectModeTemp <- selectMode
  tempStartTemp <- tempStart
  tempCoolFactTemp <- tempCoolFact
  tempCheckOrdersTemp <- tempCheckOrders
  maxIterTemp <- maxIter
  orderShufflesTemp <- orderShuffles
  stopDiffTemp <- stopDiff

  ## clusterApply insists for a local copy of the search function
  searchFunc <-.searchSA

  ## for fixed .Random.seed let all child R processes have fixed, but different, seed status
  saved.seed <- .Random.seed
  len <- length(saved.seed)
  seeds <- sapply(1:clusterNodes, function(i) saved.seed[as.integer(1+(len-1)*runif(1,0,1))])
 
  res <- clusterApply(cl, 1:clusterNodes, function(i)
                      return(searchFunc(eval,
                                        numnodes, numsamples,
                                        data, perturbations, 
                                        categories, maxCategories,
                                        maxParentSetTemp, maxComplexityTemp, startorder[[i]], 
                                        parentsPool=parentsPoolTemp,
					fixedParentsPool=fixedParentsPoolTemp,
                                        selectModeTemp, 
                                        tempStartTemp, tempCoolFactTemp, tempCheckOrdersTemp, 
                                        maxIterTemp, orderShufflesTemp, stopDiffTemp,
                                        bFastTemp, bEchoTemp,
                                        seeds[[i]]))
                      )
  
  maxeval <- eval
  maxprob <- -Inf
  for(cc in 1:length(res)) {
    nnets <- res[[cc]]@nets
    if(length(nnets) < 1)
      next
        
    if(selectMode == "AIC") {
      selnet <- cnFindAIC(nnets)
    }
    else {
      if(selectMode == "BIC")
        selnet <- cnFindBIC(nnets, numsamples)
      else 
        selnet <- nnets[[length(nnets)]]
    }
      
    if(is(selnet, "catNetwork") &&
       maxprob < selnet@likelihood) {
      maxeval <- res[[cc]]
      maxprob <- selnet@likelihood
    }
  }

  stopCluster(cl)

  for(i in 1:length(maxeval@nets)) {
    if(is.null(maxeval@nets[[i]]))
      next
    maxeval@nets[[i]]@categories <- categories
  }
  
  t2 <- proc.time()
  maxeval@time <- t2[3] - t1[3]
  
  return(maxeval)
}


cnSearchHist <- function(data, perturbations,  
                         maxParentSet, maxComplexity,
                         parentsPool = NULL, fixedParentsPool = NULL, 
                         niter = 32, echo=FALSE) {

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
  
  if(maxParentSet < 1)
    maxParentSet <- 1
  if(numnodes < 1 || numsamples < 1)
    stop("No valid sample is specified.")

  if(length(nodenames) < numnodes) {
    nodenames <- seq(1, numnodes)
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

  if(fast) {
    ## call the C-function
    .Call("ccnReleaseCache", PACKAGE="catnet")
  }

  norders <- niter
  res <- vector("list", norders)
  
  for(k in 1:norders) {

    order <- sample(1:numnodes)
    if(echo)
      cat("\n Order: ", order, "\n")

    t1 <- proc.time()
    
    res[[k]] <- optimalNetsForOrder(data, perturbations, 
                                    categories, as.integer(maxCategories),
                                    as.integer(maxParentSet), as.integer(maxComplexity), order, 
                                    parentsPool, fixedParentsPool, 
                                    fast, echo, useCache=FALSE)
    
    t2 <- proc.time()
    mins <- floor((t2[1]-t1[1])/60)
    secs <- t2[1]-t1[1] - 60*mins
    if(echo)
      cat(k, "\\", norders, "search time: ", mins, "min, ", secs, "sec\n")

    if(k == 1)
      mhisto <- parHisto(res[[k]], order)
    else
      mhisto <- mhisto + parHisto(res[[k]], order)
 
  }

  rownames(mhisto)<-nodenames
  colnames(mhisto)<-nodenames
  
  return(mhisto)
}


cnSearchHistCluster <- function(data, perturbations,  
                                maxParentSet, maxComplexity,
                                parentsPool = NULL, fixedParentsPool = NULL, 
                                niter = 32,
                                clusterNodes = 2,
                                clusterHost = "localhost",
                                echo = FALSE) {

  if(!require("snow")) {
    stop("Snow not found.")
  }

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame")

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
  
  if(maxParentSet < 1)
    maxParentSet <- 1
  if(numnodes < 1 || numsamples < 1)
    stop("No valid sample is specified.")

  if(length(nodenames) < numnodes) {
    nodenames <- seq(1, numnodes)
  }

  r <- .categorizeSample(data, perturbations)
  data <- r$data
  perturbations <- r$perturbations
  categories <- r$categories
  maxCategories <- r$maxCategories

  if(is.null(perturbations)) {
    dims <- dim(data)
    perturbations <- matrix(rep(0,dims[1]*dims[2]), dims[1], dims[2])
  }

  norders <- niter
  res <- vector("list", norders)

  cl <- makeSOCKcluster(rep(clusterHost, clusterNodes))

  ## create new instances
  maxParentSetTemp <- as.integer(maxParentSet)
  maxComplexityTemp <- as.integer(maxComplexity)
  parentsPoolTemp <- parentsPool
  fixedParentsPoolTemp <- fixedParentsPool
  bFastTemp <- TRUE
  bEchoTemp <- echo
  bUseCache <- FALSE
  
  if(bFastTemp) {
    ## call the C-function
    .Call("ccnReleaseCache", PACKAGE="catnet")
  }

  k <- 0
  while(k < norders) {

    order <- vector("list", clusterNodes)
    for(cc in 1:clusterNodes)
        order[[cc]] <- sample(1:numnodes)

    t1 <- proc.time()
    
    clres <- clusterApply(cl, order, function(ord) {
      return(optimalNetsForOrder(data, perturbations, 
                                 categories=categories,
                                 maxCategories=as.integer(maxCategories),
                                 maxParentSet=maxParentSetTemp,
                                 maxComplexity=maxComplexityTemp,
                                 nodeOrder=ord, 
                                 parentsPool=parentsPoolTemp, fixedParentsPool=fixedParentsPoolTemp, 
                                 bFastTemp, bEchoTemp, bUseCache))
    })

    t2 <- proc.time()
    mins <- floor((t2[3]-t1[3])/60)
    secs <- t2[3]-t1[3] - 60*mins
    if(echo)
      cat(k+clusterNodes, "\\", norders, "search time: ", mins, "min, ", secs, "sec\n")
    
    for(cc in 1:clusterNodes) {
      k <- k + 1
      res[[k]] <- clres[[cc]]
                    
      if(k == 1)
        mhisto <- parHisto(res[[k]], order[[cc]])
      else
        mhisto <- mhisto + parHisto(res[[k]], order[[cc]])
    }
      
  }

  stopCluster(cl)
  
  rownames(mhisto)<-nodenames
  colnames(mhisto)<-nodenames
  
  return(mhisto)
}

