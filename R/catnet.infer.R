#########################################################################
# Categorical Network Class Methods
# Inference

setMethod("addNetworkNode", "catNetwork",
          function(object, newnode, newparents, loglik, problist, nodecats) {
            object@nodes <- c(object@nodes, newnode)
            object@categories[[length(object@nodes)]] <- nodecats
            object@parents[[length(object@nodes)]] <- newparents
            object@probabilities[[length(object@nodes)]] <- problist
            if(object@maxParents < length(newparents))
              object@maxParents <- length(newparents)
            if(object@maxCategories < length(nodecats))
              object@maxCategories <- length(nodecats)
            object@numnodes <- length(object@nodes)
            c <- 1
            pc <- sapply(newparents, function(j) length(object@categories[[j]]))
            if(length(pc) > 0)
              c <- prod(pc)
            nodeComplx <- as.integer((length(nodecats)-1)*c)
            object@nodeComplexity <- c(object@nodeComplexity, nodeComplx)
            object@complexity <- object@complexity + nodeComplx
            object@nodeLikelihood <- c(object@nodeLikelihood, loglik)
            object@likelihood <- object@likelihood + loglik
            return(object)
            })

setMethod("replaceNetworkNode", "catNetwork",
          function(object, newnode, newparents, loglik, problist) {
            if(length(which(object@nodes == newnode)) < 1)
              return(object)
            object@parents[[newnode]] <- newparents
            object@probabilities[[newnode]] <- problist
            if(object@maxParents < length(newparents))
              object@maxParents <- length(newparents)
            excomplexity <- object@complexity - object@nodeComplexity[newnode]
            object@likelihood <- object@likelihood - object@nodeLikelihood[newnode]
            c <- 1
            pc <- sapply(newparents, function(j) length(object@categories[[j]]))
            if(length(pc) > 0)
              c <- prod(pc)
            object@nodeComplexity[newnode] <- as.integer((length(object@categories[[newnode]])-1)*c)
            object@complexity <- as.integer(excomplexity + object@nodeComplexity[newnode])
            object@nodeLikelihood[newnode] <- loglik
            object@likelihood <- object@likelihood + loglik
            return(object)
            })

nodeMatParents <- function(object, nodeorder) {
            n <- length(object@nodes)
            mat <- matrix(rep(0, n*n), n, n)
            for(j in 1:length(object@parents)) {
              plist <- object@parents[[j]]
              if(length(plist) == 0)
                next
              for(i in plist) {
                mat[nodeorder[j], nodeorder[i]] <- 1
              }
            }
            rownames(mat) <- object@nodes
            colnames(mat) <- object@nodes
            return(mat)
            }

setMethod("cnMatParents", c("catNetwork", "missing"), 
          function(object) {
            n <- length(object@nodes)
            nodeorder <- seq(1,n)
            return(nodeMatParents(object, nodeorder))
            })

setMethod("cnMatParents", c("catNetwork", "vector"), 
          function(object, nodeorder) {
            n <- length(object@nodes)
            if(missing(nodeorder) || is.null(nodeorder))
              nodeorder <- seq(1,n)
            return(nodeMatParents(object, nodeorder))
            })

# in the list of networks [listnet] 
# for each complexity keep one network,
# the network with the highest likelihood 
setMethod("updateNetworkNode", "catNetwork",
          function(newnet, listnet) {
            # assume node-nets are kept in listnet and indexed by their complexity
            c <- newnet@complexity
            if(length(listnet) < c)
              return(listnet)
            #cat(newnet@complexity,"\n")
            pnode <- listnet[[c]]
            
            if(is.null(pnode) || pnode@likelihood < newnet@likelihood) {
              listnet[[c]] <- newnet
            }
            
            return(listnet)
          })

.categorizeSample <- function(data, perturbations, object = NULL) {

  if(is.matrix(data)) {
    numnodes <- nrow(data)
    numsamples <- ncol(data)

    ## save row/column names
    samplenames <- dimnames(data)
    if(is.numeric(data)) {
      data <- matrix(as.integer(data), nrow=dim(data)[1])
    }
    
    maxCategories <- 1
    categories <- vector("list", numnodes)
    for(i in 1:numnodes) {
      cats <- NULL
      for(c in data[i,]) {
        if(is.na(c))
          next
        if(sum(cats==c)==0)
          cats <- c(cats,c)
      }
      ## include a dummy
      if(length(cats) <= 1)
        cats <- c(cats, -1)
      cats <- sort(cats)
      categories[[i]] <- cats
      if(maxCategories < length(cats))
        maxCategories <- length(cats)
      
      for(j in 1:numsamples) {
        c <- which(cats == data[i,j])
        data[i,j] <- c
      }
    }

    data <- matrix(as.integer(data), nrow=numnodes)
    dimnames(data) <- samplenames

    if(!is.null(object) && object@numnodes == numnodes) {
      ## make sure data and object's nodes have the same order
      rnames <- rownames(data)
      tsamp <- data
      if(!is.null(rnames)) {
        for(nnode in 1:object@numnodes) {
          i <- which(object@nodes == rnames[nnode])
          if(length(i) == 1) {
            if(nnode != i)
              data[i,] <- tsamp[nnode,]
          }
          else
            stop("Object's nodes and sample nodes are different")
        }
        data <- matrix(as.integer(data), nrow=numnodes)
        rownames(data) <- object@nodes
      }
    }
  } ## is.matrix
  else {
    ## data is data.frame format
    numnodes <- ncol(data)
    numsamples <- nrow(data)
    
    fdata <- data
    data <- matrix(rep(NA, numnodes*numsamples), nrow=numnodes)
    
    maxCategories <- 1
    categories <- vector("list", numnodes)
    for(i in 1:numnodes) {
      cats <- levels(fdata[,i])
      lencat <- length(cats)
      if(lencat > 0) {
        data[i,] <- as.integer(fdata[,i])
      }
      else {
        cats <- NULL
        for(c in fdata[,i]) {
          if(is.na(c))
            next
          if(sum(cats==c)==0)
          cats <- c(cats,c)
        }
        cats <- sort(cats)
        lencat <- length(cats)
        for(j in 1:numsamples) {
          c <- which(cats == fdata[j,i])
          if(length(c) == 1)
            data[i,j] <- c
        }
      }
      categories[[i]] <- cats
      
      if(maxCategories < lencat)
        maxCategories <- lencat
      
    }

    rownames(data) <- colnames(fdata)
    colnames(data) <- rownames(fdata)
  }

  if(length(rownames(data)) < numnodes)
    rownames(data) <- 1:numnodes

  if(!is.null(perturbations) && !is.matrix(perturbations)) {
    perturbations <- as.matrix(t(perturbations))
  }

  if(maxCategories*maxCategories > numsamples) {
    cat("The small sample. Continue? ('y' or 'n')\n")
    if(scan("", what="character", nmax=1, quiet=TRUE) != "y" ) 
      stop("Operation canceled")
  }

  if(maxCategories > 16) {
    cat("The data seems to have too many categories. The operation can be very long and memory consuming. Continue? ('y' or 'n')\n")
    if(scan("", what="character", nmax=1, quiet=TRUE) != "y" ) 
      stop("Operation canceled")
  }
  
  return(list(data=data, categories=categories, maxCategories=maxCategories, perturbations=perturbations))
}

# lists all (k=parsize)-combinations from a parent set (parset)
combinationSets <- function(outlist, curlist, parset, allowedset, parsize) {
  if(parsize < 1)
    return(outlist)
  n <- length(curlist)
  if(n > 0)
    ancestor <- curlist[[n]]
  else
    ancestor <- 0
  if(parsize == 1) {
    nout <- length(outlist)
    if(nout == 0)
      outlist <- vector("list", 1)
    nout <- nout + 1
    for(i in (1:length(parset))) {
      if(length(which(allowedset == parset[i])) <= 0)
        next
      if(parset[i] >= ancestor) {
        outlist[[nout]] <- c(curlist, parset[i])
        nout <- nout + 1
      }
    }
  }
  else {
    for(i in (1:length(parset))) {
      if(length(which(allowedset == parset[i])) <= 0)
        next
      if(parset[i] >= ancestor)
        outlist <- combinationSets(outlist, c(curlist, parset[i]), parset[-i], allowedset, parsize-1)
    }
  }
  return(outlist)
}

# internal functions
optimalNetsForOrder <- function(
           data, perturbations, 
           categories, maxCategories, 
           maxParentSet, maxComplexity, nodeOrder,
           parentsPool = NULL, fixedParentsPool = NULL,  
           fast = FALSE, echo=TRUE, useCache=TRUE) {
    
  ## [perturbations] is a binary matrix of the dimensions of data, 
  ## value 1 means the sample node is perturbed, so ignore its likelihood

  numnodes <- dim(data)[1]
  numsamples <- dim(data)[2] 
  nodenames <- seq(1,numnodes)
  if(length(dim(data)) == 2) {
    nodenames <- dimnames(data)[[1]]
  }

  if(fast) {
    bestnets <- .Call("ccnOptimalNetsForOrder", 
                      data, perturbations, 
                      maxParentSet, maxComplexity, nodeOrder,
                      parentsPool, fixedParentsPool, useCache, 
                      PACKAGE="catnet")
    if(length(nodenames) == numnodes && length(bestnets) > 0) {
      for(i in 1:length(bestnets))
	if(is.null(bestnets[[i]]))
		stop("Failed to find a network")
	else
        	bestnets[[i]]@nodes <- nodenames[nodeOrder]
    }
    return(bestnets)
  }
  
  if(!is.null(fixedParentsPool)) {
    for(fixpars in fixedParentsPool) {
      if(length(fixpars) > maxParentSet) {
        maxParentSet <- length(fixpars)
      }
    } 
  }

  if(echo)
    cat("nodeOrder: ", nodeOrder,"\n")

  ## Add the first node in order
  problist <- initSampleProb(nodeOrder[1], NULL, categories, NULL)
  if(!is.null(perturbations)) {
    psubsamples <- data[,perturbations[nodeOrder[1],]==0]
    problist <- setNodeSampleProb(nodeOrder[1], NULL, categories,
                                  NULL, problist, psubsamples)
  }
  else
    problist <- setNodeSampleProb(nodeOrder[1], NULL, categories,
                                  NULL, problist, data)
  loglik <- nodeCondLoglik(nodeOrder[1], NULL, categories, NULL, problist)

  netnode <- new("catNetwork", numnodes=1, maxParents=0, maxCategories)
  netnode@nodes[[1]] <- nodeOrder[1]
  ##netnode@parents <- vector("list", 1)
  netnode@categories[1] <- categories[nodeOrder[1]]

  netnode@probabilities[[1]] <- problist
  netnode@nodeLikelihood[1] <- loglik
  netnode@likelihood <- loglik
  
   ## the list of current optimal networks indexed by their complexity
  listnet <- vector("list", 1)
  listnet[[netnode@complexity]] <- netnode
  
  for(i in (2:length(nodeOrder))) {

    nnode <- nodeOrder[i]
    
    nodidx <- 1:(i-1)
    if(i > 1 && !is.null(parentsPool)) {
      exclude <- NULL
      for(nn in 1:(i-1))
        if(length(which(parentsPool[[nnode]]==nodeOrder[nn])) <= 0)
          exclude <- c(exclude, nn)
      if(!is.null(exclude))
        nodidx <- nodidx[-exclude] 
    }

    ## Next only works if fixedParentsPool is subset of parentsPool
    fixedIdPars <- NULL
    if(!is.null(fixedParentsPool) && !is.null(fixedParentsPool[[nnode]])) {
      exclude <- NULL
      for(nn in 1:length(nodidx)) {
        if(length(which(fixedParentsPool[[nnode]]==nodeOrder[nodidx[nn]])) > 0) {
            exclude <- c(exclude, nn)
          }
      }
      if(!is.null(exclude)) {
        fixedIdPars <- nodidx[exclude]
        nodidx <- nodidx[-exclude]
      }
    }
    
    ## keep best networks with nodes nodidx[1:(i-1)] in oldlistnet
    oldlistnet <-  listnet
    listnet <- vector("list", maxComplexity)

    for(d in (0:min(length(nodidx), maxParentSet - length(fixedIdPars)))) {
  
      if(d == 0) {
        comblist <- vector("list", 1)
        if(!is.null(fixedIdPars))
          comblist[[1]] <- fixedIdPars
      }
      else {
        comblist <- combinationSets(NULL, NULL, 1:(i-1), nodidx, d)
        for(nn in 1:length(comblist)) {
          comblist[[nn]] <- c(fixedIdPars, comblist[[nn]])
        }
      }
      
      if(echo)
        cat(paste("Node ", nnode, " parents = ", d + length(fixedIdPars),
                ##" nets = ", length(oldlistnet),
                "; Parse ", length(comblist), " combination(s).\n", sep=""))

      sprogress <- 1
      progress <- 0
      
      ## all parent set in [comblist] impose the same complexity on [nnode]
      ## thus we find the one with maximal likelihood
      opt.loglik <- -Inf
      opt.idparlist <- NULL
      opt.problist <- NULL
      for(idparlist in comblist) {

        if(is.null(idparlist))
          parlist <- NULL
        else 
          parlist <- nodeOrder[idparlist]


        problist <-
          initSampleProb(nnode, parlist, categories, seq(1, length(parlist)))
        if(!is.null(perturbations)) {
          psubsamples <- data[, perturbations[nnode,]==0]
          if(is.null(dim(psubsamples)) && length(psubsamples) > 0) {
            psubsamples <- as.matrix(psubsamples, rows=length(psubsamples))
          }
          problist <- setNodeSampleProb(nnode, parlist, categories,
                                        seq(1,length(parlist)), problist, psubsamples)
        }
        else
          problist <- setNodeSampleProb(nnode, parlist, categories,
                                        seq(1,length(parlist)), problist, data)
        
        loglik <- nodeCondLoglik(nnode, parlist, categories, seq(1,length(parlist)), problist)          
        
        if(loglik > opt.loglik) {
          opt.loglik <- loglik
          opt.idparlist <- idparlist
          opt.problist <- problist
        }        
      }

      for(pnet in oldlistnet) {
        if(is.null(pnet))
          next
        newnet <- addNetworkNode(pnet, nnode, opt.idparlist, opt.loglik, opt.problist, categories[[nnode]])
        if(newnet@complexity > maxComplexity)
          next
        listnet <- updateNetworkNode(newnet, listnet)
      } # for(pnet in oldlistnet)
    } # for d
  } # for i

  ## compactify listnet by removing empty slots
  j <- 0
  for(i in (1:length(listnet)))
    if(!is.null(listnet[[i]]))
      j <- j + 1
  newlist <- vector("list", j)
  j <- 1
  for(i in (1:length(listnet))) {
    if(!is.null(listnet[[i]])) {
      newlist[[j]] <- listnet[[i]]
      j <- j + 1
    }
  }

  # normalize probabilities
  if(echo)
    cat("Normalization...")
  bestnets <- newlist
  if(length(bestnets) > 0)
    for(i in 1:length(bestnets))
      bestnets[[i]] <- normalizeProb(bestnets[[i]])
  if(echo) {
    cat(" Done\n")
    cat(length(bestnets), " networks found. \n")
  }

  if(length(nodenames) == numnodes && length(bestnets) > 0)
    for(i in 1:length(bestnets))
      bestnets[[i]]@nodes <- nodenames[nodeOrder]
  
  return(bestnets)
}


findOptimalNetworks <- function(object, data, perturbations = NULL, maxParentSet = 0, maxComplexity = 0, fast = FALSE, echo=FALSE) {

  r <- .categorizeSample(data, perturbations, object)
  data <- r$data

  dims <- dim(data)
  numnodes <- dims[1]
  
  if(is.null(perturbations))
    perturbations <- matrix(rep(0,dims[1]*dims[2]), dims[1], dims[2])
  if(maxParentSet <= 0)
    maxParentSet <- object@maxParents
  if(maxComplexity < numnodes)
    maxComplexity <- cnComplexity(object)
  
  pertdims <- dim(perturbations)
  if(dims[1] != object@numnodes)
    stop("Incompatible sample dimensions.\n")
  if(dims[1] != pertdims[1] ||
     dims[2] != pertdims[2])
    stop("Incompatible perturbation dimensions.\n")
  numsamples <- dims[2]
  if(numsamples < object@maxCategories)
    stop(paste("No enough samples. At least ", object@maxCategories," has to be given.\n"))   

  t1 <- proc.time()
  if(echo) {
    cat("Proc time: ", t1,"\n")  
  }

  nodidx <- cnOrder(object)
  ## reorder the nodes so object can be compared later
  object <- cnReorderNodes(object, nodidx)

  bestnets <- optimalNetsForOrder(data, perturbations,
                                  object@categories, object@maxCategories,
                                  maxParentSet, maxComplexity, nodidx, 
                                  parentsPool = NULL, fixedParentsPool = NULL, 
                                  fast, echo, useCache=FALSE)

  t2 <- proc.time()
  if(echo) {
    cat("Proc time: ", t2,"\n")
    mins <- floor((t2[1]-t1[1])/60)
    secs <- t2[1]-t1[1] - 60*mins
    cat("Time elapsed: ", mins, "min, ", secs, "sec\n")
  }
  
  out <- findNetworkDistances(object, numsamples, bestnets, extended = FALSE)

  out@time <- t2[1] - t1[1]

  return(out)
}

setMethod("cnEvaluate", c("catNetwork"), 
          function(object, data, perturbations = NULL, maxParentSet = 0, maxComplexity = 0) {

	    if(!is.matrix(data) && !is.data.frame(data))
              stop("data should be a matrix or data frame of node categories")

            return(findOptimalNetworks(object, data, perturbations, maxParentSet, maxComplexity, fast = TRUE, echo = FALSE))
            })

