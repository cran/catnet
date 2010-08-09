
#########################################################################
# Categorical Network Class Methods
# Random Sampling

##setMethod("cnSamples", c("catNetwork"),   
##	function(object, numsamples=1, output="frame", as.index=FALSE) {
##	  if(!is.numeric(numsamples) && !is.integer(numsamples))
##            stop("The number of samples should be integer")
##	  numsamples <- as.integer(numsamples)
##	  data <- sapply(1:numsamples, function(x) genRandomCats(object))
##	  data <- matrix(data, ncol=numsamples)
##          if(!as.index) {
##            for(i in 1:object@numnodes) {
##              cats <- object@categories[[i]]
##              for(j in 1:numsamples) {
##                data[i, j] <- cats[as.integer(data[i,j])]
##              }
##            }
##          }
##	  rownames(data)<-object@nodes
##	  if(!missing(output) && output=="matrix")
##            return(data)
##          return(as.data.frame(t(data)))
##	})

setMethod("cnSamples", c("catNetwork"),
	function(object, numsamples=1, perturbations=NULL, output="frame", as.index=FALSE) {
	  if(!is.numeric(numsamples) && !is.integer(numsamples))
            stop("The number of samples should be integer")
	  numsamples <- as.integer(numsamples)
          if(!is.null(perturbations) && is.vector(perturbations))
            data <- sapply(1:numsamples, function(x) x<-genRandomCatsPert(object, perturbations))
          else
            data <- sapply(1:numsamples, function(x) genRandomCats(object))
	  data <- matrix(data, ncol=numsamples)
          rownames(data)<-object@nodes
          if(!as.index) {
            for(i in 1:object@numnodes) {
              cats <- object@categories[[i]]
              for(j in 1:numsamples) {
                data[i, j] <- cats[as.integer(data[i,j])]
              }
            }
            if(missing(output) || output!="matrix")
              data <- as.data.frame(t(data))
          }
          else {
            if(missing(output) || output!="matrix") {
              data <- as.data.frame(t(data))
              for(i in 1:object@numnodes)
                data[,i] <- as.integer(data[,i])
            }
            else {
              for(i in 1:object@numnodes)
                data[i, ] <- as.integer(data[i,])
            }
          }
          return(data)
	})

# recursive probability search
findProbSlot <- function(idroot, ppars, pcatlist, idx, problist, psample) {
  if(is.null(ppars) || length(idx) < 1) {
    return(problist)
  }
  idnode <- ppars[idx[1]]
  if(psample[idnode] < 1 || psample[idnode] > length(pcatlist[[idnode]]))
    stop("Wrong category\n")
  return(findProbSlot(idnode, ppars, pcatlist, idx[-1], problist[[psample[idnode]]], psample))
}

genRandomCats <- function(object) {  
  data <- rep(-1, length(object@nodes))  
  rpath <- orderNodesDescend(object@parents)
  for(i in (1:length(rpath))) {
    ##cat(i,", ")
    nnode <- rpath[i]
    if(length(object@categories[[nnode]]) < 2){
      data[nnode] <- -1
      next
    }
    if(is.null(object@parents[[nnode]])) {
      r <- runif(1,0,1)
      icat <- 1
      rcat <- object@probabilities[[nnode]][icat]
      while(rcat < r) {
        icat <- icat + 1
        rcat <- rcat + object@probabilities[[nnode]][icat]
      }
      data[nnode] <- icat
      next
    }
    plist <- findProbSlot(nnode, object@parents[[nnode]], object@categories,
                          seq(1,length(object@parents[[nnode]])), object@probabilities[[nnode]], data)
    r <- runif(1,0,1)
    icat <- 1
    rcat <- plist[icat]
    while(rcat < r) {
      icat <- icat + 1
      rcat <- rcat + plist[icat]
    }
    data[nnode] <- icat
    ##cat("\n")
  }
  return(data) 
}

genRandomCatsPert <- function(object, perturbations) {
  data <- rep(-1, length(object@nodes))  
  rpath <- orderNodesDescend(object@parents)
  for(i in (1:length(rpath))) {
    ##cat(i,", ")
    nnode <- rpath[i]
    if(perturbations[nnode] >= 1 && perturbations[nnode] <= length(object@categories[[nnode]])) {
      data[nnode] <- perturbations[nnode]
      next
    }
    if(length(object@categories[[nnode]]) < 2){
        data[nnode] <- -1
        next
      }
    if(is.null(object@parents[[nnode]])) {      
      r <- runif(1,0,1)
      icat <- 1
      rcat <- object@probabilities[[nnode]][icat]
      while(rcat < r) {
        icat <- icat + 1
        rcat <- rcat + object@probabilities[[nnode]][icat]
      }
      data[nnode] <- icat
      next
    }
    plist <- findProbSlot(nnode, object@parents[[nnode]], object@categories,
                   seq(1,length(object@parents[[nnode]])), object@probabilities[[nnode]], data)
    r <- runif(1,0,1)
    icat <- 1
    rcat <- plist[icat]
    while(rcat < r) {
      icat <- icat + 1
      rcat <- rcat + plist[icat]
    }
    data[nnode] <- icat
    ##cat("\n")
  }
  return(data) 
}

# creates a tree-list of zeroes
# for a node [idroot] with vector [ppars] of parents and [pcatlist] of categories
# idx is a recursion controlling index set
initSampleProb <- function(idroot, ppars, pcatlist, idx) {
  if(is.null(ppars) || length(idx) < 1) {
    return(rep(0, length(pcatlist[[idroot]])))
  }
  idnode <- ppars[idx[1]]
  lapply(seq(1, length(pcatlist[[idnode]])),
         function(cat, idroot, ppars, pcatlist, idx)
         initSampleProb(idroot, ppars, pcatlist, idx),
         idroot, ppars, pcatlist, idx[-1])
}

# for a node [idroot] with parents [ppars] and
# probability tree-list [problist], 
# update the conditional probability 
# based on a sample [psample] by incrementing the coresponding frequency value
updateSampleProb <- function(idroot, ppars, pcatlist, idx, problist, psample) {
  if(is.null(ppars) || length(idx) < 1) {
    if(psample[idroot] > length(pcatlist[[idroot]]))
      stop("Wrong category\n")
    problist[psample[idroot]] <- problist[psample[idroot]] + 1
    return(problist)
  }
  idnode <- ppars[idx[1]]
  poutlist <- lapply(seq(1,length(pcatlist[[idnode]])),
           function(cat, idroot, ppars, pcatlist, idx, problist, psample) {
             if(psample[idnode] != cat)
               return(problist[[cat]])
             else
               return(updateSampleProb(idroot, ppars, pcatlist, idx[-1], problist[[cat]], psample))
             },
           idroot, ppars, pcatlist, idx, problist, psample)
  return(poutlist)
}

# here [psample] contains the samples categories for the nodes listed in [ppars]
# while [pcounts] gives corresponding frequences
# update is for all categories simultaneously
setNodeSampleProb <- function(idroot, ppars, pcatlist, idx, problist, data) {
  if(is.null(ppars) || length(idx) < 1) {
    if(!is.null(dim(data)) && dim(data)[2] > 0)
      problist <- sapply(1:length(pcatlist[[idroot]]),
                         function(ncat) {                           
                           return(sum(data[idroot,]==ncat))
                         })
    else
      problist <- sapply(1:length(pcatlist[[idroot]]), function(cat) return(0))
    return(problist)
  }
  idnode <- ppars[idx[1]]
  poutlist <- lapply(1:length(pcatlist[[idnode]]),
      function(catidx, idroot, ppars, pcatlist, idx, problist, data, pcounts) {
      if(!is.null(dim(data)) && dim(data)[1] >= idnode){
        psubsample <- data[,data[idnode,] == catidx]
        ## beware: psubsample may bacome a vector
        if(is.null(dim(psubsample)) && length(psubsample) > 0)
          psubsample <- as.matrix(psubsample, rows=length(psubsample))
      }
      else
        psubsample <- NULL
      return(setNodeSampleProb(idroot, ppars, pcatlist, idx[-1], problist[[catidx]], psubsample))
    },
                     idroot, ppars, pcatlist, idx, problist, data)
  return(poutlist)
}

normalizeProb <- function(net) {
  if(!is(net,"catNetwork"))
    stop("net should be catNetwork.")
  if(length(net@parents) < net@numnodes)
      net@parents <- c(net@parents, vector("list", net@numnodes - length(net@parents)))
  if(length(net@parents) != length(net@probabilities))
     stop("length(net@parents) != length(net@probabilities)")
  for(i in 1:length(net@probabilities)) {
    idx <- NULL
    if(length(net@parents[[i]]) > 0)
      idx <- 1:length(net@parents[[i]])
    net@probabilities[[i]] <- normalizeProbSlot(i, net@parents[[i]], net@categories, idx, net@probabilities[[i]])
  }
  return(net)
}

# based on a sample [psample] with parents [ppars] and
# probability tree-list [problist], 
# normalize the conditional probability tree-list at [idroot] so that it sums to 1
normalizeProbSlot <- function(idroot, ppars, pcatlist, idx, problist) {
  if(is.null(ppars) || length(idx) < 1) {
    if(length(problist) != length(pcatlist[[idroot]]) || length(pcatlist[[idroot]]) < 2) {
       return(problist)
    }
    ps <- sum(problist)
    if(!is.na(ps) && is.numeric(ps)) {
      if(ps > 0)
        problist <- problist/ps
      else {
        ## set neutral probability
        nn <- length(problist)
        avg <- 1 / nn
        problist <- rep(avg, nn)
        problist[nn] <- 1 - sum(problist[1:(nn-1)])
      }
    }
    return(problist)
  }
  idnode <- ppars[idx[1]]
  poutlist <- lapply(seq(1,length(pcatlist[[idnode]])),
                     function(cat) {
                       normalizeProbSlot(idroot, ppars, pcatlist, idx[-1], problist[[cat]])
                     })
  return(poutlist)
}

# calculates the likelihood at node [idroot] with parents [ppars] and
# conditional probability tree-list [problist]
# which is assumed to contain non-normalized frequences
# multinomial model is assumed 
nodeCondLoglik <- function(idroot, ppars, pcatlist, idx, problist) {
  if(is.null(ppars) || length(idx) < 1) {
    if(length(pcatlist[[idroot]]) < 2) {
       return(0)
    }
    probs <- problist
    ps <- sum(probs)

    ####################################################
    ## Next condition determines whether parent sets with
    ## non-sample-populated slots should be discarded
    ##if(ps <= 0)
    ##  return(-Inf)
    if(ps > 0) {
      ps <- 1 / ps
      probs <- probs * ps
    }
    ####################################################
    probs[probs==0] <- 1
    return(sum(problist*log(probs)))
  }
  idnode <- ppars[idx[1]]
  poutlist <- sapply(seq(1,length(pcatlist[[idnode]])),
                     function(cat, idroot, ppars, pcatlist, idx, problist)
                     nodeCondLoglik(idroot, ppars, pcatlist, idx[-1], problist[[cat]]),
                     idnode, ppars, pcatlist, idx, problist
                     )
  return(sum(poutlist))
}

.setSampleProb <- function(object, data) {
  if(!is(object, "catNetwork"))
    stop("Object should be catNetwork.")
  
  if(dim(data)[1] != object@numnodes)
    stop("Incompatible sample dimensions.\n")
  numsamples <- dim(data)[2]
  if(numsamples < 1)
    stop("No samples\n")

  for(nnode in (1:object@numnodes)) {
    object@probabilities[[nnode]] <- initSampleProb(nnode, 
                                                object@parents[[nnode]],
                                                object@categories,
                                                seq(1,length(object@parents[[nnode]])))
  }
  
  for(j in (1:numsamples)) {
    ps <- data[,j]
    for(nnode in 1:object@numnodes) {
      ## increment frequency      
      object@probabilities[[nnode]] <- updateSampleProb(nnode, object@parents[[nnode]], object@categories,
                                                        seq(1,length(object@parents[[nnode]])), object@probabilities[[nnode]], ps)
    }
    
    for(nnode in (1:object@numnodes)) {
      object@probabilities[[nnode]] <- normalizeProbSlot(nnode, object@parents[[nnode]], object@categories,
                                                         seq(1,length(object@parents[[nnode]])), object@probabilities[[nnode]])
    }
  }
  
  return(object)
}


setMethod("cnSetProb", "catNetwork",
          function(object, data) {
            
            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame of categories")
            
            r <- .categorizeSample(data, NULL, object)
            data <- r$data
            object@categories <- r$categories
            object@maxCategories <- r$maxCategories
            
            if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
	      stop("The number of nodes in the object and data should be equal.")

	    rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows should be named after the nodes of the object.")

            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              data <- data[norder,]
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }

            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The row names should correspond to the object nodes.")

	    fast <- TRUE

            numnodes <- dim(data)[1]
            numsamples <- dim(data)[2]
            
	    if(fast) {
		## needs to create a probability list
		object@probabilities <- lapply(seq(1, numnodes),
                  function(parid) {
                      setDefaultProb(parid, object@parents[[parid]], object@categories,
                                    seq(1, length(object@parents[[parid]])))
                  })

		newobject <- .Call("ccnSetProb", 
                      	object, data, NULL, 
                      	PACKAGE="catnet")

                ## awkward but necessary
		newobject@nodes <- object@nodes
	    }
	    else {
              newobject <- .setSampleProb(object, data)
            }
            
            newobject@categories <- object@categories
            newobject@maxCategories <- object@maxCategories
            
	    return(newobject)
          }
)

.nodeLikelihood <- function(idroot, ppars, pcatlist, idx, problist, psample) {
  if(is.null(ppars) || length(idx) < 1) {
    ##cat(idroot," ", psample[idroot], "\n")
    if(psample[idroot] > length(pcatlist[[idroot]]))
      stop("Wrong category\n")
    cat <- psample[idroot]
    ##cat(idroot, ": nodeLikelihood ", cat, ", ", problist[cat], "\n")
    return(problist[cat])
  }
  idnode <- ppars[idx[1]]
  if(idnode > length(psample))
    stop("Wrong sample")
  cat <- psample[idnode]
  ##cat(idroot," ", idnode, " cat = ", cat)
  return(.nodeLikelihood(idroot, ppars, pcatlist, idx[-1], problist[[cat]], psample))
}

.networkLikelihood <- function(object, data) {
  if(dim(data)[1] != object@numnodes)
    stop("Incompatible sample dimensions.\n")
  
  numsamples <- dim(data)[2]
  if(numsamples < 1)
    stop("No samples\n")

  r <- .categorizeSample(data, NULL, object)
  data <- r$data
  
  floglik <- 0
  for(j in (1:numsamples)) {
    ps <- data[,j]
    liklist <- sapply(seq(1,object@numnodes), function(nnode) {
      .nodeLikelihood(nnode, object@parents[[nnode]], object@categories, 
                      seq(1,length(object@parents[[nnode]])),
                      object@probabilities[[nnode]], ps)
    })
    ##floglik <- floglik + sum(sapply(liklist, function(x) if(x > 0) log(x) else -Inf))
    floglik <- floglik + sum(sapply(liklist, function(x) if(x > 0) log(x) else 0))
  }

  return(floglik)
}


setMethod("cnLoglik", c("catNetwork"), 
          function(object, data, bysample=FALSE) {

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame of categories")

            if(is.data.frame(data))
              data <- as.matrix(t(data))

            if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
              stop("The number of nodes in  the object and data should be equal")
            
            rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows/cols should be named after the nodes of the object.")
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              ## keep the matrix format !!
              data <- as.matrix(data[norder,])
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The data names should correspond to the object nodes.")

            r <- .categorizeSample(data, NULL, object, ask=FALSE)
            data <- r$data
            if(object@maxCategories < r$maxCategories)
              stop("Data has more categories than the object")
            
            ##if(missing(fast) || is.null(fast)) 
            fast <- TRUE
            
            numnodes <- dim(data)[1]
            numsamples <- dim(data)[2]
            
            if(fast) {
              loglik <- .Call("ccnLoglik", 
                              object, data, NULL, 
                              PACKAGE="catnet")
              if(bysample)
                return(loglik)
              return(sum(loglik))
            }
            else
              return(.networkLikelihood(object, data))
          })


.nodeSampleLoglik <- function(nnode, parentSet, data, categories, maxCategories) {

  numnodes <- dim(data)[1]
  numsamples <- dim(data)[2]

  nodeprob <- initSampleProb(nnode, parentSet, categories, seq(1,length(parentSet)))
  for(j in (1:numsamples)) {
    ps <- data[,j]
    ## increment frequency      
    nodeprob <- updateSampleProb(nnode, parentSet, categories,
                                 seq(1,length(parentSet)), nodeprob, ps)
  }

  nodeprob <- normalizeProbSlot(nnode, parentSet, categories,
                   seq(1,length(parentSet)), nodeprob)

  # calculate likelihood (a kind of redundancy)
  floglik <- 0
  for(j in (1:numsamples)) {
    ps <- data[,j]
    lik <- .nodeLikelihood(nnode, parentSet, categories,
                      seq(1,length(parentSet)),
                      nodeprob, ps)
    floglik <- floglik + log(lik)
  }
  return(floglik)
}

cnNodeSampleLoglik <- function(node, parents, data, perturbations = NULL) {

  if(!is.matrix(data) && !is.data.frame(data))
    stop("data should be a matrix or data frame of node categories")

  if(is.matrix(data)) {
    numnodes <- dim(data)[1]
    numsamples <- dim(data)[2]  
  }
  else {
    numnodes <- dim(data)[2]
    numsamples <- dim(data)[1]  
  }
  
  if(numsamples < 1)
    stop("No samples\n")

  r <- .categorizeSample(data, perturbations)
  data <- r$data
  categories <- r$categories
  maxCategories <- r$maxCategories
  
  if(!is.list(node) && !is.list(parents)) {
    if(node < 1 || node > numnodes)
      stop("Incompatible sample dimensions.\n")
    if(!is.null(perturbations)) {
      data <- data[, perturbations[node,]==0]
    }  
    return(.nodeSampleLoglik(node, parents, data, categories, maxCategories))
  }

  if(is.list(node) && is.list(parents)) {
    if(length(node) != length(parents))
      stop("'node' and 'parents' lists should be of equal length")
    loglik <- rep(0, length(node))
    for(i in 1:length(node)) {
      nod <- node[[i]]
      par <- parents[[i]]
      ##cat("nod=", nod, ", par=", par, "\n")
      if(nod < 1 || nod > numnodes)
        stop("Incompatible sample dimensions.\n")
      if(!is.null(perturbations)) {
        subdata <- data[, perturbations[node,]==0]
        loglik[i] <- .nodeSampleLoglik(nod, par, subdata, categories, maxCategories)
      }
      else {
        loglik[i] <- .nodeSampleLoglik(nod, par, data, categories, maxCategories)
      }
    }
    return(loglik)
  }

  return(0)
}


setMethod("cnNodeLoglik", c("catNetwork"), 
          function(object, node, data, perturbations=NULL) {

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame of categories")

            if(is.data.frame(data))
              data <- as.matrix(t(data))

            if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
              stop("The number of nodes in  the object and data should be equal")
            
            rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows/cols should be named after the nodes of the object.")
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              ## keep the matrix format !!
              data <- as.matrix(data[norder,])
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The data names should correspond to the object nodes.")

            if(!is.numeric(node) || node < 1 || node > object@nodes)
              stop("Wrong node")
            node <- as.integer(node)
              
            r <- .categorizeSample(data, NULL, object, ask=FALSE)
            data <- r$data
            if(object@maxCategories < r$maxCategories)
              stop("Data has more categories than the object")
            
            fast <- TRUE
            
            numnodes <- dim(data)[1]
            numsamples <- dim(data)[2]
            
            if(fast) {
              loglik <- .Call("ccnNodeLoglik", 
                              object, node, data, NULL, 
                              PACKAGE="catnet")
              return(loglik)
            }
            else
              return(.nodeLikelihood(node,
                                     object@parents[[node]], object@categories, 
                                     seq(1,length(object@parents[[node]])),
                                     object@probabilities[[node]], data))
          })


setMethod("cnNodeLoglikError", c("catNetwork"), 
          function(object, node, data, perturbations=NULL) {

            if(!is.matrix(data) && !is.data.frame(data))
              stop("'data' should be a matrix or data frame of categories")

            if(is.data.frame(data))
              data <- as.matrix(t(data))

            if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
              stop("The number of nodes in  the object and data should be equal")
            
            rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows/cols should be named after the nodes of the object.")
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              ## keep the matrix format !!
              data <- as.matrix(data[norder,])
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }
            
            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The data names should correspond to the object nodes.")

            if(!is.numeric(node) || node < 1 || node > object@nodes)
              stop("Wrong node")
            node <- as.integer(node)
              
            r <- .categorizeSample(data, NULL, object, ask=FALSE)
            data <- r$data
            if(object@maxCategories < r$maxCategories)
              stop("Data has more categories than the object")
            
            numnodes <- dim(data)[1]
            numsamples <- dim(data)[2]
            
            loglik <- .Call("ccnNodeLoglikError", 
                            object, node, data, NULL, 
                            PACKAGE="catnet")
            return(loglik)
          })
