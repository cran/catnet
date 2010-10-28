#########################################################################
# Categorical Network Class Methods  
 
## generate a catNetwork object 
setMethod("initialize", "catNetwork",  
          function(.Object, ...) { 
            if(!validObject(.Object)) 
              stop("Error creating the object.") 
            .Object@objectName <- "catNetwork" 
            return(.Object) 
          }) 
 
setMethod("initialize", "catNetwork",  
          function(.Object, numnodes, maxParents = 0, numCategories = 2, p.default=FALSE, p.delta1=0.01, p.delta2=0.01, ... ) { 
 
            .Object@objectName <- "catNetwork" 
 
            if(nargs() >= 2 && is.numeric(numnodes) && is.numeric(maxParents) && is.numeric(numCategories)) 
              return(genRandomCatnet(.Object, numnodes, maxParents, numCategories, p.default, p.delta1, p.delta2)) 
 
            return(edges2catnet(.Object, c(1), edges = vector("list", 1))) 
          }) 
 
cnRandomCatnet <- function(numnodes, maxParents, numCategories, p.delta1=0.01, p.delta2=0.01) { 
  return(new("catNetwork", as.integer(numnodes), as.integer(maxParents), as.integer(numCategories),
             as.logical(FALSE), p.delta1, p.delta2 )) 
}
 
cnCatnetFromGraph <- function(graph, numCategories = 2, cats = NULL, probs = NULL) { 
  maxCategories <- numCategories 
  for(cat in cats) 
    if(maxCategories < length(cat)) 
      maxCategories <- length(cat) 
  object <- new("catNetwork", numnodes=1, maxParents=0, maxCategories=2) 
  object <- graph2catnet(object, graph, as.integer(maxCategories), cats, probs) 
  return(object) 
} 
 
cnCatnetFromEdges <- function(nodes, edges, numCategories = 2) { 
  numnodes <- length(nodes) 
  if(numnodes < 1) 
    stop("'At least one node should be specified") 
  if(length(edges) != numnodes) 
    stop("'edges' should be a list with the length of 'nodes'") 
  if(numCategories < 2) { 
    numCategories <- 2 
    warning("'numCategories' is set to 2") 
  } 
   
  maxParents <- 0 
  for(nodeedges in edges) 
    if(maxParents < length(nodeedges)) 
      maxParents <- length(nodeedges) 
   
  object <- new("catNetwork", numnodes=numnodes, as.numeric(maxParents), as.integer(numCategories)) 
  object <- edges2catnet(object, nodes, edges, as.integer(numCategories)) 
  return(object) 
} 
 
cnCatnetFromSif <- function(file, numCategories = 2) { 
  fsif <- read.table(file) 
  nodes <- NULL 
  for(s in fsif[,1]) 
    if(length(which(nodes==s)) > 0) 
      next 
    else 
      nodes <- c(nodes, s) 
  for(s in fsif[,3]) 
    if(length(which(nodes==s)) > 0) 
      next 
    else 
      nodes <- c(nodes, s) 
  numnodes <- length(nodes) 
  edges <- vector("list", numnodes) 
  for(n in 1:dim(fsif)[1]) { 
    k <- which(nodes==fsif[n,1]) 
    ##m <- which(nodes==fsif[n,3]) 
    m <- as.character(fsif[n,3]) 
    if(length(which(edges[[k]]==m)) > 0) 
      next 
    edges[[k]] <- c(edges[[k]], m) 
  } 
 
  object <- new("catNetwork", numnodes=length(nodes)) 
  object <- edges2catnet(object, nodes, edges, as.integer(numCategories)) 
  return(object) 
} 
 
cnNew <- function(nodes, cats, parents, probs = NULL, p.delta1=0.01, p.delta2=0.01) {
 
  if(length(nodes) < 1) 
    stop("At least 1 node is needed") 
 
  if(length(nodes) != length(parents)) 
    stop("Incompatible parent list") 
 
  if(length(nodes) != length(cats)) 
    stop("Incompatible category list") 
   
  maxCategories <- 0 
  for(cat in cats) if(!is.null(cat) && maxCategories < length(cat)) maxCategories <- length(cat) 
 
  maxParents <- 0
  i <- 1
  for(par in parents) {
    if(!is.null(par) && maxParents < length(par)) {
      maxParents <- length(par)
      parents[[i]] <- as.integer(par)
    }
    i <- i + 1
  }
 
  object <- new("catNetwork", length(nodes), maxParents, maxCategories) 
  object@nodes <- nodes 
  object@parents <- parents 
  object@categories <- cats 
  if(!is.null(probs)) { 
    if(length(nodes) != length(probs)) 
      stop("Incompatible probability list") 
    object@probabilities <- probs 
  }
  else {
    object@probabilities <-
      lapply(seq(1, object@numnodes), function(parid) {
        if(length(object@parents[[parid]]) > 0)
          setRandomProb(parid, object@parents[[parid]], object@categories,
                      seq(1, length(object@parents[[parid]])), p.delta1, p.delta2)
        else
          setRandomProb(parid, object@parents[[parid]], object@categories, NULL, p.delta1, p.delta2)
      })
  }

  object@nodeComplexity <- sapply(1:object@numnodes, function(x) nodeComplexity(object, x)) 
  object@complexity <- as.integer(sum(object@nodeComplexity)) 
   
  object@likelihood <- 0 
  object@nodeLikelihood <- NA
  
  if(!validCatNetwork(object, TRUE)) 
    stop("Incompatible parameters") 
   
  return(object)   
} 
 
setMethod("cnAddNode", c("catNetwork", "character", "character", "ANY", "ANY"), 
          function(object, node, cats, parents = NULL, prob = NULL) { 
 
            object@numnodes <- as.integer(object@numnodes + 1) 
             
            if(length(which(object@nodes == node))>0) { 
              node <- paste("N", object@numnodes ,sep="") 
            } 
            object@nodes <- c(object@nodes, node) 
             
            if(length(cats) < 2) { 
              warning("The node should have at least two categories") 
              cats <- paste("C", 1:2 ,sep="") 
            } 
            object@categories[[object@numnodes]] <- cats 
            if(length(cats) > object@maxCategories) 
              object@maxCategories <-  length(cats) 
 
            if(is.null(parents)) { 
              warning("Set random parents") 
              idx <- sample(1:(object@numnodes-1)) 
              parents <- idx[1:object@maxParents] 
            } 
            object@parents[[object@numnodes]] <- parents 
            if(length(parents) > object@maxParents) 
               object@maxParents <-  length(parents) 
 
            if(is.null(prob)) { 
              warning("Set default conditional probability") 
              prob <- setDefaultProb(object@numnodes, parents, object@categories, 1:length(parents)) 
            } 
            object@probabilities[[object@numnodes]] <- prob 
 
            return(object) 
            }) 
 
setMethod("show", "catNetwork", 
          function(object) { 
            if(is(object, "catNetwork")) 
              cat("A catNetwork object with ", object@numnodes, " nodes, ", 
                  object@maxParents, " parents, ", object@maxCategories, 
                  " categories,\n Likelihood = ", object@likelihood, 
                  ", Complexity = ", object@complexity, ".\n") 
            }) 
 
setMethod("cnNumNodes", "catNetwork", function(object) length(object@nodes)) 
 
genRandomCatnet <- function(.Object, numnodes, maxparents, maxcats, defaultprob, p.delta1, p.delta2) { 
 
  nodes <- sapply(seq(1,numnodes), function(i) paste("N", i ,sep=""))             
  parents <- genRandomParents(numnodes, maxparents) 
   
  cats <- lapply(seq(1:numnodes), function(i, maxcats) 
                 paste("C", seq(1:maxcats), sep=""), 
                 maxcats) 

  if(p.delta1<0) p.delta1 <- 0
  if(p.delta2<0) p.delta2 <- 0
  while(p.delta1+p.delta2>=0.5) {
    p.delta1 <- p.delta1/2
    p.delta2 <- p.delta2/2
  }
  
  probs <- lapply(seq(1, numnodes), 
                  function(parid) { 
                    if(defaultprob) 
                      setDefaultProb(parid, parents[[parid]], cats, 
                                    seq(1, length(parents[[parid]]))) 
                    else 
                      setRandomProb(parid, parents[[parid]], cats, 
                                    seq(1, length(parents[[parid]])), p.delta1, p.delta2) 
                  }) 
 
  .Object@numnodes <- as.integer(numnodes) 
  .Object@nodes <- nodes 
  .Object@meta <- "" 
  .Object@maxParents <- as.integer(maxparents) 
  .Object@parents <- parents 
  .Object@categories <- cats 
  .Object@maxCategories <- as.integer(maxcats) 
  .Object@probabilities <- probs 
 
  .Object@nodeComplexity <- sapply(1:as.integer(numnodes), function(x) nodeComplexity(.Object, x)) 
  .Object@complexity <- as.integer(sum(.Object@nodeComplexity)) 
   
  .Object@likelihood <- 0 
  .Object@nodeLikelihood <- NA 
   
  return(.Object)   
} 
 
# returns vector   
setMethod("cnNodes", c("catNetwork", "missing"),  
          function(object) { 
            which <- seq(1:object@numnodes) 
            cnNodes(object, which) 
          }) 

setMethod("cnNodes", c("catNetwork", "vector"),  
          function(object, which) { 
            if(is.null(which)) 
              which <- seq(1:object@numnodes) 
            as.character(object@nodes[which]) 
          }) 
 
setMethod("cnEdges", c("catNetwork", "missing"),  
          function(object) { 
            which <- seq(1:object@numnodes) 
            cnEdges(object, which) 
          }) 

setMethod("cnEdges", c("catNetwork", "character"),  
          function(object, which) { 
            id <- sapply(which, function(node) which(object@nodes == node))
            if(length(id) < 1)
              return(NULL)
            cnEdges(object, id) 
          })

setMethod("cnEdges", c("catNetwork", "vector"),  
	function(object, which) { 
          if(object@numnodes < 1) 
            return(NULL) 
          if(is.null(which)) 
            which <- seq(1, object@numnodes) 
	  onodes <- cnNodes(object)   
	  plist <- vector("list", length(which)) 
          for(j in (1:length(which))) 
            names(plist)[[j]] <- onodes[which[j]] 
	  for(j in (1:object@numnodes)) { 
	    parlist <- object@parents[[j]]   
	    if(length(parlist) < 1) next 
	    for(i in (1:length(parlist))) {   
		idx <- which(which == parlist[i]) 
                if(length(idx) < 1) next 
                for(jj in 1:length(idx)) 
                  plist[[idx[jj]]] <- c(plist[[idx[jj]]], onodes[j]) # letter nodes    
	    }  
	  } 
          if(length(plist)==0) 
            return(plist) 
          i <- 1 
          while(i<=length(plist)) 
            if(is.null(plist[[i]])) 
              plist <- plist[-i] 
            else 
              i <- i + 1 
          return(plist) 
	})   
   
 
setMethod("cnMatEdges", "catNetwork",  
          function(object) { 
            mout <- NULL 
            if(object@numnodes < 1) 
              return(mout) 
            for(n in 1:object@numnodes) { 
              if(length(object@parents[[n]]) < 1) 
                next 
              for(j in object@parents[[n]]) { 
                mout <- c(mout, object@nodes[j], object@nodes[n]) 
              } 
            } 
            if(is.null(mout)) 
              return(NULL) 
            return(matrix(mout, nc=2, byrow=TRUE)) 
          }) 
 
# returns list 
setMethod("cnParents", c("catNetwork", "missing"),  
          function(object) { 
            which <- seq(1:object@numnodes) 
            cnParents(object, which) 
          }) 

setMethod("cnParents", c("catNetwork", "character"),  
          function(object, which) { 
            id <- sapply(which, function(node) which(object@nodes == node))
            if(length(id) < 1)
              return(NULL)
            cnParents(object, id) 
          })
 
setMethod("cnParents", c("catNetwork", "vector"),  
          function(object, which) { 
            if(is.null(which)) 
              which <- seq(1, object@numnodes) 
            plist <- lapply(which, function(n) { 
              if(is.null(object@parents[[n]])) 
                return(NULL) 
              as.character(object@nodes[object@parents[[n]]]) 
            })
            if(length(plist)==0) 
              return(plist)            
            if(!is.list(plist))
              return(as.vector(plist))
            pnames <- object@nodes[which]
            plist <- setNames(plist, pnames) 
            i <- 1 
            while(i<=length(plist)) 
              if(is.null(plist[[i]])) 
                plist <- plist[-i] 
              else 
                i <- i + 1 
            return(plist) 
          }) 
 
matParents <- function(object, nodeorder) { 
            n <- object@numnodes 
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
 
## generates the binary matrix of pairwise node directed connections 
setMethod("cnMatParents", c("catNetwork", "missing"),  
          function(object) { 
            n <- object@numnodes 
            nodeorder <- seq(1,n) 
            return(matParents(object, nodeorder)) 
          })

setMethod("cnMatParents", c("catNetwork", "vector"),  
          function(object, nodeorder) { 
            n <- object@numnodes 
            if(missing(nodeorder) || is.null(nodeorder)) 
              nodeorder <- seq(1,n) 
            return(matParents(object, nodeorder)) 
          }) 
 
listProbSet <- function(idroot, ppars, pcatlist, idx, problist, strin) { 
  if(is.null(ppars) || length(idx) < 1) { 
    if(length(pcatlist[[idroot]]) != length(problist)) {
      ##cat(idroot, ",   ", length(pcatlist[[idroot]]), ", ", length(problist), "\n")
      warning("Wrong probability slot") 
      return("") 
    } 
    strout <- sapply(seq(1, length(pcatlist[[idroot]])), function(j, strin, pcatlist, problist) { 
      paste("[", strin, "]", pcatlist[j], "  ", problist[j], "\n", sep="") 
    }, strin, pcatlist[[idroot]], problist) 
    return(paste(strout, sep="", collapse="")) 
  } 
  idnode <- ppars[idx[1]] 
  idx <- idx[-1] 
  strout <- sapply(seq(1,length(pcatlist[[idnode]])), 
    function(cat, strin) { 
      listProbSet(idroot, ppars, pcatlist, idx, problist[[cat]], paste(strin, pcatlist[[idnode]][cat], sep=" ")) 
    }, strin) 
  return(paste(strout, sep="", collapse="")) 
} 
 
 
setMethod("cnProb", c("catNetwork", "missing"),  
          function(object) { 
            which <- seq(1:object@numnodes) 
            cnProb(object, which) 
          }) 

setMethod("cnProb", c("catNetwork", "character"),  
          function(object, which) {
            id <- sapply(which, function(node) which(object@nodes == node))
            if(length(id) < 1)
              return(NULL)
            cnProb(object, id) 
          })

setMethod("cnProb", c("catNetwork", "vector"),  
          function(object, which) { 
            if(is.null(which)) 
              which <- seq(1, object@numnodes) 
            str <- sapply(which, function(n) { 
              ##paste(paste("Node[", n, "], Parents: ", sep=""), 
              ##             paste(object@parents[[n]], collapse=" "), "\n",  
              ##             listProbSet(n, object@parents[[n]], object@categories, 
              ##                           seq(1,length(object@parents[[n]])), 
              ##                           object@probabilities[[n]], ""), 
              ##             collapse="", sep="") 
              paste(paste("Node[", object@nodes[n], "], Parents: ", sep=""), 
                    paste(object@nodes[object@parents[[n]]], sep=",", collapse=" "), "\n",  
                    listProbSet(n, object@parents[[n]], object@categories, 
                                seq(1,length(object@parents[[n]])), 
                                object@probabilities[[n]], ""), 
                    collapse="", sep="") 
            }) 
            cat(str)
          }) 
 
setMethod("cnOrder", "catNetwork",    
	function(object) { 
          cnOrderNodes(object@parents) 
	})  
 
setMethod("cnOrder", "list",    
	function(object) { 
          cnOrderNodes(object) 
	}) 
 
validCatNetwork <- function(obj, quietly=FALSE) {
  res = TRUE
  onodes <- obj@nodes   
  nnodes <- length(onodes)   
  opars <- obj@parents   
  omaxpars <- obj@maxParents   
  ocats <- obj@categories   
  omaxcats <- obj@maxCategories   
  oprob <- obj@probabilities   
  if(nnodes < 1 || omaxcats < 2) {   
    if(nnodes < 1 && !quietly)   
      cat("At least one node is needed.\n")   
    if(omaxcats < 2 && !quietly)   
      cat("At least two categories are needed.\n")   
    res <- FALSE   
    return(res)   
  }   
  for(i in (1:length(onodes))) {   
    if(length(opars[[i]]) > omaxpars){   
      omaxpars <- length(opars[[i]])   
    }   
    if(length(ocats[[i]]) > omaxcats){   
      omaxcats <- length(ocats[[i]])   
    }   
    while(length(opars[[i]]) == 0 && i < length(opars))   
      i <- i + 1   
    if(length(opars[[i]]) == 0 && i >= length(opars))   
      break 
    nn <- 1 
    for(j in (1:length(opars[[i]]))) {   
      #cat(i, length(opars[[i]]),"\n")   
      ipar <- opars[[i]][j]   
      if(ipar < 1 || ipar > nnodes) {   
         if(!quietly)   
           cat("Wrong parent set: node=", i,", parent=", j, ".\n")   
        res <- FALSE   
        return(res)   
      }   
      #Check DAG condition; Note that the parent set is an increasing sequence   
      if(ipar == i)  {   
        if(!quietly)   
          cat("Wrong network, not a DAG: self-reference for node ", i, ".\n")   
        res <- FALSE   
        return(res)   
      }   
      if(ipar < i && length(opars[[ipar]])) {   
        for(ii in (1:length(opars[[ipar]])) ) {   
          jj <- opars[[ipar]][ii]   
          if(jj == i) {   
            if(!quietly)   
              cat("Wrong network, not a DAG: nodes ", i, " and ", ipar, " are mutually connected.\n")   
            res <- FALSE   
            return(res)   
          }   
          if(jj > i)   
            break;   
        }   
      }   
      nn <- nn*length(ocats[[j]]) 
    } 
  } 
 
  for(j in 1:nnodes) { 
    if(!checkProbSet(j, opars[[j]], ocats, seq(1,length(opars[[j]])), oprob[[j]])) { 
      if(!quietly)   
        cat("Wrong probability for node ", j, "\n")   
      res <- FALSE   
      return(res)   
    }   
  } 
 
  # update 
  obj@maxParents <- omaxpars   
  obj@maxCategories <- omaxcats 

  if(!isDAG(onodes, opars)) {   
    if(!quietly)   
      cat("Wrong network: not a DAG.\n")   
    res <- FALSE   
    return(res)   
  } 
   
  return(res)    
}   
 
nodeComplexity <- function(object, nnode) { 
  ll <- sapply(object@parents[[nnode]], function(i) length(object@categories[[i]])) 
  if(length(ll)>0) 
    return(prod(ll)*(length(object@categories[[nnode]])-1)) 
  else 
    return(length(object@categories[[nnode]])-1) 
} 
 
parentSetComplexity <- function(numnodecats, parents, categories) { 
  ll <- sapply(parents, function(i) length(categories[[i]])) 
  if(length(ll)>0) 
    return(prod(ll)*(numnodecats-1)) 
  else 
    return(numnodecats-1) 
} 
 
setMethod("cnComplexity", signature("catNetwork"), function(object, node) {
  if(missing(node)) { 
    pc <- sapply(1:object@numnodes, function(x) nodeComplexity(object, x)) 
    return(as.integer(sum(pc))) 
  }
  if(is.character(node))
    node <- which(object@nodes == node)
  if(!is.numeric(node)) { 
    pc <- sapply(1:object@numnodes, function(x) nodeComplexity(object, x)) 
    return(as.integer(sum(pc))) 
  } 
  return(as.integer(nodeComplexity(object, as.integer(node))))   
})  

setMethod("cnSubNetwork", "catNetwork",  
function(object, nodeIndices, indirectEdges = FALSE) { 

  if(is.character(nodeIndices)) {
    nodeIndices <- sapply(nodeIndices, function(node) which(object@nodes == node))
    if(length(nodeIndices) < 1)
      return(object)
  }
  
  nodeIndices <- nodeIndices[nodeIndices<=object@numnodes] 
  nodes <- object@nodes[nodeIndices] 
  numnodes <- length(nodes) 
  if(numnodes < 1) 
    stop("Can't create a network without nodes.") 
 
  parents <- vector("list", numnodes)   
  categories <- vector("list", numnodes) 
   
  maxParents <- 0 
  maxCategories <- 0 
  for(j in 1:numnodes) { 
    i <- nodeIndices[j] 
    categories[[j]] <- object@categories[[i]] 
    if(indirectEdges==TRUE) { 
      nodepars <- object@parents[[i]] 
      while(!is.null(nodepars)) { 
        pool <- NULL 
        ##cat(nodepars, "\n") 
        for(par in nodepars) { 
          id <-  which(nodeIndices==par) 
          if(length(id) > 0 && 
             ## prevent repeating 
             length(which(parents[[j]] == id)) <= 0) 
            parents[[j]] <- c(parents[[j]], id) 
          else 
            pool <- c(pool, object@parents[[par]]) 
        } 
        nodepars <- pool 
      } 
    } 
    else { 
      for(par in object@parents[[i]]) { 
        id <- which(nodeIndices==par) 
        if(length(id) > 0) 
          parents[[j]] <- c(parents[[j]], id) 
      } 
    } 
    if(maxParents < length(parents[[j]])) 
      maxParents <- length(parents[[j]]) 
    if(maxCategories < length(categories[[j]])) 
      maxCategories <- length(categories[[j]])
  } 
   
  newnet <- new("catNetwork", numnodes, as.integer(maxParents), as.integer(maxCategories)) 
  newnet@nodes <- nodes 
  newnet@parents <- parents 
  newnet@categories <- categories 
  problist <- lapply(1:numnodes, function(parid, obj) { 
    setDefaultProb(parid, parents[[parid]], categories, 1:length(obj@parents[[parid]])) 
  }, newnet) 
  newnet@probabilities <- problist 
  
  pc <- sapply(1:newnet@numnodes, function(x) nodeComplexity(newnet, x)) 
  newnet@complexity <- as.integer(sum(pc))
      
  return(newnet) 
}) 
 
 
setProbSlot <- function(idroot, ppars, pcatlist, idx, problist, catvec, probslot) { 
  if(is.null(ppars) || length(idx) < 1) { 
    if(length(problist) != length(probslot)) 
      stop("length(problist) != length(probslot)") 
    problist <- probslot 
    return(problist) 
  } 
  idnode <- ppars[idx[1]] 
  cat <- catvec[idx[1]] 
  if(cat < 1 || cat > length(pcatlist[[idnode]])) 
    stop("Wrong category\n") 
  return(setProbSlot(idnode, ppars, pcatlist, idx[-1], problist[[cat]], catvec, probslot)) 
} 
 
 
reorderNodeProb <- function(idroot, ppars, categories, catvec, idx, prob, 
                            nodeIndices, newparents, newcategories) { 
  if(is.null(ppars) || length(idx) < 1) { 
    newroot <- nodeIndices[idroot] 
    cat("newroot = ", newroot, "\n") 
    newidx <- NULL 
    if(length(newparents) > 0) 
      newidx <- seq(1, length(newparents)) 
    ##cat("newidx = ", newidx, "\n") 
    newprob <- initSampleProb(newroot, newparents, newcategories, newidx) 
    ##catvec <- catvec[nodeIndices] 
    ##cat("catvec = ", catvec, "\n") 
    ##cat("len(prob) = ", length(prob), "\n") 
    newprob <- setProbSlot(newroot, newparents, newcategories, 
                           newidx, newprob, catvec, prob) 
    return(newprob) 
  } 
  idnode <- ppars[idx[1]] 
  lapply(seq(1, length(categories[[idnode]])), 
         function(cat, idx) 
         reorderNodeProb(idroot, ppars, categories, c(catvec, cat), idx, prob[[cat]], 
                         nodeIndices, newparents, newcategories), 
         idx[-1]) 
} 
 
setMethod("cnReorderNodes", c("catNetwork", "vector"),  
function(object, nodeIndices) {

  if(is.character(nodeIndices)) {
    nodeIndices <- sapply(nodeIndices, function(node) which(object@nodes == node))
    if(length(nodeIndices) < 1)
      return(object)
  }

  nodeIndices <- nodeIndices[nodeIndices<=object@numnodes] 
  if(length(nodeIndices) != object@numnodes) { 
    warning("length(nodeIndices) != object@numnodes") 
    return(object) 
  } 
  ##cat(nodeIndices, "\n") 
  nodeIndicesInvert <- nodeIndices 
  for(i in 1:object@numnodes) { 
    nodeIndicesInvert[i] = which(nodeIndices == i) 
  } 
 
  ##cat(nodeIndicesInvert, "\n") 
   
  newnodes <- object@nodes[nodeIndices] 
  parents <- vector("list", object@numnodes) 
  categories <- vector("list", object@numnodes) 
  probabilities <- vector("list", object@numnodes) 
 
  for(i in 1:object@numnodes) { 
    if(length(object@parents[[nodeIndices[i]]]) > 0) {
      parents[[i]] <- nodeIndicesInvert[object@parents[[nodeIndices[i]]]]
      names(parents[[i]]) <- newnodes[parents[[i]]]
    }
    categories[[i]] <- object@categories[[nodeIndices[i]]] 
  }

  for(i in 1:object@numnodes) { 
    ##probabilities[[nodeIndices[i]]] <- reorderNodeProb(i, object@parents[[i]], object@categories, 
    ##                                      catvec=NULL, 
    ##                                      seq(1, length(object@parents[[i]])), 
    ##                                      object@probabilities[[i]], 
    ##                                      nodeIndices, parents, categories) 
    probabilities[[i]] <- object@probabilities[[nodeIndices[i]]] 
  } 
 
  object@nodes <- newnodes 
  object@categories <- categories 
  object@parents <- parents 
  object@probabilities <- probabilities 
   
  return(object) 
}) 

cnSetSeed <- function(seed) {
  ## set the seed both in R and in the C-library
  .Call("ccnSetSeed", seed, PACKAGE="catnet")
  set.seed(seed)
}
