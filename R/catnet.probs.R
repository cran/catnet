#########################################################################
# Categorical Network Class Methods
# Probability Calculations

checkProbSet <- function(idroot, ppars, pcatlist, idx, problist) {
  if(length(ppars) == 0 || length(ppars) < 1 || length(idx) < 1) {
    if(is.null(problist) || length(pcatlist[[idroot]]) != length(problist)) {
      return(FALSE)
    }
    if(!is.nan(sum(problist)) || sum(problist) == -Inf)
      return(TRUE)
    if(sum(problist<0) > 0) {
      warning("Probability slot with negative values", "\n")
      return(FALSE)
    }
    if(as.integer(sum(problist)) > 1) {
      warning("Probability slot with sum ", sum(problist), "\n")
      return(FALSE)
    }
    return(TRUE)
  }
  idnode <- ppars[idx[1]]
  res <- sapply(seq(1,length(pcatlist[[idnode]])),
                function(cat)
    checkProbSet(idroot, ppars, pcatlist, idx[-1], problist[[cat]]))
  if(sum(res) < length(pcatlist[[idnode]]))
    return(FALSE)
  return(TRUE)
}

## Attach a new leaf-node with a parent set [leafpars] and conditional probability [leafproblist]
## to a tree-list [ptree] that includes all nodes less in order than the new node
## THE PARENTS IN [leafpars] HAS THE SAME ORDER AS THAT IN [ptree],
## the latter is assured if the catNetwork is created by genRandomParents

# recursive probability assignment
setRandomProb <- function(idroot, ppars, pcatlist, idx) {
  if(is.null(ppars) || length(ppars) < 1 || length(idx) < 1) {
    if(length(pcatlist[[idroot]]) < 2) {
      return(NULL)
    }
    plist <- sapply(pcatlist[[idroot]], function(x) runif(1,0.01,0.99))
    plist <- plist/sum(plist)
    plist <- sapply(seq(1,length(plist)), function(n, plist){
      plist[n] <- floor(100*plist[n])/100
      }, plist)
    #cat(sum(plist),"\n")
    plist[1] <- 1 - sum(plist[-1])
    #poutlist <- c(poutlist, plist)
    return(as.vector(plist))
  }
  else {
    id <- ppars[idx[1]]
    #cat(ppars[idx], id, "\n")
    poutlist <- lapply(pcatlist[[id]],
           function(cat, idroot, ppars, pcatlist, idx)
           setRandomProb(idroot, ppars, pcatlist, idx),
           idroot, ppars, pcatlist, idx[-1])
  }
}

# recursive probability assignment
setDefaultProb <- function(idroot, ppars, pcatlist, idx) {
  if(is.null(ppars) || length(ppars) < 1 || length(idx) < 1) {
    if(length(pcatlist[[idroot]]) < 2) {
      return(NULL)
    }
    plist <- rep(1, length(pcatlist[[idroot]]))
    plist <- plist/sum(plist)
    plist <- sapply(seq(1,length(plist)), function(n, plist){
      plist[n] <- floor(100*plist[n])/100
      }, plist)
    plist[1] <- 1 - sum(plist[-1])
    return(as.vector(plist))
  }
  else {
    id <- ppars[idx[1]]
    ##cat(id, ": ",idx, "\n")
    poutlist <- lapply(pcatlist[[id]],
           function(cat, idroot, ppars, pcatlist, idx)
           setDefaultProb(idroot, ppars, pcatlist, idx),
           idroot, ppars, pcatlist, idx[-1])
  }
}
