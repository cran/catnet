#########################################################################
# Categorical Network Class Methods
# Joint Probability Calculations

probTreeAddLeaf <- function(ptree, leafproblist, leafpars, idtree, idleafpars, pcatlist) {
  if(length(idtree) < 1) {
    ##cat("idtree: ", idtree,"\n")
    ##cat("idleafpars: ", idleafpars,"\n")
    if(length(idleafpars)>0){
      cat(idleafpars)
      stop("Length(idleafpars) should be zero.")
    }
    ## tree is a scalar, leafproblist is a vector
    if(is.null(ptree))
       return(leafproblist)
    return(ptree*leafproblist)
  }
  treenode <- idtree[1]
  if(length(idleafpars) > 0)
    leafnode <- leafpars[idleafpars[1]]
  else
    leafnode <- 0
  ##cat("Nodes: ", idtree,"; ", leafnode, idleafpars, "\n")
  if(length(idtree) > 0 && treenode == leafnode)
    poutlist <- lapply(seq(1, length(pcatlist[[leafnode]])),
                       function(cat, ptree, leapproblist, leafpars, idtree, idleafpars, pcatlist)
                       probTreeAddLeaf(ptree[[cat]], leafproblist[[cat]], leafpars, idtree, idleafpars, pcatlist), 
                       ptree, leafproblist, leafpars, idtree[-1], idleafpars[-1], pcatlist
                       )
  else
    poutlist <- lapply(seq(1, length(pcatlist[[treenode]])),
                       function(cat, ptree, leapproblist, leafpars, idtree, idleafpars, pcatlist)
                       probTreeAddLeaf(ptree[[cat]], leafproblist, leafpars, idtree, idleafpars, pcatlist), 
                       ptree, leafproblist, leafpars, idtree[-1], idleafpars, pcatlist
                       )
  return(poutlist)
}


## 
probTreeToMatrix <- function(ptree, idx, pcatlist, prob, offset) {
  if(length(idx) < 1) {
    prob[offset] <- ptree
    return(prob)
  }
  idnode <- idx[1]
  nodeoff <- 1
  for(i in 1:length(pcatlist)) ## number of nodes == length(pcatlist)
    if(i > idnode)
      nodeoff <- nodeoff*length(pcatlist[[i]])
  ##cat(idnode, nodeoff, length(pcatlist), "\n")
  for(cat in 1:length(pcatlist[[idnode]])) {
    off <- offset + nodeoff*(cat-1)
    prob <- probTreeToMatrix(ptree[[cat]], idx[-1], pcatlist, prob, off)
  }
  return(prob)
}


.jointProb <- function(numnodes, nodes, parents, probabilities, categories) {
  ##cat(nodes, "\n")
  nodesOrder <- orderNodesDescend(parents)
  while(length(nodesOrder) < numnodes) {
    i <- 1
    while(sum(nodesOrder==i)>0)
      i <- i + 1
    nodesOrder <- c(nodesOrder, i)
  }
  if(length(parents) < numnodes) {
    parents <- c(parents, vector("list", numnodes-length(parents)))
  }
  ##cat("\nJointProb" , numnodes, length(parents), nodesOrder, "\n")

  proc.time()
  
  for(i in 1:numnodes) {
    nnode <- nodesOrder[i]
    ptree <- NULL
    idtree <- NULL
    if(i > 1)
      idtree <- nodesOrder[1:(i-1)]
    idleafpars <- NULL
    if(length(parents[[nnode]]) > 0)
      idleafpars <- 1:length(parents[[nnode]])
    ##cat("add ", i, ":", nnode, " " ,parents[[nnode]], "\n")
    ptree <- probTreeAddLeaf(ptree,
                              probabilities[[nnode]],
                              parents[[nnode]],
                              idtree,
                              idleafpars, 
                              categories)
  }

  ## list the joint prob in a matrix
  n <- 1
  for(i in 1:numnodes)
    n <- n*length(categories[[i]])
  prob <- rep(0, n)
  prob <- probTreeToMatrix(ptree, nodesOrder, categories, prob, 1)
  proc.time()
  return(prob)
}


setMethod("cnJointProb", "catNetwork",
          function(object) {
            if(!is(object, "catNetwork"))
              stop("catNetwork object is required.")
            return(.jointProb(object@numnodes, object@nodes, object@parents, object@probabilities, object@categories))
          })

setMethod("cnJointKLdist", "catNetwork",
          function(object1, object2,...) {
            if(!is(object1, "catNetwork") || !is(object2, "catNetwork"))
              stop("catNetwork object is required.")
            if(object1@numnodes != object2@numnodes)
              stop("Number of nodes should be equal.")
            p1 <- .jointProb(object1)
            p2 <- .jointProb(object2)
            probs <- p1
            probs[p2==0] <- 0
            p2[p2==0] <- 1
            p1[p1==0] <- 1
            return(sum(probs*log(p1/p2)))
          })
