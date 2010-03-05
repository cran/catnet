#########################################################################
# Categorical Network Class Methods
# Find-members

setMethod("cnFind", "catNetworkEvaluate",
          function(object, complexity) {
            if(!is(object, "catNetworkEvaluate"))
               stop("catNetworkEvaluate object should be specified.")
            if(length(object@nets) < 1)
              return(NULL)
            complx <- sapply(object@nets, function(pnet) pnet@complexity)
            idx <- which(complx == complexity)
            if(length(idx) == 0) {
              for(i in 1:length(complx))
                if(complx[i] > complexity)
                  break
              idx <- i
            }
            return(object@nets[[idx]])
            })

setMethod("cnFind", "list",
          function(object, complexity) {
            idx <- NULL
            idx <- sapply(object, function(pnet) {
              if(is(pnet, "catNetwork"))
                return(pnet@complexity==complexity)
              else
                return(FALSE)
            })
            if(is.null(idx))
              return(NULL)
            if(sum(idx) == 0) {
              idx <- sapply(object, function(pnet) {
                if(is(pnet, "catNetwork"))
                  return(abs(cnComplexity(pnet)-complexity))
              else
                return(Inf)
              })
              cc <- idx
              idx[cc==min(cc)] <- TRUE
              idx[cc!=min(cc)] <- FALSE              
            }
            idx <- which(idx==max(idx))
            id <- max(idx)
            return(object[[id]])
            })

setMethod("cnFindAIC", "catNetworkEvaluate", function(object) {
  objectlist <- object@nets
  liststr <- ""
  maxobj <- NULL
  maxaic <- -Inf
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curaic <- object@likelihood - object@complexity
    if(maxaic < curaic) {
      maxaic <- curaic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("cnFindAIC", "list", function(object) {
  objectlist <- object
  liststr <- ""
  maxobj <- NULL
  maxaic <- -Inf
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curaic <- object@likelihood - object@complexity
    if(maxaic < curaic) {
      maxaic <- curaic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("cnFindBIC", "catNetworkEvaluate", function(object, numsamples) {
  objectlist <- object@nets
  liststr <- ""
  maxobj <- NULL
  maxbic <- -Inf
  numsamples <- object@numsamples
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curbic <- object@likelihood - 0.5*object@complexity*log(numsamples)
    if(maxbic < curbic) {
      maxbic <- curbic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("cnFindBIC", "list", function(object, numsamples) {
  objectlist <- object
  liststr <- ""
  maxobj <- NULL
  maxbic <- -Inf
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curbic <- object@likelihood - 0.5*object@complexity*log(numsamples)
    if(maxbic < curbic) {
      maxbic <- curbic
      maxobj <- object
    }
  }
  return(maxobj)
})
