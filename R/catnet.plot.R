#########################################################################
# Plotting Methods  

setMethod("cnDot", "catNetwork", function(object, file="", format="ps") { 
  if(length(object@meta)>0) 
    str <- sprintf("\"%s, \\nComplexity %d, \\nLogLikelihood %5.3f\"[shape=plaintext]", 
                   as.character(object@meta), object@complexity, object@likelihood) 
  else 
    str <- sprintf("\"catNetwork with \\nComplexity %d, \\nLogLikelihood %5.3f\"[shape=plaintext]", 
                   object@complexity, object@likelihood)
  noedges <- TRUE
  strout <- sapply(seq(1, length(object@parents)), function(n) { 
    if(is.null(object@parents[[n]])) 
      return("") 
    else{
      noedges <- FALSE
      #cat(n, length(object@parents[[n]]), "\n") 
      paste(sapply(object@parents[[n]], function(j) { 
        #cat(j) 
        if(length(object@parents[[j]]) > 0 && length(which(object@parents[[j]] == n)) > 0) { 
          ## double-edge in both directions 
          ##paste("\"", object@nodes[j], "\" -> \"", object@nodes[n], "\" [dir=both, style=dashed];\n", collapse="", sep="")
          paste("\"", object@nodes[j], "\" -> \"", object@nodes[n], "\" [style=dashed];\n", collapse="", sep="") 
        } 
        else 
          paste("\"", object@nodes[j], "\" -> \"", object@nodes[n], "\";\n", collapse="", sep="") 
      }), collapse="", sep="") 
    } 
  }) 
  strout <- paste(str, paste(strout, collapse="", sep="")) 
  str <- paste("digraph G {\n", strout, "}\n", collapse="", sep="") 
  if(!missing(file) && !is.null(file)) { 
 
    ## get the full path to the file 
    file <- paste(getwd(), "/", file, sep="") 
 
    write(str, file=paste(file,".dot",sep="")) 
 
    dotviewer <- as.character(Sys.getenv("R_DOTVIEWER")) 
    if(dotviewer != "") {
      if(format == "ps") 
        strdotcall<-paste(dotviewer, " -Tps \"", file, ".dot\"", " -o \"", file, ".ps\"", sep="")
      else  
        strdotcall<-paste(dotviewer, " -Tpdf \"", file, ".dot\"", " -o \"", file, ".pdf\"", sep="") 
      try(system(strdotcall, intern=TRUE, ignore.stderr=TRUE), silent = TRUE) 
 
      pdfviewer <- as.character(Sys.getenv("R_PDFVIEWER")) 
      if(pdfviewer != "") {
        if(format == "ps") 
          strevincecall<-paste(pdfviewer, " \"", file, ".ps\"", sep="")
        else
          strevincecall<-paste(pdfviewer, " \"", file, ".pdf\"", sep="") 
        try(system(strevincecall, intern=FALSE, wait=FALSE, ignore.stderr=TRUE), silent = TRUE) 
      } 
    } 
  } 
  else 
    cat(str) 
  }) 
 
 
setMethod("cnDot", "list", function(object, file="", format="ps") { 
  if(!is.list(object)) 
    return("") 
  objectlist <- object 
  liststr <- "" 
  i <- 1 
  for(object in objectlist) {

    if(is(object, "catNetwork")) {
      
      str <- sprintf("\"%s, \\nComplexity %d, \\nLogLikelihood %5.3f\"[shape=plaintext]", 
                     as.character(object@meta), object@complexity, object@likelihood) 
      strout <- sapply(seq(1, length(object@parents)), function(n) { 
        if(is.null(object@parents[[n]])) {
          warning("network without edges")
          return("")
        }
        else{ 
          paste(sapply(object@parents[[n]], function(j) { 
            if(length(object@parents[[j]]) > 0 && length(which(object@parents[[j]] == n)) > 0) { 
              paste("\"", object@nodes[j], "\" -> \"", object@nodes[n], "\" [style=dashed];\n", collapse="", sep="") 
            } 
            else 
              paste("\"", object@nodes[j], "\" -> \"", object@nodes[n], "\";\n", collapse="", sep="") 
          }), collapse="", sep="") 
        } 
      }) 
      strout <- paste(str, paste(strout, collapse="", sep="")) 
      str <- paste("digraph ", sprintf("G%d", i), "{\n", strout, "};\n", collapse="", sep="")
    } ## catNetwork

    if(is.matrix(object)) {
      
      medges <- as.matrix(object) 
      if(dim(medges)[1] != dim(medges)[2] || dim(medges)[1] < 2) {
        warning("Wrong matrix")
        next
      }
      rnames <- rownames(medges) 
      if(is.null(rnames)) 
        rnames <- 1:dim(medges)[1] 
      nnodes <- dim(medges)[1] 
      strout <- "" 
      for(row in 1:nnodes) { 
        for(col in 1:nnodes) { 
          if(medges[row,col] <= 0) 
            next 
          if(medges[col,row] > 0)  
            strout <- paste(strout, rnames[col], " -> ", rnames[row], " [style=dashed];\n", collapse="", sep="") 
          else 
            strout <- paste(strout, rnames[col], " -> ", rnames[row], ";\n", collapse="", sep="") 
        } 
      }
      str <- paste("digraph G {\n", strout, "}\n", collapse="", sep="")
    }
    
    liststr <- paste(liststr, str, "", sep="") 
    i <- i + 1 
  }
  
  if(!missing(file) && !is.null(file)) { 
 
    ## get the full path to the file 
    file <- paste(getwd(), "/", file, sep="") 
 
    write(liststr, file=paste(file,".dot",sep="")) 
 
    dotviewer <- as.character(Sys.getenv("R_DOTVIEWER")) 
    if(dotviewer != "") {
      if(format == "ps") 
        strdotcall<-paste(dotviewer, " -Tps \"", file, ".dot\"", " -o \"", file, ".ps\"", sep="")
      else
        strdotcall<-paste(dotviewer, " -Tpdf \"", file, ".dot\"", " -o \"", file, ".pdf\"", sep="") 
      try(system(strdotcall, intern=TRUE, ignore.stderr=TRUE), silent = TRUE) 
     
      pdfviewer <- as.character(Sys.getenv("R_PDFVIEWER")) 
      if(pdfviewer != "") {
        if(format == "ps") 
          strevincecall<-paste(pdfviewer, " \"", file, ".ps\"", sep="")
        else
          strevincecall<-paste(pdfviewer, " \"", file, ".pdf\"", sep="") 
        try(system(strevincecall, intern=TRUE, ignore.stderr=TRUE), silent = TRUE) 
      } 
    } 
  } 
  else 
    cat(liststr) 
}) 
 
setMethod("cnDot", "matrix", function(object, file="", format="ps") { 
  if(!is(object, "matrix")) 
    stop("Specify a valid square matrix.") 
  medges <- as.matrix(object) 
  if(dim(medges)[1] != dim(medges)[2] || dim(medges)[1] < 2) 
    stop("Specify a valid square matrix.") 
  rnames <- rownames(medges) 
  if(is.null(rnames)) 
    rnames <- 1:dim(medges)[1] 
  nnodes <- dim(medges)[1] 
  strout <- "" 
  for(row in 1:nnodes) { 
    for(col in 1:nnodes) { 
      if(medges[row,col] <= 0) 
        next
      ##cat(rnames[col], " -> ", rnames[row], "\n") 
      if(medges[col,row] > 0)  
        ## double-edge in both directions 
        ##strout <- paste(strout, rnames[row], " -> ", rnames[col], " [dir=both, style=dashed];\n", collapse="", sep="")
        strout <- paste(strout, rnames[col], " -> ", rnames[row], " [style=dashed];\n", collapse="", sep="") 
      else 
        strout <- paste(strout, rnames[col], " -> ", rnames[row], ";\n", collapse="", sep="") 
    } 
  }
  str <- paste("digraph G {\n", strout, "}\n", collapse="", sep="")
  
  if(!missing(file) && !is.null(file)) { 
 
    ## get the full path to the file 
    file <- paste(getwd(), "/", file, sep="") 
 
    write(str, file=paste(file,".dot",sep="")) 
 
    dotviewer <- as.character(Sys.getenv("R_DOTVIEWER")) 
    if(dotviewer != "") { 
      if(format == "ps")
        strdotcall<-paste(dotviewer, " -Tps \"", file, ".dot\"", " -o \"", file, ".ps\"", sep="")
      else
        strdotcall<-paste(dotviewer, " -Tpdf \"", file, ".dot\"", " -o \"", file, ".pdf\"", sep="")
      try(system(strdotcall, intern=TRUE, ignore.stderr=TRUE), silent = TRUE) 
     
      pdfviewer <- as.character(Sys.getenv("R_PDFVIEWER")) 
      if(pdfviewer != "") {
        if(format == "ps") 
          strevincecall<-paste(pdfviewer, " \"", file, ".ps\"", sep="")
        else
          strevincecall<-paste(pdfviewer, " \"", file, ".pdf\"", sep="")
        try(system(strevincecall, intern=TRUE, ignore.stderr=TRUE)) 
      } 
    } 
  } 
  else 
    cat(str) 
  }) 
 
setMethod("cnPlot", "catNetwork", 
        function(object, file = NULL) { 
	    err <- FALSE  
            err <- as.logical(Sys.getenv("R_CATNET_USE_IGRAPH"))
            if(err) { 
              err <- FALSE 
              try(err <- require(igraph), TRUE)
              if(err) { 
                medges <- cnMatEdges(object) 
                if(is.null(medges)) 
                  err <- FALSE 
                else { 
                  igr <- graph.edgelist(medges)
                  nodenames <- get.vertex.attribute(igr, "name")

                  caps <- capabilities()
                  id <- which(names(caps)=="tcltk")
                  usetcltk <- FALSE
                  if(length(id)>0)
                    usetcltk <- caps[id]
                  if(usetcltk)
                    tkplot(igr, vertex.label=nodenames)
                  else
                    plot(igr, vertex.label=nodenames) 
                } 
              }
            }
            if(!err) { 
              if(is.null(file) || file == "") 
                return(cnDot(object, "unknown", "pdf")) 
              else 
                return(cnDot(object, file, "pdf")) 
            } 
          }) 

