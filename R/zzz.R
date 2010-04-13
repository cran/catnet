.onLoad <- function(lib, pkg, where) {
  dotpath <- as.character(Sys.getenv("R_DOTVIEWER"))
  if(dotpath == "") {
    try(err <- require(igraph), TRUE)
    if(err)
      Sys.setenv(R_CATNET_USE_IGRAPH=TRUE)
    else {
      Sys.setenv(R_CATNET_USE_IGRAPH=FALSE)
      cat("No plotting capabilities detected. Type 'help(catnet)'\n")
    }
  }
  else {
    Sys.setenv(R_CATNET_USE_IGRAPH=FALSE)
  }
  set.seed(02081969)
}

