.onLoad <- function(lib, pkg) {
  dotpath <- as.character(Sys.getenv("R_DOTVIEWER"))
  if(dotpath == "") {
    pl <- .packages(all.available=TRUE)
    err <- FALSE
    for(pp in pl)
      if(pp=="igraph") {
        err <- TRUE
        break
      }
    if(err)
      Sys.setenv(R_CATNET_USE_IGRAPH=TRUE)
    else {
      Sys.setenv(R_CATNET_USE_IGRAPH=FALSE)
      packageStartupMessage("No plotting capabilities detected. Type 'help(catnet)'\n")
    }
  }
  else {
    Sys.setenv(R_CATNET_USE_IGRAPH=FALSE)
  }
  set.seed(02081969)
}

