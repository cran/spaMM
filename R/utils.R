# ~ try() but catches warnings : 
##' Catch *and* save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
tryCatch.W.E <- function(expr) {
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
   				     warning = w.handler),
       warning = W)
}

## meeting the width.cutoff in deparse :-(
DEPARSE <- function(expr) {
  paste(deparse(expr),collapse="")
}

overcat <- function(msg, prevmsglength) {
  msglength <- nchar(msg)
  if (prevmsglength>0) {cat("\r")}    
  cat(msg)
  return(msglength)
}

