# function from demo(error.catching) which silently catches a warning and includes it as a member of a returned list
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