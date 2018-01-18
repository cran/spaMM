# ~ try() but catches warnings : 
## Catch *and* save both errors and warnings, and in the case of
## a warning, also keep the computed result.
.tryCatch_W_E <- function(expr) {
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
   				     warning = w.handler),
       warning = W)
}

.minimalErrorHandler <- function(cond) invisible(structure("See 'condition' attribute", class = "try-error", condition = cond))

## meeting the width.cutoff in deparse :-(
.DEPARSE <- function(expr) { paste(deparse(expr),collapse="") }

overcat <- function(msg, prevmsglength) {
  if (prevmsglength>0) {cat("\r")}    
  cat(msg)
  return(nchar(msg))
}

# removed all 'pastefrom' code in version 2.2.44

## quite ad hoc
.prettify_num <- function(x,nsmall=2L) {
  if (abs(x)>1) {
    Ldigits <- floor(log10(abs(x)))+nsmall+1L
    format(round(x,Ldigits), digits=nsmall, nsmall=nsmall) ## character !
  } else signif(x,6) ## double !
}
