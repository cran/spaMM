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
# Details of ?deparse imply that parse/deparse operations are not completely reversible. 
.DEPARSE <- function(expr,collapse="") { paste(deparse(expr),collapse=collapse) }
# parse(text=paste(deparse(binomial()), collapse="")) fails
# parse(text=paste(deparse(binomial()), collapse="\n")) works
# "\n" not OK for formulas; 
# From v 2.5.7, .stripFormula() (now .preprocess_formula()) no longer uses .DEPARSE() [and I neglected stats::reformulate]

overcat <- function(msg, prevmsglength) {
  if (prevmsglength) {cat("\r")}   # \b effect is terminal-dependent so cannot be relied upon.  
  cat(msg)
  return(nchar(msg))
}

# removed all 'pastefrom' code in version 2.2.44

.diagfast <- function(x, nc=ncol(x)) {
  diagPos <- seq.int(1L,nc^2,nc+1L)
  x[diagPos]
}

## quite ad hoc
.prettify_num <- function(x,nsmall=2L) {
  if (abs(x)>1) {
    Ldigits <- floor(log10(abs(x)))+nsmall+1L
    format(round(x,Ldigits), digits=nsmall, nsmall=nsmall) ## character !
  } else signif(x,6) ## double !
}

.prompt <- function() {
  cat("\nPause: 'c' to continue, 'b' for browsing, 's' to stop")
  res <- readLines(n = 1)
  if (res=='s') stop("User-requested stop().")
  if (res=='b') browser()
}

# must work on S4 objects but we avoid defining a new S4 class... hence e don't manipulate class
.subcol_wAttr <- function(X, j, drop) {
  Xattr <- attributes(X)
  X <- X[,j=j,drop=drop]
  names_lostattrs <- setdiff(names(Xattr), names(attributes(X)))
  attributes(X)[names_lostattrs] <- Xattr[names_lostattrs] ## not mostattributes hich messes S4 objects ?!
  return(X)
}
