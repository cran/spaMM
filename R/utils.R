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
.pretty_summ_lambda <- function(lambda_est,processed) {
  cum_n_u_h <- processed$cum_n_u_h
  nrand <- length(cum_n_u_h)-1L
  ulam <- numeric(nrand)
  for (rd in seq_len(nrand)) ulam[rd] <- unique(lambda_est[(cum_n_u_h[rd]+1L):(cum_n_u_h[rd+1L])])
  return(ulam)
}

## quite ad hoc
.prettify_num <- function(x,nsmall=2L) {
  if (is.na(x) || is.null(x)) {
    return(x)
  } else if (abs(x)>1) {
    Ldigits <- floor(log10(abs(x)))+nsmall+1L
    format(round(x,Ldigits), digits=nsmall, nsmall=nsmall) ## character !
  } else signif(x,6) ## double !
}

.prettify_list <- function(li,nsmall=2L) {
  lapply(li,.prettify_num,nsmall=nsmall)
}

.prompt <- function() { # levite qq lourdeurs de browser(); 
  cat("\nPause: 'c' to continue, 'b' for browsing, 's' to stop") # other ahrs equivalent to c
  res <- readLines(n = 1)
  if (res=='s') stop("User-requested stop().")
  if (res=='b') browser("To debug, run options(error=recover); stop()" )
}

# must work on S4 objects but we avoid defining a new S4 class... hence e don't manipulate class
.subcol_wAttr <- function(X, j, drop) {
  Xattr <- attributes(X)
  X <- X[,j=j,drop=drop]
  names_lostattrs <- setdiff(names(Xattr), names(attributes(X)))
  attributes(X)[names_lostattrs] <- Xattr[names_lostattrs] ## not mostattributes hich messes S4 objects ?!
  return(X)
}

## F I X M E not yet used, to be tested next time I see an inline call from do.call(...)
print.spcall <- function (x, digits = NULL, quote = TRUE, na.print = NULL, print.gap = NULL, 
          right = FALSE, max = NULL, useSource = TRUE, ...) {
  ## print.default: 
  # args <- pairlist(digits = digits, quote = quote, na.print = na.print, 
  #                  print.gap = print.gap, right = right, max = max, useSource = useSource, 
  #                  ...)
  # missings <- c(missing(digits), missing(quote), missing(na.print), 
  #               missing(print.gap), missing(right), missing(max), missing(useSource))
  # .Internal(print.default(x, args, missings))
  ##
  abyss <- lapply(x, str) 
}

.all.subsets <- function(set) {
  if (is.null(set)) return(list(list())) # so that resu is always a list of lists, with at least one element.
  len <- length(set)
  TFsubsets <- expand.grid(rep(list(c(TRUE, FALSE)),len))
  resu <- apply(TFsubsets, 1L, function(TF) { set[TF] })
  if (! length(resu)) resu <- list(list()) # the case for input set=list()
  return(resu)
}

# See https://stackoverflow.com/questions/32806974/detecting-whether-shiny-runs-the-r-code
.check_frames <- function (which) {
  if (is.null(which)) return(0L)
  # Look for `which` call somewhere in the call stack.
  frames <-  sys.frames()
  calls <- lapply(sys.calls(), `[[`, 1)
  call_name <- function (call)
    if (is.function(call)) '<closure>' else deparse(call)
  call_names <- vapply(calls, call_name, character(1))
  target_call <- grep(which, call_names)
  return(length(target_call))
}

.mvrnorm <- function(n = 1, mu, Sigma, tol = 1e-6, empirical = FALSE, tcross_Sigma=NULL) {
  if (is.null(tcross_Sigma)) {
    return(mvrnorm(n, mu, Sigma, tol, empirical))
  } else { # simplification of MASS::mvrnorm code
    p <- length(mu)
    nr <- nrow(tcross_Sigma)
    nc <- ncol(tcross_Sigma)
    if (nr != p) stop("incompatible arguments (nrow(tcross_Sigma) != length(mu))")
    X <- matrix(rnorm(nc * n), n)
    if (empirical) { # not tested
      X <- scale(X, TRUE, FALSE)
      X <- X %*% svd(X, nu = 0)$v
      X <- scale(X, FALSE, TRUE)
    }
    X <- drop(mu) + tcrossprod(tcross_Sigma, X) #eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(tcross_Sigma))) 
      nm <- dn[[1L]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1) 
      drop(X)
    else as.matrix(t(X))
  }
}
