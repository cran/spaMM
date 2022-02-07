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

.prompt <- function() { # evite qq lourdeurs de browser(); 
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

.call4print <- function (x, max.print=10L) {
  x <- as.list(x)
  for (it in seq_along(x)) {
    xx <- x[[it]]
    if (is.numeric(xx) & length(xx)> max.print) {
      x[[it]] <- capture.output(str(xx))
      names(x)[it] <- paste0("str(long numeric '",names(x[it]),"')")
    }
  }
  as.call(x)
}

print.bootci4print <- function(x, ...) {
  x$call <- .call4print(x$call)
  class(x) <- "bootci"
  print(x)
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

.inRstudio <- function(silent=FALSE, bool=TRUE) {
  trygetfn <- try(get("getActiveProject",envir = asNamespace("rstudioapi")), silent=silent) # fails if rstudioapi not installed
  if ( tryres <- ( ! inherits(trygetfn,"try-error"))) {
    tryres <- try(trygetfn(), silent=silent) # fails if rstudioapi not running ie not an Rstudio session
    # otherwise NULL if no active project, or the path if there is an active project.
  }
  if (bool) tryres <- ! inherits(tryres,"try-error") # to return a boolean : whether we are in an Rstudio session or not
  # otherwise retuen value is (try-error, or NULL, or path).
  tryres 
}

projpath <- local({
  pp <- NULL
  function() {
    if (is.null(ppp <- .spaMM.data$options$projpath)) {
      if (is.null(pp)) {
        projpathinRstudio <- .inRstudio(silent=FALSE, bool=FALSE)
        if (inherits(projpathinRstudio,"try-error") || is.null(projpathinRstudio)) { # not an Rstudio session || no active project
          if (interactive()) {
            message('Need to give the project path, say "C:/home/francois/travail/stats/spaMMplus/spaMM":')
            pp <<- readline(prompt="Enter path: ")
          } else {
            message('Need to start in the projpath, say "C:/home/francois/travail/stats/spaMMplus/spaMM", so that getwd() finds it.')
            pp <<- getwd()
          }
        } else pp <<- projpathinRstudio
      }
      return(pp)
    } else return(ppp)
  }
})

.bracket <- function(vec, min, max) {
  vec[vec < min] <- min
  vec[vec > max] <- max
  vec
}

.get_bare_fnname <- function(fun) {
  if (is.function(fun)) { # from do.call(spaMM::fitme, args = args) => it's a closure
    fnname <- names(which(sapply(list(fitme=spaMM::fitme,
                                      HLfit=spaMM::HLfit,
                                      HLCor=spaMM::HLCor,
                                      corrHLfit=spaMM::corrHLfit), identical, y=fun)))
    # the fact that we can rename functions shows how this can fail... 
    # solution is to have fnname in object$how$fnname, called before .get_bare_fnname()
    # but this is not backward-compatible, in which case...
    if ( ! length(fnname)) {
      funtext <- capture.output(fun)
      if (length(grep("fitme_body",funtext))) {
        fnname <- "fitme"
      } else if (length(grep("corrHLfit_body",funtext))) {
        fnname <- "fitme"
      } else if (length(grep("HLCor_body",funtext))) {
        fnname <- "HLCor"
      } else if (length(grep("HLfit_body",funtext))) {
        fnname <- "HLfit"
      } else stop("Unable to find the function that created the fit object. Retry after assigning its name to <fit object>$how$fnname .")
    }
    
  } else { # assuming it's a 'name' (direct call or do.call("fitme", args = args)) or fn got by by get(...)
    fnname <- sub("^.*?::","",deparse(fun)) ## remove any "spaMM::" in the name, from which paste() would return c("::","spaMM",<>)
  }
  return(fnname)
}

.get_bare_fnname.HLfit <- function(object, call.=getCall(object)) {
  if (is.null(fnname <- object$how$fnname)) fnname <- .get_bare_fnname(fun=call.[[1L]])
  return(fnname)
}

.unlist <- function(x) unlist(x, recursive=FALSE, use.names = FALSE)

..trDiagonal <- function(n) .trDiagonal(n=n, unitri = FALSE)

.rawsolve <- function(A) Matrix::solve(A, ..trDiagonal(n=ncol(A))) 

## specializations of Matrix::bdiag() and ::.bdiag():
.bdiag_dtC <- function(..., uplo="L") {  # if not all matrices are not consistently lower triangular (resp upper triangular) we should transpose some
  # if ((nA <- nargs()) == 0) return(new("dgCMatrix"))
  if ((nA <- nargs()) == 1L && !is.list(...)) 
    return(as(..., "CsparseMatrix"))
  alis <- if (nA == 1 && is.list(..1)) 
    ..1
  else list(...)
  if (length(alis) == 1L) 
    return(as(alis[[1]], "CsparseMatrix"))
  #
  nl <- length(alis)
  for (it in seq_len(nl)) alis[[it]] <- Matrix::.diagU2N(alis[[it]], cl="CsparseMatrix") # using Matrix::.diagU2N() as public version of Matrix:::as_Csp2()
  i_off <- c(0L, cumsum(vapply(alis, nrow, 1L)))
  i <- p <- x <- vector("list", nl)
  cumsum_p <- 0L
  for (it in seq_len(nl)) {
    i[[it]] <- alis[[it]]@i + i_off[it]
    p[[it]] <- cumsum_p+alis[[it]]@p[-1]
    x[[it]] <- alis[[it]]@x
    cumsum_p <- cumsum_p+length(i[[it]])
  }
  nc <- i_off[nl+1L]
  r <- new("dtCMatrix")
  r@x <- .unlist(x) 
  r@i <- .unlist(i)
  r@p <- c(0L,.unlist(p))
  r@Dim <- c(nc,nc)
  r@uplo <- uplo
  r
}

.bdiag_dsC <- function(...) {  
  # if ((nA <- nargs()) == 0) return(new("dgCMatrix"))
  if ((nA <- nargs()) == 1L && !is.list(...)) 
    return(as(..., "CsparseMatrix"))
  alis <- if (nA == 1 && is.list(..1)) 
    ..1
  else list(...)
  if (length(alis) == 1L) 
    return(as(alis[[1]], "CsparseMatrix"))
  #
  nl <- length(alis)
  Tlst <- vector("list", nl)
  for (it in seq_len(nl)) alis[[it]] <- Matrix::.diagU2N(alis[[it]], cl="CsparseMatrix") # using Matrix::.diagU2N() as public version of Matrix:::as_Csp2()
  #
  uplos <- vapply(alis, slot, ".", "uplo")
  tLU <- table(uplos)
  if (length(tLU) == 1) {
    useU <- uplos[1] == "U"
  } else {
    useU <- diff(tLU) >= 0L
    if (useU && (hasL <- tLU[1] > 0L)) {
      for (it in which(hasL)) Tlst[[hasL]] <- t(Tlst[[hasL]])
    } else if (!useU && (hasU <- tLU[2] > 0L)) {
      for (it in which(hasU)) Tlst[[hasU]] <- t(Tlst[[hasU]])
    }
  }
  i_off <- c(0L, cumsum(vapply(alis, nrow, 1L)))
  i <- p <- x <- vector("list", nl)
  cumsum_p <- 0L
  for (it in seq_len(nl)) {
    i[[it]] <- alis[[it]]@i + i_off[it]
    p[[it]] <- cumsum_p+alis[[it]]@p[-1]
    x[[it]] <- alis[[it]]@x
    cumsum_p <- cumsum_p+length(i[[it]])
  }
  nc <- i_off[nl+1L]
  r <- new("dsCMatrix")
  r@x <- .unlist(x) 
  r@i <- .unlist(i)
  r@p <- c(0L,.unlist(p))
  r@Dim <- c(nc,nc)
  r@uplo <-   if (useU) "U" else "L"
  r
}

# Not used. Most efficient of several syntaxes tried in devel/code_optimizations/structure.R
.structure_list  <- function(.Data, attr_list) {
  attr_list$names <- names(.Data)
  attributes(.Data) <- attr_list
  .Data
} 
