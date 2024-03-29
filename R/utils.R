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

# ~ suppressWarnings(try(, silent=TRUE)) ; try() doc recommends tryCatch() when silent=TRUE
.silent_W_E <- function(expr) {
  w.handler <- function(w){ # warning handler
    invokeRestart("muffleWarning")
  }
  withCallingHandlers(tryCatch(expr, error = function(e) e),
                      warning = w.handler)
}

.silent_M_E <- function(expr) {
  m.handler <- function(m){ 
    invokeRestart("muffleMessage")
  }
  withCallingHandlers(tryCatch(expr, error = function(e) e),
                      message = m.handler)
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

# must work on S4 objects but we avoid defining a new S4 class... hence we don't manipulate class
# This is called on 'A'matrices and on X.pv. In the latter case for an mv fit colnames MUST be present.
.subcol_wAttr <- function(X, j, drop) {
  Xattr <- attributes(X)
  X <- X[,j=j,drop=drop]
  names_lostattrs <- setdiff(names(Xattr), names(attributes(X)))
  if ( ! is.null(col_ranges <- Xattr$col_ranges)) { # mv model; tested but numInfo(zut0) in test-mv-extra which is NOT part of long tests.
    n_subm <- length(col_ranges)
    ncol_vec <- integer(n_subm)
    colmatches <- match(Xattr$dimnames[[2]], colnames(X))
    for (mv_it in seq_len(n_subm)) {
      colrange <- col_ranges[[mv_it]]  
      col_ranges[[mv_it]] <- na.omit(colmatches[colrange])  
    }
    Xattr$col_ranges <- col_ranges
  }
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

# if bool = TRUE, returns test whether in Rstudio session
# if bool = FALSE  (useful only for devel), returns path of project, or NULL, or a try-error
.inRstudio <- function(silent=FALSE, bool=TRUE) {
  tryres <- exists("RStudio.Version", envir = globalenv()) # from rstudioapi:::callLauncherFun 
  if (tryres && ! bool) { # we're in Rstudio
    trygetfn <- try(get("getActiveProject",envir = asNamespace("rstudioapi")), silent=silent) # fails if rstudioapi not installed
    if ( tryres <- ( ! inherits(trygetfn,"try-error"))) {
      tryres <- try(trygetfn(), silent=silent) # NULL if no active project, or the path if there is an active project.
    }
    tryres 
  } else tryres
}

projpath <- local({
  pp <- NULL
  function() {
    if (is.null(ppp <- .spaMM.data$options$projpath)) {
      if (is.null(pp)) {
        projpathinRstudio <- .inRstudio(silent=FALSE, bool=FALSE)
        if (inherits(projpathinRstudio,"try-error") || is.null(projpathinRstudio)) { # not an Rstudio session || no active project
          if (interactive()) {
            message('Need to give the project path, say "D:/home/francois/travail/stats/spaMMplus/spaMM":')
            pp <<- readline(prompt="Enter path: ")
          } else {
            message('Need to start in the projpath, say "D:/home/francois/travail/stats/spaMMplus/spaMM", so that getwd() finds it.')
            pp <<- getwd()                  # would reach here in an R CMD check session (but in principle projpath() is never called in this context, 
                                            # except perhaps wrapped in a if (file_exists())
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

..get_bare_fnname <- function(fun) { #Back compat. function => can ignore spaMM::fitmv
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
  if (is.null(fnname <- object$how$fnname)) fnname <- ..get_bare_fnname(fun=call.[[1L]])
  return(fnname)
}

.unlist <- function(x) unlist(x, recursive=FALSE, use.names = FALSE)

.relist_rep <- function(x, skeleton) relist(rep(x,length(unlist(skeleton,use.names = FALSE))),skeleton)

..trDiagonal <- function(n) .trDiagonal(n=n, unitri = FALSE)

# .rawsolve <- function(A) Matrix::solve(A, ..trDiagonal(n=ncol(A))) # wrapper useful only when ..trDiagonal() was memoised.

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
  for (it in seq_len(nl)) if (alis[[it]]@diag=="U") alis[[it]] <- Matrix::.diagU2N(alis[[it]], cl="CsparseMatrix") # using Matrix::.diagU2N() as public version of Matrix:::as_Csp2()
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
  # 2022/12/04: Removed a faulty line  here, copied from .bdiag_dtC, exposed as a spaMM bug with Matrix 1.5-3, as could be expected from the .diagU2N doc
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

.shermanMstep_de <- function(Ainv,u) { # ad-hoc code for perturbation of inverse of Qt.Q, using symmetry of Ainv, and the ad-hoc -2 factor
  # Ainv %*% kronecker(u,tu) %*% Ainv is the kronecker product of the folowing rowSums:
  rSAu <- rowSums(sweep(Ainv, MARGIN=2L, u, `*`)) # '.matrix_times_Dvec'
  num <- kronecker(t(rSAu),rSAu)  
  denom <- 1 -2*sum(u*(Ainv %*% u)[,1]) 
  Ainv +num/(denom/2)
}

.shermanM_de <- function(Qt, indic) { # ad-hoc code for perturbation of inverse of (Qt.Q=I) 
  # This thing is neither PD nor ND
  steps <- which(indic !=0)
  Ainv <- diag(nrow(Qt))
  for (it in steps) Ainv <- .shermanMstep_de(Ainv=Ainv,u=Qt[,it])
  Ainv
}

## The sparse version is pure Rcpp stuff, see .Rcpp_adhoc_shermanM_sp

### Utility for devel, produces the sXaug matrix for spcorr algo from an spprec sXaug matrix
### conceived for default values shown in formals, may need non-default values otherwise (sure for ZAL)
.spprec2spcorr <- function(spprec, obsInfo=TRUE, ZAL=spprec$AUGI0_ZX$ZAfix, zInfo=NULL) {
  AUGI0_ZX <- spprec$AUGI0_ZX
  H_w.resid <- .BLOB(spprec)$H_w.resid
  H_global_scale <- .calc_H_global_scale(H_w.resid)
  weight_X <- .calc_weight_X(H_w.resid, H_global_scale, obsInfo=obsInfo) ## sqrt(s^2 W.resid)  # -> .... sqrt(w.resid * H_global_scale)
  w.ranef <- attr(spprec,"w.ranef")
  ZAL_scaling <- 1/sqrt(w.ranef*H_global_scale) ## Q^{-1/2}/s
  ZAL <- .Matrix_times_Dvec(ZAL,ZAL_scaling) # should use another input ZAL more generally
  Zero_X <- rbind2(AUGI0_ZX$ZeroBlock, AUGI0_ZX$X.pv)
  I_ZAL <- rbind2(AUGI0_ZX$I, ZAL)
  Xscal <- cbind2(I_ZAL, Zero_X)
  attr(Xscal,"AUGI0_ZX") <- AUGI0_ZX # environment => cheap access to its 'envir$updateable' variable or anything else 
  spcorr <- do.call(def_sXaug_Matrix_CHM_H_scaled,
                    list(Xaug=Xscal, weight_X=weight_X, w.ranef=w.ranef, H_global_scale=H_global_scale))
  if (is.null(zInfo)) {
    spcorr
  } else {
    BLOBp <- .BLOB(spprec)
    BLOBc <- .BLOB(spcorr)
    typematch <- c(BLOBp$nonSPD, BLOBc$nonSPD, BLOBp$signs_in_WLS_mat, BLOBc$signs_in_WLS_mat)
    if ( ! all(typematch==c(F,T,T,F))) {
      # exclude case where  spprec G is SPD but spcorr H is nonSPD => WLS_mat differ
      wzAug <- c(zInfo$y2_sscaled/ZAL_scaling, (zInfo$z1_sscaled)*weight_X)
      Vscaled_beta <- get_from_MME(spcorr,szAug=wzAug) 
      seq_n_u_h <- seq_len(length(ZAL_scaling))
      v_h_beta <- list(v_h=Vscaled_beta[seq_n_u_h]*ZAL_scaling,
           beta=Vscaled_beta[-seq_n_u_h])
    } else v_h_beta <- NULL
    list(spcorr=spcorr, typematch=typematch, v_h_beta=v_h_beta)
  }
}

.safe_true <- function(cond) {( identical(cond,TRUE) ||       (( ! is.na(cond)) && is.numeric(cond) && cond>0) )   } # must exclude NA_real_ and allow 1 or 1L...

# .pMatrix_perm <- function(pMat) {
#   if (.hasSlot(pMat, "margin") && pMat@margin==2) {
#     sort.list(pMat@perm)
#   } else pMat@perm
# }
# 
# .tpMatrix_perm <- function(tpMat) {
#   if (.hasSlot(tpMat, "margin") && pMat@margin==2) {
#     sort.list(pMat@perm)
#   } else pMat@perm
# }

.dist_fn <- proxy::dist   #  see Zvariants/parallelDist.R for alternatives ...

# "someone suggested on the R-devel mailing list that maybe NCOL(NULL) would better give 0 instead of 1 as it currently does."
# This function is designed to reproduce the old NCOL() behaviour whatever the current NCOL() does.
.old_NCOL <- function(x) if (is.null(x)) {1L} else NCOL(x)
