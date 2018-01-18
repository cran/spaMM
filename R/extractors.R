# getPar extract values from a list of lists, controlling that there is no redundancies between the lists => useful t merge lists 
# 'which' can be any way of indexing a list
.getPar <- function(parlist,name,which=NULL) {
  if ( ! is.null(which)) parlist <- parlist[[which]] 
  val <- parlist[[name]] 
  if (is.null(val)) { ## ie name not found a topmost level; scan sublists:
    vallist <- lapply(parlist, function(sublist) {
      if (is.list(sublist)) {sublist[[name]]} else {NULL}
    })
    ll <- sapply(vallist,length)
    ll <- which(ll>0)
    if (length(ll)>1) {
      stop(paste("Found several instances of element '",name,"' in nested list: use 'which' to resolve this.",sep=""))
    } else if (length(ll)==0L) {
      val <- NULL
    } else val <- vallist[[ll]]
  }
  val
}

# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"b") ## 2
# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"c") ## 4
# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"a") ## error
# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"a",which=1) ## 1
# .getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"d") ## NULL

# LMatrix assumed to be dense
.get_invL <- function(object, regul.threshold=1e-7) { ## computes inv(L) [not inv(Corr): see calc_invColdoldList]
  strucList <- object$strucList
  if ( ! is.null(strucList)) {
    if (is.null(object$envir$invL)) {
      cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      resu <- diag(cum_n_u_h[length(cum_n_u_h)])
      if (object$spaMM.version < "2.2.116") {
        ranefs <- attr(object$ZAlist,"ranefs") 
      } else ranefs <- attr(object$ZAlist,"exp_ranef_strings") 
      for (Lit in seq_len(length(strucList))) {
        lmatrix <- strucList[[Lit]]
        if ( ! is.null(lmatrix)) {
          affecteds <- which(ranefs %in% attr(lmatrix,"ranefs"))
          type <-  attr(lmatrix,"type")
          condnum <- kappa(lmatrix,norm="1")
          invlmatrix <- NULL
          if ( ! is.null(latentL_blob <- attr(lmatrix,"latentL_blob"))) { ## from .process_ranCoefs
            compactchol_Q <- latentL_blob$compactchol_Q 
            if (is.null(compactchol_Q)) {
              invlmatrix <- .makelong(solve(latentL_blob$u),longsize=ncol(lmatrix),as_matrix=TRUE)
            } else invlmatrix <- .makelong(t(compactchol_Q),longsize=ncol(lmatrix),as_matrix=TRUE) ## L=Q^{-T} => invL=Q^T
            # as_matrix necessary for resu[u.range, u.range] <- invlmatrix
          } else if (type == "from_Q_CHMfactor")  {
            invlmatrix <- t(as(attr(lmatrix,"Q_CHMfactor"),"sparseMatrix")) ## L=Q^{-T} => invL=Q^T
            invlmatrix <- as.matrix(invlmatrix) ## for [<-.matrix
          } else if (type == "cholL_LLt")  {
            if (condnum<1/regul.threshold) {
              invlmatrix <- try(forwardsolve(lmatrix,diag(ncol(lmatrix))),silent=TRUE)
              if (inherits(invlmatrix,"try-error")) invlmatrix <- NULL
            }
            if (is.null(invlmatrix)) Rmatrix <- t(lmatrix)
          } else { ## Rcpp's symSVD, or R's eigen() => LDL (also possible bad use of R's svd, not normally used)
            if (condnum<1/regul.threshold) {
              decomp <- attr(lmatrix,attr(lmatrix,"type")) ## of corr matrix !
              if ( all(abs(decomp$d) > regul.threshold) ) {
                invlmatrix <-  try(.ZWZt(decomp$u,sqrt(1/decomp$d)),silent=TRUE)
                if (inherits(invlmatrix,"try-error")) invlmatrix <- NULL
              }
            }
            if (is.null(invlmatrix)) Rmatrix <- qr.R(qr(t(lmatrix))) 
          }
          if (is.null(invlmatrix)){
            # chol2inv is quite robust in the sens of not stopping, even without any regularization.
            # Nevertheless (1) lmatrix %*% invlmatrix may be only roughly = I:
            #   if we don't regularize we expect departures from I due to numerical precision;
            #   if we regularize we expect departures from I even with exact arithmetic...
            #              (2) unregul. chol2inv result may still cause problems in later computations ?
            singular <- which(abs(diag(Rmatrix))<regul.threshold) 
            if (length(singular)) {
              if (spaMM.getOption("wRegularization")) warning("regularization required.")
              nc <- ncol(Rmatrix)
              diagPos <- seq.int(1L,nc^2,nc+1L)[singular]
              Rmatrix[diagPos] <- sign(Rmatrix[diagPos])* regul.threshold
            }
            invLLt <- chol2inv(Rmatrix) ## 
            invlmatrix <- crossprod(lmatrix, invLLt) ## regularized (or not) solve(lmatrix)
          }
          for (it in affecteds) {
            u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
            resu[u.range,u.range] <- invlmatrix   ## FR->FR it would be more eficcent to maintain a block structure here ?
            # but this allows an lmatrix to affect several rand effects
          }  
        }
      }
      object$envir$invL <- resu
    }
    return(object$envir$invL)
  } else return(NULL)
}




## fitted= X.beta + ZLv where we want to be able to write Lv as Cw = L.L'.w 
# => w = inv(L').v
.calc_invL_coeffs <- function(object,newcoeffs) { ##replaces dv/dloglan  by dw/dloglam
  strucList <- object$strucList
  if ( ! is.null(strucList)) {
    cum_n_u_h <- attr(object$lambda,"cum_n_u_h") ## FIXME: avoid using object$lambda
    if (object$spaMM.version < "2.2.116") {
      ranefs <- attr(object$ZAlist,"ranefs") 
    } else ranefs <- attr(object$ZAlist,"exp_ranef_strings") 
    for (Lit in seq_len(length(strucList))) {
      lmatrix <- strucList[[Lit]]
      if (! is.null(lmatrix)) { ## spatial or random-coef
        lmatrixranefs <- attr(lmatrix,"ranefs")
        if (object$spaMM.version>"1.11.56" && is.null(lmatrixranefs)) stop('attr(lmatrix,"ranefs") missing')
        affecteds <- which(ranefs %in% lmatrixranefs)
        for (it in affecteds) {
          u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
          ## dense calculation on not necess triangular lmatrix (!?). 
          ## Could surely be optimized, but by redfining lamtrix à la sXaug 
          newcoeffs[u.range] <- solve(t(lmatrix),newcoeffs[u.range])   ## transpose
        }  
      }
    }
  }
  return(newcoeffs)
}


.getHLfit <- function(fitobject) {
  if (inherits(fitobject,"HLfit")) {
    fitobject    
  } else if (inherits(fitobject,"HLCor")) {
    fitobject$hlfit    
  } else if (inherits(fitobject,"corrHLfit")) {
    fitobject$hlfit    
  }
} 

fitted.HLfit <- function(object,...) {
  object <- .getHLfit(object)
  object$fv
}

.get_BinomialDen <- function(object) {
  if (object$spaMM.version<"2.0.17") {
    return(object$weights)
  } else return(object$BinomialDen)
}

residuals.HLfit <- function(object, type = c("deviance", "pearson", "response"), ...) {
  object <- .getHLfit(object)
  type <- match.arg(type)
  BinomialDen <- .get_BinomialDen(object) 
  if (is.null(BinomialDen)) BinomialDen <- 1L
  y <- object$y / BinomialDen ## le y 
  mu <- object$fv # on 0-1 probability scale for binomial models
  wts <- object$prior.weights * BinomialDen ## the BinomialDen are in the $prior.weights of a glm object, but not of an HLfit one
  res <- switch(type, deviance = #if (object$df.residual > 0) 
                  { d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, 
                                                  wts), 0))
    ifelse(y > mu, d.res, -d.res)
  } #else rep.int(0, length(mu))
  , pearson = (y - mu) * sqrt(wts)/sqrt(object$family$variance(mu)), 
  response = y - mu)
  #if (!is.null(object$na.action))  res <- naresid(object$na.action, res)
  res
}

ranef.HLfit <- function(object,type="correlated",...) {
  object <- .getHLfit(object)
  lambda.object <- object$lambda.object
  print_namesTerms <- lambda.object$print_namesTerms
  repGroupNames <- unlist(lapply(seq_len(length(print_namesTerms)), function(it) {
    names(print_namesTerms[[it]]) <- rep(names(print_namesTerms)[it],length(print_namesTerms[[it]]))
  })) ## makes group identifiers unique (names of coeffs are unchanged)
  if (type=="correlated") {
    uv_h <- object$v_h 
  } else uv_h <- object$ranef #random effects \eqn{u}
  cum_n_u_h <- attr(object$ranef,"cum_n_u_h")
  if (object$spaMM.version < "2.2.116") {
    ranefs <- attr(object$ZAlist,"ranefs") 
  } else ranefs <- attr(object$ZAlist,"exp_ranef_strings") 
  colNames <- lapply(object$ZAlist,"colnames")
  # compute Lv from v:
  strucList <- object$strucList
  RESU <- vector("list", length(ranefs))
  for (it in seq_along(ranefs)) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    res <- uv_h[u.range] # vector
    if (length(print_namesTerms[[it]])>1L) { # random-coef term
      n_cols <- length(colNames[[it]])/length(print_namesTerms[[it]])
      if (type == "correlated") {
        res <- strucList[[it]] %*% res ## matrix
        res <- structure(res, dimnames=list(colNames[[it]][1:n_cols], print_namesTerms[[it]]))
      } else {
        res <- structure(matrix(res, ncol=n_cols,byrow=TRUE), dimnames= list(NULL, colNames[[it]][1:n_cols])) ## matrix
      }
    } else {
      if (type == "correlated" && ! is.null(strucList[[it]])) {
        res <- as.vector(strucList[[it]] %*% res) ## vector
      } 
      names(res) <- colNames[[it]]
    }
    RESU[[it]] <- res
  }
  names(RESU) <- ranefs
  RESU
}


fixef.HLfit <- function(object,...) {
  object <- .getHLfit(object)
  object$fixef    
}

logLik.HLfit <- function(object, which=NULL, ...) {
  object <- .getHLfit(object)
  if (is.null(which)) {
    mess <- .REMLmess(object)
    which <- switch(mess, 
                    "by stochastic EM."= "logLapp",
                    "by Laplace ML approximation (p_v)."= "p_v",
                    "by h-likelihood approximation."= "p_v",
                    "by ML."= "p_v",
                    "by Laplace REML approximation (p_bv)."= "p_bv",
                    "by REML."= "p_bv",
                    "by non-standard REML"= "p_bv",
                    stop(paste("No default '",which,"' value for '",mess,"' estimation method.",sep=""))
                    ) 
  }
  resu  <- object$APHLs[[which]]
  names(resu) <- which
  return(resu)
}

vcov.HLfit <- function(object,...) {
  object <- .getHLfit(object)
  beta_cov <- .get_beta_cov_any_version(object)
  class(beta_cov) <- c("vcov.HLfit",class(beta_cov))
  return(beta_cov)
}

.get_beta_cov_any_version <- function(object) {
  beta_cov <- object$beta_cov ## set by HLfit using get_from_MME(, which="beta_cov")
  if (is.null(beta_cov)) {
    return(.get_beta_cov(object))
  } else return(beta_cov)
}

# addition post 1.4.4
Corr <- function(object,...) { ## compare ?VarCorr
  trivial <- "No non-trivial correlation matrix for this random effect"
  strucList <- object$strucList
  locfn <- function(it) {
    resu <- object$cov.mats[[it]]
    if (is.null(resu)) resu <- .tcrossprod(strucList[[it]],NULL)
    if (is.null(resu)) resu <- trivial 
    return(resu) ## may still be NULL
  }
  if ( ! is.null(strucList)) {
    resu <- lapply(seq_len(length(strucList)), locfn) ## list for the different ranefs   
  } else {
    message(trivial)
    resu <- list(`1`=trivial) 
  }
  return(resu)
}

dev_resids <- function(object,...) {
  mu <- predict(object)
  BinomialDen <- .get_BinomialDen(object) 
  if (is.null(BinomialDen)) BinomialDen <- 1
  object$family$dev.resids(object$y/BinomialDen,mu,BinomialDen)
}

deviance.HLfit <- function(object,...) {
  dev_res2 <- dev_resids(object=object,...)
  return(sum(dev_res2))
}  

get_fixefVar <- function(...) {
  mc <- match.call(expand.dots = TRUE)
  mc$variances$fixefVar <- TRUE
  mc[[1L]] <- quote(spaMM::predict.HLfit)
  attr(eval(mc,parent.frame()),"fixefVar")
}

get_predVar <- function(...) {
  mc <- match.call(expand.dots = TRUE)
  mc$variances$predVar <- TRUE
  mc[[1L]] <- quote(spaMM::predict.HLfit)
  attr(eval(mc,parent.frame()),"predVar")
}

get_residVar <- function(...) {
  mc <- match.call(expand.dots = TRUE)
  mc$variances$residVar <- TRUE
  mc[[1L]] <- quote(spaMM::predict.HLfit)
  attr(eval(mc,parent.frame()),"residVar")
}

get_respVar <- function(...) {
  mc <- match.call(expand.dots = TRUE)
  mc$variances$respVar <- TRUE
  mc[[1L]] <- quote(spaMM::predict.HLfit)
  attr(eval(mc,parent.frame()),"respVar")
}

get_intervals <- function(...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- quote(spaMM::predict.HLfit)
  if (is.null(mc$intervals)) mc$intervals <- "respVar"
  attr(eval(mc,parent.frame()),"intervals") ## intervals by gaussian approx (possibly student'), not LR 
}


get_any_IC <- function(object,...,verbose=interactive()) {
  info_crits <- .get_info_crits(object)
  likelihoods <- numeric(0)
  if (!is.null(info_crits$mAIC))  likelihoods <- c(likelihoods,"       marginal AIC:"=info_crits$mAIC)
  if (!is.null(info_crits$cAIC))  likelihoods <- c(likelihoods,"    conditional AIC:"=info_crits$cAIC)
  if (!is.null(info_crits$dAIC))  likelihoods <- c(likelihoods,"     dispersion AIC:"=info_crits$dAIC)
  if (!is.null(info_crits$GoFdf)) likelihoods <- c(likelihoods,"       effective df:"=info_crits$GoFdf)
  if (verbose) {
    astable <- as.matrix(likelihoods)
    write.table(format(astable, justify="right"), col.names=FALSE, quote=FALSE) 
  }
  invisible(likelihoods)
}

AIC.HLfit <- function(object, ..., k,verbose=interactive()) {
  get_any_IC(object,...,verbose=verbose)
}

get_RLRTSim_args <- function(object,...) {
  if (object$family$family !="gaussian" 
      || (object$models[["eta"]]=="etaHGLM" && any(attr(object$rand.families,"lcrandfamfam")!="gaussian"))) {
    warning("'object' is not the fit of a LMM, while RLRTSim() methods were conceived for LMMs.")
  }
  X.pv <- as.matrix(object$X.pv)
  qrX.pv <- qr(X.pv)
  ZAL <- get_ZALMatrix(object)
  if(inherits(ZAL,"Matrix")) ZAL <- as.matrix(ZAL)
  sqrt.s <- diag(ncol(ZAL))
  return(list(X=X.pv, Z=ZAL, qrX = qrX.pv, sqrt.Sigma=sqrt.s))
}
# get_dispVar : dispVar in not a returned attribute

get_rankinfo <- function(object) return(attr(object$X.pv,"rankinfo")) 
