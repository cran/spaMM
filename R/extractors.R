# 'which' can be any way of indexing a list
getPar <- function(parlist,name,which=NULL) {
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

# getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"b") ## 2
# getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"c") ## 4
# getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"a") ## error
# getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"a",which=1) ## 1
# getPar(list("1"=list(a=1,b=2),"2"=list(a=3,c=4)),"d") ## NULL


calc_invL <- function(object,regul.threshold=1e-7) { ## computes inv(L) [not inv(Corr): see calc_invColdoldList]
  LMatrix <- attr(object$predictor,"LMatrix") 
  if ( ! is.null(LMatrix)) {
    cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
    resu <- diag(cum_n_u_h[length(cum_n_u_h)])
    ranefs <- attr(object$ZAlist,"ranefs") 
    if ( ! is.list(LMatrix)) LMatrix <- list(dummyid=LMatrix) ## so that next loop works
    for (Lit in seq_len(length(LMatrix))) {
      lmatrix <- LMatrix[[Lit]]
      affecteds <- which(ranefs %in% attr(lmatrix,"ranefs"))
      type <-  attr(lmatrix,"type")
      condnum <- kappa(lmatrix,norm="1")
      invlmatrix <- NULL
      if (type == "cholL_LLt")  {
        if (condnum<1/regul.threshold) {
          invlmatrix <- try(forwardsolve(lmatrix,diag(ncol(lmatrix))),silent=TRUE)
          if (inherits(invlmatrix,"try-error")) invlmatrix <- NULL
        }
        if (is.null(invlmatrix)) Rmatrix <- t(lmatrix)
      } else { ## Rcpp's symSVD, or R's eigen() => LDL (also possible bad use of R's svd, not normally used)
        if (condnum<1/regul.threshold) {
          decomp <- attr(lmatrix,attr(lmatrix,"type")) ## of corr matrix !
          if ( all(abs(decomp$d) > regul.threshold) ) {
            invlmatrix <-  try(ZWZt(decomp$u,sqrt(1/decomp$d)),silent=TRUE)
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
        if (length(singular)>0L) {
          if (spaMM.getOption("wRegularization")) warning("regularization required.")
          nc <- ncol(Rmatrix)
          diagPos <- seq.int(1L,nc^2,nc+1L)[singular]
          Rmatrix[diagPos] <- sign(Rmatrix[diagPos])* regul.threshold
        }
        invLLt <- chol2inv(Rmatrix) ## 
        invlmatrix <- t(lmatrix) %*% invLLt ## regularized (or not) solve(lmatrix)
      }
      for (it in affecteds) {
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        resu[u.range,u.range] <- invlmatrix   ## FR->FR it would be more eficcent to maintain a block structure here ?
        # but this allows an lmatrix to affect several rand effects
      }  
    }
    return(resu)
  } else return(NULL)
}

calc_invColdoldList <- function(object,regul.threshold=1e-7) { ## returns a list of inv(Corr) from the LMatrix
  LMatrix <- attr(object$predictor,"LMatrix") 
  if ( ! is.null(LMatrix)) {
    cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
    resu <- list()
    ranefs <- attr(object$ZAlist,"ranefs") 
    vec_n_u_h <- attr(attr(object$lambda,"cum_n_u_h"),"vec_n_u_h") ## !
    if ( ! is.list(LMatrix)) LMatrix <- list(dummyid=LMatrix) ## so that next loop works
    for (it in seq(ranefs)) resu[[it]] <- Diagonal(vec_n_u_h[it]) 
    for (Lit in seq_len(length(LMatrix))) {
      lmatrix <- LMatrix[[Lit]]
      affecteds <- which(ranefs %in% attr(lmatrix,"ranefs"))
      ## end of designL.from.corr implies either type is cholL_LLt or we have decomp $u and $d 
      type <-  attr(lmatrix,"type")
      invCoo <- NULL
      if (type == "cholL_LLt")  {
        Rmatrix <- t(lmatrix)
      } else { ## Rcpp's symSVD, or R's eigen() => LDL (also possible bad use of R's svd, not normally used)
        condnum <- kappa(lmatrix,norm="1")
        if (condnum<1/regul.threshold) {
          decomp <- attr(lmatrix,attr(lmatrix,"type")) ## of corr matrix !
          invCoo <-  try(ZWZt(decomp$u,1/decomp$d),silent=TRUE)
          if (inherits(invCoo,"try-error")) invCoo <- NULL
        }
        if (is.null(invCoo)) Rmatrix <- qr.R(qr(t(lmatrix))) 
      }
      if (is.null(invCoo)){ ## see comments on chol2inv in calc_invL()
        singular <- which(abs(diag(Rmatrix))<regul.threshold) 
        if (length(singular)>0L) {
          if (spaMM.getOption("wRegularization")) warning("regularization required.")
          nc <- ncol(Rmatrix)
          diagPos <- seq.int(1L,nc^2,nc+1L)[singular]
          Rmatrix[diagPos] <- sign(Rmatrix[diagPos])* regul.threshold
        }
        invCoo <- chol2inv(Rmatrix) ## 
      }
      for (aff in affecteds) resu[[aff]] <- invCoo
    }
    return(resu)
  } else return(NULL)
}


## fitted= X.beta + ZLv where we want to be able to write Lv as Cw = L.L'.w 
# => w = inv(L').v
calc_invL_coeffs <- function(object,newcoeffs) {
  LMatrix <- attr(object$predictor,"LMatrix") ## note transpose below
  if ( ! is.null(LMatrix)) {
    cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
    ranefs <- attr(object$ZAlist,"ranefs") 
    if ( ! is.list(LMatrix)) LMatrix <- list(dummyid=LMatrix)
    for (Lit in seq_len(length(LMatrix))) {
      lmatrix <- LMatrix[[Lit]]
      affecteds <- which(ranefs %in% attr(lmatrix,"ranefs"))
      for (it in affecteds) {
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        newcoeffs[u.range] <- solve(t(lmatrix),newcoeffs[u.range])   ## transpose
      }  
    }
  }
  return(newcoeffs)
}


getHLfit <- function(fitobject) {
  if (inherits(fitobject,"HLfit")) {
    fitobject    
  } else if (inherits(fitobject,"HLCor")) {
    fitobject$hlfit    
  } else if (inherits(fitobject,"corrHLfit")) {
    fitobject$hlfit    
  }
} 

fitted.HLfit <- function(object,...) {
  object <- getHLfit(object)
  object$fv
}

#ranef.HLfit <- function(object,ranef.class=TRUE,...) {
ranef.HLfit <- function(object,...) {
  object <- getHLfit(object)
  lambda.object <- object$lambda.object
  namesTerms <- lambda.object$namesTerms
  repGroupNames <- unlist(lapply(seq_len(length(namesTerms)),function(it) {
    names(namesTerms[[it]]) <- rep(names(namesTerms)[it],length(namesTerms[[it]]))
  })) ## makes group identifiers unique (names of coeffs are unchanged)
  coefficients <- unlist(lambda.object$namesTerms) ## FR->FR store this one for all in lambda.object ?
  ## cf Group and Term columns in output generated by summary.HL()
  res <- object$ranef #random effects \eqn{u}
  attr(res,"nams") <- paste(repGroupNames,coefficients) ## one name for each lambda coefficient. Cf bug detector in print.ranef
  class(res) <- c("ranef",class(res)) ## for print.ranef
  res
}

print.ranef <- function(x,...) {
  cum_n_u_h <- attr(x,"cum_n_u_h")
  nams <- attr(x,"nams")
  if (length(cum_n_u_h) != length(nams)+1L) { ## bug detector (hopefully no bug to be detected)
    message("length(cum_n_u_h) != length(nams)+1L in print.ranef: minor bug, but don't trust the output of print.ranef")
  }
  lapply(seq_len(length(nams)), function(it) {
    #cat(paste(nams[it], " (", cum_n_u_h[it + 1]-cum_n_u_h[it], " levels)\n",sep=""))    
    cat(paste(nams[it], " :\n",sep=""))    
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    print(x[u.range])
  })
  invisible(x)
}

fixef.HLfit <- function(object,...) {
  object <- getHLfit(object)
  object$fixef    
}

logLik.HLfit <- function(object, which=NULL, ...) {
  object <- getHLfit(object)
  if (is.null(which)) {
    mess <- REMLmess(object)
    which <- switch(mess, 
                    "by stochastic EM."= "logLapp",
                    "by ML approximation (p_v)."= "p_v",
                    "by h-likelihood approximation."= "p_v",
                    "by ML."= "p_v",
                    "by REML approximation (p_bv)."= "p_bv",
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
  object <- getHLfit(object)
  beta_cov <- get_beta_cov_any_version(object)
  class(beta_cov) <- c("vcov.HLfit",class(beta_cov))
  return(beta_cov)
}

get_beta_cov_any_version <- function(object) {
  if (object$spaMM.version>"1.9.22") {
    return(object$get_beta_cov(object))
  } else {
    return(object$beta_cov)
  }
}

# addition post 1.4.4
Corr <- function(object,...) { ## compare ?VarCorr
  if ( ! is.null(LMatrix <- attr(object$predictor,"LMatrix"))) {
    Corr <- tcrossprodCpp(LMatrix)    
  } else {
    message("No 'non-trivial' correlation matrix of random effects")
    Corr <- NA
  }
  return(Corr)
}

dev_resids <- function(object,...) {
  mu <- predict(object)
  weights <- object$weights
  if (is.null(weights)) weights <- 1
  object$family$dev.resids(object$y/weights,mu,weights)
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
  attr(eval(mc,parent.frame()),"intervals")
}


get_any_IC <- function(object,...,verbose=interactive()) {
  info_crits <- object$get_info_crits(object)
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
  ZAL <- as.matrix(object$get_ZALMatrix(object))
  sqrt.s <- diag(ncol(ZAL))
  return(list(X=X.pv, Z=ZAL, qrX = qrX.pv, sqrt.Sigma=sqrt.s))
}
# get_dispVar : dispVar in not a returned attribute