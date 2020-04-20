# LMatrix assumed to be dense except in the trivial case of identity matrix
.get_invL_HLfit <- function(object, regul.threshold=1e-7) { ## computes inv(L) [not inv(Corr): see calc_invColdoldList]
  strucList <- object$strucList
  if (is.null(strucList)) {
    return(NULL) ## no ranefs
  } else if (is.null(unlist(strucList))) { ## ranefs with trivial strucList
    if (is.null(object$envir$invL)) {
      cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      object$envir$invL <- Diagonal(n=cum_n_u_h[length(cum_n_u_h)]) 
    }
    return(object$envir$invL)
  } else {
    if (is.null(object$envir$invL)) {
      cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      resu <- diag(cum_n_u_h[length(cum_n_u_h)]) # has to be dense since we will assign blocks
      if (object$spaMM.version < "2.2.116") {
        ranefs <- attr(object$ZAlist,"ranefs") 
      } else ranefs <- attr(object$ZAlist,"exp_ranef_strings") 
      for (Lit in seq_len(length(strucList))) {
        lmatrix <- strucList[[Lit]]
        if ( ! is.null(lmatrix)) {
          type <-  attr(lmatrix,"type")
          invlmatrix <- NULL
          if ( ! is.null(latentL_blob <- attr(lmatrix,"latentL_blob"))) { ## from .process_ranCoefs
            compactchol_Q <- latentL_blob$compactchol_Q 
            if (is.null(compactchol_Q)) {
              invlmatrix <- .makelong(solve(latentL_blob$design_u),longsize=ncol(lmatrix),as_matrix=TRUE)
            } else invlmatrix <- .makelong(t(compactchol_Q),longsize=ncol(lmatrix),as_matrix=TRUE) ## L=Q^{-T} => invL=Q^T
            # as_matrix necessary for resu[u.range, u.range] <- invlmatrix
          } else if (inherits(lmatrix,"dCHMsimpl")) { # before any test on type...
            invlmatrix <- t(as(lmatrix, "sparseMatrix"))
          } else if (type == "from_AR1_specific_code")  {
            invlmatrix <- solve(lmatrix) # cost of solve sparse triangular matrix
            invlmatrix <- as.matrix(invlmatrix) ## for [<-.matrix
          } else if (type == "from_Q_CHMfactor")  {
            invlmatrix <- t(as(attr(lmatrix,"Q_CHMfactor"),"sparseMatrix")) ## L=Q^{-T} => invL=Q^T ## correct but requires the attribute => numerical issues in computing Q_CHMfactor
            invlmatrix <- as.matrix(invlmatrix) ## for [<-.matrix
          } else if (type == "cholL_LLt")  {
            condnum <- kappa(lmatrix,norm="1")
            if (condnum<1/regul.threshold) {
              invlmatrix <- try(forwardsolve(lmatrix,diag(ncol(lmatrix))),silent=TRUE)
              if (inherits(invlmatrix,"try-error")) invlmatrix <- NULL
            }
            if (is.null(invlmatrix)) Rmatrix <- t(lmatrix)
          } else { ## Rcpp's symSVD, or R's eigen() => LDL (also possible bad use of R's svd, not normally used)
            condnum <- kappa(lmatrix,norm="1")
            if (condnum<1/regul.threshold) {
              decomp <- attr(lmatrix,attr(lmatrix,"type")) ## of corr matrix !
              if ( all(abs(decomp$d) > regul.threshold) ) {
                invlmatrix <-  try(.ZWZt(decomp$u,sqrt(1/decomp$d)),silent=TRUE) ## try() still allowing for no (0) regul.threshold; not useful ?
                if (inherits(invlmatrix,"try-error")) invlmatrix <- NULL
              }
            }
            if (is.null(invlmatrix)) Rmatrix <- .lmwithQR(t(lmatrix),yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled # no pivoting compared to qr.R(qr(t(lmatrix))) 
          }
          if (is.null(invlmatrix)){
            # chol2inv is quite robust in the sense of not stopping, even without any regularization.
            # Nevertheless (1) lmatrix %*% invlmatrix may be only roughly = I:
            #   if we don't regularize we expect departures from I due to numerical precision;
            #   if we regularize we expect departures from I even with exact arithmetic...
            #
            # But regul. chol2inv result still causes problems in later computations!
            #
            # singular <- which(abs(diag(Rmatrix))<regul.threshold) 
            # if (length(singular)) {
            #   if (spaMM.getOption("wRegularization")) warning("regularization required.")
            #   nc <- ncol(Rmatrix)
            #   diagPos <- seq.int(1L,nc^2,nc+1L)[singular]
            #   Rmatrix[diagPos] <- sign(Rmatrix[diagPos])* regul.threshold
            # }
            # invLLt <- chol2inv(Rmatrix) ## 
            #
            invLLt <- try(chol2inv(Rmatrix),silent=TRUE)
            if (inherits(invLLt,"try-error") || max(abs(range(invLLt)))> 1e12) {
              invLLt <- ginv(crossprod(Rmatrix))
            }
            invlmatrix <- .crossprod(lmatrix, invLLt) ## regularized (or not) solve(lmatrix)
          }
          u.range <- (cum_n_u_h[Lit]+1L):(cum_n_u_h[Lit+1L])
          resu[u.range,u.range] <- invlmatrix   
        }
      }
      resu <- as(resu,"sparseMatrix") # previously broke test twolambda vs onelambda by effect on nearly singular matrix:
      # ...by changing the twolambda result.
      # # Without it I haD two messages : 
      # 1: In .calc_logdisp_cov(object, dvdloglamMat = dvdloglamMat, dvdlogphiMat = dvdlogphiMat,  :
      # Numerical precision issue in computation of the information matrix for dispersion parameters:
      #   the prediction variance may be inaccurate.
      # 2: In .force_solve(logdispInfo) : The matrix looks exactly singular.
      # # While with it I haD only the first message, and a result less consistent with onelambda
      # Alternatively to the present conversion, 
      # I could reproduce these two symptoms (only the first message, and a result less consistent) 
      # by calling .crossprodCpp in .calc_invV_factors() -> .crossprod(ZAfix, wrZ)
      # (as tested by hacking the return value of  the single call to .crossprod() in this get_predVar() test)
      # and the only impact of the .crossprod() is here to change the numerical precision of the $r_x_n element in the 
      # return value of .calc_invV_factors() by effects of order 1e-13 (this element being dgeMatrix whether .crossprodCpp was called or not).
      object$envir$invL <- resu
    }
    return(object$envir$invL) # May be Diagonal() => calling code must handle that case.
  }
}

.useless_get_sp_invL_HLfit <- function(object, regul.threshold=1e-7) { ## computes inv(L) [not inv(Corr): see calc_invColdoldList]
  strucList <- object$strucList
  if (is.null(strucList)) {
    return(NULL) ## no ranefs
  } else if (is.null(unlist(strucList))) { ## ranefs with trivial strucList
    if (is.null(object$envir$invL)) {
      cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      object$envir$invL <- Diagonal(n=cum_n_u_h[length(cum_n_u_h)])
    }
    return(object$envir$invL)
  } else {
    if (is.null(object$envir$invL)) {
      cum_n_u_h <- attr(object$lambda,"cum_n_u_h")
      seq_n_u_h <- diff(cum_n_u_h)
      resu <- vector("list", length(strucList))
      if (object$spaMM.version < "2.2.116") {
        ranefs <- attr(object$ZAlist,"ranefs") 
      } else ranefs <- attr(object$ZAlist,"exp_ranef_strings") 
      for (Lit in seq_len(length(strucList))) {
        lmatrix <- strucList[[Lit]]
        if ( is.null(lmatrix)) {
          resu[[Lit]] <- Diagonal(n=seq_n_u_h[Lit])
        } else if (inherits(lmatrix,"dCHMsimpl")) { # before any test on type...
          resu[[Lit]] <- t(as(lmatrix, "sparseMatrix"))
        } else {
          type <-  attr(lmatrix,"type")
          if ( ! is.null(latentL_blob <- attr(lmatrix,"latentL_blob"))) { ## from .process_ranCoefs
            compactchol_Q <- latentL_blob$compactchol_Q 
            if (is.null(compactchol_Q)) {
              resu[[Lit]] <- .makelong(solve(latentL_blob$design_u),longsize=ncol(lmatrix),as_matrix=FALSE)
            } else resu[[Lit]] <- .makelong(t(compactchol_Q),longsize=ncol(lmatrix),as_matrix=FALSE) ## L=Q^{-T} => invL=Q^T
            # as_matrix necessary for resu[u.range, u.range] <- invlmatrix
          } else if (type == "from_AR1_specific_code")  {
            resu[[Lit]] <- solve(lmatrix) # cost of solve sparse triangular matrix
          } else if (type == "from_Q_CHMfactor")  {
            resu[[Lit]] <- t(as(attr(lmatrix,"Q_CHMfactor"),"sparseMatrix")) ## L=Q^{-T} => invL=Q^T ## correct but requires the attribute => numerical issues in computing Q_CHMfactor
          } else if (type == "cholL_LLt")  {
            condnum <- kappa(lmatrix,norm="1")
            if (condnum<1/regul.threshold) {
              resu[[Lit]] <- try(forwardsolve(lmatrix,diag(ncol(lmatrix))),silent=TRUE)
              if (inherits(resu[[Lit]],"try-error")) invlmatrix <- NULL
            }
            if (is.null(resu[[Lit]])) {
              Rmatrix <- t(lmatrix)
            } else resu[[Lit]] <- as(resu[[Lit]],"sparseMatrix")
          } else { ## Rcpp's symSVD, or R's eigen() => LDL (also possible bad use of R's svd, not normally used)
            condnum <- kappa(lmatrix,norm="1")
            if (condnum<1/regul.threshold) {
              decomp <- attr(lmatrix,attr(lmatrix,"type")) ## of corr matrix !
              if ( all(abs(decomp$d) > regul.threshold) ) {
                resu[[Lit]] <-  try(.ZWZt(decomp$u,sqrt(1/decomp$d)),silent=TRUE) ## try() still allowing for no (0) regul.threshold; not useful ?
                if (inherits(resu[[Lit]],"try-error")) resu[[Lit]] <- NULL
              }
            }
            if (is.null(resu[[Lit]])) {
              Rmatrix <- .lmwithQR(t(lmatrix),yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled # no pivoting compared to qr.R(qr(t(lmatrix))) 
            } else resu[[Lit]] <- as(resu[[Lit]],"sparseMatrix")
          }
          if (is.null(resu[[Lit]])){
            # chol2inv is quite robust in the sense of not stopping, even without any regularization.
            # Nevertheless (1) lmatrix %*% invlmatrix may be only roughly = I:
            #   if we don't regularize we expect departures from I due to numerical precision;
            #   if we regularize we expect departures from I even with exact arithmetic...
            #
            # But regul. chol2inv result still causes problems in later computations!
            #
            # singular <- which(abs(diag(Rmatrix))<regul.threshold) 
            # if (length(singular)) {
            #   if (spaMM.getOption("wRegularization")) warning("regularization required.")
            #   nc <- ncol(Rmatrix)
            #   diagPos <- seq.int(1L,nc^2,nc+1L)[singular]
            #   Rmatrix[diagPos] <- sign(Rmatrix[diagPos])* regul.threshold
            # }
            # invLLt <- chol2inv(Rmatrix) ## 
            #
            invLLt <- try(chol2inv(Rmatrix),silent=TRUE)
            if (inherits(invLLt,"try-error") || max(abs(range(invLLt)))> 1e12) {
              invLLt <- ginv(crossprod(Rmatrix))
            }
            resu[[Lit]] <- .crossprod(lmatrix, invLLt) ## regularized (or not) solve(lmatrix)
            resu[[Lit]] <- as(resu[[Lit]],"sparseMatrix")
          }
        }
        resu <- do.call(Matrix::bdiag,resu)
      }
      object$envir$invL <- resu
    }
    return(object$envir$invL)
  }
}




## fitted= X.beta + ZLv where we want to be able to write Lv as Cw = L.L'.w 
# => w = inv(L').v
.calc_invL_coeffs <- function(object,newcoeffs) { 
  strucList <- object$strucList
  if ( ! is.null(strucList)) {
    cum_n_u_h <- attr(object$lambda,"cum_n_u_h") ## FIXME: avoid using object$lambda
    for (Lit in seq_len(length(strucList))) {
      lmatrix <- strucList[[Lit]]
      if (inherits(lmatrix,"dCHMsimpl")) {
        u.range <- (cum_n_u_h[Lit]+1L):(cum_n_u_h[Lit+1L])
        newcoeffs[u.range] <- as(lmatrix,"sparseMatrix") %*% newcoeffs[u.range]
      } else if (! is.null(lmatrix)) { ## spatial or random-coef
        u.range <- (cum_n_u_h[Lit]+1L):(cum_n_u_h[Lit+1L])
        ## dense calculation on not necess triangular lmatrix (!?). 
        ## solve( _t_(lmatrix)) may not allow the efficient use of solveWrap. 
        ## But this is a one-time calculation whose results are saved. No optimization attempted.
        newcoeffs[u.range] <- solve(t(lmatrix),newcoeffs[u.range])   ## newcoeffs must be a _vector_
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
  #colNames <- lapply(object$ZAlist,"colnames") without a slow lapply():
  colNames <- vector("list", length(object$ZAlist))
  names(colNames) <- names(object$ZAlist)
  for (it in seq_along(object$ZAlist)) {
    verif <- colnames(object$ZAlist[[it]])
    if ( ! is.null(verif)) colNames[[it]] <- verif # do not assign <- NULL to list member
    # the case NULL: I could stop(paste0("Misformed object: colnames(object$ZAlist[[",it,"]]) is NULL.")) 
    # but this occurrence may depend on user input through AMatrices...
  }
  # compute Lv from v:
  strucList <- object$strucList
  RESU <- vector("list", length(ranefs))
  for (it in seq_along(ranefs)) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    res <- uv_h[u.range] # vector
    if ((nr <- length(print_namesTerms[[it]]))>1L) { # random-coef term
      n_cols <- length(colNames[[it]])/length(print_namesTerms[[it]])
      if (type == "correlated") {
        res <- strucList[[it]] %*% res ## matrix
        res <- structure(matrix(res, ncol=n_cols,byrow=TRUE), dimnames=list(print_namesTerms[[it]],colNames[[it]][1:n_cols]))
        res <- t(res) # despite the t() it makes it ~ to the vector in the alternative case (both operate as n_u_h-row matrices) 
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
  class(RESU) <- c("ranef", "list")
  RESU ## TODO: ~lme4:::ranef.merMod(mod, condVar = TRUE) & ajouter des arguments "variances" et "intervals" à ta fonction ranef() (Alex 7/5/2018)
}

print.ranef <- function(x, max.print=40L, ...) {
  oldopt <- options(max.print=max.print)
  print.default(x)
  options(oldopt)
}

fixef.HLfit <- function(object,...) {
  object <- .getHLfit(object)
  object$fixef    
}

.get_phiW <- function(object, newdata, 
                      mu, # mu is a list (or data frame) whose length() is the # of replicates needed
                      prior.weights=object$prior.weights, phi_type, needed) { # phi/prior.weights
  phi_model <-  object$models[["phi"]]
  is_phiW_fix_btwn_sims <- FALSE
  if (phi_model == "") { ## the count families, or user-given phi
    if (length(object$phi)>1L) {
      if (is.null(newdata)) {
        message(paste0("simulate.HLfit() called on an original fit where phi was given but not constant.\n",
                       "This phi will be used, but is that relevant?"))
      } else stop("I do not know what to simulate when 'newdata' is not NULL and the original fit's phi was given but not constant.")
    }  
    newphiMat <- matrix(object$phi,ncol=length(mu),nrow=length(mu[[1L]]))
    if (identical(attr(prior.weights,"unique"),TRUE)) {
      phiW <- newphiMat/prior.weights[1L]
      is_phiW_fix_btwn_sims <- TRUE
    } else phiW <- .Dvec_times_matrix(1/prior.weights,newphiMat)  ## warnings or errors if something suspect
  } else if (phi_type=="predict") {
    newphiVec <- switch(phi_model,
                        "phiGLM" = predict(object$phi.object$glm_phi, 
                                           newdata=newdata, type="response"), ## vector
                        "phiHGLM" = predict(object$phi_model, newdata=newdata, type=phi_type)[ ,1L],
                        "phiScal" = rep(object$phi,length(mu[[1L]])),
                        stop('Unhandled object$models[["phi"]]')
    ) ## VECTOR in all cases, becomes matrix later
    if (identical(attr(prior.weights,"unique"),TRUE)) {
      phiW <- newphiVec/prior.weights[1L]
      if (phi_model %in% c("phiScal","phiGLM")) is_phiW_fix_btwn_sims <- TRUE
    } else phiW <- newphiVec/prior.weights  ## warnings or errors if something suspect
    phiW <- matrix(phiW,nrow=length(phiW), ncol=length(mu))  # vector -> matrix
  } else { # any other phi_type 
    newphiMat <- switch(phi_model,
                        "phiGLM" = as.matrix(simulate(object$phi.object$glm_phi, 
                                                      newdata=newdata, nsim=needed)), ## data frame -> matrix
                        "phiHGLM" = simulate(object$phi_model, newdata=newdata, type=phi_type, nsim=needed),
                        "phiScal" = matrix(object$phi,ncol=length(mu),nrow=length(mu[[1L]])),
                        stop('Unhandled object$models[["phi"]]')
    ) ## already MATRIX in all cases
    if (identical(attr(prior.weights,"unique"),TRUE)) {
      phiW <- newphiMat/prior.weights[1L]
      if (phi_model=="phiScal") is_phiW_fix_btwn_sims <- TRUE
    } else phiW <- .Dvec_times_matrix(1/prior.weights,newphiMat)  ## warnings or errors if something suspect
    # F I X M E add diagnostics ?
  } # phiW is always a matrix
  attr(phiW,"is_phiW_fix_btwn_sims") <- is_phiW_fix_btwn_sims
  return(phiW) ## always MATRIX
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
                    stop(paste0("No default '",which,"' value for '",mess,"' estimation method."))
                    ) 
  }
  if (which=="logL_Lap") {
    if (all(unlist(object$family[c("family","link")]==c("Gamma","log")))) {
      ZAL <- .compute_ZAL(XMatrix=object$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg(object)) 
      w.obs <- structure(object$w.resid * (object$y/object$fv)[,1],unique=FALSE)
      d2hdv2 <- .calcD2hDv2(ZAL,w.obs,object$w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
      hlik <- object$APHLs$hlik
      resu <- hlik -determinant(d2hdv2)$modulus[1L]/2 + ncol(d2hdv2)* log(2*pi)/2
    } else if (.is_link_canonical(object$family)) {
      message("logL_Lap = p_v because the response-family link is canonical")
      resu  <- object$APHLs[["p_v"]]
    } else stop("logL_Lap computation not yet implemented for this (family,link) combination.")
  } else if (which=="cliks") { ## undoc
    if (object$family$family=="binomial") {
      muFREQS <- predict(object)
      mu <- muFREQS * object$BinomialDen
    }
    clik_fn <- .get_clik_fn(object$family)
    phi_est <- .get_phiW(object=object, newdata=NULL, mu=as.data.frame(mu), phi_type="predict", needed=1L,...)[,1L] # ... allows to overcome the default prior.weights
    reweiclik <- .calc_clik(mu=mu,phi_est=phi_est,object, clik_fn=clik_fn, summand=TRUE, ...) # ... allows to overcome the default prior.weights
    colnames(reweiclik) <- "clik"
    return(reweiclik) ## currently 1-col matrix; avoids names(resu) <- which !
  } else resu  <- object$APHLs[[which]]
  names(resu) <- which
  return(resu)
}

vcov.HLfit <- function(object, ...) {
  object <- .getHLfit(object)
  beta_cov <- .get_beta_cov_any_version(object)
  class(beta_cov) <- c("vcov.HLfit",class(beta_cov))
  return(beta_cov)
}

.get_beta_cov_any_version <- function(object) {
  if (object$spaMM.version > "2.5.20") {
    beta_cov <- object$envir$beta_cov_info$beta_cov
  } else {
    beta_cov <- object$beta_cov
    attr(beta_cov,"beta_v_cov") <- NULL 
  } 
  if (is.null(beta_cov)) beta_cov <- .get_beta_cov_info(object)$beta_cov ## notably for GLM or if was stripped from envir by stripHLfit()
  return(beta_cov)
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

# for a lme4::VarCorr() equivalent; generic is nlme::VarCorr 
VarCorr.HLfit <- function(x, sigma=1, message.=TRUE, ...) {
  if (message.) message("The design of VarCorr.HLfit() is not yet fully specified.") 
  lambda.object <- x$lambda.object
  namesTerms <- lambda.object$print_namesTerms ## list of vectors of variable length
  linklam_coeff_list <- lambda.object$coefficients_lambdaS ## used beyond the next line
  lamtable <- .lambda_table_fn(namesTerms, x, lambda.object,linklam_coeff_list)
  loctable <- data.frame(Group=lamtable[,"Group"],Term=lamtable[,"Term"],Variance=lamtable[,"Var."],"Std.Dev."=sqrt(lamtable[,"Var."]))
  if ("Corr." %in% colnames(lamtable)) loctable <- cbind(loctable,Corr=lamtable[,"Corr."])
  rownames(loctable) <- lamtable[,"Term"]
  if (x$family$family %in% c("gaussian","Gamma")) {
    phi.object <- x$phi.object
    if ( ! is.null(phi_outer <- phi.object$phi_outer)) {
      phi_line <- data.frame(Group="Residual",Term="(Intercept)",Variance=phi_outer, "Std.Dev."=sqrt(phi_outer))
      if ("Corr." %in% colnames(lamtable)) phi_line <- cbind(phi_line,Corr=NA)
      loctable <- rbind(loctable,phi_line)
    } else {
      if (x$models[["phi"]]=="phiHGLM") {
        #cat("Residual dispersion model includes random effects:\n  use summary(<fit object>$resid_fit) to display results.\n")       
      } else if ((loc_p_phi <- length(phi.object$fixef))==1L) {
        namesX_disp <- names(phi.object$fixef)
        dispOffset <- attr(x$resid.predictor,"off")
        if (!is.null(dispOffset)) dispOffset <- unique(dispOffset)
        if (length(namesX_disp)==1 && namesX_disp[1]=="(Intercept)" && length(dispOffset)<2L) {
          # constant phi: we can display it
          phi_est <- (phi.object$fixef)
          if (length(dispOffset)==1L) phi_est <- phi_est+dispOffset
          resid.family <- eval(x$resid.family)
          phi_est <- resid.family$linkinv(phi_est)
        } ## else phi not constant; We don't try to display it
        phi_line <- data.frame(Group="Residual",Term="(Intercept)",Variance=phi_est, "Std.Dev."=sqrt(phi_est))
        if ("Corr." %in% colnames(lamtable)) phi_line <- cbind(phi_line,Corr=NA)
        loctable <- rbind(loctable,phi_line)
      }                                                 
    }
  }
  rownames(loctable) <- NULL
  return(loctable)
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

get_predVar <- function(..., variances=list(), which="predVar") {
  mc <- match.call(expand.dots = TRUE)
  if (which=="naive") {
    variances$naive <- TRUE
    variances$predVar <- FALSE
    variances$disp <- FALSE # we only want the "naive" component
  } else if (which=="respVar") {
    variances$respVar <- TRUE
  } else if (which=="residVar") {
    variances$residVar <- TRUE
  } else if (which=="fixefVar") {
    variances$fixefVar <- TRUE
  } else if (which=="intervals") {
    if (is.null(mc$intervals)) mc$intervals <- "predVar" # but other intervals can be obtained if mc$intervals is not NULL
  } else variances$predVar <- TRUE
  mc$variances <- variances
  mc[[1L]] <- get("predict.HLfit", asNamespace("spaMM")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  attr(eval(mc,parent.frame()),which)
}

get_residVar <- function(...) {
  get_predVar(...,which="residVar")
}

get_respVar <- function(...) {
  get_predVar(...,which="respVar")
}

get_intervals <- function(...) {
  get_predVar(...,which="intervals")
}

get_fixefVar <- function(...) {
  get_predVar(...,which="fixefVar")
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

get_ranPars <- function(object, which=NULL, ...) {
  CorrEst_and_RanFix <- object$CorrEst_and_RanFix
  if (is.null(which)) {
    return(CorrEst_and_RanFix)
  } else if (which=="corrPars") {
    resu <- CorrEst_and_RanFix$corrPars
    if ( ! is.null(resu)) resu <- structure(resu, type=attr(CorrEst_and_RanFix,"type")$corrPars)
    return(resu)
  } else if (which=="lambda") {
    resu <- CorrEst_and_RanFix$lambda
    if ( ! is.null(resu)) resu <- structure(resu, type=attr(CorrEst_and_RanFix,"type")$lambda)
    return(resu)
  } 
}

.get_from_ranef_info <- function(object,which="sub_corr_info") {
  if (is.null(object$ranef_info)) { ## occurs for < v2.6.15
    if (is.null(object[[which]])) { ## occurs for < v2.4.57
      return(object$corr_info) 
    } else return(object[[which]])
  } else return(object$ranef_info[[which]])
}

# older version of formula.HLfit removed from  [v2.6.59
formula.HLfit <- function(x, which="hyper", ...) {
  ## stats:::formula.default looks for x$formula then x$call$formula. 
  # So formula(object) should be enough, EXCEPT that if it finds neither (no explicitly named formula in the call), 
  # it evaluates the call, in which case print(formula(<HLfit>)) leads to an infinite recursion
  # since form is then an HLfit object so print(form...) will call summary.HLfit()...
  form <- NULL
  if (x$spaMM.version> "2.6.20") form <- x$predictor
  # Older object or objects that are not strictly HLfit objects, as detailed below:
  if (is.null(form) && x$spaMM.version> "2.4.35") form <- x$call$formula
  # Objects that are not strictly hLfit objects, as follows:
  ## finds something for HLfit member object in an HLfitlist  (they have a call with $processed and no $formula 
  #          bc it's not clear where to get a non-processed HLCor or HLfit call if the oricall was an outer optimizing fn)
  #    or    
  ## should not occur on a finished HLfit object but may on a call with $processed accessed in a debugging session
  #    or
  ## resid_fit !
  if (is.null(form)) form <- getCall.HLfit(x)$processed$predictor # from $call$processed$predictor
  #
  ## Now we must have a non-null form
  #
  if (which=="hyper") { # suitable for displays and for any function that reprocesses (update, refit, ... 
    #               and MFSDR -> stats::step() hence default="hyper")
    if ( ! is.null(resu <- attr(form,"hyper_info")$formula) ) return(resu)
  } else if (which=="no_offset") {
    if ( ! is.null(resu <- attr(form,"no_offset")) ) return(resu)
  }
  return(form) ## which = "" or NULL attributes
} ## extends stats::formula generic

nobs.HLfit <- function(object, ...) {length(object$y)}

get_any_IC <- function(object,...,verbose=interactive(),also_cAIC=TRUE) {
  info_crits <- .get_info_crits(object,also_cAIC=also_cAIC)
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

AIC.HLfit <- function(object, ..., k,verbose=interactive(),also_cAIC=TRUE) {
  get_any_IC(object,...,verbose=verbose,also_cAIC=also_cAIC)
}

extractAIC.HLfit <- function(fit, scale, k=2L, ..., verbose=FALSE) { ## stats::extractAIC generic
  df <- fit$dfs[["pforpv"]]
  aic <- AIC(object=fit, ..., verbose = verbose,also_AIC=FALSE)[["       marginal AIC:"]] # does not use k
  if (k !=2L) aic <- aic + (k - 2)*df
  c(edf=df, AIC=aic) 
}

.get_XZ_0I <- function(object) { ## there is a .calc_XZ_0I from the processed AUGI0_ZX
  pforpv <- ncol(object$X.pv)
  nrd <- length(object$w.ranef)
  nobs <- nrow(object$X.pv)
  ZAL <- get_ZALMatrix(object, force_bind=TRUE) # force bind for rbind2()  
  if (is.null(ZAL)) {
    XZ_0I <- object$X.pv ## and general code below works and gives the same result for augmented or not
  } else {
    XZ_0I <- cbind2(
      rbind2(object$X.pv, matrix(0,nrow=nrd,ncol=pforpv)), 
      rbind2(ZAL, diag(nrow=nrd))
    ) 
  }
  return(XZ_0I)
}

get_ZALMatrix <- function(object,as_matrix, force_bind=FALSE) {
  if ( ! missing(as_matrix)) stop("'as_matrix' is deprecated")
  if (length(ZAlist <- object$ZAlist)) { ## ou tester if (object$models[["eta"]]=="etaGLM")
    if (is.null(object$envir$ZALMatrix) ||
        (force_bind && inherits(object$envir$ZALMatrix,"ZAXlist")) ) {
      object$envir$ZALMatrix <- .compute_ZAL(XMatrix=object$strucList, ZAlist=ZAlist,as_matrix=FALSE, force_bind=force_bind) 
    }
    return(object$envir$ZALMatrix)
  } else return(NULL) 
}

.get_ZAfix <- function(object, as_matrix) {
  if (is.null(object$envir$ZAfix)) object$envir$ZAfix <- .ad_hoc_cbind(object$ZAlist, as_matrix=as_matrix)  
  return(object$envir$ZAfix)
}

# left pseudo inverse of (augmented) design matrix, (X_a' W X_a)^{-1} X_a' W
.get_WLS_ginv <- function(object, augmented, XZ_0I=NULL) { 
  ## gets inv(tX_a invSig_a X_a).tX_a invSig_a that gives hat(beta,v_h)
  if (is.null(XZ_0I))  XZ_0I <- .get_XZ_0I(object)
  ww <- c(object$w.resid, object$w.ranef) ## NOT sqrt()
  Wei_XZ_0I <- .calc_wAugX(XZ_0I=XZ_0I, sqrt.ww=ww) ## 2nd argument name misleading
  if (is.null(tcrossfac_beta_v_cov <- object$envir$beta_cov_info$tcrossfac_beta_v_cov)) {
    tcrossfac_beta_v_cov <- .get_beta_cov_info(object)$tcrossfac_beta_v_cov
  }
  beta_v_cov <- .tcrossprod(tcrossfac_beta_v_cov)
  augXWXXW <- tcrossprod(beta_v_cov, Wei_XZ_0I) ## = solve(crossprod(wAugX)) %*% crossprod(wAugX, diag(x=sqrt.ww))
  if (augmented) {
    return(augXWXXW)
  } else {
    return(augXWXXW[seq_len(ncol(object$X.pv)),seq_len(nrow(object$X.pv)),drop=FALSE])
  }
}

.get_control_dist <- function(fitobject, char_rd) {
  if (is.null(optimInfo <- attr(fitobject,"optimInfo"))) {
    # HLCor, or fitme with nothing to optimize... (corrHLfit with nothing to optimize has an optimInfo...)
    outer_call <- getCall(fitobject) 
    return(outer_call$dist.method[[char_rd]]) ## (from a preprocessed version copied back to the call) ## imposes char_rd here
  } else return(optimInfo$LUarglist$moreargs[[char_rd]]$control.dist)
} # F I X M E but the two are not equivalent since the moreargs version has gone through .provide_rho_mapping() which converts NULL rho.mapping into explicit rho_mappings

get_matrix <- function(object, which="model.matrix", augmented=TRUE, ...) {
  switch(which,
         "model.matrix"= object$X.pv, #model.matrix(object, ...), ## X
         "ZAL"=get_ZALMatrix(object),                             ## ZAL
         "AugX"=.get_XZ_0I(object),                               ## X_a
         "wAugX"={
           XZ_0I <- .get_XZ_0I(object)
           ww <- c(object$w.resid, object$w.ranef)
           .calc_wAugX(XZ_0I=XZ_0I, sqrt.ww=sqrt(ww)) 
         },                               
         "wei_AugX"={
           XZ_0I <- .get_XZ_0I(object)
           ww <- c(object$w.resid, object$w.ranef) ## NOT sqrt()
           Wei_XZ_0I <- .calc_wAugX(XZ_0I=XZ_0I, sqrt.ww=ww) ## 2nd argument name misleading
         },                               
         "left_ginv"= .get_WLS_ginv(object, augmented=augmented), ## X_a^- = (X_a' W X_a)^{-1} X_a' W
         "hat_matrix"= { ## hat projection matrix                 ## P_a = X_a X_a^- = X_a (X_a' W X_a)^{-1} X_a' W
           XZ_0I <- .get_XZ_0I(object) 
           WLS_ginv <- .get_WLS_ginv(object, augmented=TRUE, XZ_0I=XZ_0I)
           XZ_0I %*% WLS_ginv 
          },
         stop("Unhandled 'which' value in get_matrix()")
  )
}

model.matrix.HLfit <- function(object, ...) object$X.pv

"how" <- function(object, ...) UseMethod("how")

how.default <- function(object, ...) {message(paste("No 'how' method defined for objects of class",class(object)))} 

how.HLfit <- function(object, devel=FALSE, verbose=TRUE, format=print, ...) {
  info <- object$how
  if (is.null(info)) {
    info <- list(MME_method=setdiff(object$MME_method,c("list")),
                 fit_time=object$fit_time,
                 "spaMM.version"=object$spaMM.version)
  }
  if (verbose) {
    fun <- getCall(object)[[1L]]
    if (is.function(fun)) { # from do.call(spaMM::fitme, args = args) => it's a closure
      fnname <- names(which(sapply(list(fitme=spaMM::fitme,
                                        HLfit=spaMM::HLfit,
                                        HLCor=spaMM::HLCor,
                                        corrHLfit=spaMM::corrHLfit), identical, y=fun)))
    } else { # assuming it's a 'name' (direct call or do.call("fitme", args = args)) or fn got by by get(...)
      fnname <- paste(sub("^.*?::","",fun)) ## remove any "spaMM::" in the name, from which paste() would return c("::","spaMM",<>)
    }
    if ( ! length(fnname)) {
      format(paste0("Model fitted by spaMM, version ",info[["spaMM.version"]],
                   ", in ",info$fit_time,"s using method: ",paste(info$MME_method,collapse=","),"."))
      if (packageVersion("spaMM")==info[["spaMM.version"]]) {
        message("(Fitting function could not be identified)") # and it's not clear why
      } else { # object fitted with different version of spaMM, identical() may not work
        # function could not be identified, but we know why
      }
    } else format(paste0("Model fitted by spaMM::", paste(fnname),", version ",info[["spaMM.version"]],
                 ", in ",info$fit_time,"s using method: ",paste(info$MME_method,collapse=","),"."))
  }
  if  (!  is.null(resid_fit <- object$resid_fit)) {
    resid_info <- object$resid_fit$how
    if (is.null(resid_info)) {
      resid_info <- list( MME_method=setdiff(object$resid_fit$MME_method,c("list")) )
    }
    if (verbose) format(paste0("Residual model fitted using method: ",paste(resid_info$MME_method, collapse=","),"."))
    info$resid_info <- resid_info
  }
  if (devel && ! is.null(switches <- info$switches)) {
    if (verbose) format(paste(paste(names(switches),"=",switches), collapse = ", "))
  }  
  invisible(info)
}

response <- function(object, ...) object$y[,1L]

family.HLfit <- function(object, ...) object$family

model.frame.HLfit <- function(formula, ...) {
  object <- formula
  form <- formula(object) 
  form <- .asNoCorrFormula(form) ## strips out the spatial information, retaining the variables
  frame.form <- .subbarsMM(form) ## this comes from lme4 and converts (...|...) terms to some "+" form 
  environment(frame.form) <- environment(form)
  mf_call <- call("model.frame", data=object$data, formula=frame.form,  drop.unused.levels=TRUE) # language object reused later
  eval(mf_call) ## data.frame with all vars required for fixef and for ranef, and a "terms" attribute
} ## model frame -> data.frame; model.matrix -> matrix...

## Might eventually try to match the merMod version:
#   str(model.frame(fittedModel))
#   'data.frame':	200 obs. of  3 variables:
#     $ observedResponse: int  1 1 1 1 0 0 1 0 1 1 ...
#   $ Environment1    : num  -0.776 0.119 -0.29 0.527 0.726 ...
#   $ group           : Factor w/ 10 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
#   - attr(*, "terms")=Classes 'terms', 'formula'  language observedResponse ~ Environment1 + (1 + group)
#   .. ..- attr(*, "variables")= language list(observedResponse, Environment1, group)
#   .. ..- attr(*, "factors")= int [1:3, 1:2] 0 1 0 0 0 1
#   .. .. ..- attr(*, "dimnames")=List of 2
#   .. .. .. ..$ : chr [1:3] "observedResponse" "Environment1" "group"
#   .. .. .. ..$ : chr [1:2] "Environment1" "group"
#   .. ..- attr(*, "term.labels")= chr [1:2] "Environment1" "group"
#   .. ..- attr(*, "order")= int [1:2] 1 1
#   .. ..- attr(*, "intercept")= int 1
#   .. ..- attr(*, "response")= int 1
#   .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#     .. ..- attr(*, "predvars")= language list(observedResponse, Environment1, group)
#   .. ..- attr(*, "dataClasses")= Named chr [1:3] "numeric" "numeric" "factor"
#   .. .. ..- attr(*, "names")= chr [1:3] "observedResponse" "Environment1" "group"
#   .. ..- attr(*, "predvars.fixed")= language list(observedResponse, Environment1)
#   .. ..- attr(*, "varnames.fixed")= chr [1:2] "observedResponse" "Environment1"
#   .. ..- attr(*, "predvars.random")= language list(observedResponse, group)
#   - attr(*, "formula")=Class 'formula'  language observedResponse ~ Environment1 + (1 | group)
#   .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#     
#     
# starting from 
#   'data.frame':	200 obs. of  3 variables:
#     $ observedResponse: int  1 1 1 1 0 0 1 0 1 1 ...
#   $ Environment1    : num  -0.776 0.119 -0.29 0.527 0.726 ...
#   $ group           : Factor w/ 10 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
#   - attr(*, "terms")=Classes 'terms', 'formula'  language observedResponse ~ Environment1 + (1 + group)
#   .. ..- attr(*, "variables")= language list(observedResponse, Environment1, group)
#   .. ..- attr(*, "factors")= int [1:3, 1:2] 0 1 0 0 0 1
#   .. .. ..- attr(*, "dimnames")=List of 2
#   .. .. .. ..$ : chr [1:3] "observedResponse" "Environment1" "group"
#   .. .. .. ..$ : chr [1:2] "Environment1" "group"
#   .. ..- attr(*, "term.labels")= chr [1:2] "Environment1" "group"
#   .. ..- attr(*, "order")= int [1:2] 1 1
#   .. ..- attr(*, "intercept")= int 1
#   .. ..- attr(*, "response")= int 1
#   .. ..- attr(*, ".Environment")=<environment: 0x000000002871d6e8> 
#     .. ..- attr(*, "predvars")= language list(observedResponse, Environment1, group)
#   .. ..- attr(*, "dataClasses")= Named chr [1:3] "numeric" "numeric" "factor"
#   .. .. ..- attr(*, "names")= chr [1:3] "observedResponse" "Environment1" "group"

