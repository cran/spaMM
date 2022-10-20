.get_u_h <- function(object) {
  if (is.null(u_h <- object$u_h)) u_h <- attr(object$v_h, "u_h") # non-null $ranef is before v3.8.34
  u_h    
}

.inv_Lmatrix <- function(lmatrix, type=attr(lmatrix,"type"), regul.threshold) {
  invlmatrix <- NULL
  if (inherits(lmatrix,"dCHMsimpl")) { # before any test on type, because dCHMsimpl has a @type slot
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
  #
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
  invlmatrix
}

## Old comment: computes inv(L) [not inv(Corr): see .get_invColdoldList];
## => not so clear now that latent d's are removed
.get_invL_HLfit <- function(object, regul.threshold=1e-7) { 
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
          if ( ! is.null(latentL_blob <- attr(lmatrix,"latentL_blob"))) { ## from .process_ranCoefs
            if (is.null(compactchol_Q_w <- latentL_blob$compactchol_Q_w)) { # object$spaMM.version < "3.8.33"
              latent_d <- latentL_blob[["d"]] # will be needed for back compat whatever future devel
              compactchol_Q_w <- .Matrix_times_Dvec(latentL_blob$compactchol_Q, 1/sqrt(latent_d)) # tcrossfac of full precision matrix
            }
            if ( ! is.null(kron_Y <- object$ranef_info$sub_corr_info$kron_Y_LMatrices[[Lit]])) {
              # longL = solve(t(compactchol_Q_w)) \otimes Lunique
              # invL =solve(A=t(compactchol_Q_w) \otimes B=Lunique) = solve(A) \otimes solve(B)
              # see comments about in kron_Y_LMatrices in .get_invColdoldList()
              kron_Y <- .inv_Lmatrix(kron_Y, regul.threshold=regul.threshold)
            }
            invlmatrix <- .makelong(t(compactchol_Q_w),longsize=ncol(lmatrix),as_matrix=TRUE, kron_Y=kron_Y) # ___FIXME___ allow kron_long=FALSE ?
            # as_matrix=TRUE necessary for resu[u.range, u.range] <- invlmatrix
          } else invlmatrix <- .inv_Lmatrix(lmatrix, regul.threshold=regul.threshold)
          u.range <- (cum_n_u_h[Lit]+1L):(cum_n_u_h[Lit+1L])
          resu[u.range,u.range] <- invlmatrix   
        }
      }
      resu <- as(resu,"sparseMatrix") # previously broke test twolambda vs onelambda by effect on nearly singular matrix:
      # ...by changing the twolambda result.
      # # Without it I had two messages : 
      # 1: In .calc_logdisp_cov(object, dvdloglamMat = dvdloglamMat, dvdlogphiMat = dvdlogphiMat,  :
      # Numerical precision issue in computation of the information matrix for dispersion parameters:
      #   the prediction variance may be inaccurate.
      # 2: In .force_solve(logdispInfo) : The matrix looks exactly singular.
      # # While with it I had only the first message, and a result less consistent with onelambda
      # Alternatively to the present conversion, 
      # I could reproduce these two symptoms (only the first message, and a result less consistent) 
      # by calling .crossprodCpp_d in .calc_invV_factors() -> .crossprod(ZAfix, wrZ)
      # (as tested by hacking the return value of  the single call to .crossprod() in this get_predVar() test)
      # and the only impact of the .crossprod() is here to change the numerical precision of the $r_x_n element in the 
      # return value of .calc_invV_factors() by effects of order 1e-13 (this element being dgeMatrix whether .crossprodCpp_d was called or not).
      object$envir$invL <- resu
    }
    return(object$envir$invL) # May be Diagonal() => calling code must handle that case.
  }
}

## possibly useful comment from a removed function:
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
        # Ugly, BUT only back-compat code since now all latent_d's are nullified for post-fit computations: 
        if (! is.null(latent_d <- attr(lmatrix,"latentL_blob")[["d"]])) { # should be FALSE in current spaMM version
          # if I changed the meaning of strucList[[]] between the fit and the postfit (in .post_process_v_h_LMatrices)
          # I must adjust accordingly here. This is uglier.
          # and here we used the compact 'd' for testing and use the already expanded one next:
          newcoeffs[u.range] <- newcoeffs[u.range]/sqrt(object$envir$sXaug$AUGI0_ZX$envir$latent_d_list[[Lit]])   
        }
        newcoeffs[u.range] <- solve(t(lmatrix),newcoeffs[u.range])   ## newcoeffs must be a _vector_
      }
    }
  }
  return(newcoeffs)
}

# this must be for back compat since currently all fit object meet the first condition.
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

.residuals_bare <- function(type, y, mu, family, wts) {
  if (type=="deviance") {
    res <- sign(y-mu)*sqrt(pmax((family$dev.resids)(y, mu, wts), 0))
  } else if (type=="pearson") res <- (y - mu) * sqrt(wts)/sqrt(family$variance(mu)) # also for LLMs, cf e.g. Cribary-Neta & Zeileis
  drop(res)
}

residuals.HLfit <- function(object, type = c("deviance", "pearson", "response", "working", "std_dev_res"), force=FALSE, ...) {
  object <- .getHLfit(object)
  type <- match.arg(type)
  BinomialDen <- .get_BinomialDen(object) 
  if (is.null(BinomialDen)) BinomialDen <- 1L
  y <- object$y / BinomialDen ## le y 
  mu <- object$fv # on 0-1 probability scale for binomial models
  if (type=="response") {
    return(drop(y-mu))
  } else if (type=="working") { # residuals in glm.fit
    family <- object$family
    return(drop(y - mu)/family$mu.eta(family$linkfun(mu))) 
  } else if (type=="std_dev_res") {
    if (force || is.null(res <- object$std_dev_res[,1])) {
      std_dev_res <- .std_dev_resids(object, phi_est=residVar(object, which="phi"), 
                                     lev_phi=hatvalues(object, type="std"))$std_dev_res
      res <- (sign(y-mu) * std_dev_res)[,1]
    }
  } else { # deviance and pearson residuals.
    pw <- object$prior.weights
    family <- object$family
    if (is.null(family)) {
      wts <- .unlist(pw) * BinomialDen ## the BinomialDen are in the $prior.weights of a glm object, but not of an HLfit one
      families <- object$families
      res <- vector("list", length(families))
      cum_nobs <- attr(families,"cum_nobs")
      for (mv_it in seq_along(families)) {
        resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
        res[[mv_it]] <- .residuals_bare(type, y[resp_range], mu[resp_range], families[[mv_it]], wts[resp_range])
      }
      res <- unlist(res, recursive = FALSE, use.names=TRUE)
    } else {
      wts <- as.vector(pw) * BinomialDen ## the BinomialDen are in the $prior.weights of a glm object, but not of an HLfit one
      res <- .residuals_bare(type, y, mu, family, wts)
      #if (!is.null(object$na.action))  res <- naresid(object$na.action, res)
    }
  }
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
  } else uv_h <- .get_u_h(object) #random effects \eqn{u}
  cum_n_u_h <- attr(.get_u_h(object),"cum_n_u_h")
  if (object$spaMM.version < "2.2.116") {
    ranefs <- attr(object$ZAlist,"ranefs") 
  } else ranefs <- attr(object$ZAlist,"exp_ranef_strings") 
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
        if (inherits(strucList[[it]], "dCHMsimpl")) {
          res <- as.vector(solve(strucList[[it]], res, system="Lt"))
        } else res <- as.vector(strucList[[it]] %*% res) ## vector
      } 
      names(res) <- colNames[[it]]
    }
    RESU[[it]] <- res
  }
  names(RESU) <- ranefs
  class(RESU) <- c("ranef", "list")
  RESU ## __F I X M E__ TODO: ~lme4:::ranef.merMod(mod, condVar = TRUE) & ajouter des arguments "variances" et "intervals" à ta fonction ranef() (Alex 7/5/2018)
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

.get_phi_fit <- function(object, mv_it=NULL) {
  phi_model <- object$models[["phi"]]
  if (is.null(mv_it)) {
    phi_fit <- switch(phi_model,
                      "phiGLM" = object$phi.object$glm_phi, ## glm
                      "phiHGLM" = object$resid_fit, ## hlfit
                      "phiScal" = object$phi, ## scalar
                      # "" = object$phi, ## scalar for the count families, or user-given phi  # but switch does not handle ""
                      # stop('Unhandled object$models[["phi"]]')
                      object$phi ## scalar for the count families, or user-given phi
    )  
  } else {
    phi_fit <- switch(phi_model[[mv_it]],
                      "phiGLM" = object$phi.object[[mv_it]]$glm_phi, ## glm
                      "phiHGLM" = object$resid_fits[[mv_it]], ## hlfit        
                      "phiScal" = object$phi[[mv_it]], ## scalar 
                      object$phi[[mv_it]] ## scalar for the count families, or user-given phi...
    )  
  }
  return(phi_fit)
}

#------------------------------       get_residVar vs residVar (or simulate)      -------------------------
# get_residVar() (anything that leads to variances$residVar being TRUE)
# |-> .predict_body()
#     |-> .add_residVar()
#         |-> .warn_pw() + **.calcResidVar()**: no pw computations
#                           |-> .get_glm_phi() + predict( ,newdata) + Gamma-specif code
#                                     ^ 
#                                     |
#                                    v
#        |-> uses hard-coded pw + [.get_glm_phi() or .get_phi_fit()] + predict(, newdata allowed)
# |-> **.get_phiW(, newdata)** + Gamma-specif code
# residVar() or simulate()
#
# The source code of .calcResidVar_phiW() (not actually used) shows how we could substitute .get_phiW() to some of its code,
# and illustrates with this is not such a good idea.
#__________________________________________________________________________________________________________

residVar <- function(object, which="var", submodel=NULL, newdata=NULL) {
  phi_model <- object$models[["phi"]]
  if (which=="formula") {
    if (length(phi_model)>1L && is.null(submodel)) stop("'submodel' index required to extract residual-model formula.") 
    .get_phiform(object, mv_it=submodel)
  } else if (which=="family") {
    if (length(phi_model)>1L && is.null(submodel)) stop("'submodel' index required to extract residual-model family") 
    eval(.get_phifam(object, mv_it=submodel))
  } else if (which=="fam_parm") {
    family <- object$family
    if (is.null(family)) {
      if (is.null(submodel)) stop("'submodel' index required to extract residual-model family parameter") 
      family <- object$families[[submodel]]
    }
    if (family$family=="COMPoisson") {
      return(environment(family$aic)$nu)
    } else if (family$family %in% c("beta_resp")) {
      return(environment(family$aic)$prec)
    } else if (family$family %in% c("negbin","negbin1","negbin2")) {
      return(environment(family$aic)$shape)
    } else return(NA)
  } else if (which=="fit") {
    if (length(phi_model)>1L && is.null(submodel)) stop("'submodel' index required to extract residual-model fit") 
    #                                          # We could return a list of objects but too heterogeneous to be convenient at user level.
    .get_phi_fit(object, mv_it=submodel)
  } else if (which %in% c("phi","var")) {
    mu <- predict(object, newdata=newdata) # 1-col matrix
    cum_nobs <- attr(object$families,"cum_nobs")
    if ( ! is.null(submodel)) {
      resp_range <- .subrange(cumul=cum_nobs, it=submodel)
      phi <- .get_phiW(object=object, newdata=newdata, dims=c(length(resp_range), 1L),
                       phi_type="predict",prior.weights=object$prior.weights[[submodel]],
                       phi_mod_class=phi_model[submodel], mv_it=submodel)
      if (which=="var" && object$families[[submodel]]$family=="Gamma") {
        return(as.vector(phi * mu[resp_range]^2)) # var
      } else return(as.vector(phi)) # phi
    } else {
      phi <- .get_phiW(object=object, newdata=newdata, 
                            dims=dim(mu), phi_type="predict")
      if (which=="var") {
        if ( ! is.null(cum_nobs)) {
          var <- phi
          for (mv_it in seq_len(length(cum_nobs)-1L)) {
            if (object$families[[mv_it]]$family=="Gamma") {
              resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
              var[resp_range] <- var[resp_range] * mu[resp_range]^2
            }
          }
          return(as.vector(var)) # var
        } else {
          if (object$family$family=="Gamma") {
            return(phi * mu^2) # var
          } else return(as.vector(phi)) # var
        }
      } else return(as.vector(phi)) 
    }
  } else stop("Unknown 'which' type")
} 

# __F I X M E__: from R v3.3.0: New S3 generic function sigma() with methods for extracting the estimated standard deviation aka “residual standard deviation” from a fitted model. 

.get_phiW <- function(object, newdata, 
                      dims, 
                      prior.weights=object$prior.weights, phi_type, needed,
                      phi_mod_class=object$models[["phi"]],
                      mv_it=NULL
                      ) { # returns phi/prior.weights
  if (length(phi_mod_class)>1L) {
    cum_nobs <- attr(object$families,"cum_nobs")
    phiWs <- vector("list", length(phi_mod_class)) # will be a list of MATRICES
    for (mv_it in seq_along(phiWs)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      phiWs[[mv_it]] <- .get_phiW(object, newdata, 
                                  dims=c(length(resp_range),dims[2]), 
                                  prior.weights=prior.weights[[mv_it]], phi_type, needed,
                                  phi_mod_class=phi_mod_class[mv_it],
                                  mv_it=mv_it)
    }
    vec_is_phiW_fix_btwn_sims <- lapply(phiWs, attr, which="is_phiW_fix_btwn_sims")
    vec_is_phiW_fix_btwn_sims <- unlist(vec_is_phiW_fix_btwn_sims, recursive=FALSE, use.names = FALSE)
    return(structure(do.call(rbind,phiWs), # the phiWs should be 1- or nsim-col matrices
                     is_phiW_fix_btwn_sims=vec_is_phiW_fix_btwn_sims))
  }
  is_phiW_fix_btwn_sims <- FALSE
  if (phi_mod_class == "") { ## the count families, or user-given phi
    phi_fit <- .get_phi_fit(object, mv_it=mv_it)
    if (length(phi_fit)>1L) {
      if (is.null(newdata)) {
        message(paste0("simulate.HLfit() called on an original fit where phi was given but not constant.\n",
                       "This phi will be used, but is that relevant?"))
      } else stop("I do not know what to simulate when 'newdata' is not NULL and the original fit's phi was given but not constant.")
    }  
    newphiMat <- matrix(phi_fit, nrow=dims[1], ncol=dims[2]) # n_repl is typically 1, or nsim
    if (identical(attr(prior.weights,"unique"),TRUE)) {
      phiW <- newphiMat/prior.weights[1L]
      is_phiW_fix_btwn_sims <- TRUE
    } else phiW <- .Dvec_times_matrix(1/prior.weights, newphiMat)  ## warnings or errors if something suspect
  } else {
    if (phi_mod_class=="phiGLM") {
      phi_fit <- .get_glm_phi(object, mv_it=mv_it) # object of class "glm"
    } else phi_fit <- .get_phi_fit(object, mv_it=mv_it) # more diverse object
    if (phi_type=="predict") {
      newphiVec <- switch(phi_mod_class,
                          "phiGLM" = predict(phi_fit, newdata=newdata, type="response"), ## vector
                          "phiHGLM" = predict(phi_fit, newdata=newdata, type="response")[ ,1L],
                          "phiScal" = rep(phi_fit, dims[1]),
                          stop('Unhandled object$models[["phi"]]')
      ) ## VECTOR in all cases, becomes matrix later
      if (identical(attr(prior.weights,"unique"),TRUE)) {
        phiW <- newphiVec/prior.weights[1L]
        if (phi_mod_class %in% c("phiScal","phiGLM")) is_phiW_fix_btwn_sims <- TRUE
      } else phiW <- newphiVec/prior.weights  ## warnings or errors if something suspect
      phiW <- matrix(phiW,nrow=length(phiW), ncol=dims[2])  # vector -> matrix
    } else { # any other phi_type 
      newphiMat <- switch(phi_mod_class,
                          "phiGLM" = as.matrix(simulate(phi_fit, newdata=newdata, nsim=needed)), ## data frame -> matrix
                          "phiHGLM" = simulate(phi_fit, newdata=newdata, type=phi_type, nsim=needed),
                          "phiScal" = matrix(phi_fit,nrow=dims[1],ncol=dims[2]),
                          stop('Unhandled object$models[["phi"]]')
      ) ## already MATRIX in all cases
      if (identical(attr(prior.weights,"unique"),TRUE)) {
        phiW <- newphiMat/prior.weights[1L]
        if (phi_mod_class=="phiScal") is_phiW_fix_btwn_sims <- TRUE
      } else phiW <- .Dvec_times_matrix(1/prior.weights,newphiMat)  ## warnings or errors if something suspect
      # F I X M E add diagnostics ?
    } # phiW is always a matrix
  }
  attr(phiW,"is_phiW_fix_btwn_sims") <- is_phiW_fix_btwn_sims
  return(phiW) ## always MATRIX
}

.get_objLik <- function(object) { # conceived to work on both 'processed' and HLfit objects
  mess <- .REMLmess(object)
  which <- switch(mess, 
                  "by stochastic EM."= "logLapp",
                  "by Laplace ML approximation (p_v)."= "p_v", # old objects before distinction p|p^\circ
                  "by ML (p_v approximation of logL)."= "p_v",
                  "by ML (P_v approximation of logL)."= "P_v",
                  "by h-likelihood approximation."= "p_v",
                  "by ML."= "p_v",
                  "by Laplace REML approximation (p_bv)."= "p_bv", # old objects before distinction p|^\circ
                  "by REML (p_bv approximation of restricted logL)."= "p_bv",
                  "by REML (P_bv approximation of restricted logL)."= "P_bv",
                  "by REML."= "p_bv",
                  "by non-standard REML"= "p_bv",
                  stop(paste0("No default '",which,"' value for '",mess,"' estimation method."))
  ) 
}


logLik.HLfit <- function(object, which=NULL, ...) {
  object <- .getHLfit(object)
  if (is.null(which)) which <- .get_objLik(object)
  if (which=="cliks") { ## *undoc* and the summand=TRUE non-default has no effect when clik_fn uses family()$aic (__F I X M E__) 
    if (object$family$family=="binomial") {
      muFREQS <- predict(object)
      mu <- muFREQS * object$BinomialDen
    }
    clik_fn <- .get_clik_fn(object$family)
    phi_est <- .get_phiW(object=object, newdata=NULL, dims=c(length(mu),1L), phi_type="predict", needed=1L,...)[,1L] # ... allows to overcome the default prior.weights
    reweiclik <- .calc_clik(mu=mu,phi_est=phi_est,object, clik_fn=clik_fn, summand=TRUE, ...) # ... allows to overcome the default prior.weights
    colnames(reweiclik) <- "clik"
    return(reweiclik) ## currently 1-col matrix; avoids names(resu) <- which !
  } else resu  <- object$APHLs[[tolower(which)]]
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
    beta_cov <- object$envir$"beta_cov" # special version assigned by .init_promises_spprec(., nullify_X_ones =TRUE) 
    if (is.null(beta_cov)) beta_cov <- object$envir$beta_cov_info$beta_cov
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
    if (is.null(resu)) resu <- .tcrossprod(strucList[[it]],NULL, perm=TRUE)
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


.add_varCorr_phi_lines <- function(x, loctable, corrFill, families=x$families, 
                         family=x$family, phi.object=x$phi.object, phimodel=x$models[["phi"]], mv_it=NULL) {
  if ( ! is.null(families)) {
    for (mv_it in seq_along(families)) {
      loctable <- .add_varCorr_phi_lines(x=x, loctable=loctable, corrFill=corrFill, families=NULL, 
                               family=families[[mv_it]], phi.object=phi.object[[mv_it]], phimodel=phimodel[mv_it],
                               mv_it=mv_it)
    }
    return(loctable)
  }
  if (family$family %in% c("gaussian","Gamma")) {
    if ( ! is.null(phi_outer <- phi.object$phi_outer)) { 
      phi_line <- data.frame(Group="Residual",Term="(Intercept)",Variance=phi_outer, "Std.Dev."=sqrt(phi_outer))
      if ("Corr." %in% colnames(loctable)) phi_line <- cbind(phi_line,corrFill, row.names=NULL)
      loctable <- rbind(loctable,phi_line)
    } else {
      if (phimodel=="phiHGLM") { 
        #cat("Residual dispersion model includes random effects:\n  use summary(<fit object>$resid_fit) to display results.\n")       
      } else if ((loc_p_phi <- length(phi.object$fixef))==1L) {
        namesX_disp <- names(phi.object$fixef)
        phiform <- .get_phiform(x, mv_it) 
        dispOffset <- attr(phiform,"off") 
        if (!is.null(dispOffset)) dispOffset <- unique(dispOffset)
        if (length(namesX_disp)==1 && namesX_disp[1]=="(Intercept)" && length(dispOffset)<2L) {
          # constant phi: we can display it
          phi_est <- (phi.object$fixef)
          if (length(dispOffset)==1L) phi_est <- phi_est+dispOffset
          resid.family <- eval(.get_phifam(x, mv_it))  
          phi_est <- resid.family$linkinv(phi_est)
        } ## else phi not constant; We don't try to display it
        if (is.null(mv_it)) { 
          grptxt <- "Residual"
        } else grptxt <- paste0("Residual_",mv_it)
        phi_line <- data.frame(Group=grptxt,Term="(Intercept)",Variance=phi_est, "Std.Dev."=sqrt(phi_est))
        if ("Corr." %in% colnames(loctable)) phi_line <- cbind(phi_line, corrFill, row.names=NULL)
        loctable <- rbind(loctable,phi_line)
      }                                                 
    }
  }
  return(loctable)
}

# for a lme4::VarCorr() equivalent; generic is nlme::VarCorr 
VarCorr.HLfit <- function(x, sigma=1, add_residVars=TRUE, verbose=TRUE, ...) {
  loctable <- NULL
  if ( ! is.null(lambda.object <- x$lambda.object)) {
    #.legend_lambda(object, type = "family")
    namesTerms <- lambda.object$print_namesTerms ## list of vectors of variable length
    linklam_coeff_list <- lambda.object$coefficients_lambdaS ## used beyond the next line
    lamtable <- .lambda_table_fn(namesTerms, x, lambda.object,linklam_coeff_list)
    nonunique_colnames <- colnames(lamtable)
    loctable <- lamtable[,seq_len(ncol(lamtable))] # subsetting automatically generates unique names for the Corr. columns, but attributes are dropped
    for (it in seq_len(nrow(loctable))) {
      # That's not good bc the two lines for Group 'gridcode' are filled with the Intercept ($lambda_list never contains the adjd...)
      # if (is.na(lamtable[[it,"Var."]])) lamtable[[it,"Var."]] <- lambda.object$lambda_list[[lamtable[[it,"Group"]]]][[lamtable[[it,"Term"]]]]
      if (is.na(loctable[[it,"Var."]])) {
        if ((term. <- loctable[[it,"Term"]])=="(Intercept)") {
          loctable[[it,"Var."]] <- lambda.object$lambda_list[[loctable[[it,"Group"]]]][[loctable[[it,"Term"]]]] # not a glm coef (cf inverse link)
        } # else do nothing bc it's not clear what to do
      }
    }
    loctable <- data.frame(Group=loctable[,"Group"],Term=loctable[,"Term"],Variance=loctable[,"Var."],"Std.Dev."=sqrt(loctable[,"Var."]))
    if ("Corr." %in% colnames(lamtable)) {
      corrFill <- lamtable[, nonunique_colnames=="Corr.", drop=FALSE]
      loctable <- cbind(loctable, corrFill, row.names=NULL)
      corrFill <- corrFill[1, , drop=FALSE]
      corrFill[] <- NA
    }
  } 
  if (add_residVars) loctable <- .add_varCorr_phi_lines(x, loctable, corrFill)
  rownames(loctable) <- NULL
  if (is.null(loctable)) {
    if (verbose) message("VarCorr() found no variance to report.")
  } else {
    attr(loctable,"random_slope_ncol_geq_1_pos") <- attr(lamtable,"random_slope_ncol_geq_1_pos") # indices are rows in table
    attr(loctable,"random_slope_ncol_geq_1_rows") <- attr(lamtable,"random_slope_ncol_geq_1_rows") # indices are rows in table
    attr(loctable,"row_map") <- attr(lamtable,"row_map") # map from ranef indices to rows
  }
  return(loctable)
} 


.dev_resids <- function(object, fv=object$fv, y=object$y, BinomialDen=object$BinomialDen, family=object$family, 
                        families=object$families, phi_est=NULL, lev_phi,...) {
  if ( ! is.null(families)) { # mv case, list of families
    cum_nobs <- attr(families,"cum_nobs")
    dev_res <- vector("list",length(families))
    fvs <- attr(fv,"mv")
    for (mv_it in seq_along(families)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      if (is.null(phi_est)) {
        dev_res[[mv_it]] <- .dev_resids(fv=fvs[[mv_it]], y=y[resp_range], BinomialDen=BinomialDen[resp_range], family=families[[mv_it]], 
                                        families=NULL, ...)
      } else  dev_res[[mv_it]] <- .dev_resids(fv=fvs[[mv_it]], y=y[resp_range], BinomialDen=BinomialDen[resp_range], family=families[[mv_it]], 
                                                   families=NULL, phi_est=phi_est[[mv_it]], lev_phi=lev_phi[resp_range],...)
    }
    if ( ! is.null(phi_est)) {
      matlist <- do.call("rbind", dev_res)
      dev_res <- list(std_dev_res=unlist(matlist[,"std_dev_res"], recursive = FALSE, use.names = FALSE), 
                      dev_res=unlist(matlist[,"dev_res"], recursive = FALSE, use.names = FALSE))
      
    } else dev_res <- unlist(dev_res, recursive = FALSE, use.names = FALSE)
  } else {
    if (family$family=="binomial") {
      dev_res <- family$dev.resids(y/BinomialDen, mu=fv, BinomialDen) # mu is proba for binomial
    } else {
      mu <- attr(fv,"mu_U") # mu of untruncated latent variable if it exists 
      if (is.null(mu)) mu <- fv # otherwise expectation of response
      dev_res <- family$dev.resids(y,mu=fv,wt=rep(1,length(fv)))
    }
    if ( ! is.null(phi_est)) dev_res <- list(std_dev_res=dev_res/(phi_est*(1-lev_phi)),
                                             dev_res=dev_res)
  }
  dev_res
}

dev_resids <- function(object, ...) .dev_resids(object, ...) # hides the default argumentt of .dev_resids()

.std_dev_resids <- function(object, phi_est, lev_phi, ...) .dev_resids(object, phi_est=phi_est,lev_phi=lev_phi, ...) # idem and more


deviance.HLfit <- function(object,...) {
  dev_res2 <- dev_resids(object=object,...)
  return(sum(dev_res2))
}  

.check_predVar_type_confusion <- local({
  predVar_type_warned <- FALSE
  function(variances, type="link", intervals=NULL) {
    if ( ! environment(.check_predVar_type_confusion)$predVar_type_warned &&
         identical(variances$predVar,TRUE) && type!="link" && is.null(intervals)) {
      # combination of options suggesting that users tried to get a predVar on response scale 
      warning(paste("Checking possible confusion in get_predVar() call:\n",
                    'variances$predVar is requested with explicit type != "link" but is always returned on linear predictor scale.\n',
                    "You may use the 'intervals' argument or get_intervals(., variances=list(predVar=TRUE))\n if you wish to translate linear predictor uncertainty on response scale."),
              call.=FALSE)
      predVar_type_warned <<- TRUE
    }
  }
})

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
    # if (is.null(mc$intervals)) mc$intervals <- "predVar" # but other intervals can be obtained if mc$intervals is not NULL
  } else variances$predVar <- TRUE
  mc$variances <- variances
  subcall <- mc[c(1L,which(names(mc) %in% c("variances","type","intervals")))]
  subcall[[1L]] <- get(".check_predVar_type_confusion", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  resu <- eval(subcall,parent.frame())
  mc[[1L]] <- get("predict.HLfit", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  resu <- eval(mc,parent.frame())
  structure(attr(resu,which), respnames=attr(resu,"respnames"))
}

get_residVar <- function(...) {
  get_predVar(...,which="residVar")
}

get_respVar <- function(...) {
  get_predVar(...,which="respVar")
}

get_intervals <- function(..., intervals="predVar") {
  get_predVar(...,intervals=intervals, which="intervals")
}

get_fixefVar <- function(...) {
  get_predVar(...,which="fixefVar")
}



.get_RLRTSim_args <- function(object) { 
  if (object$family$family !="gaussian" 
      || (object$models[["eta"]]=="etaHGLM" && any(attr(object$rand.families,"lcrandfamfam")!="gaussian"))) {
    warning("'object' is not the fit of a LMM, while RLRTSim() methods were conceived for LMMs.")
  }
  X.pv <- as.matrix(object$X.pv)
  qrX.pv <- qr(X.pv)
  ZA <- .get_ZAfix(object, as_matrix=TRUE) 
  Lmat <- .get_LMatrix(object)
  corrmat <- tcrossprod(as.matrix(Lmat))
  sqrt_Sigma <- chol(corrmat) # RLRsim expects "The upper triangular Cholesky factor of the correlation matrix of the random effect"
  if(inherits(sqrt_Sigma,"Matrix")) sqrt_Sigma <- as.matrix(sqrt_Sigma)
  return(list(X=X.pv, Z=ZA, qrX = qrX.pv, sqrt.Sigma=sqrt_Sigma))
}

.get_LRTSim_args <- function(fullfit,nullfit) { 
  if (fullfit$family$family !="gaussian" 
      || (fullfit$models[["eta"]]=="etaHGLM" && any(attr(fullfit$rand.families,"lcrandfamfam")!="gaussian"))) {
    warning("'object' is not the fit of a LMM, while RLRTSim() methods were conceived for LMMs.")
  }
  X.pv <- as.matrix(fullfit$X.pv)
  ZA <- .get_ZAfix(fullfit, as_matrix=TRUE) 
  Lmat <- .get_LMatrix(fullfit)
  corrmat <- tcrossprod(as.matrix(Lmat))
  sqrt_Sigma <- chol(corrmat) # RLRsim expects "The upper triangular Cholesky factor of the correlation matrix of the random effect"
  if(inherits(sqrt_Sigma,"Matrix")) sqrt_Sigma <- as.matrix(sqrt_Sigma)
  return(list(X=X.pv, Z=ZA, q=fullfit$dfs[["pforpv"]]-nullfit$dfs[["pforpv"]], sqrt.Sigma=sqrt_Sigma))
}

get_RLRsim_args <- function(fullfit, nullfit, verbose=TRUE, REML=NA, ...) {
  if ( ! inherits(fullfit,"HLfit")) stop("'fullfit' should inherit from class 'HLfit'")
  is_REML <-  .REMLmess(fullfit,return_message=FALSE)
  if ( ! is.na(REML)) if (REML!=is_REML) warning("'object' is not a restricted likelihood fit. LRTSim() should be used next.")
  if (is_REML) {
    args <- .get_RLRTSim_args(fullfit)
    if (verbose) message("input for LRTSim() from package 'RLRsim' produced.")
  } else {
    if ( ! inherits(nullfit,"HLfit")) stop("'nullfit' should inherit from class 'HLfit'")
    args <- .get_LRTSim_args(fullfit,nullfit)
    if (verbose) message("input for RLRTSim() from package 'RLRsim' produced.")
  }
  return(args)
}

get_RLRTSim_args <- function(object, verbose=TRUE, ...) get_RLRsim_args(fullfit=object, verbose=verbose, REML=TRUE, ...)

# get_dispVar : dispVar in not a returned attribute

get_rankinfo <- function(object) return(attr(object$X.pv,"rankinfo")) 

.rm_Xi_ncol_geq_1 <- function(varcorr) {
  if ( ! is.null(row_map <- attr(varcorr,"row_map"))) {
    ranef_indices <- seq_along(row_map)
    if (  length(corr_rC_rows <- attr(varcorr,"random_slope_ncol_geq_1_rows"))) {
      varcorr <- varcorr[ - corr_rC_rows,,drop=FALSE] # previously this line was missing and the next code line
      ranef_indices <- ranef_indices[ - attr(varcorr,"random_slope_ncol_geq_1_pos")]
      ranef_map <- row_map[ranef_indices]
    } else   ranef_map <- row_map
    names(ranef_map) <- ranef_indices
    for (rd in ranef_indices) ranef_map[[rd]] <- rep(rd, length(ranef_map[[rd]]))
    # removed rows with a NA "Corr." value , which was nonsense. Now all Xi_ncol>1-ranCoefs are removed before the next code line.
    varcorr$ranef <- .unlist(ranef_map)
  } else varcorr$ranef <- seq(nrow(varcorr))
  varcorr
}

.get_fitted_ranPars <- function(object) {
  CorrEst_and_RanFix <- object$CorrEst_and_RanFix
  resu <- CorrEst_and_RanFix$corrPars
  resu <- .remove_from_cP(resu, u_names=names(which(unlist(attr(CorrEst_and_RanFix,"type")$corrPars)=="fix")))
  # Overwrite any fixed lambda or phi in all cases:
  if (length(lambda_list <- object$lambda.object$lambda_list)) {
    which_fitted_simple_lam <- which(sapply(lambda_list, length)==1L &  # removes (most) ranCoefs params [separate component of return value]
                                       object$lambda.object$type!="fixed") 
    resu$lambda <- .unlist(lambda_list[which_fitted_simple_lam]) 
    if (length(resu$lambda)) names(resu$lambda) <- which_fitted_simple_lam
  }
  resu$ranCoefs <-  .get_rC_inits_from_hlfit(object,type=c("inner","outer_ranCoefs")) # overlap with Xi_ncol ranCoefs?
  # ___F I X M E___ part of the pb is that of recovering info about the different types of rC from an HLfit object.
  resu
} 

.get_fittedPars_ran_phi <- function(object) { ## ___F I X M E___ ultimately a get_fittedPars extractor; 
                                              ## main pb may be heterogeneous structure of mv phi
  resu <- .get_fitted_ranPars(object)
  resu <- .add_famPars_outer(resu, fitobject=object)
  if ( is.null(object$families)) {
    phi <- residVar(object, which="fit") # the syntax to get a single scalar when it's appropriate; 
                                              #  object$phi.object$fittedPars may contain Gamma GLM coefs.
    if ( ! identical(attr(phi,"constr_phi"),TRUE)) resu$phi <- phi
  } else {
    phi <- list()
    for (mv_it in seq_along(object$families)) {
      phi_it <- residVar(object, which="fit", submodel = mv_it)
      if ( ! identical(attr(phi_it,"constr_phi"),TRUE)) phi[as.character(mv_it)] <- list(phi_it)
    }
    resu$phi <- phi
  }
  resu
} 

get_ranPars <- function(object, which=NULL, 
                        verbose=TRUE, # does not control all verbosity, only that which is TRUE by default.
                        lambda_names="Group.Term",
                        ...) {
  CorrEst_and_RanFix <- object$CorrEst_and_RanFix
  if (is.null(which)) {
    resu <- CorrEst_and_RanFix
    ## #resu <- .add_famPars_outer(parlist=resu, fitobject=object, type_attr=FALSE)  # adds only estimated ones. Are fixed ones in CorrEst_and_RanFix?
    varcorr <- VarCorr(object, add_residVars=FALSE, verbose=FALSE) 
    if ( ! is.null(varcorr)) {
      varcorr <- .rm_Xi_ncol_geq_1(varcorr) # removes (most) ranCoefs params [not part of return value]
      lambdatable <- na.omit(varcorr)
      if ( ! is.null(lambdatable)) {
        resu$lambda <- lambdatable[,"Variance"] 
        if (lambda_names=="Group.Term") {
          names(resu$lambda) <- paste(lambdatable[,"Group"],lambdatable[,"Term"],sep=".")
        } else names(resu$lambda) <- lambdatable$ranef
      }
    }
    return(resu)
  } else if (which=="fitted") {
    resu <- .get_fitted_ranPars(object)
    resu$phi <- NULL
    return(resu) # no type attribute
  } else if (which=="corrPars") {
    resu <- CorrEst_and_RanFix$corrPars
    if ( ! is.null(resu)) resu <- structure(resu, type=attr(CorrEst_and_RanFix,"type")$corrPars)
    return(resu)
  } else if (which=="lambda") { 
    varcorr <- VarCorr(object, add_residVars=FALSE, verbose=FALSE)
    if ( ! is.null(varcorr)) {
      if ( length(corr_rC_rows <- attr(varcorr,"random_slope_ncol_geq_1_rows"))) varcorr <- varcorr[-corr_rC_rows,,drop=FALSE] # previously this line was missing and the next code line
      # removed rows with a NA "Corr." value , which was nonsense. Now all Xi_ncol>1-ranCoefs are removed before the next code line.  
      varcorr <- na.omit(varcorr)
      varcorr <- .rm_Xi_ncol_geq_1(varcorr)
      resu <- structure(varcorr[,"Variance"], names=paste(varcorr[,"Group"],varcorr[,"Term"],sep="."))
      return(resu)
    }   # else returns NULL
  } else if (which=="outer_lambda") { 
    resu <- CorrEst_and_RanFix$lambda
    if ( ! is.null(resu)) resu <- structure(resu, type=attr(CorrEst_and_RanFix,"type")$lambda)
    return(resu)
  } else if (which=="ranef_var") {
    varcorr <- VarCorr(object, add_residVars=FALSE, verbose=verbose)
    if ( ! is.null(varcorr)) {
      varcorr <- .rm_Xi_ncol_geq_1(varcorr)
      lambdatable <- na.omit(varcorr)
      if ( ! is.null(lambdatable)) {
        lambda <- lambdatable[,"Variance"] 
        if (lambda_names=="Group.Term") {
          names(lambda) <- paste(lambdatable[,"Group"],lambdatable[,"Term"],sep=".")
        } else names(lambda) <- lambdatable$ranef
        resu <- list(Var=lambda,
                     outer=CorrEst_and_RanFix$lambda,
                     lambda_est=object$lambda.object$lambda_est, # long version: one value for each level of each ranef
                     lambda_list=object$lambda.object$lambda_list)
        attr(resu, "type") <- attr(CorrEst_and_RanFix,"type")$lambda
        return(resu)
      }   # else returns NULL
    } # else return NULL
  } else warning("'which' value not handled.")
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

nobs.HLfit <- function(object, ...) {length(object$y)} # see Details in help("nobs") for what it is useful for.

get_any_IC <- function(object, nsim=0L, ...,verbose=interactive(), also_cAIC=TRUE, short.names=NULL) {
  info_crits <- .get_info_crits(object,also_cAIC=also_cAIC, nsim=nsim, ...)
  info_crits <- info_crits[intersect(c("mAIC","cAIC","b_cAIC","dAIC","GoFdf"),names(info_crits))] # standard order better for display
  descriptive <- c(mAIC="       marginal AIC:",
                   dAIC="     dispersion AIC:",
                   GoFdf="       effective df:",
                   cAIC="    conditional AIC:",
                   b_cAIC="    conditional AIC:")
  comments <- c(mAIC  ="           ",
                dAIC  ="           ",
                GoFdf ="           ",
                cAIC  ="(plug-in)  ",
                b_cAIC="(bootstrap)")
  likelihoods <- uIC <- unlist(info_crits)
  names(likelihoods) <- descriptive[names(info_crits)]
  comments <- comments[names(info_crits)]
  if (verbose) {
    astable <- as.matrix(likelihoods)
    astable <- format(astable, justify="right")
    astable <- cbind(astable, comments, names(uIC))
    firstcolname <- paste0("         criterion", paste(substr("           ",1,nchar(astable[1,1])-3L)),"value")
    astable <- rbind(c(firstcolname, "  method   ", "short name"), astable)
    write.table(astable, col.names=FALSE, quote=FALSE) 
  }
  if (is.null(short.names)) short.names <-  ! is.null(info_crits$b_cAIC)
  if (short.names) {
    invisible(uIC)
  } else invisible(likelihoods)
}

AIC.HLfit <- function(object, ..., nsim=0L, k, verbose=interactive(), also_cAIC=TRUE, short.names=NULL) {
  dotlist <- list(...)
  ndots <- length(dotlist)
  if (ndots) { # there is an unnamed second argument or further unmatched arguments
    is_HLfit <- logical(ndots)
    for (it in seq_len(ndots)) is_HLfit[it] <- inherits(dotlist[[it]],"HLfit")
    which_is_HLfit <- which(is_HLfit)
    fitlist <- dotlist[which_is_HLfit]
    dotlist <- dotlist[which( ! is_HLfit)]
    if (nHLfits <- length(which_is_HLfit)) {
      dotexps <- as.list(substitute(...()))
      dotfitnames <- dotexps[which_is_HLfit]
      ICs <- vector("list", nHLfits+1L) 
      mc <- match.call(expand.dots = FALSE)
      mc[[1L]] <- get("get_any_IC", asNamespace("spaMM"), inherits=FALSE)
      mc["..."] <- NULL
      mc[["verbose"]] <- FALSE
      # With the new formals, all arguments beyond the fits should be named. But we make an exception for back compat with old formals:
      if (is.null(names(dotlist)[1L])) names(dotlist)[1L] <- "nsim" 
      mc[names(dotlist)] <- dotlist
      ICs[[1L]] <- eval(mc)
      colnams <- names(ICs[[1L]])
      objectname <- paste(mc[[2L]])
      for (it in seq_len(nHLfits)) {
        mc[["object"]] <- fitlist[[it]]
        ICs[[it+1L]] <- eval(mc)
      }
      
      ICs <- do.call("rbind.data.frame",ICs)
      rownames(ICs) <- c(objectname, dotfitnames)
      colnames(ICs) <- colnams
      if (verbose) {return(ICs)} else {invisible(ICs)}
    } else { # back -compatibility fix for case were second argument was unnamed, which was previously matched to the nsim argument
      mc <- match.call()
      names(mc)[[3L]] <- "nsim"
      mc[[1L]] <- get("get_any_IC", asNamespace("spaMM"), inherits=FALSE)
      eval(mc)
    }
  } else get_any_IC(object, nsim=nsim, ..., verbose=verbose, also_cAIC=also_cAIC, short.names=short.names) # no dots => no unnames second argument
}

extractAIC.HLfit <- function(fit, scale, k=2L, ..., verbose=FALSE) { ## stats::extractAIC generic
  df <- fit$dfs[["pforpv"]] # cf Value and Examples of extractAIC.HLfit showing in which sense this is the correct value.
  aic <- AIC(object=fit, ..., verbose = verbose, also_cAIC=FALSE, short.names=TRUE)[["mAIC"]] # does not use k
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

get_ZALMatrix <- function(object, force_bind=TRUE) {
  if (length(ZAlist <- object$ZAlist)) { ## ou tester if (object$models[["eta"]]=="etaGLM")
    if (is.null(object$envir$ZALMatrix) ||
        (force_bind && inherits(object$envir$ZALMatrix,"ZAXlist")) ) {
      object$envir$ZALMatrix <- .compute_ZAL(XMatrix=object$strucList, ZAlist=ZAlist,as_matrix=FALSE, 
                                             bind. = TRUE, force_bindable=force_bind) 
    }
    return(object$envir$ZALMatrix)
  } else return(NULL) 
}

.get_LMatrix <- function(object) {
  if (length(ZAlist <- object$ZAlist)) { ## ou tester if (object$models[["eta"]]=="etaGLM")
    if (is.null(object$envir$LMatrix)) {
      for (rd in seq_along(ZAlist)) ZAlist[[rd]] <- .symDiagonal(n = ncol(ZAlist[[rd]])) #  strucList elements may be null and then we need this.
      object$envir$LMatrix <- .compute_ZAL(XMatrix=object$strucList, ZAlist=ZAlist, as_matrix=FALSE, 
                                           bind.=TRUE, force_bindable=FALSE) # reproduces old defaults; __F I X M E__ set force_bindable=TRUE? 
    }
    return(object$envir$LMatrix)
  } else return(NULL) 
}


.get_ZAfix <- function(object, as_matrix) {
  if (.is_spprec_fit(object)) { 
    return(object$envir$sXaug$AUGI0_ZX$ZAfix) # using new (v3.2.19) spprec object$envir$sXaug
  }  else { # (( !spprec) || older spaMM) 
    if (is.null(object$envir$ZAfix)) object$envir$ZAfix <- .ad_hoc_cbind(object$ZAlist, as_matrix=as_matrix)  
    return(object$envir$ZAfix)
  }
}

# left pseudo inverse of (augmented) design matrix, (X_a' W X_a)^{-1} X_a' W
.get_WLS_ginv <- function(object, augmented, XZ_0I=NULL) { 
  ## gets inv(tX_a invSig_a X_a).tX_a invSig_a that gives hat(beta,v_h)
  if (is.null(XZ_0I))  XZ_0I <- .get_XZ_0I(object)
  ww <- c(.get_H_w.resid(object), object$w.ranef) ## NOT sqrt()
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

.get_fixef_WLS_ginv <- function(object, X.pv=object$X.pv) { 
  sqrt_ww <- sqrt(.get_H_w.resid(object))
  wX <- .Dvec_times_m_Matrix(sqrt_ww, X.pv) 
  if (inherits(wX,"sparseMatrix")) wX <- as.matrix(wX) # avoids suspect sparse-X code
  qrwX <- qr(wX)
  
  # if(inherits(wX,"sparseMatrix")) {
  #   solve(qrwX,diag(x=sqrt_ww))
  # } else {
  #   # .m_Matrix_times_Dvec(backsolve(qr.R(qrwX),t(qr.Q(qrwX))),sqrt_ww)
  qr.solve(qrwX,diag(x=sqrt_ww))
  # }
}


if (FALSE) { # permuted QR example - but singular too, so notexactly what we want
  X <- cbind(int = 1,
             b1=rep(1:0, each=3), b2=rep(0:1, each=3),
             c1=rep(c(1,0,0), 2), c2=rep(c(0,1,0), 2), c3=rep(c(0,0,1),2)
  )
  X <- as(X, "sparseMatrix")
  qX <- qr(X)
  drop0(R. <- qr.R(qX), tol=1e-13) # columns *permuted*: c3 b1 ..
  Q. <- qr.Q(qX)
  qI <- sort.list(qX@q) # the inverse 'q' permutation
  (X. <- drop0(Q. %*% R.[, qI], tol=1e-13))## 
  X-X.
  X %*% qr.coef(qX, diag(6)) %*% X # bad, not pseudoinverse
  solve(qX,diag(6)) %*% X # bad
  solve(R.[, qI] %*% t(Q.)) %*% X # bad
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
         "ZAL"=get_ZALMatrix(object, force_bind=TRUE),            ## ZAL
         "L"= .get_LMatrix(object),                               ## L
         "ZA"= .get_ZAfix(object,as_matrix=FALSE),                ## ZA   
         "AugX"=.get_XZ_0I(object),                               ## X_a
         "wAugX"={
           XZ_0I <- .get_XZ_0I(object)
           ww <- c(.get_H_w.resid(object), object$w.ranef)
           .calc_wAugX(XZ_0I=XZ_0I, sqrt.ww=sqrt(ww)) 
         },                               
         "wei_AugX"={
           XZ_0I <- .get_XZ_0I(object)
           ww <- c(.get_H_w.resid(object), object$w.ranef) ## NOT sqrt()
           Wei_XZ_0I <- .calc_wAugX(XZ_0I=XZ_0I, sqrt.ww=ww) ## 2nd argument name misleading
         },                               
         "left_ginv"= .get_WLS_ginv(object, augmented=augmented), ## X_a^- = (X_a' W X_a)^{-1} X_a' W
         "hat_matrix"= { ## hat projection matrix                 ## P_a = X_a X_a^- = X_a (X_a' W X_a)^{-1} X_a' W
           XZ_0I <- .get_XZ_0I(object) 
           WLS_ginv <- .get_WLS_ginv(object, augmented=TRUE, XZ_0I=XZ_0I)
           XZ_0I %*% WLS_ginv 
          },
         "fixef_left_ginv"= .get_fixef_WLS_ginv(object, ...), ## X^- = (X' W X)^{-1} X' W    # use the dots to pass an alternative X.pv
         stop("Unhandled 'which' value in get_matrix()")
  )
}

model.matrix.HLfit <- function(object, ...) object$X.pv

.prettify_method <- function(MME_method, by_y_augm) {
  if (by_y_augm) {
    for (it in seq_along(MME_method)) {
      MME_method[it] <- switch(MME_method[it],
                               AUGI0_ZX_spprec = "sparse-precision method for y-augmented matrix",
                               sXaug_Matrix_QRP_CHM_scaled = "sparse-correlation method for y-augmented matrix (Cholesky)", # get_absdiagR_blocks, not sXaug method
                               #sXaug_Matrix_cholP_scaled = "sparse-correlation method for y-augmented matrix (Cholesky)", # get_absdiagR_blocks, not sXaug method
                               AUGI0_ZX_sparsePrecision = "sparse-precision method for y-augmented matrix",
                               sXaug_EigenDense_QRP_Chol_scaled = "dense-correlation method for y-augmented matrix (chol with QR fallback)",
                               dgCMatrix = NA,
                               matrix = NA,
                               array = NA,
                               MME_method[it])
    }
  } else for (it in seq_along(MME_method)) {
    MME_method[it] <- switch(MME_method[it],
                             sXaug_Matrix_QRP_CHM_scaled = "sparse-correlation (QR) methods",
                             sXaug_Matrix_CHM_H_scaled = "sparse-correlation (Cholesky) methods",
                             AUGI0_ZX_spprec = "sparse-precision methods",
                             AUGI0_ZX_sparsePrecision = "sparse-precision methods",
                             sXaug_EigenDense_QRP_Chol_scaled = "dense-correlation (QR) methods",
                             #sXaug_Matrix_cholP_scaled = "sparse-correlation (Cholesky) methods (experimental)",
                             "(G)LM" = "methods for (G)LMs",
                             LLF = "methods for LLMs",
                             dgCMatrix = NA,
                             matrix = NA,
                             array = NA,
                             MME_method[it])
  }
  MME_method
}

"how" <- function(object, ...) UseMethod("how")

how.default <- function(object, ...) {message(paste("No 'how' method defined for objects of class",class(object)))} 

how.HLfit <- function(object, devel=FALSE, verbose=TRUE, format=print, ...) {
  info <- object$how # has all info, incl. e.g. fnname. Most of the code below is for back compt
  if (is.null(info)) {
    info <- list(MME_method=setdiff(object$MME_method,c("list")),
                 fit_time=object$fit_time,
                 "spaMM.version"=object$spaMM.version)
  }
  if (verbose) {
    fnname <- .get_bare_fnname.HLfit(object) # checks first object$how !
    if ( ! length(fnname)) {
      if (packageVersion("spaMM")==info[["spaMM.version"]]) {
        message("(Fitting function could not be identified)") # and it's not clear why
      } else { # object fitted with different version of spaMM, identical() may not work
        # function could not be identified, but we know why
      }
      mess <- "Model fitted by spaMM"
    } else mess <- paste0("Model fitted by spaMM::", paste(fnname))
    pretty_method <- .prettify_method(info$MME_method, by_y_augm=identical(info$switches["augZXy_cond"][[1]], TRUE))
    mess <- paste0(mess,", version ",info[["spaMM.version"]],
                   ", in ",info$fit_time,"s using ",paste(na.omit(pretty_method),collapse=","))
    if (object$models[["eta"]]=="etaHGLM") {
      if (identical(info[["obsInfo"]], TRUE)) {
        mess <- paste0(mess," (with obs. info. matrix)")
      } else if (identical(info[["obsInfo"]], FALSE)) {
        mess <- paste0(mess," (with exp. info. matrix)")
      } # else obsInfo may be OL, which results in no parenthetical detail here (canonical link GLM)
    } # otherwise for GLMs, the objective function does not depend on the use of obs vs exp weights . 
    mess <- paste0(mess,".")
    format(mess)  
  }
  if  (!  is.null(resid_fits <- object$resid_fits)) {
    if (verbose) {
      for (mv_it in seq_along(resid_fits)) if ( ! is.null(resid_fits[[mv_it]])) format(
        paste0("Residual model for sub-model ",mv_it," fitted using method: ",
               paste(resid_fits[[mv_it]]$how$MME_method, collapse=","),"."))
    }
    info$resid_info <- lapply(resid_fits,`[[`, x="how")
  } else if  (!  is.null(resid_fit <- object$resid_fit)) {
    resid_info <- object$resid_fit$how
    if (is.null(resid_info)) { # back compat for old spaMM objects
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

how.HLfitlist <- function(object, devel=FALSE, verbose=TRUE, format=print, ...) {
  info <- attr(object,"how")
  if (verbose) {
    fnname <- info$fnname
    if ( ! length(fnname)) {
      format(paste0("Model fitted by spaMM, version ",info[["spaMM.version"]],
                    ", in ",info$fit_time,"s."))
      if (packageVersion("spaMM")==info[["spaMM.version"]]) {
        message("(Fitting function could not be identified)") # and it's not clear why
      } else { # object fitted with different version of spaMM, identical() may not work
        # function could not be identified, but we know why
      }
    } else format(paste0("Model fitted by spaMM::", paste(fnname),", version ",info[["spaMM.version"]],
                         ", in ",info$fit_time,"s."))
  }
  hows <- lapply(object, how, devel=devel, verbose=FALSE, format=format, ...)
  #info <- t(structure(unlist(hows, recursive=FALSE),dim=c(length(hows[[1]]),length(hows)))) # assuming its square...
  invisible(c(info, list(hows=hows)))
}

response <- function(object, ...) object$y[,1L]

family.HLfit <- function(object, ...) {
  family <- object$family
  if (is.null(family)) {
    if ( ! is.null(object$families)) stop(paste("post-fit functions looking for 'family' in a multivariate-response fit should fail.\n", 
                                                "Contact the package maintainer for these post-fit functions if you wish extended functionality"))
  }
  family
}

.model.frame <- function(formula., object, 
                         is_framed= ! is.null(attr(object$data, "terms")), ...) {
  if (is_framed) return(object$data)
  # ELSE
  frame.form <- .subbarsMM(formula.) ## this comes from lme4 and converts (...|...) terms to some "+" form 
  environment(frame.form) <- environment(formula.)
  mf_call <- call("model.frame", data=object$data, formula=frame.form, drop.unused.levels=TRUE, 
                  weights=getCall(object)$prior.weights) # language object reused later
  eval(mf_call) ## data.frame with all vars required for fixef and for ranef, and a "terms" attribute
}

# (1) multcomp::glht() calls model.frame() so this is not purely internal and thus arguments such as with_resid must have defaults
# (2) https://developer.r-project.org/model-fitting-functions.html provide some thoughts
model.frame.HLfit <- function(formula, # the generic was poorly conceived... formula is not a formula
                              ...) {
  object <- formula # .../... so we rename for clarity
  formula. <- formula(object) 
  if (inherits(formula.,"list")) { 
    # in mv case the merged case should still have a single data frame without attributes, even if submodels's processed stored the model frames
    resu <- vector("list", length(formula.))
    for (mv_it in seq_along(formula.)) {
      form <- formula.[[mv_it]] 
      resu[[mv_it]] <- .model.frame(formula.=form, object, is_framed=FALSE, ...)
    }
    return(resu)
  }
  .model.frame(formula., object, ...)
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
#   .. ..- attr(*, "predvars.fixed")= language list(observedResponse, Environment1)   # <==============================
#   .. ..- attr(*, "varnames.fixed")= chr [1:2] "observedResponse" "Environment1"     # <==============================
#   .. ..- attr(*, "predvars.random")= language list(observedResponse, group)         # <==============================
#   - attr(*, "formula")=Class 'formula'  language observedResponse ~ Environment1 + (1 | group)     # <==============================
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

.get_old_info_uniqueGeo <- function(fitobject, char_rd=NULL) {
  if (fitobject$spaMM.version < "3.10.35") {
    info <- attr(fitobject,"info.uniqueGeo")
  } else info <- fitobject$ranef_info$info_oldUniqueGeo
  if ( ! is.null(char_rd)) {
    if (fitobject$spaMM.version > "2.3.18") { ## test TRUE for version > 2.3.18:
      info <- info[[char_rd]]
    } # else object was a single matrix 
  }
  info
}

df.residual.HLfit <- function(object, ...) {
  dfs <- object$dfs
  dfs$p_fixef_phi <- NULL # For dhglms the dfs of the residual model is more generally given by .calc_p_phi().... 
  length(object$y) - sum(.unlist(dfs)) # ... But the aim of the present fn is precisely to remove dfs of other parameters.
}

weights.HLfit <- function(object, type, ...) {
  if (type=="prior") {
    object$prior.weights
  } else object$envir$H_w.resid*residVar(object,which="phi")
}


