# The 'good' procedure using hlcorcall
.numInfo_objfn <- function(x, 
                           hlcorcall, # Conceived for a processed call (call from numInfo()); but .p4m_by_outer_optim() may use an unprocessed call.
                           # ie a call with a $processed argument
                           # which is not yet the case for .p4m_by_outer_optim() fits using a call to .p4m_by_iters().
                           # such a call may allow a $processed_call arg (not yet successfully used) rather than a $processed arg.
                           skeleton, 
                           transf, # signals transformed input. FALSE when called from numInfo(), 
                                   # TRUE when called from as_LMLT <- function(., transf=TRUE) as the argument is passed all the way down to here
                                   # If transf is true then 'x' is assumed to be in transf space AND the skeleton must match that.
                           objective,
                           moreargs,
                           full_beta,
                           objfn.extras=list()
                           ) {
  parlist <- relist(x, skeleton) # loses the keepInREML attribute, but this does not matter bc this attr is effective only in preprocessing..
  # Here if the original fit was an HLCor call there was no outer optim hence no moreargs computed and the 'moreargs' attr is NULL
  # But it is required only if there are hyper parameters, which are not supposed to be handled by HLCor
  # although they have to be handled by HLCor.obj)
  if (length(.unlist(parlist$trRanCoefs)) &&
      (length(objfn.extras[["user.lower"]]$ranCoefs) || length(objfn.extras[["user.upper"]]$ranCoefs))
  ) parlist  <- .apply_transformed_box_constr(parlist, skeleton=NULL, 
                                              user.lower=objfn.extras[["user.lower"]], 
                                              user.upper=objfn.extras[["user.upper"]], transf=transf)
  if ("etaFix" %in% names(parlist)) {
    # LUarglist contains any beta info (for outer beta optim) in $beta elements
    # this is converted to etaFix$beta by HLCor.obj() or HLfit.obj()
    # but .numInfo_objfn() expects a skeleton with $etaFix$beta....
    beta <- NULL
    if (transf && ! is.null(trBeta <- parlist$etaFix$trBeta)) {
    beta <- .spaMM.data$options$.betaInv(trBeta) 
    } else {
      beta <- parlist$etaFix$beta 
    }
    full_beta[names(beta)] <- beta
    hlcorcall$etaFix$beta <- full_beta
    parlist$etaFix <- NULL
  }

  # For PQL/L, if beta is in "which", the beta values are fixed to the fixef() when the grad of 
  # other, outer-estimated parameters is computed, and the resulting gradient is not zero!
  # Only when the beta are reestimated by h-lik does the gradient for the parameters is zero. 
  # So for PQL/L it might be better to check gradients with beta excluded 
  # (so they are not fixed to fixef() here).
  # But still, what will the full numInfo means?
  if ( ! is.null(hlcorcall[["processed"]])) {
    parlist <- .expand_hyper(parlist, hlcorcall$processed$hyper_info, moreargs=moreargs) ## input ranPars contains both unconstrained ranPars and $hyper
    parlist <- .canonizeRanPars(parlist, corr_info=hlcorcall$processed$corr_info, checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)
  } else {
    # transformed parameter in ~ user-level 'fixed' argument may be ignored.
    parlist <- .canonizeRanPars(parlist, corr_info=NULL, checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)
  }
  hlcorcall$fixed <- .modify_list(hlcorcall$fixed, parlist) 
  refit <- eval(hlcorcall)
  if (is.null(objective)) { # call from .get_covbeta() : there must not be an etaFix in the skeleton
    refit
  } else ( - refit$APHLs[[objective]])
}

.ad_hoc_jac_transf <- function(parlist, bdiag.=TRUE, moreargs) { # acobian: rows for vector-valued function, cols for elements of its argument
  gr <- parlist
  for (st in names(parlist)) {
    gr[[st]] <- switch(st,
                       "lambda"= {
                         dLam <- parlist[[st]]
                         for (it in seq_along(dLam)) {
                           if (dLam[it]<2e-4) {side <- 1} else side <- NA # ___F I X M E____ rethink
                           dLam[it] <- grad(.dispFn, dLam[it], side=side)
                         } 
                         if (bdiag. && length(dLam)>1L) dLam <- diag(x=dLam)
                         dLam
                       },
                       "ranCoefs"= {
                         drC <- parlist[[st]]
                         for (it in seq_along(drC)) drC[[it]] <- jacobian(.ranCoefsFn, drC[[it]], rC_transf=.spaMM.data$options$rC_transf) 
                         if (bdiag.) drC <- Matrix::bdiag(drC)
                         drC
                       },
                       "phi"= {
                         if (inherits(parlist[[st]],"list")) { # mv fit
                           dphilist <- parlist[[st]]
                           for (it in seq_along(dphilist)) dphilist[[it]] <- grad(.dispFn, dphilist[[it]])
                           if (bdiag.) dphilist <- diag(x=.unlist(dphilist))
                           dphilist
                         } else grad(.dispFn, parlist[[st]])
                       },
                       "rdisPars"= {
                         if (inherits(parlist[[st]],"list")) { # mv fit
                           dphilist <- parlist[[st]]
                           for (it in seq_along(dphilist)) dphilist[[it]] <- rep(1, length(dphilist[[it]]))
                           if (bdiag.) dphilist <- diag(x=.unlist(dphilist))
                           dphilist
                         } else rep(1, length(parlist[[st]]))
                       },
                       #
                       "etaFix" = {
                         dbeta <- grad(.spaMM.data$options$.betaFn, parlist[[st]][["beta"]])
                         if (bdiag. && length(dbeta)>1L) dbeta <- diag(x=dbeta)
                         dbeta
                        },
                       "nu" = grad(.nuFn, parlist[[st]], NUMAX=moreargs$NUMAX),
                       "rho" = grad(.rhoFn, parlist[[st]], RHOMAX=moreargs$RHOMAX),
                       "longdep" = grad(.longdepFn, parlist[[st]], LDMAX=moreargs$LDMAX),
                       "kappa" = grad(.kappaFn, parlist[[st]], KAPPAMAX=moreargs$KAPPAMAX),
                       "corrPars" = {
                         dcorrlist <- parlist[[st]]
                         char_rds <- names(dcorrlist)
                         for (it in seq_along(dcorrlist)) {
                           dcorrlist_it <- .ad_hoc_jac_transf(dcorrlist[[it]], moreargs=moreargs[[char_rds[[it]]]])
                           if (bdiag.  && length(dcorrlist_it)>1L) dcorrlist_it <- diag(x=.unlist(dcorrlist_it))
                           dcorrlist[[it]] <- dcorrlist_it
                         }
                         if (bdiag.) dcorrlist <- Matrix::bdiag(dcorrlist)
                         dcorrlist
                       }, # to identically process the "corrPars" sublist
                       # no transfo for other params =>
                       if (inherits(parlist[[st]],"list")) { 
                         dvec <- rep(1, length(.unlist(parlist[[st]])))
                         if ( ! bdiag.) {
                           dvec
                         } else relist(dvec, parlist[[st]])
                       } else rep(1, length(parlist[[st]]))
                       #                       stop(paste("Parameter transformation for",st,"not yet handled in .ad_hoc_grad()")) # COMP_nu would be most problematic to implement; but NB_sape and beta_prec easy ?
    )
  }
  if (bdiag.) gr <- Matrix::bdiag(gr) # "inefficient"
  gr
}

.ad_hoc_grXhessians_transf <- function(parlist, grad_list) {
  resu <- parlist
  for (st in names(parlist)) {
    hess <- switch(st,
                   "lambda"= {
                     ghlambdas <- parlist[[st]]
                     for (it in seq_along(ghlambdas)) ghlambdas[it] <- hessian(.dispFn, ghlambdas[it])
                     ghlambdas <- grad_list$trLambda * ghlambdas
                     if (length(ghlambdas)>1L) ghlambdas <- diag(x=ghlambdas)
                     ghlambdas
                   },
                   "ranCoefs"= {
                     rancoefs_list <- parlist[[st]]
                     for (it in seq_along(rancoefs_list)) {
                       ranCoefs_it <- rancoefs_list[[it]]
                       gend <- numDeriv::genD(.ranCoefsFn, ranCoefs_it, rC_transf=.spaMM.data$options$rC_transf)
                       hessians <- gend$D[,-seq_along(ranCoefs_it)] # each ith line represent the lower triangle of the hessian of rC_transf[i]
                       grXhessians <- grad_list$trRanCoefs[[it]] %*% hessians
                       grXhess <- diag(length(ranCoefs_it))
                       .lower.tri(grXhess,diag = TRUE) <- colSums(grXhessians)
                       .upper.tri(grXhess,diag = FALSE) <- grXhess[lower.tri(grXhess,diag = FALSE)]
                       rancoefs_list[[it]] <- grXhess
                     }
                     Matrix::bdiag(rancoefs_list)
                   },
                   "phi"= {
                     if (inherits(parlist[[st]],"list")) { # mv fit
                       d2philist <- parlist[[st]]
                       for (it in seq_along(d2philist)) d2philist[[it]] <- hessian(.dispFn, d2philist[[it]])
                       h <- .unlist(d2philist)*.unlist(grad_list$trPhi)
                       diag(x=h)
                     } else grad_list$trPhi*hessian(.dispFn, parlist[[st]])
                   },
                   "rdisPars"= {
                     if (inherits(parlist[[st]],"list")) { # mv fit
                       diag(0, nrow=length(.unlist(parlist[[st]])))
                     } else diag(0,nrow=length(parlist[[st]]))
                   },
                   #
                   "etaFix" = {
                     ghbetas <- parlist[[st]]$beta
                     for (it in seq_along(ghbetas)) ghbetas[it] <- hessian(.spaMM.data$options$.betaFn, ghbetas[it])
                     ghbetas <- grad_list$etaFix$trBeta * ghbetas
                     if (length(ghbetas)>1L) ghbetas <- diag(x=ghbetas)
                     ghbetas
                   },
                   #
                   # processing corrPars structured sublist
                   "nu" =  grad_list$trNu * hessian(.nuFn, parlist[[st]]),
                   "rho" = {
                     ghrhovec <- parlist[[st]]
                     for (it in seq_along(ghrhovec)) ghrhovec[it] <- hessian(.rhoFn, ghrhovec[it])
                     ghrhovec <- grad_list$trRho * ghrhovec
                     if (length(ghrhovec)>1L) ghrhovec <- diag(x=ghrhovec)
                     ghrhovec
                   },
                   "longdep" = grad_list$trLongDep * hessian(.longdepFn, parlist[[st]]),
                   "kappa" =  grad_list$trKappa * hessian(.kappaFn, parlist[[st]]),
                   "corrPars" = {
                     ghcorrlist <- parlist[[st]]
                     for (it in seq_along(ghcorrlist)) {
                       ghcorrlist_it <- .ad_hoc_grXhessians_transf(ghcorrlist[[it]], grad_list=grad_list$corrPars[it])
                       ghcorrlist[[it]] <- .unlist(ghcorrlist_it)
                     }
                     ghcorrlist <- .unlist(ghcorrlist)
                     if (length(ghcorrlist)>1L) ghcorrlist <- diag(x=ghcorrlist)     
                     ghcorrlist
                   }, 
                   # no transfo for other params =>
                   if (inherits(parlist[[st]],"list")) {
                     diag(0, nrow=length(.unlist(parlist[[st]])))
                   } else diag(0,nrow=length(parlist[[st]]))
                   #                       stop(paste("Parameter transformation for",st,"not yet handled in .ad_hoc_grad()")) # COMP_nu would be most problematic to implement; but NB_sape and beta_prec easy ?
    )
    resu[[st]] <- hess 
  }
  Matrix::bdiag(resu) # "inefficient"
}

## (yet partial) inverse of .canonizeRanpars. This work bc the latter fn do not assume all parameters are transformed 
.ad_hoc_trRanpars <- function(ranPars,
                              # moreargs, # might be needed later
                              rC_transf=.spaMM.data$options$rC_transf) {
  trRanpars <- list()
  for (st in names(ranPars)) { # MUST keep parameter order
    switch(st, # here working through its side effects, not its return value
           "lambda" = {trRanpars$trLambda <-.dispFn(ranPars$lambda)},
           "phi" = {
             if (is.list(phi <- ranPars$phi)) {
             trRanpars$trPhi <- lapply(phi, .dispFn)
             } else trRanpars$trPhi <-.dispFn(phi)
           },
           "ranCoefs" = {trRanpars$trRanCoefs <- lapply(ranPars$ranCoefs, .ranCoefsFn, rC_transf=rC_transf)},
           {trRanpars[st] <- ranPars[st]} 
    )
  }
  trRanpars
}


.post_process_hlcorcall <- function(hlcorcall, 
                                    ranpars, # beware canonical/non canonical in later extensions of this fn.
                                    # optional for beta numDerivs:
                                    beta_eta=NULL, fitobject,
                                    ori_off=model.offset.HLfit(fitobject)) {
  processed <- hlcorcall$processed
  if (is.list(processed)) {
    proc1 <- processed[[1L]]
  } else proc1 <- processed
  .assignWrapper(processed, paste0("return_only <- \"",proc1$objective,"APHLs\""))
  # I must clean the preprocessed info for fixed ranCoefs...
  if (! is.null(ranpars$trRanCoefs)) {
    rancoefs <- .ranCoefsInv(ranpars$trRanCoefs, rC_transf= .spaMM.data$options$rC_transf)
  } else rancoefs <- ranpars$ranCoefs
  if (! is.null(rancoefs)) {
    for (char_rd in names(rancoefs)) {
      rd <- as.numeric(char_rd)
      processed$ranCoefs_blob$is_set[rd] <- FALSE
      processed$ranCoefs_blob$LMatrices[rd] <- list(NULL)
    }
  }
  if (length(beta_eta)) { # same as for outer optim of beta in fitme()/fitmv()
    ## cf explanations on code simular to this block in preprocessing functions. 
    # Here for numInfo computation we mix features of inner and outer optim. Not the most lucid block of code...
    # A comment in HLfit_body() says "AUGI0_ZX$X.pv must correspondingly have been reduced by .preprocess()"
    # Indeed. There was an etaFix$beta in preprocessing, which allowed as call to 
    # .process_betaFix() before merged_X was put into AUGI0_ZX.
    betanames <- names(fixef(fitobject))
    X.pv <- model.matrix(fitobject)
    X_off <-.subcol_wAttr(X.pv, j=betanames, drop=FALSE)
    X.pv <- .subcol_wAttr(X.pv, j=setdiff(colnames(X.pv),betanames), drop=FALSE)
    #     ori_off <- model.offset.HLfit(fitobject) # processed$off differs from it as get_HLCorcall(., etaFix) -> .preprocess(...etaFix) adds the etaFix-derived offset.
    processed$X_off_fn <- .def_off_fn(X_off, ori_off=ori_off) 
    processed[["vecdisneeded_ori"]] <-  processed[["vecdisneeded"]]
    processed[["vecdisneeded"]] <- processed[["vecdisneeded"]] & ncol(X.pv) 
  }
  info <- list(objective=proc1$objective)
  info # not the primary effect of the fn.
}

.calc_grad_thr <- function(skeleton, fitobject, beta_eta=skeleton$etaFix$beta) {
  thr <- relist(rep(0.1, length(unlist(skeleton, recursive = TRUE, use.names = FALSE))), skeleton)
  if (fitobject$how$spaMM.version>"4.1.58") { # => made scale_info available afterwards 
    if (length(beta_eta)) thr$etaFix$beta <- 0.1* attr(model.matrix(fitobject),"scale_info")
    if (length(rdisPars <- skeleton$rdisPars)) {
      if (is.list(rdisPars)) {
        for (mv_it in names(rdisPars)) {
          thr$rdisPars[[mv_it]] <- 0.1*attr(fitobject$families[[mv_it]]$resid.model$X,"scale_info")
        }
      } else thr$rdisPars <- 0.1*attr(fitobject$family$resid.model$X,"scale_info")
    }
  }
  thr # structured list...
}

.noisy_genD <- function(func, x, method.args = list(), 
                        ...) {
  args <- list(eps = 1e-04, d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07), 
               r = 4, v = 2)
  args[names(method.args)] <- method.args
  d <- args$d
  r <- args$r
  v <- args$v
  if (v != 2) 
    stop("The current code assumes v is 2 (the default).")
  f0 <- func(x, ...)
  n <- length(x)
  h0 <- abs(d * x) + args$eps * (abs(x) < args$zero.tol)
  D <- matrix(0, length(f0), (n * (n + 3L))/2L)
  Daprox <- matrix(0, length(f0), r)
  Hdiag <- matrix(0, length(f0), n)
  Haprox <- matrix(0, length(f0), r)
  df <- list()
  it <- 1L
  df[[it]] <- c(x*0, f0, 0, 0)
  it <- it+1L
  for (i in 1:n) {
    h <- h0
    for (k in 1:r) {
      xx <- x + (i == (1:n)) * h
      f1 <- func(xx, ...)
      df[[it]] <- c(xx-x, f1, (h/h0)[i], (h/h0)[i])
      it <- it+1L
      xx <- x - (i == (1:n)) * h
      f2 <- func(xx, ...)
      df[[it]] <- c(xx-x, f2, (h/h0)[i], (h/h0)[i])
      it <- it+1L
      Daprox[, k] <- (f1 - f2)/(2 * h[i])
      Haprox[, k] <- (f1 - 2 * f0 + f2)/h[i]^2
      h <- h/v
    }
    for (m in 1:(r - 1)) for (k in 1:(r - m)) {
      Daprox[, k] <- (Daprox[, k + 1] * (4^m) - Daprox[, k])/(4^m - 1)
      Haprox[, k] <- (Haprox[, k + 1] * (4^m) - Haprox[, k])/(4^m - 1)
    }
    D[, i] <- Daprox[, 1]
    Hdiag[, i] <- Haprox[, 1]
  }
  u <- n
  for (i in 1:n) {
    for (j in 1:i) {
      u <- u + 1L
      if (i == j) 
        D[, u] <- Hdiag[, i]
      else {
        h <- h0
        for (k in 1:r) {
          xx <- x + (i == (1:n)) * h + (j == (1:n)) * h
          f1 <- func(xx, ...)
          df[[it]] <- c(xx-x, f1, (h/h0)[i], (h/h0)[j])
          it <- it+1L
          xx <- x - (i == (1:n)) * h - (j == (1:n)) * h
          f2 <- func(xx, ...)
          df[[it]] <- c(xx-x, f2, (h/h0)[i], (h/h0)[j])
          it <- it+1L
          Daprox[, k] <- (f1 - 2 * f0 + f2 - Hdiag[, i] * h[i]^2 - Hdiag[, j] * h[j]^2)/(2 * h[i] * h[j])
          h <- h/v
        }
        for (m in 1:(r - 1)) for (k in 1:(r - m)) Daprox[, k] <- 
            (Daprox[, k + 1] * (4^m) - Daprox[,k])/(4^m - 1)
        D[, u] <- Daprox[, 1]
      }
    }
  }
  df <- as.data.frame(do.call(rbind,df))
  colnames(df) <- c(paste0("x",seq_along(x)),"y","i","j")
  D <- list(D = D, p = length(x), f0 = f0, func = func, x = x, 
            d = d, method = "Richardson", method.args = args, df=df)
  class(D) <- "Darray"
  invisible(D)
}

# wrapper for .noisy_genD()
.spaMM_hessian <- function (func, x, method.args = list(), ...) {
  args <- list(eps = 1e-04, d = 0.1, zero.tol = sqrt(.Machine$double.eps/7e-07), 
               r = 4, v = 2, show.details = FALSE)
  args[names(method.args)] <- method.args
  
  DLM <- .noisy_genD(func, x, method.args = args, ...)
  
  coefs <- DLM$D[-seq(length(x))]
  H <- matrix(nrow=length(x), ncol=length(x))
  H[upper.tri(H,diag=TRUE)] <- coefs
  H <- t(H)
  H[upper.tri(H,diag=TRUE)] <- coefs
  
  attr(H,"df") <- DLM$df
  H
}

.post_process_hessian <- function(resu, parnames, removand_rC, sing, verbose) {
  dimnames(resu) <- list(parnames, parnames)
  
  # remove partially-fixed ranCoefs and those identified by check deriv
  if ( length(removand_rC)) {
    parnames <- setdiff(parnames, removand_rC) 
    resu <- resu[parnames, parnames, drop=FALSE]
  }
  
  if (sing) {
    ev <- eigen(resu, only.values = TRUE)$values
    if (any(ev < sing)) {
      attr(ev,"sing") <- sing
      class(ev) <- c(class(ev),"singeigs")
      if (verbose) message("Information matrix has suspiciously small eigenvalues.")
    }
    attr(resu,"eigvals") <- ev
  }
  resu
}

.smooth_hessian <- function(df, x) {
  form <- as.formula(
    paste0("y ~ poly(cbind(", 
           paste0("x",seq_along(x), collapse=","),
           "),degree=2,raw=TRUE)")
  )
  lmfit <- fitme(form, data=df, resid.model = ~ i*j)
  coefs <- fixef(lmfit)[-c(1L, 1L+cumsum(seq_along(x)))]
  hess <- matrix(0,nrow=length(x), ncol=length(x))
  mask <- upper.tri(hess,diag=TRUE)
  hess[mask] <- coefs
  hess <- hess+ t(hess) # doubling on the diag: OK here
  hess
}



# ___F I X M E____ I should enable numInfo() on resid models, [but then see comment on FIXME in ..calcPHI()]
# working on GLMs ($phi.object$glm_phi) and HGLMs ($resid_fit)
# while not a joint numInfo it would not be worse than other implementations.
# ___F I X M E____ store the result in the fitobject along which 'which' and perhaps a few other attrs, 
# so that it can be reused conditionnally on such attrs...
numInfo <- function(fitobject, 
                    transf=FALSE, 
                    which=NULL,
                    check_deriv=TRUE,       # outer estim without refit = no leverages checks... => it seems better to always check them
                    sing=1e-05,
                    verbose=FALSE,
                    refit_hacks=list(),
                    attrs=NULL, # "df", "smoothed"
                    return.="",
                    # method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE),
                    ...) {
  ## We need an X_off_fn so that the etaFix is used to build an offset. 
  ## IRLS function do not really handle etaFix. We need an etaFix at preprocessing stage so that columns are suppressed from AUGI0_ZX$X.pv
  ## => => get_HLCorcall(fitobject, ... etaFix=list(beta=fixef(fitobject)))
  ## Currently X_off_fn is set by .preprocess() only given an init beta, not a given beta, => we set up it in this function
  
  ### REML: the resulting SEs are consistent with those from the beta table... (with keepInREML used in numInfo())
  is_REML <- .REMLmess(fitobject,return_message=FALSE)
  is_PQL_sl <- fitobject$HL[1L]==0L # fixed effects estimated by h-lik
  if (is.null(which)) {
    if (is_REML || is_PQL_sl) {
      which <-      c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "COMP_nu", "beta_prec", "rdisPars")
    } else which <- c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "COMP_nu", "beta_prec", "rdisPars", 
                      "beta")
  } else {
    if (is_PQL_sl &&
        "beta" %in% which) warning("'beta'in 'which' argument may give confusing results when PQL approx. has been used.",
                                   immediate. = TRUE)
  }
  # if (is_REML && "beta" %in% which) {
  #   REMLformula <- formula(fitobject)
  # } else REMLformula <- NULL

  # but the skeleton understood by hlcorcall is already a re-merging of fixed and oprimized values, 
  # so it also needs a full ranCoefs, the hessian must first be computed on full ranCoefs, 
  # and fixed columns be removed afterwards 

  where_rP <- .get_fittedPars(fitobject, fixef=FALSE, verbose=verbose, partial_rC="keep", phifits=FALSE, phiPars=FALSE) 
  for_which_rP <- .get_fittedPars(fitobject, which=which, fixef=FALSE, verbose=verbose, partial_rC="keep", phifits=FALSE, phiPars=FALSE) # may be zero-length
  if ("beta" %in% which) {beta_eta <- structure(fixef(fitobject), keepInREML=TRUE)} else beta_eta <- NULL
  outer_call <- getCall(fitobject)
  if ( ! is.null(refit_verbose <- refit_hacks$verbose)) {
    refit_verbose <- .modify_list(outer_call$verbose, refit_verbose)
  } else refit_verbose <- outer_call$verbose
  if (.get_bare_fnname.HLfit(fitobject, call.=outer_call)=="pois4mlogit") {
    outer_call$data <- fitobject$data # with "good" .dynoffset
    outer_call$init <- get_inits_from_fit(fitobject)$init # !! get inits from the fit, not the call
    outer_call["initfn"] <- NULL 
    hlcorcall <- outer_call
    objective <- .get_objective(fitobject)
    proc_info <- list(objective=objective) 
  } else {
    hlcorcall <- get_HLCorcall(fitobject, 
                               fixed=where_rP, # If one use 'for_which_rP' here, parameters still get fixed cf (hlcorcall$fixed),
                               # but to what is usually  default initial values (of not interest here)
                               # This is a "feature", not a design decision....
                               # fitmv() does not handle REMLformula *outside* the submodels, so either fitmv should be modified, or
                               # the correct REMLformula must be automatically set as side effect of the keepInREML attribute,
                               # or (chosen solution) the required effect must be achieved by keepInREML => cf .preprocess_X_XRe_off(), and related code for fitmv()
                               #                           REMLformula=REMLformula,
                               etaFix=list(beta=beta_eta),  # => the etaFix argument means that processed$off is modified by it
                               verbose=refit_verbose,
                               init=NULL # not the call's$init (and get_inits_from_fit(fitobject)$init may conflict with fixed values)
    ) 
    proc_info <- .post_process_hlcorcall(hlcorcall, ranpars=for_which_rP, beta_eta=beta_eta, fitobject=fitobject)
    if (is.null(check_deriv)) check_deriv <- (
      length(for_which_rP$lambda) && any(for_which_rP$lambda<1e-6) ||
        ( ! is.null(fitobject$warnings$allLeveLam1))
    ) 
  }
  #
  skeleton <- for_which_rP
  if (length(beta_eta)) skeleton$etaFix$beta <- beta_eta
  if (verbose) print(skeleton)
  
  # cannot remove partially-fixed ranCoefs too early. They must be kept in skeleton in all cases
  ufixed <- na.omit(unlist(outer_call$fixed))
  fixednames <- names(ufixed) # .preprocess_fixed() removed any fancy names given by users.
  uskeleton <- unlist(skeleton)
  removand_rC <- intersect(fixednames,names(uskeleton))
  
  if (return.=="grad") {
    side <- .calc_grad_side_arg(skeleton)
    uside <- unlist(side)
    gr_neg_APHL <- grad(func = .numInfo_objfn, x =unlist(skeleton), side=uside, skeleton=skeleton, 
                        hlcorcall=hlcorcall, transf=transf, full_beta=fixef(fitobject), 
                        objective=proc_info$objective, moreargs=.get_moreargs(fitobject), ...)
    names(gr_neg_APHL) <- c(names(unlist(skeleton[setdiff(names(skeleton), "etaFix")])), 
                            names(skeleton$etaFix$beta))
    return( - gr_neg_APHL)
  }
  
  if (check_deriv) {
    thr <- .calc_grad_thr(skeleton, fitobject, beta_eta) # effect visible on numInfo(fitme(cases~1+(1|id),family=negbin1(), data=scotlip, resid.model = ~ population))
    removand <- .check_numDeriv_task(skeleton, .numInfo_objfn, hlcorcall, transf, proc_info, 
                                     moreargs=.get_moreargs(fitobject), full_beta=fixef(fitobject), 
                                     thr=unlist(thr), ...) # assumes untransformed param
    tmp <- uskeleton
    tmp[removand] <- NaN
    # (1) check on finally retained parameters
    if ( ! length(unlist(.rmNaN(relist(tmp,skeleton))))) {
      warning("No fitted (co-)variance parameters whose information matrix could be evaluated.", immediate. = TRUE)
      return(NULL)
    }
    # (2) but put back rC parameters than cannot yet be removed:
    # Actually we cannot yet remove other elements of ranCoefs:
    # the only ambiguous case would be the case where .check_numDeriv_task() and partial_rC together flag
    # a full ranCoef for removal. In that case a more efficient numInfo computation would be possible by removing the ranCoef.
    # That does not seem worth the effort.
    check_rC <- intersect(names(tmp[which(removand)]), names(unlist(skeleton["ranCoefs"])))
    removand_rC <- unique(c(removand_rC, check_rC))
    tmp[removand_rC] <- uskeleton[removand_rC] 
    tmp <- relist(tmp,skeleton)
    skeleton <- .rmNaN(tmp)
  } else if ( ! length(skeleton)) {
    warning("No fitted (co-)variance parameters whose information matrix could be evaluated.", immediate. = TRUE)
    return(NULL)
  } 
  
  # Final value always refer to untransformed params:
  parnames <- c(names(unlist(skeleton[setdiff(names(skeleton), "etaFix")])), names(skeleton$etaFix$beta))
  if (transf) { # from CANON to TRANSF
    canon_skeleton <- skeleton
    skeleton <- .ad_hoc_trRanpars(skeleton)
    if (length(beta_eta)) {
      skeleton$etaFix$trBeta <- .spaMM.data$options$.betaFn(beta_eta)
      skeleton$etaFix$beta <- NULL
    }
    #  and the .numInfo_objfn must back-transform parameters
  }
  #
  
  resu <- .spaMM_hessian(func = .numInfo_objfn, x = unlist(skeleton), 
                  skeleton=skeleton, hlcorcall=hlcorcall, 
                  transf=transf, # whether x is on untransformed scale
                  objective=proc_info$objective, full_beta=fixef(fitobject), 
                  moreargs=.get_moreargs(fitobject), ...)
  
  if ("smoothed" %in% attrs) {
    df <- attr(resu,"df")
    if ( ! is.null(df)) {
      smoothed <- .smooth_hessian(df, x=parnames)
      smoothed <- .post_process_hessian(resu=smoothed, parnames, removand_rC, sing, verbose)
      attr(resu,"smoothed") <- smoothed
      attr(resu,"df") <- NULL
    }
  }
  if ( ! ("df" %in% attrs)) attr(resu,"df") <- NULL

  # .assignWrapper(processed, paste0("return_only <- NULL")) # Not necess bc $processed is created by the local get_HLCorcall() and freed when this fn is exited. 
  if (transf) {
    jacTransf <- .ad_hoc_jac_transf(canon_skeleton, moreargs=.get_moreargs(fitobject))
    resu <- as.matrix(crossprod(jacTransf, resu %*% jacTransf))
    if (is_REML) { # grad_obj should be ~ 0 , for ML fits at least, so this computation should not be necessary
      grad_obj <- grad(func = .numInfo_objfn, x = unlist(skeleton), # on transformed scale
                       skeleton=skeleton, hlcorcall=hlcorcall, 
                       transf=transf, full_beta=fixef(fitobject), 
                       objective=proc_info$objective, 
                       moreargs=.get_moreargs(fitobject), ...)
      grad_list <- relist(grad_obj, skeleton) 
      gr_obj_X_hess_transf <- .ad_hoc_grXhessians_transf(canon_skeleton, grad_list) # 2nd deriv of each parameter transform.
      resu <- resu+gr_obj_X_hess_transf
    }
  }
  resu <- .post_process_hessian(resu, parnames, removand_rC, sing, verbose)
  resu
}


print.singeigs <- function(x, ...) {
  sing <- attr(x,"sing")
  xx <- sapply(x, function(v) {
    ifelse(v<sing, cli::style_underline(signif(v)), signif(v))
  })
  cat(paste(xx))
}

