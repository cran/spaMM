# The 'good' procedure using hlcorcall
.numInfo_objfn <- function(x, hlcorcall, skeleton, 
                           transf, # signals transformed input. FALSE when called from numInfo(), 
                                   # TRUE when called from as_LMLT <- function(., transf=TRUE) as the argument is passed all the way down to here
                           objective,
                           moreargs
                           ) {
  # browser()
  #print(x)
  parlist <- relist(x, skeleton)
  # Here if the original fit was an HLCor call there was no outer optim hence no moreargs computed and the 'moreargs' attr is NULL
  # But it is required only if there are hyper parameters, which are not supposed to be handled by HLCor
  # although they have to be handled by HLCor.obj)
  parlist <- .expand_hyper(parlist, hlcorcall$processed$hyper_info, moreargs=moreargs) ## input ranPars contains both unconstrained ranPars and $hyper
  if (transf) {
    if ( ! is.null(trBeta <- parlist$etaFix$trBeta)) {
      hlcorcall$etaFix$beta <- .betaInv(trBeta) 
      parlist$etaFix <- NULL
    }
    # print(hlcorcall$etaFix$beta)
    parlist <- .canonizeRanPars(parlist, corr_info=hlcorcall$processed$corr_info, checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)
  } else {
      hlcorcall$etaFix <- parlist$etaFix 
      parlist$etaFix <- NULL
  }
  hlcorcall$fixed <- .modify_list(hlcorcall$fixed, parlist) 
# str(hlcorcall$fixed)
#str(hlcorcall$etaFix)
  refit <- eval(hlcorcall)
#print(refit)
  if (is.null(objective)) {
    refit
  } else ( - refit$APHLs[[objective]])
}


.ad_hoc_grad_transf <- function(parlist, unlist.=TRUE) {
  gr <- parlist
  for (st in names(parlist)) {
    gr[[st]] <- switch(st,
                       "lambda"= grad(.dispFn, parlist[[st]]), 
                       "ranCoefs"= {
                         drC <- parlist[[st]]
                         for (it in seq_along(drC)) drC[[it]] <- grad(.dispFn, drC[[it]])
                         drC
                       },
                       "phi"= {
                         if (inherits(parlist[[st]],"list")) { # mv fit
                           dphilist <- parlist[[st]]
                           for (it in seq_along(dphilist)) dphilist[[it]] <- grad(.dispFn, dphilist[[it]])
                           dphilist
                         } else grad(.dispFn, parlist[[st]])
                       },
                       #
                       "beta" = grad(.betaFn, parlist[[st]]),
                       "etaFix" = .ad_hoc_grad_transf(parlist[[st]]), # to identically process the "etaFix" sublist (its "beta" element)
                       #
                       "nu" = grad(.nuFn, parlist[[st]]),
                       "rho" = grad(.rhoFn, parlist[[st]]),
                       "longdep" = grad(.longdepFn, parlist[[st]]),
                       "kappa" = grad(.kappaFn, parlist[[st]]),
                       "corrPars" = {
                         dcorrlist <- parlist[[st]]
                         for (it in seq_along(dcorrlist)) dcorrlist[[it]] <- .ad_hoc_grad_transf(.dispFn, dcorrlist[[it]])
                         dcorrlist
                       }, # to identically process the "corrPars" sublist
                       stop("element not yet handled in .ad_hoc_grad()")
    )
  }
  if (unlist.) gr <- unlist(gr, recursive=TRUE, use.names = FALSE)
  gr
}

.ad_hoc_diag_hess_transf <- function(parlist, unlist.=TRUE) {
  gr <- parlist
    for (st in names(parlist)) {
      hess <- switch(st,
                     "lambda"= {
                       d2lambdas <- parlist[[st]]
                       for (it in seq_along(d2lambdas)) d2lambdas[it] <- hessian(.dispFn, d2lambdas[it])
                       d2lambdas
                     },
                     "ranCoefs"= {
                       rancoefs_list <- parlist[[st]]
                       for (it in seq_along(rancoefs_list)) rancoefs_list[[it]] <- hessian(.ranCoefsFn, rancoefs_list[[it]])
                       rancoefs_list
                     },
                     "phi"= {
                       if (inherits(parlist[[st]],"list")) { # mv fit
                         d2philist <- parlist[[st]]
                         for (it in seq_along(d2philist)) d2philist[[it]] <- hessian(.dispFn, d2philist[[it]])
                         d2philist
                       } else hessian(.dispFn, parlist[[st]])
                     },
                     #
                     # processing etaFix structured sublist 
                     "beta" = {
                       d2betas <- parlist[[st]]
                       for (it in seq_along(d2betas)) d2betas[it] <- hessian(.betaFn, d2betas[it])
                       d2betas
                     },
                     "etaFix" = .ad_hoc_diag_hess_transf(parlist[[st]]),
                     #
                     # processing corrPars structured sublist
                     "nu" = hessian(.nuFn, parlist[[st]]),
                     "rho" = {
                       d2rhovec <- parlist[[st]]
                       for (it in seq_along(d2rhovec)) d2rhovec[it] <- hessian(.rhoFn, d2rhovec[it])
                       d2rhovec
                     },
                     "longdep" = hessian(.longdepFn, parlist[[st]]),
                     "kappa" = hessian(.kappaFn, parlist[[st]]),
                     "corrPars" = {
                       d2corrlist <- parlist[[st]]
                       for (it in seq_along(d2corrlist)) d2corrlist[[it]] <- .ad_hoc_diag_hess_transf(.dispFn, d2corrlist[[it]])
                       d2corrlist
                     }, 
                     stop("element not yet handled in .ad_hoc_diag_hess()")
      )
      gr[[st]] <- hess 
    }
  if (unlist.) gr <- unlist(gr, recursive=TRUE, use.names = FALSE)
  gr
}

## (yet partial) inverse of .canonizeRanpars 
.ad_hoc_trRanpars <- function(ranPars,
                              moreargs=NULL, rC_transf=.spaMM.data$options$rC_transf) {
  trRanpars <- ranPars
  if ( ! is.null(ranPars$lambda)) {
    trRanpars$trLambda <-.dispFn(ranPars$lambda)
    trRanpars$lambda <- NULL
  }
  if ( ! is.null(phi <- ranPars$phi)) {
    if (is.list(phi)) {
      trRanpars$trPhi <- lapply(phi, .dispFn)
    } else trRanpars$trPhi <-.dispFn(ranPars$phi)
    trRanpars$phi <- NULL
  }
  if ( ! is.null(ranPars$ranCoefs)) {
    trRanpars$trRanCoefs <- lapply(ranPars$ranCoefs, .ranCoefsFn, rC_transf=rC_transf)
    trRanpars$ranCoefs <- NULL
  }
  trRanpars
}


.post_process_hlcorcall <- function(hlcorcall, ranpars,
                                    # optional for beta numDerivs
                                    beta_eta=NULL, fitobject) {
  processed <- hlcorcall$processed
  if (is.list(processed)) {
    proc1 <- processed[[1L]]
  } else proc1 <- processed
  .assignWrapper(processed, paste0("return_only <- \"",proc1$objective,"APHLs\""))
  # I must clean the preprocessed info for fixed ranCoefs...
  if (! is.null(ranpars$ranCoefs)) {
    for (char_rd in names(ranpars$ranCoefs)) {
      rd <- as.numeric(char_rd)
      processed$ranCoefs_blob$is_set[rd] <- FALSE
      processed$ranCoefs_blob$LMatrices[rd] <- list(NULL)
    }
  }
  if (length(beta_eta)) {
    betanames <- names(fixef(fitobject))
    X.pv <- fitobject$X.pv
    X_off <-.subcol_wAttr(X.pv, j=betanames, drop=FALSE)
    X.pv <- .subcol_wAttr(X.pv, j=setdiff(colnames(X.pv),betanames), drop=FALSE)
    ori_off <- model.offset.HLfit(fitobject) # processed$off differs from it as get_HLCorcall(., etaFix) -> .preprocess(...etaFix) adds the etaFix-derived offset.
    processed$X_off_fn <- .def_off_fn(X_off, ori_off=ori_off)
  }
  info <- list(objective=proc1$objective)
  info
}

numInfo <- function(fitobject, 
                    transf=FALSE, 
                    which=c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "beta_prec", "beta"),
                    verbose=FALSE,
                    ...) {
  ## We need an X_off_fn so that the etaFix is used to build an offset. 
  ## IRLS function do not really handle etaFix. We need an etaFix at preprocessing stage so that columns are suppressed from AUGI0_ZX$X.pv
  ## => => get_HLCorcall(fitobject, ... etaFix=list(beta=fixef(fitobject)))
  ## Currently X_off_fn is set by .preprocess() only given an init beta, not a given beta, => we set up it in this function
  
  is_REML <- .REMLmess(fitobject,return_message = FALSE) # used anyway, say gradient comput.
  ### REML: the resulting SEs are consistent with those from the beta table... (with keepInREML used in numInfo())
  # 
  # if (is_REML) { # perhaps add an 'objective' argument ?
  #   if ("beta" %in% which) {
  #     warning("numInfo() computation on REML fits:\n derivatives of log restricted likelihood wrt 'beta' coefficients may not be meaningful.", immediate.=TRUE)
  #     # I could implement computation of separate blocks for p_bv vs p_v-maximizing parameters. What about cross-blocks of second derivatives ?
  #   } else if (verbose) message("numInfo() computation on REML fits: derivatives are those of log restricted likelihood.", immediate.=TRUE)
  # }
  # 
  
  where_rP <- .get_fittedPars(fitobject, fixef=FALSE, verbose=verbose) 
  for_which_rP <- .get_fittedPars(fitobject, which=which, fixef=FALSE, verbose=verbose) 
  if ("beta" %in% which) {beta_eta <- structure(fixef(fitobject), keepInREML=TRUE)} else beta_eta <- NULL
  hlcorcall <- get_HLCorcall(fitobject, 
                             fixed=where_rP, # If one use 'for_which_rP' here, parameters still get fixed cf (hlcorcall$fixed),
                                             # but to what is usually  default initial values (of not interest here)
                                             # This is a "feature", not a design decision....
                             etaFix=list(beta=beta_eta)) 
  # => the etaFix argument means that processd$off is modified by it 
  proc_info <- .post_process_hlcorcall(hlcorcall, ranpars=for_which_rP, beta_eta=beta_eta, fitobject=fitobject)
  #
  skeleton <- for_which_rP
  if (length(beta_eta)) skeleton$etaFix$beta <- beta_eta
  if (verbose) print(skeleton)
  if (transf) {
    canon_skeleton <- skeleton
    skeleton <- .ad_hoc_trRanpars(for_which_rP)
    if (length(beta_eta)) {
      skeleton$etaFix$trBeta <- .betaFn(beta_eta)
    }
    #  and the .numInfo_objfn must back-transform parameters
  }
  moreargs <- .get_moreargs(fitobject)
  #
  resu <- hessian(func = .numInfo_objfn, x = unlist(skeleton), 
                  skeleton=skeleton, hlcorcall=hlcorcall, 
                  transf=transf, # whether x is on untransformed scale
                  objective=proc_info$objective, 
                  moreargs=moreargs, ...)
  # .assignWrapper(processed, paste0("return_only <- NULL")) # Not necess bc $processed is created by the local get_HLCorcall() and freed when this fn is exited. 
  if (transf) {
    grTransf <- .ad_hoc_grad_transf(canon_skeleton)
    grTransf <- diag(x=grTransf, ncol=length(grTransf))
    resu <- grTransf %*% resu %*% grTransf
    if (is_REML) { # grad_obj should be ~ 0 , for ML fits at least, so this computation should not be necessary
      grad_obj <- grad(func = .numInfo_objfn, x = unlist(skeleton), # on transformed scale
                       skeleton=skeleton, hlcorcall=hlcorcall, 
                       transf=transf, 
                       objective=proc_info$objective, 
                       moreargs=moreargs, ...)
      diag_hess_transf <- .ad_hoc_diag_hess_transf(canon_skeleton)
      resu <- resu+diag(x=diag_hess_transf*grad_obj)
    }
  }
  parnames <- c(names(unlist(for_which_rP)), names(beta_eta))
  dimnames(resu) <- list(parnames, parnames)
  resu
}
