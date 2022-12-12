# The 'good' procedure using hlcorcall
.numInfo_objfn <- function(x, hlcorcall, skeleton, 
                           transf, # signals transformed input. FALSE when called from numInfo(), 
                                   # TRUE when called from as_LMLT <- function(., transf=TRUE) as the argument is passed all the way down to here
                           objective,
                           moreargs
                           ) {
  # browser()
  #print(x)
  parlist <- relist(x, skeleton) # loses the keepInREML attribute, but this does not matter bc this attr is effective only in preprocessing..
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
  if (is.null(objective)) { # call from .get_covbeta() : there must not be an etaFix in the skeleton
    refit
  } else ( - refit$APHLs[[objective]])
}

.ad_hoc_jac_transf <- function(parlist, bdiag.=TRUE) { # acobian: rows for vector-valued function, cols for elements of its argument
  gr <- parlist
  for (st in names(parlist)) {
    gr[[st]] <- switch(st,
                       "lambda"= {
                         dLam <- parlist[[st]]
                         for (it in seq_along(dLam)) {
                           if (dLam[it]<2e-4) {side <- 1} else side <- NA # ____F I X M E____ rethink
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
                       #
                       "etaFix" = {
                         dbeta <- grad(.betaFn, parlist[[st]][["beta"]])
                         if (bdiag. && length(dbeta)>1L) dbeta <- diag(x=dbeta)
                         dbeta
                        },
                       "nu" = grad(.nuFn, parlist[[st]]),
                       "rho" = grad(.rhoFn, parlist[[st]]),
                       "longdep" = grad(.longdepFn, parlist[[st]]),
                       "kappa" = grad(.kappaFn, parlist[[st]]),
                       "corrPars" = {
                         dcorrlist <- parlist[[st]]
                         for (it in seq_along(dcorrlist)) {
                           dcorrlist_it <- .ad_hoc_jac_transf(.dispFn, dcorrlist[[it]])
                           if (bdiag.  && length(dcorrlist_it)>1L) dcorrlist_it <- diag(x=.unlist(dcorrlist_it))
                           dcorrlist[[it]] <- dcorrlist_it
                         }
                         if (bdiag.) dcorrlist <- Matrix::bdiag(dcorrlist)
                         dcorrlist
                       }, # to identically process the "corrPars" sublist
                       stop("element not yet handled in .ad_hoc_grad()")
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
                   #
                   "etaFix" = {
                     ghbetas <- parlist[[st]]$beta
                     for (it in seq_along(ghbetas)) ghbetas[it] <- hessian(.betaFn, ghbetas[it])
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
                   stop("element not yet handled in .ad_hoc_diag_hess()")
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
                    which=NULL,
                    check_deriv=TRUE,       # outer estim without refit = no leverages checks... => it seems better to always check them
                    verbose=FALSE,
                    # method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE),
                    ...) {
  ## We need an X_off_fn so that the etaFix is used to build an offset. 
  ## IRLS function do not really handle etaFix. We need an etaFix at preprocessing stage so that columns are suppressed from AUGI0_ZX$X.pv
  ## => => get_HLCorcall(fitobject, ... etaFix=list(beta=fixef(fitobject)))
  ## Currently X_off_fn is set by .preprocess() only given an init beta, not a given beta, => we set up it in this function
  
  ### REML: the resulting SEs are consistent with those from the beta table... (with keepInREML used in numInfo())
  # 
  # if (is_REML) { # perhaps add an 'objective' argument ?
  #   if ("beta" %in% which) {
  #     warning("numInfo() computation on REML fits:\n derivatives of log restricted likelihood wrt 'beta' coefficients may not be meaningful.", immediate.=TRUE)
  #     # I could implement computation of separate blocks for p_bv vs p_v-maximizing parameters. What about cross-blocks of second derivatives ?
  #   } else if (verbose) message("numInfo() computation on REML fits: derivatives are those of log restricted likelihood.", immediate.=TRUE)
  # }
  # 
  is_REML <- .REMLmess(fitobject,return_message=FALSE)
  if (is.null(which)) {
    if (is_REML) {
      which <- c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "beta_prec")
    } else which <- c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "beta_prec", "beta")
  }
  # if (is_REML && "beta" %in% which) {
  #   REMLformula <- formula(fitobject)
  # } else REMLformula <- NULL
  
  where_rP <- .get_fittedPars(fitobject, fixef=FALSE, verbose=verbose) 
  for_which_rP <- .get_fittedPars(fitobject, which=which, fixef=FALSE, verbose=verbose) # may be zero-length
  if ("beta" %in% which) {beta_eta <- structure(fixef(fitobject), keepInREML=TRUE)} else beta_eta <- NULL
  hlcorcall <- get_HLCorcall(fitobject, 
                             fixed=where_rP, # If one use 'for_which_rP' here, parameters still get fixed cf (hlcorcall$fixed),
                                             # but to what is usually  default initial values (of not interest here)
                                             # This is a "feature", not a design decision....
 # fitmv() does not handle REMLformula *outside* the submodels, so either fitmv should be modified, or
 # the correct REMLformula must be automatically set as side effect of the keepInREML attribute,
 # or (chosen solution) the required effect must be achieved by keepInREML => cf .preprocess_X_XRe_off(), and related code for fitmv()
 #                           REMLformula=REMLformula,
                             etaFix=list(beta=beta_eta)) # => the etaFix argument means that processed$off is modified by it 
  proc_info <- .post_process_hlcorcall(hlcorcall, ranpars=for_which_rP, beta_eta=beta_eta, fitobject=fitobject)
  if (is.null(check_deriv)) check_deriv <- (
      length(for_which_rP$lambda) && any(for_which_rP$lambda<1e-6) ||
      ( ! is.null(fitobject$warnings$allLeveLam1))
    ) 
  #
  skeleton <- for_which_rP
  if (length(beta_eta)) skeleton$etaFix$beta <- beta_eta
  if (verbose) print(skeleton)
  if (check_deriv) skeleton <- .check_numDeriv_task(skeleton, .numInfo_objfn, hlcorcall, transf, proc_info, 
                                                        moreargs=.get_moreargs(fitobject), ...) # assumes untransformed param
  partially_fixed <- attr(skeleton, "partially_fixed") # may be zero-length (___F I X M E use this info to more elegantly fix thos pars in the numDeriv:: calls?)
  # Final value always refer to untransformed params:
  parnames <- c(names(unlist(skeleton[setdiff(names(skeleton), "etaFix")])), names(skeleton$etaFix$beta))
  if (transf) {
    canon_skeleton <- skeleton
    skeleton <- .ad_hoc_trRanpars(for_which_rP)
    if (length(beta_eta)) {
      skeleton$etaFix$trBeta <- .betaFn(beta_eta)
    }
    #  and the .numInfo_objfn must back-transform parameters
  }
  #
  
  resu <- hessian(func = .numInfo_objfn, x = unlist(skeleton), 
                  skeleton=skeleton, hlcorcall=hlcorcall, 
                  transf=transf, # whether x is on untransformed scale
                  objective=proc_info$objective, 
                  moreargs=.get_moreargs(fitobject), ...)
  # .assignWrapper(processed, paste0("return_only <- NULL")) # Not necess bc $processed is created by the local get_HLCorcall() and freed when this fn is exited. 
  if (transf) {
    jacTransf <- .ad_hoc_jac_transf(canon_skeleton)
    resu <- crossprod(jacTransf, resu %*% jacTransf)
    is_REML <- .REMLmess(fitobject,return_message = FALSE) # used anyway, say gradient comput.
    if (is_REML) { # grad_obj should be ~ 0 , for ML fits at least, so this computation should not be necessary
      grad_obj <- grad(func = .numInfo_objfn, x = unlist(skeleton), # on transformed scale
                       skeleton=skeleton, hlcorcall=hlcorcall, 
                       transf=transf, 
                       objective=proc_info$objective, 
                       moreargs=.get_moreargs(fitobject), ...)
      grad_list <- relist(grad_obj, skeleton) 
      gr_obj_X_hess_transf <- .ad_hoc_grXhessians_transf(canon_skeleton, grad_list) # 2nd deriv of each parameter transform.
      resu <- resu+gr_obj_X_hess_transf
    }
  }
  dimnames(resu) <- list(parnames, parnames)
  if ( length(partially_fixed)) { 
    parnames <- setdiff(parnames, partially_fixed)
    resu <- resu[parnames, parnames, drop=FALSE]
  }
  resu
}


# call by (private) .score_test and should presumably notbe called on REML fits; otherwise  (___F I X M E___) rethink the REML case
# currently tested by tests-private/test-score-test.R
.score <- function(fitobject, 
                   fixed_H0=list(),
                   etaFix_H0=list(),
                   refit_H0=update(fitobject, fixed=fixed_H0, etaFix=etaFix_H0),
                   transf, 
                   which=c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "beta_prec", "beta"),
                   verbose=FALSE,
                   ...) {
  where_rP <- .get_fittedPars(refit_H0, fixef=FALSE, verbose=verbose) # ma be NULL if all are fixed
  where_rP <- .modify_list(where_rP,fixed_H0) # merge fixed and fitted ones
  if ("beta" %in% which) {beta_eta <- structure(fixef(refit_H0), keepInREML=TRUE)} else beta_eta <- NULL 
  hlcorcall <- get_HLCorcall(fitobject, 
                             fixed=where_rP, 
                             etaFix=list(beta=beta_eta)) 
  for_which_rP <- where_rP[intersect(names(where_rP),which)]
  proc_info <- .post_process_hlcorcall(hlcorcall, ranpars=for_which_rP, beta_eta=beta_eta, fitobject=fitobject)
  
  skeleton <- for_which_rP
  if (length(beta_eta)) skeleton$etaFix$beta <- beta_eta
  if (verbose) print(skeleton)
  side <- .calc_grad_side_arg(skeleton)

  if (transf) {
    canon_skeleton <- skeleton
    skeleton <- .ad_hoc_trRanpars(for_which_rP)
    if (length(beta_eta)) {
      skeleton$etaFix$trBeta <- .betaFn(beta_eta)
    }
    #  and the .numInfo_objfn must back-transform parameters
  }

  resu <- grad(func = .numInfo_objfn, x = unlist(skeleton), side=unlist(side), skeleton=skeleton, hlcorcall=hlcorcall, 
             transf=transf, objective=proc_info$objective, 
             moreargs= .get_moreargs(fitobject), ...)
  
  if (transf) {
    jacTransf <- .ad_hoc_jac_transf(canon_skeleton)
    resu <- resu %*% jacTransf
  }
  names(resu) <- c(names(unlist(for_which_rP)), names(beta_eta))
  resu
}

if (FALSE) {
  # using a 'singfit' from test-rank.R
  spaMM:::.score(singfit,fixed_H0=list(lambda=c("1"=0.01)), etaFix_H0 = list(beta=fixef(singfit)+1),transf=FALSE)
  spaMM:::.score(singfit,fixed_H0=list(lambda=c("1"=0.01)), etaFix_H0 = list(beta=fixef(singfit)+1),transf=TRUE)
}

# Currently not used:
.numInfo_at <- function(fitobject, 
                   fixed_H0=list(),
                   etaFix_H0=list(),
                   transf, 
                   which=NULL,
                   verbose=FALSE,
                   ...) {
  refit_H0 <- update(fitobject, fixed=fixed_H0, etaFix=etaFix_H0)
  where_rP <- .get_fittedPars(refit_H0, fixef=FALSE, verbose=verbose) # ma be NULL if all are fixed
  where_rP <- .modify_list(where_rP,fixed_H0) # merge fixed and fitted ones
  if ("beta" %in% which) {beta_eta <- structure(fixef(refit_H0), keepInREML=TRUE)} else beta_eta <- NULL
  is_REML <- .REMLmess(fitobject,return_message=FALSE)
  if (is.null(which)) {
    if (is_REML) {
      which <- c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "beta_prec")
    } else which <- c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "beta_prec", "beta")
  }
  if (is_REML && "beta" %in% which) {
    REMLformula <- formula(fitobject)
  } else REMLformula <- NULL
  
  hlcorcall <- get_HLCorcall(fitobject, 
                             fixed=where_rP,
                             REMLformula=REMLformula,
                             etaFix=list(beta=beta_eta)) 
  for_which_rP <- where_rP[intersect(names(where_rP),which)]
  proc_info <- .post_process_hlcorcall(hlcorcall, ranpars=for_which_rP, beta_eta=beta_eta, fitobject=fitobject)
  
  skeleton <- for_which_rP
  if (length(beta_eta)) skeleton$etaFix$beta <- beta_eta
  if (verbose) print(skeleton)
  
  side <- relist(rep(NA,length(.unlist(skeleton))),skeleton) 
  side$lambda[skeleton$lambda<5e-5] <- 1
  
  if (transf) {
    canon_skeleton <- skeleton
    skeleton <- .ad_hoc_trRanpars(for_which_rP)
    if (length(beta_eta)) {
      skeleton$etaFix$trBeta <- .betaFn(beta_eta)
    }
    #  and the .numInfo_objfn must back-transform parameters
  }
  
  ###### finishes as numInfo... except test on is_REML removed...
  
  resu <- hessian(func = .numInfo_objfn, x = unlist(skeleton), 
                  skeleton=skeleton, hlcorcall=hlcorcall, 
                  transf=transf, # whether x is on untransformed scale
                  objective=proc_info$objective, 
                  moreargs=.get_moreargs(fitobject), ...)
  # .assignWrapper(processed, paste0("return_only <- NULL")) # Not necess bc $processed is created by the local get_HLCorcall() and freed when this fn is exited. 
  if (transf) {
    # Assuming each tr_i is function of a cn_i only:
    #   d2_cn1,cn2 logL(cn) = (d2_tr logL(tr)) (d_cn1 tr1) (d_cn2 tr2) + (d_tr1 logL(tr)) (d2_cn1,cn2 tr1)  
    #                  =                  .                            + nonzero only for cn1 ~ cn2  =>   (d_tr1 logL(tr)) (d2_cn1 tr1)  
    # otherwise second term is rather d_cn2 ((grad_tr logL).(jac_cn1 tr)) => overall (grad_tr logL) %*% hessian(transfo)
    jacTransf <- .ad_hoc_jac_transf(canon_skeleton)
    resu <- crossprod(jacTransf, resu %*% jacTransf) # (d2_tr logL(tr)) (d_cn1 tr1) (d_cn2 tr2)
    grad_obj <- grad(func = .numInfo_objfn, x = unlist(skeleton), # on transformed scale
                     skeleton=skeleton, hlcorcall=hlcorcall, 
                     transf=transf, 
                     objective=proc_info$objective, 
                     moreargs=.get_moreargs(fitobject), ...)
    grad_list <- relist(grad_obj, skeleton) 
    diag_hess_transf <- .ad_hoc_grXhessians_transf(canon_skeleton, grad_list) # 2nd deriv of each parameter transform.
    resu <- resu+diag(x=diag_hess_transf*grad_obj) 
  }
  parnames <- c(names(unlist(for_which_rP)), names(beta_eta))
  dimnames(resu) <- list(parnames, parnames)
  resu
  
}

if (FALSE) {
  (null_numInfo <-   spaMM:::.numInfo_at(spfit,fixed_H0 = list(lambda=c(0.2084315, 0.0342933), phi=spfit$phi), etaFix_H0 = list(beta=fixef(spfit)), transf=TRUE))
  (null_numInfo <-   spaMM:::.numInfo_at(spfit,fixed_H0 = list(lambda=spfit$lambda, phi=spfit$phi), etaFix_H0 = list(beta=fixef(spfit)), transf=TRUE))
}

# Compare to solve(merDeriv::vcov(., full=TRUE, information="expected"))
.expInfo <- function (object) {
  beta_expI <- solve(vcov(object)) # well there's nothing here when beta is fixed...
  logdisp_cov <- .get_logdispObject(object)$logdisp_cov
  logdisp_expI <- solve(logdisp_cov)
  ranpars_phi <- .get_ranPars_phi(object, wo_fixed=FALSE)
  precs <- 1/c(ranpars_phi$lambda,ranpars_phi$phi) # for ranCoefs, this fails bc the logdisp_cov hass elements fro the 'lambda' but there is no $lambda
  precmat <- diag(x=precs, ncol=length(precs))
  disp_expI <- precmat %*% logdisp_expI %*% precmat
  Matrix::bdiag(beta_expI,disp_expI)
}

.expInfo_at <- function(fitobject, 
                        fixed_H0=list(),
                        etaFix_H0=list(),
                        transf, 
                        which=c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "beta_prec", "beta"),
                        verbose=FALSE,
                        ...) {
  refit4beta_expI <- update(fitobject, fixed=fixed_H0)
  beta_expI <- solve(vcov(refit4beta_expI))
  refit_H0 <- update(fitobject, fixed=fixed_H0, etaFix=etaFix_H0)
  logdisp_cov <- .calc_logdispObject(refit_H0, force_fixed=TRUE)$logdisp_cov
  logdisp_expI <- solve(logdisp_cov)
  ranpars_phi <- .get_ranPars_phi(refit_H0, wo_fixed=FALSE)
  precs <- 1/unlist(ranpars_phi[c("lambda","phi")]) # for ranCoefs, this fails bc the logdisp_cov hass elements fro the 'lambda' but there is no $lambda
  precmat <- diag(x=precs, ncol=length(precs))
  disp_expI <- precmat %*% logdisp_expI %*% precmat
  resu <- Matrix::bdiag(beta_expI,disp_expI)
  rownames(resu) <- colnames(resu) <- c(colnames(beta_expI), names(precs))
  resu
}

.score_test <- function(fitobject, # ___F I X M E___ API through LRT option ? protect against calling it on REML fits?
                        y=response(fitobject), 
                        X=fitobject$X.pv, 
                        Z=get_ZALMatrix(fitobject), 
                        fixed_H0=list(),
                        etaFix_H0=list(), ...){
  fixed_H0$lambda <- pmax(1e-6,fixed_H0$lambda)
  fixed_H0$lambda <- .reformat_lambda(fixed_H0$lambda, nrand=length(fixed_H0$lambda), full_lambda=FALSE)
  nullfit <- update(fitobject, fixed=fixed_H0, etaFix=etaFix_H0)
  null_beta <- fixef(nullfit)
  full_lam_phi <- .get_ranPars_phi(fitobject, wo_fixed=FALSE)[c("lambda","phi")] # order lambda,phi reversed // lmmstest !
  triv_scale <- rep(c(1,2,2), c(length(null_beta), length(full_lam_phi$lambda), length(full_lam_phi$phi) )) # order lambda,phi !
  
  varnames <- c(names(null_beta),names(unlist(full_lam_phi)))
  
  score <-  .score(fitobject, refit_H0=nullfit, fixed_H0 = fixed_H0, etaFix_H0 = list(beta=null_beta), transf=FALSE)[varnames]
  
  null_expInfo <- .expInfo_at(fitobject,fixed_H0 = fixed_H0, etaFix_H0 = list(beta=null_beta), transf=FALSE)[varnames,varnames]
  diag_scale <- diag(triv_scale)
  
  # lmmstest's 'loglik' component is logLik(nullfit); it is not used here.
  scaled_score <- - score * triv_scale
  scaled_inf <- diag_scale %*% null_expInfo %*% diag_scale
  
  np <- length(varnames)
  testnames <- (varnames %in% c(names(unlist(fixed_H0)), names(etaFix_H0$beta)))
  test_idx <- which(testnames)
  fitted_lam_phi <- .get_ranPars_phi(fitobject, wo_fixed=TRUE)[c("lambda","phi")] # order lambda,phi !
  
  fixedpars <- setdiff(names(unlist(full_lam_phi)),names(unlist(fitted_lam_phi)))
  fix_idx <- which(varnames %in% fixedpars)
  if (length(pb <- intersect(fix_idx,test_idx))) {
    stop(paste("Parameters",paste(varnames[pb],collapse=", "),"fixed both by 'fixed_H0' and in the full model."))
  }
  no_idx <- (1L:np)[-c(test_idx, fix_idx)]
  s <- scaled_score[test_idx]
  I_block <- scaled_inf[test_idx, test_idx]
  if(length(no_idx)){
    I_block <- I_block -
      scaled_inf[test_idx, no_idx] %*% solve(scaled_inf[no_idx, no_idx], scaled_inf[no_idx, test_idx])
  }
  test_stat <- sum(s * qr.solve(I_block, s))
  df <- length(test_idx)
  p_val <- stats::pchisq(q = test_stat, df = df, lower.tail = FALSE)
  c("chi_sq" = test_stat, "df" = df, "p_val" = p_val)
}
  
