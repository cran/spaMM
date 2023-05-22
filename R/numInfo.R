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
      hlcorcall$etaFix$beta <- .spaMM.data$options$.betaInv(trBeta) 
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
    X.pv <- model.matrix(fitobject)
    X_off <-.subcol_wAttr(X.pv, j=betanames, drop=FALSE)
    X.pv <- .subcol_wAttr(X.pv, j=setdiff(colnames(X.pv),betanames), drop=FALSE)
    ori_off <- model.offset.HLfit(fitobject) # processed$off differs from it as get_HLCorcall(., etaFix) -> .preprocess(...etaFix) adds the etaFix-derived offset.
    processed$X_off_fn <- .def_off_fn(X_off, ori_off=ori_off)
  }
  info <- list(objective=proc1$objective)
  info
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

numInfo <- function(fitobject, 
                    transf=FALSE, 
                    which=NULL,
                    check_deriv=TRUE,       # outer estim without refit = no leverages checks... => it seems better to always check them
                    sing=1e-05,
                    verbose=FALSE,
                    refit_hacks=list(),
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
      which <-      c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "COMP_nu", "beta_prec", "rdisPars")
    } else which <- c("lambda", "ranCoefs", "corrPars", "hyper", "phi", "NB_shape", "COMP_nu", "beta_prec", "rdisPars", 
                      "beta")
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
  
  hlcorcall <- get_HLCorcall(fitobject, 
                             fixed=where_rP, # If one use 'for_which_rP' here, parameters still get fixed cf (hlcorcall$fixed),
                                             # but to what is usually  default initial values (of not interest here)
                                             # This is a "feature", not a design decision....
                             # fitmv() does not handle REMLformula *outside* the submodels, so either fitmv should be modified, or
                             # the correct REMLformula must be automatically set as side effect of the keepInREML attribute,
                             # or (chosen solution) the required effect must be achieved by keepInREML => cf .preprocess_X_XRe_off(), and related code for fitmv()
                             #                           REMLformula=REMLformula,
                             etaFix=list(beta=beta_eta),  # => the etaFix argument means that processed$off is modified by it
                             verbose=refit_verbose
  ) 
  proc_info <- .post_process_hlcorcall(hlcorcall, ranpars=for_which_rP, beta_eta=beta_eta, fitobject=fitobject)
  if (is.null(check_deriv)) check_deriv <- (
      length(for_which_rP$lambda) && any(for_which_rP$lambda<1e-6) ||
      ( ! is.null(fitobject$warnings$allLeveLam1))
    ) 
  #
  skeleton <- for_which_rP
  if (length(beta_eta)) skeleton$etaFix$beta <- beta_eta
  if (verbose) print(skeleton)
  
  # cannot remove partially-fixed ranCoefs too early. They must be kept in skeleton in all cases
  ufixed <- na.omit(unlist(outer_call$fixed))
  fixednames <- names(ufixed) # .preprocess_fixed() removed any fancy names given by users.
  uskeleton <- unlist(skeleton)
  removand_rC <- intersect(fixednames,names(uskeleton))
  
  if (check_deriv) {
    thr <- .calc_grad_thr(skeleton, fitobject, beta_eta) # effect visible on numInfo(fitme(cases~1+(1|id),family=negbin1(), data=scotlip, resid.model = ~ population))
    removand <- .check_numDeriv_task(skeleton, .numInfo_objfn, hlcorcall, transf, proc_info, 
                                     moreargs=.get_moreargs(fitobject), 
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
  
  resu <- hessian(func = .numInfo_objfn, x = unlist(skeleton), 
                  skeleton=skeleton, hlcorcall=hlcorcall, 
                  transf=transf, # whether x is on untransformed scale
                  objective=proc_info$objective, 
                  moreargs=.get_moreargs(fitobject), ...)
  # .assignWrapper(processed, paste0("return_only <- NULL")) # Not necess bc $processed is created by the local get_HLCorcall() and freed when this fn is exited. 
  if (transf) {
    jacTransf <- .ad_hoc_jac_transf(canon_skeleton)
    resu <- as.matrix(crossprod(jacTransf, resu %*% jacTransf))
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

print.singeigs <- function(x, ...) {
  sing <- attr(x,"sing")
  xx <- sapply(x, function(v) {
    ifelse(v<sing, crayon::underline(signif(v)), signif(v))
  })
  cat(paste(xx))
}

