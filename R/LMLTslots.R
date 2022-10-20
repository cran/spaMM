.infofn_fixed <- function(x, fitobject, skeleton, transf) {
  ranpars <- relist(x, skeleton)
  if (transf) ranpars <- .canonizeRanPars(ranpars, corr_info=fitobject$ranef_info$sub_corr_info, checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)
  # a merge parlists may be necessary
  # print(ranpars)
  refit <- update(fitobject, fixed=ranpars ) # ___F I X M E___ a way to save time by update returning only a logLik?
  - logLik(refit)  
}

.infofn_ranFix <- function(x, fitobject, skeleton, transf) {
  ranpars <- relist(x, skeleton)
  if (transf) ranpars <- .canonizeRanPars(ranpars, corr_info=fitobject$ranef_info$sub_corr_info, checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)
  # a merge parlists may be necessary
  refit <- update(fitobject, ranFix=ranpars ) 
  - logLik(refit)  
}

.numInfo <- function(fitobject, skeleton=.get_fittedPars_ran_phi(fitobject), transf) {
  fittingfn <- .get_bare_fnname.HLfit(fitobject)
  if ("fixed" %in% names(formals(fittingfn))) {
    infofn <- .infofn_fixed
  } else infofn <- .infofn_ranFix
  numDeriv::hessian(func = infofn, x = unlist(skeleton), skeleton=skeleton, fitobject=fitobject, transf=transf)
}

.get_covbeta <- function(x, fitobject, skeleton, transf) {
  ranpars <- relist(x, skeleton)
  if (transf) ranpars <- .canonizeRanPars(ranpars, corr_info=fitobject$ranef_info$sub_corr_info, checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)
  # a merge parlists may be necessary?
  refit <- update(fitobject, fixed=ranpars )
  vcov(refit)
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

# reproduces the elements added to the lmerMod object by lmerTest:
.calc_numderivs_LMLT <- function(fitobject, tol = 1e-08, skeleton=NULL, transf) {
  if (is.null(skeleton)) skeleton <- .get_fittedPars_ran_phi(fitobject)
  if ( ! length(skeleton)) stop("No fitted (co-)variance parameters whose information matrix could be evaluated.")
  nuisance <- skeleton
  if (transf) skeleton <- .ad_hoc_trRanpars(skeleton)
  res <- list(skeleton=skeleton, nuisance=nuisance, # nuisance should be on canonical scale as it is displayed; skeleton is typically transformed.
              vcov_beta=vcov(fitobject))
  h <- .numInfo(fitobject, skeleton = skeleton, transf=transf)
  eig_h <- eigen(h, symmetric = TRUE)
  evals <- eig_h$values
  neg <- evals < -tol
  pos <- evals > tol
  zero <- evals > -tol & evals < tol
  if (sum(neg) > 0) {
    eval_chr <- if (sum(neg) > 1) 
      "eigenvalues"
    else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
    warning(sprintf("Model failed to converge with %d negative %s: %s", 
                    sum(neg), eval_chr, evals_num), call. = FALSE)
  }
  if (sum(zero) > 0) {
    eval_chr <- if (sum(zero) > 1) 
      "eigenvalues"
    else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
    warning(sprintf("Model may not have converged with %d %s close to zero: %s", 
                    sum(zero), eval_chr, evals_num))
  }
  pos <- eig_h$values > tol
  q <- sum(pos)
  h_inv <- with(eig_h, {
    vectors[, pos, drop = FALSE] %*% diag(1/values[pos], 
                                          nrow = q) %*% t(vectors[, pos, drop = FALSE])
  })
  res$vcov_varpar <- h_inv 
  Jac <- numDeriv::jacobian(func = .get_covbeta, x = unlist(skeleton), skeleton=skeleton, 
                            fitobject=fitobject, transf=transf)
  res$Jac_list <- lapply(1:ncol(Jac), function(i) array(Jac[, i], dim = rep(length(fixef(fitobject)), 2)))
  res
}



setClass("LMLTslots",
         representation = representation(
           nuisance="list", # so that this can be examined post-analysis. Not a slot of an lmerModLmerTest object
           beta="numeric",
           frame="data.frame", # for anova()
           X.pv="matrix", # for anova() -> model.matrix
           vcov_varpar="matrix",  
           Jac_list= "list", 
           vcov_beta= "matrix", 
           sigma = "numeric"))

# Without the next method, show() would apply the method for a merMod object (-> bug)
setMethod(f = "show",
          signature = "LMLTslots",
          definition = function(object){
            print(paste("S4 object with the following slots:", paste(slotNames("LMLTslots"), collapse=", ")))
          })

model.matrix.LMLTslots <- function(object, ...) object@X.pv

as_LMLT <- function(fitobject, nuisance=NULL, verbose=TRUE, transf=TRUE, ...) {
  numderivs <- .calc_numderivs_LMLT(fitobject, skeleton=nuisance, transf=transf) 
  vcov_beta <- numderivs$vcov_beta
  class(vcov_beta) <- "matrix"
  sigma <- sqrt(residVar(fitobject, which="fit"))
  fixef_terms <- fitobject$main_terms_info$fixef_terms
  resu <- new("LMLTslots", 
              nuisance=numderivs$nuisance,
              beta=fixef(fitobject), 
              frame=structure(model.frame(fixef_terms,fitobject$data), # contains  I(Days^2), in contrast to fitobject$data
                              formula=formula(fitobject),
                              predvars.fixed=attr(fixef_terms,"variables")), 
              X.pv=fitobject$X.pv,
              vcov_varpar=numderivs$vcov_varpar,
              Jac_list=numderivs$Jac_list,
              vcov_beta=vcov_beta,
              sigma=sigma )
  #
  # Got inspiration from https://github.com/r-dbi/odbc/commit/992d9070e51590d709483a248d119a1b70117c85 to avoid "cannot add bindings to a locked environment" stuff
  if (is.null(getClassDef("LMLT", where = .spaMM.data$class_cache, inherits = FALSE))) 
    setClass(Class="LMLT", contains = c("LMLTslots","lmerModLmerTest"), where=.spaMM.data$class_cache)
  #
  if (verbose) {
    cat("Result accounting for uncertainty for the following estimates of nuisance parameters:\n")
    print(unlist(resu@nuisance))
  }
  resu <- as(resu, "LMLT")
  invisible(resu)
}
