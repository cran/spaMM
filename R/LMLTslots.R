
.get_covbeta <- function(x, hlcorcall, skeleton, transf, fitobject) {
  if (E_S_S_A_I <- TRUE) {
    refit <- .numInfo_objfn(x, hlcorcall=hlcorcall, skeleton=skeleton, 
                            transf, # signals transformed input. T or F when called from numInfo(), 
                            # TRUE when called from as_LMLT <- function(., transf=TRUE) as the argument is passed all the way down to here
                            objective=NULL)
    vcov(refit)
  } else {
    ranpars <- relist(x, skeleton)
    if (transf) ranpars <- .canonizeRanPars(ranpars, corr_info=fitobject$ranef_info$sub_corr_info, checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf)
    # a merge parlists may be necessary?
    fittingfn <- .get_bare_fnname.HLfit(fitobject)
    refit <- update(fitobject, fixed=ranpars )    
    vcov(refit)
  }
}

# reproduces the elements added to the lmerMod object by lmerTest:
.calc_numderivs_LMLT <- function(fitobject, tol = 1e-08, skeleton=NULL, transf) {
  if (is.null(skeleton)) skeleton <- .get_fittedPars_ran_phi(fitobject)
  if ( ! length(skeleton)) stop("No fitted (co-)variance parameters whose information matrix could be evaluated.")
  #
  hlcorcall <- get_HLCorcall(fitobject, fixed=skeleton) # using skeleton on canonical scale
  proc_info <- .post_process_hlcorcall(hlcorcall, ranpars=skeleton) # modifies the $processed environment
  
  nuisance <- skeleton
  if (transf) skeleton <- .ad_hoc_trRanpars(skeleton)
  res <- list(skeleton=skeleton, nuisance=nuisance, # nuisance should be on canonical scale as it is displayed; skeleton is typically transformed.
              vcov_beta=vcov(fitobject))
  h <- hessian(func = .numInfo_objfn, x = unlist(skeleton), skeleton=skeleton, hlcorcall=hlcorcall, transf=transf, objective=proc_info$objective)
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
  .assignWrapper(hlcorcall$processed, paste0("return_only <- NULL"))  # currently we call vcov() in .get_covbeta() so we need a full return object (but this could be improved) 
  Jac <- jacobian(func = .get_covbeta, x = unlist(skeleton), hlcorcall=hlcorcall, skeleton=skeleton, 
                            # fitobject=fitobject, 
                            transf=transf)
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
  if (is.null(getClassDef("LMLT", where = .spaMM.data$class_cache, inherits = FALSE))) {
    # if (requireNamespace("lmerTest",quietly=TRUE)) { 
      # I wrote comments to the effect that this block failed if lmerTest is not loaded (hence the projected requireNamespace()).  
      # The "glitch" would occur running  "test-anova-contest.R" (now a block in test-ANOVA-&-lmerTest?) without its requireNamespace...
      # but currently unable to reproduce the problem...
      setClass(Class="LMLT", contains = c("LMLTslots","lmerModLmerTest"), where=.spaMM.data$class_cache)
      setClass(Class="LMLTinternal", contains = c("LMLTslots","lmerModLmerTest"), where=.spaMM.data$class_cache) # to avoid infinite recursions see e.g. anova.LMLT
    #} 
  }
  #
  if (verbose) {
    cat("Result accounting for uncertainty for the following estimates of nuisance parameters:\n")
    print(unlist(resu@nuisance))
  }
  resu <- as(resu, "LMLT") # (___F I X M E___) alternative would be to recode the lmerTest algo for the contrast matrix 'L', the rest seeming straightforward.
  invisible(resu)
}
