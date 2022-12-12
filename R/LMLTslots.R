
.get_covbeta <- function(x, hlcorcall, skeleton, transf, fitobject) {
  if (TRUE) {
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

.calc_grad_side_arg <- function(skeleton, tol=1e-5) { # on canonical skeleton
  side <- relist(rep(NA,length(.unlist(skeleton))),skeleton) 
  side$lambda[skeleton$lambda<tol] <- 1
  if ( ! is.null(rCs <- skeleton$ranCoefs)) {
    diagbools <- relist(rep(FALSE,length(.unlist(rCs))),rCs) 
    for (rc_it in seq_along(rCs)) {
      np <- length(rCs[[rc_it]])
      Xi_ncol <- floor(sqrt(np*2L))
      vdiagPos <- cumsum(c(1L,rev(seq(Xi_ncol-1L)+1L))) # diagpos on vector repre of half matrix, not on matrix
      diagbools[[rc_it]][vdiagPos] <- TRUE
    }
    diagbools <- .unlist(diagbools)
    urCs <- unlist(rCs)
    rC_side <- rep(NA,length(urCs))
    rC_side[diagbools & urCs <tol] <- 1
    rC_side[( ! diagbools) & urCs < -1+tol] <- 1
    rC_side[( ! diagbools) & urCs >  1-tol] <- -1 # amusingly, one can see 'improved' p_bv for corr >1 (mv zut0 test, refitted by REML)
    side$ranCoefs <- relist(rC_side,rCs)
  }
  side
}

.check_numDeriv_task <- function(skeleton, .numInfo_objfn, hlcorcall, transf, proc_info, moreargs, ...) {
  side <- .calc_grad_side_arg(skeleton)
  uside <- unlist(side)
  gr_neg_APHL <- grad(func = .numInfo_objfn, x = unlist(skeleton), side=uside, skeleton=skeleton, hlcorcall=hlcorcall, 
             transf=transf, objective=proc_info$objective, moreargs=moreargs, ...)
  
  removand <- rep(FALSE, length(uside))
  uside[is.na(uside)] <- 0
  thr <- 0.1 # ideally the threshold should depend on the internal steps of grad, which are O(1e-6) but "complicated" and not available in output
  # so a gradient of -1 for a change of 1e-6 of the argument means the APHL appears to increase by ~ 1e-6. Lower values may be negligible 
  if (any(nonzero_grad <- abs(gr_neg_APHL)>thr)) {
    max_at_bound <- (uside==-1L & gr_neg_APHL < -thr) | (uside==1L & gr_neg_APHL > thr)
    if (any(max_at_bound)) {
      warning(paste(proc_info$objective, "has nonzero gradient for",
                    paste(names(unlist(skeleton))[which(max_at_bound)],collapse=", "), 
                    "value(s) near a boundary of their range.\n  These parameter(s) may be ignored in the results."), 
              immediate. = TRUE)
      removand[which(max_at_bound)] <- TRUE # if there is no info on lambda, it is retained...
    } 
  }
  if (any(other_nonzero_grad <- (nonzero_grad & ! removand) )) {
    # nonzero gradient not diagnosed as being at bound (
    # * might still be at bound but .calc_grad_side_arg() does not know that;
    # For ranCoefs, 
    if (any(other_nonzero_grad)) {
      warning(paste(proc_info$objective, "has nonzero gradient for",
                    paste(names(unlist(skeleton))[which(other_nonzero_grad)],collapse=", "), 
                    "value(s).\n  These parameter(s) may be ignored in the results."), 
              immediate. = TRUE)
      removand[which(other_nonzero_grad)] <- TRUE
    } 
  }
  if (any(other_pb_at_bound <- (uside !=0 & ! removand))) { # eg for flat lambda at zero bound, risk of producing negative values in hessian comput
    # Only lambda is actually checked here but might need to extend this later
    warning(paste(paste(names(unlist(skeleton))[which(other_pb_at_bound)],collapse=", "), 
                  "value(s) near a boundary of their range.\n  These parameter(s) may be ignored in the results."), 
            immediate. = TRUE)
    removand[which(other_pb_at_bound)] <- TRUE
  } 

  tmp <- unlist(skeleton)
  tmp[removand] <- NaN
  tmp <- relist(tmp,skeleton)
  partially_fixed_rC <- tmp["ranCoefs"]
  if ( ! is.null(rCs <- tmp$ranCoefs)) {
    for (rc_it in seq_along(rCs)) {
      nans <- is.nan(rCs[[rc_it]])
      if (any(nans) && ! all(nans)) { # partially fixed ranCoefs: as in fits, we keep them in the variable params 
        rCs[[rc_it]] <- skeleton$ranCoefs[[rc_it]] # won't be removed by .rmNaN()
      }
    }
    tmp$ranCoefs <- rCs
  }
  skeleton <- .rmNaN(tmp)
  if ( ! length(skeleton)) stop("No fitted (co-)variance parameters whose information matrix could be evaluated.")
  # In some case we might need to regerate hlcorcall and proc_info...
  
  attr(skeleton,"partially_fixed") <- names(which(is.nan(unlist(partially_fixed_rC))))
  skeleton
}


# reproduces the elements added to the lmerMod object by lmerTest:
.calc_numderivs_LMLT <- function(fitobject, tol = 1e-08, skeleton=NULL, transf, check_deriv=NULL) {
  if (is.null(skeleton)) skeleton <- .get_ranPars_phi(fitobject, wo_fixed=TRUE)
  if ( ! length(skeleton)) stop("No fitted (co-)variance parameters whose information matrix could be evaluated.")
  hlcorcall <- get_HLCorcall(fitobject, fixed=skeleton) # using skeleton on canonical scale
  #
  proc_info <- .post_process_hlcorcall(hlcorcall, ranpars=skeleton) # modifies the $processed environment
  if (is.null(check_deriv)) check_deriv <- (length(skeleton$lambda) && any(skeleton$lambda<1e-6)) 
  if (check_deriv) skeleton <- .check_numDeriv_task(skeleton, .numInfo_objfn, hlcorcall, transf, proc_info, 
                                                    moreargs=.get_moreargs(fitobject))
  #
  
  nuisance <- skeleton
  if (transf) skeleton <- .ad_hoc_trRanpars(skeleton)
  res <- list(skeleton=skeleton, nuisance=nuisance, # nuisance should be on canonical scale as it is displayed; skeleton is typically transformed.
              vcov_beta=vcov(fitobject))
  h <- hessian(func = .numInfo_objfn, x = unlist(skeleton), skeleton=skeleton, hlcorcall=hlcorcall, transf=transf, objective=proc_info$objective)
  eig_h <- eigen(h, symmetric = TRUE)
  evals <- eig_h$values
  neg <- evals < -tol
  pos <- evals > tol
  if (sum(neg) > 0) {
    eval_chr <- if (sum(neg) > 1) 
      "eigenvalues"
    else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
    warning(sprintf("Fit failed to converge with %d negative %s of information matrix: %s", 
                    sum(neg), eval_chr, evals_num), call. = FALSE)
  }
  zero <- evals > -tol & evals < tol
  if (sum(zero) > 0) {
    eval_chr <- if (sum(zero) > 1) 
      "eigenvalues"
    else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
    warning(sprintf("Fit may not have converged, or information may be missing about some (combinations of) parameters,\n   with %d %s of information matrix close to zero: %s", 
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
  res$Jac_list <- lapply(1:ncol(Jac), function(i) array(Jac[, i], dim = rep(length(fixef(fitobject, na.rm=TRUE)), 2)))
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

as_LMLT <- function(fitobject, nuisance=NULL, verbose=TRUE, transf=TRUE, check_deriv=NULL, ...) {
  numderivs <- .calc_numderivs_LMLT(fitobject, skeleton=nuisance, transf=transf, check_deriv=check_deriv) 
  vcov_beta <- numderivs$vcov_beta
  class(vcov_beta) <- "matrix"
  sigma <- sqrt(residVar(fitobject, which="fit"))
  fixef_terms <- fitobject$main_terms_info$fixef_terms
  resu <- new("LMLTslots", 
              nuisance=numderivs$nuisance,
              beta=fixef(fitobject, na.rm=TRUE),  
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
