# y-augmented method withOUT precomputation of R_aug_ZXy
.HLfit_body_augZXy <- function(processed, fixed=list()) { 
  trace <- processed$verbose["TRACE"]
  ranFix <- .canonizeRanPars(ranPars=fixed,corr_info=NULL, checkComplete = FALSE, rC_transf=.spaMM.data$options$rC_transf)## including full-size lambda
  nobs <- length(processed$y) ## before prior.weights is evaluated
  # Updates processed$ranCoefs_blob which contains no globally fixed ranCoefs as this has been excluded by .preprocess_augZXy() 
  ranCoefs_blob <- .process_ranCoefs(processed, ranCoefs=.getPar(ranFix,"ranCoefs"), ## may be NULL, 
                                     use_tri_CORREL=.spaMM.data$options$use_tri_for_augZXy) ## *updates* *locally* a preexisting object
  LMatrices <- processed$AUGI0_ZX$envir$LMatrices
  # HLCor_body has prefilled $LMatrices for :
  #    for Matern...
  # or for CAR, when HLCor_body has updated LMatrix as fn of rho: seek adjd usage]
  ### we add the ranCoefs matrices:
  if (any(ranCoefs_blob$is_set)) {
    LMatrices[ranCoefs_blob$is_set] <- ranCoefs_blob$LMatrices[ranCoefs_blob$is_set]
    attr(LMatrices,"is_given_by")[ranCoefs_blob$is_set] <- "ranCoefs"
    # lambda_est initialized from ranCoefs_blob later !
  }
  if (processed$is_spprec) { 
    .init_AUGI0_ZX_envir_spprec_info(processed)
    .wrap_precisionFactorize_ranCoefs(processed,LMatrices)
  }
  if (processed$is_spprec) {
    ZAL <- NULL # we practically don't need it (though F I X M E: it would be nice to provide alternative info to .eval_init_lambda_guess)
  } else if ( any((attr(LMatrices,"is_given_by") !="")) ) {
    ZAL <- .compute_ZAL(XMatrix=LMatrices, ZAlist=processed$ZAlist, as_matrix=.eval_as_mat_arg(processed)) # with default bind.=TRUE: we must bind ZA_Kron objects, and the whole ZAL should be bindable for **corr algos
    ## ZAL may be modified by other call to .compute_ZAL()   
  } else { 
    ZAL <- processed$AUGI0_ZX$ZAfix ## default ZA 
  } 
  if (any(ranFix$lambda==0, na.rm=TRUE)) stop("lambda cannot be fixed to 0.")
  lam_fix_or_outer_or_NA <- processed$reserve$repNAnrand
  lam_fix_or_outer_or_NA[names(ranFix$lambda)] <- ranFix$lambda # LHS should then have length 'nrand' or else be NULL
  ###
  off <- processed$off
  ##################
  cum_n_u_h <- processed$cum_n_u_h
  n_u_h <- tail(cum_n_u_h,1L) 
  ## Initial estimate for lambda in 'compact" form
  init.lambda <- .calc_initial_init_lambda(lam_fix_or_outer_or_NA, # already the right size 'nrand' with NA's or non-fixed ones"
                                           processed=processed, ranCoefs_blob=ranCoefs_blob, 
                                           init.HLfit=NULL, fixed=ranFix)
  # expand:
  lambda_est <- .HLfit_finalize_init_lambda(models=processed$models, init.lambda, processed, ZAL=ZAL, cum_n_u_h, 
                                            vec_n_u_h=diff(cum_n_u_h), n_u_h, ranCoefs_blob)
  if (identical(processed$return_only,"p_vAPHLs")) {
    whichAPHLs <- "p_v"
  } else if (identical(processed$return_only,"p_bvAPHLs")) {
    whichAPHLs <- "p_bv" ## return value may still include p_v if it is used to compute p_bv
  } else whichAPHLs <- c("p_v","p_bv")
  ####################################################################################################
  # we don't want anything specific on u_h values:
  w.ranef <- 1/lambda_est # call to .updateW_ranefS() reduced to this for v3.6.37
  # call to m.muetablob reduced to elementary GLMweights computation for v3.6.37:
  GLMweights <- processed$prior.weights ## with attr(.,"unique")=TRUE ; .preprocess_augZXy() has checked that. 
  phi_est <- ranFix$phi 
  if (is.null(phi_est)) phi_est <- processed$phi.Fix ## not sure this is needed
  if (is.null(phi_est)) { 
    w.resid <- GLMweights 
  } else {
    w.resid <- GLMweights/phi_est ## 'weinu', must be O(n) in all cases
    attr(w.resid,"unique") <- TRUE
    attr(w.resid,"is_unit") <- FALSE 
  }
  H_global_scale <- .calc_H_global_scale(w.resid)
  if (processed$is_spprec) {
    # .HLfit_body_augZXy has called .init_AUGI0_ZX_envir_spprec_info(processed,LMatrices)...
    sXaug_arglist <- list(AUGI0_ZX=processed$AUGI0_ZX, corrPars=ranFix$corrPars,w.ranef=w.ranef,
                          cum_n_u_h=cum_n_u_h,H_w.resid=w.resid)
    if (trace) cat(".")
    sXaug <- do.call(processed$spprec_method, # ie, def_AUGI0_ZX_spprec
                     sXaug_arglist)
  } else {
    ZAL_scaling <- 1/sqrt(w.ranef*H_global_scale) ## Q^{-1/2}/s
    weight_X <- .calc_weight_X(Hobs_w.resid=w.resid, H_global_scale=H_global_scale, obsInfo=FALSE) ## sqrt(s^2 W.resid) ## should not affect the result up to precision
    if (inherits(ZAL,"sparseMatrix")) { 
      ZW <- .Dvec_times_Matrix(weight_X,.Matrix_times_Dvec(ZAL,ZAL_scaling))
      XW <- .Dvec_times_m_Matrix(weight_X,processed$AUGI0_ZX$X.pv)
      sXaug <- list(ZW=ZW,XW=XW,I=processed$AUGI0_ZX$I)
      attr(sXaug,"w.ranef") <- w.ranef
      attr(sXaug,"n_u_h") <- ncol(ZW) # mandatory for all sXaug types
      attr(sXaug,"pforpv") <- ncol(XW) # mandatory for all sXaug types
      attr(sXaug,"weight_X") <- weight_X # new mandatory 08/2018
      attr(sXaug,"H_global_scale") <- H_global_scale
      class(sXaug) <- c(class(sXaug),"sXaug_blocks") # current default for sparse ZAL; .calc_APHLs_by_augZXy_or_sXaug uses .get_absdiagR_blocks -> .updateCHMfactor, with QR fall-back 
    } else {
      Xscal <- .make_Xscal(ZAL=ZAL, ZAL_scaling = ZAL_scaling, processed=processed) # does not weights the I
      if (inherits(Xscal,"sparseMatrix")) { # 
        attr(Xscal,"AUGI0_ZX") <- processed$AUGI0_ZX # originally for .sXaug_Matrix_CHM_H_scaled. But ... this block is not currently used
      } 
      if (trace) cat(".")
      sXaug <- do.call(processed$corr_method,
                       list(Xaug=Xscal, weight_X=weight_X, w.ranef=w.ranef, H_global_scale=H_global_scale))
    }
  }
  ####################################################################################################
  augZXy_resu <- .calc_APHLs_by_augZXy_or_sXaug(sXaug=sXaug, phi_est=phi_est, # may be NULL
                                                processed=processed, which=whichAPHLs,
                                                update_info=list(allow= all(processed$AUGI0_ZX$envir$updateable))) #update_info=list(allow= (! any(unlist(ranCoefs.Fix)==0))))
  res <- list(APHLs=augZXy_resu) # may contain $phi_est
  return(res)    ########################   R E T U R N
}


