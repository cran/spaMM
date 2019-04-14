# y-augmented method withOUT precomputation of R_aug_ZXy
.HLfit_body_augZXy <- function(processed, ranFix=list()) { 
  ranFix <- .post_process_family(processed$family,ranFix) ## assign 'extra' pars and cleans ranFix
  ranFix <- .canonizeRanPars(ranPars=ranFix,corr_info=NULL, checkComplete = FALSE)## including full-size lambda
  nobs <- nrow(processed$data) ## before prior.weights is evaluated
  ### a bit of post processing
  nobs <- NROW(processed$AUGI0_ZX$X.pv)
  LMMbool <- processed$LMMbool
  models <- processed$models
  if (models[["eta"]]!="etaHGLM") stop("Not a mixed-effect model.") ## not handled
  nrand <- length(processed$ZAlist)
  cum_n_u_h <- processed$cum_n_u_h
  n_u_h <- cum_n_u_h[nrand+1L] 
  sparse_precision <- processed$sparsePrecisionBOOL
  ranCoefs.Fix <- .getPar(ranFix,"ranCoefs") ## may be NULL
  # Updates processed$ranCoefs_blob which contains no globally fixed ranCoefs as this has been excluded by .determine_augZXy() 
  ranCoefs_blob <- .process_ranCoefs(processed, ranCoefs.Fix,use_tri=.spaMM.data$options$use_tri_for_augZXy) ## UPDATES preexisting object
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
  if (processed$sparsePrecisionBOOL) { 
    .init_precision_info(processed,LMatrices) ## modifies processed$AUGI0_ZX$envir  
  }
  if (processed$sparsePrecisionBOOL) {
    ZAL <- NULL # we practically don't need it (though F I X M E: it would be nice to provide alternative info to .eval_init_lambda_guess)
  } else if ( any((attr(LMatrices,"is_given_by") !="")) ) {
    ZAL <- .compute_ZAL(XMatrix=LMatrices, ZAlist=processed$ZAlist,as_matrix=.eval_as_mat_arg(processed))
    ## ZAL may be modified by other call to .compute_ZAL()   
  } else { 
    ZAL <- processed$AUGI0_ZX$ZAfix ## default ZA 
  } 
  lambda.Fix <- ranFix$lambda # .getPar(ranFix,"lambda") ## should already have length 'nrand' or else be NULL
  if (any(lambda.Fix==0)) stop("lambda cannot be fixed to 0.")
  ###
  off <- processed$off
  ##################
  ## Initial estimate for lambda in 'compact" form
  init.lambda <- .calc_initial_init_lambda(lambda.Fix, nrand, processed, ranCoefs_blob, 
                                           init.HLfit=NULL, fixed=ranFix)
  # expand:
  lambda_est <- .HLfit_finalize_init_lambda(models, init.lambda, processed, ZAL=ZAL, cum_n_u_h, 
                                            vec_n_u_h=diff(cum_n_u_h), n_u_h, ranCoefs_blob)
  if (identical(processed$return_only,"p_vAPHLs")) {
    whichAPHLs <- "p_v"
  } else if (identical(processed$return_only,"p_bvAPHLs")) {
    whichAPHLs <- "p_bv" ## retrun value may still include p_v if it is used to compute p_bv
  } else whichAPHLs <- c("p_v","p_bv")
  ####################################################################################################
  # we don't want anything specific on u_h values:
  wranefblob <- .updateW_ranefS(processed$cum_n_u_h, processed$rand.families, lambda=lambda_est, 
                                u_h=rep(NA,n_u_h),v_h=rep(NA,n_u_h)) ## indeed
  muetablob <- .muetafn(eta=rep(NA,nobs),BinomialDen=processed$BinomialDen,processed=processed) 
  phi_est <- ranFix$phi 
  if (is.null(phi_est)) phi_est <- processed$phi.Fix ## not sure this is needed
  if (is.null(phi_est)) { 
    w.resid <- .calc_w_resid(muetablob$GLMweights,1)
  } else w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
  H_global_scale <- .calc_H_global_scale(w.resid)
  if (processed$sparsePrecisionBOOL) {
    # .HLfit_body_augZXy has called .init_precision_info(processed,LMatrices)...
    sXaug <- do.call(processed$AUGI0_ZX$envir$method, # ie, def_AUGI0_ZX_sparsePrecision
                       list(AUGI0_ZX=processed$AUGI0_ZX, corrPars=ranFix$corrPars,w.ranef=wranefblob$w.ranef,
                            cum_n_u_h=cum_n_u_h,w.resid=w.resid))
  } else {
    ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
    Xscal <- .make_Xscal(ZAL=ZAL, ZAL_scaling = ZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX)
    if (inherits(Xscal,"sparseMatrix")) { # 
      mMatrix_method <- .spaMM.data$options$Matrix_method
    } else {
      mMatrix_method <- .spaMM.data$options$matrix_method
    }
    weight_X <- .calc_weight_X(w.resid, H_global_scale) ## sqrt(s^2 W.resid) ## should not affect the result up to precision
    sXaug <- do.call(mMatrix_method,
                     list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
  }
  ####################################################################################################
  augZXy_resu <- .calc_APHLs_by_augZXy_or_sXaug(sXaug=sXaug, phi_est=phi_est, # may be NULL
                                             processed=processed, which=whichAPHLs)
  res <- list(APHLs=augZXy_resu)
  return(res)    ########################   R E T U R N
}


