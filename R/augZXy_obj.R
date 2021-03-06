.augZXy_obj <- function(ranefParsVec, skeleton, processed, ranFix=list(), ...) {
  ranefParsList <- relist(ranefParsVec, skeleton)
  ranFix <- .modify_list(ranFix, ranefParsList)
  rpType <- .modify_list(attr(ranFix, "type"), attr(skeleton, "type"))
  attr(ranFix, "type") <- rpType
  if (FALSE) {
    hlfit <- eval(call(.spaMM.data$options$augZXy_fitfn, processed=processed, ranFix=ranFix))
  } else {
    ranFix <- .canonizeRanPars(ranPars = ranFix, corr_info = NULL, 
                               checkComplete = FALSE, rC_transf = .spaMM.data$options$rC_transf)
    {
      
      trace <- processed$verbose["TRACE"]
      ranFix <- .canonizeRanPars(ranPars=ranFix,corr_info=NULL, checkComplete = FALSE, rC_transf=.spaMM.data$options$rC_transf)## including full-size lambda
      nobs <- length(processed$y) ## before prior.weights is evaluated
      nrand <- length(processed$ZAlist)
      cum_n_u_h <- processed$cum_n_u_h
      n_u_h <- cum_n_u_h[nrand+1L] 
      sparse_precision <- processed$is_spprec
      ranCoefs.Fix <- .getPar(ranFix,"ranCoefs") ## may be NULL
      # Updates processed$ranCoefs_blob which contains no globally fixed ranCoefs as this has been excluded by .preprocess_augZXy() 
      ranCoefs_blob <- .process_ranCoefs(processed, ranCoefs.Fix,use_tri_Nspprec=.spaMM.data$options$use_tri_for_augZXy) ## *updates* *locally* a preexisting object
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
        .update_AUGI0_ZX_envir_ranCoefs_info(processed,LMatrices)
      }
      if (processed$is_spprec) {
        ZAL <- NULL # we practically don't need it (though F I X M E: it would be nice to provide alternative info to .eval_init_lambda_guess)
      } else if ( any((attr(LMatrices,"is_given_by") !="")) ) {
        ZAL <- .compute_ZAL(XMatrix=LMatrices, ZAlist=processed$ZAlist,as_matrix=.eval_as_mat_arg(processed))
        ## ZAL may be modified by other call to .compute_ZAL()   
      } else { 
        ZAL <- processed$AUGI0_ZX$ZAfix ## default ZA 
      } 
      lambda.Fix <- ranFix$lambda # .getPar(ranFix,"lambda") ## should already have length 'nrand' or else be NULL
      if (any(lambda.Fix==0, na.rm=TRUE)) stop("lambda cannot be fixed to 0.")
      ###
      off <- processed$off
      ##################
      ## Initial estimate for lambda in 'compact" form
      init.lambda <- .calc_initial_init_lambda(lambda.Fix, nrand, processed, ranCoefs_blob, 
                                               init.HLfit=NULL, fixed=ranFix)
      # expand:
      lambda_est <- .HLfit_finalize_init_lambda(processed$models, init.lambda, processed, ZAL=ZAL, cum_n_u_h, 
                                                vec_n_u_h=diff(cum_n_u_h), n_u_h, ranCoefs_blob)
      if (identical(processed$return_only,"p_vAPHLs")) {
        whichAPHLs <- "p_v"
      } else if (identical(processed$return_only,"p_bvAPHLs")) {
        whichAPHLs <- "p_bv" ## return value may still include p_v if it is used to compute p_bv
      } else whichAPHLs <- c("p_v","p_bv")
      ####################################################################################################
      # we don't want anything specific on u_h values:
      w.ranef <- 1/lambda_est # call to .updateW_ranefS() reduced to this for v3.6.39
      muetablob <- .muetafn(eta=rep(NA,nobs),BinomialDen=processed$BinomialDen,processed=processed) 
      phi_est <- ranFix$phi 
      if (is.null(phi_est)) phi_est <- processed$phi.Fix ## not sure this is needed
      if (is.null(phi_est)) { 
        w.resid <- .calc_w_resid(muetablob$GLMweights,1)
      } else w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
      H_global_scale <- .calc_H_global_scale(w.resid)
      if (processed$is_spprec) {
        # .HLfit_body_augZXy has called .init_AUGI0_ZX_envir_spprec_info(processed,LMatrices)...
        if (trace) cat(".")
        sXaug <- do.call(processed$AUGI0_ZX$envir$method, # ie, def_AUGI0_ZX_sparsePrecision
                         list(AUGI0_ZX=processed$AUGI0_ZX, corrPars=ranFix$corrPars,w.ranef=w.ranef,
                              cum_n_u_h=cum_n_u_h,w.resid=w.resid))
      } else {
        ZAL_scaling <- 1/sqrt(w.ranef*H_global_scale) ## Q^{-1/2}/s
        weight_X <- .calc_weight_X(w.resid, H_global_scale) ## sqrt(s^2 W.resid) ## should not affect the result up to precision
        if (.spaMM.data$options$TRY_R_new && inherits(ZAL,"sparseMatrix")) {
          ZW <- .Dvec_times_Matrix(weight_X,.Matrix_times_Dvec(ZAL,ZAL_scaling))
          XW <- .Dvec_times_m_Matrix(weight_X,processed$AUGI0_ZX$X.pv)
          sXaug <- structure(list(ZW=ZW,XW=XW,I=processed$AUGI0_ZX$I),
                             w.ranef=w.ranef,
                             n_u_h=ncol(ZW), # mandatory for all sXaug types
                             pforpv=ncol(XW),  # mandatory for all sXaug types
                             weight_X=weight_X, # new mandatory 08/2018
                             H_global_scale=H_global_scale)
          class(sXaug) <- c(class(sXaug),"sXaug_blocks")
        } else {
          Xscal <- .make_Xscal(ZAL=ZAL, ZAL_scaling = ZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX) # does not weights the I
          if (inherits(Xscal,"sparseMatrix")) { # 
            mMatrix_method <- .spaMM.data$options$Matrix_method # does not weights the I
            # attr(Xscal,"AUGI0_ZX") <- processed$AUGI0_ZX # for experimental .sXaug_Matrix_cholP_scaled(). But useless bc 
            # .calc_APHLs_by_augZXy_or_sXaug won't use the methods associated to the Matrix_method.
            # Instead it uses .get_absdiagR_blocks -> .updateCHMfactor 
          } else {
            mMatrix_method <- .spaMM.data$options$matrix_method
          }
          if (trace) cat(".")
          sXaug <- do.call(mMatrix_method,
                           list(Xaug=Xscal, weight_X=weight_X, w.ranef=w.ranef, H_global_scale=H_global_scale))
        }
      }
      ####################################################################################################
      augZXy_resu <- .calc_APHLs_by_augZXy_or_sXaug(sXaug=sXaug, phi_est=phi_est, # may be NULL
                                                    processed=processed, which=whichAPHLs,
                                                    update_info=list(allow= all(processed$AUGI0_ZX$envir$updateable))) #update_info=list(allow= (! any(unlist(ranCoefs.Fix)==0))))
      hlfit <- list(APHLs=augZXy_resu)
    }
  }
  ###################################################
  objective <- processed$objective
  aphls <- hlfit$APHLs
  resu <- aphls[[objective]]
  if (objective=="cAIC") resu <- - resu ## for minimization of cAIC (private & experimental)
  if (resu>processed$augZXy_env$objective) {
    processed$augZXy_env$objective <- resu
    processed$augZXy_env$phi_est <- aphls[["phi_est"]]
  }
  return(resu) #
}
