# hatvalues is a generic
hatvalues.HLfit <- function(model, type="projection", which="resid", force=FALSE, ...) { 
  loctype <- switch(type,
                    marginal="hatval_Z",
                    restricted="hatval",
                    projection="hatval",
                    std="fit",
                    type)
  if ( ! (loctype %in% c("fit","hatval","hatval_Z"))) stop("Unhandled 'type' value.")
  if ( ! (which %in% c("resid","ranef","both"))) stop("Unhandled 'which' value.")
  X.Re <- model$distinctX.Re
  if (model$models[["eta"]]=="etaGLM") {
    # no sXaug, no get_from_MME
    if (loctype=="fit") {
      lev <- .get_hatvalues_FM(X.Re, augX=model$X.pv, model$w.resid)
    } else { 
      lev <- .get_hatvalues_FM(model$X.pv, augX=model$X.pv, model$w.resid)
    }
  } else {
    if (is.null(X.Re)) { # standard ReML, "hatval" was used in fit
      look_object_lev <- loctype %in% c("fit","hatval")
    } else if (ncol(X.Re)==0L) { # standard ReML, "hatval_Z" was used in fit
      look_object_lev <- loctype %in% c("fit","hatval_Z")
    } else { # non-standard REML
      look_object_lev <- loctype %in% c("fit")
    }
    lev <- NULL
    if ( ( ! force) && look_object_lev) {
      if (which=="resid") lev <- model$lev_phi 
      if (which=="ranef") lev <- model$lev_lam 
      if (which=="both") {
        lev <- list(ranef=model$lev_lam,resid=model$lev_phi) 
        if (any(sapply(lev, is.null))) lev <- NULL
      }
    }
    if (is.null(lev)) {
      sXaug <- model$envir$sXaug # NOT get_matrix() since get_matrix() provides matrices with blocks in Henderson's order,
      
      # not compatible with mMatrix_method or get_from_MME    
      if (loctype=="fit") {
        hatvals <- .get_hatvalues_MM(sXaug,X.Re=X.Re, weight_X=attr(sXaug,"weight_X"))
        n_u_h <- length(model$w.ranef)
        hatvals <- list(ranef=hatvals[seq_len(n_u_h)], resid=hatvals[-seq_len(n_u_h)])
        
        rand.families <- model$rand.families
        cum_n_u_h <- attr(model$lambda.object$lambda_list,"cum_n_u_h")
        
        v_h_bounds <- .eval_v_h_bounds(cum_n_u_h, rand.families)
        u_h <- .u_h_v_h_from_v_h(model$v_h, rand.families, cum_n_u_h=cum_n_u_h,
                                 lcrandfamfam=attr(rand.families,"lcrandfamfam"),
                                 lower.v_h=v_h_bounds$lower.v_h,
                                 upper.v_h=v_h_bounds$upper.v_h)
        # We need to identify cases where phi could vary but 
        # where the estimated projection matrix is constrained away from the true one.
        # for Poisson, etc, phi is fixed but not constrained so constr_phi_fit must be FALSE
        constr_phi_fit <- identical(attr(model$phi,"constr_fit"),TRUE)
        # Likewise for lambda, we identify cases where it is fixed by user 
        ismixed <- (model$models[[1]]=="etaHGLM")
        lambdaType <- model$lambda.object$type
        if (( ismixed && any(lambdaType %in% c("fixed","fix_ranCoefs"))) || 
            constr_phi_fit) {
          message("Some dispersion parameters were constrained: using leverages for diagnostic purposes is speculative.")
          # nevertheless we try to reproduce the HLfit behaviour in that case and so we evaluate 'fully_constr_lam_fit'
        }
        fully_constr_lam_fit <-  (ismixed && all(lambdaType %in% c("fixed","fix_ranCoefs")))
        constr_phi <- identical(attr(model$phi,"constr_phi"),TRUE)
        lev <- .calc_lev_from_hat(hatvals=hatvals, sXaug, 
                                  is_null_phi.Fix= ( ! constr_phi ) , u_h=u_h, 
                                  #object=model, # 'processed' or fit object
                                  HL=model$HL, models=model$models, need_simple_lambda=( ! fully_constr_lam_fit), 
                                  muetablob=model$muetablob, family=model$family, mu=model$muetablob$mu, 
                                  BinomialDen=model$BinomialDen, w.resid=model$w.resid, 
                                  wranefblob=.updateW_ranefS(cum_n_u_h=cum_n_u_h, rand.families=rand.families, model$lambda, u_h=u_h, model$v_h), 
                                  nobs=length(model$y), ZAL=get_ZALMatrix(model), #psi_M=rep(attr(rand.families,"unique.psi_M"),diff(cum_n_u_h)), 
                                  lambda_est=model$lambda.object$lambda_est, cum_n_u_h=cum_n_u_h, #lcrandfamfam=attr(rand.families,"lcrandfamfam"), 
                                  rand.families=rand.families, y=model$y, prior.weights=model$prior.weights, #nrand=length(lcrandfamfam), 
                                  phi_est=model$phi)
        if (which=="resid") lev <- lev$resid
        if (which=="ranef") lev <- lev$ranef
      } else {
        lev <- get_from_MME(sXaug,which=loctype, B=c("phi", "lambda"))
        if (is.list(lev)) { # depedns on mMatrix_method
          if (which=="resid") lev <- lev$lev_phi
          if (which=="ranef") lev <- lev$lev_lambda 
        } else {
          n_u_h <- length(model$w.ranef)
          if (which=="resid") lev <- lev[-seq_len(n_u_h)]
          if (which=="ranef") lev <- lev[seq_len(n_u_h)] 
          if (which=="both") lev <- list(ranef=lev[seq_len(n_u_h)], resid=lev[-seq_len(n_u_h)])
        }
      }
    }
  }
  if (which=="resid") names(lev) <- names(model$eta)
  return(lev)
} 

# Semantics: hat values: from a projection matrix, vs leverages: final standardizing coefficients
.calc_lev_from_hat <- function(hatvals, sXaug, is_null_phi.Fix, u_h,
                               #object, # all what can be deduced from the 'object'
                               HL, models, need_simple_lambda, muetablob, family, mu, 
                               BinomialDen, w.resid, wranefblob, nobs, ZAL, 
                               psi_M=rep(attr(rand.families,"unique.psi_M"),diff(cum_n_u_h)), 
                               lambda_est, cum_n_u_h, lcrandfamfam=attr(rand.families,"lcrandfamfam"), rand.families, y,  
                               prior.weights, nrand=length(lcrandfamfam), phi_est) {
  ## (HL[2]=0, HL[3]=0): previous hat matrix -> p 
  ## (HL[2]=0, HL[3]=1): notEQL -> tilde(p), (HL[2]=1 && ): full correction -> q 
  ## (HL[2]=1, HL[3]=1): full correction -> q 
  #### contribution from GLM weights
  if (HL[2L]>0L && models[[1L]]=="etaHGLM" 
      && (need_simple_lambda || is_null_phi.Fix) ) { ## LeeN01 HL(.,1) ie the + in 'EQL+'
    ## first the d log hessian / d log lambda or phi corrections
    ### For the d log hessian first the derivatives of GLM weights wrt eta 
    ##################### noter que c'est le coef2 de HL(1,.), but mu,eta may have been updated since coef2 was computed
    dlW_deta <- .calc_dlW_deta(muetablob=muetablob,family=family,
                               BinomialDen=BinomialDen,
                               w.resid=w.resid, families=family)$dlW_deta
    ### we join this with the deriv of log w.ranef wrt v_h
    dlW_deta_or_v <- c(dlW_deta, wranefblob$dlogWran_dv_h )  ## vector with n+'r' elements
    # dlogWran_dv_h is 0 gaussian ranef; d2mudeta2 is 0 for identity link => vector is 0 for LMM
    ## else we continue the computation of the d log hessian term d2 log dens u/ dv dloglambda
    ## where we ignore the contribution of the log Jacobian, log(dth/du), to log dens u since it is not fn of lambda
    ## hence this is d2 log dens th(u)/ dv dloglambda
    if (any(dlW_deta_or_v!=0L)) {
      lev_phi_range <- 1L:nobs
      leve__dlW_deta_or_v <- c(hatvals$resid,hatvals$ranef) * dlW_deta_or_v
      leve__dlW_deta_or_v__ZALI <-  leve__dlW_deta_or_v[lev_phi_range] %*% ZAL +  leve__dlW_deta_or_v[-(lev_phi_range)]
      
      if (need_simple_lambda) {
        neg.d2f_dv_dloglam <- .calc_neg_d2f_dv_dloglam(dlogfthdth=(psi_M - u_h)/lambda_est, ## the d log density of th(u)
                                                       cum_n_u_h, lcrandfamfam, rand.families, u_h)
        dvdloglamMat <-  .calc_dvdloglamMat_new(neg.d2f_dv_dloglam,
                                                sXaug=sXaug) # not d2hdv2_info
        dleve <- as.vector(leve__dlW_deta_or_v__ZALI %*% dvdloglamMat) # (r+n).(r+n)Xr.rXr = r (each element is a sum over r+n terms= a trace)
        hatvals$ranef <- hatvals$ranef - dleve  
      } 
      ## 
      if (is_null_phi.Fix) {
        dh0deta <- ( w.resid *(y-mu)/muetablob$dmudeta ) ## 12/2013 supp BinomialDen (soit Bin -> phi fixe=1, soit BinomialDen=1)
        dvdlogphiMat <- .calc_dvdlogphiMat_new(dh0deta=dh0deta,ZAL=ZAL, 
                                               sXaug=sXaug)  # not d2hdv2_info
        dleve <- as.vector(leve__dlW_deta_or_v__ZALI %*% dvdlogphiMat) # (r+n) . (r+n)Xr . rXn = n (each element is a sum over r+n terms= a trace)
        hatvals$resid <- hatvals$resid - dleve  
      } 
    }
  }
  if (HL[2L]>1) {stop("Need a_i correction in Table 7 of NohL07 ie derivatives of second order correction wrt disp param.")}
  #### contribution from exact likelihood function instead of EQL
  if (HL[3L]!=0 ) {## HL(.,.,1) ie , p_bv(h), not EQL p_bv(q+), LeeNP p89; distinction does not arise for PQL <=> Gaussian ranefs...  
    # lambda
    if (models[[1L]]=="etaHGLM" && need_simple_lambda) ## d h/ d !log! lambda correction     
      hatvals$ranef <- hatvals$ranef + .corr_notEQL_lambda(nrand,cum_n_u_h,lambda_est,lcrandfamfam) 
    # phi hence not poiss,binom:
    if (inherits(family,"family")) {
      if (family$family=="Gamma" && is_null_phi.Fix ) { ## d h/ d !log! phi correction (0 for gauss. resid. error). Not tied to REML
        phiscaled <- phi_est/eval(prior.weights) ## 08/2014 ## bug "*" corrected -> "/" 2015/03/05
        hatvals$resid <- hatvals$resid +  1+2*(log(phiscaled)+digamma(1/phiscaled))/phiscaled ## LNP p. 89 and as in HGLMMM IWLS_Gamma
      }    
    } else { # mv case, list of families
      cum_nobs <- attr(family,"cum_nobs")
      for (mv_it in seq_along(family)) {
        fam <- family[[mv_it]]
        if (fam$family=="Gamma" && is_null_phi.Fix ) { ## d h/ d !log! phi correction (0 for gauss. resid. error). Not tied to REML
          resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
          phiscaled <- phi_est[[mv_it]]/eval(prior.weights[[mv_it]]) 
          hatvals$resid[resp_range] <- hatvals$resid[resp_range] +  1+2*(log(phiscaled)+digamma(1/phiscaled))/phiscaled ## LNP p. 89 and as in HGLMMM IWLS_Gamma
        }    
      }
    }
  }
  hatvals
}

