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
    H_w.resid <- .get_H_w.resid(object=model)
    if (loctype=="fit") {
      lev <- .get_hatvalues_FM(X.Re, augX=model$X.pv, H_w.resid)
    } else { 
      lev <- .get_hatvalues_FM(model$X.pv, augX=model$X.pv, H_w.resid)
    }
  } else {
    #### Check available info:
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
    #### end of checks of available info.
    if (is.null(lev)) {
      sXaug <- model$envir$sXaug # NOT get_matrix() since get_matrix() provides matrices with blocks in Henderson's order,
      if (loctype=="fit") { 
        #### First gte the hat values, strito sensu
        hatvals <- .get_hatvalues_MM(sXaug,X.Re=X.Re, weight_X=attr(sXaug,"weight_X")) 
        # => typically calls get_from_MME(sXaug,which= < "hatval"|"hatval_Z" >, B=c("phi", "lambda"))
        #### Next correct them into standardizing leverages:
        n_u_h <- length(model$w.ranef)
        hatvals <- list(ranef=hatvals[seq_len(n_u_h)], resid=hatvals[-seq_len(n_u_h)])
        
        rand.families <- model$rand.families
        cum_n_u_h <- attr(model$lambda.object$lambda_list,"cum_n_u_h")
        
        v_h_bounds <- .eval_v_h_bounds(cum_n_u_h, rand.families)
        u_h <- .u_h_v_h_from_v_h(model$v_h, rand.families, cum_n_u_h=cum_n_u_h,
                                 lower.v_h=v_h_bounds$lower.v_h,
                                 upper.v_h=v_h_bounds$upper.v_h)
        # We need to identify cases where phi could vary but 
        # where the estimated projection matrix is constrained away from the true one.
        # for Poisson, etc, phi is fixed but not constrained so constr_phi_fit must be FALSE
        if (is.list(model$phi)) { # mv fit
          any_constr_phi_fit <- any(sapply(model$phi, function(v) identical(attr(v,"constr_phi"),TRUE)))
        } else any_constr_phi_fit <- identical(attr(model$phi,"constr_fit"),TRUE)
        # Likewise for lambda, we identify cases where it is fixed by user 
        ismixed <- (model$models[[1]]=="etaHGLM")
        lambdaType <- model$lambda.object$type
        if (( ismixed && any(lambdaType %in% c("fixed","fix_ranCoefs"))) || 
            any_constr_phi_fit) {
          message("Some dispersion parameters were constrained: using leverages for diagnostic purposes is speculative.")
          # nevertheless we try to reproduce the HLfit behaviour in that case and so we evaluate 'fully_constr_lam_fit'
        }
        fully_constr_lam_fit <-  (ismixed && all(lambdaType %in% c("fixed","fix_ranCoefs")))
        if (is.list(model$phi)) { # mv fit
          constr_phi <- all(sapply(model$phi, function(v) identical(attr(v,"constr_phi"),TRUE)))
        } else constr_phi <- identical(attr(model$phi,"constr_phi"),TRUE)
        H_w.resid <- .get_H_w.resid(object=model)
        # Call to .hatvals2std_lev() is a painful hack... but hatvalues() is not the most important extractor in spaMM...
        lev <- .hatvals2std_lev( # tested e.g. by h2s <- hatvalues(adjfit,type="std")[rnge]; more generally hatvalues(.,type="std", force=TRUE)
          hatvals=hatvals, sXaug=sXaug,
          processed=model, # .calc_dlW_deta() may use $y
          anynull_phi.Fix= ( ! constr_phi ) , # should replicate .anyNULL(phi.Fix) in the fit
          u_h=u_h, 
          need_simple_lambda=( ! fully_constr_lam_fit), muetablob=model$muetablob, 
          w.resid=attr(H_w.resid, "w.resid"),
          H_w.resid=H_w.resid,
          wranefblob = .updateW_ranefS(cum_n_u_h=cum_n_u_h, rand.families=model$rand.families, model$lambda, u_h=u_h, model$v_h), 
          nobs=length(model$y), ZAL=get_ZALMatrix(model, force_bind=FALSE), 
          # with default psi_M=rep(attr(rand.families,"unique.psi_M"),diff(cum_n_u_h)), 
          lambda_est=model$lambda.object$lambda_est, cum_n_u_h=cum_n_u_h, #lcrandfamfam=attr(rand.families,"lcrandfamfam"), 
          phi_est=model$phi
        )
        if (which=="resid") lev <- lev$resid
        if (which=="ranef") lev <- lev$ranef
      } else {
        lev <- get_from_MME(sXaug,which=loctype, B=c("phi", "lambda")) # loctype=hatval or hatval_Z
        if (is.list(lev)) { # depends on corr_method
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
.hatvals2std_lev <- function(hatvals, sXaug, anynull_phi.Fix, u_h,
                             processed, # (fit object in post-fit calls)
                             HL=processed$HL, models=processed$models, need_simple_lambda, muetablob, 
                             mu=muetablob$mu, 
                             BinomialDen=processed$BinomialDen, 
                             w.resid,
                             H_w.resid=.BLOB(sXaug)$H_w.resid, # for obsInfo, we need it, generally not for H_w.resid itself but for its attributes
                             wranefblob, nobs, ZAL, 
                             psi_M=rep(attr(rand.families,"unique.psi_M"),diff(cum_n_u_h)), 
                             lambda_est, cum_n_u_h, lcrandfamfam=attr(rand.families,"lcrandfamfam"), 
                             rand.families=processed$rand.families, y=processed$y,  
                             prior.weights=processed$prior.weights, nrand=length(lcrandfamfam), phi_est) {
  ## (HL[2]=0, HL[3]=0): previous hat matrix -> p 
  ## (HL[2]=0, HL[3]=1): notEQL -> tilde(p), (HL[2]=1 && ): full correction -> q 
  ## (HL[2]=1, HL[3]=1): full correction -> q 
  #### contribution from GLM weights
  if (HL[2L]>0L && models[["eta"]]=="etaHGLM" && 
      (need_simple_lambda || anynull_phi.Fix) ) { ## LeeN01 HL(.,1) ie the + in 'EQL+'
    ## first the d log hessian / d log lambda or phi corrections
    ### For the d log hessian first the derivatives of GLM weights wrt eta 
    ##################### noter que c'est le coef2 de HL(1,.), but mu,eta may have been updated since coef2 was computed
    dlW_deta <- .calc_dlW_deta(muetablob=muetablob, 
                               processed=processed, 
                               # using .calc_dlW_deta()" default $family and $families from processed
                               BinomialDen=BinomialDen, 
                               w.resid=w.resid, # w.resid and not H_w.resid # potentially the list with $w_resid element, etc.  
                               # in most cases the result .calc_dlW_deta() simply don't use w.resid to compute only $dlW_deta, although
                               #      * .calc_dlW_deta()'s dlW_deta needs the structured w.resid for some truncated models
                               #      * .calc_dlW_deta() needs w.resid in other calls for $coef1
                               Hratio_factors=attr(H_w.resid,"Hratio_factors") # This is still needed...
                              )$dlW_deta
    ### we join this with the deriv of log w.ranef wrt v_h
    dlW_deta_or_v <- c(dlW_deta, wranefblob$dlogWran_dv_h )  ## vector with n+'r' elements
    # dlogWran_dv_h is 0 gaussian ranef; d2mudeta2 is 0 for identity link => vector is 0 for LMM
    ## else we continue the computation of the d log hessian term d2 log dens u/ dv dloglambda
    ## where we ignore the contribution of the log Jacobian, log(dth/du), to log dens u since it is not fn of lambda
    ## hence this is d2 log dens th(u)/ dv dloglambda
    #
    # dleve_phi is the vector as.vector(leve__dlW_deta_or_v__ZALI %*% invd2hdv2 %*% t(ZAL) %*% diag(dh0deta))
    # which has been computed as 
    # dvdlogphiMat <- .calc_dvdlogphiMat_new(dh0deta=dh0deta,ZAL=ZAL, 
    #                                        sXaug=sXaug)  # not d2hdv2_info
    # dleve <- as.vector(leve__dlW_deta_or_v__ZALI %*% dvdlogphiMat) # (r+n) . (r+n)Xr . rXn = n (each element is a sum over r+n terms= a trace)
    # But it's best computed backward, that is by looking at the transposed expression diag(dh0deta) %*% ZAL %*% invd2hdv2 %*% leve__dlW_deta_or_v__ZALI
    #  and computing successive matrix product of the transpose leftwards so that its rhs ('t_lhs' in the code) is always a vector.
    # Same logic for dleve_lam 
    #
    if (any(dlW_deta_or_v!=0L)) {
      lev_phi_range <- 1L:nobs
      leve__dlW_deta_or_v <- c(hatvals$resid,hatvals$ranef) * dlW_deta_or_v
      leve__dlW_deta_or_v__ZALI <-  drop(leve__dlW_deta_or_v[lev_phi_range] %*% ZAL) +  leve__dlW_deta_or_v[-(lev_phi_range)]
      t_lhs <- .ad_hoc_solve_d2hdv2_Bvec(B=leve__dlW_deta_or_v__ZALI, sXaug=sXaug) # as.vector() result
      if (need_simple_lambda) {
        neg.d2f_dv_dloglam <- .calc_neg_d2f_dv_dloglam(dlogfthdth=(psi_M - u_h)/lambda_est, ## the d log density of th(u)
                                                       cum_n_u_h, lcrandfamfam, rand.families, u_h)
        dleve_lam <- neg.d2f_dv_dloglam * t_lhs # vector
        hatvals$ranef <- hatvals$ranef - dleve_lam  
      } 
      ## 
      if (anynull_phi.Fix) { # => *there is some Gamma or gaussian *GLM* hence not binomial hence BinomialDen (which is in w.resid anyway) is 1. No further BinomialDen here.
        if (processed$how$obsInfo) { 
          if (is.list(w.resid)) { # mv obsInfo case. Can still be GLM but no special code for that case => always as if no GLM
            dh0deta <- .unlist(lapply(muetablob$mv, getElement, name="dlogcLdeta")) # =  drop(.unlist(w.resid$mvlist)*(y-mu)/muetablob$dmudeta) when only GLMs
          } else {
            dh0deta <- muetablob$dlogcLdeta # test code: gaussian(inverse) fit of wafers data in test-devel-LLM
            # it happens that this is a single Gamma or gaussian *GLM* so drop(.unlist(w.resid$mvlist)*(y-mu)/muetablob$dmudeta) would be correct if GLMweights and w.resid where computed.
          }
        } else { # => single *there is some Gamma or gaussian GLM*. Whether obsInfo or not, we still have w.resid vector (but in obsInfo case, we could use muetablob$dlogcLdeta)
          dh0deta <- ( w.resid *(y-mu)/muetablob$dmudeta ) #  muetablob$dlogcLdeta # 
        }
        dleve_phi <- dh0deta * as.vector(ZAL %*% t_lhs) # vector
        hatvals$resid <- hatvals$resid - dleve_phi  
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
    if (is.null(families <- processed$families)) {
      if (processed$family$family=="Gamma" && anynull_phi.Fix ) { ## d h/ d !log! phi correction (0 for gauss. resid. error). Not tied to REML
        phiscaled <- phi_est/eval(prior.weights) ## 08/2014 ## bug "*" corrected -> "/" 2015/03/05
        hatvals$resid <- hatvals$resid +  1+2*(log(phiscaled)+digamma(1/phiscaled))/phiscaled ## LNP p. 89 and as in HGLMMM IWLS_Gamma
      }    
    } else { # mv case, list of families
      cum_nobs <- attr(families,"cum_nobs")
      for (mv_it in seq_along(families)) {
        fam <- families[[mv_it]]
        if (fam$family=="Gamma" && anynull_phi.Fix ) { ## d h/ d !log! phi correction (0 for gauss. resid. error). Not tied to REML
          resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
          phiscaled <- phi_est[[mv_it]]/eval(prior.weights[[mv_it]]) 
          hatvals$resid[resp_range] <- hatvals$resid[resp_range] +  1+2*(log(phiscaled)+digamma(1/phiscaled))/phiscaled ## LNP p. 89 and as in HGLMMM IWLS_Gamma
        }    
      }
    }
  }
  hatvals
}

