HLfit_body <- function(processed, 
         control.HLfit=list(), ## used both by preprocess and HLfit_body
         init.HLfit = list(), ## not from processed: this is affected by HLCor_body -> .canonizeRanPars(ranPars) post .preprocess()ing
         #                       Thus if in a HLCor call we can expect a $corrPars in sp3 code.
         ranFix=list(), ##  possibly trLambda but not necessarily.
         etaFix=list() ## beta, v_h (or even u_h)
) {
  processed$envir$ranFix <- ranFix # for diagnostics reported by div_info() (seek '$ranFixes') [_F I X M E_ rethink] 
  family <- processed$family # may be null in mv case
  ranFix <- .post_process_family(family, ranFix, processed$families) ## assign 'extra' COMPoisson or negbin pars and cleans ranFix of them
  # modifications in environment(family$aic) ARE modifications in environment(processed$family$aic); there is no 'local' copy of that envir. 
  # next line to be called before we extract anything (lambda, ranCoefs... ) from ranFix:
  ranFix <- .canonizeRanPars(ranPars=ranFix,corr_info=NULL, checkComplete = FALSE, rC_transf=.spaMM.data$options$rC_transf)## including full-size lambda
  #data <- processed$data
  verbose <- processed$verbose
  
  predictor <- attr(processed$predictor,"no_offset") 
  prior.weights <- processed$prior.weights
  
  warningList <- list()
  ## when adding verbose elements, remind that these might be lost through corrHLfit -> HLCor cf dotlist$verbose <- verbose[intersect(...]
  ##
  y <- processed$y
  nobs <- length(y) ## before prior.weights is evaluated
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  HL <- processed$HL
  spaMM_tol <- processed$spaMM_tol
  iter_mean_dispFix <- processed$iter_mean_dispFix
  iter_mean_dispVar <- processed$iter_mean_dispVar
  max.iter <- processed$max.iter
  BinomialDen <- processed$BinomialDen
  #X.pv <- processed$AUGI0_ZX$X.pv
  X.Re <- processed$`X.Re` ## should be NULL _xor_ distinct from X.pv
  ### a bit of post processing
  pforpv <- ncol(processed$AUGI0_ZX$X.pv)
  models <- processed$models
  LMMbool <- attr(models,"LMMbool")
  #### Note that HLCor modifies the L matrix (=> ZAL cannot be preprocessed by corrHLfit and must be recomputed each time 
  ZAlist <- processed$ZAlist ## : ZAlist is a list of design matrices 
  nrand <- length(ZAlist)
  cum_n_u_h <- processed$cum_n_u_h
  n_u_h <- cum_n_u_h[nrand+1L] 
  vec_n_u_h <- diff(cum_n_u_h)
  if (models[["eta"]]=="etaHGLM") {
    sparse_precision <- processed$is_spprec
    ranCoefs.Fix <- .getPar(ranFix,"ranCoefs") ## may be NULL
    ranCoefs_blob <- .process_ranCoefs(processed, ranCoefs.Fix, use_tri_CORREL=TRUE) ## UPDATES preexisting object # no augZXy pb here
    LMatrices <- processed$AUGI0_ZX$envir$LMatrices
    # HLCor_body has prefilled $LMatrices for :
    #    for Matern...
    # or for CAR, when HLCor_body has updated LMatrix as fn of rho: seek adjd usage]
    # or was copied from attr(predictor,"LMatrix")
    ### we add the ranCoefs matrices:
    if (any(ranCoefs_blob$is_set)) {
      LMatrices[ranCoefs_blob$is_set] <- ranCoefs_blob$LMatrices[ranCoefs_blob$is_set]
      attr(LMatrices,"is_given_by")[ranCoefs_blob$is_set] <- "ranCoefs"
      # lambda_est initialized from ranCoefs_blob later !
    }
    if (processed$is_spprec) { 
      which_inner_ranCoefs <- which(ranCoefs_blob$isRandomSlope & ( ! ranCoefs_blob$is_set))
      attr(LMatrices,"is_given_by")[which_inner_ranCoefs] <- "inner_ranCoefs"
      .init_AUGI0_ZX_envir_spprec_info(processed)
      .wrap_precisionFactorize_ranCoefs(processed,LMatrices)
      ZAL <- .compute_ZAL(XMatrix=LMatrices, ZAlist=ZAlist,as_matrix=.eval_as_mat_arg(processed), bind.=FALSE)
    } else if ( any((attr(LMatrices,"is_given_by") !="")) ) { # if there are non-trivial L matrices
      ZAL <- .compute_ZAL(XMatrix=LMatrices, ZAlist=ZAlist,as_matrix=.eval_as_mat_arg(processed), 
                          bind.=(.spaMM.data$options$bind_ZAL || HL[1L]=="SEM")) # always TRUE in practice except for devel experiments
      ## ZAL may be modified by other call to .compute_ZAL()   
    } else { 
      ZAL <- processed$AUGI0_ZX$ZAfix ## default ZA 
    } 
  } else if (models[[1]]=="etaGLM") {
    ZAL <- NULL ## 
    u_h <- v_h <- lev_lambda <- numeric(0)
  }
  ranFix$lambda <- .reformat_lambda(ranFix$lambda, nrand, namesTerms=attr(ZAlist,"namesTerms"), full_lambda=TRUE) # necessary to standardize names before next line
  if (any(ranFix$lambda==0,na.rm=TRUE)) stop("lambda cannot be fixed to 0.")
  lam_fix_or_outer_or_NA <- processed$reserve$repNAnrand
  lam_fix_or_outer_or_NA[names(ranFix$lambda)] <- ranFix$lambda # .getPar(ranFix,"lambda") ## should already have length 'nrand' or else be NULL
  ###
  off <- processed$off
  ##################
  valid_inits <- c("fixef","phi","lambda","v_h","rho","nu","Nugget","ARphi","corrPars","ranCoefs")
  unknowns <- setdiff(names(init.HLfit),valid_inits) ## in spaMM 3.0 several names should disappear   
  if (length(unknowns)) {
    if ("beta" %in% unknowns) message("  Use 'fixef' rather than 'beta' in 'init.HLfit'.")
    if ("beta_eta" %in% unknowns) message("  Use 'fixef' rather than 'beta_eta' in 'init.HLfit'.")
    stop(paste0("unhandled elements in 'init.HLfit'. Allowed ones are: ",paste(valid_inits, collapse=",") ))
  }
  ###################
  HLfit_corrPars <- init.HLfit$corrPars 
  corr_est <- .get_cP_stuff(init.HLfit,"rho") ## not spaMM 3.0## init.HLfit[intersect(c("nu","rho","Nugget","ARphi"),names(init.HLfit))]
  if (! is.null(corr_est)) {
    corr_est <- list(rho=corr_est) ## not yet spaMM 3.0 
    ## nothing is done of that; and then of processed$REMLformula anywhere in HLfit_body! (except passed to post-fit fn).
    # corrEstBlob <- .eval_corrEst_args(family=family,rand.families=rand.families,predictor=predictor,data=processed$data,X.Re=X.Re,
    #                                  REMLformula=processed$REMLformula, ranFix=ranFix,
    #                                  Optimizer=control.HLfit$optimizer)
    # corrEst.args <- corrEstBlob$corrEst.args ## but corrEstBlob also has $corrEst.form which will stay there for later use
  }
  if (need_simple_lambda <- need_ranefPars_estim <- models[[1]]=="etaHGLM") {
    need_simple_lambda <- any(is.na(lam_fix_or_outer_or_NA) & ! ranCoefs_blob$is_set)
    need_ranefPars_estim <-  (need_simple_lambda || ! is.null(corr_est))
  } 
  #
  whichadj <- which(attr(ZAlist,"exp_ranef_types")=="adjacency") ## bug presumably corrected here 30/12/2107
  test_inner_estim_rho <- (length(whichadj) && ! is.null(adjd <- attr(LMatrices[[whichadj]],"symsvd")$adjd))
  phi.Fix <- processed$phi.Fix
  if (.anyNULL(phi.Fix)) phi.Fix <- .modify_list(phi.Fix,.getPar(ranFix,"phi")) 
  ## => initial value is preprocessed value. If the latter is NULL, this remains NULL, 
  # except when  RHS was set in final call of outer estimation 
  # (it would be misleading to compute leverages in such a final call)
  # ..modify_list() is necessary for fitmv: phi.Fix <- .getPar(ranFix,"phi") may replace a partially NULL by a full NULL.
  # This however means that (internally-set) ranFix phi component should be a *named list*.
  nothing_to_fit <-  ((! need_ranefPars_estim) && pforpv==0L && (! .anyNULL(phi.Fix)) 
                      && (models[[1]]=="etaHGLM" && (! is.null(etaFix$v_h))) )
  
  ### case where nothing to fit #############################################
  if ( nothing_to_fit ) { 
    if (test_inner_estim_rho) {  
      u.range <- (cum_n_u_h[whichadj]+1L):cum_n_u_h[whichadj+1L]
      adj_rho <- corr_est$rho
      if (is.null(adj_rho)) adj_rho <- .getPar(ranFix,"rho") ## could this occur with test_inner_estim_rho ?
      if (is.null(adj_rho)) adj_rho <- init.HLfit$rho
      if ( ! is.null(adj_rho)) fixed_adjacency_info <- list(whichadj=whichadj, u.range=u.range, coeffs=1/(1-adj_rho*adjd))
    } else fixed_adjacency_info <- NULL
    return(.nothing_to_fit(phi.Fix, off, models, etaFix, rand.families, cum_n_u_h, lcrandfamfam, 
                           lam_fix_or_outer_or_NA, vec_n_u_h, n_u_h, fixed_adjacency_info, ZAL, BinomialDen, processed)) 
  }   ### RETURN !! ## FR->FR but p_bv is not returned. (fixme?: not of class HLfit)
  
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  
  ### Initial values  for lambda, phi and beta from lam_fix_or_outer_or_NA, phi.Fix, or init.HLfit ##### 
  ## Initial estimate for beta  (etaFix acts directly in auglinmodfit)
  beta_eta <- processed$port_env$port_fit_values$fixef # scaled
  if (is.null(beta_eta)) { # then we look at user value and scale it <=> user values must be unscaled.
    if (pforpv) {
      beta_eta <- init.HLfit$fixef
      if ( ! is.null(beta_eta) && ! is.null(attr(processed$AUGI0_ZX$X.pv,"scaled:scale"))) {
        beta_eta <- .scale(beta=beta_eta,X=processed$AUGI0_ZX$X.pv)
      }
    } else beta_eta <- numeric(0L) # don't leave it NULL so that we don't try to get it from inits_by_glm.
  } 
  ## Initial estimate for phi 
  vec_nobs <- processed$vec_nobs
  phi_est <- phi.Fix
  phi_est <- .denullify(phi_est, modifier=processed$port_env$port_fit_values$phi_est, vec_nobs=vec_nobs)
  phi_est <- .denullify(phi_est, modifier=init.HLfit$phi, vec_nobs=vec_nobs)
  if (models[[1]]=="etaHGLM") { 
    ## Initial estimate for u_h, v_h 
    psi_M <- rep(attr(rand.families,"unique.psi_M"),diff(cum_n_u_h))
    v_h <- processed$intervalInfo$init_v_h
    if (is.null(v_h)) v_h <- .initialize_v_h(psi_M=psi_M, etaFix=etaFix, 
                                             init.HLfit=init.HLfit, ## checks init.HLfit$v_h
                                             cum_n_u_h=cum_n_u_h, rand.families=rand.families, port_env=processed$port_env,
                                             v_h_list=processed$reserve$rand_list
                                             )
    u_h <- .u_h_v_h_from_v_h(v_h, rand.families=rand.families, cum_n_u_h=cum_n_u_h,
                             lcrandfamfam=lcrandfamfam, lower.v_h=NULL, upper.v_h=NULL, u_list=processed$envir$u_list)
    ## Initial estimate for lambda in 'compact" form
    init.lambda <- .calc_initial_init_lambda(lam_fix_or_outer_or_NA, nrand, processed, ranCoefs_blob, init.HLfit, ranFix)
  } else init.lambda <- NULL
  ###
  ######### missing Initial estimates for mu, phi, lambda by GLM ####################
  if (  is.null(beta_eta) || .anyNULL(phi_est) || anyNA(init.lambda) ) { 
    inits_by_glm <- .get_inits_by_glm(processed, reset=quote(family$family %in% c("COMPoisson","negbin")) ) # quoted to be applied to each family
    ## : uses processed$y, $BinomialDen, [["control.glm"]]
  }
  # delayedAssign("inits_by_glm", {
  #   .get_inits_by_glm(processed, reset=quote(family$family %in% c("COMPoisson","negbin")) )
  #   }
  # )
  ## Finalize initial values for beta_eta
  if (is.null(beta_eta) ) beta_eta <- inits_by_glm$beta_eta # from .lm.fit or lm.fit using scaled X.pv, hennce must be scaled.
  intervalInfo <- processed$intervalInfo
  if (!is.null(intervalInfo)) {
    parmcol <- attr(intervalInfo$parm,"col")
    beta_eta[parmcol] <- intervalInfo$init ## already appropriately scaled if X.pv has been scaled
  }  
  ## Finalize initial values for phi 
  phi_est <- .denullify(phi_est, modifier=inits_by_glm$phi_est, vec_nobs=vec_nobs)
  ## Finalize initial values for lambda
  if (models[[1]]=="etaHGLM") {
    lambda_est <- .HLfit_finalize_init_lambda(models, init.lambda, processed, ZAL, cum_n_u_h, 
                                              vec_n_u_h, n_u_h, ranCoefs_blob) # mv __FIXME__ ->.calc_fam_correctect_guess() uses total nrand rather than nrand for submodels that contain the ranef
  }
  ###
  ## predictor from initial values
  if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
    eta <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
  } else  eta <- off +drop(processed$AUGI0_ZX$X.pv %*% beta_eta) ## no iteration hence no updating  ## FREQS
  ## conversion to mean of response variable (COUNTS for binomial)
  muetablob <- .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
  mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  w.resid <- .calc_w_resid(muetablob$GLMweights, phi_est) ## 'weinu', must be O(n) in all cases
  conv.phi <- conv.lambda <- conv.corr <- FALSE
  if (models[[1]]=="etaHGLM") {
    if (! is.null(intervalInfo)) intervalInfo$ranFix <- ranFix
    W_ranefS_constant_args <- processed$reserve$W_ranefS_constant_args
    wranefblob <- do.call(".updateW_ranefS",c(W_ranefS_constant_args, list(u_h=u_h,v_h=v_h,lambda=lambda_est))) ## initialization !
    if (pforpv==0L && !is.null(etaFix$v_h)) {
      maxit.mean <- 0 ## used in test near the end...
    } else if ( LMMbool && is.null(intervalInfo) ) {
      maxit.mean <- 1 ## 1 is sufficient for LMM as Hessian does not vary with beta_eta  => quadratic function;
      # ./ E X C E P T for calculation of confidence intervals: at least two intervalSteps are required. Quite visible when 
      # dispersion params are outer-estimated (fitme) in which case there isno outer iteration to compensate a small maxit.mean  
    } else { ## even h maximization in *G*LMMs 
      if ( ! .anyNULL(phi.Fix) && ! need_simple_lambda) { ## allFix hence no true outer iteration 
        maxit.mean <- iter_mean_dispFix 
      } else maxit.mean <- iter_mean_dispVar # If phi.Fix and lam_fix_or_outer_or_NA, the only way to have 'outer' convergence is to have 'inner' convergence
    } 
  } else if (models[[1]]=="etaGLM") {
    if (! is.null(intervalInfo)) {
      intervalInfo$parmcol_X <- parmcol 
    }
    if ( ! .anyNULL(phi.Fix)) { ## 
      maxit.mean <- iter_mean_dispFix 
    } else maxit.mean <- iter_mean_dispVar # If phi.Fix and lam_fix_or_outer_or_NA, the only way to have 'outer' convergence is to have 'inner' convergence
  }
  prev_lik <- -Inf
  conv_logL <- NA
  dvdlogphiMat <- NULL
  dvdloglamMat <- NULL
  penul_lambda <- NULL
  hatvals <- NULL
  ########################################
  ######### Main loop ####################
  ########################################
  iter <- 0L
  if ( ! is.null(vec_nobs)) {
    PHIblob <- list(multiPHI=rep(list(list(prevmsglength=0L)), length(vec_nobs))) # consistently with .calcmultiPHIs() output 
  } else PHIblob <- list(prevmsglength=0L)
  if (HL[1]=="SEM") { ## specif probit
    processed$SEMargs$qr_X <- qr(processed$AUGI0_ZX$X.pv) 
    locarglist <- list(processed=processed, ZAL=ZAL, beta_eta=beta_eta,
                       off=off, corr_est=corr_est, init.lambda=attr(lambda_est,"init.lambda"),
                       lambda.Fix=lam_fix_or_outer_or_NA, LMatrices=LMatrices, verbose=verbose)
    SEMblob <- .probitgemWrap("SEMwrap",arglist=locarglist, pack="probitgem") # eval(as.call(c(quote(SEMwrap),logarglist)))
    beta_eta <- SEMblob$beta_eta
    corr_est["rho"] <- .get_rho_from_SEM_output(SEMblob, lam_fix_or_outer_or_NA)
    if (is.null(SEMblob$glm_lambda)) {
      lambda_est <- lam_fix_or_outer_or_NA
    } else {
      lambda_est <- predict(SEMblob$glm_lambda,type="response")
    }
    u_h <- v_h <- SEMblob$v_h
    logLapp <- SEMblob$logLapp
    attr(logLapp,"seInt") <- SEMblob$seInt ## may be NULL
    ## for summary.HLfit -> beta_cov (_info ?)
    eta <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
    muetablob <- .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
    w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    wranefblob <- do.call(".updateW_ranefS",c(W_ranefS_constant_args, list(u_h=u_h,v_h=v_h,lambda=lambda_est))) ## no fit, likelihood computation
    sqrt.ww <- sqrt(c(w.resid,wranefblob$w.ranef))
    ##
  } else while ( TRUE ) { ##the main loop with steps: new linear predictor, new leverages, new phi, new w.resid, new lambda, new fn(lambda)
    if (models[[1]]=="etaHGLM") {
      if (! is.null(intervalInfo)) { # There is a distinct interval procedure in .calc_etaGLMblob()
        intervalInfo$corr_est <- corr_est ## currently needs to be added ex tempo... apparently only for warn_intervalStep()
        intervalInfo$beta_eta <- beta_eta ## ex tempo: the current estimate of the CI bound
        # to compute a likelihood for theta:=(beta, lambda, phi_pars) The predictions of any ranef must be the ones under theta.
        # This means that if these predictions change with beta the lik does not really be computed 
        # unless, beta,V have reached equilibrium (the present end of the loop).
        # However, in case of mixed model for phi, we would need to *predict* the new phi values as function of given phi parameters
        # and of new deviance residuals implied by the (beta, lambda) and new (beta,V) 
        # (for a phi GLM, no such problem: the phi pars determine the phi values irrespective of the deviance residuals)
        # Hence we can warn about suspiciously high likelihood when:
        intervalInfo$phi_pred_OK <- (intervalInfo$no_phi_pred || conv.phi) 
        #where no_phi_pred is a perhaps too strict condition and || conv.phi possible a too loose one.
      } ## else intervalInfo remains NULL
      H_global_scale <- .calc_H_global_scale(w.resid)
      if ( processed$is_spprec ) { 
        adj_rho <- corr_est$rho
        if (is.null(adj_rho)) adj_rho <- .getPar(ranFix,"rho") ## could this occur with test_inner_estim_rho ?
        if (is.null(adj_rho)) adj_rho <- init.HLfit$rho
        auglinmodblob <- .solve_IRLS_as_spprec(ZAL=ZAL, y=y, n_u_h=cum_n_u_h[nrand+1L],
                                               #H_global_scale=H_global_scale, # de facto ignored through trivial tests
                                               lambda_est=lambda_est, # homoscedastic for CAR! confusing!
                                               muetablob=muetablob,
                                               off=off, maxit.mean=maxit.mean, etaFix=etaFix,
                                               ## supplement for LevenbergM
                                               beta_eta=beta_eta,
                                               ## supplement for ! GLMM
                                               wranefblob=wranefblob, u_h=u_h, v_h=v_h, w.resid=w.resid, phi_est=phi_est,
                                               for_init_z_args=list(nrand=nrand, psi_M=psi_M), 
                                               for_intervals=intervalInfo,
                                               ##
                                               verbose=verbose,
                                               processed=processed,corrPars=list(rho=adj_rho))
      } else {
        auglinmodblob <-  .solve_IRLS_as_ZX(processed$AUGI0_ZX$X.pv, ZAL, y, n_u_h=cum_n_u_h[nrand+1L], 
                                            H_global_scale=H_global_scale, 
                                            lambda_est=lambda_est, muetablob=muetablob,
                                            off=off, maxit.mean=maxit.mean, etaFix=etaFix,
                                            ## supplement for LevenbergM
                                            beta_eta=beta_eta,
                                            ## supplement for ! GLMM
                                            wranefblob=wranefblob, u_h=u_h, v_h=v_h, w.resid=w.resid, phi_est=phi_est,
                                            for_init_z_args=list(nrand=nrand, psi_M=psi_M), 
                                            for_intervals=intervalInfo,
                                            ##
                                            verbose=verbose,
                                            processed=processed)
      }
      ##############################
      beta_eta <- auglinmodblob$beta_eta
      wranefblob <- auglinmodblob$wranefblob # updated only if !LMM
      muetablob <- auglinmodblob$muetablob 
      mu <- muetablob$mu ## testé par accident, necess dans test COMPoisson HLfit...
      w.resid <- auglinmodblob$w.resid # updated only if !LMM
      v_h <- auglinmodblob$v_h
      u_h <- auglinmodblob$u_h
      #eta <- auglinmodblob$eta
      innerj <- auglinmodblob$innerj
    } else if (models[[1]]=="etaGLM") {
      if (pforpv>0L  && maxit.mean) {
        ## resultat nomme' auglinmodblob pour avancer vers simplif du code, mais je ne peux suppser qu'auglinmodblob est toujours calculé
        auglinmodblob <- .calc_etaGLMblob(processed=processed,  mu=mu, eta=muetablob$sane_eta, muetablob=muetablob, old_beta_eta=beta_eta, 
                                          w.resid=w.resid, phi_est=phi_est, off=off, maxit.mean=maxit.mean, 
                                          verbose=verbose, for_intervals=intervalInfo, Xtol_rel=spaMM_tol$Xtol_rel)
        #eta <- auglinmodblob$eta 
        muetablob <- auglinmodblob$muetablob 
        mu <- muetablob$mu
        beta_eta <- auglinmodblob$beta_eta 
        w.resid <- auglinmodblob$w.resid
        innerj <- auglinmodblob$innerj
      } else auglinmodblob <- NULL 
    } # end etaGLM...
    if(inherits(mu,"Matrix")) mu <- drop(mu) ## pb calcul deviance_residual 
    ########## LEVERAGES
    #### base from hat matrix
    if (models[[1]]=="etaHGLM") {
      if (need_ranefPars_estim || .anyNULL(phi.Fix)) {
        if (maxit.mean==0L) {
          stop("(!) Computation of leverages with maxit.mean=0: check that this is meaningful.")
        } # ELSE rWW was updated in the inner loop for betaV
        hatval <- .get_hatvalues_MM(auglinmodblob$sXaug,X.Re=X.Re, auglinmodblob$weight_X)
        hatvals <- list(
          ranef=hatval[seq(n_u_h)],  ## for the ranef residuals (lambda)
          resid=hatval[(n_u_h+1L):(n_u_h+nobs)] ## for the error residuals (phi)
        )
      }
    } else if (.anyNULL(phi.Fix)) { ## phi estim in GLM fitted by ML. 
      hatvals <- list(resid=.get_hatvalues_FM(X.Re, augX=processed$AUGI0_ZX$X.pv, w.resid))
    }  
    ## (HL[2]=0, HL[3]=0): previous hat matrix -> p 
    ## (HL[2]=0, HL[3]=1): notEQL -> tilde(p), (HL[2]=1 && ): full correction -> q 
    ## (HL[2]=1, HL[3]=1): full correction -> q 
    if ( ! is.null(processed$vec_nobs)) {
      leverages <- .calc_lev_from_hat(hatvals, sXaug=auglinmodblob$sXaug, 
                                      is_null_phi.Fix=.anyNULL(phi.Fix), u_h=u_h, 
                                      #object=processed,
                                      HL=HL, models=models, need_simple_lambda=need_simple_lambda, muetablob=muetablob, 
                                      family=processed$families, 
                                      mu=mu, 
                                      BinomialDen=BinomialDen, w.resid=w.resid, wranefblob=wranefblob, nobs=nobs, ZAL=ZAL, 
                                      psi_M=psi_M, lambda_est=lambda_est, cum_n_u_h=cum_n_u_h, lcrandfamfam=lcrandfamfam, rand.families=rand.families, 
                                      y=y, prior.weights=processed$prior.weights, nrand=length(lcrandfamfam), phi_est=phi_est)
    } else leverages <- .calc_lev_from_hat(hatvals, sXaug=auglinmodblob$sXaug, 
                                           is_null_phi.Fix=.anyNULL(phi.Fix), u_h=u_h, 
                                           #object=processed,
                                           HL=HL, models=models, need_simple_lambda=need_simple_lambda, muetablob=muetablob, family=family, mu=mu, 
                                           BinomialDen=BinomialDen, w.resid=w.resid, wranefblob=wranefblob, nobs=nobs, ZAL=ZAL, 
                                           psi_M=psi_M, lambda_est=lambda_est, cum_n_u_h=cum_n_u_h, lcrandfamfam=lcrandfamfam, rand.families=rand.families, 
                                           y=y, prior.weights=processed$prior.weights, nrand=length(lcrandfamfam), phi_est=phi_est)    
    lev_phi <- leverages$resid
    lev_lambda <- leverages$ranef
    
    ######### Dispersion Estimates for phi #####################
    if (.anyNULL(phi.Fix)) { ## if phi is estimated (phi.Fix set to 1 for Poisson, Binomial)
      ## leverages have been computed before the  inner loop, which did not change the design matrices 
      lev_phi[lev_phi>1-1e-8] <- 1-1e-8
      ## updated residuals from updated mu must be used (LeeNP p.161) 
      ## Once in Gamma GLM, y=25.75563, y-mu=-5.996906e-08, yielded negative dev.resid
      
      if ( ! is.null(vec_nobs)) { # fitmv case
        # models sharing random effects have no reason to share phi parameters.
        PHIblob <- .calcmultiPHIs(processed=processed, y=y,mu=mu, wt= eval(prior.weights), lev_phi=lev_phi, phimodel=models[["phi"]], 
                                  #
                                  verbose=verbose, control.HLfit=control.HLfit, phi.Fix=phi.Fix,
                                  # hglm:
                                  iter=iter, prev_multiPHIblob=PHIblob)
        next_phi_est <- PHIblob$next_phi_est # value of *phi* (not phi_i:= phi/prior.weights as pw are used in GLMweights, not here)  
        # phi_est <- unlist(phi_est) # for the comparison to unlist(next_phi_est). Assumes that phi_est is not used later... 
        #                                                         BUT is it (-> port_fit_value -> comparison by .modify_list)
        for (mv_it in seq_along(PHIblob)) {
          if (models[["phi"]][mv_it]=="phiHGLM") {
          } else {
            if (! is.null(locw <- PHIblob[[mv_it]]$glm_phi$warnmess)) warningList$innerPhiGLM <- locw
          }
        }
      } else {
        # to obtain the phi estimate given by summary.glm(), one must use
        #  dev.res <- wt * ((y-mu)/family$mu.eta(eta))^2 ## wt * EQL residuals
        # but the logLik differs from that given by logLik.glm(). See Details in ?HLfit
        PHIblob <- .calcPHI(processed,
                            dev.res= family$dev.resids(y,mu,wt= eval(prior.weights)), # eta-family !
                            #: times pw to be an estimate of same phi across level of response
                            # but not same phi as when there is no pw !
                            # double pw => double phi_est so that phi_est_i :=phi_est/pw_i is unchanged
                            lev_phi=lev_phi,
                            phimodel=models[["phi"]],
                            verbose=verbose, 
                            # glm:
                            control.HLfit=control.HLfit,
                            # hglm:
                            iter=iter, prev_PHIblob=PHIblob)
        next_phi_est <- PHIblob$next_phi_est # value of *phi* (not phi_i:= phi/prior.weights as pw are used in GLMweights, not here)  
        if (models[["phi"]]=="phiHGLM") {
        } else {
          if (! is.null(locw <- PHIblob$glm_phi$warnmess)) warningList$innerPhiGLM <- locw
        }
      }
      # phi_est may be (a single value|a long vector| a list of such values, for mv).
      if (is.list(phi_est)) {
        conv.phi <- TRUE
        for (mv_it in seq_along(phi_est)) {
          conv.phi <- all(abs(next_phi_est[[mv_it]]-phi_est[[mv_it]]) < 
                            spaMM_tol$Xtol_rel*phi_est[[mv_it]] + processed$spaMM_tol$Xtol_abs) 
          if ( ! conv.phi) break
        }
      } else conv.phi <- all(abs(next_phi_est-phi_est) < spaMM_tol$Xtol_rel*phi_est + processed$spaMM_tol$Xtol_abs) ## 'weak convergence'... 
      # :           # may substract scalar (initial value) to vector (predictions from phi model) 
    } else {conv.phi <- TRUE} ## there is a phi.Fix, already -> phi_est
    ###############################################################
    ######### Dispersion Estimates for lambda #####################
    ###############################################################
    
    if (need_ranefPars_estim) { ## lambda must be estimated 
      levlam_bound <- 1 - .spaMM.data$options$regul_lev_lambda
      if (any(abs(lev_lambda) > levlam_bound)) { ## abs... not commented when written...
        lev_lambda[lev_lambda>levlam_bound] <- levlam_bound
        warningList$leveLam1 <- TRUE
      }      
      ################## ranefEstargs mustcontain arguments for makeCoveEst1 => its names are constrained
      ranefEstargs <- list(u_h=u_h,ZAlist=processed$ZAlist,cum_n_u_h=cum_n_u_h,
                           prev_LMatrices=LMatrices,
                           #lcrandfamfam=lcrandfamfam,
                           processed=processed,
                           init_ranCoefs=init.HLfit$ranCoefs,
                           prev_lambda_est=lambda_est,
                           w.resid=w.resid)
      if (any(ranCoefs_blob$isRandomSlope)) { ## if random-slope model
        ranefEstargs <- c(ranefEstargs,list(phi_est=phi_est,
                                            as_matrix=( ! inherits(ZAL,"Matrix")),v_h=v_h))
        H_global_scale <- .calc_H_global_scale(w.resid)
        ## MakeCovEst defines à local ZAL and the eta,mu, w.resid must generally be recomputed locally for this ZAL
        ranefEstargs$MakeCovEst_pars_not_ZAL_or_lambda <- list(X.pv=processed$AUGI0_ZX$X.pv, y=y, n_u_h=cum_n_u_h[nrand+1L], 
                                                               H_global_scale=H_global_scale, 
                                                               muetablob=NULL,
                                                               off=off, maxit.mean=maxit.mean, etaFix=etaFix,
                                                               ## supplement for LevenbergM
                                                               beta_eta=beta_eta,
                                                               ## supplement for ! GLMM
                                                               u_h=u_h, v_h=v_h, w.resid=NULL, phi_est=phi_est,
                                                               for_init_z_args=list(nrand=nrand, psi_M=psi_M), 
                                                               for_intervals=intervalInfo,
                                                               verbose=c(TRACE=FALSE,trace=FALSE),
                                                               ##
                                                               processed=processed)
      }
      calcRanefPars_blob <- .calcRanefPars(HLfit_corrPars=HLfit_corrPars,
                                           lev_lambda=lev_lambda,
                                           ranefEstargs=ranefEstargs,
                                           ranCoefs_blob=ranCoefs_blob,
                                           lam_fix_or_outer_or_NA=lam_fix_or_outer_or_NA,
                                           rand.families=rand.families,
                                           psi_M=psi_M,
                                           verbose=verbose,
                                           iter=iter,
                                           control=processed[["control.glm"]],
                                           maxLambda=processed$maxLambda
      )
      HLfit_corrPars <- calcRanefPars_blob$HLfit_corrPars
      next_info_for_conv_rC <- calcRanefPars_blob$HLfit_corrPars$info_for_conv_rC ## vector
      next_lambda_est <- calcRanefPars_blob$next_lambda_est  ## T H I S is what affects next beta,v estimates when rho is inner estimated in CAR.
      ########################################
      # => variation of log(u^2/lambda) = simplified likRanU convergence  (from 1.9.31) (see below use of verbose["trace"] to print diagnostics)
      # conv lambda code depending on gaussian_u_ranges had been ineffective for a long time and commented-out in v. 3.6.38 
      # then fully removed in v.3.6.39. Also ***conv.lambda no longer depends on iter>1L***.
      if ( .anyNULL(phi.Fix) && max(.unlist(phi_est))<1e-4) { ## default $Xtol_abs a bit lax for maximizing lik of quasi deterministic response (phi and lambda -> 0)
        loc_Xtol_abs <- min(spaMM_tol$Xtol_abs, 1e-7 )
      } else loc_Xtol_abs <- spaMM_tol$Xtol_abs
      conv.lambda <- all( abs(next_lambda_est-lambda_est) < spaMM_tol$Xtol_rel *lambda_est + loc_Xtol_abs)  
      if (!is.null(next_info_for_conv_rC)) {
        if (conv.corr <- (iter>1L)) {
          dlogL_rC <- abs(info_for_conv_rC$obj-next_info_for_conv_rC$obj) 
          if ( ! (conv.corr <- (dlogL_rC<1e-07))) { # strict bu sufficient condition on logL convergence, else:
            # the good conve crit on the parameters is difficult to find. Cf comments on HLfit6 in .makeCovEst1()
            rel_dcov <- abs(info_for_conv_rC$ranCoefs-next_info_for_conv_rC$ranCoefs)/(1+abs(info_for_conv_rC$ranCoefs))
            rel_dcov <- abs(info_for_conv_rC$ranCoefs-next_info_for_conv_rC$ranCoefs)/(1+abs(info_for_conv_rC$ranCoefs))
            conv.corr <- (max(rel_dcov) < spaMM_tol$Xtol_rel ) 
            #print(rel_dcov)
          }  
        } 
        info_for_conv_rC <- next_info_for_conv_rC
      } else conv.corr <- TRUE
    } else { conv.lambda <- conv.corr <- TRUE } ## end if need_ranefPars_estim else...
    #
    iter <- iter+1L ## here first from 0 to 1
    ###### convergence: 
    if ( conv.phi && conv.lambda && conv.corr) {
      break 
    } else if (iter>=max.iter) { ## did not converge...
      break 
    } else { ## update and continue
      if (.anyNULL(phi.Fix)) phi_est <- next_phi_est
      if (models[[1]]=="etaHGLM") {
        if ( ! is.null(corr_est) ) {
          corr_est <- list(rho=.getPar(HLfit_corrPars,"rho")) ## not spaMM 3.0 : single inner-optimized rho
          #LMatrix est constant!= decomp$u
        }
        if (need_ranefPars_estim) { ## next_lambda_est/ next ranCoefs/ next rho adjacency available:
          if (any(ranCoefs_blob$isRandomSlope)) {
            LMatrices <- calcRanefPars_blob$next_LMatrices ## keep L factor of corr mats for all models 
            if (processed$is_spprec) {
              for (rt in which_inner_ranCoefs) {
                processed$AUGI0_ZX$envir$precisionFactorList[[rt]] <- 
                  .precisionFactorize(latentL_blob=attr(LMatrices[[rt]],"latentL_blob"),
                                      rt=rt, longsize=ncol(LMatrices[[rt]]), processed=processed,
                                      cov_info_mat=processed$corr_info$cov_info_mats[[rt]])
              }
            }
            ZAL <- .compute_ZAL(XMatrix=LMatrices, ZAlist=processed$ZAlist,as_matrix=( ! inherits(ZAL,"Matrix")),
                                bind.= ! processed$is_spprec) 
            if ( ! LMMbool ) {
              ## ZAL is modified hence wranefblob must be modified (below) but also eta-> mu->GLMweights
              ## .makeCovEst1 may have reestimated beta but we do not take this into account nor any resulting change in the 'blobs'
              eta <- off + drop(processed$AUGI0_ZX$X.pv %*% beta_eta) + drop(ZAL %id*% v_h) 
              muetablob <- .muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
            }
          }
          ###################################################################################################
          if ( FALSE && NCOL(processed$X_lambda)==1L) {## FR->FR this almost worked ./.
            ## ./. mais voir premier exemple dans ?spaMMqui boucle...comprends pas                ..... currently not $X_lambda 
            if (iter>2L) {
              loglam0 <- log(antepenul_lambda)  
              loglam1 <- log(penul_lambda)  
              loglam2 <- log(lambda_est[1L])  
              loglam3 <- log(next_lambda_est[1L])  
              slope1 <- (loglam2-loglam1)/(loglam1-loglam0)
              slope2 <- (loglam3-loglam2)/(loglam2-loglam1) ## FR->FR need to handle 0 denoms
              if ((abs(log(abs(slope2)))>log(1.05)) ##  <0.95 || abs(slope2)>1.05) ## estimate will not explode
                  && (logarg <- slope1/slope2)>0   
                  && abs(log(logarg))<0.1 # we are in geom phase
              ) {               
                geom_est_loglam <- loglam2+(loglam3-loglam2)/(1-slope2)
                print(c(iter,loglam0,loglam1,loglam2,loglam3,geom_est_loglam))
                if ( ! is.na(geom_est_loglam) && ! is.infinite(geom_est_loglam) ) {
                  geom_est_lam <- exp(geom_est_loglam)
                  next_lambda_est <- rep(geom_est_lam,length(next_lambda_est))
                } 
              } else {
                print(c(iter,loglam0,loglam1,loglam2,loglam3))
              }
            }
            antepenul_lambda <- penul_lambda
            penul_lambda <- lambda_est[1L]
          }
          ###################################################################################################
          # UPDATE:
          ##      in particular, also in the adjacency case if rho was updated but not the lambda param.
          lambda_est <- next_lambda_est
          wranefblob <- do.call(".updateW_ranefS",c(W_ranefS_constant_args, list(u_h=u_h,v_h=v_h,lambda=lambda_est))) ## bc lambda was modified
        } 
      }
      if ( .anyNULL(phi.Fix) || ( models[[1]]=="etaHGLM" && any(ranCoefs_blob$isRandomSlope) && ! LMMbool) ) { ## phi or (ZAL -> mu) modified
        w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est) ## bc phi was updated. 'weinu', must be O(n) in all cases 
      }
      ## conv_logL either used to break the loop, Xor required only in last two iters for diagnostics 
      if (processed$break_conv_logL || verbose["trace"] ) {
        next_lik <- .calc_APHLs_from_ZX(auglinmodblob=auglinmodblob, which="p_v", processed=processed)$p_v
        # this does not depend on the latest ranPars estimates, since sXaug was not updated afterwards... 
        conv_logL <- abs(next_lik - prev_lik)/(0.1 + abs(next_lik)) < 1e-8 # ~ glm.fit convergence
        if (processed$break_conv_logL && conv_logL) break 
        prev_lik <- next_lik
      } else conv_logL <- NA
      ##
      if (verbose["trace"]) {
        print(paste("iteration ",iter,
                    #"; convergence criteria for phi, lambda, corr pars , conv_lambda_vs_u, conv_rel_lambda: ",
                    "; convergence criteria for phi, lambda, corr pars: ",
                    paste0(c( conv.phi , conv.lambda, conv.corr), #,conv_lambda_vs_u, conv_rel_lambda),
                           collapse = " "),
                    "; logL & conv_logL:", next_lik, conv_logL))
        if (need_simple_lambda) { 
          #print(range(logrel_crit))
          #print(range(reldlam_crit))
        }
      } 
    } 
    ##### end convergence block
  } ## end main loop while ( TRUE )
  ########################################
  ######### END main loop ################
  ########################################
  if (verbose["trace"]) {
    if (iter==max.iter) {
      mess <- paste("(beta,v)/lambda/phi iterations failed to converge in",max.iter,"iterations")
      message(mess)
    } else {
      message(paste("(beta,v)/lambda/phi iterations in HLfit() converged in",iter,"iterations"))
    }
  }
  #
  if (HL[1]=="SEM") {
    APHLs <- list(logLapp=logLapp) ## keeps attributes
    APHLs$clik <- .calc_clik(mu,phi_est,processed) ## useful for .get_info_crits()
  } else {
    if (models[[1]]=="etaHGLM") {
      if (identical(processed$return_only,"p_vAPHLs")) {
        whichAPHLs <- "p_v"
      } else if (identical(processed$return_only,"p_bvAPHLs")) {
        whichAPHLs <- "p_bv" ## return value may still include p_v if it is used to compute p_bv
      } else whichAPHLs <- c("p_v","p_bv")
      APHLs <- .calc_APHLs_from_ZX(auglinmodblob,which=whichAPHLs,processed)
    } else { 
      APHLs <- .calc_APHLs_from_ZX(auglinmodblob=NULL,processed, which="p_v",
                                   sXaug=NULL, phi_est, lambda_est=NULL, dvdu=NULL, u_h=NULL, muetablob=muetablob ) 
    }
  }
  ######################### potential R E T U R N here: cases without p_bv
  if ( identical(processed$return_only,"p_vAPHLs")) {
    # Following comment no longer clear, but this may have referred to the non-existence of 'processed'  in optimthrousmooth code
    # a bit of ugly coding, but optimthroughsmooth calls HLCor, not HLCor.obj, thus it cannot directly control return_only. So either leave as is, or move the test to HLCor, or modify optimthroughsmooth to call HLCor.obj  
    if (HL[1]=="SEM") { # lambda used for smoothing.
      res <- list(APHLs=APHLs,lambda=SEMblob$lambda) 
    } else {
      res <- list(APHLs=APHLs) 
      if ( ! is.null(oldp_v <- processed$port_env$objective)) {
        .update_port_fit_values(old_obj=oldp_v,new_obj=APHLs$p_v, 
                                port_fit_values=list(fixef=beta_eta,v_h=v_h,phi_est=phi_est), 
                                models=models, processed=processed, control.HLfit=control.HLfit,
                                lambda_est=lambda_est, 
                                PHIblob=PHIblob)
      } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
        assign("objective",APHLs$p_v,envir=processed$port_env)
      }
    }
    return(res)   ########################   R E T U R N
  } 
  ## ELSE continue: make sur p_bv is included
  if (HL[1] != "SEM") {
    if (models[[1]]=="etaHGLM") {
      ## nothing to do as APHLs already contain p_bv.
      ## cf notes 19/08/2016 pour calcul APHLs et IC's for phiHGLM 
    } else { ## G L M
      ## ML: X.Re non NULL mais ncol(X.Re)=0
      X.REML <- X.Re
      if (is.null(X.Re)) {X.REML <- processed$AUGI0_ZX$X.pv} ## REML standard
      if ( ncol(X.REML) ) { ## REML standard || REML non standard
        Md2hdb2 <- .ZtWZwrapper(X.REML,w.resid)
        ladb <- .LogAbsDetWrap(Md2hdb2,logfac=-log(2*pi))
      } else ladb <- 0 ## fit ML standard : ncol(X.Re)=0, p_bv=p_v hence d2hdpbv reduces to d2hdv2
      APHLs <- list(p_v=APHLs$clik, p_bv=APHLs$clik - ladb/2) ## 07/2016 for inference about phi in GLMs 
    }
  }
  ######################### potential R E T U R N here: with p_bv
  if ( identical(processed$return_only,"p_bvAPHLs")) {
    if ( ! is.null(oldp_bv <- processed$port_env$objective)) {
      .update_port_fit_values(old_obj=oldp_bv,new_obj=APHLs$p_bv, 
                              port_fit_values=list(fixef=beta_eta,v_h=v_h,phi_est=phi_est), 
                              models=models, processed=processed, control.HLfit=control.HLfit,
                              lambda_est=lambda_est, 
                              PHIblob=PHIblob)
    } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
      assign("objective",APHLs$p_bv,envir=processed$port_env)
    }
    res <- list(APHLs=APHLs)
    return(res)    ########################   R E T U R N
  } 
  # beta_cov code removed from here in v1.9.24
  ######################
  ######################
  ######################
  ##### LAMBDA and other RANEF PARS
  ## (1) we count (most) inner-estimated ranef pars (but we add some later...)
  if (need_ranefPars_estim) {
    bloc_lambda_args <- list(models=models, 
                             processed=processed, lam_fix_or_outer_or_NA=lam_fix_or_outer_or_NA, 
                             cum_n_u_h=cum_n_u_h, next_LMatrices=LMatrices)
    if (HL[1]=="SEM") {
      bloc_lambda_args$SEMblob <- SEMblob
    } else {
      bloc_lambda_args$calcRanefPars_blob <- calcRanefPars_blob
      bloc_lambda_args$lev_lambda <- lev_lambda
    }
    #
    process_resglm_blob <- do.call(".bloc_lambda",bloc_lambda_args)
    # for a 3-par ranCoef estimated internally, counted in p_lambda, 
    # 2 pars are counted from the coefficients_lambdaS, and one by p_corr added to p_lambda
    # For outer-estimated raneCoefs, see     p_lambda <- p_lambda + p_ranCoefs ... below
    coefficients_lambdaS <- process_resglm_blob$coefficients_lambdaS # list
    p_lambda <- length(.unlist(coefficients_lambdaS)) - length( which(attr(init.lambda,"type") %in% c("fixed", "outer_hyper", "fix_hyper"))) 
    p_adjd <- numeric(length(coefficients_lambdaS) )
    for (it in seq_len(length(coefficients_lambdaS))) p_adjd[it] <- length( which(names(coefficients_lambdaS[[it]])=="adjd")) 
    p_adjd <- sum(p_adjd) 
    p_lambda <- p_lambda-p_adjd
    var_ranCoefs <- with(ranCoefs_blob, (isRandomSlope & ! is_set))
    if (any(var_ranCoefs)) {
      p_corr <- vector("list", length(LMatrices))
      for (it in which(var_ranCoefs)) {
        dimL <- NROW(attr(LMatrices[[it]],"latentL_blob")$compactcovmat)
        if (dimL==0L) stop("'latentL_blob' attribute missing to LMatrices[[it]].")
        p_corr[[it]] <- (dimL-1)*dimL/2
      }
      p_corr <- sum(.unlist(p_corr))
      p_lambda <- p_lambda + p_corr
    }
  } else {
    p_lambda <- p_adjd <- 0
    process_resglm_blob <- list(lambda_pred_list=as.list(rep(NA,nrand)))
  }
  ## (2) we count outer estimated ones, except hy_lam ones
  if (! is.null(preproFix <- processed$lambda.Fix)) {
    #p_lambda <- p_lambda + length(which(is.na(preproFix[ ! is.na(lambda.Fix)])))
    #p_lambda <- p_lambda - length(which(attr(init.lambda,"type")=="outer_hyper"))
    p_lambda <- p_lambda + length(which(attr(init.lambda,"type")=="outer")) 
  }##### PHI: 
  # if (models[["phi"]]=="phiHGLM") {
  #   ## nothing to be done since we have a full fitme'd object. We could same some time by hacking _fitme_ 
  #   ##  so that it does not perform its final call but still returns the phi_est (ie its mu). Not urgent.
  # } # else nothing to do here, as the phi_GLM is built one request from summary() if missing 
  ######################################
  ## BUILD full RETURN VALUE
  ######################################
  #
  ###################
  ## LIKELIHOODS
  ###################
  res <- list(APHLs=APHLs)
  ###################
  ## FIXEF, ETA, ... 
  ###################
  res$X.pv <- processed$AUGI0_ZX$X.pv
  if ( ! is.null(attr(res$X.pv,"scaled:scale"))) {
    beta_eta <- .unscale(beta=beta_eta, X=res$X.pv)
    res$X.pv <- .unscale(res$X.pv) ## usefully not in an environment
  } 
  ## Assuming beta_eta is a vector, not a matrix
  #  if ( ! is.null(namesOri <- attr(processed$AUGI0_ZX$X.pv,"namesOri"))) { ## including NA's names (and etaFix$beta names) # Condition always TRUE ?
  namesOri <- attr(processed$AUGI0_ZX$X.pv,"namesOri")
  nc <- length(namesOri)
  beta_etaOri <- rep(NA,nc)
  names(beta_etaOri) <- namesOri
  beta_etaOri[names(beta_eta)] <- beta_eta ## keeps the original NA's
  beta_etaOri[names(etaFix$beta)] <- etaFix$beta  ## no longer in X.pv 2015/03
  res$fixef <- beta_etaOri ## fixme I should keep out the fixed ones for summary ? newetaFix code assumes the opposite
  # } else {
  #   names(beta_eta) <- colnames(processed$AUGI0_ZX$X.pv)
  #   res$fixef <- beta_eta 
  # } 
  if (identical(processed$return_only,"confint_bound")) {
    return(res)    ########################   R E T U R N fixef + APHLs
  }
  res$eta <- muetablob$sane_eta ## convenient for defining starting values... and also sometimes used by predict() 
  etanames <- .unlist(attr(processed$data,"validrownames"))
  if (is.null(etanames)) { # must be all case except main (ie not residModel) fit of a fitmv (=> code to make sure of post-fit residProcessed$data in fitmv())
    names(res$eta) <- rownames(processed$data)
  } else names(res$eta) <- etanames
  res$muetablob <- list(mu=muetablob$mu,dmudeta=muetablob$dmudeta) # for get_logdispObject, added 11/2016
  ###################
  ## DATA
  ###################
  res$data <- processed$data
  if (is.null(family) || family$family=="binomial") { # null for mv case
    res$BinomialDen <- BinomialDen # we could put it in all cases...
  }
  res$y <- y ## counts for Pois/bin
  res$prior.weights <- prior.weights # with attrs, and possibly a call
  res$prior.weights[] <- eval(prior.weights) ## see Gamma()$simulate # eval the call but keep the attributes.
  ###################
  ## MODEL info
  ###################
  res$family <- family
  res$families <- processed$families
  res$ranFix <- ranFix ## currently as a uniform template consistent with projected changes ; except that lamFix, phiFix info is now in lambda.object, etc
  ## If an outer optimizer has been called,
  #  "fix" and "outer" parameters are given these types by .get_refit_args() after the optimization call, 
  # then HLfit called again and we reach this point.
  #  This means ranFix gets its type from there *if* properly retained by .canonizeRanPars() 
  #  Then we add inner-optimized parameters, with "var" type added by .get_CorrEst_and_RanFix()
  if ( ! is.null(corr_est) && ! is.null(init.HLfit$corrPars)) corr_est <- list(corrPars=relist(corr_est$rho,init.HLfit$corrPars)) ## not yet spaMM 3.0
  # Canonical, and inherits all info about outer-optimized corrPars through HLfit's ranFix argument:
  res$CorrEst_and_RanFix <- .get_CorrEst_and_RanFix(ranFix, corr_est) # corr_est parameters are inner-estimated and of type "var"
  #
  if ( ! is.null(res$CorrEst_and_RanFix$corrPars)) {
    ## corrPars types are "outer", "fix", or "var" (for inner optimized rho)
    p_corrPars <- length(which(unlist(attr(res$CorrEst_and_RanFix,"type")$corrPars, use.names = FALSE) == "outer")) 
    p_var_corrPars <- length(which(unlist(attr(res$CorrEst_and_RanFix,"type")$corrPars, use.names = FALSE) == "var")) 
    if (p_var_corrPars != p_adjd) {
      warning("Suspect mismatch in df count (p_var_corrPars != p_adjd). AIC computations may be incorrect.")
    } else {
      p_corrPars <- p_corrPars + p_var_corrPars # consistent with df count in the case of outer optimization
    }
    # we add outer ranCoefs to p_lambda, same as inner ranCoefs added to p_lambda previously
    p_lambda <- p_lambda + .dfs_ranCoefs(types=attr(res$CorrEst_and_RanFix,"type"))
    # at this point, the _type_ of the corrPars controlled by hyperparameters have been removed, and that for lambda is "fix" (!)
    hy_var <- which(unlist(attr(res$CorrEst_and_RanFix,"type")$hyper, use.names = FALSE) != "fix")
    if (length(hy_var)) {
      p_corrPars <- p_corrPars + length(grep("hy_trK", names(hy_var)))
      p_lambda <- p_lambda + length(grep("hy_trL", names(hy_var)))
    }
    res$corrPars <- structure(res$CorrEst_and_RanFix$corrPars, # ## subset of the above: F I X M E (?) redundancy but convenient when examining fits
                              type=attr(res$CorrEst_and_RanFix,"type")$corrPars,
                              message='Use get_ranPars(.,which="corrPars") to extract "corrPars" cleanly from fit object.')
  } else p_corrPars <- 0  
  res$dfs <- list(pforpv=pforpv, p_lambda=p_lambda, p_fixef_phi=processed$p_fixef_phi, p_corrPars=p_corrPars)
  res$models <- structure(models, LMMbool=LMMbool)
  res$main_terms_info <- processed$main_terms_info ## used by predict
  res$predictor <- processed$predictor ##  all post fitting functions expect PROCESSED predictor
  res$vec_nobs <- processed$vec_nobs ## non-null for fitmv
  #
  if (models[[1]] == "etaHGLM") res$ZAlist <- processed$ZAlist ## needed for prediction variance
  res$REMLformula <- processed$REMLformula  # only for .REMLmess()... but it's still a simple way to pass the info. Perhaps put it elsewhere in res?
  ###################
  ## OBJECTIVE and ALGORITHMs
  ###################
  res$HL <- HL ## info on fitting objective
  res$how <- list(spaMM.version=packageVersion("spaMM"),
                  MME_method=.get_MME_method(auglinmodblob, HL),
                  switches=c(augZXy_cond=processed$augZXy_cond,
                             use_spprec_QR=.spaMM.data$options$use_spprec_QR)
  )
  # !!! res$MME_method used later in this fn !!!
  res$MME_method <- structure(res$how$MME_method,
                              message="Please use how(<fit object>)[['MME_method']] to extract this information cleanly.")                 
  res$spaMM.version <- structure(res$how$spaMM.version, ## this is NOT a string and comparison with a string is suitably def'ed (as detailed in ?package_version)
                                 message="Please use how(<fit object>)[['spaMM.version']] to extract this information cleanly.")                 
  if (HL[1]=="SEM") res$SEM_info <- SEMblob$SEM_info ## info
  ###################
  ## FITTED VALUES
  ###################
  res$fv <- .mu_U2fv(mu, # For truncated fams, input mu of latent untruncated  variable; output fv is different but mu is kept as an attribute
                     # For binomial: input mu is count, output fv is proba
                     BinomialDen, muetablob=muetablob, processed)
  ###################
  ## LEVERAGES and REML (ie either phi OR lambda was estimated)
  ###################
  if (HL[1]!="SEM") { ## both lev_phi and deviance_residual missing otherwise
    if (.anyNULL(phi.Fix) || need_simple_lambda) { ## in either case all leverages are computed and it makes sense to consider the residuals
      res$lev_phi <- lev_phi
      dev_res_blob <- .std_dev_resids(res, phi_est=phi_est, lev_phi=lev_phi)
      res$std_dev_res <- sign(y-mu) * dev_res_blob$std_dev_res
      dev_res <- dev_res_blob$dev_res # needed below
    }
    if (need_simple_lambda) res$lev_lambda <- lev_lambda
  }  
  res$distinctX.Re <- X.Re ## NULL if not distinct from X.pv
  ###################
  ## ALL other LAMBDA returns
  ###################
  res$rand.families <- rand.families 
  ##
  res$QRmethod <- processed$QRmethod
  #
  ## $w.ranef and $w.resid not doc'ed, as there is no mention of the augmented mode lin the doc.
  if (is.list(w.resid)) { ## truncated 'family', or 'families' with some truncated one(s), but not all mv cases
    res$w.resid <- w.resid$w_resid ## useful for .get_info_crits() and get_LSmatrix
  } else res$w.resid <- w.resid ## useful for .get_info_crits() and get_LSmatrix()
  # res$w.resid is always a vector
  if (is.null(attr(res$w.resid,"unique"))) attr(res$w.resid,"unique") <- length(unique(res$w.resid))==1L # is.null() => for mv    
  #
  if (models[["eta"]]=="etaHGLM") {
    #
    sub_corr_info <- mget(c("corr_families","corr_types", "AMatrices"),  processed$corr_info) 
    kron_Y_LMatrices <- vector("list", nrand)
    for (it in seq_len(nrand)) kron_Y_LMatrices[it] <- list(attr(processed$corr_info$cov_info_mats[[it]],"blob")$Lunique)
    sub_corr_info$kron_Y_LMatrices <- kron_Y_LMatrices
    # 
    res$ranef_info <- list(sub_corr_info=sub_corr_info, hyper_info=processed$hyper_info,
                           vec_normIMRF=processed$AUGI0_ZX$vec_normIMRF)
    res$w.ranef <- wranefblob$w.ranef ## useful for .get_info_crits() and get_LSmatrix()
    #
    res$lambda.object <- .make_lambda_object(nrand, lambda_models=models[["lambda"]], cum_n_u_h, lambda_est, 
                                             process_resglm_blob, rand.families=processed$rand.families,
                                             ZAlist, LMatrices, lambdaType=attr(init.lambda,"type"))
    
    strucBlob <- .post_process_v_h_LMatrices(next_LMatrices=LMatrices, v_h=v_h, u_h=u_h,
                                             processed=processed, ranCoefs_blob=ranCoefs_blob) 
    res$strucList <- strucBlob$strucList
    res$v_h <- strucBlob$v_h
    # $lambda is a vector that NO LONGER contains lambdas for ranCoefs (too confusing, the more so as they are often 1)
    res[["lambda"]] <- structure(.get_lambdas_notrC_from_hlfit(res, type="adhoc"), cum_n_u_h=cum_n_u_h) ## redundant but very convenient except for programming
  } ## else all these res$ elements are NULL
  ###################
  ## ALL other PHI returns
  ###################
  if ( ! is.null(vec_nobs) ) {
    res$residModels <- vector("list", length(vec_nobs))
    cum_nobs <- attr(res$families,"cum_nobs")
    for (mv_it in seq_along(vec_nobs)) {
      res$residModels[[mv_it]] <- list(formula=processed$residModels[[mv_it]]$formula,
                                       family=attr(processed$residModels[[mv_it]]$family,"quoted"))
      if (models[["phi"]][mv_it]=="phiHGLM") {
        if (is.null(res[["resid_fits"]])) res[["resid_fits"]] <- vector("list", length(vec_nobs))
        res[["resid_fits"]][[mv_it]] <- PHIblob$multiPHI[[mv_it]]$phifit
      } else {
        if (is.null(res[["phi.object"]])) res[["phi.object"]] <- vector("list", length(vec_nobs))
        resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
        res[["phi.object"]][[mv_it]] <- .get_phi_object(phi.Fix[[mv_it]], PHIblob=PHIblob$multiPHI[[mv_it]], dev_res[resp_range], 
                                                        prior.weights=res$prior.weights[[mv_it]], 
                                                        phi.preFix=processed$phi.Fix[[mv_it]], nobs=vec_nobs[mv_it], 
                                                        control=processed[["control.glm"]])
      }
    }
    res[["phi"]] <- phi_est 
  } else {
    res$residModel <- list(formula=processed$residModel$formula,
                           family=attr(processed$residModel$family,"quoted"))
    if (models[["phi"]]=="phiHGLM") {
      res[["resid_fit"]] <- PHIblob$phifit
      res[["phi"]] <- phi_est 
    } else {
      res[["phi.object"]] <- .get_phi_object(phi.Fix, PHIblob, dev_res, prior.weights=res$prior.weights, 
                                             phi.preFix=processed$phi.Fix, nobs, control=processed[["control.glm"]])
      # phi_est comes from PHIblob$next_phi_est, not from final glm,  hence is in minimal form
      if (models[["phi"]]=="phiScal") {res[["phi"]] <- phi_est[1L]} else res[["phi"]] <- phi_est 
    }
  }
  ################### the magic environment
  if (models[[1L]]=="etaHGLM") {
    if (HL[1]=="SEM") {
      res$envir <- list2env(list(dvdloglamMat=dvdloglamMat, dvdlogphiMat=dvdlogphiMat), ## provided if available
                            parent=environment(HLfit_body))
      # bla
    } else {
      res$diagnostics$m_grad_obj <- auglinmodblob$m_grad_obj # typically NULL for LMM
      ## for beta_v_cov and cAIC's p_d
      # These elements are protected from deletion in stripHLfit() by explicit setdiff(stripnames,c("G_CHMfactor",...)):
      if ("AUGI0_ZX_sparsePrecision" %in% res$MME_method) {
        sXaug <- auglinmodblob$sXaug # 
        #sXaug$AUGI0_ZX$X.pv <- res[["X.pv"]] # NO that would affect the confint -> int_sXaug
        .init_promises_spprec(sXaug, non_X_ones=FALSE, nullify_X_ones =TRUE, intervalInfo=processed$intervalInfo) # initialize no promise, only removes X ones...
        #                                             (but confint might fully reinit them)
        # qrXa is generally an unevaluated promise => .calc_Md2hdvb2_info_spprec() tests its status. NULLified here
        # factor_inv_Md2hdv2 is used by spprec code for .calc_d2hdv2_info(). Kept as promise here
        # $envir aims to provide both the BLOB functionality and the $sXaug itself:
        envir <- sXaug$BLOB
        envir$sXaug <- sXaug # => BLOB within itself; memory-cheap 
        envir$dvdloglamMat <- dvdloglamMat
        envir$dvdlogphiMat <- dvdlogphiMat
        res$envir <- envir 
      } else {
        res$envir <- list2env(list(dvdloglamMat=dvdloglamMat, dvdlogphiMat=dvdlogphiMat,## provided if available
                                   sXaug=auglinmodblob$sXaug), ## F I X M E definitely useful, but could we remove some elements?
                              # big-ranefs.R is a good test
                              parent=environment(HLfit_body))
        attr(res$envir$sXaug,"scaled:scale") <- attr(processed$AUGI0_ZX$X.pv,"scaled:scale")
      }
    }
  } else res$envir <- list2env(list(dvdlogphiMat=dvdlogphiMat), ## provided if available
                               parent=environment(HLfit_body)) # scale not used post-fit for GLMs, otherwise it might be worth to provide it as above.
  ###################
  ## WARNINGS
  ###################
  .hack_options_error(message=NULL)
  ## translation of warnings in user-more friendly form 
  if ( ! is.null(warningList$resLam0) && warningList$resLam0) { 
    warningList$resLam0 <- "lambda residuals numerically 0 were replaced by 1e-6"
  }
  if ( ! is.null(warningList$resLamInf) && warningList$resLamInf) { 
    warningList$resLamInf <- "lambda residuals numerically >1e10 were replaced by 1e10"
  }
  if (! is.null(warningList$leveLam1) && warningList$leveLam1) {
    warningList$leveLam1 <- paste("lambda leverages numerically 1 were replaced by 1-",
                                  .spaMM.data$options$regul_lev_lambda,"(as controlled by option 'regul_lev_lambda')")
  }
  if ( ! is.null(warningList$resPhi0) && warningList$resPhi0) { 
    warningList$resPhi0 <- "phi residuals numerically 0 were replaced by 1e-6"
  }
  if ( ! is.null(warningList$resPhiInf) && warningList$resPhiInf) { 
    warningList$resPhiInf <- "phi residuals numerically >1e10 were replaced by 1e10"
  }
  if (! is.null(warningList$levePhi1) && warningList$levePhi1) {
    warningList$levePhi1 <- "phi leverages numerically 1 were replaced by 1 - 1e-8"
  }
  if (! is.null(warningList$negLevLam) && warningList$negLevLam) {
    warningList$negLevLam <- "Negative leverages for lambda were replaced by 1e-8"
  }
  if (! is.null(locw <- warningList$innerPhiGLM)) {
    warningList$innerPhiGLM <- paste0("'",locw,"' in some sub-final iteration(s) of phi estimation;")
  }
  if (! is.null(locw <- warningList$innerLamGLM)) {
    warningList$innerLamGLM <- paste0("'",locw,"' in some sub-final iteration(s) of lambda estimation;")
  }
  if (! is.null(locw <- processed$envir$inits_by_glm$conv_info)) {
    warningList$inits_by_glm <- locw ## added 2019/04/03
  }
  if ( HL[1]!="SEM" && maxit.mean>1 
       && ( models[[1]]=="etaHGLM" || pforpv>0L) ## cases where iterations are needed : including pforpv>0L for GLM
       && innerj==maxit.mean ) {
    warningList$innerNotConv <- paste0("linear predictor estimation did not converge;",
                                       if ( models[[1]]=="etaHGLM" && ! LMMbool && ! processed$LevenbergM["LM_start"] ) {
                                         " try control.HLfit=list(LevenbergM=TRUE), or"
                                       },
                                       " increase 'control.HLfit$max.iter.mean' above ",maxit.mean)
  }
  if ( (! is.na(conv_logL)) && iter==max.iter) {
    maxitmess <- paste0("Estimates did not converge;",
                        if ( ! LMMbool && ! processed$LevenbergM["LM_start"] ) {
                          " try control.HLfit=list(LevenbergM=TRUE), or"
                        },
                        " increase 'max.iter' above ",max.iter,
                        "\n (see help('HLfit') for details about 'max.iter')")
    if (models[["eta"]]=="etaHGLM") {
      if (conv_logL  && ! conv.lambda) {
        mainNotConv <- paste0("p_v apparently converged but lambda estimates apparently did not.",
                              "\n This may indicate that some lambda estimate(s) should be zero.",
                              "\n Otherwise try increasing 'max.iter' above ",max.iter,
                              "\n (see help(HLfit) for details about 'max.iter')")          
      } else mainNotConv <- maxitmess        
    } else mainNotConv <- maxitmess        
    attr(mainNotConv,"diagnostics") <- c( conv.phi=conv.phi , conv.lambda=conv.lambda, 
                                          conv.corr=conv.corr )
    warningList$mainNotConv <- mainNotConv
  }
  res$warnings <- warningList
  ###
  ### experimental cAIC minimization (completely experimental)
  if ( identical(processed$return_only,"cAICAPHLs")) {
    APHLs <- .get_info_crits(res)["cAIC"]
    if ( ! is.null(oldcAIC <- processed$port_env$objective)) {
      .update_port_fit_values(old_obj= - oldcAIC,new_obj= - APHLs[["cAIC"]], 
                              port_fit_values=list(fixef=beta_eta,v_h=v_h,phi_est=phi_est), 
                              models=models, processed=processed, control.HLfit=control.HLfit,
                              lambda_est=lambda_est, 
                              PHIblob=PHIblob)
    } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
      assign("objective",APHLs[["cAIC"]],envir=processed$port_env)
    }
    res <- list(APHLs=APHLs)
    return(res)    ########################   R E T U R N
  }
  ###################
  ## SUMMARY, RETURN
  ###################
  class(res) <- c("HLfit",class(res)) 
  if (verbose["all_objfn_calls"]) {
    seriousWarnings <- warningList[intersect(c("innerNotConv","mainNotConv"),names(warningList))]
    if (length(seriousWarnings) ) { 
      warningsss <- paste0("In HLfit :\n",unlist(seriousWarnings))
      abyss <- sapply(warningsss, warning, call.=FALSE) 
      warningList[setdiff(names(warningList),c("innerNotConv","mainNotCov"))] <- NULL
    }
  } 
  if (verbose["trace"]) {
    if (length(warningList) ) {
      warningsss <- paste0(unlist(warningList),"\n")
      abyss <- sapply(warningsss,cat) 
    }
  }
  # cleanup: for diagnostic, use
  # sort(sapply(ls(<object>$envir), function(x)
  # +             object.size(get(x, envir = <object>$envir))),decreasing=TRUE)
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"res")) ## empties the whole local envir except the return value
  return(res)
}
