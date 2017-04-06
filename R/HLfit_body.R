HLfit_body <- function(processed, #resid.model= ~ 1, 
                       verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),
                       control.HLfit=list(), ## used both by preprocess and HLfit_body
                       init.HLfit = list(), ## apparently not used by preprocess
                       ranFix=list(), ## phi, lambda, possibly nu, rho if not in init.HLfit
                       etaFix=list() ## beta, v_h (or even u_h)
) {
  data <- processed$data
  family <- processed$family
  ranFix <- post_process_family(family,ranFix) ## assign 'extra' pars and cleans ranFix
  if (spaMM.getOption("wDEVEL2")) {  ## must ensure match with future object
    ## rPparlist NULLin direct call of HLfit
    if ( is.null(rPparlist <- attr(ranFix,"parlist"))) {
      #attr(ranFix,"types") <- typelist ## used for:
      #parlist <- merge_parlist(,new=ranFix,types=NULL) ## will use attr(ranFix,"types") --> "fix"
      parlist <- merge_parlist(,new=ranFix,types="fix") ## will use attr(ranFix,"types") --> "fix"
      parlist <- merge_parlist(parlist,new=init.HLfit,types="var")
      attr(ranFix,"parlist") <- parlist ## pallist has fix+var contrary to ranFix[]
      #str(ranFix)
    } else {
      ## ranPars may have a preesisting parlist but eg with trLambda instead of lambda
      if ( ! identical(rpType <- attr(ranFix,"type"), rpTypes <- lapply(attr(rPparlist,"types"),tolower))) {
        #stop(" ! identical(attr(ranPars,\"type\"),attr(attr(ranPars,\"parlist\"),\"types\"))")
        # can be different bc rpType handles named vectors (rho,lambda) differently...
        # and also bc attr(ranFix, "type") may be empty list() contrary to attr(ranFix, "types")  
        #str(rpType)
        #str(rpTypes)
      }
    }
  }
  rpblob <- canonizeRanPars(ranPars=ranFix,corr.model="",checkComplete = FALSE) 
  ranFix <- rpblob$ranPars ## including full-size lambda
  
  formula <- processed$predictor
  prior.weights <- processed$prior.weights
  
  warningList <- list()
  ## when adding verbose elements, remind that these might be lost through corrHLfit -> HLCor cf dotlist$verbose <- verbose[intersect(...]
  #  verbose <- setControl(verbose=verbose)
  ##
  phi.Fix <- processed$phi.Fix
  if (is.null(phi.Fix)) phi.Fix <- getPar(ranFix,"phi") ## if set in final call of outer estimation 
  #
  nobs <- nrow(data) ## before prior.weights is evaluated
  predictor <- processed$predictor
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  HL <- processed$HL
  stop.on.error <- processed$stop.on.error ## to control issues with matrix computations; F by default
  conv.threshold <- processed$conv.threshold
  iter.mean.dispFix <- processed$iter.mean.dispFix
  iter.mean.dispVar <- processed$iter.mean.dispVar
  max.iter <- processed$max.iter
  resid.predictor <- processed$residModel$predictor 
  BinomialDen <- processed$BinomialDen
  y <- processed$y
  REMLformula <- processed$REMLformula ## should no longer be modified
  X.pv <- processed$`X.pv`
  X.Re <- processed$`X.Re` ## should be NULL _xor_ distinct from X.pv
  ### a bit of post processing
  nobs <- NROW(X.pv)
  pforpv <- ncol(X.pv)
  LMMbool <- processed$LMMbool
  models <- processed$models
  #### Note that HLCor modifies the L matrix (inprocessed$predictor if required) => ZAL cannot be preprocessed by corHLfit and must be recomputed each time 
  ZAlist <- processed$ZAlist ## : ZAlist is a list of design matrices 
  nrand <- length(ZAlist)
  cum_n_u_h <- processed$cum_n_u_h
  n_u_h <- cum_n_u_h[nrand+1L] 
  vec_n_u_h <- diff(cum_n_u_h)
  if (models[["eta"]]=="etaHGLM") {
    LMatrix <- processed$AUGI0_ZX$envir$LMatrix
    if ( ! is.null(LMatrix) ) { ## then LMatrix is presumably variable 
      ## [for CAR, HLCor_body must have updated LMatrix as fn of rho: seek adjd usage]
      ZALlist <- .compute_ZAXlist(XMatrix=LMatrix,ZAlist=processed$ZAlist)
      userLfixeds <- attr(ZALlist,"userLfixeds") # for MakeCovEst, set to (a vector of) FALSE
      ZAL <- .post_process_ZALlist(ZALlist,as_matrix=.eval_as_mat_arg(processed)) 
      ## ZAL may be modified by other call to .post_process_ZALlist   
    } else { ## ZAL was  accessible to preprocess; it must be in attr(predictor,"ZALMatrix")
      ZAL <- processed$AUGI0_ZX$ZAfix
      userLfixeds <- rep(FALSE,nrand) ## length matters if random-slope  as second effect 
    }
  } else if (models[[1]]=="etaGLM") {
    ZAL <- NULL ## 
    ZALlist <- NULL
    u_h <- v_h <- lev_lambda <- numeric(0)
  }
  lambda.Fix <- getPar(ranFix,"lambda") ## should already have length 'nrand' or else be NULL
  if (is.null(lambda.Fix)) lambda.Fix <- rep(NA,nrand) ## FIXME: put in canonizeRanpars ? depends onother uses of this fn.
  if (any(lambda.Fix[!is.na(lambda.Fix)]==0)) stop(pastefrom("lambda cannot be fixed to 0.",prefix="(!) From "))
  ###
  off <- attr(processed$predictor,"offsetObj")$total
  next_cov_est_vec <- NULL ## will be tested
  ##################
  unknowns <- setdiff(names(init.HLfit),c("fixef","phi","lambda","v_h","rho","nu","Nugget","ARphi"))   
    if (length(unknowns)>0) {
      mess <- pastefrom("unhandled elements in 'init.HLfit'.",prefix="(!) From ")
      message(mess)
      if ("beta" %in% unknowns) message("  Use 'fixef' rather than 'beta' in 'init.HLfit'.")
      stop()
    }
  ###################
  corr_est <- init.HLfit[intersect(c("nu","rho","Nugget","ARphi"),names(init.HLfit))]
  if (length(corr_est)==0L) {
    corr_est <- NULL
  } else {
    corrEstBlob <- eval.corrEst.args(family=family,rand.families=rand.families,predictor=predictor,data=data,X.Re=X.Re,
                                     REMLformula=REMLformula,ranFix=ranFix,
                                     Optimizer=control.HLfit$Optimizer)
    corrEst.args <- corrEstBlob$corrEst.args ## but corrEstBlob also has $corrEst.form which will stay there for later use
  }
  need_ranefPars_estim <-  (models[[1]]=="etaHGLM" && (anyNA(lambda.Fix) || ! is.null(corr_est)))
  #
  whichadj <- which(attr(attr(ZAlist,"ranefs"),"type")=="adjacency")
  test_inner_estim_rho <- (length(whichadj)>0L && ! is.null(adjd <- attr(LMatrix,"symsvd")$adjd))
  nothing_to_fit <-  ((! need_ranefPars_estim) && ncol(X.pv)==0L && (! is.null(phi.Fix)) 
                      && (models[[1]]=="etaHGLM" && (! is.null(etaFix$v_h))) )
  init_adjacency_info <- NULL
  if (nothing_to_fit && test_inner_estim_rho) {  
    u.range <- (cum_n_u_h[whichadj]+1L):cum_n_u_h[whichadj+1L]
    adj_rho <- corr_est$rho
    if (is.null(adj_rho)) adj_rho <- getPar(ranFix,"rho") ## could this occur with test_inner_estim_rho ?
    if (is.null(adj_rho)) adj_rho <- init.HLfit$rho
    if ( ! is.null(adj_rho)) init_adjacency_info <- list(whichadj=whichadj, u.range=u.range, coeffs=1/(1-adj_rho*adjd))
  } 
  
  ### case where nothing to fit #############################################
  if ( nothing_to_fit ) { ## nothing to fit. We just want a likelihood
    ### a bit the same as max.iter<1 ... ?
    phi_est <- phi.Fix
    eta <- off
    if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
      ## we need u_h in calc_APHLS...() and v_h here for eta...
      v_h <- etaFix$v_h
      u_h <- etaFix$u_h
      if (is.null(u_h)) {u_h <- .u_h_v_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,
                                                  lcrandfamfam=lcrandfamfam,lower.v_h=NULL,upper.v_h=NULL)}
      lambda_est <- resize_lambda(lambda.Fix,vec_n_u_h,n_u_h, adjacency_info=init_adjacency_info)
      eta <- eta + drop(ZAL  %id*id%  etaFix$v_h) ## updated at each iteration
    } ## FREQS
    ## conversion to mean of response variable (COUNTS for binomial)
    muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
    w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
      wranefblob <- .updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## no fit, likelihood computation
      if (models[["phi"]]=="phiScal") {
        H_global_scale <- phi_est[1L]
      } else H_global_scale <- exp(mean(log(phi_est)))
      ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
      Xscal <- make_Xscal(ZAL, ZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX)
      weight_X <- sqrt(H_global_scale*w.resid) ## sqrt(s^2 W.resid)
      if (inherits(Xscal,"Matrix")) {
        mMatrix_method <- .spaMM.data$options$Matrix_method
      } else mMatrix_method <- .spaMM.data$options$matrix_method
      sXaug <- do.call(mMatrix_method,
                       list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
    } else sXaug <- NULL 
    return(list(APHLs=calc_APHLs_from_ZX(processed=processed, which="p_v", sXaug=sXaug, phi_est=phi_est, 
                                         lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu)))
    ### RETURN !! ## FR->FR but p_bv is not returned.
  } 
  
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################

  ### Initial values  for lambda, phi and beta from lambda.Fix, phi.Fix, or init.HLfit ##### 
  ## Initial estimate for beta  (etaFix acts directly in auglinmodfit)
  beta_eta <- processed$port_env$port_fit_values$fixef
  if (is.null(beta_eta)) {
    beta_eta <- numeric(pforpv)
    if (pforpv>0) beta_eta <- init.HLfit$fixef
  }
  ## Initial estimate for phi 
  phi_est <- phi.Fix
  if (is.null(phi_est)) phi_est <- processed$port_env$port_fit_values$phi_est
  if (is.null(phi_est)) phi_est <- init.HLfit$phi ## must be 'response' value(s)
  if (models[[1]]=="etaHGLM") { 
    ## Initial estimate for u_h, v_h 
    gaussian_u_ranges <- processed$gaussian_u_ranges 
    psi_M <- rep(attr(rand.families,"unique.psi_M"),diff(cum_n_u_h))
    v_h <- initialize_v_h(psi_M=psi_M, etaFix=etaFix, init.HLfit=init.HLfit,
                          cum_n_u_h=cum_n_u_h, rand.families=rand.families, port_env=processed$port_env)
    u_h <- .u_h_v_h_from_v_h(v_h, rand.families=rand.families, cum_n_u_h=cum_n_u_h,
                            lcrandfamfam=lcrandfamfam, lower.v_h=NULL, upper.v_h=NULL)
    ## Initial estimate for lambda in 'compact" form
    init.lambda <- lambda.Fix ## already the right size 'nrand' with NA's or non-fixed ones
    lambdaType <- rep("",nrand) 
    whichInner <- which(is.na(lambda.Fix))
    lambdaType[whichInner] <- "inner" ## sure
    lambdaType[lambdaType=="" & is.na(processed$lambda.Fix)] <- "outer" ## sure
    lambdaType[lambdaType==""] <- "fix" ## sure
    if (! is.null(init.HLfit$lambda)) init.lambda[whichInner] <- init.HLfit$lambda[whichInner]
  } else init.lambda <- NULL
  ###
  ######### missing Initial estimates for mu, phi, lambda by GLM ####################
  if (  is.null(beta_eta) || is.null(phi_est) || anyNA(init.lambda) ) { 
    reset <- (family$family %in% c("COMPoisson","negbin"))
    inits_by_glm <- .get_inits_by_glm(processed,family=family,reset=reset) # using possibly postprocessed family
    ## : uses processed$y, $BinomialDen, attr($predictor,"offsetObj"), $control.glm.
  }
  ## Finalize initial values for beta_eta
  if (is.null(beta_eta) ) beta_eta <- inits_by_glm$beta_eta
  intervalInfo <- control.HLfit$intervalInfo
  if (!is.null(intervalInfo)) {
    parmcol <- attr(intervalInfo$parm,"col")
    beta_eta[parmcol] <- intervalInfo$init 
  }  
  ## Finalize initial values for phi 
  if (is.null(phi_est) ) {
    ## there are many cases that .get_init_phi does not handle
    if ( models[["phi"]] == "phiScal") {
      if (is.call(processed$prior.weights)) {
        phi_est <- .get_init_phi(processed,weights=NULL)
      } else phi_est <- .get_init_phi(processed,weights=processed$prior.weights)
    }
    if ( is.null(phi_est) ) phi_est <- inits_by_glm$phi_est
    if (models[["phi"]] != "phiScal") {
      phi_est <- rep(phi_est,nobs) ## moche ## FR->FR why is this necess ?
    }
  }
  ## Finalize initial values for lambda
  if (models[[1]]=="etaHGLM") { ## the basic case (LMM, GLMM...)
    ## we typically need  both init.lambda and lambda_est...
    if ( anyNA(init.lambda) ) { ## some inner may remain undetermined
      NA_in_compact <- is.na(init.lambda)
      prev_lambda_est <- processed$port_env$port_fit_values$lambda_est
      if (is.null(prev_lambda_est)) {
        stillNAs <- which(NA_in_compact)
        ## if reset is TRUE init.lambda is recomputed. Otherwise it is computed only once and then 'got'
        default.init.lambda <- .get_init_lambda(processed, reset=reset, stillNAs=stillNAs,
                                                         init_lambda_by_glm=inits_by_glm$lambda)
        init.lambda[stillNAs] <- default.init.lambda[stillNAs]
        lambda_est <- resize_lambda(init.lambda,vec_n_u_h,n_u_h, adjacency_info=init_adjacency_info)
      } else { ## do not replace current fixed values !
        lambda_est <- resize_lambda(init.lambda,vec_n_u_h,n_u_h, adjacency_info=init_adjacency_info)
        NA_in_expanded <- resize_lambda(NA_in_compact,vec_n_u_h,n_u_h, adjacency_info=init_adjacency_info)
        lambda_est[NA_in_expanded] <- prev_lambda_est[NA_in_expanded]
      } 
    } else lambda_est <- resize_lambda(init.lambda,vec_n_u_h,n_u_h, adjacency_info=init_adjacency_info)
    attr(init.lambda,"type") <- lambdaType
  }
  ###
  ## predictor from initial values
  if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
    eta <- off + drop(X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
  } else  eta <- off +drop(X.pv %*% beta_eta) ## no iteration hence no updating  ## FREQS
  ## conversion to mean of response variable (COUNTS for binomial)
  muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
  mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
  conv.phi <- conv.lambda <- conv.corr <- FALSE
  if (models[[1]]=="etaHGLM") {
    if (! is.null(intervalInfo)) {
      intervalInfo$parmcol_ZX <- n_u_h+parmcol 
      intervalInfo$parmcol_X <- parmcol 
      intervalInfo$ranFix <- ranFix
      if (processed$HL[1]==0L) { intervalInfo$likfn <- "hlik" } else {intervalInfo$likfn <- "p_v"} ##  use h in PQL/L (> v1.5.59) 
    }
    wranefblob <- .updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## initilization !
    if (ncol(X.pv)==0L && !is.null(etaFix$v_h)) {
      maxit.mean <- 0 ## used in test near the end...
    } else if ( LMMbool && is.null(intervalInfo) ) {
      maxit.mean <- 1 ## 1 is sufficient for LMM as Hessian does not vary with beta_eta  => quadratic function;
      # ./ E X C E P T for calculation of confidence intervals: at least two intervalSteps are required. Quite visible when 
      # dispersion params areouter-estimated (fitme) in which case there isno outer iteration to compensate a small maxit.mean  
    } else { ## even h maximization in *G*LMMs 
      if ( ! is.null(phi.Fix) && ! anyNA(lambda.Fix)) { ## allFix hence no true outer iteration 
        maxit.mean <- iter.mean.dispFix 
      } else maxit.mean <- iter.mean.dispVar # If phi.Fix and lambda.Fix, the only way to have 'outer' convergence is to have 'inner' convergence
    } 
    # next_LMatrices <- LMatrix ## originally <- NULL, modified <- LMatrix 2015/06/03, cf notes of that day
    next_LMatrices <- replicate(length(ZAlist), list(NULL))
    if (!is.null(LMatrix)) next_LMatrices[[which(attr(ZAlist,"ranefs") %in% attr(LMatrix,"ranefs"))]] <- LMatrix
  } else if (models[[1]]=="etaGLM") {
    ## FR->FR on doit pouvoir mettre maxit.mean = 0 dans pas mal de cas ?
    ## attention toutefois à COMPoisson negbin...
    if ( ! is.null(phi.Fix)) { ## 
      maxit.mean <- iter.mean.dispFix 
    } else maxit.mean <- iter.mean.dispVar # If phi.Fix and lambda.Fix, the only way to have 'outer' convergence is to have 'inner' convergence
  }
  prev_lik <- -Inf
  conv_logL <- NA
  dvdlogphiMat <- NULL
  dvdloglamMat <- NULL
  penul_lambda <- NULL
  residProcessed <- processed$residProcessed ## currently NULL except for phiHGLM
  ########################################
  ######### Main loop ####################
  ########################################
  iter <- 0L
  prevmsglength <- 0L
  if (HL[1]=="SEM") { ## specif probit
    #SEMblob <- probitgem::SEMwrap(processed, ZAL, beta_eta, off, corr_est, init.lambda, lambda.Fix, LMatrix, verbose) ## passes CHECK if CRAN knows it
    SEMblob <- eval(as.call(c(quote(SEMwrap),list(processed=processed, ZAL=ZAL, beta_eta=beta_eta, 
                                                             off=off, corr_est=corr_est, init.lambda=init.lambda, 
                                                             lambda.Fix=lambda.Fix, LMatrix=LMatrix, verbose=verbose))))
    beta_eta <- SEMblob$beta_eta
    corr_est["rho"] <- SEMblob$corr_est["rho"] ## may again be NULL
    lambda_est <- predict(SEMblob$glm_lambda,type="response")
    u_h <- v_h <- SEMblob$v_h
    logLapp <- SEMblob$logLapp
    attr(logLapp,"seInt") <- SEMblob$seInt ## may be NULL
    ## for summary.HLfit -> .calc_beta_cov:
    eta <- off + drop(X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
    muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
    w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    wranefblob <- .updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## no fit, likelihood computation
    sqrt.ww <- sqrt(c(w.resid,wranefblob$w.ranef))
    ##
  } else while ( TRUE ) { ##the main loop with steps: new linear predictor, new leverages, new phi, new w.resid, new lambda, new fn(lambda)
    if (models[[1]]=="etaHGLM") {
      if (! is.null(intervalInfo)) {
        intervalInfo$corr_est <- corr_est ## currently needs to be added ex tempo... apparently only for warn_intervalStep()
        intervalInfo$beta_eta <- beta_eta ## ex tempo: the current estimate of the CI bound
      } ## else intervalInfo remains NULL
      if (models[["phi"]]=="phiScal") {
        H_global_scale <- phi_est[1L]
      } else H_global_scale <- exp(mean(log(phi_est)))
      if (spaMM.getOption("wDEVEL")) {
        adj_rho <- corr_est$rho
        if (is.null(adj_rho)) adj_rho <- getPar(ranFix,"rho") ## could this occur with test_inner_estim_rho ?
        if (is.null(adj_rho)) adj_rho <- init.HLfit$rho
        auglinmodblob <- fit_as_sparsePrecision(X.pv, ZAL, y, n_u_h=cum_n_u_h[nrand+1L],
                               H_global_scale=1,
                               lambda_est=lambda_est, muetablob=muetablob,
                               off=off, maxit.mean=maxit.mean, etaFix=etaFix,
                               ## for ! LMM
                               eta=eta,
                               ## supplement for LevenbergM
                               beta_eta=beta_eta,
                               ## supplement for ! GLMM
                               wranefblob, u_h, v_h, w.resid, phi_est,
                               for_init_z_args=list(nrand=nrand, psi_M=psi_M, ranFix=ranFix), # corr_est=corr_est),
                               for_intervals=intervalInfo,
                               ##
                               processed=processed,rho=adj_rho)
      } else auglinmodblob <- fit_as_ZX(X.pv, ZAL, y, n_u_h=cum_n_u_h[nrand+1L], 
                                 H_global_scale=H_global_scale, 
                                 lambda_est=lambda_est, muetablob=muetablob,
                                 off=off, maxit.mean=maxit.mean, etaFix=etaFix,
                                 ## for ! LMM
                                 eta=eta, 
                                 ## supplement for LevenbergM
                                 beta_eta=beta_eta,
                                 ## supplement for ! GLMM
                                 wranefblob, u_h, v_h, w.resid, phi_est,
                                 for_init_z_args=list(nrand=nrand, psi_M=psi_M, ranFix=ranFix), # corr_est=corr_est),
                                 for_intervals=intervalInfo,
                                 ##
                                 processed=processed)
      ##############################
      beta_eta <- auglinmodblob$beta_eta
      wranefblob <- auglinmodblob$wranefblob # updated only if !LMM
      muetablob <- auglinmodblob$muetablob # updated only if !LMM
      mu <- muetablob$mu ## testé par accident, necess dans test COMPoisson HLfit...
      w.resid <- auglinmodblob$w.resid # updated only if !LMM
      d2hdv2 <- auglinmodblob$d2hdv2 # updated only if !LMM; and NULL when fit_as_ZX is called
      v_h <- auglinmodblob$v_h
      u_h <- auglinmodblob$u_h
      eta <- auglinmodblob$eta
      innerj <- auglinmodblob$innerj
    } else if (models[[1]]=="etaGLM") {
      if (pforpv>0  && maxit.mean>0L) {
        ## nomme auglinmodblob pour avancer vers simplif du code, mais je ne peux suppser qu'auglinmodblob est toujours calculé
        auglinmodblob <- calc_etaGLMblob(processed=processed,  mu=mu, eta=eta, muetablob=muetablob, beta_eta=beta_eta, 
                                      w.resid=w.resid, phi_est=phi_est, off=off, maxit.mean=maxit.mean, 
                                      verbose=verbose, conv.threshold=conv.threshold)
        eta <- auglinmodblob$eta 
        muetablob <- auglinmodblob$muetablob 
        mu <- muetablob$mu
        beta_eta <- auglinmodblob$beta_eta 
        w.resid <- auglinmodblob$w.resid
        innerj <- auglinmodblob$innerj
      } else innerj <- 0L
    } # end etaGLM...
    if(inherits(mu,"Matrix")) mu <- drop(mu) ## pb calcul deviance_residual 
    # if (verbose["trace"]) {print(paste("beta=",paste(signif(beta_eta,4),collapse=", ")),quote=F)}
    ########## LEVERAGES
    #### base from hat matrix
    if (models[[1]]=="etaHGLM") {
      if (need_ranefPars_estim || is.null(phi.Fix)) {
        if (maxit.mean==0L) {
          stop("(!) Computation of leverages with maxit.mean=0: check that this is meaningful.")
        } # ELSE rWW was updated in the inner loop for betaV
        hatval <- .get_hatvalues(auglinmodblob$sXaug,X.Re=X.Re, auglinmodblob$weight_X)
        lev_lambda <- hatval[seq(n_u_h)]  ## for the ranef residuals (lambda)
        lev_phi <- hatval[(n_u_h+1L):(n_u_h+nobs)] ## for the error residuals (phi)
        hatval <- c(lev_phi,lev_lambda) ## old order; patch for further hatval corr 
      }
    } else if (is.null(phi.Fix)) { ## phi estim in GLM fitted by ML. Here I had more ad hoc code up to v1.7.42
      if ( ! is.null(X.Re) ) { 
        wAugXleve <- calc_wAugX(augX=X.Re,sqrt.ww=sqrt(w.resid)) # rWW%*%X.Re 
        lev_phi <- leverages(wAugXleve)
      } else { ## basic REML, leverages from the same matrix used for estimation of beta
        wAugX <- calc_wAugX(augX=X.pv,sqrt.ww=sqrt(w.resid)) # rWW %*% X.pv 
        lev_phi <- leverages(wAugX)
      }
    }
    ## (HL[2]=0, HL[3]=0): previous hat matrix -> p 
    ## (HL[2]=0, HL[3]=1): notEQL -> tilde(p), (HL[2]=1 && ): full correction -> q 
    ## (HL[2]=1, HL[3]=1): full correction -> q 
    #### contribution from GLM weights
    if (HL[2L]>0 && models[[1L]]=="etaHGLM" 
        && (anyNA(lambda.Fix) || is.null(phi.Fix)) ) { ## LeeN01 HL(.,1) ie the + in 'EQL+'
      ## first the d log hessian / d log lambda or phi corrections
      ### For the d log hessian first the derivatives of GLM weights wrt eta 
      ##################### noter que c'est le coef2 de HL(1,.), but mu,eta may have been updated since coef2 was computed
      dlW_deta <- .calc_dlW_deta(dmudeta=drop(muetablob$dmudeta),family=family,mu=drop(mu),eta=drop(eta),
                                BinomialDen=BinomialDen,canonicalLink=processed$canonicalLink)$dlW_deta
      ### we join this with the deriv of log w.ranef wrt v_h
      dlW_deta_or_v <- c(dlW_deta, wranefblob$dlogWran_dv_h )  ## vector with n+'r' elements
      # dlogWran_dv_h is 0 gaussian ranef; d2mudeta2 is 0 for identity link => vector is 0 for LMM
      ## else we continue the computation of the d log hessian term d2 log dens u/ dv dloglambda
      ## where we ignore the contribution of the log Jacobian, log(dth/du), to log dens u since it is not fn of lambda
      ## hence this is d2 log dens th(u)/ dv dloglambda
      if (any(dlW_deta_or_v!=0L)) {
        lev_phi_range <- 1L:nobs
        leve__dlW_deta_or_v <- hatval * dlW_deta_or_v
        leve__dlW_deta_or_v__ZALI <-  leve__dlW_deta_or_v[lev_phi_range] %*% ZAL +  leve__dlW_deta_or_v[-(lev_phi_range)]
        
        if (anyNA(lambda.Fix)) {
          sXaug <- auglinmodblob$sXaug
          dvdloglamMat <-  .calc_dvdloglamMat_new(dlogfthdth=(psi_M - u_h)/lambda_est, ## the d log density of th(u)
                                                 cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,rand.families=rand.families,
                                                 u_h=u_h,
                                                 sXaug=sXaug,
                                                 stop.on.error=stop.on.error)
          dleve <- as.vector(leve__dlW_deta_or_v__ZALI %*% dvdloglamMat) # (r+n).(r+n)Xr.rXr = r (each element is a sum over r+n terms= a trace)
          lev_lambda <- lev_lambda - dleve  
        } 
        ## 
        if (is.null(phi.Fix)) {
          dh0deta <- ( w.resid *(y-mu)/muetablob$dmudeta ) ## 12/2013 supp BinomialDen (soit Bin -> phi fixe=1, soit BinomialDen=1)
          sXaug <- auglinmodblob$sXaug
          dvdlogphiMat <- .calc_dvdlogphiMat_new(dh0deta=dh0deta,ZAL=ZAL, sXaug=sXaug, stop.on.error=stop.on.error)
          dleve <- as.vector(leve__dlW_deta_or_v__ZALI %*% dvdlogphiMat) # (r+n) . (r+n)Xr . rXn = n (each element is a sum over r+n terms= a trace)
          lev_phi <- lev_phi - dleve  
        } 
      }
    }
    if (HL[2L]>1) {stop("Need a_i correction in Table 7 of NohL07 ie derivatives of second order correction wrt dips param.")}
    #### contribution from exact likelihood function instead of EQL
    if (HL[3L]!=0 ) {## HL(.,.,1) ie , p_bv(h), not EQL p_bv(q+), LeeNP p89; distinction does not arise for PQL <=> Gaussian ranefs...  
      # lambda
      if (models[[1L]]=="etaHGLM" && anyNA(lambda.Fix)) ## d h/ d !log! lambda correction     
        lev_lambda <- lev_lambda + corr.notEQL.lambda(nrand,cum_n_u_h,lambda_est,lcrandfamfam) 
      # phi hence not poiss,binom:
      if (family$family=="Gamma" && is.null(phi.Fix) ) { ## d h/ d !log! phi correction (0 for gauss. resid. error). Not tied to REML
        phiscaled <- phi_est/eval(prior.weights) ## 08/2014 ## bug "*" corrected -> "/" 2015/03/05
        lev_phi <- lev_phi +  1+2*(log(phiscaled)+digamma(1/phiscaled))/phiscaled ## LNP p. 89 and as in HGLMMM IWLS_Gamma
      }    
    }
    ######### Dispersion Estimates for phi #####################
    if (is.null(phi.Fix)) { ## if phi is estimated (phi.Fix set to 1 for Poisson, Binomial)
      ## leverages have been computed before the  inner loop, which did not change the design matrices 
      lev_phi[lev_phi>1-1e-8] <- 1-1e-8
      ## updated residuals from updated mu must be used (LeeNP p.161) [not so in dhglmfit !!]
      ## Once in Gamma GLM, y=25.75563, y-mu=-5.996906e-08, yielded negative dev.resid
      if (models[["phi"]]=="phiHGLM") {
        residProcessed$prior.weights <- structure((1-lev_phi)/2,unique=FALSE) # expected structure in processed.
        # uses of prior weights matches thatin input of calcPHI -> dispGammaGLM 
        residProcessed$data$.phi <- family$dev.resids(y,mu,wt= eval(prior.weights))/(1-lev_phi) 
        residProcessed$y <- residProcessed$data$.phi
        #cat("\nbefore")
        phifitarglist <- processed$residModel
        ## which may include processed$residModel$fixed, which is a *list* created by preprocess
        phifitarglist$processed <- residProcessed 
        ## A resid.model's fixed lambda, phi have been preprocessed 
        ##   and info in processed$residModel$fixed must match that in residProcessed
        if (iter==0) {
          ## may be NULL, but allows a previous 'close' [si condition on port_fit_values$phifit_init] 
          ##   outer Hfit to be used to initialize the 1st inner fitme of the new HLfit 
          phifitarglist$processed$port_env$port_fit_values <- processed$port_env$port_fit_values$phifit_init
          phifit_init <- list()
          if (is.null(processed$residModel$fixed$phi)) phifit_init$phi <- 1 
        } else {
          phifit_init <- phifit$corrPars[attr(phifit$corrPars,"type")=="outer"] ## may be an empty list
          if (length(phifit_init)>0L) {
            ## FR->FR init must be in user scale, optr output is 'unreliable' => complex code
            LUarglist <- attr(phifit,"optimInfo")$LUarglist ## might be NULL 
            if (! is.null(LUarglist) && LUarglist$optim.scale=="transformed") {
              LowUp <- do.call("makeLowerUpper",LUarglist)
              if ( !is.null(phifit_init$nu) ) {
                NUMAX <- LUarglist$NUMAX
                initnu <- nuFn(phifit_init$nu,NUMAX=NUMAX) 
                initnu <- initnu +1e-4 *( (LowUp$lower$trNu+LowUp$upper$trNu)/2 -initnu)
                phifit_init$nu <- nuInv(initnu,NUMAX=NUMAX)
              }
              if ( !is.null(phifit_init$rho) ) {
                RHOMAX <- LUarglist$RHOMAX
                initrho <- rhoFn(phifit_init$rho,RHOMAX=RHOMAX) 
                initrho <- initrho +1e-4 *( (LowUp$lower$trRho+LowUp$upper$trRho)/2 -initrho)
                phifit_init$rho <- rhoInv(initrho,RHOMAX=RHOMAX)
              }
            }                      
          }
          if (all(is.na(residProcessed$lambda.Fix))) phifit_init$lambda <- phifit$lambda
          if (is.null(processed$residModel$fixed$phi)) phifit_init$phi <- 1 
          if (is.null(processed$residModel$fixed$fixef)) phifitarglist$init.HLfit$fixef <- fixef(phifit)
          phifitarglist$init.HLfit$v_h <- phifit$v_h
        }
        phifitarglist$init <- phifit_init
        phifitarglist$verbose <- .reformat_verbose(NULL,For="corrHLfit") ## constrained
        # if ( length(phifit_init)>0L) cat(paste("\nINIT phi fit in iter=",iter+1L,"     ",
        #           paste(names(phifit_init),"=", signif(unlist(phifit_init),6),
        #                 collapse=", ",sep=""),"\n",sep=""))
        phifit <- do.call("fitme_body",phifitarglist)
        #summary(phifit)
        prevmsglength <- overcat(paste("phi fit in HLfit's iter=",iter+1L,
                                       ", .phi[1]=",signif(phifit$y[1],5),", ",
                                       paste(c(names(phifit$corrPars),"lambda"),"=",
                                               signif(c(unlist(phifit["corrPars"]),phifit$lambda),6),
                                               collapse=", ",sep=""),
                                       ";           ",sep=""),prevmsglength
            )
        ## the final printing after iterations comes from HLfit.obj/HLCor.obj -> if (mc$verbose["objective"]) ...
        next_phi_est <- phifit$fv
      } else {
        # to obtain the phi estimate given by summary.glm(), one must use
        #  dev.res <- wt * ((y-mu)/family$mu.eta(eta))^2 ## wt * EQL residuals
        # but the logLik differs from that given by logLik.glm(). See Details in ?HLfit
        calcPHIblob <- calcPHI(oriFormula=attr(resid.predictor,"oriFormula"),
                               dev.res= family$dev.resids(y,mu,wt= eval(prior.weights)),
                               #: times pw to be an estimate of same phi accross level of response
                               # but not same phi as when there is no pw !
                               # double pw => double phi_est so that phi_est_i :=phi_est/pw_i is unchanged
                               data=data,
                               family=processed$residModel$family,
                               lev_phi=lev_phi,
                               control=processed$control.glm,
                               phimodel=models[["phi"]],
                               verbose=verbose,
                               control.phi=control.HLfit$`control.phi`)
        if (! is.null(locw <- calcPHIblob$glm_phi$warnmess)) warningList$innerPhiGLM <- locw
        next_phi_est <- calcPHIblob$next_phi_est # value of *phi* (not phi_i:= phi/prior.weights as pw are usd inGLMweights, not here)  
      }
      if (all(abs(next_phi_est-phi_est) < conv.threshold* (phi_est+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.phi <- TRUE ## 'weak convergence'... 
      } else conv.phi <- FALSE
    } else {conv.phi <- TRUE} ## there is a phi.Fix, already -> phi_est
    ###############################################################
    ######### Dispersion Estimates for lambda #####################
    ###############################################################
    
    if (need_ranefPars_estim) { ## lambda must be estimated 
      if (any(abs(lev_lambda) > 1 - 1e-8)) { ## abs... not commented when written...
        lev_lambda[lev_lambda>1-1e-8] <- 1-1e-8
        warningList$leveLam1 <- TRUE
      }      
      cov_est_vec <- next_cov_est_vec
      ################## ranefEstargs mustcontain arguments for makeCoveEst1 => its names are constrained
      ranefEstargs <- list(u_h=u_h,ZAlist=processed$ZAlist,cum_n_u_h=cum_n_u_h,
                           prev_LMatrices=next_LMatrices,
                           userLfixeds=userLfixeds,
                           processed=processed,
                           w.resid=w.resid)
      if (attr(processed$ZAlist,"anyRandomSlope")) { ## if random-slope model
        covEstmethod <- .spaMM.data$options$covEstmethod ## note same call in calcRanefPars
        if (is.null(covEstmethod)) stop("spaMM.getOption('covEstmethod') should not be NULL")
        if (covEstmethod == "makeCovEst1") {## le plus bourrin
          ranefEstargs <- c(ranefEstargs,list(phi_est=phi_est,
                                              as_matrix=( ! inherits(ZAL,"Matrix")),v_h=v_h))
          if (models[["phi"]]=="phiScal") {
            H_global_scale <- phi_est[1L]
          } else H_global_scale <- exp(mean(log(phi_est)))
          ranefEstargs$auglinfixedpars <- list(X.pv=X.pv, y=y, n_u_h=cum_n_u_h[nrand+1L], 
                                               H_global_scale=H_global_scale, 
                                               muetablob=muetablob,
                                               off=off, maxit.mean=maxit.mean, etaFix=etaFix,
                                               ## for ! LMM
                                               eta=eta, 
                                               ## supplement for LevenbergM
                                               beta_eta=beta_eta,
                                               ## supplement for ! GLMM
                                               u_h=u_h, v_h=v_h, w.resid=w.resid, phi_est=phi_est,
                                               for_init_z_args=list(nrand=nrand, psi_M=psi_M, ranFix=ranFix), # corr_est=corr_est),
                                               for_intervals=intervalInfo,
                                               ##
                                               processed=processed)
        } else { ## pour makeCovEst2
          covEstarglist$clik <- sum(processed$loglfn.fix(mu,y,eval(prior.weights)/phi_est)) ## constant over optim in cov estim
          covEstarglist$prevZAL <- ZAL
        }
      }
      calcRanefPars_blob <- calcRanefPars(corrEstList=list(corr_est=corr_est), 
                                          lev_lambda=lev_lambda,
                                          ranefEstargs=ranefEstargs,
                                          lambda.Fix=lambda.Fix,
                                          rand.families=rand.families,
                                          lcrandfamfam=lcrandfamfam,
                                          psi_M=psi_M,
                                          verbose=verbose,
                                          iter=iter,
                                          control=processed$control.glm
      )
      next_corr_est <- calcRanefPars_blob$next_corrEstList$corr_est
      next_cov_est_vec <- calcRanefPars_blob$next_corrEstList$cov_est_vec
      next_LMatrices <- calcRanefPars_blob$next_LMatrices ## random slope: must converge to the L factor of corr mat
      next_lambda_est <- calcRanefPars_blob$next_lambda_est  ## T H I S is what affects next beta,v estimates when rho is inner estimated in CAR.
      ########################################
      # => variation of log(u^2/lamdba) = simplified likRanU convergence  (from 1.9.31)
      next_u_vs_lambda <- 2*log(abs(u_h[gaussian_u_ranges])+1e-12)-log(next_lambda_est[gaussian_u_ranges])    
      conv_lambda_vs_u <- (iter> 1L && 
                             all( abs(next_u_vs_lambda-prev_u_vs_lambda) < 500*conv.threshold) )        
      prev_u_vs_lambda <- next_u_vs_lambda
      ## Absolu ou relatif selon valeur de lambda:
      conv_rel_lambda <- all( abs(next_lambda_est-lambda_est)/(lambda_est+0.1) < conv.threshold ) 
      # with 0.1 in denom, lambda converges at 1e-6 if conv.threshold is 1e-5; ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
      conv.lambda <- ( conv_lambda_vs_u && conv_rel_lambda )  
    } else { conv.lambda <- TRUE } ## end if anyNA lambda.Fix else ...
    #
    if (!is.null(next_cov_est_vec)) {
      if (iter>1 && abs(cov_est_vec-next_cov_est_vec) < conv.threshold ) { 
        conv.corr <- TRUE 
      } else conv.corr - FALSE       
    } else conv.corr <- TRUE
    #
    iter <- iter+1L ## here first from 0 to 1
    ###### convergence: 
    if ( conv.phi && conv.lambda && conv.corr) {
      break 
    } else if (iter>=max.iter) { ## did not converge...
      break 
    } else { ## update and continue
      if ( is.null(phi.Fix)) {
        phi_est <- next_phi_est
        w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## bc phi was updated. 'weinu', must be O(n) in all cases 
      }
      if (attr(processed$ZAlist,"anyRandomSlope") || ! is.null(corr_est) ) {
        ## FR->FR incompat entre randomslope code using XMatrix=next_LMatrices and code using XMatrix=LMatrix: 
        if (attr(processed$ZAlist,"anyRandomSlope")) {
          ZALlist <- .compute_ZAXlist(XMatrix=next_LMatrices, ZAlist=processed$ZAlist)
        } else if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="adjacency") {
          corr_est <- next_corr_est 
          #LMatrix est constant!= decomp$u
          #ZALlist <- .compute_ZAXlist(XMatrix=LMatrix,ZAlist=processed$ZAlist)
        } else if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="SAR_WWt") {
          corr_est <- next_corr_est ## 
          ## 
        }
        ZAL <- .post_process_ZALlist(ZALlist,as_matrix=( ! inherits(ZAL,"Matrix"))) 
      }
      if (models[[1]]=="etaHGLM" && (anyNA(lambda.Fix) || test_inner_estim_rho)) { ## next_lambda_est available:
        ##      in particular, also in the adjacency case if corr_est was updated but not the lambda param.
        ###################################################################################################
        if ( FALSE && NCOL(processed$X_lambda)==1L) {## FR->FR this almost works ./.
          ## ./. mais voir premier exemple dans ?spaMMqui boucle...comprends pas  
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
        lambda_est <- next_lambda_est
      }
      if (models[[1]]=="etaHGLM" && (anyNA(lambda.Fix) || ! is.null(corr_est))) { ## lambda or u_h were modified
        wranefblob <- .updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h)
      } 
      ## conv_logL either used to break the loop, Xor required only in last two iters for diagnostics 
      if (processed$break_conv_logL 
          || ( verbose["trace"] && iter>=max.iter-1L )) {
        next_lik <- calc_APHLs_from_ZX(auglinmodblob=auglinmodblob, which="p_v", processed=processed)
        # this does not depend on the latest ranPars estimates, sicne sXaug was not updated afterwards... 
        conv_logL <- abs(next_lik - prev_lik)/(0.1 + abs(next_lik)) < 1e-8 # ~ glm.fit convergence
        if (processed$break_conv_logL && conv_logL) break 
        prev_lik <- next_lik
      } else conv_logL <- NA
      ##
      if (verbose["trace"]) {
        print(paste("iteration ",iter,"; convergence criteria (phi, lambda, corr [, conv_lambda_vs_u, conv_rel_lambda]): ",
                    paste(c( conv.phi , conv.lambda, conv.corr, 
                             conv_lambda_vs_u, conv_rel_lambda),collapse = " "),sep=""))
        if (models[[1]]=="etaHGLM" && anyNA(lambda.Fix)) { 
          #print(range(logrel_crit))
          #print(range(reldlam_crit))
        }
        print("================================================")
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
      mess <- pastefrom(mess,prefix="(!) From ")
      message(mess)
    } else {
      message(paste("(beta,v)/lambda/phi iterations in HLfit() converged in",iter,"iterations"))
    }
  }
  #
  if (HL[1]=="SEM") {
    APHLs <- list(logLapp=logLapp) ## keeps attributes
  } else {
    if (models[[1]]=="etaHGLM") {
      # if (processed$HL[1]==0L) {p_v_obj <-"hlik"} else p_v_obj <-"p_v" # meme pour PQL/ l'objective n'est pas hlik
      if (identical(processed$return_only,"p_vAPHLs")) {
        whichAPHLs <- "p_v"
      } else if (identical(processed$return_only,"p_bvAPHLs")) {
        whichAPHLs <- "p_bv" ## retrun value may still include p_v if it is used to compute p_bv
      } else whichAPHLs <- c("p_v","p_bv")
      APHLs <- calc_APHLs_from_ZX(auglinmodblob,which=whichAPHLs,processed)
    } else { 
      APHLs <- calc_APHLs_from_ZX(auglinmodblob=NULL,processed, which="p_v",
                                          sXaug=NULL, phi_est, lambda_est=NULL, dvdu=NULL, u_h=NULL, mu ) 
    }
  }
  ######################### potential R E T U R N here: cases without p_bv
  if ( identical(processed$return_only,"p_vAPHLs")) {
    # a bit of ugly coding, but optimthroughsmooth calls HLCor, not HLCor.obj, thus it cannot directly control return_only. So either leave as is, or move the test to HLCor, or modify optimthroughsmooth to call HLCor.obj  
    if (HL[1]=="SEM") { # lambda used for smoothing.
      res <- list(APHLs=APHLs,lambda=SEMblob$lambda) 
    } else {
      res <- list(APHLs=APHLs) 
      ### voir verbose["objective"] pour des affichages dans les fonctions .obj
      if ( ! is.null(oldp_v <- processed$port_env$objective)) {
        if ( APHLs$p_v> oldp_v ) {
          port_fit_values <- list(fixef=beta_eta,v_h=v_h,phi_est=phi_est)
          if (models[[1]]=="etaHGLM") port_fit_values$lambda_est <- lambda_est
          ## and use residProcessed$port_env$port_fit_values (if any) to initiate (iter=0) 
          ## the next fitme call for phi_fit, in the next parent call to HLfit
          port_fit_values$phifit_init <- processed$residProcessed$port_env$port_fit_values
          ## Use FALSE to inhibit all port_env usage:
          assign("port_fit_values",port_fit_values,envir=processed$port_env)
          if ( ! identical(control.HLfit$write_port_env,FALSE)) assign("objective",APHLs$p_v,envir=processed$port_env) 
        } else { ## keep objective value; two cases for init values:
          if (abs(oldp_v-APHLs$p_v)<1) {
            ## small decrease: do nothing, keep old port_env_values
          } else assign("port_fit_values",NULL,envir=processed$port_env) ## remove starting value, no useful for large variations 
        }
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
      if (is.null(X.Re)) {X.REML <- X.pv} ## REML standard
      if ( ncol(X.REML)>0 ) { ## REML standard || REML non standard
        Md2hdb2 <- .ZtWZwrapper(X.REML,w.resid)
        ladb <- LogAbsDetWrap(Md2hdb2,logfac=-log(2*pi))
      } else ladb <- 0 ## fit ML standard : ncol(X.Re)=0, p_bv=p_v hence d2hdpbv reduces to d2hdv2
      APHLs <- list(p_v=APHLs$clik, p_bv=APHLs$clik - ladb/2) ## 07/2016 for inference about phi in GLMs 
    }
  }
  ######################### potential R E T U R N here: with p_bv
  if ( identical(processed$return_only,"p_bvAPHLs")) {
    ### voir verbose["objective"] pour affichages (.obj)
    if (FALSE) {
      if ( ! is.null(oldp_bv <- processed$port_env$objective)) {
        if (abs(oldp_bv-APHLs$p_bv)<1) {
          #cat(APHLs$p_bv,"")
          port_fit_values <- list(fixef=beta_eta,v_h=v_h,phi_est=phi_est)
          if (models[[1]]=="etaHGLM") port_fit_values$lambda_est <- lambda_est
          ## and use residProcessed$port_env$port_fit_values (if any) to initiate (iter=0)
          ## the next fitme call for phi_fit, in the next parent call to HLfit
          port_fit_values$phifit_init <- processed$residProcessed$port_env$port_fit_values
        } else {
          #cat("n ")
          port_fit_values <- NULL
        }
        #cat(beta_eta[1],phi_est[1],lambda_est[1],"\n")
        assign("port_fit_values",port_fit_values,envir=processed$port_env)
      }
      if ( ! identical(control.HLfit$write_port_env,FALSE)) assign("objective",APHLs$p_bv,envir=processed$port_env) ## use FALSE to inhibit all port_env usage
    } else {
      if ( ! is.null(oldp_bv <- processed$port_env$objective)) {
        if ( APHLs$p_bv> oldp_bv ) {
          port_fit_values <- list(fixef=beta_eta,v_h=v_h,phi_est=phi_est)
          if (models[[1]]=="etaHGLM") port_fit_values$lambda_est <- lambda_est
          ## and use residProcessed$port_env$port_fit_values (if any) to initiate (iter=0) 
          ## the next fitme call for phi_fit, in the next parent call to HLfit
          port_fit_values$phifit_init <- processed$residProcessed$port_env$port_fit_values
          ## Use FALSE to inhibit all port_env usage:
          assign("port_fit_values",port_fit_values,envir=processed$port_env)
          if ( ! identical(control.HLfit$write_port_env,FALSE)) assign("objective",APHLs$p_bv,envir=processed$port_env) 
        } else { ## keep objective value; two cases for init values:
          if (abs(oldp_bv-APHLs$p_bv)<1) {
            ## small decrease: do nothing, keep old port_env_values
          } else assign("port_fit_values",NULL,envir=processed$port_env) ## remove starting value, no useful for large variations 
        }
      } else if ( ! identical(control.HLfit$write_port_env,FALSE)) {
        assign("objective",APHLs$p_bv,envir=processed$port_env)
      }
    }
    return(list(APHLs=APHLs))    ########################   R E T U R N
  }
  # beta_cov code removed from here in v1.9.24
  ######################
  ######################
  ######################
  ##### LAMBDA and other RANEF PARS
  # there is one lambda 
  ## FR->FR the code defines lambda.Fix as a vector with one element per ranef
  ## hence its not clear how lambda.Fix  is handled in this case
  ## but coefficients lambda (list) below may have a two-element vector in this case ?
  ## (1) we count inner-estimated ranef pars
  if (need_ranefPars_estim) {
    bloc_lambda_args <- list(models=models, #init.lambda=init.lambda, 
                             processed=processed, lambda.Fix=lambda.Fix, 
                             cum_n_u_h=cum_n_u_h, next_LMatrices=next_LMatrices)
    if (HL[1]=="SEM") {
      bloc_lambda_args$SEMblob <- SEMblob
    } else {
      bloc_lambda_args$calcRanefPars_blob <- calcRanefPars_blob
      bloc_lambda_args$lev_lambda <- lev_lambda
    }
    #
    process_resglm_blob <- do.call(".bloc_lambda",bloc_lambda_args)
    #
    coefficients_lambdaS <- process_resglm_blob$coefficients_lambdaS # list
    p_lambda <- length( unlist(process_resglm_blob$coefficients_lambdaS)) - length( which(attr(init.lambda,"type")=="fix")) 
    p_adjd <- sum(unlist(lapply(process_resglm_blob$coefficients_lambdaS,function(v) which(names(v)=="adjd"))))
    p_lambda <- p_lambda + p_adjd
    if (attr(processed$ZAlist,"anyRandomSlope")) {
      ## lignes suiv supposent que L_matrix decrit random slope model
      p_corr <- sum(unlist(lapply(next_LMatrices,function(mat) {
        dimL <- nrow(attr(mat,"Lcompact"))
        (dimL-1)*dimL/2
      })))
      p_lambda <- p_lambda+p_corr
    }
  } else {
    p_lambda <- 0
    process_resglm_blob <- list(print_lambdas=as.list(rep(NA,nrand)))
  }
  ## (2) we count outer estimated ones
  if (! is.null(preproFix <- processed$lambda.Fix))
  p_lambda <- p_lambda + length(which(is.na(preproFix[ ! is.na(lambda.Fix)])))
  ##### PHI: 
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
  res$dfs <- c(pforpv=pforpv, p_lambda=p_lambda, p_phi=processed$p_phi)
  ###################
  ## DATA
  ###################
  res$data <- data ## very useful for simulate...
  if (family$family=="binomial") {
    res$weights <- BinomialDen
  }
  res$y <- y ## counts for Pois/bin
  #res$prior.weights <- structure(eval(prior.weights), unique=identical(attr(prior.weights,"unique"),TRUE)) ## see Gamma()$simulate
  res$prior.weights <- structure(eval(prior.weights), unique=attr(prior.weights,"unique")) ## see Gamma()$simulate
  ###################
  ## MODEL info
  ###################
  res$family <- family
  res$X.pv <- X.pv
  res$ranFix <- ranFix ## currently as a uniform template consistent with projected changes ; excpt that lamFix, phiFix info is now in lambda.object, etc
  #if (spaMM.getOption("wDEVEL2")) str(ranFix)
  correst_and_ranfix <- .get_CorrEst_and_RanFix(ranFix, corr_est)
  res$corrPars <- correst_and_ranfix$corrPars ## subset of the following:
  res$CorrEst_and_RanFix <- correst_and_ranfix$CorrEst_and_RanFix 
  res$models <- models
  res$fixef_terms <- processed$fixef_terms ## added 2015/12/09 for predict
  res$fixef_levels <- processed$fixef_levels ## added 2015/12/09 for predict
  res$predictor <- predictor ##  all post fitting functions expect PROCESSED predictor
  ## noter que le 'call'$formula peut encore avoir attr LMatrix 
  attr(res$predictor,"LMatrix") <- NULL ## processed$AUGI0_ZX$envir$LMatrix ## ?FIXME? at least temporary fix
  #
  if (models[[1]] == "etaHGLM") res$ZAlist <- processed$ZAlist ## needed for prediction variance
  res$REMLformula <- REMLformula ## copy without modif of processed$REMformula, given that 'processed' is not returned
  ###################
  ## ALGORITHM
  ###################
  res$HL <- HL ## info on fitting method
  ###################
  ## FITTED VALUES
  ###################
  if (family$family=="binomial") {
    res$fv <- mu/BinomialDen ## cf glm(binomial): fitted values are frequencies 
  } else {res$fv <- mu} ## fitted values may be counts (cf poisson), or reals
  ###################
  ## FIXEF, ETA, ... 
  ###################
  if ( ! is.null(namesOri <- attr(X.pv,"namesOri"))) { ## includins NA's names (and etaFix$beta names)
    nc <- length(namesOri)
    beta_etaOri <- rep(NA,nc)
    names(beta_etaOri) <- namesOri
    beta_etaOri[names(beta_eta)] <- beta_eta ## keeps the original NA's
    beta_etaOri[names(etaFix$beta)] <- etaFix$beta  ## no longer in X.pv 2015/03
    res$fixef <- beta_etaOri ## FR->FR I should keep out the fixed ones for summary ? newetaFix code assumes the opposite
  } else {
    names(beta_eta) <- colnames(X.pv)
    res$fixef <- beta_eta
  }
  res$eta <- eta ## convenient for defining starting values...
  res$muetablob <- muetablob # for get_logdispObject, added 11/2016
  if (models[[1L]]=="etaHGLM" && HL[1]!="SEM") res$beta_cov <- get_from_MME(auglinmodblob$sXaug,"beta_cov")
  ###################
  ## LEVERAGES and REML (ie either phi OR lambda was estimated)
  ###################
  if (HL[1]!="SEM") { ## both lev_phi and deviance_residual missing otherwise
    if (is.null(phi.Fix) || anyNA(lambda.Fix)) { ## in either case all leverages are computed and it makes sense to consider the residuals
      res$lev_phi <- lev_phi
      dev_res <- family$dev.resids(y,mu,wt=1) * res$prior.weights
      res$std_dev_res <- sign(y-mu) * dev_res/(phi_est*(1-lev_phi)) ## should all have variance 1
    }
    if (anyNA(lambda.Fix)) res$lev_lambda <- lev_lambda
  }  
  res$distinctX.Re <- X.Re ## NULL if not distinct from X.pv
  ###################
  ## ALL other LAMBDA returns
  ###################
  res$rand.families <- rand.families 
  ##
  res$ranef <- structure(u_h,cum_n_u_h=cum_n_u_h) ## FR->FR added cum_n_u_h attribute 11/2014: slightly duplicates info in lambda object
  res$v_h <- v_h
  res$QRmethod <- processed$QRmethod
  ## FR->FR $w.ranef and $w.resid not doc'ed, as there is no mention of the augmented mode lin the doc.
  res$w.resid <- w.resid ## useful for .get_info_crits() and get_LSmatrix()
  if (models[["eta"]]=="etaHGLM") {
    res$w.ranef <- wranefblob$w.ranef ## useful for .get_info_crits() and get_LSmatrix()
    #
    res$lambda.object <- .make_lambda_object(nrand, lambda_models=models[["lambda"]], cum_n_u_h, lambda_est, 
                                            init.lambda, ## for attr(.,"type")
                                            process_resglm_blob, 
                                            coefficients_lambdaS, ## may come for process_resglm_blob OR from elsewhere
                                            ZAlist, next_LMatrices)
    res$"lambda" <- structure(unlist(res$lambda.object$lambda),cum_n_u_h=cum_n_u_h) ## redundant but very convenient
    if (attr(ZAlist,"anyRandomSlope")) {
      ## weird coding as next_LMatrices may have NULL elements;
      res$cov.mats <- lapply(next_LMatrices,function(mat) {
        if (is.null(mat)) {
          return(NULL)
        } else return(ZWZt(attr(mat,"Lcompact")
                           ,exp(coefficients_lambdaS[[which(attr(ZAlist,"Xi_cols")>1L)]]))) 
      })
    }
    ## building a more comprehensive list of descriptors of the structure of random effects
    strucList <- replicate(length(ZAlist), list(NULL)) 
    if (is.list(next_LMatrices)) { ## excludes primitive case of a single Matern matrix
      for (it in seq_len(length(ZAlist))) {
        lmatrix <- next_LMatrices[[it]]
        if ( ! is.null(lmatrix)) {
          corr.model <- attr(lmatrix, "corr.model") 
          if (is.null(corr.model)) warning('attr(next_LMatrix, "corr.model") is NULL')
          if (corr.model=="random-coef") {
            longLmatrix <- .makelong(attr(lmatrix,"Lcompact"),longsize=ncol(ZAlist[[it]]))
            ## extra attributes are hidden from str() because lmatrix is an S4 object...
            ## ranefs itself has attribute type="(.|.)"
            ## attr(lmatrix,"Lcompact") doit etre ce que .calc_latentL()$u retrouve
            lmatrix <- do.call(structure,c(list(.Data=longLmatrix),attributes(lmatrix)[c("Lcompact","par","ranefs","corr.model")]))
            attr(lmatrix,"cov.mat") <- ZWZt(attr(lmatrix,"Lcompact"),
                                            exp(coefficients_lambdaS[[which(attr(ZAlist,"Xi_cols")>1L)]])) 
          } 
          strucList[[it]] <- lmatrix  
        } 
      }
    }
    lmatrix <- processed$AUGI0_ZX$envir$LMatrix
    ## typically with attributes type="cholL_LLt", corr.model="Matern", ranefs
    if (!is.null(lmatrix)) {
      ## find ZAlist elements affected by LMatrix element
      affecteds <- which(attr(ZAlist,"ranefs") %in% attr(lmatrix,"ranefs"))
      for (it in affecteds) { strucList[[it]] <- lmatrix }
    }
    res$strucList <- strucList
  } ## else all these res$ elements are NULL
  ###################
  ## ALL other PHI returns
  ###################
  res$resid.predictor <- resid.predictor ## even if phi.Fix (04/2016), expected in summary of final hlcor call
  res$resid.family <- attr(processed$residModel$family,"quoted")  ## attribute used only for compact return
  # phi_est comes from calcPHIblob$next_phi_est, not from final glm,  hence is in minimal form
  if (models[["phi"]]=="phiScal") {res$phi <- phi_est[1]} else res$phi <- phi_est
  if (is.null(phi.Fix)) {
    if (models[["phi"]]=="phiHGLM") {
      res$resid_fit <- phifit
    } else {
      beta_phi <- calcPHIblob$beta_phi 
      names(beta_phi) <- unlist(lapply(names(beta_phi),function(st) {
        if (substr(st,1,1)=="X") {return(substring(st,2))} else {return(st)}
      })) ## removes "X" without guessing any order or length
      # FR->FR redundant info for summary, a nettoyer 
      phi.object <- list(fixef=beta_phi)
      phi.object$glm_phi <- calcPHIblob$glm_phi
      if (is.null(phi.object[["glm_phi"]])) {
        # delays computation of glm_phi
        glm_phi_args <- list(dev.res=dev_res*res$prior.weights,
                             control=processed$control.glm,
                             etastart=rep(calcPHIblob$beta_phi,nobs)) ## no glm <=> formula was ~1
        phi.object <- c(phi.object, list(glm_phi_args=glm_phi_args ) )
      } 
      res$phi.object <- phi.object
    }
  } else {
    ## important distinction for (summary, df of LRTs:
    if (is.null(processed$phi.Fix)) { ## absent from original call
      res$phi.object <- list(phi_outer=structure(phi.Fix,type="var")) ## hlcor call of corrHLfit / HLfit call post fitme ?
    } else res$phi.object <- list(phi_outer=structure(phi.Fix,type="fix"))
  }
  ################### the magic environment
  res$envir <- list2env(list(dvdloglamMat=dvdloglamMat, dvdlogphiMat=dvdlogphiMat), ## provided if available
                        parent=environment(HLfit_body))
  ###################
  ## WARNINGS
  ###################
  ## translation of warnings in user-more friendly form ##FR -> FR  a revoir
  if ( ! is.null(warningList$resLam0) && warningList$resLam0) { 
    warningList$resLam0 <- "lambda residuals numerically 0 were replaced by 1e-6"
  }
  if ( ! is.null(warningList$resLamInf) && warningList$resLamInf) { 
    warningList$resLamInf <- "lambda residuals numerically >1e10 were replaced by 1e10"
  }
  if (! is.null(warningList$leveLam1) && warningList$leveLam1) {
    warningList$leveLam1 <- "lambda leverages numerically 1 were replaced by 1 - 1e-8"
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
    warningList$innerPhiGLM <- paste("'",locw,"' in some sub-final iteration(s) of phi estimation;", sep="")
  }
  if (! is.null(locw <- warningList$innerLamGLM)) {
    warningList$innerLamGLM <- paste("'",locw,"' in some sub-final iteration(s) of lambda estimation;", sep="")
  }
  if ( HL[1]!="SEM" && maxit.mean>1 ## cases where iterations are needed 
      && ( ( models[[1]]=="etaHGLM"  && innerj==maxit.mean) 
           || 
           ( models[[1]]=="etaGLM" && pforpv>0 && innerj==maxit.mean)
        )) {
    warningList$innerNotConv <- paste("linear predictor estimation did not converge.\n",
                                      "Try increasing 'max.iter.mean' above ",maxit.mean,sep="")
  }
  if ( (! is.na(conv_logL)) && iter==max.iter) {
    maxitmess <- paste("Estimates did not converge. Try 'control.HLfit=list(LevenbergM=TRUE)';",
                       "\n else increase 'max.iter' above ",max.iter,
                       "\n (see help('HLfit') for details about 'max.iter')",sep="")
    if (models[["eta"]]=="etaHGLM") {
      if (conv_logL  && ! conv.lambda) {
        mainNotConv <- paste("p_v apparently converged but lambda estimates apparently did not.",
                             "\n This may indicate that some lambda estimate(s) should be zero.",
                             "\n Otherwise try increasing 'max.iter' above ",max.iter,
                             "\n (see help(HLfit) for details about 'max.iter')",sep="")          
      } else mainNotConv <- maxitmess        
      attr(mainNotConv,"diagnostics") <- c( conv.phi=conv.phi , conv.lambda=conv.lambda, 
                                            conv.corr=conv.corr, conv_lambda_vs_u=conv_lambda_vs_u,
                                            conv_rel_lambda=conv_rel_lambda )
    } else {
      mainNotConv <- maxitmess        
      attr(mainNotConv,"diagnostics") <- c( conv.phi=conv.phi , conv.lambda=conv.lambda, 
                                            conv.corr=conv.corr )
    }
    warningList$mainNotConv <- mainNotConv
  }
  res$warnings <- warningList
  res$spaMM.version <- packageVersion("spaMM")
  ###
  ###################
  ## SUMMARY, RETURN
  ###################
  class(res) <- c("HLfit",class(res)) 
  if (verbose["summary"]) {
    summary(res) 
  }
  if (verbose["warn"]) {
    seriousWarnings <- warningList[intersect(c("innerNotConv","mainNotConv"),names(warningList))]
    if (length(seriousWarnings)>0 ) { 
      abyss <- sapply(length(seriousWarnings),function(i) {
        warning(paste("In HLfit :\n",seriousWarnings[[i]],sep=""),call.=FALSE)}) 
      warningList[setdiff(names(warningList),c("innerNotConv","mainNotCov"))] <- NULL
    }
  } 
  if (verbose["trace"]) {
    if (length(warningList)>0 ) {
      abyss <- sapply(length(warningList),function(i) {cat(warningList[[i]]);cat("\n")}) 
    }
  }
  # cleanup: for diagnostic, use
  # sort(sapply(ls(<object>$envir), function(x)
  # +             object.size(get(x, envir = <object>$envir))),decreasing=TRUE)
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"res")) ## empties the whole local envir except the return value
  return(res)
}