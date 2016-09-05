HLfit_body <- function(processed, #resid.model= ~ 1, 
                       verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),
                       control.HLfit=list(), ## used both by preprocess and HLfit_body
                       init.HLfit = list(), ## apparently not used by preprocess
                       ranFix=list(), ## phi, lambda, possibly nu, rho if not in init.HLfit
                       etaFix=list() ## beta, v_h (or even u_h)
) {
  #mc <- match.call() 
  data <- processed$data
  family <- processed$family
  if (family$family=="COMPoisson") {
    if ( ! is.null(ranFix$COMP_nu)) { ## optimisation call
      assign("nu",ranFix$COMP_nu,envir=environment(family$aic))
      ranFix$COMP_nu <- NULL
    } else {
      checknu <- try(environment(family$aic)$nu,silent=TRUE)
      if (inherits(checknu,"try-error")) stop(attr(checknu,"condition")$message)
    }
  } else if (family$family=="negbin") {
    if ( ! is.null(ranFix$NB_shape)) { ## optimisation call
      assign("shape",ranFix$NB_shape,envir=environment(family$aic))
      ranFix$NB_shape <- NULL
    } else {
      checktheta <- try(environment(family$aic)$shape,silent=TRUE)
      if (inherits(checktheta,"try-error")) stop(attr(checktheta,"condition")$message)
    }
  }
  
  formula <- processed$predictor
  prior.weights <- processed$prior.weights
  
  rpblob <- canonizeRanPars(ranPars=ranFix,corr.model="",checkComplete = FALSE) 
  ranFix <- rpblob$ranPars
  
  corr_est <- init.HLfit[intersect(c("nu","rho","Nugget","ARphi"),names(init.HLfit))]
  if (length(corr_est)==0L) corr_est <- NULL
  
  warningList <- list()
  ## whene addingverbose elements, remind that these might be lost through corrHLfit -> HLCor cf dotlist$verbose <- verbose[intersect(...]
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
  if (HL[1]=="SEM") SEMargs <- processed$SEMargs 
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
  vec_n_u_h <- attr(cum_n_u_h,"vec_n_u_h")
  if (models[["eta"]]=="etaHGLM") {
    LMatrix <- attr(predictor,"LMatrix")
    TT <- processed$AUGXZ_0I 
    if ( ! is.null(LMatrix) ) { ## then LMatrix is presumably variable and TT was precomputed without it
      ZALlist <- computeZAXlist(XMatrix=LMatrix,ZAlist=processed$ZAlist)
      userLfixeds <-attr(ZALlist,"userLfixeds") # for MakeCovEst, set to (a list of) FALSE
      ZAL <- post.process.ZALlist(ZALlist,as_matrix=(processed$QRmethod=="Matrix::qr")) ## may be modified by other call to post.process.ZALlist   
      TT[1:nobs,(pforpv+1L):ncol(TT)] <- ZAL # TT <- calcTT(X001=Xpv001,ZAL) 
      # now we have TT
    } else { ## ZAL was  accessible to preprocess; it must be in attr(predictor,"ZALMatrix")
      ZAL <- processed$ZA
      ## I could have stored ZAL <- attr(predictor,"ZALMatrix") in processed$ZA ut let's keep the distinction so far
      if (is.null(ZAL)) ZAL <- attr(predictor,"ZALMatrix")
      userLfixeds <- FALSE 
    }
    tZAL <- as.matrix(t(ZAL))## used by auglinmodfit
  } else if (models[[1]]=="etaGLM") {
    ZAL <- NULL ## 
    ZALlist <- NULL
    u_h <- v_h <- lev_lambda <- numeric(0)
    TT <- X.pv
  }
  ### a bit of post processing // repeat of code in preprocess...
  lambda.Fix <- getPar(ranFix,"lambda") 
  if (is.null(lambda.Fix)) lambda.Fix <- rep(NA,nrand)
  if (any(lambda.Fix[!is.na(lambda.Fix)]==0)) stop(pastefrom("lambda cannot be fixed to 0.",prefix="(!) From "))
  ###
  off <- attr(processed$predictor,"offsetObj")$total
  next_cov12_est <- NULL ## will be tested
  ##################
  if (is.character(init.HLfit)) { ## at this point init.HLfit is a string or not. Elsewhere it can be a list
    spaMM.options(INIT.HLFITNAME=init.HLfit) ## if a string, copied in...
  } else {
    spaMM.options(INIT.HLFITNAME=NA)  
    # init.HLfitName <- NULL
    unknowns <- names(init.HLfit)[!names(init.HLfit) %in% c("fixef","phi","lambda","v_h","rho","nu","Nugget","ARphi")] 
    if (length(unknowns)>0) {
      mess <- pastefrom("unhandled elements in 'init.HLfit'.",prefix="(!) From ")
      message(mess)
      if ("beta" %in% unknowns) message("  Use 'fixef' rather than 'beta' in 'init.HLfit'.")
      stop()
    }
  }
  ###################
  if ( ! is.null(corr_est)) {
    corrEstBlob <- eval.corrEst.args(family=family,rand.families=rand.families,predictor=predictor,data=data,X.Re=X.Re,
                                     REMLformula=REMLformula,ranFix=ranFix,
                                     Optimizer=control.HLfit$Optimizer)
    corrEst.args <- corrEstBlob$corrEst.args ## but corrEstBlob also has $corrEst.form which will stay there for later use
  }

  ### case where nothing to fit #############################################
  if (is.null(corr_est) && 
      ncol(X.pv)==0L &&
      !is.null(phi.Fix) &&
      (models[[1]]=="etaGLM" || (!is.null(etaFix$v_h) &&  ! anyNA(lambda.Fix))) 
  ) { ## nothing to fit. We just want a likelihood
    ### a bit the same as max.iter<1 ... ?
    phi_est <- phi.Fix
    eta <- off
    if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
      ## we need u_h in calc.p_v() and v_h here for eta...
      v_h <- etaFix$v_h
      u_h <- etaFix$u_h
      if (is.null(u_h)) {u_h <- u_h_v_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,lower.v_h=NULL,upper.v_h=NULL)}
      lambda_est <- resize_lambda(lambda.Fix,vec_n_u_h,n_u_h)
      eta <- eta + drop(ZAL  %id*id%  etaFix$v_h) ## updated at each iteration
    } ## FREQS
    ## conversion to mean of response variable (COUNTS for binomial)
    muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
    w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
      wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## no fit, likelihood computation
      d2hdv2 <- calcD2hDv2(ZAL,w.resid,wranefblob$w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    }
    return(list(APHLs=calc.p_v(mu=muetablob$mu, u_h=u_h, dvdu=wranefblob$dvdu,
                               lambda_est=lambda_est, phi_est=phi_est, d2hdv2=d2hdv2, processed=processed))) 
    ### RETURN !! ## FR->FR but p_bv is not returned.
  } 
  
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################

  ### Initial values  for lambda, phi and beta from lambda.Fix, phi.Fix, or init.HLfit ##### 
  ## Initial estimate for beta
  beta_eta <- numeric(pforpv)
  if (pforpv>0) beta_eta <- init.HLfit$fixef
  ## Initial estimate for phi 
  phi_est <- phi.Fix
  if (is.null(phi.Fix)) phi_est <- init.HLfit$phi ## must be a vector of 'response' values
  ## Initial estimate for lambda 
  if (models[[1]]=="etaHGLM") { 
    gaussian_u_ranges <- processed$gaussian_u_ranges 
    psi_M <- rep(attr(rand.families,"unique.psi_M"),attr(cum_n_u_h,"vec_n_u_h"))
    v_h <- initialize.v_h(psi_M=psi_M,etaFix=etaFix,init.HLfit=init.HLfit,cum_n_u_h=cum_n_u_h,rand.families=rand.families)
    u_h <- u_h_v_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,lower.v_h=NULL,upper.v_h=NULL)
    init.lambda <- lambda.Fix ## already the right size 
    type <- rep("",nrand) ## actually fix or outer
    whichInner <- which(is.na(lambda.Fix))
    type[whichInner] <- "inner" ## sure
    type[type=="" & is.na(processed$lambda.Fix)] <- "outer" ## sure
    type[type==""] <- "fix" ## sure
    if (! is.null(init.HLfit$lambda)) init.lambda[whichInner] <- init.HLfit$lambda[whichInner]
  }
  ###
  ######### missing Initial estimates for mu, phi, lambda by GLM ####################
  if (  is.null(beta_eta) || is.null(phi_est) || anyNA(init.lambda) ) { 
    reset <- (family$family %in% c("COMPoisson","negbin"))
    inits_by_glm <- processed$get_inits_by_glm(processed,family=family,reset=reset) # using possibly postprocessed family
    ## : uses processed$y, $BinomialDen, attr($predictor,"offsetObj"), $control.glm.
  }
  ## Finalize initial values for beta_eta
  if (is.null(beta_eta) ) beta_eta <- inits_by_glm$beta_eta 
  if (!is.null(control.HLfit$intervalInfo)) {
    parmcol <- attr(control.HLfit$intervalInfo$parm,"col")
    beta_eta[parmcol] <- control.HLfit$intervalInfo$init 
  }  
  ## Finalize initial values for phi 
  if (is.null(phi_est) ) {
    ## there are many cases that get_init_phi does not handle and then it may have not been created
    if ( ! is.null(processed$get_init_phi)) {
      if (is.call(processed$prior.weights)) {
        phi_est <- processed$get_init_phi(processed,weights=NULL)
      } else phi_est <- processed$get_init_phi(processed,weights=processed$prior.weights)
    }
    if ( is.null(phi_est) ) phi_est <- inits_by_glm$phi_est
    if (models[["phi"]] != "phiScal") {
      phi_est <- rep(phi_est,nobs) ## moche ## FR->FR why is this necess ?
    }
  }
  ## Finalize initial values for lambda
  if (models[[1]]=="etaHGLM") { ## the basic case (LMM, GLMM...)
    if ( anyNA(init.lambda) ) { ## some inner may remain undetermined
      stillNAs <- which(is.na(init.lambda))
      ## if reset is TRUE init.lambda is recomputed. Otherwise it is'get' ie recomputed only if not previously computed 
      default.init.lambda <- processed$get_init_lambda(processed, reset=reset, stillNAs=stillNAs,
                                                       init_lambda_by_glm=inits_by_glm$lambda)
      init.lambda[stillNAs] <- default.init.lambda[stillNAs]
    } 
    attr(init.lambda,"type") <- type
    lambda_est <- resize_lambda(init.lambda,vec_n_u_h,n_u_h)
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
  conv.phi <- FALSE; conv.lambda <- FALSE; conv.corr <- FALSE
  if (models[[1]]=="etaHGLM") {
    wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## initilization !
    ddi_or_matrix_ZAL <- post_process_ZAL(ZAL,attr(w.resid,"unique"))
    d2hdv2 <- calcD2hDv2(ddi_or_matrix_ZAL,w.resid,wranefblob$w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    if (ncol(X.pv)==0L && !is.null(etaFix$v_h)) {
      maxit.mean <- 0 ## used in test near the end...
    } else if ( LMMbool && is.null(control.HLfit$intervalInfo) ) {
      maxit.mean <- 1 ## 1 is sufficient for LMM as Hessian does not vary with beta_eta  => quadratic function;
      # ./ E X C E P T for calculation of confidence intervals: at least two intervalSteps are required. Quite visible when 
      # dispersion params areouter-estimated (fitme) in which case there isno outer iteration to compensate a small maxit.mean  
    } else { ## even h maximization in *G*LMMs 
      if ( ! is.null(phi.Fix) && ! anyNA(lambda.Fix)) { ## allFix hence no true outer iteration 
        maxit.mean <- iter.mean.dispFix 
      } else maxit.mean <- iter.mean.dispVar # If phi.Fix and lambda.Fix, the only way to have 'outer' convergence is to have 'inner' convergence
    } 
    next_LMatrices <- LMatrix ## next_LMatrices originally for random slopes, originally <- NULL, modified 2015/06/03, cf notes of that day
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
    SEMargs$qr.XtX <- QRwrap(crossprodCpp(X.pv),useEigen=FALSE) ## qr(t(X.pv)%*%X.pv) ## pas sur que FALSE gagne du temps
    SEMargs$beta_eta <- beta_eta
    SEMargs$corr_est <- corr_est["rho"] ## may be NULL depending in particular on init.HLfit
    SEMargs$ZA <- processed$ZAlist[[1]]
    SEMargs$lambda <- init.lambda ## unique(lambda_est)
    SEMargs$lambda.Fix <- lambda.Fix ## may be NA
    if (is.null(LMatrix)) {
      locdim <- ncol(SEMargs$ZA)
      SEMargs$symSVD <- list(corr.model="identity",
                             symsvd=list(u=Diagonal(n=locdim),
                                         ## diagonal matrix (ddiMatrix) with @diag = "U"
                                         d=rep(1,locdim)), 
                             dim=rep(locdim,2)
      )
    } else SEMargs$symSVD <- attributes(LMatrix) ## includes dim(LMAtrix)
    SEMargs$ZAL <- ZAL ## FR->FR should be ddi_or_matrix_ZAL ?
    SEMargs$off <- off
    #if (SEMargs$SEMlogL=="p_v") SEMargs$mc <- mc ## pass HLfit call args
    SEMargs$stop.on.error <- stop.on.error
    ## following two lines may not go in preprocess if X.pv is modified by HLfit
    SEMargs$X.pv <- X.pv
    SEMargs$X_lamres <- processed$X_lamres
    SEMargs$whichy1 <- (y==1) ##FR->FR in preprocess ?
    SEMargs$verbose <- verbose["SEM"]
    SEMargs$control.glm <- processed$control.glm
    ##
    SEMblob <- do.call("SEMbetalambda",SEMargs)  ########## CALL
    beta_eta <- SEMblob$beta_eta
    lambda_est <- resize_lambda(SEMblob$lambda,vec_n_u_h,n_u_h)
    corr_est["rho"] <- SEMblob$corr_est["rho"] ## may again be NULL
    u_h <- v_h <- SEMblob$v_h
    logLapp <- SEMblob$logLapp
    attr(logLapp,"seInt") <- SEMblob$seInt ## may be NULL
    ## for calc_beta_cov:
    eta <- off + drop(X.pv %*% beta_eta) + drop(ZAL %id*% v_h)
    muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
    w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## no fit, likelihood computation
    sqrt.ww <- sqrt(c(w.resid,wranefblob$w.ranef))
    wAugX <- calc_wAugX(augX=TT,sqrt.ww=sqrt.ww)
    ##
  } else while ( TRUE ) { ##the main loop with steps: new linear predictor, new leverages, new phi, new w.resid, new lambda, new fn(lambda)
    if (models[[1]]=="etaHGLM") {
      ##############################
      auglinmodblob <- auglinmodfit(TT=TT,ZAL=ZAL,
                                    tZAL=tZAL,
                                    lambda_est=lambda_est,wranefblob=wranefblob,
                                    d2hdv2=d2hdv2,w.resid=w.resid,beta_eta=beta_eta,
                                    maxit.mean=maxit.mean,eta=eta,u_h=u_h,v_h=v_h,
                                    control.HLfit=control.HLfit,
                                    X.pv=X.pv,etaFix=etaFix,
                                    cum_n_u_h=cum_n_u_h,psi_M=psi_M,
                                    muetablob=muetablob,
                                    phi_est=phi_est,verbose=verbose,
                                    ranFix=ranFix,
                                    corr_est=corr_est, ## only for messages
                                    processed=processed,
                                    ZALtZAL=NULL,
                                    ddi_or_matrix_ZAL=ddi_or_matrix_ZAL
      ) ## HL(.,.) estim of beta, v for given lambda,phi
      ##############################
      beta_eta <- auglinmodblob$beta_eta
      v_h <- auglinmodblob$v_h
      u_h <- auglinmodblob$u_h
      eta <- auglinmodblob$eta
      wranefblob <- auglinmodblob$wranefblob
      muetablob <- auglinmodblob$muetablob
      mu <- muetablob$mu ## testé par accident, necess dans test COMPoisson HLfit...
      w.resid <- auglinmodblob$w.resid
      d2hdv2 <- auglinmodblob$d2hdv2
      wAugX <- auglinmodblob$wAugX
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
      if (anyNA(lambda.Fix) || is.null(phi.Fix)) {
        if (maxit.mean==0L) {
          stop("(!) Computation of leverages with maxit.mean=0: check that this is meaningful.")
        } # ELSE rWW was updated in the inner loop for betaV
        hatval <- calc_hatval(X.Re=X.Re,
                              sqrt.ww=auglinmodblob$sqrt.ww,
                              wAugX=as.matrix(wAugX), ## promise...
                              levQ=auglinmodblob$levQ)
        if (any(abs(hatval) > 1 - 1e-8)) {
          hatval <- pmin(hatval,1-1e-8) # ifelse(abs(hatval) > 1 - 1e-8, 1 - 1e-8,hatval)
          warningList$leveLam1 <-TRUE
        }
        lev_phi <- hatval[1L:nobs] ## for the error residuals (phi)
        lev_lambda <- hatval[(nobs+1L):(nobs+cum_n_u_h[nrand+1L])]  ## for the ranef residuals (lambda)
      }
    } else { ## GLM fitted by ML. Here I had more ad hoc code up to v1.7.42
      if ( ! is.null(X.Re) ) { 
        wAugXleve <- calc_wAugX(augX=X.Re,sqrt.ww=sqrt(w.resid)) # rWW%*%X.Re 
        if (is.null(phi.Fix)) lev_phi <- leverages(wAugXleve)
      } else { ## basic REML, leverages from the same matrix used for estimation of beta
        wAugX <- calc_wAugX(augX=X.pv,sqrt.ww=sqrt(w.resid)) # rWW %*% X.pv 
        if (is.null(phi.Fix)) lev_phi <- leverages(wAugX)
      }
    }
    #### contribution from GLM weights
    if (HL[2]>0 && models[[1]]=="etaHGLM" 
        && (anyNA(lambda.Fix) || is.null(phi.Fix)) ) { ## LeeN01 HL(.,1) ie the + in 'EQL+'
      ## (0): previous hat matrix -> p, notEQL -> tilde(p), (1): full correction -> q 
      ## first the d log hessian / d log lambda or phi corrections then, IF HL[3]>0, the notEQL correction
      ## For the d log hessian first the derivatives of GLM weights wrt eta 
      ##################### noter que c'est le coef2 de HL(1,.), but mu,eta may have been updated since coef2 was computed
      dlW_deta <- calc.dlW_deta(dmudeta=muetablob$dmudeta,family=family,mu=mu,eta=eta,
                                BinomialDen=BinomialDen,canonicalLink=processed$canonicalLink)$dlW_deta
      ## we join this with the deriv of log w.ranef wrt v_h
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
          dvdloglamMat <-  calc.dvdloglamMat(dlogfthdth=(psi_M - u_h)/lambda_est, ## the d log density of th(u)
                                             cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,rand.families=rand.families,
                                             u_h=u_h,d2hdv2=d2hdv2,stop.on.error=stop.on.error)
          dleve <- leve__dlW_deta_or_v__ZALI %*% dvdloglamMat # (r+n).(r+n)Xr.rXr = r (each element is a sum over r+n terms= a trace)
          lev_lambda <- lev_lambda - as.vector(dleve)  
        } 
        ## 
        if (is.null(phi.Fix)) {
          dh0deta <- ( w.resid *(y-mu)/muetablob$dmudeta ) ## 12/2013 supp BinomialDen (soit Bin -> phi fixe=1, soit BinomialDen=1)
          dvdlogphiMat <- calc_dvdlogphiMat(dh0deta=dh0deta,ZAL=ZAL,d2hdv2=d2hdv2,stop.on.error=stop.on.error)
          dleve <- leve__dlW_deta_or_v__ZALI %*% dvdlogphiMat # (r+n) . (r+n)Xr . rXn = n (each element is a sum over r+n terms= a trace)
          lev_phi <- lev_phi - as.vector(dleve)  
        } 
      }
    }
    if (HL[2]>1) {stop("Need a_i correction in Table 7 of NohL07 ie derivatives of second order correction wrt dips param.")}
    #### contribution from exact likelihood function instead of EQL
    if (HL[3]!=0 ) {## HL(.,.,1) ie , p_bv(h), not EQL p_bv(q+), LeeNP p89; distinction does not arise for PQL <=> Gaussian ranefs...  
      # lambda
      if (models[[1]]=="etaHGLM" && anyNA(lambda.Fix)) ## d h/ d !log! lambda correction     
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
      lev_phi <- pmin(lev_phi, 1 - 1e-8)
      ## updated residuals from updated mu must be used (LeeNP p.161) [not so in dhglmfit !!]
      ## Once in Gamma GLM, y=25.75563, y-mu=-5.996906e-08, yielded negative dev.resid
      if (models[["phi"]]=="phiHGLM") {
        residProcessed$prior.weights <- structure((1-lev_phi)/2,unique=FALSE) # expected structure in processed.
        # uses of prior weights matches thatin input of calcPHI -> dispGammaGLM 
        residProcessed$data$.phi <- family$dev.resids(y,mu,wt= eval(prior.weights))/(1-lev_phi) 
        residProcessed$y <- residProcessed$data$.phi
        #cat("\nbefore")
        #browser()
        if (iter==0) {
          initphifit <- list()
          if (is.null(processed$residModel$fixed$phi)) initphifit$phi <- 1 
          ## FR->FR Fixing phi makes a difference in crack example !
        } else {
          initphifit <- phifit$corrPars[attr(phifit$corrPars,"type")=="outer"] ## may be an empty list
          if (length(initphifit)>0L) {
            ## FR->FR init must be in user scale, optr output is 'unreliable' => complex code
            LUarglist <- attr(phifit,"optimInfo")$LUarglist ## might be NULL 
            if (! is.null(LUarglist) && LUarglist$optim.scale=="transformed") {
              LowUp <- do.call("makeLowerUpper",LUarglist)
              if ( !is.null(initphifit$nu) ) {
                NUMAX <- LUarglist$NUMAX
                initnu <- nuFn(initphifit$nu,NUMAX=NUMAX) 
                initnu <- initnu +1e-4 *( (LowUp$lower$trNu+LowUp$upper$trNu)/2 -initnu)
                initphifit$nu <- nuInv(initnu,NUMAX=NUMAX)
              }
              if ( !is.null(initphifit$rho) ) {
                RHOMAX <- LUarglist$RHOMAX
                initrho <- rhoFn(initphifit$rho,RHOMAX=RHOMAX) 
                initrho <- initrho +1e-4 *( (LowUp$lower$trRho+LowUp$upper$trRho)/2 -initrho)
                initphifit$rho <- rhoInv(initrho,RHOMAX=RHOMAX)
              }
            }                      
          }
          if (all(is.na(residProcessed$lambda.Fix))) initphifit$lambda <- phifit$lambda
          if (is.null(processed$residModel$fixed$phi)) initphifit$phi <- 1 
        }
        phifitarglist <- processed$residModel
        ## which may include processed$residModel$fixed, which is a *list* created by preprocess
        phifitarglist$processed <- residProcessed 
        ## A resid.model's fixed lambda, phi have been preprocessed 
        ##   and info in processed$residModel$fixed must match that in residProcessed
        phifitarglist$init <- initphifit
        phifitarglist$verbose <- reformat_verbose(NULL,For="corrHLfit") ## constrained
        phifit <- do.call("fitme_body",phifitarglist)
        #summary(phifit)
        prevmsglength <- overcat(paste("phi fit in iter=",iter+1L,
                                       ", .phi[1]=",signif(phifit$y[1],5),", ",
                                       paste(c(names(phifit$corrPars),"lambda"),"=",
                                               signif(c(unlist(phifit["corrPars"]),phifit$lambda),6),
                                               collapse=", ",sep=""),
                                       "; parent: ",sep=""),
                                 prevmsglength)
        ## the final printing after iterations comes from HLfit.obj/HLCor.obj -> if (mc$verbose["objective"]) ...
        next_phi_est <- phifit$fv
      } else {
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
    ######### Dispersion Estimates for lambda #####################
    if (models[[1]]=="etaHGLM" && anyNA(lambda.Fix)) { ## lambda must be estimated 
      if (any(abs(lev_lambda) > 1 - 1e-8)) { ## abs... not commented when written...
        lev_lambda <- pmin(lev_lambda,1-1e-8) #ifelse(abs(lev_lambda) > 1 - 1e-8, 1 - 1e-8,lev_lambda)
        warningList$leveLam1 <- TRUE
      }      
      cov12_est <- next_cov12_est
      ##################
      ranefEstargs <- list(u_h=u_h,ZAlist=processed$ZAlist,cum_n_u_h=cum_n_u_h,
                           prev_LMatrices=next_LMatrices,
                           userLfixeds=userLfixeds,
                           w.resid=w.resid,
                           processed=processed)
      
      if (attr(processed$ZAlist,"anyRandomSlope")) { ## if random-slope model
        covEstmethod <- .spaMM.data$options$covEstmethod ## note same call in calcRanefPars
        if (is.null(covEstmethod)) stop("spaMM.getOption('covEstmethod') should not be NULL")
        if (covEstmethod == "makeCovEst1") {## le plus bourrin
          ranefEstargs <- c(ranefEstargs,list(phi_est=phi_est,
                                              locTT=TT,v_h=v_h))
          ranefEstargs$auglinfixedpars <- list(d2hdv2=d2hdv2,w.resid=w.resid,beta_eta=beta_eta,
                                               maxit.mean=maxit.mean,eta=eta,u_h=u_h,v_h=v_h,
                                               control.HLfit=control.HLfit,
                                               X.pv=X.pv,etaFix=etaFix,
                                               cum_n_u_h=cum_n_u_h,psi_M=psi_M,
                                               muetablob=muetablob,
                                               phi_est=phi_est,verbose=verbose,
                                               ranFix=ranFix,
                                               corr_est=corr_est, 
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
      next_cov12_est <- calcRanefPars_blob$next_corrEstList$cov12_est
      next_LMatrices <- calcRanefPars_blob$next_LMatrices ## random slope: must converge to the L factor of corr mat
      next_lambda_est <- calcRanefPars_blob$next_lambda_est 
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
    ############# experimental spatial model estimation
    if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="Matern") { ## this code does not apply for the random slope model
      corrEst.args$formula <- Predictor(formula=corrEstBlob$corrEst.form,offset=off + X.pv %*% beta_eta) ## FR->FR ugly input for offset
      corrEst.args$init.corrHLfit <- corr_est ## this entails use of optim() (or another Optimizer) on these parameters
      if (nrand>1) stop("code needed for corr Estim within HLfit with multiple lambda parameters") ## FR->FR
      corrEst.args$ranFix$lambda <- unique(lambda_est)
      corrEst.args$ranFix$phi <- phi_est 
      corrEst.args$init.HLfit$v_h <- v_h  
      ## corrEst.args$HLmethod <- .... ## default REML  ~ ML here
      #       if (FALSE) { ## seems to work...
      #       locprocessed <- preprocess(control.HLfit=control.HLfit,HLmethod=HLmethod,
      #                                  predictor=Predictor(formula=corrEst.form,offset=off + X.pv %*% beta_eta),phi.Fix=phi_est,                 
      #                                  resid.predictor=resid.formula, ## must be ignored, but no default... =>preprocess could be improved
      #                                  REMLformula=corrEst.args$REMLformula,data=data,
      #                                  family=family,BinomialDen=BinomialDen,rand.family=rand.family)
      #       corrEst.args$processed <- locprocessed ## risky
      #       }
      pff <- do.call("corrHLfit",corrEst.args)
      next_corr_est <- pff$corrPars[names(corr_est)] ## rho,nu,  pas trRho, trNu 
      #FR->FR maybe conv_threshold a bit strict here...
      if (all(abs(log(unlist(next_corr_est)/unlist(corr_est))) < conv.threshold) ) { ## 
        conv.corr <- TRUE ## this is the simplest, best case. ## but if slow geometric decrease to 0, this is never true 
      } else if (all(abs(unlist(next_corr_est)-unlist(corr_est)) < conv.threshold* (unlist(corr_est)+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.corr <- TRUE ## 'weak convergence'... 
      } else conv.corr <- FALSE
    } else {
      if (!is.null(next_cov12_est)) {
        if (iter>1 && abs(cov12_est-next_cov12_est) < conv.threshold ) { 
          conv.corr <- TRUE 
        } else conv.corr - FALSE       
      } else conv.corr <- TRUE
    }
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
          ZALlist <- computeZAXlist(XMatrix=next_LMatrices,ZAlist=processed$ZAlist)
        } else if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="Matern") {
          corr_est <- next_corr_est 
          LMatrix <- attr(pff$predictor,"LMatrix")
          ZALlist <- computeZAXlist(XMatrix=LMatrix,ZAlist=processed$ZAlist)
        } else if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="adjacency") {
          corr_est <- next_corr_est 
          #LMatrix <- attr(pff$predictor,"LMatrix") ## FR->FR mais LMatrix devrait simplement être decomp$u
          #ZALlist <- computeZAXlist(XMatrix=LMatrix,ZAlist=processed$ZAlist)
        } else if (! is.null(corr_est) && attr(LMatrix,"corr.model")=="SAR_WWt") {
          corr_est <- next_corr_est ## 
          ## 
        }
        ZAL <- post.process.ZALlist(ZALlist,as_matrix=(processed$QRmethod=="Matrix::qr")) 
        ddi_or_matrix_ZAL <- post_process_ZAL(ZAL,attr(w.resid,"unique"))
        TT[1:nobs,(pforpv+1L):ncol(TT)] <- ZAL ## TT <- calcTT(X001=Xpv001,ZAL) 
      }
      if (models[[1]]=="etaHGLM" && anyNA(lambda.Fix)) { ## lambda was modified
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
        wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h)
      } 
      if (models[[1]]=="etaHGLM") {
        if (anyNA(lambda.Fix) || is.null(phi.Fix) || ! is.null(corr_est)) { ## w.ranef or w.resid or ZAL were modified 
          d2hdv2 <- calcD2hDv2(ddi_or_matrix_ZAL, w.resid, wranefblob$w.ranef) ## - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
        }
      } 
      ## conv_logL either used to break the loop, Xor required only in last two iters for diagnostics 
      if (processed$break_conv_logL 
          || ( verbose["trace"] && iter>=max.iter-1L )) {
        next_lik <- calc.p_v(mu=mu, u_h=u_h, dvdu=wranefblob$dvdu, 
                             lambda_est=lambda_est, phi_est=phi_est, d2hdv2=d2hdv2,processed=processed)$p_v ## no p_bv...
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
    ## FR->FR tries to save time by using attr(d2hdv2,"APHLs"), but hisis messy
    if (models[[1]]=="etaHGLM") {
      if (pforpv==0L) { 
        d2hdv2 <- calcD2hDv2(ddi_or_matrix_ZAL, w.resid, wranefblob$w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
#       APHLs <- calc.p_v(mu=mu, u_h=u_h, dvdu=wranefblob$dvdu, 
#                          lambda_est=lambda_est, phi_est=phi_est, d2hdv2=d2hdv2,processed=processed)
      } #else APHLs <- attr(d2hdv2,"APHLs") ## non-null if auglinmodfit provided it
#      if (is.null(APHLs$lad)) { ## either APHLs is null or it contains hlik but not p_v or it has both p_v and lad
        APHLs <- calc.p_v(mu=mu, u_h=u_h, dvdu=wranefblob$dvdu, 
                          lambda_est=lambda_est, phi_est=phi_est, d2hdv2=d2hdv2,processed=processed)
#      }
    } else APHLs <- calc.p_v(mu=mu, u_h=u_h, dvdu=wranefblob$dvdu, 
                        lambda_est=lambda_est, phi_est=phi_est, d2hdv2=d2hdv2,processed=processed) ## GLM, d2hdv2 missing
  }
  ######################### potential R E T U R N here: cases without p_bv
  if ( identical(processed$return_only,"p_vAPHLs")) {
    # a bit of ugly coding, but optimthroughsmooth calls HLCor, not HLCor.obj, thus it cannot directly control return_only. So either leave as is, or move the test to HLCor, or modify optimthroughsmooth to call HLCor.obj  
    if (HL[1]=="SEM") { # lambda used for smoothing.
      RESU <- list(APHLs=APHLs,lambda=SEMblob$lambda)    ########################   R E T U R N
    } else RESU <- list(APHLs=APHLs)    ########################   R E T U R N
    ### voir verbose["objective"]
    # if (verbose["hack"]) {
    #   print(paste(paste(unlist(RESU$APHLs),collapse=", "), 
    #               paste(unique(lambda_est),collapse=", ")))
    # }
    return(RESU)   ########################   R E T U R N
  } ## ELSE continue: compute p_bv
  #
  ## computation p_bv ## cf notes 19/08/2016 pour alcul APHLs et IC's for phiHGLM 
  if (HL[1] != "SEM") {
    ## ML: X.Re non NULL mais ncol(X.Re)=0
    if ( is.null(X.Re) ||  ncol(X.Re)>0 ) {
      ## REML standard || REML non standard
      if (is.null(X.Re)) {X.REML <- X.pv} else X.REML <- X.Re
      Md2hdb2 <- ZtWZwrapper(X.REML,w.resid)
      if (models[[1]]=="etaHGLM") {
        hessnondiag <- suppressWarnings(crossprod(ZAL, sweep(X.REML, MARGIN = 1, w.resid, `*`))) ## Matrix or matrix depending on ZAL
        Md2hdbv2 <- rbind2(cbind2(Md2hdb2, t(hessnondiag)),
                           cbind2(hessnondiag, - d2hdv2)) 
        ladbv <- LogAbsDetWrap(Md2hdbv2,logfac=-log(2*pi))
      } else { ## GLM
        ladb <- LogAbsDetWrap(Md2hdb2,logfac=-log(2*pi))
      }
    } else { 
      ## fit ML standard : ncol(X.Re)=0, p_bv=p_v hence d2hdpbv reduces to d2hdv2
      if (models[[1]]=="etaHGLM") {
        Md2hdbv2 <- - d2hdv2 
        ladbv <- APHLs$lad
      } else ladb <- 0
    }
    if (models[[1]]=="etaHGLM") {
      APHLs$p_bv <- APHLs$hlik-ladbv/2  
    } else  APHLs <- list(p_v=APHLs$clik, 
                          p_bv=APHLs$clik - ladb/2) ## 07/2016 for inference about phi in GLMs 
    if ( ! is.null(APHLs$second.corr)) APHLs$p_bv <- APHLs$p_bv + APHLs$second.corr
  }
  ######################### potential R E T U R N here: with p_bv
  if ( identical(processed$return_only,"p_bvAPHLs")) {
    ### voir verbose["objective"]
    # if (verbose["hack"]) {
    #   print(paste(paste(unlist(APHLs),collapse=", "), 
    #               paste(unique(lambda_est),collapse=", ")))
    # }
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
  if (models[[1]]=="etaHGLM" && anyNA(lambda.Fix)) {
    bloc_lambda_args <- list(models=models, init.lambda=init.lambda, processed=processed, lambda.Fix=lambda.Fix, 
                             cum_n_u_h=cum_n_u_h, next_LMatrices=next_LMatrices)
    if (HL[1]=="SEM") {
      bloc_lambda_args$SEMblob <- SEMblob
    } else {
      bloc_lambda_args$calcRanefPars_blob <- calcRanefPars_blob
      bloc_lambda_args$lev_lambda <- lev_lambda
    }
    #
    process_resglm_blob <- do.call("bloc_lambda",bloc_lambda_args)
    #
    coefficients_lambdaS <- process_resglm_blob$coefficients_lambdaS # list
    p_lambda <- length( unlist(coefficients_lambdaS[ ! attr(init.lambda,"type")=="fix"])) 
    if (attr(processed$ZAlist,"anyRandomSlope")) {
      ## lignes suiv supposent que L_matrix decrit random slope model
      p_corr <- sum(unlist(lapply(next_LMatrices,function(mat) {
        dimL <- nrow(attr(mat,"Lcompact"))
        (dimL-1)*dimL/2
      })))
      p_lambda <- p_lambda+p_corr
    }
  } else p_lambda <- 0    
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
  res$corrPars <- get_DispCorrPars(ranFix, corr_est) 
  ## FR->FR il serait logique ? de regrouper $ranFix et $corrPars dans la sortie ? Diffcile car corrPars inclut fixed and variable corr pars
  res$models <- models
  res$fixef_terms <- processed$fixef_terms ## added 2015/12/09 for predict
  res$fixef_levels <- processed$fixef_levels ## added 2015/12/09 for predict
  res$predictor <- predictor ##  all post fitting functions expect PROCESSED predictor
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
  ## FR->FR $w.ranef and $w.resid not doc'ed, as there is no mention of the augmented mode lin the doc.
  res$w.resid <- w.resid ## useful for calc_info_crits() (+ get_LSmatrix)
  if (models[["eta"]]=="etaHGLM") {
    res$w.ranef <- wranefblob$w.ranef ## useful for calc_info_crits() (+ get_LSmatrix)
    #
    res$lambda.object <- make_lambda_object(nrand, lambda_models=models[["lambda"]], cum_n_u_h, lambda_est, 
                                            init.lambda, ## for attr(.,"type")
                                            coefficients_lambdaS, 
                                            process_resglm_blob, ZAlist, next_LMatrices)
    res$"lambda" <- res$lambda.object$lambda ## redundant but very convenient
    if (attr(ZAlist,"anyRandomSlope")) {
      ## weird coding as next_LMatrices may have NULL elements;
      res$cov.mats <- lapply(next_LMatrices,function(mat) {
        ZWZt(attr(mat,"Lcompact")
             ,exp(coefficients_lambdaS[[which(attr(ZAlist,"Xi_cols")>1L)]]))
      })
    }
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
        phi.object <- c(phi.object, list(glm_phi_args=glm_phi_args, get_glm_phi=create_get_glm_phi() ) )
      } 
      res$phi.object <- phi.object
    }
  } else {
    ## important distinction for (summary, df of LRTs:
    if (is.null(processed$phi.Fix)) { ## absent from original call
      res$phi.object <- list(phi_outer=structure(phi.Fix,type="var")) ## hlcor call of corrHLfit / HLfit call post fitme ?
    } else res$phi.object <- list(phi_outer=structure(phi.Fix,type="fix"))
  }
  ###################
  ## private hack
  ###################
  #    if ( ! is.null(init.HLfitName)) {
  if ( ! is.na(spaMM.getOption("INIT.HLFITNAME"))) {
    nextinit.HLfit <- list()
    nextinit.HLfit$fixef <- beta_eta
    nextinit.HLfit$v_h <- v_h
    if (models[["eta"]]=="etaHGLM" && anyNA(lambda.Fix)) nextinit.HLfit$lambda <- lambda_est
    spaMM.options(INIT.HLFITNAME=nextinit.HLfit)
    ##assign(init.HLfitName, nextinit.HLfit,pos=".GlobalEnv")
  }  
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
    if (models[["eta"]]=="etaHGLM") {
      if (conv_logL  && ! conv.lambda) {
        mainNotConv <- paste("p_v apparently converged but lambda estimates apparently did not.",
                             "\n This may indicate that some lambda estimate(s) should be zero.",
                             "\n Otherwise try increasing 'max.iter' above ",max.iter,
                             "\n (see help(HLfit) for details about 'max.iter')",sep="")          
      } else mainNotConv <- paste("Estimates did not converge. Try increasing 'max.iter' above ",max.iter,
                                  "\n (see help(HLfit) for details about 'max.iter')",sep="")        
      attr(mainNotConv,"diagnostics") <- c( conv.phi=conv.phi , conv.lambda=conv.lambda, 
                                            conv.corr=conv.corr, conv_lambda_vs_u=conv_lambda_vs_u,
                                            conv_rel_lambda=conv_rel_lambda )
    } else {
      mainNotConv <- paste("Estimates did not converge. Try increasing 'max.iter' above ",max.iter,
                           "\n (see help(HLfit) for details about 'max.iter')",sep="")        
      attr(mainNotConv,"diagnostics") <- c( conv.phi=conv.phi , conv.lambda=conv.lambda, 
                                            conv.corr=conv.corr )
    }
    warningList$mainNotConv <- mainNotConv
  }
  res$warnings <- warningList
  res$spaMM.version <- packageVersion("spaMM")
  ###################
  ## Create fns with their own local environment (Chambers, p. 127)
  # The create_... functions have only arguments not in res
  # The created function only has the containing  object are argument 
  ###################
  res$get_w_h_coeffs <- create_get_w_h_coeffs()
  res$get_ZALMatrix <- create_get_ZALMatrix()
  res$get_beta_cov <- create_get_beta_cov()
  res$get_beta_w_cov <- create_get_beta_w_cov()
  res$get_invColdoldList <- create_get_invColdoldList()
  ## Sig may bestored in the envir of $get_logdispObject:
  res$get_logdispObject <- create_get_logdispObject(dvdloglamMat, dvdlogphiMat, 
                                                    muetablob=muetablob, stop.on.error=stop.on.error)
  ###
  res$get_info_crits <- create_get_info_crits(pforpv, p_lambda, processed$p_phi)
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
  # sort(sapply(ls(environment(<object>$get_logdispObject)), function(x)
  # +             object.size(get(x, envir = environment(<object>$get_logdispObject)))),decreasing=TRUE)
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"res")) ## empties the whole local envir except the return value
  return(res)
}