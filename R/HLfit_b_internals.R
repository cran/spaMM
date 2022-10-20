.post_process_family_it <- function(family, ranFix, char_mv_it) {
  if (family$family=="COMPoisson") {
    if ( ! is.null(ranFix$COMP_nu) && ! is.na(COMP_nu <- ranFix$COMP_nu[char_mv_it])) { ## optimisation call 
      #  COMP_nu[char_mv_it] to get a NA where [[char_mv_it]] generates an error
      #  But it is important to drop the name that, if present, would mess names of return values of .COMPxxx() fns
      assign("nu",COMP_nu[[1]],envir=environment(family$aic))
      ranFix$COMP_nu[char_mv_it] <- NA
      ranFix$COMP_nu <- na.omit(ranFix$COMP_nu)
    } else {
      checknu <- substitute(nu, env=environment(family$aic)) 
      if (inherits(checknu,"call")) eval(checknu)
    }
  } else if (family$family=="beta_resp") {
    if ( ! is.null(ranFix$beta_prec) && ! is.na(beta_prec <- ranFix$beta_prec[char_mv_it])) { ## optimisation call 
      #  COMP_nu[char_mv_it] to get a NA where [[char_mv_it]] generates an error
      #  But it is important to drop the name that, if present, would mess names of return values of .COMPxxx() fns
      assign("prec",beta_prec[[1]],envir=environment(family$aic))
      ranFix$beta_prec[char_mv_it] <- NA
      ranFix$beta_prec <- na.omit(ranFix$beta_prec)
    } else {
      checkprec <- substitute(prec, env=environment(family$aic)) 
      if (inherits(checkprec,"call")) eval(checkprec)
    }
  } else if (family$family %in% c("negbin1","negbin2")) {
    if ( ! is.null(ranFix$NB_shape) && ! is.na(NB_shape <- ranFix$NB_shape[char_mv_it])) { ## fitme -> HLCor -> HLfit
      assign("shape",NB_shape,envir=environment(family$aic))
      ranFix$NB_shape[char_mv_it] <- NA
      ranFix$NB_shape <- na.omit(ranFix$NB_shape)
    } else if ( ! is.null(ranFix$trNB_shape) && ! is.na(trNB_shape <- ranFix$trNB_shape[char_mv_it])) { ## fitme -> HLCor -> HLfit
      assign("shape",.NB_shapeInv(trNB_shape),envir=environment(family$aic))
      ranFix$trNB_shape[char_mv_it] <- NA
      ranFix$trNB_shape <- na.omit(ranFix$trNB_shape)
    } else if ( ! is.null(ranFix$trbeta_prec) && ! is.na(trbeta_prec <- ranFix$trbeta_prec[char_mv_it])) { ## fitme -> HLCor -> HLfit
      assign("prec",.beta_precInv(trbeta_prec),envir=environment(family$aic))
      ranFix$trbeta_prec[char_mv_it] <- NA
      ranFix$trbeta_prec <- na.omit(ranFix$trbeta_prec)
    } else {
      checktheta <- substitute(shape, env=environment(family$aic)) 
      if (inherits(checktheta,"call")) eval(checktheta)
    }
  }
  return(ranFix)
}

.post_process_respfamilies <- function(family, ranFix, families=NULL) {
  if ( ! is.null(families)) {
    for (mv_it in seq_along(families)) {
      ranFix <- .post_process_family_it(family=families[[mv_it]], ranFix, char_mv_it=as.character(mv_it)) 
    }
    return(ranFix)
  }
  if (family$family=="COMPoisson") {
    if ( ! is.null(ranFix$COMP_nu)) { ## optimisation call
      assign("nu",ranFix$COMP_nu,envir=environment(family$aic))
      ranFix$COMP_nu <- attr(ranFix,"type")$COMP_nu <- NULL
    } else {
      checknu <- substitute(nu, env=environment(family$aic)) 
      if (inherits(checknu,"call")) eval(checknu)
    }
  } else if (family$family=="beta_resp") {
    if ( ! is.null(ranFix$beta_prec)) { ## optimisation call
      assign("prec",ranFix$beta_prec,envir=environment(family$aic))
      ranFix$beta_prec <- attr(ranFix,"type")$beta_prec <- NULL
    } else if ( ! is.null(ranFix$trbeta_prec)) { ## fitme -> HLfit directly (FIXME: unify both cases ?)
      assign("prec",.beta_precInv(ranFix$trbeta_prec),envir=environment(family$aic))
      ranFix$trbeta_prec <- attr(ranFix,"type")$trbeta_prec <- NULL
    } else {
      checkprec <- substitute(prec, env=environment(family$aic)) 
      if (inherits(checkprec,"call")) eval(checkprec)
    }
  } else if (family$family  %in% c("negbin1","negbin2")) {
    if ( ! is.null(ranFix$NB_shape)) { ## fitme -> HLCor -> HLfit
      assign("shape",ranFix$NB_shape,envir=environment(family$aic))
      ranFix$NB_shape <- attr(ranFix,"type")$NB_shape <- NULL
    } else if ( ! is.null(ranFix$trNB_shape)) { ## fitme -> HLfit directly (FIXME: unify both cases ?)
      assign("shape",.NB_shapeInv(ranFix$trNB_shape),envir=environment(family$aic))
      ranFix$trNB_shape <- attr(ranFix,"type")$trNB_shape <- NULL
    } else {
      checktheta <- substitute(shape, env=environment(family$aic)) 
      if (inherits(checktheta,"call")) eval(checktheta)
    }
  }
  return(ranFix)
}

.resize_lambda <- function(lambda,vec_n_u_h,n_u_h,adjacency_info=NULL) {
  if  (length(lambda)<length(vec_n_u_h)) {
    newlambda <- rep(NA,length(vec_n_u_h))
    names(newlambda) <- seq_along(vec_n_u_h)
    newlambda[names(lambda)] <- lambda
    lambda <- newlambda
  } 
  if (length(lambda)==length(vec_n_u_h)) {
    lambda_est <- rep(lambda,vec_n_u_h)
    if ( ! is.null(adjacency_info)) { 
      lambda_est[adjacency_info$u.range] <- lambda[adjacency_info$whichadj]*adjacency_info$coeffs 
    }
  } else if (length(lambda)==n_u_h) {
    lambda_est <- lambda
  } else {stop("Initial lambda cannot be mapped to levels of the random effect(s).")}
  lambda_est
}

.calc_clik <- function(mu, # as given by .muetafn(), ie, prediction of counts in binomial case...
                       phi_est,processed, clik_fn=processed$clik_fn, 
                       prior.weights=processed$prior.weights,
                       summand=FALSE, muetaenv=NULL) { 
  BinomialDen <- processed$BinomialDen
  y <- processed$y
  if ( ! is.null(vec_nobs <- processed$vec_nobs)) { # mv case, list of families
    cum_nobs <- attr(processed$families,"cum_nobs")
    cliks <- vector("list",length(vec_nobs)) # clik_fn returns a single value when it uses family()$aic  (e.g., poisson) (so summand=FALSE fails) and a vector otherwise
    for (mv_it in seq_along(vec_nobs)) {
      fam <- processed$families[[mv_it]]
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      theta <- .theta.mu.canonical(mu[resp_range]/BinomialDen[resp_range],fam)
      if (fam$family=="binomial") {
        cliks[[mv_it]] <- clik_fn[[mv_it]](theta, y[resp_range]/BinomialDen[resp_range], BinomialDen[resp_range], eval(prior.weights[[mv_it]])/phi_est[[mv_it]])
      } else {
        phi <- phi_est[[mv_it]]
        phi[phi<1e-12] <- 1e-10
        cliks[[mv_it]] <- clik_fn[[mv_it]](theta, y[resp_range], eval(prior.weights[[mv_it]])/phi)
      }
    }
    cliks <- unlist(cliks, recursive = FALSE, use.names = FALSE)
  } else {
    # clik_fn has theta as argument (and may indeed may sometimes use the theta value rather than mu) hence uses the canonical link to deduce mu from theta (or gets the mu attribute of theta) 
    # If clik_fn was written in terms of mu the .theta.mu.canonical() call could be avoided. 
    # Indeed for 'LLM's theta=mu (identity pseudo-canonical lonk).
    # Only other use of clik_fn is for cAIC bootstrap correction, only for GLM families...
    family <- processed$family
    theta <- .theta.mu.canonical(mu/BinomialDen,family)  # Could be a bottleneck for CMP if attr(mu,"lambda") were missing.
    if (family$family=="binomial") {
      cliks <- clik_fn(theta, y/BinomialDen, BinomialDen, eval(prior.weights)/phi_est)
    } else if (family$family=="COMPoisson") {
      cliks <- clik_fn(theta, y, eval(prior.weights)/phi_est, muetaenv=muetaenv)
    } else if (family$family=="beta_resp") {
      cliks <- clik_fn(theta, y, pw=eval(prior.weights))
    } else {
      phi_est[phi_est<1e-12] <- 1e-10 ## 2014/09/04 local correction, has to be finer than any test for convergence 
      ## creates upper bias on clik but should be more than compensated by the lad
      ## correcting the lad makes an overall upper bias for small (y-theta) at "constant" corrected phi 
      ## this can be compensated by correcting the lad LESS.
      cliks <- clik_fn(theta, y, eval(prior.weights)/phi_est)
    }
  }
  if (summand) {
    attr(cliks,"unique") <- NULL
    return(cliks)
  } else return(sum(cliks))
}

# this is used for LLM too
.eval_gain_clik_LevM <- function(LevenbergMstep_result,family, X.pv ,coefold,clikold, phi_est, processed, offset) {  
  dbeta <- LevenbergMstep_result$dbetaV
  beta <- coefold + dbeta
  eta <- drop(X.pv %*% beta) + offset
  eta <- .sanitize_eta(eta, y=processed$y, family=family, max=40) 
  ## Here I can use 
  muetablob <- .muetafn(eta=eta, BinomialDen=processed$BinomialDen, processed=processed, phi_est=phi_est) 
  ## which returns a $mu=muCOUNT in all cases.   
  # if (family$family=="binomial") {
  #   muFREQS <- family$linkinv(eta)
  #   mu <- muFREQS * processed$BinomialDen
  # } else mu <- family$linkinv(eta)
  clik <- .calc_clik(mu=muetablob$mu, phi_est=phi_est,processed=processed) 
  if (is.infinite(clik) || is.na(clik)) {  
    gainratio <- -1
    conv_clik <- Inf
  } else {
    summand <- dbeta*(LevenbergMstep_result$rhs+ LevenbergMstep_result$dampDpD * dbeta) 
    ## In the summand, all terms should be positive. conv_dbetaV*rhs should be positive. 
    # However, numerical error may lead to <0 or even -Inf
    #  Further, if there are both -Inf and +Inf elements the sum is NaN and the fit fails.
    summand[summand<0] <- 0
    denomGainratio <- sum(summand)
    #cat("eval_gain_LM ");print(c(devold,dev,denomGainratio))
    dclik <- clik-clikold
    conv_clik <- abs(dclik)/(1+abs(clik))
    gainratio <- 2*dclik/denomGainratio 
  }
  return(list(gainratio=gainratio,clik=clik,beta=beta,eta=eta,
              muetablob=muetablob, # with $mu=muCOUNT
              conv_clik=conv_clik))
}  

# fixed-effect main response
# HLfit_body -> .calc_etaGLMblob -> this fn
.do_damped_WLS_glm <- function(wX, LM_wz, damping, X.pv, clik, family, old_beta_eta, phi_est, off, processed, verbose) {
  restarted <- FALSE
  dampingfactor <- 2
  while(TRUE) { 
    if (inherits(wX,"Matrix")) {
      LevenbergMstep_result <- .LevenbergMsolve_Matrix(wAugX=wX,LM_wAugz=LM_wz,damping=damping) 
    } else LevenbergMstep_result <- .LevenbergMstepCallingCpp(wAugX=wX,LM_wAugz=LM_wz,damping=damping) 
    levMblob <- .eval_gain_clik_LevM(LevenbergMstep_result=LevenbergMstep_result,
                                       X.pv=X.pv, clikold=clik, family=family,
                                       coefold=old_beta_eta,
                                       phi_est=phi_est, offset=off,
                                       processed=processed)
    gainratio <- levMblob$gainratio
    conv_crit <- max(levMblob$conv_clik, 
                     abs(LevenbergMstep_result$rhs)/(1+clik))
    if (is.nan(gainratio)) {
      break
    } else if (gainratio>0) { ## success
      damping <- damping * max(1/3,1-(2*gainratio-1)^3)  
      dampingfactor <- 2
      break
    } else if (dampingfactor>4 ## iter not accessible for test
               && gainratio==0) { # apparently flat deviance
      if (conv_crit < 1e-8) { # we were at optimum
        damping <- 1e-7
        if (verbose["trace"]) cat("#")
        break ## apparently flat dev
      } else { ## measurable gradient, but we don't move => too high damping ? (if damping from previous LevM step was too high)
        if ( ! restarted) { # condition to avoid infinite loop
          damping <- 1e-7
          dampingfactor <- 2
          restarted <- TRUE
          if (verbose["trace"]) cat("-")
          # and continue # but this allows an  infinite loop
        } else {
          # hmm; well, break and diagnose...
          if (verbose["trace"]) cat("?")
          break
        }
      }
    } else { ## failure: increase damping and continue iterations
      damping <- dampingfactor*damping
      dampingfactor <- dampingfactor*2
    } 
    if (damping>1e10) break # stop("reached damping=1e10")
  } ## while TRUE
  RESU <- levMblob
  RESU$damping <- damping
  RESU$dbetaV <- LevenbergMstep_result$dbetaV
  return(RESU)
}


# called by HLfit_body for models[[1]]=="etaGLM"; tested eg by test-COMPoisson:
.calc_etaGLMblob <- function(processed, muetablob, 
                         mu=muetablob$mu, eta=muetablob$sane_eta, 
                         old_beta_eta, ## scaled since X.pv is scaled; same for return beta_eta. An init.HLfit$fixef would be (.)/attr(spaMM:::.scale(zut$X.pv),"scaled:scale")
                         w.resid,
                         phi_est, 
                         off=processed$off, 
                         maxit.mean, 
                         verbose, 
                         for_intervals=NULL,
                         Xtol_rel=processed$spaMM_tol$Xtol_rel) {
    BinomialDen <- processed$BinomialDen
    X.pv <- processed$AUGI0_ZX$X.pv
    y <- processed$y
    family <- processed$family
    LM_called <- FALSE
    qr_X <- NA
    damping <- 1e-7 ## as suggested by Madsen-Nielsen-Tingleff... # Smyth uses abs(mean(diag(XtWX)))/nvars
    newclik <- .calc_clik(mu=mu,phi_est=phi_est,processed=processed) ## handles the prior.weights from processed
    for (innerj in seq_len(maxit.mean)) {
      ## breaks when Xtol_rel is reached
      clik <- newclik
      # Historical oddity: the fit has worked with code which was OK for solving, but not for CI as the CI code suppresses 
      # a column of the design matrix, which is not sufficient on the premultiplied (scaled X) system.
      z1 <- .calc_z1(muetablob, w.resid, y, off, cum_nobs=attr(processed$families,"cum_nobs"))
      if (is.list(w.resid)) {
        sqrtW <- sqrt(w.resid$w_resid)        
      } else sqrtW <- sqrt(w.resid)
      wX <- .Dvec_times_m_Matrix(sqrtW, X.pv) ## keeps colnames: important for intervalStep_glm
      szAug <- sqrtW * z1  
      # names(szAug) <- colnames(X.pv) ## also important for intervalStep_glm
      ## simple QR solve with LevM fallback
      if ( ! is.null(for_intervals) || ! LM_called) {
        if ( ! is.null(for_intervals)) {
          currentDy <- (for_intervals$fixeflik-newclik)
          if (currentDy < -1e-4 && 
              (is.null(bestlik <- processed$envir$confint_best$lik) || newclik > bestlik)) {
            if (is.null(bestlik)) {
              locmess <- paste("A higher",names(for_intervals$fixeflik),"was found than for the original fit.",
                               "\nThis suggests the original fit did not fully maximize",names(for_intervals$fixeflik),
                               "\nExpect more information at end of computation.")
              message(locmess)
            }
            processed$envir$confint_best$lik <- newclik
            processed$envir$confint_best$beta_eta <- .unscale(X.pv, old_beta_eta)
          }
          intervalBlob <- .intervalStep_glm(old_beta=old_beta_eta,
                                            sXaug=wX,
                                            szAug=szAug,
                                            for_intervals=for_intervals,
                                            currentlik=newclik,currentDy=currentDy)
          beta_eta <- intervalBlob$beta
        } else if ( ! LM_called)  {
          qr_X <- qr(wX,tol=spaMM.getOption("qrTolerance")) 
          beta_eta <- .safesolve_qr_vector(qr_X, szAug)
          beta_eta <- drop(beta_eta)
        }
        names(beta_eta) <- colnames(X.pv)
        # # PROBLEM is that NaN/Inf test does not catch all divergence cases so we need this :
        eta <- off + drop(X.pv %*% beta_eta) ## updated at each inner iteration
        muetablob <- .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed, phi_est=phi_est) 
        newclik <- .calc_clik(mu=muetablob$mu, phi_est=phi_est,processed=processed) 
      }  
      if ( is.null(for_intervals) &&
        (LM_called || # always use LevM when it has been called before: motivated by spaMM_glm example, cf email Alex, 01/06/2022, 14:13
         newclik < clik-1e-5 || anyNA(beta_eta) || any(is.infinite(beta_eta))) ) { 
        ## more robust LevM
        LM_called <- TRUE
        LM_wz <- z1*sqrtW - (wX %*% old_beta_eta)
        damped_WLS_blob <- .do_damped_WLS_glm(wX, LM_wz, damping, X.pv, clik, family, old_beta_eta, phi_est, off, processed, verbose)
        beta_eta <- damped_WLS_blob$beta 
        eta <- damped_WLS_blob$eta #off + drop(X.pv %*% beta_eta) ## updated at each inner iteration
        muetablob <- damped_WLS_blob$muetablob # .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
        newclik <- damped_WLS_blob$clik
        damping <- damped_WLS_blob$damping
        dbetaV <- damped_WLS_blob$dbetaV
      } else dbetaV <- beta_eta - old_beta_eta
      mu <- muetablob$mu ## needed to update z1
      w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est, obsInfo=processed$how$obsInfo) ## 'weinu', must be O(n) in all cases
      if (verbose["trace"]) {
        print(paste0("Inner iteration ",innerj))
        print_err <- c(beta_eta=beta_eta)
        if (innerj>1L) print_err <- c(norm.dbetaV=sqrt(sum(dbetaV^2)),print_err)
        print(print_err)
        print("================================================")
      } 
      #if(innerj>600) browser()
      if (maxit.mean>1 && (damping>1e100 || mean(abs(dbetaV)) < Xtol_rel)) break
      #### done with one inner iteration
      old_beta_eta <- beta_eta
    } ## end for (innerj in 1:maxit.mean)
    names(beta_eta) <- colnames(X.pv)
    return(list(eta=muetablob$sane_eta, muetablob=muetablob, beta_eta=beta_eta, w.resid=w.resid, innerj=innerj,
                sXaug=structure(NA,class="(G)LM"), qr_X=qr_X))
  }

.calc_std_leverages <- function(models, need_ranefPars_estim, phi.Fix, auglinmodblob, n_u_h, nobs, processed, w.resid, u_h, 
                                need_simple_lambda, muetablob, wranefblob, ZAL, lambda_est, cum_n_u_h, 
                                lcrandfamfam=attr(processed$rand.families,"lcrandfamfam"), phi_est) {
  hatvals <- NULL 
  if (models[[1]]=="etaHGLM") {
    if (need_ranefPars_estim || .anyNULL(phi.Fix)) {
      hatval <- .get_hatvalues_MM(auglinmodblob$sXaug,X.Re=processed$X.Re, auglinmodblob$weight_X)
      hatvals <- list(
        ranef=hatval[seq(n_u_h)],  ## for the ranef residuals (lambda)
        resid=hatval[(n_u_h+1L):(n_u_h+nobs)] ## for the error residuals (phi)
      )
    }
  } else if (.anyNULL(phi.Fix)) { ## phi estim in GLM fitted by ML. 
    hatvals <- list(resid=.get_hatvalues_FM(processed$X.Re, augX=processed$AUGI0_ZX$X.pv, w.resid))
  }
  ## (HL[2]=0, HL[3]=0): previous hat matrix -> p 
  ## (HL[2]=0, HL[3]=1): notEQL -> tilde(p), (HL[2]=1 && ): full correction -> q 
  ## (HL[2]=1, HL[3]=1): full correction -> q 
  leverages <- .hatvals2std_lev(
    hatvals, #### base from hat matrix
    sXaug=auglinmodblob$sXaug, anynull_phi.Fix=.anyNULL(phi.Fix), u_h=u_h, 
    need_simple_lambda=need_simple_lambda, muetablob=muetablob, wranefblob=wranefblob, nobs=nobs, ZAL=ZAL, 
    w.resid=w.resid, 
    psi_M=processed$psi_M, lambda_est=lambda_est, cum_n_u_h=cum_n_u_h, lcrandfamfam=lcrandfamfam, 
    nrand=length(lcrandfamfam), phi_est=phi_est,
    processed=processed)
  leverages
}

.calc_APHLs_GLM <- function(processed, w.resid, clik) {
  ## ML: X.Re non NULL mais ncol(X.Re)=0
  X.REML <- processed$X.Re
  if (is.null(X.REML)) {X.REML <- processed$AUGI0_ZX$X.pv} ## REML standard
  if ( ncol(X.REML) ) { ## REML standard || REML non standard
    Md2hdb2 <- .ZtWZwrapper(X.REML,w.resid)
    ladb <- .LogAbsDetWrap(Md2hdb2,logfac=-log(2*pi))
    p_bv <- clik - ladb/2
    if ( ! is.null(X_scale <- attr(processed$AUGI0_ZX$X.pv,"scaled:scale"))) p_bv <- p_bv -sum(log(X_scale))
    APHLs <- list(p_v=clik, p_bv=p_bv)
  } else APHLs <- list(p_v=clik, p_bv=clik)
  APHLs
}

.get_fixed_adjacency_info <- function(whichadj, LMatrices, cum_n_u_h, corr_est, ranFix, init.HLfit) {
  test_inner_estim_rho <- (length(whichadj) && ! is.null(adjd <- attr(LMatrices[[whichadj]],"symsvd")$adjd))
  if (test_inner_estim_rho) {  
    u.range <- (cum_n_u_h[whichadj]+1L):cum_n_u_h[whichadj+1L]
    adj_rho <- corr_est$rho
    if (is.null(adj_rho)) adj_rho <- .getPar(ranFix,"rho") ## could this occur with test_inner_estim_rho ?
    if (is.null(adj_rho)) adj_rho <- init.HLfit$rho
    if ( ! is.null(adj_rho)) fixed_adjacency_info <- list(whichadj=whichadj, u.range=u.range, coeffs=1/(1-adj_rho*adjd))
  } else fixed_adjacency_info <- NULL
  fixed_adjacency_info
}

.get_init_beta <- function(processed, pforpv, init.HLfit, X.pv=processed$AUGI0_ZX$X.pv) {
  beta_eta <- processed$port_env$port_fit_values$fixef # scaled
  if (is.null(beta_eta)) { # then we look at user value and scale it <=> user values must be unscaled.
    if (pforpv) {
      beta_eta <- init.HLfit$fixef
      if ( ! is.null(beta_eta) && ! is.null(attr(X.pv,"scaled:scale"))) {
        beta_eta <- .scale(beta=beta_eta,X=X.pv)
      }
    } else beta_eta <- numeric(0L) # don't leave it NULL so that we don't try to get it from inits_by_glm.
  } 
  beta_eta
}

.p_corr_inner_rC <- function(ranCoefs_blob, LMatrices) {
  var_ranCoefs <- with(ranCoefs_blob, (isRandomSlope & ! is_set)) # inner ranCoefs
  if (any(var_ranCoefs)) {
    p_corr <- vector("list", length(LMatrices))
    for (it in which(var_ranCoefs)) {
      dimL <- NROW(attr(LMatrices[[it]],"latentL_blob")$compactcovmat)
      if (dimL==0L) stop("'latentL_blob' attribute missing to LMatrices[[it]].")
      p_corr[[it]] <- (dimL-1)*dimL/2
    }
    sum(.unlist(p_corr))
  } else 0L
}
