post_process_family <- function(family, ranFix) {
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
  return(ranFix)
}

  resize_lambda <- function(lambda,vec_n_u_h,n_u_h,adjacency_info=NULL) {
  if  (length(lambda)==length(vec_n_u_h)) {
    lambda_est <- rep(lambda,vec_n_u_h)
    if ( ! is.null(adjacency_info)) { 
      lambda_est[adjacency_info$u.range] <- lambda[adjacency_info$whichadj]*adjacency_info$coeffs 
    }
  } else if (length(lambda)==n_u_h) {
    lambda_est <- lambda
  } else {stop("Initial lambda cannot be mapped to levels of the random effect(s).")}
  lambda_est
}

.calc_clik <- function(mu,phi_est,processed) { 
  BinomialDen <- processed$BinomialDen
  loglfn.fix <- processed$loglfn.fix
  y <- processed$y
  family <- processed$family
  theta <- .theta.mu.canonical(mu/BinomialDen,family)  
  if (family$family=="binomial") {
    clik <- sum(loglfn.fix(theta,y/BinomialDen,BinomialDen,1/(phi_est))) ## freq a revoir
  } else {
    phi_est[phi_est<1e-12] <- 1e-10 ## 2014/09/04 local correction, has to be finer than any test for convergence 
    ## creates upper bias on clik but should be more than compensated by the lad
    ## correcting the lad makes an overall upper bias for small (y-theta) at "constant" corrected phi 
    ## this can be compensated by correcting the lad LESS.
    clik <- sum(loglfn.fix(theta,y,eval(processed$prior.weights)/phi_est)) ## note (prior) weights meaningful only for gauss/ Gamma 
    clik
  }
}

.eval_gain_LevM_GLM <- function(LevenbergMstep_result,family, X.pv ,coefold,clikold,phi_est,processed, offset) {  
  dbeta <- LevenbergMstep_result$dbetaV
  beta <- coefold + dbeta
  eta <- drop(X.pv %*% beta) + offset
  if (family$link=="log") {eta[eta>30] <-30} ## cf similar code in muetafn
  if (family$link=="inverse" && family$family=="Gamma") {
    etamin <- sqrt(.Machine$double.eps)
    eta[eta<etamin] <- etamin ## eta must be > 0
  } ## cf similar code in muetafn
  mu <- family$linkinv(eta)
  clik <- .calc_clik(mu=mu, phi_est=phi_est,processed=processed) 
  if (is.infinite(clik) || is.na(clik)) {  
    gainratio <- -1
  } else {
    summand <- dbeta*(LevenbergMstep_result$rhs+ LevenbergMstep_result$dampDpD * dbeta) 
    ## In the summand, all terms should be positive. conv_dbetaV*rhs should be positive. 
    # However, numerical error may lead to <0 or even -Inf
    #  Further, if there are both -Inf and +Inf elements the sum is NaN and the fit fails.
    summand[summand<0] <- 0
    denomGainratio <- sum(summand)
    #cat("eval_gain_LM ");print(c(devold,dev,denomGainratio))
    gainratio <- 2*(clik-clikold)/(1e-8+denomGainratio)
  }
  return(list(gainratio=gainratio,clik=clik,beta=beta,eta=eta,mu=mu))
}  

.calc_etaGLMblob <- function(processed, 
                         mu, eta, muetablob, beta_eta, w.resid, ## those in output
                         phi_est, 
                         off, 
                         maxit.mean, 
                         verbose, 
                         for_intervals=NULL,
                         conv.threshold) {
    BinomialDen <- processed$BinomialDen
    X.pv <- processed$X.pv
    y <- processed$y
    family <- processed$family
    stop.on.error <- processed$stop.on.error
    damping <- 1e-7 ## as suggested by Madsen-Nielsen-Tingleff... # Smyth uses abs(mean(diag(XtWX)))/nvars
    dampingfactor <- 2
    newclik <- .calc_clik(mu=mu,phi_est=phi_est,processed=processed) ## handles the prior.weights from processed
    for (innerj in seq_len(maxit.mean)) {
      ## breaks when conv.threshold is reached
      clik <- newclik
      z1 <- eta+(y-mu)/muetablob$dmudeta-off ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
      ## simple QR solve with LevM fallback
      old_beta_eta <- beta_eta
      tXinvS <- sweep(t(X.pv),2L,w.resid,`*`)  
      rhs <-  tXinvS %*% z1
      if ( ! is.null(for_intervals)) {
        currentDy <- (for_intervals$fitlik-newclik)
        if (currentDy <0) .warn_intervalStep(newclik,for_intervals)
        intervalBlob <- .intervalStep_glm(old_beta=beta_eta,
                                         sXaug=tXinvS%*%X.pv,
                                         szAug=rhs,
                                         for_intervals=for_intervals,
                                         currentlik=newclik,currentDy=currentDy)
        beta_eta <- intervalBlob$beta
      } else {
        qr.XtinvSX <- QRwrap(tXinvS%*%X.pv,useEigen=FALSE) ## Cholwrap tested  ## pas sur que FALSE gagne du temps
        beta_eta <- solveWrap.vector( qr.XtinvSX , rhs ,stop.on.error=stop.on.error)
      }
                                           
      # # PROBLEM is that NaN/Inf test does not catch all divergence cases so we need this :
      eta <- off + drop(X.pv %*% beta_eta) ## updated at each inner iteration
      muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
      newclik <- .calc_clik(mu=muetablob$mu, phi_est=phi_est,processed=processed) 
      if ( is.null(for_intervals) &&
        (newclik < clik-1e-5 || anyNA(beta_eta) || any(is.infinite(beta_eta))) ) { 
        ## more robust LevM
        sqrt.ww <- sqrt(w.resid)
        wX <- calc_wAugX(augX=X.pv,sqrt.ww=sqrt.ww)
        LM_wz <- z1*sqrt.ww - (wX %*% old_beta_eta)
        while(TRUE) { 
          if (inherits(wX,"Matrix")) {
            LevenbergMstep_result <- LevenbergMsolve_Matrix(wAugX=wX,LM_wAugz=LM_wz,damping=damping)
          } else LevenbergMstep_result <- LevenbergMstepCallingCpp(wAugX=wX,LM_wAugz=LM_wz,damping=damping) 
          levMblob <- .eval_gain_LevM_GLM(LevenbergMstep_result=LevenbergMstep_result,
                                         X.pv=X.pv, clikold=clik, family=family,
                                         coefold=old_beta_eta,
                                         phi_est=phi_est, offset=off,
                                         processed=processed)
          if (levMblob$gainratio<0) { ## failure: increase damping and continue iterations
            damping <- dampingfactor*damping
            dampingfactor <- dampingfactor*2
          } else break ## success
        } ## while TRUE
        ## on success :
        damping <- damping * max(1/3,1-(2*levMblob$gainratio-1)^3)  
        dampingfactor <- 2
        dbetaV <- LevenbergMstep_result$dbetaV
        beta_eta <- levMblob$beta 
        eta <- off + drop(X.pv %*% beta_eta) ## updated at each inner iteration
        muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
        newclik <- levMblob$clik
      } else dbetaV <- beta_eta - old_beta_eta
      mu <- muetablob$mu ## needed to update z1
      w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
      if (verbose["trace"]) {
        print(paste("Inner iteration ",innerj,sep=""))
        print_err <- c(beta_eta=beta_eta)
        if (innerj>1) print_err <- c(norm.dbetaV=sqrt(sum(dbetaV^2)),print_err)
        print(print_err)
        print("================================================")
      } 
      if (maxit.mean>1) {
        if (mean(abs(dbetaV)) < conv.threshold) break; ## FR->FR mean(abs) is not standard ?  
      }
      #### done with one inner iteration
    } ## end for (innerj in 1:maxit.mean)
    names(beta_eta) <- colnames(X.pv)
    return(list(eta=eta, muetablob=muetablob, beta_eta=beta_eta, w.resid=w.resid, innerj=innerj))
  }
