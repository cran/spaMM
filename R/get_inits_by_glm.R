.safe_llm <- function(X.pv, Yy, prior.weights, off, family, processed, resu) {
  resllm <- llm.fit(x=X.pv,  y=Yy,  weights = eval(prior.weights),  offset = off, 
                    family = family,  control = processed[["control.glm"]])
  if ( (! resllm$converged) && family$link !="inverse") {
    .spaMM.data$options$xLM_conv_crit <- list(max=-Inf) # So that no-convergence in later .calc_dispGammaGLM) call can be independently assessed
    ## divergence more likely than from spaMM_glm.fit, and fewer options then...  a simple control of 'start' may work well for log link:
    altresllm <- llm.fit(x=X.pv,  y=Yy,  weights = eval(prior.weights),  offset = off, 
                         start = rep(0,ncol(X.pv)),
                         family = family,  control = processed[["control.glm"]])
    if (altresllm$deviance < resllm$deviance) resllm <- altresllm
  }
  if ( (! resllm$converged ) && (conv_crit <- .spaMM.data$options$xLM_conv_crit$max>0)) { 
    if ( ! .spaMM.data$options$xLM_conv_silent) resu$conv_info <- 
        paste("In initialization, llm.fit() did not converge within",
              resllm$iter,"iterations (criterion:",
              paste(names(conv_crit),"=",signif(unlist(conv_crit),3),collapse=", "),").\n",
              "Use <fitting function>(., control.glm=list(maxit=<larger value>,..)) to maybe fix this.")  # but in practice this has occurred bc there were bugs in the family member fns. 
    .spaMM.data$options$xLM_conv_crit <- list(max=-Inf) # So that no-convergence in later .calc_dispGammaGLM) call can be independently assessed
  }
  resllm
}

# (1) This handles truncated fams (since glm -> glm.fit handles Tpoisson() etc)
# (2) This ignores obsInfo so results will differ from the equivalent LLF family
# (3) note possible overdisp re-assignment.
.safe_glm <- function(locfamfam, resu, X.pv, Y, prior.weights, off, family, processed) {
  if (locfamfam=="COMPoisson") {
    glm.fit <- glm.nodev.fit
    # Redef overdisp to account for the fact that deviance is not stored in resu$resxlm:
    resxlm <- NULL # only bc R CMD check is unable to diagnose the delayedAssign properly.
    delayedAssign("overdisp", {
      muetaenv <- resxlm$muetaenv
      unsc_dev <- sum(resxlm$family$dev.resids(y=resxlm$y, mu=resxlm$fitted.values, muetaenv=muetaenv)) # muetaenv is in 'resu' environment 
      # => prior weights ignored in this calculation bc they are not scaling factor for logl. 
      as.numeric(unsc_dev/df.residual(resxlm))
    }, eval.env = resu, assign.env = resu)
  } 
  
  tryglm <- .tryCatch_W_E(glm.fit(x=X.pv, 
                                  y=Y, 
                                  weights = eval(prior.weights), 
                                  offset = off, family = family, 
                                  control = processed[["control.glm"]]))
  if (inherits((resglm <- tryglm$value),"error") || 
      ( ! resglm$converged ### && any(fitted(resglm)>1e20)
      ) # this occurred in Gamma(log) models or negbin(log)
  ) { # not often... hzut in test-COMPoisson-difficult (glm not converged).
    # use more robust method, and for log link, additionally set an initial value that will override the default mustart <- y of family()$initialize..
    if (family$link=="log") {
      start <- rep(0, ncol(X.pv))
    } else start <- NULL
    resglm <- spaMM_glm.fit(x=X.pv, # -> LevM -> COMP bottleneck
                            y=Y, 
                            weights = eval(prior.weights), 
                            offset = off, family = family, start=start,
                            control = processed[["control.glm"]])
    if ( (! resglm$converged ) && (conv_crit <- .spaMM.data$options$xLM_conv_crit$max>0)) { 
      if ( ! .spaMM.data$options$xLM_conv_silent) resu$conv_info <- 
          paste("In initialization, spaMM_glm.fit() did not converge within",
                resglm$iter,"iterations (criterion:",
                paste(names(conv_crit),"=",signif(unlist(conv_crit),3),collapse=", "),").\n",
                "Use <fitting function>(., control.glm=list(maxit=<larger value>,..)) to control this.")
      # This info will then be copied as part of the whole 'inits_by_xLM' into processed$envir$inits_by_xLM;
      #  a later .calc_inits_by_xLM() call may overwrite the whole 'inits_by_xLM' so that only the last initialization's $conv_info can be reported,
      #  by a message in a fit summary.
      # .calc_dispGammaGLM() has a distinct tracking mechanism: a *warning* may be emitted each time .check_conv_dispGammaGLM_reinit() is called,
      #  that is, within each call of the main fitting functions, incl. repetitive 'inner' calls of (say) HLfit() within fitme(). 
      #  Nevertheless, the .calc_dispGammaGLM() tracking mechanism also reads convergence info in .spaMM.data$options$xLM_conv_crit 
      #  since this is where xLM fitting functions write it. Hence the need to reset it:
      .spaMM.data$options$xLM_conv_crit <- list(max=-Inf) # So that no-convergence in later .calc_dispGammaGLM) call can be independently assessed
    }
  } else if ( ! is.null(tryglm$warning) && 
              substring(tryglm$warning$message,0,14)=="maxn truncated") .spaMM.data$options$COMP_maxn_warned <- FALSE
  ## -> see explanations in .COMP_maxn()
  resglm
}



.calc_inits_by_xLM <- function(processed, X.pv=processed$AUGI0_ZX$X.pv, family=processed$family,
                               y=processed$y, ## requested by the formula
                               Y=.get_from_terms_info(terms_info=processed$main_terms_info, which="Y"),
                               prior.weights=processed$prior.weights,   ## do not try to eval() it outside of the .wfit function call; else nasty crashes may occur.
                               off=processed$off
) {
  overdisp <- NULL # only for R CMD check
  ## Prior to optimization the family parameters may not be assigned, so: 
  orifamfam <- family$family
  if (orifamfam %in% c("negbin1","negbin2")) {
    if (variable_rdispar <- inherits(substitute(shape, env=environment(family$aic)),"call")) family <- 
        Poisson(family$link, trunc=environment(family$aic)$trunc) # equiv to negbin with minimal dispersion
  } else if (orifamfam=="beta_resp") {
    if (variable_rdispar <- inherits(substitute(prec, env=environment(family$aic)),"call")) { # only if beta_prec not assigned
      loc_link <- family$link
      family <- beta_resp(link=loc_link, prec=1e4) 
    }
  } else if (orifamfam=="betabin") {
    if (variable_rdispar <- inherits(substitute(prec, env=environment(family$aic)),"call")) { # only if beta_prec not assigned
      loc_link <- family$link
      family <- betabin(link=loc_link, prec=1e4) # well that's ~binomial except that deviance() returns an number for it while it returns NULL for a binomial fit 
    } # else use original family with fixed prec. The test case is bbin_llmm_fix. 
  } else if (orifamfam=="COMPoisson") {
    loc_link <- family$link
    if (loc_link=="loglambda") loc_link <- "log"
    if (variable_rdispar <- inherits(substitute(nu, env=environment(family$aic)),"call")) family <- poisson(loc_link) # only if COMP_nu not assigned
    # old _F I X M E__ : much time is spent on inits for glm when nu<1
    # if (inherits(nuchk <- substitute(nu, env=environment(family$aic)),"call") || nuchk >0.25) family <- poisson(loc_link)
    # : This helps for he GLM but not much otherwise. Other ideas ?
  }
  # Given that overdisp is relative to 'prec'=1e4 above, 
  # some of the code below is designed to return e4 when overdisp=0.
  # However, that makes less sense with the additional bracketing.
  
  locfamfam <- family$family # family may have been modified locally!
  if (locfamfam  %in% c("binomial","betabin") && NCOL(y)==1L) { 
    BinomialDen <- processed$BinomialDen ## requested by the formula
    begform <-"cbind(y,BinomialDen-y)~"  
  } else {begform <-"y~"}
  ###################################################if (pforpv==0) {endform <-"0"} else 
  pforpv <- ncol(X.pv)
  if( ! pforpv && # no predictor variable
      (locfamfam %in% c("binomial","poisson","beta_resp", "betabin") || 
       # 1st condition => includes cases where original family was (T)negbin with free disp param (and pforpv=0)L, a *poisson* GLM *with* an Intercept is fitted (the resulting beta is ignored);
       # 2nd condition => if original family was negbin with FIXED disp param, a negbin GLM is fitted; it must have an intercept otherwise 
       #                                                  eta=0=> mu untruncated=1 => y>1 is impossible (and returns negative deviance=> negative init lambda)
       (locfamfam  %in% c("negbin1","negbin2") && # with FIXED shape
        family$zero_truncated) 
      )) X.pv <- matrix(1,ncol=1, nrow=nrow(X.pv))
  n_lambda <- sum(attr(processed$ZAlist,"Xi_cols"))
  if (locfamfam %in% c("Gamma","gaussian")) {lam_fac <- 1/(n_lambda+1L)} else lam_fac <- 1/n_lambda
  
  resu <- new.env(parent=environment(.calc_inits_by_xLM)) 
  delayedAssign("overdisp", {
    # When resxlm is not a glm nor an HLfit object, this uses the stats:::<>.default method of the two extractors.
    as.numeric(deviance(resxlm)/df.residual(resxlm))
  }, eval.env = resu, assign.env = resu)
  # Using "prior-weighted phi-unscaled" deviance residuals in the GLM-family case.
  # which are extracted by [deviance(), taking prior weights into account for families where it is a scaling factor for logl]
  # or recomputed by with(resxlm, family$dev.resids(y, mu=fitted.values, weights)) #    note 'weights' 
  ## => For non-GLM families, prior weights have a distinct meaning and are handled differently.
  # In any case, this is distinct from resid(resxlm) which gives resxlm$residuals, which are the "working residuals" ie (y-mu)/mu.eta, 
  
  if (locfamfam=="gaussian" && family$link=="identity") {
    if (inherits(X.pv,"sparseMatrix")) {
      resxlm <- .spaMM_lm.wfit(x=X.pv,y=y,offset=off,w=eval(prior.weights))
    } else resxlm <- lm.wfit(x=X.pv,y=y,offset=off,w=eval(prior.weights))
  } else if ( locfamfam == "betabin") { 
    resxlm <- .safe_llm(X.pv, Yy=Y, prior.weights, off, family, processed, resu)
  } else if ( locfamfam %in% c("negbin1", "beta_resp")) { # with variable dispersion param; I cannot use $flags bc e.g. for COMPoisson, the local family may be stats::poisson
    resxlm <- .safe_llm(X.pv, Yy=drop(Y), prior.weights, off, family, processed, resu)
  } else { ## family is GLM 
    resxlm <- .safe_glm(locfamfam, resu, X.pv, Y, prior.weights, off, family, processed)
  }
  # 
  # I cannot easily use family $flags to distinguish cases bc e.g. for COMPoisson, the local family may be stats::poisson
  
  if (locfamfam=="gaussian" && family$link=="identity") {
    fv <- fitted(resxlm)
    dev <- resid(resxlm)^2
    if (! is.null(resxlm$weights)) dev <- dev/resxlm$weights
    guess <- sum(dev)/resxlm$df.residual
    if (is.nan(guess)) { # resxlm$df.residual=0, possibly wider issue with requested fit, cannot be resolved from here.
      resu$lambda <- resu$phi_est <- 1e-04
    } else resu$lambda <- resu$phi_est <- sum(dev)/resxlm$df.residual # /lam_fac
  } else if (orifamfam == "betabin") {
    # overdisp <- deviance(resxlm)/df.residual(resxlm)
    delayedAssign("lambda", {
      0.025*log(1+overdisp) # 
    }, eval.env = resu, assign.env = resu)
    if (variable_rdispar) delayedAssign("beta_prec", {
      min(100, max(0.15, 100/(1e-2+overdisp))) # bbin_llmm has init=45.33... and constraining <10 is a bad idea...
    }, eval.env = resu, assign.env = resu)
  } else if (orifamfam=="beta_resp") { # Now, 'overdisp is 1 if lambda=0' We need lambda low in that case, but not too low (optim issues)
    delayedAssign("lambda", {
      #max(0.1,min(log(overdisp),1)) #  1st beta_resp example of the doc is good test.
      0.025*log(1+overdisp) # 
    }, eval.env = resu, assign.env = resu)
    if (variable_rdispar) delayedAssign("beta_prec", {
      min(10, max(0.15, 100/(1e-2+overdisp))) 
    }, eval.env = resu, assign.env = resu)
  } else if ( orifamfam =="negbin1") {
    delayedAssign("lambda", {
      #max(0.1,log(overdisp)) 
      0.0125*log(1+overdisp) # 
    }, eval.env = resu, assign.env = resu)
    if (variable_rdispar) delayedAssign("NB_shape", {
      min(10, max(0.15, 2/(2e-4+overdisp))) # test: negbin1 in test-LLM, notably tnb1
    }, eval.env = resu, assign.env = resu)
    ## I tried weird things in version 4.0.14. In any case don't forget there is also a distinct .calc_fam_corrected_guess() function.
  } else { # Other cases, where locfamfam is a GLM family, but orifamfam was possibly non-GLM
    # lambda (before corrections beyond this function, such as ceiling the initial value by 0.2 etc.)
    if (orifamfam=="binomial") {
      if (max(resxlm$prior.weights)==1L) { ## binary response
        resu$lambda <- 1
        # alternatives below technically works but does not seem appropriate.
      } else {
        delayedAssign("lambda", {
          # 0.025*log(1+overdisp) # unconvincing attempt 03/2023
          overdisp
        }, eval.env = resu, assign.env = resu)
      }
    } else {
      delayedAssign("lambda", {
        overdisp
      }, eval.env = resu, assign.env = resu)
    } 

    # fam parm and phi
    # if (orifamfam == "negbin1") {
    #   delayedAssign("NB_shape", {
    #     min(10, max(0.15, 2/overdisp)) # test: negbin1 in test-LLM, notably tnb1
    #   }, eval.env = resu, assign.env = resu)
    # } else 
    if (orifamfam == "negbin2") {
      delayedAssign("NB_shape", {
        # __F I X M E___? try  something like       min(10, max(0.15, 2/(2e-4+overdisp))) ?
        min(20, max(0.15, 4/overdisp)) # test: spaMMintro's HLnegbin case; fit_05 of twinR back-compat check
      }, eval.env = resu, assign.env = resu)
    } else  delayedAssign("phi_est", {
      overdisp
    }, eval.env = resu, assign.env = resu)
    
  }
  
  if (pforpv) {
    ## Two potential problems (1) NA's pour param non estimables (cas normal); 
    ## (2) "glm.fit: fitted probabilities numerically 0 or 1 occurred" which implies separation or large offset
    if (max(abs(c(coefficients(resxlm))),na.rm=TRUE)>1e10) { ## na.rm v1.2 
      warning(paste0("(!) Apparent divergence of estimates in a fixed-effect pre-fit.\n",
                     "    Check your data for extreme values, separation or bad offset values."))
    } 
    beta_eta <- c(coefficients(resxlm)) ## this may include NA's. Testcase: HLfit(Strength ~ Material*Preheating+Method,data=weld)
    #    "c is sometimes used for its side effect of removing attributes except names..."
    if (all(names(beta_eta)=="X.pv")) { ## si la formula etait y ~X.pv-1
      names(beta_eta) <- colnames(resxlm$model$X.pv)
    } else names(beta_eta) <- gsub("^X\\.pv","",names(beta_eta)) ## removes initial "X.pv" without guessing any order or length
    resu$beta_eta <- beta_eta
  } else resu$beta_eta <- numeric(0) # for combining it in initialization of mv fit
  resu$resxlm <- resxlm
  resu$converged <- resxlm$converged
  return(resu)
}

# function called within HLfit for missing inits... 
# but also sometimes in fitme_body prior to optimization (cf complex condition for call of .eval_init_lambda_guess())
.get_inits_by_xLM <- function(processed, reset=FALSE,
                              # args used only as promises:
                              X.pv=processed$AUGI0_ZX$X.pv, family=processed$family) {
  if ( is.null(inits_by_xLM <- processed$envir$inits_by_xLM)) {
    if ( ! is.null(vec_nobs <- processed$vec_nobs) ) {
      do_eval <- rep(TRUE, length(vec_nobs))
      new_mvlist <- new_phi_ests <- vector("list", length(vec_nobs))
      beta_eta <- numeric(ncol(X.pv))
      lambdas <- numeric(length(vec_nobs))
    } else do_eval <- TRUE
  } else if ( ! is.null(vec_nobs <- processed$vec_nobs) ) {
    do_eval <- logical(length(vec_nobs))
    for (mv_it in seq_along(vec_nobs)) {
      family <- processed$families[[mv_it]] # necessary bc reset=quote(...) uses it!
      do_eval[mv_it] <- eval(reset)
    }
    new_mvlist <- inits_by_xLM$mvlist
    new_phi_ests <- inits_by_xLM$phi_ests 
    beta_eta <- inits_by_xLM$beta_eta
    lambdas <- inits_by_xLM$lambdas
  } else {
    family <- processed$family
    do_eval <- eval(reset)
  }
  #
  if (any(do_eval)) {
    if ( ! is.null(vec_nobs)) {
      cum_nobs <- attr(processed$families,"cum_nobs")
      col_ranges <-  attr(X.pv,"col_ranges")
      for (mv_it in which(do_eval)) {
        resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
        col_range <-  col_ranges[[mv_it]]
        new_mvlist[[mv_it]] <- .calc_inits_by_xLM(processed, X.pv=X.pv[resp_range,col_range, drop=FALSE],                                                   
                                                  family=processed$families[[mv_it]],
                                                  y=processed$y[resp_range], ## requested by the formula
                                                  Y=.get_from_terms_info(terms_info=processed$main_terms_info, which="Y", mv_it=mv_it), 
                                                  prior.weights=processed$prior.weights[[mv_it]], ## do not try to eval() it outside of the .wfit function call; else nasty crashes may occur.
                                                  off=processed$off[resp_range] )
        # in mv cases the benefits of promises in .calc_inits_by_xLM() return value are lost (__F I X M E___):
        beta_eta[col_range] <- new_mvlist[[mv_it]]$beta_eta
        lambdas[mv_it] <- new_mvlist[[mv_it]]$lambda
        if (processed$families[[mv_it]]$family %in% c("gaussian","Gamma")) {
          new_phi_ests[[mv_it]] <- new_mvlist[[mv_it]]$phi_est
        } else new_phi_ests[[mv_it]] <- 1 # no impact anyway?
      }
      names(new_mvlist) <- seq_along(do_eval)
      names(new_phi_ests) <- which(do_eval)
      processed$envir$inits_by_xLM <- list(mvlist=.modify_list(inits_by_xLM$mvlist,new_mvlist), # info. But maybe we won't need it ?
                                           beta_eta=beta_eta, # vector
                                           phi_est=.modify_list(inits_by_xLM$phi_est,new_phi_ests), # list
                                           lambdas=lambdas, # vector
                                           lambda=exp(mean(log(lambdas))) # scalar
      )
    } else {
      if (! is.null(processed$X_off_fn)) {
        X.pv <- environment(processed$X_off_fn)$X_off
      } else X.pv <- processed$AUGI0_ZX$X.pv
      processed$envir$inits_by_xLM <- .calc_inits_by_xLM(processed, X.pv=X.pv) # univariate response case: 
                                      # returns an environment with promises for deviance-dependent elements
    }
  }
  return(processed$envir$inits_by_xLM)
} 

