.calc_inits_by_xLM <- function(processed, X.pv=processed$AUGI0_ZX$X.pv, family=processed$family,
                               y=processed$y, ## requested by the formula
                               Y=.get_from_terms_info(terms_info=processed$main_terms_info, which="Y"),
                               prior.weights=processed$prior.weights,   ## do not try to eval() it outside of the .wfit function call; else nasty crashes may occur.
                               off=processed$off
                               ) {
  ## Prior to optimization the family parameters may not be assigned, so: 
  orifamfam <- family$family
  if (orifamfam %in% c("negbin1","negbin2")) {
    if (variable_rdispar <- inherits(substitute(shape, env=environment(family$aic)),"call")) family <- Poisson(family$link, trunc=environment(family$aic)$trunc) 
  } else if (orifamfam=="beta_resp") {
    if (variable_rdispar <- inherits(substitute(prec, env=environment(family$aic)),"call")) { # only if beta_prec not assigned
      loc_link <- family$link
      family <- beta_resp(link=loc_link, prec=1) 
    }
  } else if (orifamfam=="COMPoisson") {
    loc_link <- family$link
    if (loc_link=="loglambda") loc_link <- "log"
    if (variable_rdispar <- inherits(substitute(nu, env=environment(family$aic)),"call")) family <- poisson(loc_link) # only if COMP_nu not assigned
    # problem with test glmmTMB: most of the time is spent on inits for glm when nu<1
    # if (inherits(nuchk <- substitute(nu, env=environment(family$aic)),"call") || nuchk >0.25) family <- poisson(loc_link)
    # : This helps for he GLM but not much otherwise. Other ideas ? __F I X M E__
  }
  locfamfam <- family$family
  if (locfamfam=="binomial" && NCOL(y)==1L) { 
    BinomialDen <- processed$BinomialDen ## requested by the formula
    begform <-"cbind(y,BinomialDen-y)~"  
  } else {begform <-"y~"}
  ###################################################if (pforpv==0) {endform <-"0"} else 
  pforpv <- ncol(X.pv)
  if( ! pforpv && # no predictor variable
      (locfamfam %in% c("binomial","poisson","beta_resp") || 
       # 1st condition => includes cases where original family was (T)negbin with free disp param (and pforpv=0)L, a *poisson* GLM *with* an Intercept is fitted (the resulting beta is ignored);
       # 2nd condition => if original family was negbin with FIXED disp param, a negbin GLM is fitted; it must have an intercept otherwise 
       #                                                  eta=0=> mu untruncated=1 => y>1 is impossible (and returns negative deviance=> negative init lambda)
       (locfamfam  %in% c("negbin1","negbin2") && # with FIXED shape
        family$zero_truncated) 
      )) X.pv <- matrix(1,ncol=1, nrow=nrow(X.pv))
  n_lambda <- sum(attr(processed$ZAlist,"Xi_cols"))
  if (locfamfam %in% c("Gamma","gaussian")) {lam_fac <- 1/(n_lambda+1L)} else lam_fac <- 1/n_lambda
  resu <- list() 
  
  if (locfamfam=="gaussian" && family$link=="identity") {
    
    if (inherits(X.pv,"sparseMatrix")) {
      resglm <- .spaMM_lm.wfit(x=X.pv,y=y,offset=off,w=eval(prior.weights))
    } else resglm <- lm.wfit(x=X.pv,y=y,offset=off,w=eval(prior.weights))
    fv <- fitted(resglm)
    dev <- resid(resglm)^2
    if (! is.null(resglm$weights)) dev <- dev/resglm$weights
    guess <- sum(dev)/resglm$df.residual
    if (is.nan(guess)) { # resglm$df.residual=0, possibly wider issue with requested fit, cannot be resolved from here.
      resu$lambda <- resu$phi_est <- 1e-04
    } else resu$lambda <- resu$phi_est <- sum(dev)/resglm$df.residual # /lam_fac
    
  } else if ( locfamfam %in% c("negbin1", "beta_resp")) { # with variable dispersion param; I cannot use $flags bc e.g. for COMPoisson, the local family may be stats::poisson
    
    resglm <- llm.fit(x=X.pv,  y=drop(Y),  weights = eval(prior.weights),  offset = off, 
                      family = family,  control = processed[["control.glm"]])
    if ( ! resglm$converged && family$link !="inverse") {
      ## divergence more likely than from spaMM_glm.fit, and fewer options then...  a simple control of 'start' may work well for log link:
      resglm <- llm.fit(x=X.pv,  y=drop(Y),  weights = eval(prior.weights),  offset = off, 
                        start = rep(0,ncol(X.pv)),
                        family = family,  control = processed[["control.glm"]])
    }
    if ( (! resglm$converged ) && ( ! .spaMM.data$options$xLM_conv_silent)
         && (conv_crit <- .spaMM.data$options$xLM_conv_crit$max>0)) { 
      resu$conv_info <- paste("In initialization, llm.fit() did not converge within",
                              resglm$iter,"iterations (criterion:",
                              paste(names(conv_crit),"=",signif(unlist(conv_crit),3),collapse=", "),").\n",
                              "Use <fitting function>(., control.glm=list(maxit=<larger value>,..)) to control this.")
      # No dispGammaGLM for residual variation in these families => no need to reset .spaMM.data$options$xLM_conv_crit here ? (see explanation in glm case)
      # BUT: cross effects in mv fits ? => play safe
      .spaMM.data$options$xLM_conv_crit <- list(max=-Inf) # So that no-convergence in later .calc_dispGammaGLM) call can be independently assessed
    }
    #
    overdisp <- as.numeric(deviance(resglm)/df.residual(resglm)) # that's the way the residual dispersion param is estimated in GLMs;
    # if (overdisp<1e-10) browser()
    # since resglm is not a glm nor an HLfit object, this uses the stats:::<>.default method of the two extractors.
    if (locfamfam=="beta_resp") { # Now, 'overdisp is 1 if lambda=0' We need lambda low in that case, but not too low (optim issues)
      resu$lambda <- max(0.1,min(log(overdisp),1)) #  1st beta_resp example of the doc is good test. 
      if (variable_rdispar) resu$beta_prec <- min(10, max(0.15, 1/(2*overdisp))) # will be ignored by rdisPars controls for *structured* beta_prec
    } else resu$lambda <- max(0.1,log(overdisp))
    
    ## I tried weird things here in version 4.0.14. In any case don't forget there is also a distinct .calc_fam_corrected_guess() function.

  } else { ## LOCfamily is GLM, even when fitted by obsInfo
    # (1) This handles truncated fams (since glm -> glm.fit handles Tpoisson() etc)
    # (2) This ignores obsInfo so results will differ from the equivalent LLF family
    if (locfamfam=="COMPoisson") glm.fit <- glm.nodev.fit 
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
      if ( ( ! resglm$converged) && ( ! .spaMM.data$options$xLM_conv_silent)
           && (conv_crit <- .spaMM.data$options$xLM_conv_crit$max>0)) {
        resu$conv_info <- paste("In initialization, spaMM_glm.fit() did not converge within",
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
    #
    overdisp <- as.numeric(deviance(resglm)/resglm$df.residual)
    if (orifamfam == "negbin1") {
      resu$NB_shape <- min(10, max(0.15, 2/overdisp)) # test: negbin1 in test-LLM, notably tnb1
    } else if (orifamfam == "negbin2") {
      resu$NB_shape <- min(20, max(0.15, 4/overdisp)) # test: spaMMintro's HLnegbin case; fit_05 of twinR back-compat check
    }
    resu$phi_est <- overdisp #/n_lambda 
    if (locfamfam=="binomial" && max(resglm$prior.weights)==1L) { ## binary response
      resu$lambda <- 1
    } else {
      # fv <- fitted(resglm)
      # resu$lambda <- sum(resid(resglm)^2/(resglm$prior.weights*family$variance(fv)))/resglm$df.residual # was a birth pang...
      # resid(resglm) gives resglm$residuals, which are the "working residuals" ie (y-mu)/mu.eta one the RHS of the IRLS equation
      # hence they are not an appropriate way to reconstruct the deviance residuals
      # which are given by with(resglm, family$dev.resids(y, mu=fitted.values, weights)), already taking prior weights into account
      # but this leads to 
      resu$lambda <- resu$phi_est
    }
  } ## end else (LOCfamily is 'GLM')
  if (pforpv) {
    ## Two potential problems (1) NA's pour param non estimables (cas normal); 
    ## (2) "glm.fit: fitted probabilities numerically 0 or 1 occurred" which implies separation or large offset
    if (max(abs(c(coefficients(resglm))),na.rm=TRUE)>1e10) { ## na.rm v1.2 
      warning(paste0("(!) Apparent divergence of estimates in a fixed-effect pre-fit.\n",
                     "    Check your data for extreme values, separation or bad offset values."))
    } 
    beta_eta <- c(coefficients(resglm)) ## this may include NA's. Testcase: HLfit(Strength ~ Material*Preheating+Method,data=weld)
    if (all(names(beta_eta)=="X.pv")) { ## si la formula etait y ~X.pv-1
      names(beta_eta) <- colnames(resglm$model$X.pv)
    } else names(beta_eta) <- gsub("^X\\.pv","",names(beta_eta)) ## removes initial "X.pv" without guessing any order or length
    resu$beta_eta <- beta_eta
  } else resu$beta_eta <- numeric(0) # for combining it in initialization of mv fit
  resu$converged <- resglm$converged
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
        new_mvlist[[mv_it]] <- .calc_inits_by_xLM(processed, X.pv=X.pv[resp_range,col_range, drop=FALSE],                                                   family=processed$families[[mv_it]],
                                                  y=processed$y[resp_range], ## requested by the formula
                                                  Y=.get_from_terms_info(terms_info=processed$main_terms_info, which="Y", mv_it=mv_it), 
                                                  prior.weights=processed$prior.weights[[mv_it]], ## do not try to eval() it outside of the .wfit function call; else nasty crashes may occur.
                                                  off=processed$off[resp_range] )
        beta_eta[col_range] <- new_mvlist[[mv_it]]$beta_eta
        lambdas[mv_it] <- new_mvlist[[mv_it]]$lambda
        if (processed$families[[mv_it]]$family %in% c("gaussian","Gamma")) {
          new_phi_ests[[mv_it]] <- new_mvlist[[mv_it]]$phi_est
        } else new_phi_ests[[mv_it]] <- 1 # no impact anyway?
      }
      names(new_mvlist) <- names(new_phi_ests) <- which(do_eval)
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
      processed$envir$inits_by_xLM <- .calc_inits_by_xLM(processed, X.pv=X.pv) # univariate response case
    }
  }
  return(processed$envir$inits_by_xLM)
} 

