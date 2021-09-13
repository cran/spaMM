.calc_inits_by_glm <- function(processed, X.pv=processed$AUGI0_ZX$X.pv, family=processed$family,
                               y=processed$y, ## requested by the formula
                               Y=.get_from_terms_info(terms_info=processed$main_terms_info, which="Y"),
                               prior.weights=processed$prior.weights,   ## do not try to eval() it outside of the .wfit function call; else nasty crashes may occur.
                               off=processed$off
                               ) {
  ## if .get_inits_by_glm is called prior to optimization the family parameters may not be assigned, so: 
  if (family$family=="negbin") {
    if (inherits(substitute(shape, env=environment(family$aic)),"call")) family <- Poisson(family$link, trunc=environment(family$aic)$trunc) 
  } else if (family$family=="COMPoisson") {
    loc_link <- family$link
    if (loc_link=="loglambda") loc_link <- "log"
    if (inherits(substitute(nu, env=environment(family$aic)),"call")) family <- poisson(loc_link)
    # problm with test glmmTMB: most of the time is spent on inits for glm when nu<1
    # if (inherits(nuchk <- substitute(nu, env=environment(family$aic)),"call") || nuchk >0.25) family <- poisson(loc_link)
    # : This helps for he GLM but not much otherwise. Other ideas ? __F I X M E__
  }
  
  if (family$family=="binomial" && NCOL(y)==1L) { 
    BinomialDen <- processed$BinomialDen ## requested by the formula
    begform <-"cbind(y,BinomialDen-y)~"  
  } else {begform <-"y~"}
  ###################################################if (pforpv==0) {endform <-"0"} else 
  pforpv <- ncol(X.pv)
  if( ! pforpv && # no predictor variable
      (family$family %in% c("binomial","poisson") || 
        # 1st condition => if original family was (T)negbin with free disp param (and pforpv=0)L, a *poisson* GLM *with* an Intercept is fitted (the resulting beta is ignored);
        # 2nd condition => if original family was negbin with FIXED disp param, a negbin GLM is fitted; it must have an intercept otherwise 
        #                                                  eta=0=> mu untruncated=1 => y>1 is impossible (and returns negative deviance=> negative init lambda)
        (family$family == "negbin" && # with FIXED shape
         family$zero_truncated)
    )) X.pv <- matrix(1,ncol=1, nrow=nrow(X.pv))
  n_lambda <- sum(attr(processed$ZAlist,"Xi_cols"))
  if (family$family %in% c("Gamma","gaussian")) {lam_fac <- 1/(n_lambda+1L)} else lam_fac <- 1/n_lambda
  resu <- list() 
  if (family$family=="gaussian" && family$link=="identity") {
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
  } else { ## GLM
    #
    if (family$family=="COMPoisson") glm.fit <- glm.nodev.fit 
    tryglm <- .tryCatch_W_E(glm.fit(x=X.pv, 
                                    y=Y, 
                                    weights = eval(prior.weights), 
                                    offset = off, family = family, 
                                    control = processed[["control.glm"]]))
    if (inherits((resglm <- tryglm$value),"error") || 
        ( ! resglm$converged && any(fitted(resglm)>1e20)) # this occurred in Gamma(log) models or negbin(log)
    ) {
      if (family$family=="gaussian" && family$link %in% c("log","inverse") && any(Y==0)) {
        ## gaussian(<"log"|"inverse">)$initialize is deficient
        family$initialize <- expression({
          if (any(y < 0)) 
            stop("negative values not allowed for the 'gaussian' family with such link")
          n <- rep.int(1, nobs)
          mustart <- y + 0.1
        })
      }
      resglm <- spaMM_glm.fit(x=X.pv, # ->LevM -> COMP bottleneck
                              y=Y, 
                              weights = eval(prior.weights), 
                              offset = off, family = family, 
                              control = processed[["control.glm"]])
      if ( ( ! .spaMM.data$options$spaMM_glm_conv_silent)
           && (conv_crit <- environment(spaMM_glm.fit)$spaMM_glm_conv_crit$max>0)) {
        resu$conv_info <- paste(".get_inits_by_glm() -> spaMM_glm.fit did not yet converge at iteration",
                                resglm$iter,"(criterion:",
                                paste(names(conv_crit),"=",signif(unlist(conv_crit),3),collapse=", "),").\n",
                                "Use <fitting function>(., control.glm=list(maxit=<larger value>,..)) to control this.")
        assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit)) # So that no-convergence in these glms will not be warned about
      }
    } else if ( ! is.null(tryglm$warning) && 
                substring(tryglm$warning$message,0,14)=="maxn truncated") .spaMM.data$options$COMP_maxn_warned <- FALSE
    ## -> see explanations in .COMP_maxn()
    #
    resu$phi_est <- as.numeric(deviance(resglm)/resglm$df.residual) #/n_lambda 
    if (family$family=="binomial" && max(resglm$prior.weights)==1L) { ## binary response
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
  } ## end else ('GLM')
  if (pforpv) {
    ## Two potential problems (1) NA's pour param non estimables (cas normal); 
    ## (2) "glm.fit: fitted probabilities numerically 0 or 1 occurred" which implies separation or large offset
    if (max(abs(c(coefficients(resglm))),na.rm=TRUE)>1e10) { ## na.rm v1.2 
      warning(paste0("(!) Apparent divergence of estimates in a *GLM* analysis of the data.\n",
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
.get_inits_by_glm <- function(processed, X.pv=processed$AUGI0_ZX$X.pv, family=processed$family, reset=FALSE) {
  if ( is.null(inits_by_glm <- processed$envir$inits_by_glm)) {
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
    new_mvlist <- inits_by_glm$mvlist
    new_phi_ests <- inits_by_glm$phi_ests 
    beta_eta <- inits_by_glm$beta_eta
    lambdas <- inits_by_glm$lambdas
  } else {
    family <- processed$family
    do_eval <- eval(reset)
  }
  #
  if (any(do_eval)) {
    if ( ! is.null(vec_nobs)) {
      cum_nobs <- attr(processed$families,"cum_nobs")
      cum_ncol_X <- attr(X.pv,"cum_ncol")
      for (mv_it in which(do_eval)) {
        resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
        col_range <- cum_ncol_X[mv_it]+seq_len(cum_ncol_X[mv_it+1L]-cum_ncol_X[mv_it]) # avoid (n+1:n) problem 
        new_mvlist[[mv_it]] <- .calc_inits_by_glm(processed, X.pv=X.pv[resp_range,col_range, drop=FALSE], 
                                                  family=processed$families[[mv_it]],
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
      processed$envir$inits_by_glm <- list(mvlist=.modify_list(inits_by_glm$mvlist,new_mvlist), # info. But maybe we won't need it ?
                                           beta_eta=beta_eta, # vector
                                           phi_est=.modify_list(inits_by_glm$phi_est,new_phi_ests), # list
                                           lambdas=lambdas, # vector
                                           lambda=exp(mean(log(lambdas))) # scalar
      )
    } else processed$envir$inits_by_glm <- .calc_inits_by_glm(processed)
  }
  return(processed$envir$inits_by_glm)
} 

