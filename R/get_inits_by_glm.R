# function called within HLfit for missing inits... 
# but also sometimes in fitme_body prior to optimization (cf complex condition for call of .eval_init_lambda_guess())
.get_inits_by_glm <- function(processed, family=processed$family, reset=FALSE) {
  if (reset || is.null(processed$envir$inits_by_glm)) { 
    ## if .get_inits_by_glm is called prior to optimization the family parameters may not be assigned, so: 
    if (family$family=="negbin") {
      checktheta <- suppressWarnings(try(environment(family$aic)$shape,silent=TRUE))
      if (inherits(checktheta,"try-error")) family <- Poisson(family$link, trunc=environment(family$aic)$trunc) 
    } else if (family$family=="COMPoisson") {
      checknu <- suppressWarnings(try(environment(family$aic)$nu,silent=TRUE))
      if (inherits(checknu,"try-error")) family <- poisson() ## do not use the loglambda link of the COMPoisson!
    }
    y <- processed$y ## requested by the formula
    if (family$family=="binomial" && NCOL(y)==1L) { 
      BinomialDen <- processed$BinomialDen ## requested by the formula
      begform <-"cbind(y,BinomialDen-y)~"  
    } else {begform <-"y~"}
    ###################################################if (pforpv==0) {endform <-"0"} else 
    X.pv <- processed$AUGI0_ZX$X.pv ## possibly requested by the formula
    pforpv <- ncol(X.pv)
    if(pforpv) {
      endform <-"X.pv-1" ## pas besoin de rajouter une constante vue qu'elle est deja dans X
    } else {
      if (family$family %in% c("binomial","poisson")) {
        endform <- "1" ## no meaningful glm without fixed effect in this case !
      } else {endform <- "0"}
    }
    locform <- as.formula(paste(begform, endform))
    prior.weights <- processed$prior.weights   ## do not try to eval() it outside of the .wfit function call; else nasty crashes may occur.
    n_lambda <- sum(attr(processed$ZAlist,"Xi_cols"))
    if (family$family %in% c("Gamma","gaussian")) {lam_fac <- 1/(n_lambda+1L)} else lam_fac <- 1/n_lambda
    resu <- list() 
    if (family$family=="gaussian" && family$link=="identity") {
      if (inherits(X.pv,"sparseMatrix")) {
        resglm <- .spaMM_lm.wfit(x=X.pv,y=y,offset=processed$off,w=eval(prior.weights))
      } else resglm <- lm.wfit(x=X.pv,y=y,offset=processed$off,w=eval(prior.weights))
      fv <- fitted(resglm)
      dev <- resid(resglm)^2
      if (! is.null(resglm$weights)) dev <- dev/resglm$weights
      guess <- sum(dev)/resglm$df.residual
      if (is.nan(guess)) { # resglm$df.residual=0, possibly wider issue with requested fit, cannot be resolved from here.
        resu$lambda <- resu$phi_est <- 1e-04
      } else resu$lambda <- resu$phi_est <- sum(dev)/resglm$df.residual # /lam_fac
    } else { ## GLM
      #
      resglm <- .tryCatch_W_E(glm.fit(x=X.pv, 
                                      y=processed$HLframes$Y, 
                                      weights = eval(prior.weights), 
                                      offset = processed$off, family = family, 
                                      control = processed[["control.glm"]]))$value 
      if (inherits(resglm,"error") || 
          ( ! resglm$converged && any(fitted(resglm)>1e20)) # this occurred in Gamma(log) models or negbin(log)
      ) {
        if (TRUE) {
          #control.glm <- processed[["control.glm"]]
          #control.glm$maxit <- control.glm$maxit*2
          resglm <- spaMM_glm.fit(x=X.pv, 
                                  y=processed$HLframes$Y, 
                                  weights = eval(prior.weights), 
                                  offset = processed$off, family = family, 
                                  control = processed[["control.glm"]])
          # F I X M E but this may still be poor; default control is still maxit=25. Think about (1) better diagnotic than any(fitted(resglm)>1e20) ? (2) check on non-convergence at this second resglm ? 
        } else {
          resglm <- withCallingHandlers(
            {
              #control.glm <- processed[["control.glm"]]
              #control.glm$maxit <- control.glm$maxit*2
              spaMM_glm.fit(x=X.pv, 
                            y=processed$HLframes$Y, 
                            weights = eval(prior.weights), 
                            offset = processed$off, family = family, 
                            control = processed[["control.glm"]])
            },
            warning = function(w){
              ## Ignore convergence diagnostics from *this* spaMM_glm.fit() call (only, as it is the first call within HLfit)
              #if(w$message != "glm.fit: algorithm did not converge" ) # the message is never this !
              warning(w$message)
              invokeRestart("muffleWarning")
            } 
          )
        }
        if ( ( ! .spaMM.data$options$spaMM_glm_conv_silent)
             && (conv_crit <- environment(spaMM_glm.fit)$spaMM_glm_conv_crit$max>0)) {
          resu$conv_info <- paste(".get_inits_by_glm() -> spaMM_glm.fit did not yet converge at iteration",
                                  resglm$iter,"(criterion:",
                                  paste(names(conv_crit),"=",signif(unlist(conv_crit),3),collapse=", "),").\n",
                                  "Use <fitting function>(., control.glm=list(maxit=<larger value>,..)) to control this.")
          assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit)) # So that no-convergence in these glms will not be warned about
        }
      }
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
    } 
    resu$converged <- resglm$converged
    #parent.env(environment(processed$get_inits_by_glm)) <- environment(stats::glm)
    processed$envir$inits_by_glm <- resu
  } 
  return(processed$envir$inits_by_glm)
} 

