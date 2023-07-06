# called by .calc_inits_by_xLM() -> llm.fit() ; or by .calc_etaXLMblob()
# BinomialDen is distinctly handled between these two cases. 
.calc_dlogL_blob <- function(eta, mu, y, weights, family, phi, muetaenv, BinomialDen) {
  famfam <- family$family
  if (famfam=="COMPoisson")  {
    dlogLdmu <- family$DlogLDmu(mu=mu, y=y, thetaMuDerivs=muetaenv$thetaMuDerivs_2) 
    d2logLdmu2 <- family$D2logLDmu2(mu=mu, y=y, thetaMuDerivs=muetaenv$thetaMuDerivs_2)
  } else if (famfam=="betabin")  {
    if (missing(BinomialDen)) { # The missingness currently means that this has been called through llm.fit, 
      # AND that mu is muFREQS 
      BinomialDen <- environment(family$aic)$.BinomialDen 
      muFREQS <- mu
    } else muFREQS <- mu/BinomialDen # call from .calc_etaXLMblob: mu is muCOUNT, BinomialDen provided, $.BinomialDen not initialized
    dlogLdmu <- family$DlogLDmu(muFREQS=muFREQS, y, BinomialDen=BinomialDen, wt=weights) 
    d2logLdmu2 <- family$D2logLDmu2(muFREQS=muFREQS, y, BinomialDen=BinomialDen, wt=weights) 
  } else if (famfam=="binomial")  { 
    # if (missing(BinomialDen)) { # The missingness currently means that this has been called through llm.fit 
    #   ## AND that mu is muFREQS.
    # 
    #   ## glm(., family=spaMM:::.statsfam2LLF(stats::binomial), method="llm.fit") could lead here.
    #   ## fitme(., family=spaMM:::.statsfam2LLF(stats::binomial)) -> .calc_inits_by_xLM presumably 
    #   ## does not lead here bc .calc_inits_by_xLM should still call glm.fit() rather that llm.fit() for binomial models
    #
    #   BinomialDen <- environment(family$aic)$.BinomialDen #  It evaluates to NULL currently
    #   dlogLdmu <- family$DlogLDmu(muFREQS=mu, y, BinomialDen=BinomialDen, muCOUNT=mu*BinomialDen) 
    #   d2logLdmu2 <- family$D2logLDmu2(muFREQS=mu, y, BinomialDen=BinomialDen) 
    # } else { # call from .calc_etaXLMblob: mu is muCOUNT, BinomialDen provided, ...$.BinomialDen not initialized
      muFREQS <- mu/BinomialDen 
      dlogLdmu <- family$DlogLDmu(muFREQS=muFREQS, y, BinomialDen=BinomialDen, muCOUNT=mu) 
      d2logLdmu2 <- family$D2logLDmu2(muFREQS=muFREQS, y, BinomialDen=BinomialDen) 
    # }
  } else {
    dlogLdmu <- family$DlogLDmu(mu=mu, y=y, wt=weights, phi=phi) 
    d2logLdmu2 <- family$D2logLDmu2(mu=mu, y=y, wt=weights, phi=phi) 
  }
  dmudeta <- family$mu.eta(eta)
  d2mudeta2 <- family$D2muDeta2(eta)
  dlogcLdeta <- dlogLdmu*dmudeta
  d2logcLdeta2 <- d2logLdmu2*dmudeta^2+dlogLdmu*d2mudeta2
  list(dlogcLdeta=as.vector(dlogcLdeta), d2logcLdeta2=as.vector(d2logcLdeta2)) # as.vector() to drop attributes inherited from y or weights  
}


## Solve off weighted least square equation, for solution if damping=0, and else for change in solution of damped system (for LevM)
.llm_step_solver <- function(dlogL_blob, eta, offset, X, start, control, damping=0L) {
  dlogcLdeta <- dlogL_blob$dlogcLdeta
  d2logcLdeta2 <- dlogL_blob$d2logcLdeta2
  grad <- drop(crossprod(X, dlogL_blob$dlogcLdeta))
  if (any(d2logcLdeta2> -.Machine$double.eps )) { # positive d2logcLdeta2 but also exactly zero (Tnegbin(sqrt), with y=mu=1...)
    
    negHess <- crossprod(X, .Dvec_times_m_Matrix( - d2logcLdeta2, X))
    # it would be nice to devise a QR based solution here (using signs) EXCEPT that the solution of the exact system by QR may not be useful (no logL maxim)  
    cholH <- tryCatch(chol(negHess),error=function(e) e)
    if (inherits(cholH, "simpleError")) { 
      signs <- sign( - d2logcLdeta2)
      sqrtw_hess <- sqrt(abs(d2logcLdeta2))
      wX <- X * sqrtw_hess # x[good, , drop = FALSE] * sqrtw_hess
      if (damping) {
        LevM_result <- .LevenbergMsolveCpp(wX, grad, damping)
        LevM_result$rhs <- grad
        return(LevM_result)
      } else {
        start <-  get_from_MME.default(wX , szAug= sqrtw_hess* (eta-offset) + dlogcLdeta/sqrtw_hess, tol=min(1e-07, control$epsilon/1000))
        return(start)
      }
    } else if (damping) { # tested by # tested by fitme(cases~I(prop.ag/10)+offset(log(expec)),family=negbin1(shape=1/2.94), data=scotlip)
      dampDpD <- diag(negHess)*damping
      diag(negHess) <- diag(negHess)+dampDpD
      cholH <- chol(negHess)
      dbeta <- backsolve(cholH, backsolve(cholH, grad, transpose = TRUE)) # solve(mH, mH beta + grad) = beta + mH grad = beta - D2logLDbeta2 . DlogLDbeta
      return(list(dbetaV=drop(dbeta), dampDpD=dampDpD, rhs=grad)) 
    } else {
      Hbeta_grad <- crossprod(X,  - d2logcLdeta2 * (eta-offset))+ grad # (H beta = X' W X beta= X' W (eta-off)) + (grad = X .dlogcLdeta)
      beta <- backsolve(cholH, backsolve(cholH, Hbeta_grad, transpose = TRUE)) # solve(mH, mH beta + grad) = beta + mH grad = beta - D2logLDbeta2 . DlogLDbeta
      return(drop(beta)) # new beta, or dbeta, depending on rhs.
    }
    
  } else { # standard QR is safe
    sqrtw_hess <- sqrt( - d2logcLdeta2)
    wX <- X * sqrtw_hess # x[good, , drop = FALSE] * sqrtw_hess
    if (damping) {
      if (is.matrix(wX)) {
        LevM_result <- .LevenbergMsolveCpp(wX, grad, damping)
      } else LevM_result <- .LevenbergMsolve_Matrix(wX, grad, damping)
      LevM_result$rhs <- grad
      return(LevM_result)
    } else {
      # fit <-  .lm.fit(wX , sqrtw_hess* (eta-offset) + dlogcLdeta/sqrtw_hess, tol=min(1e-07, control$epsilon/1000))
      # if (any(!is.finite(fit$coefficients))) {
      #   conv <- FALSE
      #   warning(gettextf("non-finite coefficients at iteration %d", 
      #                    iter), domain = NA)
      #   break
      # }
      # if (nobs < fit$rank) 
      #   stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation", 
      #                         "X matrix has rank %d, but only %d observations"), 
      #                fit$rank, nobs), domain = NA)
      # if (!singular.ok && fit$rank < nvars) stop("singular fit encountered")
      ## calculate updated values of eta and mu with the new coef:
      ## get_from_MME_default with matrix and Matrix methods
      start <- get_from_MME_default(wX , szAug=sqrtw_hess* (eta-offset) + dlogcLdeta/sqrtw_hess, tol=min(1e-07, control$epsilon/1000))
      return(start)
    }
  }
}



llm.fit <- function (x, 
            y, # for glm.fit, y is "a vector of observations of length n" and here it is important this is a vector otherwise I may need a few more drop()"s.
            # But then the source code of glm.fit contains ynames <- if (is.matrix(y)) rownames(y) .. 
            weights = rep(1, nobs), 
            start = NULL, 
            etastart = NULL, mustart = NULL, offset = rep(0, nobs), family, 
            control = list(maxit=200), intercept = TRUE, singular.ok = TRUE) {
  control <- do.call("glm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) # cf comment on y argument... messy. (_F I X M E_)
    rownames(y)
  else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  
  if (is.null(weights)) weights <- rep.int(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  aic <- family$aic
  valideta <- family$valideta
  validmu <- family$validmu
  
  DlogLDmu <- family$DlogLDmu
  D2logLDmu2 <- family$D2logLDmu2
  D2muDeta2 <- family$D2muDeta2
  dev.resids <- family$dev.resids
  #
  # for binomial and betabin, family$initialize changes y 2 col -> 1 col; redefines y as a frequency and the weights as binomial weights
  # for betabin, this also redefines sat_logL to handle its BinomialDen argument 
  if (is.null(mustart)) {
    eval(family$initialize) 
  } else {
    mukeep <- mustart
    eval(family$initialize) 
    mustart <- mukeep
  }
  
  if (fam_is_COMP <- (family$family=="COMPoisson")) { # This *can* occur e.g.
    # glm(eggs ~ humans_eaten, family = COMPoisson(link=log, nu = 1), data = AliensMix, method="llm.fit")
    # Not that loglambda link is prevented bc .D2muDeta2() (etc) does not handle this link.
    muetaenv <- NULL
    dev.resids <- function(y, mu, wt) {family$dev.resids(y, mu, wt, muetaenv=muetaenv)}
    if (family$link=="loglambda") {
      linkinv <- function(eta) {family$linkinv(eta,muetaenv=muetaenv)}
      mu.eta <- function(eta) {family$mu.eta(eta,muetaenv=muetaenv)}
    } else {
      linkinv <- family$linkinv
      mu.eta <- family$mu.eta
    }
  } else {
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
  }
  if ( !is.function(linkinv)) 
    stop("'family' argument seems not to be a valid family object", 
         call. = FALSE)
  # Check whether beta will be constrained by the need for positive eta
  beta_bounded <- ( ( ! valideta(-1e-8) )|| # test not sufficient! next test needed for Gamma(inverse of identity)
                      ( ! validmu(linkinv(-1e-8)) ))
  n <- NULL ## to avoid an R CMD check NOTE which cannot see that n will be set by eval(family$initialize)
  
  coefold <- NULL
  
  ## Provide provisional eta; but see 1st iteration of loop which avoids the need for things like .gauss_initialize_in_Xbeta_image() used ones in spaMM_glm.fit.
  if (!is.null(etastart)) {
    eta <- etastart
  } else if (!is.null(start)) {
    if (length(start) != nvars) {
      stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                    nvars, paste(deparse(xnames), collapse = ", ")), domain = NA)
    } else {
      coefold <- start
      eta <- offset + drop(x %*% start)
    }
  } else eta <- drop(family$linkfun(mustart))
  eta <- .sanitize_eta(eta, y=y, family=family)   
  ##
  
  if (fam_is_COMP) muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
  mu <- linkinv(eta)
  if (!(validmu(mu) && valideta(eta))) 
    stop("cannot find valid starting values: please specify some", 
         call. = FALSE)
  
  devold <- sum(dev.resids(y, mu=mu, wt=weights))
  boundary <- conv <- FALSE
  damping <- 1e-7 ## as suggested by Madsen-Nielsen-Tingleff... # Smyth uses abs(mean(diag(XtWX)))/nvars
  dampingfactor <- 2
  ## the orginal glm.fit algo allows dmudeta=0 for 'individuals' with (prior)weight=0, and removes the other.
  ## This allows for changes in the dimensions of the design matrix over loop iterations... We remove this.
  goodinit <- weights > 0 
  for (iter in 1L:control$maxit) {
    ## this uses WLS in the first iteration and Levenberg Marquardt in later ones. # __F I X M E__ less systematic LevM ? LevM remains cheap for fied-effect models.
    if (iter==1L) {
      if (all(!goodinit)) {
        conv <- FALSE
        warning("no observations informative at iteration 1", domain = NA)
        break
      }  
      
      ## The aim here is to provide starting values of eta compatible with the model, rather than an unconstrained eta(mu=y),
      # so that initial dev is >= the ML deviance under the ftted model, suitable for LM iterations; 
      # while initialize() is typically mustart <- y, not  compatible with the model. => low control$epsilon necessary
      # Conversely, if all y's are quite small, a low control$epsilon may raise spurious warnings => inmportant to control it in .calc_dispGammaGLM()
      dlogL_blob <- .calc_dlogL_blob(eta, mu, y, weights, family, phi=1, muetaenv) # phi=1: either the dispersion param is distinct from phi and available in environment(family$aic)
      # or, for GLMs, only the relative phi values may matter for beta estimates. 
      # There should be no variation among phi values here, only possibly of prior weights.
      if (EMPTY) break # no model term except possibly an offset. return value still useful ( dev -> disperion estimate...)
      start <- .llm_step_solver(dlogL_blob, eta=eta, offset=offset, X=x, start=start, control=control, damping=0L)
      
      eta <- drop(x %*% start) + offset
      if ( ! beta_bounded) eta <- .sanitize_eta(eta, y=y, family=family, max=40) # else we use .get_valid_beta_coefs()
      if (fam_is_COMP) muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
      mu <- linkinv(eta) 
      dev <- suppressWarnings(sum(dev.resids(y, mu=mu, wt=weights)))
      boundary <- FALSE
      if (beta_bounded && 
          ( ! is.finite(dev) ||  ## NaN or Inf; 
            ! valideta(eta))             # I met the case sqrt link, eta<0 but dev OK (as function of mu=eta^2)) => stop("inner loop 2...)
          # I fixed .sanitize_eta in that case. But per se this does not protect from stop("inner loop 2; cannot correct...  below
          # so this test remains as a tempo fix (___F I X M E___)
      ) {
        
        if (is.null(coefold)) {
          if (requireNamespace("rcdd",quietly=TRUE)) {
            start <- .get_valid_beta_coefs(X=x,offset=offset,family,y,weights)
            eta <- drop(x %*% start) + offset
            if (fam_is_COMP) muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
            mu <- linkinv(eta) # could sanitize it, perhaps ?
            dev <- suppressWarnings(sum(dev.resids(y, mu=mu, wt=weights)))
          } else if ( ! identical(spaMM.getOption("rcdd_warned"),TRUE)) {
            message("llm.fit() failed. If the 'rcdd' package were installed, llm.fit() could automatically find a good starting value.")
            .spaMM.data$options$rcdd_warned <- TRUE
          }
        }
      }
      ## check for divergence
      if (!is.finite(dev)) {
        if (is.null(coefold)) 
          stop("no valid set of coefficients has been found: please supply starting values", 
               call. = FALSE)
        #warning("step size truncated due to divergence", call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit) 
            stop("inner loop 1; cannot correct step size", 
                 call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start) + offset
          if (family$family=="COMPoisson") muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
          mu <- linkinv(eta)
          dev <- suppressWarnings(sum(dev.resids(y, mu=mu, wt=weights)))
        }
        boundary <- TRUE
        if (control$trace) cat("Step halved: new deviance = ", dev, "\n", sep = "")
      }
      if (control$trace) cat("Deviance = ", dev, " (iteration ", iter, ")\n", sep = "")
      devold <- dev
      coef <- coefold <- start
    } else { # iter > 1L
      restarted <- FALSE
      
      while(TRUE) { # inner loop over damping values
        # LevenbergMstep_result <- .LevenbergMstepCallingCpp(wAugX=wX, LM_wAugz=LM_wz, damping=damping)
        ## LevenbergMstep_result includes 'rhs', the gradient of the objective, essentially  crossprod(wX,LM_wz)
        LevenbergMstep_result <- .llm_step_solver(dlogL_blob, eta=eta, offset=offset, X=x, start=start, control=control, damping=damping)
        #
        ## Use a function call to keep the change in beta->start,eta,mu,dev private:
        levMblob <- .eval_gain_dev_LevM(LevenbergMstep_result,family, x ,coefold,devold, offset,y,
                                        weights) # with default value of dev.resids argument => use the family$dev.resids function, 
        # not the llm.fit-local one => the llm.fit-local muetaenv won't be sought
        gainratio <- levMblob$gainratio
        
        if (beta_bounded) {
          conv_crit <- levMblob$conv_dev ## (dev crit in glm.fit style, vs $dbeta in auglinmodfit style)
          ## But D(dev) and dbetav may be small for non-zero gradient, part. with large damping.
          ## Hence additional check, worth when ( ! beta_bounded).
          ## It might make sense to make a stricter comparison on $conv_dev when (beta_bounded).
        } else { # we also check that the gradient vanish
          conv_crit <- max(levMblob$conv_dev, 
                           abs(LevenbergMstep_result$rhs)/(1+dev))
        }
        
        if (is.nan(gainratio)) {
          if (control$trace) cat("!")
          break
        } else if (gainratio>0) { ## success
          damping <- max(damping * max(1/3,1-(2*gainratio-1)^3), 1e-7) # lower bound as in .get_new_damping() for MMs    
          dampingfactor <- 2
          start <- levMblob$beta 
          if (fam_is_COMP) {
            muetaenv <- levMblob$muetaenv
            eta <- muetaenv$sane_eta
            mu <- muetaenv$mu
          } else {
            eta <- levMblob$eta ## FR->FR 1 -col matrix without names....
            mu <- levMblob$mu # muFREQS
          }
          mu <- levMblob$mu # muFREQS
          dev <- levMblob$dev 
          if (control$trace) cat("|")
          break
        } else if (iter>2
                   && gainratio==0) { # apparently flat deviance
          if (conv_crit < control$epsilon) { # we were at optimum
            damping <- 1e-7 
            if (control$trace) cat("#")
            break ## apparently flat dev
          } else { ## measurable gradient, but we don't move => too high damping ? (if damping from previous LevM step was too high)
            if ( ! restarted) { # condition to avoid infinite loop
              damping <- 1e-7
              dampingfactor <- 2
              restarted <- TRUE
              if (control$trace) cat("-")
              # and continue 
            } else {
              # hmm; well, break and diagnose...
              if (control$trace) cat("?")
              break
            }
          }
        } else { ## iter=1 or failure (gainratio<0): increase damping and continue iterations
          damping <- dampingfactor*damping
          dampingfactor <- dampingfactor*2
          if (control$trace) cat("+")
        } 
        if (damping>1e10) {
          break ## but an error may occur (NaN mu) (we should not reach this point except if gainratio<0?)
          #
          # A potential persistent pb is when there is a flat region of gainratio=0 for large damping, gainratio<0 for low damping,
          # and a narrow window of progress in between. No current example (a suspected one was a bug in .calc_dlogL_blob). 
          # A standard 1D optimizer would also easily miss it even when bracketing by bounds implied by previous iterations.
          # 
          # This breaks the inner loop. Alternatives could be:
          # * also break the outer loop (time gain, no drawback?)
          # * reset the damping controls? (In suspected example, they were reset in every other of the subsequent iterations of the outer loop)
          # Again, there is no current example good for testing altenatives. 
        }
      } ## end of inner loop over damping values
      
      if ( conv_crit < control$epsilon) { 
        conv <- TRUE
        coef <- start
        break
      } else {
        if (control$trace) 
          cat("Deviance = ", dev, " (iteration ", iter,", damping ", damping,", factor ", dampingfactor, ")\n", sep = "")
        devold <- dev # used first for next call to eval_gain_LM()
        coef <- coefold <- start
      }
    } ## => one WLS solve or Levenberg loop
    
    ## check for fitted values outside domain. (despite finite dev)
    if (!(valideta(eta) && validmu(mu))) { 
      if (is.null(coefold)) 
        stop("no valid set of coefficients has been found: please supply starting values", 
             call. = FALSE)
      warning("step size truncated: out of bounds", 
              call. = FALSE)
      ii <- 1
      while (!(valideta(eta) && validmu(mu))) {
        if (ii > control$maxit) stop("inner loop 2; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold)/2
        eta <- drop(x %*% start)
        eta <- eta + offset
        if (family$family=="COMPoisson") muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
        mu <- linkinv(eta)
      } ## stop()s or exits loop with valideta and mu
      boundary <- TRUE
      dev <- suppressWarnings(sum(dev.resids(y, mu=mu, wt=weights)))
      if (control$trace) cat("Step halved: new deviance = ", dev, "\n", sep = "")
    }
    
    ## update WLS system for next step 
    dlogL_blob <- .calc_dlogL_blob(eta, mu, y, weights, family, phi=1, muetaenv)
  } ## end main loop (either a break or control$maxit reached)
  
  if (! conv) {
    if (getOption("warn")==2L) { # immediate stop() ...
      warning("llm.fit() did not converge") # rather than stop() so it's tagged as converted from warning
    } else { # ... vs. *delayed* warning
      .spaMM.data$options$xLM_conv_crit <- list(max=max(.spaMM.data$options$xLM_conv_crit$max,conv_crit),latest=conv_crit)
    }
  }
  
  if (boundary) warning("spaMM_glm.fit: algorithm stopped at boundary value", call. = FALSE)
  eps <- 10 * .Machine$double.eps
  
  if (all(dlogL_blob$d2logcLdeta2<0)) {
    # regenerate the qr (etc) object.
    sqrtw_hess <- sqrt( - dlogL_blob$d2logcLdeta2)
    fit <-  .lm.fit(as.matrix(x * sqrtw_hess), sqrtw_hess* (eta-offset) + dlogL_blob$dlogcLdeta/sqrtw_hess, min(1e-07, control$epsilon/1000))
    # but the trouble is that evalGainLM may detect invalid LevM estimates, 
    # while the .lm.fit estimates (~damping=0) may be invalid (different if loop terminated with high damping)
    #start[fit$pivot] <- fit$coefficients
    #eta <- (x %*% start)[]
    #mu <- linkinv(eta <- eta + offset)
    #
    if (fit$rank < nvars) coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
    xxnames <- xnames[fit$pivot]
    # residuals <- (y - mu)/mu.eta(eta) # working residuals for a GLM
    fit$qr <- as.matrix(fit$qr)
    nr <- min(sum(goodinit), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
    }
    else Rmat <- fit$qr[1L:nvars, 1L:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
    
    # names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    names(y) <- ynames
    if (!EMPTY) 
      names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", nrow(x) - fit$rank))
    else linkinv(offset)
    rank <- if (EMPTY) {0} else {fit$rank}
    # derogative comment from glm.fit source:
    ##     ^^ is only initialize()d for "binomial" [yuck!]
    # n is in the envir of binomial() and set by binomial()$initialize() called above by eval(family$initialize) 
    fit_qr <- if ( ! EMPTY) structure(fit[c("qr", "rank", 
                                            "qraux", "pivot", "tol")], class = "qr")
    effects <- if ( ! EMPTY) fit$effects
    R <- if ( ! EMPTY) Rmat 
  } else { # tested by fitme(cases~I(prop.ag/10)+offset(log(expec)),family=negbin1(shape=1/2.94), data=scotlip)
    rank <- ncol(x)
    fit_qr <- effects <- R <- NULL
  }
  ##  residuals(., type="working"), analogous to (y - mu)/mu.eta(eta): 
  #      RHS of .lm.fit eqn for GLMs is [sqrt_]w {(etamo) + working-residuals}
  #      RHS of .lm.fit eqn is here sqrt_wHess (etamo) + dlogL_blob$dlogcLdeta/sqrt_wHess => working-res= dlogcLdeta/wHess, i.e.
  residuals <- - dlogL_blob$dlogcLdeta/dlogL_blob$d2logcLdeta2 
  wt <- - dlogL_blob$d2logcLdeta2 # speculative
  names(residuals) <- names(wt) <- names(weights) <- ynames
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  # wtdmu <- if (intercept) 
  #   sum(weights * y)/sum(weights) #   BUT generally best-fitting mu!=y
  # else linkinv(offset)
  # if (fam_is_COMP) { # The alternative code would use the existing muetaenv and wtdmu would be effectively ignored...
  #   nulldev <- sum(family$dev.resids(y, wtdmu, weights, muetaenv = NULL)) 
  # } else if (family$family %in% "betabin") { # The alternative code would use the existing muetaenv and wtdmu would be effectively ignored...
  #   nulldev <- NA # sum(family$dev.resids(y, wtdmu, weights, muetaenv = NULL)) # wtdmu should be a muFREQS... but which ? 
  # } else nulldev <- sum(family$dev.resids(y, mu=wtdmu, wt=weights))
  resdf <- n.ok - rank
  #
  if (family$family =="binomial") { ## for the sake of compatibility... base GLM code is what it is... 
    # all (base) aic functions have 'n' argument but only binomial appears to use it.
    # glm.fit calls $aic with five unnamed argument so all families that can be used with glm.fit() need these five arguments, dots being useless
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    # families that canNOT be used with glm.fit() can use a better API
  } else aic.model <- aic(y=y, mu=mu, wt=weights, dev=dev) + 2 * rank # with implicit default n=.BinomialDen for betabin
  # = hack so that .BinomialDen does not appear in the body of llm.fit...
  #
  list(coefficients = coef, 
       fitted.values = mu, 
       effects = effects, R = R, 
       rank = rank, qr = fit_qr, family = family, 
       linear.predictors = eta, deviance = sum(dev.resids(y, mu=mu, wt=weights)), aic = aic.model, 
       null.deviance = NA,   #   # Since best-fitting mu !=y, computing the null model would require a non-trivial fit...
       iter = iter, 
       residuals = residuals, # "working" weights
       weights = wt,  # 'GLM' weights
       prior.weights = weights, # as name says
       BinomialDen=environment(family$aic)$.BinomialDen, # used by simulate() for glm(., method="llm.fit")
       df.residual = resdf, df.null = nulldf, y = y, converged = conv, 
       boundary = boundary)
}


