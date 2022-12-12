.check_conv_glm_reinit <- function() {
  if ( ( ! .spaMM.data$options$spaMM_glm_conv_silent)
       && (conv_crit <- environment(spaMM_glm.fit)$spaMM_glm_conv_crit$max)>0) { 
    # This should arise from .calc_dispGammaGLM() ; $max potentially comes from any call to spaMM_glm.fit(), 
    # but .get_inits_by_glm() reinitializes it. Info from .get_inits_by_glm is tracked through processed$envir$inits_by_glm$conv_info -> HLfit's warningEnv
    warning(paste(".calc_dispGammaGLM() -> spaMM_glm.fit() did not always converge (criterion:",
                  paste(names(conv_crit),"=",signif(unlist(conv_crit),3),collapse=", "),")"))
    assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit))
  }
}

# The positive half line has vertex 0 and direction 1. vrepr would look like
# 0 1 0 # vertex (vrepr[,2]==1) 
# 0 0 1 # direction vector (vrepr[,2]==0)
# (the first col is not considered below)
# To define an initial value within it we start in (vertice=0 + epsilon*direction). 
# More generally, the Hull is defined as the set of linear combinations of the vertices or direction vectors defined by each row, 
# with weights all positive, and summing to 1 for the vertices; for direction vectors the hull goes arbitrarily far (weight large) in such directions.
.get_valid_beta_coefs <- function(X,offset,family,y,weights) {
  vrepr <- rcdd::scdd(rcdd::makeH(a1=-X,b1=offset))$output ## convex hull of feasible coefs
  if (nrow(vrepr)==0L) {
    stop("The model cannot be fitted: \n given 'x' and 'offset', only eta=0 might satisfy logical constraints.")
  } else {
    verticesRows <- (vrepr[,2]==1)
    betaS <- 0
    if (any(verticesRows)) {
      betaS <- t(vrepr[verticesRows, -c(1:2), drop = FALSE]) ## col vectors of boundary beta values
      betaS <- betaS + 0.001 * rowMeans(betaS) ## moving a bit inside
    }
    if (any( ! verticesRows)) {
      ## generate combinations of vertices and directions:
      dirS <- t(vrepr[ ! verticesRows,-c(1:2),drop=FALSE])
      dirS <- apply(dirS,2L, function(v) v/norm(matrix(v),type="2")) ## normalize direction vectors
      if (NCOL(dirS)>1L) {
        dirS <- 0.001* (dirS + rowMeans(dirS))
      } else dirS <- 0.002* dirS ## ~1+offset(...)
      if (is.matrix(betaS)) {
        betaS.cols <- split(betaS, col(betaS))
        betaS <- do.call(cbind,(lapply(betaS.cols,y=dirS,`+`)))
      } else betaS <- dirS + betaS ## as in Gamma(inverse) example 
    }
    if (is.array(betaS)) {
      etaS <- X %*% betaS ## col vectors of predicted etaS
      etaS <- etaS + offset
      muS <- family$linkinv(etaS)
      devS <- apply(muS,2, function(mu) sum(family$dev.resids(y, mu=mu, weights)))
      beta <- betaS[,which.min(devS)]
      return(beta)
    } else return(betaS) ## Vrepr was a single line with  0    0    1
  }
}

# For gaussian response, there is no formal constraint on model coefficients, but good starting values may be necessary to avoid overflow,
# and for the LevM method it is important that mustart belongs for the model family. So this code produces 'good' beta and deduces mu from it.
# This is distinct from .get_valid_beta_coefs() that handles formal constraint on beta, that eta must be positive.
.get_init_in_Xbeta_image <- function(X, #  design matrix
                                     offset, # true model offset, not a hack
                                     family, # conceived for gaussian(log)
                                     y, # response variable
                                     weights, # prior 'weights' of the GLM call
                                     which="mu"
) {
  link <- family$link
  linky <- switch(link,
                  "log" = log(pmax(1e-05, y)),
                  "inverse" = {
                    linky <- 1/y
                    linky[y==0] <- 1e5
                    linky
                  },
                  "logit" = pmin(1-1e-4, pmax(1e-4, y)),
                  "identity" = stop("This function should not be called for 'identity' link"), # we should not need it, as least for gaussian
                  stop("This function was not conceived for this link")
  )
  # The non-rcdd method may be sufficient but the rcdd-based one aims to reduce the  risk of overflow. Both method work on linky:='linkfun(y)' (nearly). 
  if (requireNamespace("rcdd",quietly=TRUE)) {
    # This code seeks the hull of coefficients that produces mu >= linky, 
    # uses its vertices (which are beta such that some mu=linky) as candidate values of beta, 
    # and retains the beta that has thee lowest deviance (for true y and true link)
    offset <- offset -linky 
    vrepr <- rcdd::scdd(rcdd::makeH(a1=-X,b1=offset))$output ## convex hull of feasible coefs
    verticesRows <- (vrepr[, 2] == 1)
    betaS <- t(vrepr[verticesRows, -c(1:2), drop = FALSE]) 
    etaS <- X %*% betaS ## col vectors of predicted etaS
    etaS <- etaS + offset
    muS <- family$linkinv(etaS)
    devS <- apply(muS,2, function(mu) sum(family$dev.resids(y, mu=mu, weights))) 
    if (which=="mu") {
      return(muS[,which.min(devS)])
    } else return(betaS[,which.min(devS)])
  } else {
    # if ( ! identical(spaMM.getOption("rcdd_warned"),TRUE)) {} # message not useful here? 
    fit <- lm.wfit(x=X, y=linky, weights=weights)
    if (which=="mu") {
      eta <- stats::fitted.values(fit)
      return(family$linkinv(eta))
    } else {
      return(fit$coefficients) # unaffected by pivoting (they are named so that would be easy to check)
    }
  }
  return(beta)
}

.gauss_initialize_in_Xbeta_image <- expression({
  ## original gaussian(<"log"|"inverse">)$initialize fails to produce starting values in some cases... BUT...
  n <- rep.int(1, nobs) # must be for the $aic call. stats::glm is such a weird thing. 
  # if (is.null(etastart) && is.null(start) && is.null(mustart) && 
  #     ((family$link == "inverse" && any(y == 0)) || 
  #      (family$link == "log" && any(y <= 0)))) 
  #   stop("cannot find valid starting values: please specify some")
  # mustart <- y
  # ... BUT moreover they are generally not in X.beta image, which is not OK for starting LevM algo.
  # The ok links are c("inverse", "log", "identity") and we should not need LevM for identity link in which case we ignore this constraint
  if (is.null(etastart) && is.null(start) && is.null(mustart))
    if (family$link=="identity") {
      mustart <- y
    } else if (ncol(x)) {
      mustart <- .get_init_in_Xbeta_image(X=x, offset=offset, family=family, y=y, weights=weights, which="mu")
    } else mustart <- family$linkinv(offset) # gaussian( ! identity ) with fully fixed beta coeffs
})


# same algo as .eval_gain_clik_LevM() but with different input, use of family$dev.resids(mu=muFREQS) instead of .calc_clik(mu=muCOUNT)...
# Only used to fit fixed-effect models
.eval_gain_dev_LevM <- function(LevenbergMstep_result,family, x ,coefold,devold, offset,y,weights,
                                      dev.resids=family$dev.resids) { 
  dbeta <- LevenbergMstep_result$dbetaV
  beta <- coefold + dbeta
  eta <- drop(x %*% beta) + offset
  eta <- .sanitize_eta(eta, y=y, family=family, max=40) 
  if (family$family=="COMPoisson") {
    muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
    mu <- muetaenv$mu
    dev <- suppressWarnings(sum(dev.resids(y, mu=mu, wt=weights, muetaenv=muetaenv)))
  } else {
    mu <- family$linkinv(eta) # muFREQS
    dev <- suppressWarnings(sum(dev.resids(y, mu=mu, wt=weights)))
    muetaenv <- NULL
  }
  # for binomial fam, $family(initialize) has redefined y as a frequency and the weights as binomial weights
  if (is.infinite(dev) || is.na(dev)) {  
    gainratio <- -1
    conv_dev <- Inf
  } else {
    summand <- dbeta*(LevenbergMstep_result$rhs+ LevenbergMstep_result$dampDpD * dbeta) 
    ## In the summand, all terms should be positive. conv_dbetaV*rhs should be positive. 
    # However, numerical error may lead to <0 or even -Inf
    #  Further, if there are both -Inf and +Inf elements the sum is NaN and the fit fails.
    summand[summand<0] <- 0
    denomGainratio <- sum(summand)
    ddev <- devold-dev
    conv_dev <- abs(ddev)/(1+abs(dev))
    gainratio <- 2*ddev/denomGainratio 
  }
  return(list(gainratio=gainratio,dev=dev,beta=beta,eta=eta,
              mu=mu, # muFREQS
              muetaenv=muetaenv,
              conv_dev=conv_dev))
}  


# In a fitme calls this can be called by .get_inits_from_fit(), and by get_glm_phi) -> .calc_dispGammaGLM()
spaMM_glm.fit <- local({
  spaMM_glm_conv_crit <- list(max=-Inf)
  function (x, y, weights = rep(1, nobs), 
            start = NULL, ## beta coefs... equivalent info such as mustart in Xbeta image) is needed for the LevenbergM algo: 
                          ## either provided in the call or by some self-starting function called below (such as 
                          ## .gauss_initialize_in_Xbeta_image; default default family$initialize won't be in Xbeta image generally).
            etastart = NULL, mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
            control = list(maxit=200), intercept = TRUE, singular.ok = TRUE) 
  {
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) 
      rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
      weights <- rep.int(1, nobs)
    if (is.null(offset)) 
      offset <- rep.int(0, nobs)
    aic <- family$aic
    valideta <- family$valideta
    validmu <- family$validmu
    if (family$family=="COMPoisson") {
      attr(y,"CMP") <- .CMP_attr_y(y, family)
      muetaenv <- NULL
      variance <- function(mu) {family$variance(mu,muetaenv=muetaenv)}
      dev.resids <- function(y, mu, wt) {family$dev.resids(y, mu, wt, muetaenv=muetaenv)}
      if (family$link=="loglambda") {
        linkinv <- function(eta) {family$linkinv(eta,muetaenv=muetaenv)}
        mu.eta <- function(eta) {family$mu.eta(eta,muetaenv=muetaenv)}
        beta_bounded <- FALSE
      } else {
        linkinv <- family$linkinv
        mu.eta <- family$mu.eta
        beta_bounded <- ( ( ! valideta(-1e-8) )|| # test not sufficient! next test needed for Gamma(inverse of identity)
                            ( ! validmu(linkinv(-1e-8)) ))
      }
    } else {
      variance <- family$variance
      dev.resids <- family$dev.resids
      linkinv <- family$linkinv
      mu.eta <- family$mu.eta
      if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object", 
             call. = FALSE)
      # Check whether beta will be constrained by the need for postive eta
      beta_bounded <- ( ( ! valideta(-1e-8) )|| # test not sufficient! next test needed for Gamma(inverse of identity)
                          ( ! validmu(linkinv(-1e-8)) ))
    }
    #
    # delayedAssign("positive_eta", {
    #   pos_eta <- (family$family=="Gamma" && family$link %in% c("identity","inverse"))
    #   (pos_eta || (family$family %in% c("poisson","negbin2") && family$link %in% c("identity","sqrt"))) 
    # })
    n <- NULL ## to avoid an R CMD check NOTE which cannot see that n will be set by eval(family$initialize)
    if (family$family=="gaussian" && family$link!="identity") family$initialize <- .gauss_initialize_in_Xbeta_image
    if (is.null(mustart)) {
      eval(family$initialize) ## changes y 2 col -> 1 col 
      
      ## if y=1e-200, this sets mustart=1e-200, and if later etastart is NULL, then eta <- family$linkfun(mustart)
      ##    followed by mu <- linkinv(eta) gives mu with minimal value ~1e-16
      ## This comes from Gamma(log)$linkinv being function (eta) pmax(exp(eta), .Machine$double.eps)
      ## By contrast Gamma()$linkinv is function (eta) 1/eta , and really low y yields really low mu
      
      ## Option "Gamma_min_y" may prevent extreme values except in direct call of spaMM_GLM.fit
    }
    else {
      mukeep <- mustart
      eval(family$initialize) ## or binomial fam, changes y 2 col -> 1 col; redefines y as a frequency and the weights as binomial weights
      mustart <- mukeep
    }
    if (EMPTY) {
      eta <- rep.int(0, nobs) + offset
      if (!valideta(eta)) 
        stop("invalid linear predictor values in empty model", 
             call. = FALSE)
      if (family$family=="COMPoisson") muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
      mu <- linkinv(eta)
      if (!validmu(mu)) 
        stop("invalid fitted means in empty model", call. = FALSE)
      dev <- suppressWarnings(sum(dev.resids(y, mu, weights)))
      w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
      good <- rep_len(TRUE, length(mu))
      boundary <- conv <- TRUE
      coef <- numeric()
      iter <- 0L
    } else {
      coefold <- NULL
      if (!is.null(etastart)) {
        eta <- etastart
      } else if (!is.null(start)) {
        if (length(start) != nvars) {
          stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                        nvars, paste(deparse(xnames), collapse = ", ")), domain = NA)
        } else {
          coefold <- start
          eta <- offset + as.vector(if (NCOL(x) == 1L) {x * start} else {x %*% start})
        }
      } else eta <- family$linkfun(mustart)
      # if (family$link=="log") {
      #   eta <- .sanitize_eta_log_link(eta, max=40,y=y, warn_neg_y= (family$family !="gaussian"))
      # } else if (family$link=="loglambda") {
      #   COMP_nu <- environment(family$aic)$nu 
      #   eta <- .sanitize_eta_log_link(eta, max=40, y=y, nu=COMP_nu)
      # }
      eta <- .sanitize_eta(eta, y=y, family=family)   
      if (family$family=="COMPoisson") muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
      mu <- linkinv(eta)
      if (!(validmu(mu) && valideta(eta))) 
        stop("cannot find valid starting values: please specify some", 
             call. = FALSE)
      devold <- sum(dev.resids(y, mu, weights))
      boundary <- conv <- FALSE
      damping <- 1e-7 ## as suggested by Madsen-Nielsen-Tingleff... # Smyth uses abs(mean(diag(XtWX)))/nvars
      dampingfactor <- 2
      goodinit <- weights > 0 
      for (iter in 1L:control$maxit) {
        ## this uses WLS in the first iteration and Levenberg Marquardt in later ones. 
        if (iter==1L) {
          ## The aim here is to provide starting values of eta compatible with the model, rather than an unconstrained eta(mu=y),
          # so that initial dev is >= the ML deviance under the ftted model, suitable for LM iterations; 
          # while initialize() is typically mustart <- y, not  compatible with the model. => low control$epsilon necessary
          # Conversely, if all y's are quite small, a low control$epsilon may raise spurious warnings => inmportant to control it in .calc_dispGammaGLM()
          varmu <- variance(mu)[goodinit]
          if (anyNA(varmu)) stop("NAs in V(mu): consult the package maintainer.")
          if (any(varmu == 0)) 
            stop("0s in V(mu): consult the package maintainer.")
          mu.eta.val <- mu.eta(eta)
          if (anyNA(mu.eta.val[goodinit])) 
            stop("NA/NaN in d(mu)/d(eta): consult the package maintainer.")
          good <- goodinit & (mu.eta.val != 0)
          if (all(!good)) {
            conv <- FALSE
            warning("no observations informative at iteration 1", domain = NA)
            break
          }  
          z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
          w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
          if (anyNA(z)) stop("NA/NaN in 'z': consult the package maintainer.")
          if (anyNA(w)) stop("NA/NaN in 'w': consult the package maintainer.") # suggests too large 'mu'
          # .lm.fit is a wrapper for C_Cdqrls, but with check=TRUE
          fit <-  .lm.fit(x[good, , drop = FALSE] * w, z * w, min(1e-07, control$epsilon/1000))
          if (any(!is.finite(fit$coefficients))) {
            conv <- FALSE
            warning(gettextf("non-finite coefficients at iteration %d", 
                             iter), domain = NA)
            break
          }
          if (nobs < fit$rank) 
            stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation", 
                                  "X matrix has rank %d, but only %d observations"), 
                         fit$rank, nobs), domain = NA)
          if (!singular.ok && fit$rank < nvars) stop("singular fit encountered")
          ## calculate updated values of eta and mu with the new coef:
          start[fit$pivot] <- fit$coefficients
          eta <- drop(x %*% start) + offset 
          if ( ! beta_bounded) eta <- .sanitize_eta(eta, y=y, family=family, max=40) # else we use .get_valid_beta_coefs()
          if (family$family=="COMPoisson") muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
          mu <- linkinv(eta) # automatically using the muetaenv in the COMPoisson case (cf locally defined linkinv fn)
          dev <- suppressWarnings(sum(dev.resids(y, mu, weights)))
          boundary <- FALSE
          if (beta_bounded && !is.finite(dev)) { ## NaN or Inf
            if (is.null(coefold)) {
              if (requireNamespace("rcdd",quietly=TRUE)) {
                start <- .get_valid_beta_coefs(X=x,offset=offset,family,y,weights)
                eta <- drop(x %*% start) + offset
                if (family$family=="COMPoisson") muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
                mu <- linkinv(eta) # could sanitze it, perhaps ?
                dev <- suppressWarnings(sum(dev.resids(y, mu, weights)))
              } else if ( ! identical(spaMM.getOption("rcdd_warned"),TRUE)) {
                message("glm() failed. If the 'rcdd' package were installed, spaMM_glm.fit() could automatically find a good starting value.")
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
              eta <- (x %*% start)[]
              eta <- eta + offset
              if (family$family=="COMPoisson") muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
              mu <- linkinv(eta)
              dev <- suppressWarnings(sum(dev.resids(y, mu, weights)))
            }
            boundary <- TRUE
            if (control$trace) 
              cat("Step halved: new deviance = ", dev, "\n", 
                  sep = "")
          }
          if (control$trace) 
            cat("Deviance = ", dev, " (iteration ", iter, ")\n", sep = "")
          devold <- dev
          coef <- coefold <- start
        } else {
          restarted <- FALSE
          while(TRUE) {
            if (inherits(wX,"Matrix")) { ## maybe never the case until upstream code is modified
              LevenbergMstep_result <- .LevenbergMsolve_Matrix(wAugX=wX, LM_wAugz=LM_wz, damping=damping)
            } else LevenbergMstep_result <- .LevenbergMstepCallingCpp(wAugX=wX, LM_wAugz=LM_wz, damping=damping)
            ## LevenbergMstep_result contains rhs, the gradient of the objective, essentially  crossprod(wX,LM_wz)
            #
            ## Use a function call to keep the change in beta->start,eta,mu,dev private:
            levMblob <- .eval_gain_dev_LevM(LevenbergMstep_result,family, x ,coefold,devold, offset,y,weights) # cf remark in llm.fit on using default family$dev.resids here
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
              if (family$family=="COMPoisson") {
                muetaenv <- levMblob$muetaenv
                eta <- muetaenv$sane_eta
                mu <- muetaenv$mu
              } else {
                eta <- levMblob$eta ## FR->FR 1 -col matrix without names....
                mu <- levMblob$mu # muFREQS
              }
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
            if (damping>1e10) break ## but an error may occur (NaN mu) (we should not reach this point except if gainratio<0?)
          } ## end of while(TRUE) {...}
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
        } ## .lm.fit or Levenberg step
        ## check for fitted values outside domain. (despite finite dev)
        if (!(valideta(eta) && validmu(mu))) { 
          if (is.null(coefold)) 
            stop("no valid set of coefficients has been found: please supply starting values", 
                 call. = FALSE)
          warning("step size truncated: out of bounds", 
                  call. = FALSE)
          ii <- 1
          while (!(valideta(eta) && validmu(mu))) {
            if (ii > control$maxit) 
              stop("inner loop 2; cannot correct step size", 
                   call. = FALSE)
            ii <- ii + 1
            start <- (start + coefold)/2
            eta <- drop(x %*% start)
            eta <- eta + offset
            if (family$family=="COMPoisson") muetaenv <- .CMP_muetaenv(family, pw=weights, eta)
            mu <- linkinv(eta)
          } ## stop()s or exits loop with valideta and mu
          boundary <- TRUE
          dev <- suppressWarnings(sum(dev.resids(y, mu, weights)))
          if (control$trace) 
            cat("Step halved: new deviance = ", dev, "\n", 
                sep = "")
        }
        ## 
        varmu <- variance(mu)[goodinit]
        if (anyNA(varmu)) stop("NAs in V(mu)")
        if (any(varmu == 0)) 
          stop("0s in V(mu): consult the package maintainer.")
        mu.eta.val <- mu.eta(eta)
        if (anyNA(mu.eta.val[goodinit])) 
          stop("NA/NaN in d(mu)/d(eta): consult the package maintainer.")
        good <- goodinit & (mu.eta.val != 0)
        if (all(!good)) {
          conv <- FALSE
          warning(gettextf("no observations informative at iteration %d", 
                           iter+1L), domain = NA)
          break
        }
        ##
        z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
        w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        if (anyNA(z)) stop("NA/NaN in 'z': consult the package maintainer.")
        if (anyNA(w)) stop("NA/NaN in 'w': consult the package maintainer.") # suggests too large 'mu'
        wX <- .calc_wAugX(XZ_0I=x[good, , drop = FALSE],sqrt.ww=w)
        LM_wz <- z*w - (wX %*% coefold) ## FIXME? if coefold diverges, LM_wz diverges and increasing damping may not be sufficient
      } ######### end MAIN loop (either a break or control$maxit reached)
      if (any(good) & ! conv) {
        if (getOption("warn")==2L) { # immediate stop() ...
          warning("spaMM_glm.fit() did not converge") # rather than stop() so it's tagged as converted from warning
        } else { # ... vs. *delayed* warning
          spaMM_glm_conv_crit <<- list(max=max(spaMM_glm_conv_crit$max,conv_crit),latest=conv_crit)
        }
      }
      if (boundary) 
        warning("spaMM_glm.fit: algorithm stopped at boundary value", 
                call. = FALSE)
      eps <- 10 * .Machine$double.eps
      if (family$family == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps)) 
          warning("spaMM_glm.fit: fitted probabilities numerically 0 or 1 occurred", 
                  call. = FALSE)
      }
      if (family$family == "poisson") {
        if (any(mu < eps)) 
          warning("spaMM_glm.fit: fitted rates numerically 0 occurred", 
                  call. = FALSE)
      }
      # regenerate the qr (etc) object.
      fit <-  .lm.fit(x[good, , drop = FALSE] * w, z * w, tol=min(1e-07, control$epsilon/1000))
      # but the trouble is that evalGainLM may detect invalid LevM estimates, 
      # while the .lm.fit estimates (~damping=0) may be invalid (different if loop terminated with high damping)
      #start[fit$pivot] <- fit$coefficients
      #eta <- (x %*% start)[]
      #mu <- linkinv(eta <- eta + offset)
      #
      if (fit$rank < nvars) 
        coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
      xxnames <- xnames[fit$pivot]
      fit$qr <- as.matrix(fit$qr)
      nr <- min(sum(good), nvars)
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
    }
    residuals <- (y - mu)/mu.eta(eta)
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY) 
      names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", 
                                                                  sum(good) - fit$rank))
    wtdmu <- if (intercept) 
      sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) 
      0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    # derogative comment from glm.fit source:
    ##     ^^ is only initialize()d for "binomial" [yuck!]
    # n is in the envir of binomial() and set by binomial()$initialize() called above by eval(family$initialize) 
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
         effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
         rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank", 
                                                       "qraux", "pivot", "tol")], class = "qr"), family = family, 
         linear.predictors = eta, deviance = dev, aic = aic.model, 
         null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
         df.residual = resdf, df.null = nulldf, y = y, converged = conv, 
         boundary = boundary)
  }
})

spaMM_glm <- function(formula, family = gaussian, data, weights, subset,
                      na.action, start = NULL, etastart, mustart, offset,
                      control = list(...), model = TRUE, method = c("glm.fit","spaMM_glm.fit"),
                      ## new 'method's such as spaMM_glm.fit must have the same arguments as R's default method.  
                      #  When the latter are changed the new method fails. Cf introduction of singular.ok in R-devel, Sept. 2017  
                      x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, strict=FALSE,...) {
  mc <- match.call(expand.dots=TRUE)
  mc$strict <- NULL
  ## This code should not interfere with processed$family with possibly assigned param
  family <- .checkRespFam(family)
  summaryfamily <- .as_call_family(family,get_param = TRUE) ## only for printing
  mc[[1L]] <- quote(stats::glm) # for get("glm", asNamespace("stats")), summary.glm() will print the body of glm() ! ./.
  # quote seems OK for base fns; cf usages of mc[[1L]] <- quote(stats::model.frame). othewise change the final res$call[[1L]] 
  mc$method <- method[1L] ## method must have 2 elements (doc)
  res <- .tryCatch_W_E(eval(mc,parent.frame()))
  if (inherits(res$value,"error")) {
    if ( requireNamespace("rcdd",quietly=TRUE) ) {
      mc$method <- method[2L]
      res <- eval(mc,parent.frame())
      res$call$family <- summaryfamily
      return(res)
    } else if ( ! identical(spaMM.getOption("rcdd_warned"),TRUE)) {
      stop("spaMM_glm() -> glm() failed. If the 'rcdd' package were installed, spaMM_glm() could fit the model.")
      .spaMM.data$options$rcdd_warned <- TRUE
    } 
  } else {
    if (! is.null(res$warning) &&
        ## some glm warnings are useful only to understand a failure, hence not useful here 
        # To catch specific warnings, one has to catch their translations (as its not clear how to inhibit translations cleanly)
        # If stats::lm was the called function then this seems the way: 
        res$warning$message != gettext("step size truncated due to divergence",domain="R-stats") ) {
      if (res$warning$message == gettext("glm.fit: algorithm did not converge",domain="R-stats") && strict) { # by default, strict=FALSE and non-converged result is returned.
        mc$method <- method[2L] ## changes the method arg, but the called function is still stats::glm
        res <- eval(mc,parent.frame())
        res$call$family <- summaryfamily
        return(res)
      } else warning(res$warning)
    }
    if ((conv_crit <- environment(spaMM_glm.fit)$spaMM_glm_conv_crit$max>0)) {
      warning(paste("spaMM_glm() -> spaMM_glm.fit did not yet converge at iteration",res$iter,"(criterion:",
                    paste(names(conv_crit),"=",signif(unlist(conv_crit),3),collapse=", "),")"))
      assign("spaMM_glm_conv_crit",list(max=-Inf) , envir=environment(spaMM_glm.fit))
    }
    res$value$call$family <- summaryfamily
    res$value
  }
}