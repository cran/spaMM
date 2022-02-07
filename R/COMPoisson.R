.CMP_use_asympto <- function(pow_lam_nu, nu, lambda) {
  ## Gaunt et al suggest lambda>1.5 and pow_lam_nu>1.5
  ## this ensures that nu*pow_lam_nu>1.1 but the expansion is clearly better for nu*pow_lam_nu >> 1 
  return(eval(.spaMM.data$options$CMP_asympto_cond))
}

.COMP_maxn <- function(lambda,nu) { ## this uses the simplest asympto approx to determine the discrete sum bound 
  # most pbatic case is nu->0, lambda->1^- where E -> lambda/(1-lambda) ->Inf while the first En approx <1 or <<1
  Pr_moments <- .COMP_Pr_moments(pow_lam_nu=lambda^(1/nu), lambda, nu, moments=c(1,2),use_asympto=TRUE)
  res <- app_En <- Pr_moments[1]
  opt_maxn <- .spaMM.data$options$COMP_maxn
  if (res< opt_maxn) {
    app_En2 <- Pr_moments[2]
    app_Vn <- app_En2-app_En^2 # should be positive when nu>1, see comment in source fro the asympto approx
    # BUT: app_En can be negative for lambda <<1, nu moderately > 1
    app_Vn <- max(app_Vn,lambda^(1/nu)/nu)
    app_En <- max(app_En,lambda^(1/nu))
    # the max() ensures continuity in nu=1 as the .COMP_Pr_moments equals the simpler approx there  
    res <- 4+app_En+6*sqrt(app_Vn)  ## real, to allow continuity correction 
    # 6*sqrt(app_Vn) to keep all but 1e-6 of the distrib (Feller p.245). 
    # But in fact this uses Guassian approx and fails for e.g. y=0, mu=0.09568245 => visible diff between dpois and .dCOMP
    # => add 4. 
    #   For same effect with geom maxn=-6/log10(lambda) and with log(lambda)~lambda-1 we need 6*sqrt(app_Vn)~ 6/[log(10)*(1-lambda)]
    # To ensure continuity ( & with log10() argument always <1):
    if (nu<1 && lambda<1) res <- -((1-nu)^2)*6/log10(nu+(1-nu)*lambda)+(1-(1-nu)^2)*res  ## real, to allow continuity correction 
    # => in nu->1 the correction is zero for all 1>lambda>0, in nu->0 res is -6/log10(lambda)
  }
  if (res>opt_maxn) {
    res <- opt_maxn
    if ( ! identical(.spaMM.data$options$COMP_maxn_warned,TRUE)) {
      warning(paste0("maxn truncated to ",res," for (lambda,nu)= (",lambda,",",nu,") and possibly other values."))
      .spaMM.data$options$COMP_maxn_warned <- TRUE
      ### This is first called through 
      # .get_inits_by_glm(processed, family = family, reset = quote(family$family %in% c("COMPoisson", "negbin")))
      # -> .calc_inits_by_glm(processed)
      # -> .tryCatch_W_E(glm.fit ...)
      ### so the first warning is suppressed, and no further ones will be emitted UNLESS... 
      # unless we reset .spaMM.data$options$COMP_maxn_warned <- FALSE after the .tryCatch_W_E()
      # and as .get_inits_by_glm() is called repeatedly for COMPoisson, we do so only if a warning was emitted by the glm
      # which means that COMP_maxn_warned was FALSE in entry to .get_inits_by_glm()
      # which means that COMP_maxn_warned is not set to FALSE .get_inits_by_glm() when it was TRUE in entry to it.
      if (.spaMM.data$options$need_memoise_warning) {
        warning("If the memoise package had been installed (prior to loading spaMM), faster computation could be possible.")
        .spaMM.data$options$need_memoise_warning <- FALSE
      }
    }
  }
  res
}

# debug version including R version replaced by Rcpp version
.debug_COMP_Z_integrand <- function(z,lambda,eta=log(lambda),nu,moment=0,logScaleFac=0) {
  if (missing(lambda)) {
    Cres <- .COMP_Z_integrand(z=z,eta=eta,nu=nu,moment=moment,logScaleFac=logScaleFac)
  } else Cres <- .COMP_Z_integrand(z,lambda,eta,nu,moment,logScaleFac)
  logint <- moment*log(z)+ z*eta - nu*lfactorial(z)
  res <- exp(logint-logScaleFac)
  res <- pmin(.Machine$double.xmax, res)
  res
}

## N O T the moments of the proba distr, Rather the 'num' for a moment obtained by .COMP_Z_ratio(num,denum)
# we call it repeatedly for nu=1 (and variable lambda) but this may be needed to get the continuity correction in nu=1
.COMP_Z_moment <- function(eta,nu,lambda=exp(eta), moment, 
                           pow_lam_nu = exp(eta/nu),
                           maxn=.COMP_maxn(lambda,nu),
                           use_asympto=.CMP_use_asympto(pow_lam_nu,nu,lambda)) {
  if (use_asympto) { 
    if (moment) {
      stop("Execution should not reach this point.")  
    } else { ## moment=0
      if (FALSE) {
        logScaleFac <- nu*pow_lam_nu[[1L]] ## drop any name ! otherwise return value has wrong names
        # using Gaunt et al.:
        c1 <- (nu^2-1)/24
        c2 <- (nu^2-1)*(nu^2+23)/1152
        # c3 <- (nu^2-1)*(5*nu^4-298*nu^2+11237)/414720
        # c4 <- (nu^2-1)*(5*nu^6-1887*nu^4-241041*nu^2+2482411)/39813120
        # c5 <- (nu^2-1)*(7*nu^8-7420*nu^6+1451274*nu^4-220083004*nu^2+1363929895)/6688604160
        # c6 <- (nu^2-1)*(35*nu^10 - 78295*nu^8 + 76299326*nu^6 + 25171388146*nu^4 
        #                 - 915974552561*nu^2 + 4175309343349)/4815794995200
        # c7 <- (nu^2-1)*(5*nu^12- 20190*nu^10 + 45700491*nu^8 - 19956117988*nu^6
        #                 +7134232164555*nu^4-142838662997982*nu^2 + 525035501918789)/115579079884800
        # scaled <- (1+c1/logScaleFac+c2/logScaleFac^2+c3/logScaleFac^3+c4/logScaleFac^4+
        #              c5/logScaleFac^5+c6/logScaleFac^6+c7/logScaleFac^7)
        scaled <- (1+c1/logScaleFac+c2/logScaleFac^2)
        scaled <- (scaled/((2*pi*pow_lam_nu)^((nu-1)/2)*sqrt(nu)))[[1L]] ## drop any name !
        resu <- c(logScaleFac=logScaleFac, scaled=scaled)
      } else resu <- .Rcpp_COMP_Z_asympto(nu, pow_lam_nu)
    }
  } else {
    resu <- .Rcpp_COMP_Z(moment=moment,nu=nu,lambda=lambda,maxn=maxn)
    k_maxim <- ceiling(pow_lam_nu) ## practically locates the mode 
    lfac <- lfactorial(k_maxim)
    if (is.infinite(lfac)) { ## Should not occur as the asymptotic approx for the moments of the *PDF* should have been used 
      stop(paste("Practically infinite sum for COMPoisson's nu=",nu,". The asymptotic approx for the moments of the *PDF* should have been used."))
    } 
    if (maxn>.spaMM.data$options$COMP_maxn-1) {
      logScaleFac <- (k_maxim*eta - nu*lfac)[[1L]] ## drop any name !
      ## Add approximation for tail beyond summation from zero to maxn: 
      ## integrating directly from maxn to Inf may not work, as integrate() may then miss the mode of the integrand. So:
      # (1) Always add tail beyond the highest of maxn and k_maxim
      taildef <- max(k_maxim,maxn)
      if (taildef>1e10) {
        # nintegrate behaves more poorly that what can be handled by max(0,.)
        # E.g., compare the value of 'scaled' in (spaMM:::.COMP_Z_n(lambda=exp(7.321), nu=0.25)) verss (spaMM:::.COMP_Z_n(lambda=exp(7.32), nu=0.25))
        scaled <- .do_call_wrap("quadinf", 
                                arglist=list(f=.COMP_Z_integrand, xa = taildef, xb = Inf, 
                                             eta = eta, nu = nu, moment = 1,logScaleFac = logScaleFac),
                                pack="pracma",
                                info_mess=paste0("If the 'pracma' package were available,\n",
                                                 "more accurate evaluation of an integral would be possible.")
        )$Q
        if (is.null(scaled)) scaled <- max(0,integrate(.COMP_Z_integrand,lower=taildef,upper=Inf,eta=eta,nu=nu,moment=moment,
                                                       logScaleFac=logScaleFac)$value)
      } else {
        scaled <- max(0,integrate(.COMP_Z_integrand,lower=taildef,upper=Inf,eta=eta,nu=nu,moment=moment,
                                  logScaleFac=logScaleFac)$value)
      }
      if (maxn<k_maxim) { ## then add integral from maxn to k_maxim 
        scaled <- scaled + max(0,integrate(.COMP_Z_integrand,lower=maxn,upper=k_maxim,eta=eta,nu=nu,moment=moment,
                                           logScaleFac=logScaleFac)$value)
      }
      scaled <- pmin(.Machine$double.xmax,scaled)
      resu <- .COMP_Z_sum(resu, c(logScaleFac=logScaleFac,scaled=scaled))
    }
  }
  return(resu)
}

##  the moments of the proba distr

.COMP_Pr_moments <- function(pow_lam_nu=lambda^(1/nu), lambda, nu, moments,
                             use_asympto=.CMP_use_asympto(pow_lam_nu,nu, lambda)) {
  resu <- rep(NA,length(moments))
  names(resu) <- moments
  if (use_asympto) {
    # using Gaunt et al.:
    c1 <- (nu^2-1)/24
    mu1 <- pow_lam_nu*( 1 -(nu-1)/(2*nu*pow_lam_nu) - c1/((nu*pow_lam_nu)^2)- c1/((nu*pow_lam_nu)^3))
    if ("1" %in% moments) resu["1"] <- mu1
    if ("2" %in% moments || "3" %in% moments) {
      Var <- (pow_lam_nu/nu)*( 1 + c1/((nu*pow_lam_nu)^2) + 2*c1/((nu*pow_lam_nu)^3)) # >0 when nu>1
      if ("2" %in% moments) resu["2"] <- mu1^2 + Var # supp to mu1^2 when nu>1
    }
    if ("3" %in% moments) {
      k3 <- (pow_lam_nu/nu^2)*( 1 - c1/((nu*pow_lam_nu)^2) -4*c1/((nu*pow_lam_nu)^3))
      resu["3"] <- mu1^3 + 3*mu1*Var +k3
    }
  } else {
    denum <- .COMP_Z(lambda=lambda,nu=nu) # -> calling .COMP_Pr_moments( use_asympto=TRUE) ... 
    denum_corr <- .COMP_Z(lambda=lambda,nu=1) # for continuity correction in nu=1
    if ("1" %in% moments) {
      num <- .COMP_Z_n(lambda=lambda,nu=nu)
      resu["1"] <- .COMP_Z_ratio(num,denum)
      # continuity correction wrt poisson: corr =lambda + error of the COMP_Z... functions
      corr <- .COMP_Z_ratio(.COMP_Z_n(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
      resu["1"] <- resu["1"]+(lambda-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
    }
    if ("2" %in% moments) {
      num <- .COMP_Z_n2(lambda=lambda,nu=nu)
      resu["2"] <- .COMP_Z_ratio(num,denum)
      # cotinuity correction wrt poisson: 
      corr <- .COMP_Z_ratio(.COMP_Z_n2(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
      resu["2"] <- resu["2"]+(lambda*(1+lambda)-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
    }
    if ("3" %in% moments) {
      num <- .COMP_Z_n3(lambda=lambda,nu=nu)
      resu["3"] <- .COMP_Z_ratio(num,denum)
      # cotinuity correction wrt poisson: 
      corr <- .COMP_Z_ratio(.COMP_Z_n3(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
      resu["3"] <- resu["3"]+(lambda*(1+lambda*(3+lambda))-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
    }
  }
  return(resu)
}


.COMP_Z <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  if (lambda==0) return(c(logScaleFac=0,scaled=1))  
  ## ELSE
  if (missing(eta)) eta <- log(lambda)
  resu <- .COMP_Z_moment(moment=0,eta=eta,nu=nu,lambda=lambda,maxn=maxn)
  return(resu)
  # # R version of the summation term:
  # if (nu==0) {
  #   return(c(logScaleFac=0,scaled=as.numeric(1/(1-lambda))))
  # } 
  # floorn <- floor(maxn)
  # epsn <- maxn - floorn
  # facs <- c(1,lambda/seq(floorn+1L)^nu)
  # cumprodfacs <- cumprod(facs)
  # cumprodfacs[floorn+2L] <- cumprodfacs[floorn+2L]*epsn
  # scaled <- sum(cumprodfacs)
  # if ( ! is.finite(scaled)) {
  #   logfacs <- log(facs)
  #   cumsumlogfacs <- cumsum(logfacs)
  #   refi <- which.max(cumsumlogfacs)
  #   logScaleFac <- cumsumlogfacs[refi]
  #   scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  # } else logScaleFac <- 0
  # return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z = exp(logScaleFac)*scaled
}

.COMP_Z_n <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  if (missing(eta)) eta <- log(lambda)
  resu <- .COMP_Z_moment(moment=1,eta=eta,nu=nu,lambda=lambda,maxn=maxn, use_asympto=FALSE)
  return(resu)
  # # R version:
  # if (nu==0) {
  #   return(c(logScaleFac=0,scaled=as.numeric(lambda/(1-lambda)^2)))
  # } 
  # floorn <- floor(maxn)
  # epsn <- maxn - floorn
  # seqn <- seq(2L,floorn+1L)
  # facs <- c(lambda,(seqn/(seqn-1L))*lambda/(seqn^nu))
  # cumprodfacs <- cumprod(facs)
  # cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  # scaled <- sum(cumprodfacs)
  # if ( ! is.finite(scaled)) {
  #   logfacs <- log(facs)
  #   cumsumlogfacs <- cumsum(logfacs)
  #   refi <- which.max(cumsumlogfacs)
  #   logScaleFac <- cumsumlogfacs[refi]
  #   scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  # } else logScaleFac <- 0
  # return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n = exp(logScaleFac)*scaled
}

.COMP_Z_n2 <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  if (missing(eta)) eta <- log(lambda)
  resu <- .COMP_Z_moment(moment=2,eta=eta,nu=nu,lambda=lambda,maxn=maxn, use_asympto=FALSE)
  return(resu)
  # # R version:
  # if (nu==0) {
  #   return(c(logScaleFac=0,scaled=as.numeric(lambda*(1+lambda)/(1-lambda)^3)))
  # } 
  # floorn <- floor(maxn)
  # epsn <- maxn - floorn
  # seqn <- seq(2L,floorn+1L)
  # facs <- c(lambda,((seqn/(seqn-1L))^2)*lambda/(seqn^nu))
  # cumprodfacs <- cumprod(facs)
  # cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  # scaled <- sum(cumprodfacs)
  # if ( ! is.finite(scaled)) {
  #   logfacs <- log(facs)
  #   cumsumlogfacs <- cumsum(logfacs)
  #   refi <- which.max(cumsumlogfacs)
  #   logScaleFac <- cumsumlogfacs[refi]
  #   scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  # } else logScaleFac <- 0
  # return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n2 = exp(logScaleFac)*scaled
}

.COMP_Z_n3 <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  if (missing(eta)) eta <- log(lambda)
  resu <- .COMP_Z_moment(moment=3,eta=eta,nu=nu,lambda=lambda,maxn=maxn, use_asympto=FALSE)
  return(resu)
  # # R version:
  # if (nu==0) {
  #   return(c(logScaleFac=0,scaled=as.numeric(lambda*(1+4*lambda+lambda^2)/(1-lambda)^4)))
  # } 
  # floorn <- floor(maxn)
  # epsn <- maxn - floorn
  # seqn <- seq(2L,floorn+1L)
  # facs <- c(lambda,((seqn/(seqn-1L))^3)*lambda/(seqn^nu))
  # cumprodfacs <- cumprod(facs)
  # cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  # scaled <- sum(cumprodfacs)
  # if ( ! is.finite(scaled)) {
  #   logfacs <- log(facs)
  #   cumsumlogfacs <- cumsum(logfacs)
  #   refi <- which.max(cumsumlogfacs)
  #   logScaleFac <- cumsumlogfacs[refi]
  #   scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  # } else logScaleFac <- 0
  # return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n3 = exp(logScaleFac)*scaled
}

.COMP_Z_lfacn <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,lambda/(seqn^nu)) ## n from 1 to floorn+1
  cumprodfacs <- cumprod(facs) ## sequence of lambda^n / n!^nu 
  cumsumlogn <- cumsum(c(0,log(seqn))) ## sequence of log(n!) for n from 1 to floorn+1
  cumprodfacs <- cumprodfacs*cumsumlogn ## sequence of log(n!) lambda^n / n!^nu 
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n = exp(logScaleFac)*scaled
}


.COMP_Z_ratio <- function(Z1,Z2,log=FALSE) {
  logratio <- Z1[["logScaleFac"]]-Z2[["logScaleFac"]]+log(Z1[["scaled"]]/Z2[["scaled"]]) # more overflow-robust than directly using exp()
  if (log) {
    return(logratio)
  } else return(exp(logratio))
} 

.COMP_Z_sum <- function(Z1,Z2) {
  z12sc <- Z1[["logScaleFac"]] - Z2[["logScaleFac"]]
  #keep largest logScaleFac
  if (z12sc>0) {
    logScaleFac <- Z1[["logScaleFac"]]
    Z1sc <- Z1[["scaled"]]
    Z2sc <- Z2[["scaled"]] / exp(z12sc)
  } else {
    logScaleFac <- Z2[["logScaleFac"]]
    Z1sc <- Z1[["scaled"]] * exp(z12sc)
    Z2sc <- Z2[["scaled"]]
  }
  scaled <- pmin(.Machine$double.xmax, Z1sc+Z2sc)
  return(c(logScaleFac=logScaleFac,scaled=scaled))
} 

.COMP_dnu_objfn <- function(nu,y,eta,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)) { ## not currently used
  summand <- numeric(length(y))
  for (i in seq_len(length(y))) {
    yi <- y[i]
    lambdai <- lambda[i]
    comp_z_lfacn <- .COMP_Z_lfacn(nu=nu,lambda=lambdai,maxn=maxn)
    comp_z <- .COMP_Z(nu=nu,lambda=lambdai,maxn=maxn)
    summand[i] <- lfactorial(yi)-.COMP_Z_ratio(comp_z_lfacn,comp_z,log=TRUE)
  }
  return(sum(summand))
}

.dCOMP <- function(x, mu, family_env,
                   nu=family_env$nu,
                   lambda= family_env$mu2lambda(mu), # COMPoisson(nu=nu)$linkfun(mu,log=FALSE),
                   log = FALSE, maxn=.COMP_maxn(lambda,nu)) {
  compz <- .COMP_Z(lambda=lambda,nu=nu,maxn=maxn)
  logd <- x * log(lambda) - nu* lfactorial(x) - compz[["logScaleFac"]] -log(compz[["scaled"]])
  if (log) { 
    return(logd)
  } else return(exp(logd))
}

.COMP_simulate <- function(lambda,nu,nsim=1L) { ## for scalar lambda
  maxn <- .COMP_maxn(lambda=lambda,nu=nu)
  Z <- .COMP_Z(lambda=lambda,nu=nu)
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,(seqn/(seqn-1L))*lambda/(seqn^nu))
  cumprodfacs <- cumprod(facs)
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  sample(floorn+1L,size=nsim,prob=cumprodfacs/sum(cumprodfacs), replace=TRUE)
}

.CMP_lambda2mui <- function(lambda,nu) { # mu(lambda)
  if (lambda==0) {
    return(1e-8)  
  } else {
    mu <- .COMP_Pr_moments(lambda=lambda, nu=nu, moments="1")
    return(mu) #(max(mu,1e-8)) would avoid mu=0 which fails validmu(mu) (as for poisson family)
    # but it's better to use a consistent code to sanitize the input lambda than to fix the return value.
  }
}

.CMP_lambda2mu <- function(lambda) { # mu(lambda)  
  mus <- attr(lambda,"mu") 
  if (is.null(mus)) {
    if (nu==0) {
      mus <- lambda/(1-lambda) 
    } else {
      nu <- parent.env(environment())$nu
      # mus <- sapply(lambda, .CMP_lambda2mui,nu=nu) 
      mus <- numeric(length(lambda))
      for (it in seq_along(lambda)) mus[it] <- .CMP_lambda2mui(lambda[it],nu=nu)
    }
    dim(mus) <- dim(lambda) ## may be NULL
  }
  attributes(lambda) <- NULL
  return(structure(mus,lambda=lambda))
}

# different input
.CMP_loglambda_linkinv <- function(eta,lambda=exp(eta)) { # \equiv .CMP_lambda2mu(lambda) but with nu properly passed...
  mus <- attr(lambda,"mu") 
  if (is.null(mus)) {
    if (nu==0) {
      mus <- lambda/(1-lambda) 
    } else {
      nu <- parent.env(environment())$nu # for canonical link .CMP_loglambda_linkinv has been copied to $linkinv with nu in its parent env.
      #mus <- sapply(lambda, .CMP_lambda2mui,nu=nu)
      mus <- numeric(length(lambda))
      for (it in seq_along(lambda)) mus[it] <- .CMP_lambda2mui(lambda[it],nu=nu)
    }
    dim(mus) <- dim(lambda) ## may be NULL
  }
  attributes(lambda) <- NULL
  return(structure(mus,lambda=lambda))
}



..CMP_mu2lambda <- function(mu,nu, CMP_linkfun_objfn) { 
  if (nu==1) {
    lambda <- mu ## pb du code general est qu'objfn n'est alors que l'erreur numérique de linkinv()
  } else if (mu==Inf) {
    warning(paste("Approximating lambda as 1-1e-8 in COMPoisson(nu=",nu,") for mu = Inf"))
    lambda <- 1 - 1e-8 ## 1 -> mu.eta = NaN. mu.eta being Inf would not be OK bc Inf GLMweights -> 1/w.resid=0 -> singular Sig or d2hdv2
  } else {
    app_lambda <-  max(c(0,(mu+max(0,(nu-1)/(2*nu)))))^(nu)
    lambdamin <- max(c(.Machine$double.eps, app_lambda-1)) ## (the ultimate constraint is that log(lambda)> -Inf)
    lambdamax <- max(c(.Machine$double.eps, (app_lambda+1)*c(1.01))) ## not perfect...
    # last one for low lambda,nu values
    fupper <- CMP_linkfun_objfn(lambdamax, mu=mu)
    while (fupper<0) {
      lambdamax <- lambdamax*2
      fupper <- CMP_linkfun_objfn(lambdamax, mu=mu)
      if (is.nan(fupper)) break
    }
      ## normally >0 
    if (is.nan(fupper)) {
      if ( ! identical(spaMM.getOption("COMP_geom_approx_warned"),TRUE)) {
        warning(paste("Geometric approximation tried in COMPoisson(nu=",nu,") for mu=",mu," and possibly for other nu,mu values"))
        .spaMM.data$options$COMP_geom_approx_warned <- TRUE
      }
      lambda <- mu/(1+mu) ## but it should be checked hence the warning. Ideally for large nu ?
    } else {
      flower <- CMP_linkfun_objfn(lambdamin, mu=mu) ## normally <0 
      while ((flower)>0) {
        lambdamin <- lambdamin/2
        flower <- CMP_linkfun_objfn(lambdamin, mu=mu)
        if (is.nan(flower)) break
      }
      interval <- c(lambdamin,lambdamax)
      # linkinv(0)=0 => objfn(0)= -mu
      lambda <- uniroot(CMP_linkfun_objfn,interval=interval,f.lower = flower,f.upper = fupper, mu=mu)$root
    }
  }
  return(lambda)
}

# different arguments
.CMP_mu2lambda <- function(mu) {## scalar or vector
  lambdas <- attr(mu,"lambda")
  if ( is.null(lambdas)) {
    nu <- parent.env(environment())$nu
    if (nu==0) {
      lambdas <- mu/(1+mu) 
    } else {
      #lambdas <- sapply(mu, ..CMP_mu2lambda,nu=nu, CMP_linkfun_objfn=parent.env(environment())$CMP_linkfun_objfn)
      lambdas <- numeric(length(mu))
      for (it in seq_along(mu)) lambdas[it] <- ..CMP_mu2lambda(mu[it],nu=nu, CMP_linkfun_objfn=parent.env(environment())$CMP_linkfun_objfn)
    }
  } 
  attributes(mu) <- NULL ## avoids 'mise en abime'
  return(structure(lambdas,mu=mu))
}

# different return value
.CMP_loglambda_linkfun <- function(mu) {## scalar or vector
  lambdas <- attr(mu,"lambda")
  if ( is.null(lambdas)) {
    nu <- parent.env(environment())$nu
    if (nu==0) {
      lambdas <- mu/(1+mu) 
    } else {
      #lambdas <- sapply(mu, ..CMP_mu2lambda,nu=nu, CMP_linkfun_objfn=parent.env(environment())$CMP_linkfun_objfn)
      lambdas <- numeric(length(mu))
      for (it in seq_along(mu)) lambdas[it] <- ..CMP_mu2lambda(mu[it],nu=nu, CMP_linkfun_objfn=parent.env(environment())$CMP_linkfun_objfn)
    }
  } 
  attributes(mu) <- NULL ## avoids 'mise en abime'
  return(structure(log(lambdas),mu=mu)) ## eta (ie standard linkfun value) for the Shmueli version
}

.CMP_dev.resid <- function(yi,lambdai,nu, CMP_linkfun_objfn) {
  Z2 <- .COMP_Z(lambda=lambdai,nu=nu)
  if (yi==0) { # lambda = 0,  Z1 = 1
    dev <- 2*(Z2[["logScaleFac"]]+log(Z2[["scaled"]]))
  } else {
    # eval lambda for Z1(lambda...)
    if (nu==0) {
      lambda <- yi/(1+yi)  
    } else {
      app_lambda <-  max(c(0,(yi+max(0,(nu-1)/(2*nu)))))^(nu)
      lambdamin <- max(c(.Machine$double.eps, app_lambda-1))
      while (CMP_linkfun_objfn(lambdamin, mu=yi)>0) {
        lambdamin <- lambdamin/2 ## greedy 
        if (lambdamin<1e-8) stop("lambdamin<1e-8 in COMPoisson$dev.resids()")
      }
      lambdamax <- max(c(.Machine$double.eps, (app_lambda+1)*c(1.01))) 
      while (CMP_linkfun_objfn(lambdamax, mu=yi)<0) {
        lambdamax <- lambdamax*10 ## greedy but resolves errors for low nu 
        if (lambdamax>1e100) stop("lambdamax>1e100 in COMPoisson$dev.resids()")
      }
      interval <- c(lambdamin,lambdamax)
      lambda <- uniroot(CMP_linkfun_objfn,interval=interval, mu=yi)$root
    }
    Z1 <- .COMP_Z(lambda=lambda,nu=nu)
    #
    dev <- 2*(yi*log(lambda/lambdai)-(Z1[["logScaleFac"]]-Z2[["logScaleFac"]]+log(Z1[["scaled"]]/Z2[["scaled"]])))
  }
  dev
}

.CMP_dev_resids <- function(y, mu, wt){
  # must accept, among others, vector y and scalar mu.
  lambdas <- parent.env(environment())$mu2lambda(mu=mu) # parent envir of COMPoisson()$dev.resids
  n <- length(y)
  if (length(mu)==1L) lambdas <- rep(lambdas,n)
  devs <- numeric(n)
  nu <- parent.env(environment())$nu
  CMP_linkfun_objfn <- parent.env(environment())$CMP_linkfun_objfn
  for(i in seq(n)) { devs[i] <- .CMP_dev.resid(y[i],lambdai=lambdas[i],nu=nu, CMP_linkfun_objfn=CMP_linkfun_objfn) }
  devs <- devs*wt
  devs[devs==Inf] <- .Machine$double.xmax/n ## so that total deviance may be finite 
  return(devs) 
}

.CMP_aic <- function(y, n, mu, wt, dev) {
  family_env <- parent.env(environment())
  aici <- numeric(length(y))
  for (i in seq_len(length(y))) { aici[i] <- .dCOMP(y[i], mu[i], family_env=family_env, log = TRUE) }
  -2 * sum(aici * wt)
}

.CMP_simfun <- function(object,nsim) {
  wts <- object$prior.weights
  if (any(wts != 1)) 
    warning("ignoring prior weights")
  if (inherits(object,"glm")) { # if COMPoisson() used in a glm() call simulate.lm -> family$simulate -> here
    lambdas <- attr(object$fitted.values,"lambda")
    if (is.null(lambdas)) lambdas <- sapply(mu, family$mu2lambda)
  } else { # HLfit object: 
    mu <- object$muetablob$mu
    lambdas <- attr(mu,"lambda") # exp(object$eta)
    if (is.null(lambdas)) lambdas <- sapply(mu, family$mu2lambda)
  }
  nu <- parent.env(environment())$nu
  resu <- sapply(lambdas,.COMP_simulate,nu=nu,nsim=nsim) # vector if nsim=1, nsim-row matrix otherwise
  if (nsim>1) resu <- t(resu)
  return(resu)  # vector if nsim=1, nsim-col matrix otherwise
}

.CMP_mu.eta <- function (eta,lambda=exp(eta)) {
  nu <- parent.env(environment())$nu ## not necess for R to find it, but to avoid a NOTE from CHECK
  if (nu==0) {
    resu <- lambda/(1-lambda)^2
  } else {
    resu <- numeric(length(lambda))
    for (it in seq_len(length(lambda))) {
      lambdai <- lambda[it]
      if (lambda[it]==0) {
        resu[it] <- .Machine$double.eps
      } else {
        moments <- .COMP_Pr_moments(lambda=lambda[it], nu=nu, moments=c("1","2"))
        res <- moments["2"] - moments["1"]^2
        resu[it] <- pmax(res, .Machine$double.eps)
      }
    }
  }
  resu
}

.CMP_variance <- function(mu) {
  nu <- parent.env(environment())$nu
  if (nu==0) {
    return(mu*(1+mu)) 
  } else {
    lambdas <- parent.env(environment())$mu2lambda(mu)
    resu <- numeric(length(lambdas))
    for (it in seq_len(length(lambdas))) {
      if (lambdas[it]==0) {
        resu[it] <- 1e-8
      } else {
        moments <- .COMP_Pr_moments(lambda=lambdas[it],nu=nu,moments=c("1","2"))
        resu[it] <- max(moments["2"] -moments["1"]^2,1e-8) ## pmax otherwise for low mu, Vmu=0, -> ... -> w.resid=0
      }
    }
    ## for geom V(mu) = mu(1+mu) but this cannot be a lower bound for resu.
    return(resu)
  }
}

COMPoisson <- function(nu = stop("COMPoisson's 'nu' must be specified"), 
                       link = "loglambda" # eta <-> mu link, not the eta <-> lambda log link
) { 
  .spaMM.data$options$COMP_maxn_warned <- FALSE # much better here than in .preprocess(); works with glm()
  .spaMM.data$options$COMP_geom_approx_warned <- FALSE
  if (inherits(nuch <- substitute(nu),"character") ||
      (inherits(nuch,"name") && inherits(nu, "function")) # "name" is for e.g. COMPoisson(log)
      # (but testing only "name" would catch e.g. COMPoisson(nu=nu) )
     ) { 
    if (inherits(nuch,"character")) nuch <- paste0('"',nuch,'"')
    errmess <- paste0('It looks like COMPoisson(',nuch,') was called, which absurdly means COMPoisson(nu=',nuch,
                      ').\n  Use named argument: COMPoisson(link=',nuch,') instead.')
    stop(errmess)
  }
  # When 'nu' is recognized as as call to some function ! = stop(), we eval it so it is no longer recognized as a call by .calc_optim_args()
  if (inherits(nuch,"call") && deparse(nuch[[1]])!="stop") nu <- eval(nuch) 
  
  linktemp <- substitute(link) # if link was char LHS is char ; else deparse will create a char from a language object 
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("loglambda", "log", "identity", "sqrt")
  if (inherits(link, "link-glm")) { # a make.link() object was provided
    stats <- link
    if ( ! is.null(stats$name)) linktemp <- stats$name
  } else if (linktemp=="loglambda") {
    linkfun <- .CMP_loglambda_linkfun # need to associate nu to this one (not to more classical linkfun's)
    linkinv <- .CMP_loglambda_linkinv
    mu.eta <- .CMP_mu.eta
    environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment()  ## containing nu
  } else if (linktemp %in% okLinks) {
    stats <- make.link(linktemp)
    linkfun <- stats$linkfun
    linkinv <- stats$linkinv
    mu.eta <- stats$mu.eta
  } # at this point 'stats' is always the result of make.link and 'linktemp' is always char. 
  if ( ! linktemp %in% okLinks) stop(gettextf("link \"%s\" not available for COMPoisson family; available links are %s", 
                                              linktemp, paste(sQuote(okLinks), collapse = ", ")),           domain = NA)
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0) ## from poisson()
  valideta <- function(eta) TRUE ## from poisson()
  mu2lambda <- .CMP_mu2lambda # copy to pass nu, used in in .CMP_dev_resids2()
  lambda2mu <- .CMP_lambda2mu # ## copied locally to pass nu to .CMP_lambda2mui() through the lambda2mu parent envir
  CMP_linkfun_objfn <- function(lambda, mu) {lambda2mu(lambda=lambda) -mu} # objective fn of ~.CMP_mu2lambda # called by .CMP_dev.resid() to find lambda=lambda(mu)
  variance <- .CMP_variance
  dev.resids <- .CMP_dev_resids
  aic <- .CMP_aic
  initialize <- expression({
    if (any(y < 0L)) stop("negative values not allowed for the 'COMPoisson' family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  simfun <- .CMP_simfun 
  environment(dev.resids) <- environment(aic) <- environment(simfun) <- environment(mu2lambda) <- 
    environment(variance) <- environment(lambda2mu) <- environment() ## containing nu
  ## Change the parent.env of all functions that have the same envir as aic(): 
  parent.env(environment(aic)) <- environment(.dCOMP) ## gives access to spaMM:::.dCOMP and other .COMP_ fns
  structure(list(family = structure("COMPoisson",
                                    withArgs=quote(paste0("COMPoisson(nu=",signif(nu,4),")"))), 
                 link = linktemp, linkfun = linkfun,           mu2lambda=mu2lambda,
                 linkinv = linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = valideta, simulate = simfun), 
            class = "family")
}
