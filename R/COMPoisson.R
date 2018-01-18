.CMP_use_asympto <- function(pow_lam_nu, nu) {
  ## Gaunt et al suggest lambda>1.5 and pow_lam_nu>1.5
  ## this ensures that nu*pow_lam_nu>1.1 but the expansion is clearly better for nu*pow_lam_nu >> 1 
  return(eval(.spaMM.data$options$CMP_asympto_cond))
}

.COMP_maxn <- function(lambda,nu) { ## this uses the simplest asympto approx to determine the discrete sum bound 
  pow_lam_nu <- lambda^(1/nu)
  # if (.CMP_use_asympto(pow_lam_nu,nu)) {
  #   res <- stop(".COMP_maxn() called despite .CMP_use_asympto(pow_lam_nu,nu) being TRUE") 
  # } else { 
    app_En <- lambda^(1/(nu+1e-6)) ## .COMP_asympto_P_moment(pow_lam_nu,nu+1e-6,1L) would be lower
    # and using var ~ En/nu
    res <- max(2,1+app_En+6*sqrt(app_En/(nu+1e-6)))  ## real, to allow continuity correction
    opt_maxn <- spaMM.getOption("COMP_maxn")
    if (res>opt_maxn) {
      res <- opt_maxn
      if ( ! identical(spaMM.getOption("COMP_maxn_warned"),TRUE)) {
        warning(paste("maxn truncated to",res,"for (lambda,nu)=",lambda,nu,"and possibly other values."))
        spaMM.options(COMP_maxn_warned=TRUE)
      }
    }
  # }
  res
}

.COMP_Z_integrand <- function(z,lambda,eta=log(lambda),nu,moment=0,logScaleFac=0) {
  logint <- moment*log(z)+ z*eta - nu*lfactorial(z)
  res <- exp(logint-logScaleFac)
  res <- pmin(.Machine$double.xmax, res)
  res
}

## N O T the moments of the proba distr
.COMP_Z_moment <- function(eta,nu,lambda=exp(eta), moment, 
                           pow_lam_nu = exp(eta/nu),
                           maxn=.COMP_maxn(lambda,nu),
                           use_asympto=.CMP_use_asympto(pow_lam_nu,nu)) {
  if (use_asympto) { 
    if (moment) {
      stop("Execution should not reach this point.")  
    } else { ## moment=0
      logScaleFac <- nu*pow_lam_nu[[1L]] ## drop any name ! otherwise return value has wrong names
      # using Gaunt et al.:
      c1 <- (nu^2-1)/24
      c2 <- (nu^2-1)*(nu^2+23)/1152
      c3 <- (nu^2-1)*(5*nu^4-298*nu^2+11237)/414720
      c4 <- (nu^2-1)*(5*nu^6-1887*nu^4-241041*nu^2+2482411)/39813120
      c5 <- (nu^2-1)*(7*nu^8-7420*nu^6+1451274*nu^4-220083004*nu^2+1363929895)/6688604160
      c6 <- (nu^2-1)*(35*nu^10 - 78295*nu^8 + 76299326*nu^6 + 25171388146*nu^4 
                      - 915974552561*nu^2 + 4175309343349)/4815794995200
      c7 <- (nu^2-1)*(5*nu^12- 20190*nu^10 + 45700491*nu^8 - 19956117988*nu^6
                      +7134232164555*nu^4-142838662997982*nu^2 + 525035501918789)/115579079884800
      scaled <- (1+c1/logScaleFac+c2/logScaleFac^2+c3/logScaleFac^3+c4/logScaleFac^4+c5/logScaleFac^5+c6/logScaleFac^6+c7/logScaleFac^7)
      scaled <- (scaled/((2*pi*pow_lam_nu)^((nu-1)/2)*sqrt(nu)))[[1L]] ## drop any name !
      resu <- c(logScaleFac=logScaleFac, scaled=scaled)
    }
  } else {
    resu <- .Rcpp_COMP_Z(moment=moment,nu=nu,lambda=lambda,maxn=maxn)
    k_maxim <- ceiling(pow_lam_nu)
    lfac <- lfactorial(k_maxim)
    if (is.infinite(lfac)) { ## Should not occur as the asymptotic approx for the moments of the *PDF* should have been used 
      stop(paste("Practically infinite sum for COMPoisson's nu=",nu,". The asymptotic approx for the moments of the *PDF* should have been used."))
    } 
    logScaleFac <- (k_maxim*eta - nu*lfac)[[1L]] ## drop any name !
    ## need to pinpoint the "maximum" otherwise integrate misses it. 
    scaled <- max(0,integrate(.COMP_Z_integrand,lower=k_maxim,upper=Inf,eta=eta,nu=nu,moment=moment,
                         logScaleFac=logScaleFac,stop.on.error = FALSE)$value)
    if (maxn<k_maxim) {
      scaled <- scaled + max(0,integrate(.COMP_Z_integrand,lower=maxn,upper=k_maxim,eta=eta,nu=nu,moment=moment,
                                         logScaleFac=logScaleFac,stop.on.error = FALSE)$value)
    }
    scaled <- pmin(.Machine$double.xmax,scaled)
    resu <- .COMP_Z_sum(resu, c(logScaleFac=logScaleFac,scaled=scaled))
  }
  return(resu)
}

##  the moments of the proba distr

.COMP_Pr_moments <- function(pow_lam_nu=lambda^(1/nu), lambda, nu, moments,use_asympto=.CMP_use_asympto(pow_lam_nu,nu)) {
  resu <- rep(NA,length(moments))
  names(resu) <- moments
  if (use_asympto) {
    # using Gaunt et al.:
    c1 <- (nu^2-1)/24
    mu1 <- pow_lam_nu*( 1 -(nu-1)/(2*nu*pow_lam_nu) - c1/((nu*pow_lam_nu)^2)- c1/((nu*pow_lam_nu)^3))
    if ("1" %in% moments) resu["1"] <- mu1
    if ("2" %in% moments || "3" %in% moments) {
      Var <- (pow_lam_nu/nu)*( 1 + c1/((nu*pow_lam_nu)^2) + 2*c1/((nu*pow_lam_nu)^3))
      if ("2" %in% moments) resu["2"] <- mu1^2 + Var
    }
    if ("3" %in% moments) {
      k3 <- (pow_lam_nu/nu^2)*( 1 - c1/((nu*pow_lam_nu)^2) -4*c1/((nu*pow_lam_nu)^3))
      resu["3"] <- mu1^3 + 3*mu1*Var +k3
    }
  } else {
    denum <- .COMP_Z(lambda=lambda,nu=nu)
    if ("1" %in% moments) {
      num <- .COMP_Z_n(lambda=lambda,nu=nu)
      resu["1"] <- .COMP_Z_ratio(num,denum)
    }
    if ("2" %in% moments) {
      num <- .COMP_Z_n2(lambda=lambda,nu=nu)
      resu["2"] <- .COMP_Z_ratio(num,denum)
    }
    if ("3" %in% moments) {
      num <- .COMP_Z_n3(lambda=lambda,nu=nu)
      resu["3"] <- .COMP_Z_ratio(num,denum)
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

.dCOMP <- function(x, mu, nu,
                  lambda=COMPoisson(nu=nu)$linkfun(mu,log=FALSE),
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
  sample(floorn+1L,size=nsim,prob=cumprodfacs/sum(cumprodfacs))
}

.CMP_linkinvi <- function(lambda,nu) {
  if (lambda==0) {
    return(1e-8)  
  } else {
    mu <- .COMP_Pr_moments(lambda=lambda, nu=nu, moments="1")
    return(max(mu,1e-8)) ## avoids mu=0 which fails validmu(mu) (as for poisson family)
  }
}

.CMP_linkfuni <- function(mu,nu, CMP_linkfun_objfn) {
  if (nu==1) {
    lambda <- mu ## pb du code general est qu'objfn n'est alors que l'erreur numérique de linkinv()
  } else if (mu==Inf) {
    warning(paste("Approximating lambda as 1-1e-8 in COMPoisson(nu=",nu,") for mu = Inf"))
    lambda <- 1 - 1e-8 ## 1 -> mu.eta = NaN. mu.eta being Inf would not be OK bc Inf GLMweights -> 1/w.resid=0 -> singular Sig or d2hdv2
  } else {
    app_lambda <-  max(c(0,(mu+max(0,(nu-1)/(2*nu)))))^(nu)
    lambdamin <- max(c(.Machine$double.eps, app_lambda-1)) ## (the ultimate contraint is that log(lambda)> -Inf)
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
        spaMM.options(COMP_geom_approx_warned=TRUE)
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
  lambdas <- parent.env(environment())$linkfun(mu=mu,log=FALSE)
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

.CMP_linkinv <- function(eta,lambda=exp(eta)) {
  mus <- attr(lambda,"mu") 
  if (is.null(mus)) {
    if (nu==0) {
      mus <- lambda/(1-lambda) 
    } else {
      nu <- parent.env(environment())$nu
      mus <- sapply(lambda, .CMP_linkinvi,nu=nu)
    }
    dim(mus) <- dim(lambda) ## may be NULL
  }
  attributes(lambda) <- NULL
  return(structure(mus,lambda=lambda))
}

.CMP_aic <- function(y, n, mu, wt, dev) {
  nu <- parent.env(environment())$nu
  aici <- numeric(length(y))
  for (i in seq_len(length(y))) { aici[i] <- .dCOMP(y[i], mu[i], nu=nu, log = TRUE) }
  -2 * sum(aici * wt)
}

.CMP_simfun <- function(object,nsim) {
  wts <- object$prior.weights
  if (any(wts != 1)) 
    warning("ignoring prior weights")
  lambdas <- exp(object$eta)
  nu <- parent.env(environment())$nu
  resu <- sapply(lambdas,.COMP_simulate,nu=nu,nsim=nsim)
  if (nsim>1) resu <- t(resu)
  return(resu)
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
    lambdas <- parent.env(environment())$linkfun(mu,log=FALSE)
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

.CMP_linkfun <- function(mu, ## scalar or vector
                         log=TRUE) { ## log=TRUE => returns eta; else returns lambda, with mu attribute
  lambdas <- attr(mu,"lambda")
  if ( is.null(lambdas)) {
    nu <- parent.env(environment())$nu
    if (nu==0) {
      lambdas <- mu/(1+mu) 
    } else lambdas <- sapply(mu, .CMP_linkfuni,nu=nu, CMP_linkfun_objfn=parent.env(environment())$CMP_linkfun_objfn)
  } 
  attributes(mu) <- NULL ## avoids 'mise en abime'
  if (log) {
    return(structure(log(lambdas),mu=mu)) ## eta, ie standard linkfun value
  } else {
    return(structure(lambdas,mu=mu))
  }
}

COMPoisson <- function(nu = stop("COMPoisson's 'nu' must be specified"), 
                        link = "loglambda" # eta <-> mu link, not the eta <-> lambda log link
                        ) {
  mc <- match.call()
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  okLinks <- c("loglambda")
  if (linktemp %in% okLinks) {} else {
      stop(gettextf("link \"%s\" not available for COMPoisson family; available links are %s", 
                    linktemp, paste(sQuote(okLinks), collapse = ", ")), 
           domain = NA)
  }
  linkinv <- .CMP_linkinv
  CMP_linkfun_objfn <- function(lambda, mu) {linkinv(lambda=lambda) -mu} ## def'd locally to handle the call to linkinv()
  linkfun <- .CMP_linkfun
  variance <- .CMP_variance
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0) ## from poisson()
  valideta <- function(eta) TRUE ## from poisson()
  mu.eta <- .CMP_mu.eta
  dev.resids <- .CMP_dev_resids
  aic <- .CMP_aic
  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed for the 'Poisson' family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  simfun <- .CMP_simfun
  environment(dev.resids) <- environment(linkinv) <- environment(aic) <- environment(simfun) <- 
    environment(mu.eta) <- environment(variance) <- environment(linkfun) <- environment() ## containing nu
  ## changes the parent.env of all functions: 
  parent.env(environment(aic)) <- environment(.dCOMP) ## gives access to spaMM:::.dCOMP and other .COMP_ fns
  structure(list(family = structure("COMPoisson",
                                    withArgs=quote(paste("COMPoisson(nu=",signif(nu,4),")",sep=""))), 
                 link = linktemp, linkfun = linkfun, 
                 linkinv = linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = valideta, 
                 simulate = simfun), 
            class = "family")
}

