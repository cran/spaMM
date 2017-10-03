.COMP_maxn <- function(lambda,nu) {
  app_En <- lambda^(1/(nu+1e-6))
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
  res
}

.COMP_Z <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  return(.Rcpp_COMP_Z(moment=0,nu=nu,lambda=lambda,maxn=maxn))
  # R version:
  if (nu==0) {
    return(c(logScaleFac=0,scaled=as.numeric(1/(1-lambda))))
  } 
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  facs <- c(1,lambda/seq(floorn+1L)^nu)
  cumprodfacs <- cumprod(facs)
  cumprodfacs[floorn+2L] <- cumprodfacs[floorn+2L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z = exp(logScaleFac)*scaled
}

.COMP_Z_n <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  return(.Rcpp_COMP_Z(moment=1,nu=nu,lambda=lambda,maxn=maxn))
  # R version:
  if (nu==0) {
    return(c(logScaleFac=0,scaled=as.numeric(lambda/(1-lambda)^2)))
  } 
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,(seqn/(seqn-1L))*lambda/(seqn^nu))
  cumprodfacs <- cumprod(facs)
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

.COMP_Z_n2 <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  return(.Rcpp_COMP_Z(moment=2,nu=nu,lambda=lambda,maxn=maxn))
  # R version:
  if (nu==0) {
    return(c(logScaleFac=0,scaled=as.numeric(lambda*(1+lambda)/(1-lambda)^3)))
  } 
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,((seqn/(seqn-1L))^2)*lambda/(seqn^nu))
  cumprodfacs <- cumprod(facs)
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n2 = exp(logScaleFac)*scaled
}

.COMP_Z_n3 <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  return(.Rcpp_COMP_Z(moment=3,nu=nu,lambda=lambda,maxn=maxn))
  # R version:
  if (nu==0) {
    return(c(logScaleFac=0,scaled=as.numeric(lambda*(1+4*lambda+lambda^2)/(1-lambda)^4)))
  } 
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,((seqn/(seqn-1L))^3)*lambda/(seqn^nu))
  cumprodfacs <- cumprod(facs)
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n3 = exp(logScaleFac)*scaled
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
  if (log) {
    return(Z1[["logScaleFac"]]-Z2[["logScaleFac"]]+log(Z1[["scaled"]])-log(Z2[["scaled"]]))
  } else return(exp(Z1[["logScaleFac"]]-Z2[["logScaleFac"]])*Z1[["scaled"]]/Z2[["scaled"]])
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
  num <- .COMP_Z_n(lambda=lambda,nu=nu)
  denum <- .COMP_Z(lambda=lambda,nu=nu)
  mu <- .COMP_Z_ratio(num,denum)
  if ( ! is.finite(mu) && lambda>10^nu) { ## FR->FR heuristic
    mu <- lambda^(1/nu)-(nu-1)/(2*nu)
  }  
  return(max(mu,1e-8)) ## avoids mu=0 which fails validmu(mu) (as for poisson family)
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
    lambdamax <- max(c(.Machine$double.eps, (app_lambda+1)*c(1.01))) 
    # last one for low lambda,nu values
    fupper <- CMP_linkfun_objfn(lambdamax, mu=mu)  ## normally >0 
    if (is.nan(fupper) || fupper<0) {
      if ( ! identical(spaMM.getOption("COMP_geom_approx_warned"),TRUE)) {
        warning(paste("Geometric approximation tried in COMPoisson(nu=",nu,") for mu=",mu," and possibly for other nu,mu values"))
        spaMM.options(COMP_geom_approx_warned=TRUE)
      }
      lambda <- mu/(1+mu) 
    } else {
      flower <- CMP_linkfun_objfn(lambdamin, mu=mu) ## normally <0 
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
  if (! is.null(mu <- attr(lambda,"mu"))) {
    attributes(lambda) <- NULL
    return(structure(mu,lambda=lambda))
  }
  if (nu==0) {
    mus <- lambda/(1-lambda) 
  } else {
    nu <- parent.env(environment())$nu
    mus <- sapply(lambda, .CMP_linkinvi,nu=nu)
  }
  dim(mus) <- dim(lambda) ## may be NULL
  return(mus)
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
  if (nu==0) {
    resu <- lambda/(1-lambda)^2
  } else {
    nu <- parent.env(environment())$nu
    resu <- numeric(length(lambda))
    for (it in seq_len(length(lambda))) {
      lambdai <- lambda[it]
      compz <- .COMP_Z(lambda=lambdai,nu=nu)
      compzn <- .COMP_Z_n(lambda=lambdai,nu=nu)
      compzn2 <- .COMP_Z_n2(lambda=lambdai,nu=nu)
      rn2 <- .COMP_Z_ratio(compzn2,compz)
      rn <- .COMP_Z_ratio(compzn,compz) ## mu
      res <- rn2 - rn^2 # (compz*compzn2-compzn^2)/(compz^2)
      # dmu/deta=(dmu/dlam) (dlam/deta) = lam dmu/dlam = lam (compz*compzn2-compzn^2)/(lam compz^2) = resu
      # b/c mu = compzn/compz, d compz/ dlam = compzn/lam, d compzn/ dlam = compzn2/lam
      resu[it] <- pmax(res, .Machine$double.eps)
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
    En <- En2 <- numeric(length(lambdas))
    for (it in seq_len(length(lambdas))) {
      lambda <- lambdas[it]
      En[it] <- .COMP_Z_ratio(.COMP_Z_n(lambda=lambda,nu=nu),.COMP_Z(lambda=lambda,nu=nu))
      En2[it] <- .COMP_Z_ratio(.COMP_Z_n2(lambda=lambda,nu=nu),.COMP_Z(lambda=lambda,nu=nu))
    }
    resu <- pmax(En2-En^2,1e-8) ## pmax otherwise for low mu, Vmu=0, -> ... -> w.resid=0
    ## for geom V(mu) = mu(1+mu) but this cannot be a lower bound for resu.
    return(resu)
  }
}

.CMP_linkfun <- function(mu, ## scalar or vector
                         log=TRUE) { ## log=TRUE => returns eta; else returns lambda, with mu attribute
  if ( is.null(lambdas <- attr(mu,"lambda"))) {
    nu <- parent.env(environment())$nu
    if (nu==0) {
      lambdas <- mu/(1+mu) 
    } else lambdas <- sapply(mu, .CMP_linkfuni,nu=nu, CMP_linkfun_objfn=parent.env(environment())$CMP_linkfun_objfn)
  } else attributes(mu) <- NULL ## avoids 'mise en abime'
  if (log) {
    return(log(lambdas)) ## eta, ie standard linkfun value
  } else return(structure(lambdas,mu=mu))
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
  CMP_linkfun_objfn <- function(lambda, mu) {linkinv(lambda=lambda) -mu}
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

