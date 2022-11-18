.negbin1_p0 <- function(mu,shape) {
  if (shape>1e4) {
    exp(shape*mu * log1p(-1/shape)) 
  } else (shape/(1 + shape))^(shape*mu) 
}



negbin1 <- function (shape = stop("negbin1's 'shape' must be specified"), link = "log", trunc=-1L) {
  mc <- match.call()
  if (inherits(shch <- substitute(shape),"character") ||
      (inherits(shch,"name") && inherits(shape, "function")) # "name" is for e.g. negbin(log)
      # (but testing only "name" would catch e.g. negbin(shape=shape) )
  ) { 
    if (inherits(shch,"character")) shch <- paste0('"',shch,'"')
    errmess <- paste0('It looks like negbin1(',shch,') was called, which absurdly means negbin1(shape=',shch,
                      ').\n  Use named argument: negbin1(link=',shch,') instead.')
    stop(errmess)
  }
  # When 'shape' is recognized as as call to some function ! = stop(), we eval it so it is no longer recognized as a call by .calc_optim_args()
  if (inherits(shch,"call") && deparse(shch[[1]])!="stop") shape <- eval(shch) 
  
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt")) ## all non-canonical
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name)) 
        linktemp <- stats$name
    }
    else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"", 
                       linktemp))
  }

  if ( ! is.integer(trunc)) {trunc <- round(trunc)}
  variance <- function(mu) mu + mu/shape
  validmu <- function(mu) all(mu > 0)
  # logl and derived functions  derived from negbin[2] b replacing shape by shape*mu
  # - negbin1(shape=3)$aic(y,NA,mu,wt=1)/2  = dnbinom(x=y,size=shape*mu,mu=mu,log=TRUE)
  # compared to
  # - negbin(shape=3)$aic(y,NA,mu,wt=1)/2  = dnbinom(x=y,size=shape,mu=mu,log=TRUE)
  if (trunc==0L) { # if truncated
    logl <- function(y, mu, wt) {
      mushape <- mu*shape
      term <- (y + mushape) * log(mu + mushape) - y * log(mu) + 
        lgamma(y + 1) - mushape * log(mushape) + lgamma(mushape) - lgamma(mushape + y) +
        log(1-(shape/(1 + shape))^(shape*mu)) ## log(1-p0)
      term <- drop(term)
      term[mu==0] <- y[mu==0]*log(1+shape) +log(y[mu==0]) +log( log(1+1/shape)) 
      - term * wt
    }
    DlogLDmu_0series <- function(mu,y) { # series for dlogL near mu=0
      # th PolyGamma[0, y] + 
      #   1/12 th (12 EulerGamma - 2 \[Pi]^2 th \[Mu] - th \[Mu] Log[th]^2 + 
      #              2 Log[th] (3 + th \[Mu] Log[1 + th]) - 
      #              Log[1 + th] (6 + th \[Mu] Log[1 + th]) + 
      #              6 th \[Mu] (2 PolyGamma[1, y] + 
      #                            th \[Mu] (PolyGamma[2, y] + 2 Zeta[3])))
      mushape <- mu*shape
      lsh <- log(shape)
      l1sh <- log(1+shape)
      shape*digamma(y) + shape*(12*0.57721566490-2*pi^2*mushape - mushape*lsh^2 +
                                  2*lsh*(3+mushape*l1sh) - l1sh*(6+mushape*l1sh) +
                                  6*mushape* (2*trigamma(y)+mushape*(psigamma(y, deriv=2)+2.404113806319188))  )/12
    }
    DlogLDmu <- function(mu, y, wt, n, phi) { # dlogL/dmu 
      # shape (Log[mu shape] - Log[mu (1 + shape)] - PolyGamma[0, mu shape] + PolyGamma[0, mu shape + y])
      mushape <- mu*shape
      term <- shape * (log(shape/(1 + shape)) - digamma(mushape) + digamma(mushape + y))
      term <- drop(term)
      p0 <- .negbin1_p0(mu,shape)
      Mdlog1mp0 <- shape*p0* log(shape/(1 + shape))/(1-p0)
      dlogl <- term + Mdlog1mp0 # difference of two diverging terms as mu->0
      if (any(mu<1e-4)) {
        if (length(mu)>1L) {
          if (length(y)>1L) {
            dlogl[mu<1e-4] <- DlogLDmu_0series(mu[mu<1e-4],y[mu<1e-4])
          } else dlogl[mu<1e-4] <- DlogLDmu_0series(mu[mu<1e-4],y)
        } else dlogl <- DlogLDmu_0series(mu,y)
      }
      dlogl
    }
    D2logLDmu2_0series <- function(mu,y) { # series for d2logL near mu=0
      # 1/240 \[Theta]^2 (-8 \[Pi]^2 (5 + \[Pi]^2 \[Theta]^2 \[Mu]^2) - 
      #                     20 Log[\[Theta]/(
      #                       1 + \[Theta])]^2 + \[Theta]^2 \[Mu]^2 Log[\[Theta]/(
      #                         1 + \[Theta])]^4 + 240 PolyGamma[1, y] + 
      #                     120 \[Theta] \[Mu] (2 PolyGamma[2, y] + \[Theta] \[Mu] PolyGamma[3,
      #                                                                                      y] + 4 Zeta[3]))
      mushape <- mu*shape
      lsh <- log(shape)
      l1sh <- log(1+shape)
      shape^2* (-8*pi^2* (5+(pi*mushape)^2) -
                  20*log(shape/(1 + shape))^2 + mushape^2 *log(shape/(1 + shape))^4 + 240*trigamma(y) +
                  120* mushape*(2*psigamma(y, deriv=2)+mushape*psigamma(y, deriv=3)))/240
    }
    D2logLDmu2 <- function(mu, y, wt, n, phi) { 
      # shape^2 (-PolyGamma[1, mu shape] + PolyGamma[1, mu shape + y])
      mushape <- mu*shape
      term <- shape^2 * (- trigamma(mushape) + trigamma(mushape + y))
      term <- drop(term)
      p0 <- .negbin1_p0(mu,shape)
      Md2log1mp0 <- p0* (shape*log(shape/(1 + shape))/(1-p0))^2
      d2logl <- term + Md2log1mp0 # difference of two diverging terms as mu->0 (even worse than for dlogl)
      if (any(mu<1e-2)) {
        if (length(mu)>1L) {
          if (length(y)>1L) {
            d2logl[mu<1e-2] <- D2logLDmu2_0series(mu[mu<1e-2],y[mu<1e-2])
          } else d2logl[mu<1e-2] <- D2logLDmu2_0series(mu[mu<1e-2],y)
        } else d2logl <- D2logLDmu2_0series(mu,y)
      }
      d2logl
    }
    D3logLDmu3 <- function(mu, y, wt, n, phi) { 
      # shape^2 (-PolyGamma[1, mu shape] + PolyGamma[1, mu shape + y])
      mushape <- mu*shape
      term <- shape^3 * (- psigamma(mushape, deriv=2) + psigamma(mushape + y, deriv=2))
      term <- drop(term)
      p0 <- .negbin1_p0(mu,shape)
      Md3log1mp0 <- p0*(1+p0)*(shape*log(shape/(1 + shape))/(1-p0))^3
      term + Md3log1mp0
    }
  } else {
    logl <- function(y, mu, wt) {
      mushape <- mu*shape
      term <- (y + mushape) * log(mu + mushape) - y * log(mu) + 
        lgamma(y + 1) - mushape * log(mushape) + lgamma(mushape) - lgamma(mushape + y)
      term <- drop(term)
      term[y==0L & mu==0] <- 0 # replaces NaN's with correct answer
      - term * wt
    }
    DlogLDmu <- function(mu, y, wt, n, phi) { # dlogL/dmu
      # shape (Log[mu shape] - Log[mu (1 + shape)] - PolyGamma[0, mu shape] + PolyGamma[0, mu shape + y])
      mushape <- mu*shape
      drop(shape * (log(shape/(1 + shape)) - digamma(mushape) + digamma(mushape + y)))
    }
    D2logLDmu2 <- function(mu, y, wt, n, phi) { 
      # shape^2 (-PolyGamma[1, mu shape] + PolyGamma[1, mu shape + y])
      mushape <- mu*shape
      drop(shape^2 * (- trigamma(mushape) + trigamma(mushape + y)))
    }
    D3logLDmu3 <- function(mu, y, wt, n, phi) { 
      # shape^2 (-PolyGamma[1, mu shape] + PolyGamma[1, mu shape + y])
      mushape <- mu*shape
      drop(shape^3 * (- psigamma(mushape, deriv=2) + psigamma(mushape + y, deriv=2)))
    }
  }
  
  DlogLDmu_0 <- function(y) { #  dlogL in mu=0
    shape*digamma(y) + shape*(0.57721566490+log(shape)/2 - log(1+shape)/2)
  }
  
  sat_logL <- function(y, wt, return_logL=TRUE) { 
    shapeconst <- - 1/(shape * log(shape/(1 + shape)))
    uniqy <- unique(y)
    uniqmu <- rep(NA, length(uniqy))
    uniqmu[uniqy == 0L] <- 0
    uniqmu[uniqy == 1L] <- shapeconst
    getmu <- function(y) {
      if (trunc==0L) {
        if (y==1L) return(0)
        if (DlogLDmu_0(y)<0) return(0) # DlogL < 0 in mu->0   => return mu=0
        if (shape >= 1) {
          lower <-  y-1
        } else {
          # joint series for dlogl for small mu and shape yields the following solution for mu
          # (6 (2 EulerGamma - th + Log[th] + 2 PolyGamma[0, y]))/(th (2 \[Pi]^2 + Log[th]^2 - 12 PolyGamma[1, y]))
          sol_approx <- 6*(2*0.57721566490 - shape +log(shape)+2* digamma(y) )/(shape * (2*pi^2+log(shape)^2-12*trigamma(y)))
          # seems to be either very close (low shape) or an underestimate of the solution, but let's play safe:
          lower <- 0.99 * max(0,sol_approx)
        }
      } else if (shape >= 1) {
        lower <-  y/shape
      } else lower <-  y
      interval <- c(lower, y * shapeconst)
      uniroot(DlogLDmu, interval=interval, y = y)$root
    }
    ygt1 <- (uniqy > 1L)
    uniqmu[ygt1] <- sapply(uniqy[ygt1], getmu)
    if (return_logL) {
      uniqlogl <- logl(uniqy,mu=uniqmu,wt=1)
      uniqlogl[match(y, uniqy)]
    } else uniqmu[match(y, uniqy)]
  } 
  
  dev.resids <- function(y,mu,wt) { 2*(sat_logL(y, wt=wt)-logl(y,mu=mu,wt=wt)) } # cannot use $aic() which is already a sum...
  
  aic <- function(y, n, mu, wt, dev) {
    - 2 * sum(logl(y, mu, wt))
  }
  linkfun <- function(mu,mu_truncated=FALSE) { ## mu_truncated gives type of I N put
    if (mu_truncated) { ## ie if input mu_T
      mu_U <- attr(mu,"mu_U")
      return(stats$linkfun(mu_U))
    } else return(stats$linkfun(mu)) ## mu_U -> eta_U
  }
  linkinv <- function(eta,mu_truncated=FALSE) { ## mu_truncated gives type of O U T put
    if (mu_truncated) { ## ie if return mu_T
      mu_U <- stats$linkinv(eta) ## eta is always eta_U
      p0 <- (shape/(1 + shape))^(shape*mu_U) # cf TruncatedNegBin.nb
      mu_T <- mu_U/(1-p0) 
      return(structure(mu_T,p0=p0,mu_U=mu_U))
    } else return(stats$linkinv(eta))
  }
  initialize <- expression({
    if (environment(aic)$trunc==0L) {
      if (any(y < 1L)) stop("Values <1 are not allowed for the truncated negative binomial family")
    } else if (any(y < 0L)) stop("negative values not allowed for the negative binomial family")
    n <- rep(1, nobs)
    mustart <- y + (y == 0)/6
    # mustart <- family$sat_logL(y, return_logL=FALSE) + (y == 0)/6
    # more like:
    # mustart <- y + 1/(6*environment(aic)$shape)                 + (y == 0)/6 # but still add the y==0 correction otherwie this fails badly
  })
  simfun <- function(object, nsim,zero_truncated=identical(object$family$zero_truncated,TRUE)) {
    wts <- object$prior.weights
    if (any(wts != 1)) 
      warning("ignoring prior weights")
    ftd <- fitted(object)
    .rnbinom(n=nsim * length(ftd), 
             size=ftd * shape, # c comment on logL-relate functions
             mu_str=ftd, zero_truncated=zero_truncated)
  }
  ## No ad-hoc functions dlW_deta etc. That would be ugly. The mandatory function are the D<n>logLDmu<n> functions,
  ## used to compute $Md2logcLdeta2 and $Md3logcLdeta3 in .muetafn(), and then generic code:
  ## dlW_deta <- muetablob$Md3logcLdeta3/muetablob$Md2logcLdeta2 
  ## if (calcCoef1) coef1 <- dlW_deta/muetablob$Md2logcLdeta2
  D2muDeta2 <- .D2muDeta2(linktemp)
  D3muDeta3 <- .D3muDeta3(linktemp)
  ## all closures defined here have parent.env <environment: namespace:spaMM>
  ## => Change the parent.env of all these functions (aic, dev.resids, simfun, validmu, variance): 
  # parent.env(environment(aic)) <- environment(stats::binomial) ## parent = <environment: namespace:stats>; 
  # before the change this is <environment: namespace:spaMM>, necessary to access .D2muDeta2...
  structure(list(family = structure("negbin1",
                                    withArgs=quote(paste0("negbin1(shape=",signif(shape,4),")"))), 
                 link = linktemp, linkfun = linkfun, 
                 linkinv = linkinv, variance = variance, sat_logL=sat_logL, logl=logl, dev.resids=dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun, 
                 DlogLDmu = DlogLDmu, D2logLDmu2 = D2logLDmu2, D3logLDmu3 = D3logLDmu3, 
                 D2muDeta2 = D2muDeta2, D3muDeta3 = D3muDeta3,
                 flags=list(obs=TRUE, exp=FALSE, canonicalLink=FALSE, LLgeneric=TRUE),
                 zero_truncated=(trunc==0L)), 
            class = c("LLF","family"))
}
