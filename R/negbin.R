.TNB_dev_integrand <- function(t) {} 

.rnbinom <- function(n, size, mu_str, zero_truncated=FALSE) { ## mu is either the standard mean (lambda) or a structure depending on trunc
  if (zero_truncated) {
    p0 <- attr(mu_str,"p0") ## or recompute it ? 
    Tunif <- runif(n,min=p0,max=1) ## truncated uniform
    resu <- qnbinom(Tunif, size=size, mu=attr(mu_str,"mu_U"))
  } else resu <- rnbinom(n, size=size, mu=mu_str)  
  return(resu)
}


negbin <- function (shape = stop("negbin's 'shape' must be specified"), link = "log", trunc=-1L) {
  mc <- match.call()
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
  variance <- function(mu) mu + mu^2/shape
  validmu <- function(mu) all(mu > 0)
  if (trunc==0L) {
    dev.resids <- function(y, mu, wt) {
      ## the dev.resids serves in MM to estimate phi (hence not here) or lambda in ranCoefs models (=> trunc.neg.bin with ranCoefs)
      # computation is in Mathelatica notebook.
      2 * wt * (y * log(pmax(1, y)/mu) - (y + shape) * log((y + shape)/(mu + shape))
                +log( (1-(shape/(mu+shape))^shape) / (1-(shape/(y+shape))^shape) )
      )
    }
    aic <- function(y, n, mu, wt, dev) { ## not really the aic... -2 clik
      term <- (y + shape) * log(mu + shape) - y * log(mu) + 
        lgamma(y + 1) - shape * log(shape) + lgamma(shape) - lgamma(shape + y) + 
        + log(1-(shape/(mu + shape))^shape) ## log(1-p0)
      2 * sum(term * wt)
    }
  } else {
    dev.resids <- function(y, mu, wt) {
      2 * wt * (y * log(pmax(1, y)/mu) - (y + shape) * log((y + shape)/(mu + shape)))
    }
    aic <- function(y, n, mu, wt, dev) {
      term <- (y + shape) * log(mu + shape) - y * log(mu) + 
        lgamma(y + 1) - shape * log(shape) + lgamma(shape) - lgamma(shape + y)
      2 * sum(term * wt)
    }
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
      p0 <- (shape/(mu_U + shape))^shape
      mu_T <- mu_U/(1-p0) 
      return(structure(mu_T,p0=p0,mu_U=mu_U))
    } else return(stats$linkinv(eta))
  }
  initialize <- expression({
    if (environment(aic)$trunc==0L) {
      if (any(y < 1L)) stop("Values <1 are not allowed for the truncated negative binomial family")
    } else if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
    n <- rep(1, nobs)
    mustart <- y + (y == 0)/6
  })
  simfun <- function(object, nsim,zero_truncated=identical(object$family$zero_truncated,TRUE)) {
    wts <- object$prior.weights
    if (any(wts != 1)) 
      warning("ignoring prior weights")
    ftd <- fitted(object)
    .rnbinom(n=nsim * length(ftd), size=shape, mu_str=ftd, zero_truncated=zero_truncated)
  }
  ## all closures defined here have parent.env the environment(spaMM_Gamma) ie <environment: namespace:spaMM>
  ## changes the parent.env of all these functions (aic, dev.resids, simfun, validmu, variance): 
  parent.env(environment(aic)) <- environment(stats::binomial) ## parent = <environment: namespace:stats>
  structure(list(family = structure("negbin",
                                    withArgs=quote(paste("Neg.binomial(shape=",signif(shape,4),")",sep=""))), 
                 link = linktemp, linkfun = linkfun, 
                 linkinv = linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun,
                 zero_truncated=(trunc==0L),predict=predict), 
            class = "family")
}

Tnegbin <- function(shape = stop("negbin's 'shape' must be specified"), link = "log") {
  mc <- match.call()   # avoid evaluation of promise...
  mc[[1L]] <- quote(spaMM::negbin)
  mc[["trunc"]] <- 0L
  eval(mc, parent.frame())
}

