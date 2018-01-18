.get_family_par <- function(family, par, try=FALSE, ...) {
  if (try) {
    ## not sure why I used try elsewhere...
  } else environment(family$aic)[[par]]
}

.rpois <- function(n,mu_str,zero_truncated=FALSE){ ## mu is either the standard mean (lambda) or a structure depending on trunc
  if (zero_truncated) {
    p0 <- attr(mu_str,"p0") ## or recompute it ? (p0= ppois(0, lambda= attr(ftd,"mu_U")))
    Tunif <- runif(n,min=p0,max=1) ## truncated uniform
    qpois(Tunif, lambda=attr(mu_str,"mu_U"))
  } else rpois(n=n, lambda=mu_str)
}


## hides stats::poisson (check poisson()$zero_truncated)
Poisson <- function (link = "log", trunc=-1L) {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt")
  if (linktemp %in% okLinks) 
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
    else {
      stop(gettextf("link \"%s\" not available for poisson family; available links are %s", 
                    linktemp, paste(sQuote(okLinks), collapse = ", ")), 
           domain = NA)
    }
  }
  if ( ! is.integer(trunc)) {trunc <- round(trunc)}
  variance <- function(mu) mu
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)
  if (trunc==0L) {
    dev.resids <- function(y, mu, wt) {
      r <- (wt * (y * log(y/mu) - (y - mu)  + log( (1-exp(-mu))/(1-exp(-y)) ) )) # MolasL p. 3309
      2 * r
    }
    aic <- function(y, n, mu, wt, dev) -2 * sum((dpois(y, mu, log = TRUE)-log(1-exp(-mu))) * wt) ## where (-) - log() gives + log(1-p0)
  } else {
    dev.resids <- function(y, mu, wt) {
      r <- mu * wt
      p <- which(y > 0)
      r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
      2 * r
    }
    aic <- function(y, n, mu, wt, dev) -2 * sum(dpois(y, mu, log = TRUE) * wt)
  }
  linkfun <- function(mu,mu_truncated=FALSE) { ## mu_truncated gives type of I N put
    if (mu_truncated) { ## ie if input mu_T
      mu_U <- attr(mu,"mu_U")
      return(stats$linkfun(mu_U))
    } else return(stats$linkfun(mu)) ## mu_U -> eta_U
  }
  linkinv <- function(eta,mu_truncated=FALSE) {
    if (mu_truncated) {
      mu_U <- stats$linkinv(eta) 
      p0 <- exp(-mu_U)
      mu_T <- mu_U/(1-p0) 
      return(structure(mu_T,p0=p0,mu_U=mu_U))
    } else return(stats$linkinv(eta))
  }  
  initialize <- expression({
    if (environment(aic)$trunc==0L) {
      if (any(y < 1L)) stop("Values <1 are not allowed for the truncated poisson family")
    } else if (any(y < 0)) stop("negative values not allowed for the 'Poisson' family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  simfun <- function(object, nsim,zero_truncated=identical(object$family$zero_truncated,TRUE)) { 
    wts <- object$prior.weights
    if (any(wts != 1)) 
      warning("ignoring prior weights")
    ftd <- fitted(object)
    .rpois(n=nsim * length(ftd),mu_str=ftd,zero_truncated=zero_truncated)
  }
  parent.env(environment(aic)) <- environment(stats::poisson) ## parent = <environment: namespace:stats>
  structure(list(family = "poisson", link = linktemp, linkfun = linkfun, 
                 linkinv = linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun,
                 zero_truncated=(trunc==0L)), 
            class = "family")
}

Tpoisson <- function(link="log") Poisson(link=link,trunc=0L)

#poisson <- Poisson
