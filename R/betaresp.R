beta_resp <- function (prec = stop("beta_resp's 'prec' must be specified"), link = "logit") {
  trunc <- FALSE
  mc <- match.call()
  if (inherits(shch <- substitute(prec),"character") ||
      (inherits(shch,"name") && inherits(prec, "function")) # "name" is for e.g. beta_resp(logit)
      # (but testing only "name" would catch e.g. negbin(prec=prec) )
  ) { 
    if (inherits(shch,"character")) shch <- paste0('"',shch,'"')
    errmess <- paste0('It looks like beta_resp(',shch,') was called, which absurdly means betaresp(prec=',shch,
                      ').\n  Use named argument: beta_resp(link=',shch,') instead.')
    stop(errmess)
  }
  # When 'shape' is recognized as as call to some function ! = stop(), we eval it so it is no longer recognized as a call by .calc_optim_args()
  if (inherits(shch,"call") && deparse(shch[[1]])!="stop") prec <- eval(shch) 
  
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  if (linktemp %in% c("logit", "probit", "cloglog", "cauchit")) 
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
    else stop(gettextf("\"%s\" link not available for beta response family; available links are \"logit\", \"probit\", \"cloglog\" and \"cauchit\"", 
                       linktemp))
  }

  if ( ! is.integer(trunc)) {trunc <- round(trunc)}
  variance <- function(mu) mu*(1-mu)/(1+prec)
  validmu <- function(mu) all(mu > 0 & mu < 1)
  
  logl <- function(y, mu, wt) {
    y <- drop(y)
    prec <- prec*wt 
    lgamma(prec) - lgamma(mu*prec) - 
      lgamma((1 - mu)*prec) + (mu*prec - 1)*log(y) + ((1 - mu)*prec - 1)*log(1 - y)
  }
  DlogLDmu <- function(y, mu, wt, n) { # dlogL/dmu
    y <- drop(y)
    prec <- prec*wt 
    # prec (-Log[1 - y] + Log[y] - PolyGamma[0, prec mu] + PolyGamma[0, prec - prec mu])
    prec *( log(y)-log(1-y) - digamma(prec * mu) + digamma(prec * (1-mu)))
  }
  D2logLDmu2 <- function(y, mu, wt, n) { 
    y <- drop(y)
    prec <- prec*wt 
    # prec^2 (PolyGamma[1, prec mu] + PolyGamma[1, prec - prec mu])
    - prec^2 *(trigamma(prec*mu)+trigamma(prec * (1-mu)))
  }
  D3logLDmu3 <- function(y, mu, wt, n) { 
    y <- drop(y)
    prec <- prec*wt 
    prec^3 *(- psigamma(prec*mu, deriv=2)+psigamma(prec * (1-mu), deriv=2)) # deriv=2 -> tetragamma, cf doc
  }
  
  get_sat_mu <- function(y, wt) {
    if (y>1/2) {
      interval <- c(1/2, y)
    } else interval <- c(y, 1/2)
    uniroot(DlogLDmu, interval=interval, y = y, wt=wt)$root
  }
  
  sat_logL <- function(y, wt) { 
    nr <- length(y)
    mu <- numeric(nr)
    for (it in seq_len(nr)) mu[it] <- get_sat_mu(y[it], wt=wt[it])
    logl(y,mu=mu, wt=wt)
  } 
  
  dev.resids <- function(y,mu,wt) { 2*(sat_logL(y, wt=wt)-logl(y,mu=mu,wt=wt)) } # cannot use $aic() which is already a sum...
  
  aic <- function(y, n, mu, wt, dev) {
    - 2 * sum(logl(y, mu, wt))
  }
  initialize <- expression({
    if (any(y <= 0 | y >= 1)) stop("y values must be 0 < y < 1")
    mustart <- y
  })
  simfun <- function(object, nsim) {
    wts <- object$prior.weights
    if (any(wts != 1)) 
      warning("ignoring prior weights") # NB this function is not called by spaMM anyway. See .r_resid_var() instead. 
    # So the following is more documentation by pseudocode than anything else.
    ftd <- fitted(object)
    prec <- .get_family_par(object$family)
    rbeta(n=nsim * length(ftd), 
             shape1=ftd*prec,shape2=(1-ftd)*prec)
  }
  ## all closures defined here have parent.env the environment(spaMM_Gamma) ie <environment: namespace:spaMM>
  ## changes the parent.env of all these functions (aic, dev.resids, simfun, validmu, variance): 
  D2muDeta2 <- .D2muDeta2(linktemp)
  D3muDeta3 <- .D3muDeta3(linktemp)
  parent.env(environment(aic)) <- environment(stats::binomial) ## parent = <environment: namespace:stats>; 
  # before the change this is <environment: namespace:spaMM>, necessary to access .D2muDeta2...
  structure(list(family = structure("beta_resp",
                                    withArgs=quote(paste0("beta_resp(prec=",signif(prec,4),")"))), 
                 link = linktemp, linkfun = stats$linkfun, 
                 linkinv = stats$linkinv, variance = variance,  sat_logL=sat_logL, logl=logl, dev.resids=dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun, 
                 DlogLDmu = DlogLDmu, D2logLDmu2 = D2logLDmu2, D3logLDmu3 = D3logLDmu3, 
                 D2muDeta2 = D2muDeta2, D3muDeta3 = D3muDeta3,
                 flags=list(obs=TRUE, exp=FALSE, canonicalLink=FALSE)), 
            class = c("LLF","family"))
}
