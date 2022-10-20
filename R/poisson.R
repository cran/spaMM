.get_family_par <- function(family) {
  famfam <- family$family
  if (famfam =="COMPoisson") {
    family_par <- environment(family$aic)$nu
  } else if (famfam =="beta_resp") {
    family_par <- environment(family$aic)$prec
  } else if (famfam  %in% c("negbin","negbin1","negbin2")) {
    family_par <- environment(family$aic)$shape
  }
  family_par
}

.rpois <- function(n,mu_str,zero_truncated=FALSE){ ## mu is either the standard mean (lambda) or a structure depending on trunc
  if (zero_truncated) {
    p0 <- attr(mu_str,"p0") ## or recompute it ? (p0= ppois(0, lambda= attr(ftd,"mu_U")))
    Tunif <- runif(n,min=p0,max=1) ## truncated uniform
    qpois(Tunif, lambda=attr(mu_str,"mu_U"))
  } else rpois(n=n, lambda=mu_str)
}


## hides stats::poisson (check poisson()$zero_truncated)
Poisson <- function (link = "log", trunc=-1L, LLgeneric=TRUE) {
  linktemp <- substitute(link)  # if link was char LHS is char ; else deparse will create a char from a language object 
  if ( ! is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt")
  if (linktemp %in% okLinks) {
    stats <- make.link(linktemp) # from char to make.link() return value
  } else if (is.character(link)) { # does not seem useful
    stats <- make.link(link)
    linktemp <- link
  } else if (inherits(link, "link-glm")) { # a make.link() object was provided
    stats <- link
    if (!is.null(stats$name)) linktemp <- stats$name
  } # at this point 'stats' is always the result of make.link and 'linktemp' is always char. 
  # In stats:: families, the following check is in an else statement that is never reached !!
  if ( ! linktemp %in% okLinks) stop(gettextf("link \"%s\" not available for poisson family; available links are %s", 
                                              linktemp, paste(sQuote(okLinks), collapse = ", ")),       domain = NA)
  if ( ! is.integer(trunc)) {trunc <- round(trunc)}
  variance <- function(mu) mu
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)
  
  # Derivatives of Hexp, all being the limit cases of the negbin when shape -> infty
  dlW_Hexp__detafun <- switch(  # of *Hexp* weights, used for expInfo, as well as obsInfo (as explained for the calling function .dlW_Hexp__dmu in the latter case)
    linktemp,
    "log"= function(mu) 1,
    "identity"= function(mu) -1/mu,
    "sqrt"= function(mu) rep(0, length(mu))
  )
  coef1fun <- switch(
    linktemp,
    "log"= function(mu) 1/mu,
    "identity"= function(mu) rep(-1, length(mu)),
    "sqrt"= function(mu) rep(0, length(mu))
  )
  
  if (trunc==0L) { # If truncated: 

    logl <- function(y, mu, wt) {
      resu <- dpois(y, mu, log = TRUE)-log(1-exp(-mu)) # I decided to ignore wt
      resu[y==1L & mu==0] <- 0
      resu
    }
    
    sat_logL <- function(y, wt) { 
      uniqy <- unique(y)
      uniqmu <- rep(NA, length(uniqy))
      uniqmu[uniqy == 1L] <- 0
      getmu <- function(y) {
        if (y==1L) return(0)
        uniroot(DlogLDmu, interval=c(y-1,y), y = y)$root
      }
      ygt1 <- (uniqy > 1L)
      uniqmu[ygt1] <- sapply(uniqy[ygt1], getmu)
      uniqlogl <- logl(uniqy,mu=uniqmu,wt=1)
      uniqlogl[match(y, uniqy)]
    } 
    
    ## the dev.resids is called by glm.fit...
    dev.resids <- function(y,mu,wt) { 2*(sat_logL(y, wt=wt)-logl(y,mu=mu,wt=wt)) } 
    
    ### MolasL10, pp.3309. This is weird bc their "definition" of deviance conflicts with a more fundamental def. 
    # dev.resids <- function(y, mu, wt) {
    #   r <- (wt * (y * log(y/mu) - (y - mu)  + log( (1-exp(-mu))/(1-exp(-y)) ) )) # MolasL p. 3309 (here and there, fn of latent mu, not mu of truncated response)
    #   2 * r
    # }
  } else {
    
    logl <- function(y, mu, wt) dpois(y, mu, log = TRUE) # I decided to ignore wt
    
    sat_logL <- NULL

    dev.resids <- function(y, mu, wt) {
      r <- mu * wt
      p <- which(y > 0)
      r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
      2 * r
    }
  }
  
  aic <- function(y, n, mu, wt, dev) -2 * sum(logl(y, mu)) 
  
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
    } else if (any(y < 0L)) stop("negative values not allowed for the 'Poisson' family")
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
  
  Dtheta.Dmu <- function(mu) 1/mu
  D2theta.Dmu2 <- function(mu) -1/mu^2
  
  # for all obsInfo code:
  D2muDeta2 <- .D2muDeta2(linktemp)
  D3muDeta3 <- .D3muDeta3(linktemp)
  
  
  #
  if (LLgeneric) {
    d_dlWdmu_detafun <- NULL
    if (trunc==0L) {
      DlogLDmu <- function(mu, y, wt, n, phi) { # dlogL/dmu
        term <- -1+drop(y)/mu
        p0 <- exp(-mu) # useful to avoir overflows of exp(mu)
        Mdlog1mp0 <- - p0/(1-p0)
        term + Mdlog1mp0
      }
      D2logLDmu2 <- function(mu, y, wt, n, phi) { 
        term <-  -drop(y)/mu^2
        p0 <- exp(-mu) # useful to avoir overflows of exp(mu)
        Md2log1mp0 <- p0/(1-p0)^2
        term + Md2log1mp0
      }
      D3logLDmu3 <- function(mu, y, wt, n, phi) { # element of computation of D3logLDeta3 for d logdet Hessian
        term <- 2*drop(y)/mu^3
        p0 <- exp(-mu) # useful to avoir overflows of exp(mu)
        Md3log1mp0 <- - p0*(1+p0)/(1-p0)^3
        term + Md3log1mp0
      }
    } else {
      ## link-independent function #####################
      DlogLDmu <- function(mu, y, wt, n, phi) -1+drop(y)/mu
      
      D2logLDmu2 <- function(mu, y, wt, n, phi) -drop(y)/mu^2
      
      D3logLDmu3 <- function(mu, y, wt, n, phi) 2*drop(y)/mu^3
    }
    environment(DlogLDmu) <- environment(D2logLDmu2) <- environment(D3logLDmu3) <- environment(aic) 
  } else {
    DlogLDmu <- D2logLDmu2 <- D3logLDmu3 <- NULL
    if (trunc==0L) {
      d_dlWdmu_detafun <- NULL # never available for truncated model
    } else {
      # This one is for Hobs by ad-hoc Hratio_factors method: (not available for truncated case)
      d_dlWdmu_detafun <- switch( # see Hessian_weights.nb # also W_H_exp
        linktemp,
        "log"= function(mu) -1/mu,
        "identity"= function(mu) 1/mu^2, 
        "sqrt"= function(mu) rep(0, length(mu))
      )
      
    }
  }
  
  parent.env(environment(aic)) <- environment(stats::poisson) ## parent = <environment: namespace:stats>
  structure(
    list(
      family = structure("poisson", patch="spaMM's *P*oisson"), 
      link = linktemp, linkfun = linkfun, 
      linkinv = linkinv, variance = variance, dev.resids = dev.resids, 
      aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
      validmu = validmu, valideta = stats$valideta, simulate = simfun, 
      dlW_Hexp__detafun=dlW_Hexp__detafun, coef1fun=coef1fun, # always available (needed for expInfo as well as NON-generic obsInfo) 
      d_dlWdmu_detafun=d_dlWdmu_detafun, # NULL in most cases (except NON-generic obsInfo => untruncated)
      DlogLDmu = DlogLDmu, D2logLDmu2 = D2logLDmu2, D3logLDmu3 = D3logLDmu3, D2muDeta2 = D2muDeta2, D3muDeta3 = D3muDeta3, # NULL if not LLgeneric
      flags=list(obs=TRUE, exp=TRUE, canonicalLink=(linktemp=="log"), LLgeneric=LLgeneric),
      zero_truncated=(trunc==0L)
    ), class = "family")
}

Tpoisson <- function(link="log") {
  # if we directly call  Poisson(link=link,trunc=0L) and the link was a language object (say log), the Poisson code only sees link, not log. 
  linktemp <- substitute(link)  # if link was char LHS is char ; else deparse will create a char from a language object 
  if ( ! is.character(linktemp)) linktemp <- deparse(linktemp)
  Poisson(link=linktemp,trunc=0L)
}

#poisson <- Poisson
