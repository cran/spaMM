.DlogLDmu_trunc_nb2 <- local({
  function(y, mu, wt, n, phi, shape_it=NULL) { # dlogL/dmu
    if ( ! is.null(shape_it)) shape <- shape_it
    term <- shape*(y-mu)/(mu*(mu + shape))
    term <- drop(term)
    p0 <- .negbin2_p0(mu,shape)
    Mdlog1mp0 <- - p0^(1+1/shape) /(1-p0)
    term + Mdlog1mp0
  }
})

.D2logLDmu2_trunc_nb2 <- local({
  shape <- NaN
  function(y, mu, wt, n, phi) { 
    term <-  shape *(mu^2 - y* (shape+2*mu))/(mu*(shape+mu))^2 
    term <- drop(term)
    p0 <- .negbin2_p0(mu,shape)
    Md2log1mp0 <- shape* p0 * (1+shape-p0)/ ((shape+mu)*(1-p0))^2
    term + Md2log1mp0
  }
})

.D3logLDmu3_trunc_nb2 <- local({
  shape <- NaN
  function(y, mu, wt, n, phi) { # element of computation of D3logcLdeta3 for d logdet Hessian
    term <- ( 2*shape*(-mu^3 + y *(shape^2 + 3*shape*mu + 3*mu^2)) ) / (mu^3 *(shape + mu)^3)
    term <- drop(term)
    p0 <- .negbin2_p0(mu,shape)
    Md3log1mp0 <- -shape*p0*(3*shape*(1 - p0) + 2*(1 - p0)^2 + 
                               shape^2*(1 + p0))/((shape + mu)*(1 - p0))^3
    term + Md3log1mp0
  }
})



# Far from perfect: how to substitute a warning to another ?
.rnbinom_warn <- function(cond) {
  if (grep("NA", conditionMessage(cond))) {
    enclosing_env <- parent.env(environment())
    mu_str <- get("mu_str",enclosing_env)
    if (any(is.infinite(mu_str))) { # # (n=10, size=0.1, mu=1e309)
      message("Infinite 'mu' in rnbinom() call.") 
    } else {
      size <- get("size",enclosing_env)
      if (all(mu_str>0) && all(size>0)) { # (n=10, size=0.1, mu=1e308)
        message("NAs produced by rnbinom() despite finite expectation: large 'mu' & small 'size' can result in infinite values of latent gamma draw. ")
        # pb is at level of base C sources.
      } else withRestarts({
        signalCondition(cond)
      }, muffleWarning = function() NULL)
    }
  } 
}

.rnbinom <- function(n, size, mu_str, zero_truncated=FALSE) { ## mu is either the standard mean (lambda) or a structure depending on trunc
  if (zero_truncated) {
    p0 <- attr(mu_str,"p0") ## or recompute it ? 
    Tunif <- runif(n,min=p0,max=1) ## truncated uniform
    resu <- qnbinom(Tunif, size=size, mu=attr(mu_str,"mu_U"))
  } else {
    ## rnbinom() is rpois(rgamma(shape=size,scale=mu/shape))
    ## if rgamma generates Inf [hard to control], rpois() then produces NA
    environment(.rnbinom_warn) <- environment()
    resu <- withCallingHandlers(rnbinom(n, size=size, mu=mu_str), warning=.rnbinom_warn)
  }
  return(resu)
}

.negbin2_p0 <- function(mu,shape) {
  p0 <- (shape/(mu+shape))^shape
  bad <- mu<shape/1e4
  p0[bad] <- exp(shape * log1p(-mu[bad]/shape)) # exp(-mu)
  p0
}

negbin2 <- function (shape = stop("negbin2's 'shape' must be specified"), link = "log", trunc=-1L, LLgeneric=TRUE) {
  mc <- match.call()
  resid.model <- list2env(list(off=0)) # env so that when we assign to it we don't create a new instance of the family object
  if (inherits(shch <- substitute(shape),"character") ||
      (inherits(shch,"name") && inherits(shape, "function")) # "name" is for e.g. negbin(log)
      # (but testing only "name" would catch e.g. negbin(shape=shape) )
  ) { 
    if (inherits(shch,"character")) shch <- paste0('"',shch,'"')
    errmess <- paste0('It looks like negbin[2](',shch,') was called, which absurdly means negbin[2](shape=',shch,
                      ').\n  Use named argument: negbin[2](link=',shch,') instead.')
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
  variance <- function(mu, new_fampar=NULL) {
    if ( ! is.null(new_fampar)) shape <- new_fampar # for .calcResidVar with newdata an a resid.model 
    mu + mu^2/shape
  }
  validmu <- function(mu) all(mu > 0)
  
  dlW_Hexp__detafun <- switch( # of *Hexp* weights, used for expInfo, as well as obsInfo (as explained for the calling function .dlW_Hexp__dmu in the latter case)
    linktemp,
    "log"= function(mu) shape/(mu+shape),
    "identity"= function(mu) -(2*mu + shape)/(mu^2 + mu*shape),
    "sqrt"= function(mu) -2 * sqrt(mu)/(mu + shape)
  )
  
  coef1fun <- switch(
    linktemp,
    "log"= function(mu) 1/mu,
    "identity"= function(mu) - (1+2*mu/shape),
    "sqrt"= function(mu) - sqrt(mu)/(2*shape)
  )
  
  if (trunc==0L) {
    # logl expected in all cases, with different definitions (=> single def of aic in terms of logl)
    # sat_logL not specific to obsInfo, but needed for finding saturated model when there is trucnation
    logl <- function(y, mu, wt) {
      mlogl <- (y + shape) * log(mu + shape) - y * log(mu) + 
        lgamma(y + 1) - shape * log(shape) + lgamma(shape) - lgamma(shape + y) +
        log(1-.negbin2_p0(mu,shape)) 
      mlogl[y==1L & mu==0] <- 0 
      drop(- mlogl * wt)
    }
    
    getmu <- function(y, shape_it=NULL) { # it's useful to keep shape_it=NULL fro the call to DlogLDmu() ? not sure, the latter
      if ( ! is.null(shape_it)) shape <- shape_it
      if (y==1L) return(0)
      if (shape >= 1) {
        lower <-  (y-1)*0.999
      } else lower <-  ( y-1 )*shape*0.999
      interval <- c(lower, y) # derivative always <0 in mu=y. For shape<1, perhaps a stricter upper bound could be found
      uniroot(DlogLDmu, interval=interval, y = y, shape_it=shape_it)$root
    }
    
    sat_logL <- function(y, wt) { # This is for the 0--truncated family so no y=0L case
      if (length(shape)>1L) {
        muv <- rep(NA, length(y))
        muv[y == 1L] <- 0
        for (it in which(y > 1L)) muv[it] <- getmu(y[it], shape_it=shape[it])
        logl(y,mu=muv,wt=1) # on vector y, muv and shape 
      } else {
        uniqy <- unique(y)
        uniqmu <- rep(NA, length(uniqy))
        uniqmu[uniqy == 1L] <- 0
        for (it in which(uniqy > 1L)) uniqmu[it] <- getmu(uniqy[it])
        uniqlogl <- logl(uniqy,mu=uniqmu,wt=1)
        uniqlogl[match(y, uniqy)]
      }
    } 
    
    ## the dev.resids serves in llm.fit...
    dev.resids <- function(y,mu,wt) { 2*(sat_logL(y, wt=wt)-logl(y,mu=mu,wt=wt)) } # cannot use $aic() which is already a sum...
    
  } else { # untruncated
    
    logl <- function(y, mu, wt) {
      term <- (y + shape) * log(mu + shape) - y * log(mu) + 
        lgamma(y + 1) - shape * log(shape) + lgamma(shape) - lgamma(shape + y)
      drop(- term * wt)
    }
    
    sat_logL <- NULL # its a standard GLM case where saturated mu = y
    
    dev.resids <- function(y, mu, wt) {
      2 * wt * (y * log(pmax(1, y)/mu) - (y + shape) * log((y + shape)/(mu + shape)))
    }
  }
  
  aic <- function(y, n, mu, wt, dev) { - 2 * sum(logl(y, mu, wt)) }
  
  linkfun <- function(mu,mu_truncated=FALSE) { ## mu_truncated gives type of I N put
    if (mu_truncated) { ## ie if input mu_T
      mu_U <- attr(mu,"mu_U")
      return(stats$linkfun(mu_U))
    } else return(stats$linkfun(mu)) ## mu_U -> eta_U
  }
  linkinv <- function(eta,mu_truncated=FALSE) { ## mu_truncated gives type of O U T put
    if (mu_truncated) { ## ie if return mu_T
      mu_U <- stats$linkinv(eta) ## eta is always eta_U
      p0 <- .negbin2_p0(mu_U,shape)
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
  })
  #
  simulate <- function(object, nsim,zero_truncated=identical(object$family$zero_truncated,TRUE)) { # cf comments on the beta_resp's simulate
    wts <- object$prior.weights
    if (any(wts != 1)) 
      warning("ignoring prior weights")
    ftd <- fitted(object)
    .rnbinom(n=nsim * length(ftd), size=shape, mu_str=ftd, zero_truncated=zero_truncated)
  }
  Dtheta.Dmu <- function(mu) 1/(mu*(1+mu/shape))
  D2theta.Dmu2 <- function(mu) -(1+2*mu/shape)/(mu*(1+mu/shape))^2
  #

  # for all obsInfo code:
  D2muDeta2 <- .D2muDeta2(linktemp)
  D3muDeta3 <- .D3muDeta3(linktemp)
  
  
  if (LLgeneric) {
    d_dlWdmu_detafun <- NULL
    if (trunc==0L) {
      DlogLDmu <- .DlogLDmu_trunc_nb2
      D2logLDmu2 <- .D2logLDmu2_trunc_nb2
      D3logLDmu3 <- .D3logLDmu3_trunc_nb2
    } else {
      ## link-independent function #####################
      DlogLDmu <- function(mu, y, wt, n, phi,shape_it=NULL) { 
        # The DlogLDmu() code, and higher derivatives, typically work and returns a vector, when shape is a vector taken from the function efinition environment
        # The shape_it argument is only needed when solving for DlogL=0 (uniroot) call; that's why it is not needed in functiosn for higher derivatives.
        if ( ! is.null(shape_it)) shape <- shape_it
        drop(shape*(y-mu)/(mu*(mu + shape)))
      }
      
      D2logLDmu2 <- function(mu, y, wt, n, phi) { # element of computation of Hobs weights, (-) D2logcLdeta2
        # (\[Theta] (\[Mu]^2 - y (\[Theta] + 2 \[Mu])))/(\[Mu]^2 (\[Theta] + \[Mu])^2)
        drop(shape *(mu^2 - y* (shape+2*mu))/(mu*(shape+mu))^2) 
      }
      
      D3logLDmu3 <- function(mu, y, wt, n, phi) { # element of computation of D3logcLdeta3 for d logdet Hessian
        # (2 \[Theta] (-\[Mu]^3 + y (\[Theta]^2 + 3 \[Theta] \[Mu] + 3 \[Mu]^2)))/(\[Mu]^3 (\[Theta] + \[Mu])^3)
        drop( 2*shape*(-mu^3 + y *(shape^2 + 3*shape*mu + 3*mu^2)) ) / (mu^3 *(shape + mu)^3)
      }
    }
    environment(DlogLDmu) <- environment(D2logLDmu2) <- environment(D3logLDmu3) <- environment(aic) 
  } else {
    DlogLDmu <- D2logLDmu2 <- D3logLDmu3 <- NULL
    if (trunc==0L) {
      d_dlWdmu_detafun <- NULL # never available for truncated model
    } else {
      # This one is for Hobs by ad-hoc Hratio_factors method: (not available for truncated case)
      d_dlWdmu_detafun <- switch( # see Hessian_weights.nb 
        linktemp,
        "log"= function(mu) -shape *(2*mu + shape)/(mu*(mu + shape)^2),
        "identity"= function(mu) (2 *mu^2 + 2 *mu *shape + shape^2)/(mu^2 *(mu + shape)^2), 
        "sqrt"= function(mu) 2*sqrt(mu)/(mu + shape)^2 
      )
      
    }
  }
  
  
  #
  ## all closures defined here have parent.env the environment(spaMM_Gamma) ie <environment: namespace:spaMM>
  ## changes the parent.env of all these functions (aic, dev.resids, simulate, validmu, variance): 
  # parent.env(environment(aic)) <- environment(stats::binomial) ## parent = <environment: namespace:stats>
  
  structure(list(family = structure("negbin2",
                                    withArgs=quote({
                                      if ( ! is.null(resid.model$beta)) {
                                        paste0("negbin2(shape=",paste0(signif(shape[1:min(3L,length(shape))],4),collapse=" "),"...)")
                                      } else paste0("negbin2(shape=",signif(shape,4),")")
                                    })), 
                 link = linktemp, linkfun = linkfun, 
                 linkinv = linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simulate, 
                 dlW_Hexp__detafun=dlW_Hexp__detafun, 
                 coef1fun=coef1fun, # always available (needed for expInfo as well as NON-generic obsInfo) 
                 d_dlWdmu_detafun=d_dlWdmu_detafun, # NULL in most cases (except NON-generic obsInfo => untruncated)
                 DlogLDmu = DlogLDmu, D2logLDmu2 = D2logLDmu2, D3logLDmu3 = D3logLDmu3, 
                 D2muDeta2 = D2muDeta2, D3muDeta3 = D3muDeta3, # NULL if not LLgeneric
                 resid.model=resid.model, 
                 flags=list(obs=TRUE, exp=TRUE, canonicalLink=FALSE, LLgeneric=LLgeneric),
                 zero_truncated=(trunc==0L)), 
            class = c("LLF","family"))
}

negbin <- negbin2

Tnegbin <- function(shape = stop("Tnegbin's 'shape' must be specified"), link = "log") {
  mc <- match.call()   # avoid evaluation of promise...
  mc[[1L]] <- get("negbin2", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  mc[["trunc"]] <- 0L
  eval(mc, parent.frame())
}

