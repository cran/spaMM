.TNB_dev_integrand <- function(t) {} 

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

negbin <- function (shape = stop("negbin's 'shape' must be specified"), link = "log", trunc=-1L) {
  mc <- match.call()
  if (inherits(shch <- substitute(shape),"character") ||
      (inherits(shch,"name") && inherits(shape, "function")) # "name" is for e.g. negbin(log)
      # (but testing only "name" would catch e.g. negbin(shape=shape) )
  ) { 
    if (inherits(shch,"character")) shch <- paste0('"',shch,'"')
    errmess <- paste0('It looks like negbin(',shch,') was called, which absurdly means negbin(shape=',shch,
                      ').\n  Use named argument: negbin(link=',shch,') instead.')
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
  variance <- function(mu) mu + mu^2/shape
  validmu <- function(mu) all(mu > 0)
  if (trunc==0L) {
    dev.resids <- function(y, mu, wt) {
      ## the dev.resids serves in MM to estimate phi (hence not here) or lambda in ranCoefs models (=> trunc.neg.bin with ranCoefs)
      # computation is in Mathematica notebook.
      2 * wt * (y * log(pmax(1, y)/mu) - (y + shape) * log((y + shape)/(mu + shape)) +
                log( (1-.negbin2_p0(mu,shape) )/( 1-.negbin2_p0(y,shape)) )
      )
    }
    aic <- function(y, n, mu, wt, dev) { ## not really the aic... -2 clik
      term <- (y + shape) * log(mu + shape) - y * log(mu) + 
        lgamma(y + 1) - shape * log(shape) + lgamma(shape) - lgamma(shape + y) +
        log(1-.negbin2_p0(mu,shape)) ## log(1-p0)
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
  simfun <- function(object, nsim,zero_truncated=identical(object$family$zero_truncated,TRUE)) {
    wts <- object$prior.weights
    if (any(wts != 1)) 
      warning("ignoring prior weights")
    ftd <- fitted(object)
    .rnbinom(n=nsim * length(ftd), size=shape, mu_str=ftd, zero_truncated=zero_truncated)
  }
  Dtheta.Dmu <- function(mu) 1/(mu*(1+mu/shape))
  D2theta.Dmu2 <- function(mu) -(1+2*mu/shape)/(mu*(1+mu/shape))^2
  #
  #####################################################
  # Derivatives of Hexp for computation of those of Hobs by GLM ad-hoc algo
  dlW_Hexp__detafun <- switch( # of *Hexp* weights, as explained for the calling function .dlW_Hexp__dmu
    linktemp,
    "log"= function(mu) shape/(mu+shape),
    "identity"= function(mu) -(2*mu + shape)/(mu^2 + mu*shape),
    "sqrt"= function(mu) -2 * sqrt(mu)/(mu + shape)
  )
  d_dlWdmu_detafun <- switch( # see Hessian_weights.nb # also W_H_exp
    linktemp,
    "log"= function(mu) -shape *(2*mu + shape)/(mu*(mu + shape)^2),
    "identity"= function(mu) (2 *mu^2 + 2 *mu *shape + shape^2)/(mu^2 *(mu + shape)^2), 
    "sqrt"= function(mu) 2*sqrt(mu)/(mu + shape)^2 
  )
  coef1fun <- switch(
    linktemp,
    "log"= function(mu) 1/mu,
    "identity"= function(mu) - (1+2*mu/shape),
    "sqrt"= function(mu) - sqrt(mu)/(2*shape)
  )
  # => These functions are for obsInfo. For truncated obsInfo, however, another approach will be used.
  
  if (trunc==0L) { # If truncated. 
    # One might try to correct the above fns the 1/WU_WT term... but it becomes too much of a mess
    # Dl__WT_WU__Dmu <- function(mu, shape) {
    #   p0 <- .negbin2_p0(mu=mu,shape)
    #   (p0*shape*(2*(-1 + p0)*shape + mu*(-1 + p0 + shape + p0*shape)))/
    #     ((-1 + p0)*(mu + shape) ((-1 + p0)*shape + mu*(-1 + p0 + p0*shape)))
    # }
    # dlW_Hexp__detafun <- switch( # of *Hexp* weights, as explained for the calling function .dlW_Hexp__dmu
    #   # Here adding dl__WT_WU__deta computed as Dl__WT_WU__Dmu(mu, shape) * dmu/deta
    #   linktemp,
    #   "log"= function(mu) shape/(mu+shape) + Dl__WT_WU__Dmu(mu, shape)*mu,
    #   "identity"= function(mu) -(2*mu + shape)/(mu^2 + mu*shape)+ Dl__WT_WU__Dmu(mu, shape),
    #   "sqrt"= function(mu) -2 * sqrt(mu)/(mu + shape)+ Dl__WT_WU__Dmu(mu, shape)*2*sqrt(mu)
    # )
    # and it becomes worse for d_dlWdmu_detafun()
    ## INSTEAD:
    D2muDeta2 <- .D2muDeta2(linktemp)
    D3muDeta3 <- .D3muDeta3(linktemp)
    DlogLDmu <- .DlogLDmu_trunc_nb2
    D2logLDmu2 <- .D2logLDmu2_trunc_nb2
    D3logLDmu3 <- .D3logLDmu3_trunc_nb2
    environment(DlogLDmu) <- environment(D2logLDmu2) <- environment(D3logLDmu3) <- environment(aic) 
    dHobs_trunc_aux_funs <- list(D2muDeta2=D2muDeta2,D3muDeta3=D3muDeta3,DlogLDmu=DlogLDmu,D2logLDmu2=D2logLDmu2,D3logLDmu3=D3logLDmu3)
  } else dHobs_trunc_aux_funs <- NULL
  ##########################################################
  #
  ## all closures defined here have parent.env the environment(spaMM_Gamma) ie <environment: namespace:spaMM>
  ## changes the parent.env of all these functions (aic, dev.resids, simfun, validmu, variance): 
  # parent.env(environment(aic)) <- environment(stats::binomial) ## parent = <environment: namespace:stats>
  structure(c(
    list(family = structure("negbin",
                            withArgs=quote(paste0("Neg.binomial(shape=",signif(shape,4),")"))), 
         link = linktemp, linkfun = linkfun, 
         linkinv = linkinv, variance = variance, dev.resids = dev.resids, 
         aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
         validmu = validmu, valideta = stats$valideta, simulate = simfun, 
         Dtheta.Dmu=Dtheta.Dmu, D2theta.Dmu2=D2theta.Dmu2,
         dlW_Hexp__detafun=dlW_Hexp__detafun, coef1fun=coef1fun, d_dlWdmu_detafun=d_dlWdmu_detafun,
         flags=list(obs=TRUE, exp=TRUE, canonicalLink=FALSE),
         zero_truncated=(trunc==0L)),
    dHobs_trunc_aux_funs), 
    class = "family")
}

Tnegbin <- function(shape = stop("negbin's 'shape' must be specified"), link = "log") {
  mc <- match.call()   # avoid evaluation of promise...
  mc[[1L]] <- get("negbin", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  mc[["trunc"]] <- 0L
  eval(mc, parent.frame())
}

