

.DlogLDmu_trunc_nb2 <- local({
  shape <- NaN # The function-definition environment will be replaced by the family env where this def is used.
  function(y, mu, wt, n) { # dlogL/dmu
    term <- shape*(y-mu)/(mu*(mu + shape))
    term <- drop(term)
    p0 <- .negbin2_p0(mu,shape)
    Mdlog1mp0 <- - p0^(1+1/shape) /(1-p0)
    term + Mdlog1mp0
  }
})

.D2logLDmu2_trunc_nb2 <- local({
  shape <- NaN
  function(y, mu, wt, n) { 
    term <-  shape *(mu^2 - y* (shape+2*mu))/(mu*(shape+mu))^2 
    term <- drop(term)
    p0 <- .negbin2_p0(mu,shape)
    Md2log1mp0 <- shape* p0 * (1+shape-p0)/ ((shape+mu)*(1-p0))^2
    term + Md2log1mp0
  }
})

.D3logLDmu3_trunc_nb2 <- local({
  shape <- NaN
  function(y, mu, wt, n) { # element of computation of D3logcLdeta3 for d logdet Hessian
    term <- ( 2*shape*(-mu^3 + y *(shape^2 + 3*shape*mu + 3*mu^2)) ) / (mu^3 *(shape + mu)^3)
    term <- drop(term)
    p0 <- .negbin2_p0(mu,shape)
    Md3log1mp0 <- -shape*p0*(3*shape*(1 - p0) + 2*(1 - p0)^2 + 
                               shape^2*(1 + p0))/((shape + mu)*(1 - p0))^3
    term + Md3log1mp0
  }
})



# negbin2: version of negbin() with extra functionalities "for devel purposes" -- actually documented
# SHOULD at least be kept in the Zvariants for devel purpose
negbin2 <- function (shape = stop("negbin2's 'shape' must be specified"), link = "log", trunc=-1L, LLF_only=TRUE) {
  mc <- match.call()
  if (inherits(shch <- substitute(shape),"character") ||
      (inherits(shch,"name") && inherits(shape, "function")) # "name" is for e.g. negbin(log)
      # (but testing only "name" would catch e.g. negbin(shape=shape) )
  ) { 
    if (inherits(shch,"character")) shch <- paste0('"',shch,'"')
    errmess <- paste0('It looks like negbin2(',shch,') was called, which absurdly means negbin2(shape=',shch,
                      ').\n  Use named argument: negbin2(link=',shch,') instead.')
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
  variance <- function(mu) mu + mu^2 /shape
  validmu <- function(mu) all(mu > 0)
  # logl and derived functions  derived from negbin[2] b replacing shape by shape*mu
  # - negbin1(shape=3)$aic(y,NA,mu,wt=1)/2  = dnbinom(x=y,size=shape*mu,mu=mu,log=TRUE)
  # compared to
  # - negbin(shape=3)$aic(y,NA,mu,wt=1)/2  = dnbinom(x=y,size=shape,mu=mu,log=TRUE)
  if (trunc==0L) { # if truncated
    aic <- function(y, n, mu, wt, dev) { ## not really the aic... -2 clik
      - 2 * sum(logl(y,mu,wt))
    }
    logl <- function(y, mu, wt) {
      mlogl <- (y + shape) * log(mu + shape) - y * log(mu) + 
        lgamma(y + 1) - shape * log(shape) + lgamma(shape) - lgamma(shape + y) +
        log(1-.negbin2_p0(mu,shape)) 
      mlogl[y==1L & mu==0] <- 0 
      drop(- mlogl * wt)
    }
    DlogLDmu <- .DlogLDmu_trunc_nb2
    D2logLDmu2 <- .D2logLDmu2_trunc_nb2
    D3logLDmu3 <- .D3logLDmu3_trunc_nb2
    environment(DlogLDmu) <- environment(D2logLDmu2) <- environment(D3logLDmu3) <- environment(aic) 
    
    sat_logL <- function(y, wt) { 
      uniqy <- unique(y)
      uniqmu <- rep(NA, length(uniqy))
      uniqmu[uniqy == 1L] <- 0
      getmu <- function(y) {
        if (y==1L) return(0)
        if (shape >= 1) {
          lower <-  (y-1)*0.999
        } else lower <-  ( y-1 )*shape*0.999
        interval <- c(lower, y) # derivative always <0 in mu=y. For shape<1, perhaps a stricter upper bound could be found
        uniroot(DlogLDmu, interval=interval, y = y)$root
      }
      ygt1 <- (uniqy > 1L)
      uniqmu[ygt1] <- sapply(uniqy[ygt1], getmu)
      uniqlogl <- logl(uniqy,mu=uniqmu,wt=1)
      uniqlogl[match(y, uniqy)]
    } 
    
    ## the dev.resids serves in llm.fit...
    dev.resids <- function(y,mu,wt) { 2*(sat_logL(y, wt=wt)-logl(y,mu=mu,wt=wt)) } # cannot use $aic() which is already a sum...
    
    # dev.resids <- function(y, mu, wt) {
    #   # computation is in Mathematica notebook.   # but this was nonsense ? cannot juste add a term to the untruncated case
    #   # The logic of MolasL10 is weird bc their "definition" of deviance conflicts with a more fundamental def. 
    #   2 * wt * (y * log(pmax(1, y)/mu) - (y + shape) * log((y + shape)/(mu + shape)) +
    #             log( (1-.negbin2_p0(mu,shape)) / (1-.negbin2_p0(y,shape)) )
    #   )
    # }
    
  } else {
    sat_logL <- NULL
    dev.resids <- function(y, mu, wt) {
      2 * wt * (y * log(pmax(1, y)/mu) - (y + shape) * log((y + shape)/(mu + shape)))
    }
    aic <- function(y, n, mu, wt, dev) {
      term <- (y + shape) * log(mu + shape) - y * log(mu) + 
        lgamma(y + 1) - shape * log(shape) + lgamma(shape) - lgamma(shape + y)
      2 * sum(term * wt)
    }
    
    logl <- function(y, mu, wt) {
      term <- (y + shape) * log(mu + shape) - y * log(mu) + 
        lgamma(y + 1) - shape * log(shape) + lgamma(shape) - lgamma(shape + y)
      drop(- term * wt)
    }
    
    ## link-independent function #####################
    DlogLDmu <- function(mu, y, wt, n) drop(shape*(y-mu)/(mu*(mu + shape)))
    
    D2logLDmu2 <- function(mu, y, wt, n) { # element of computation of Hobs weights, (-) D2logcLdeta2
      # (\[Theta] (\[Mu]^2 - y (\[Theta] + 2 \[Mu])))/(\[Mu]^2 (\[Theta] + \[Mu])^2)
      drop(shape *(mu^2 - y* (shape+2*mu))/(mu*(shape+mu))^2) 
    }
    
    D3logLDmu3 <- function(mu, y, wt, n) { # element of computation of D3logcLdeta3 for d logdet Hessian
      # (2 \[Theta] (-\[Mu]^3 + y (\[Theta]^2 + 3 \[Theta] \[Mu] + 3 \[Mu]^2)))/(\[Mu]^3 (\[Theta] + \[Mu])^3)
      drop( 2*shape*(-mu^3 + y *(shape^2 + 3*shape*mu + 3*mu^2)) ) / (mu^3 *(shape + mu)^3)
    }
    ###################################################
  }
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
      p0 <- .negbin2_p0(mu_U,shape) # cf TruncatedNegBin.nb
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
    .rnbinom(n=nsim * length(ftd), 
             size=shape, 
             mu_str=ftd, zero_truncated=zero_truncated)
  }
  ######################## For GLM-based algos:
  Dtheta.Dmu <- function(mu) 1/(mu*(1+mu/shape))
  D2theta.Dmu2 <- function(mu) -(1+2*mu/shape)/(mu*(1+mu/shape))^2
  if (trunc == 0L) { # If truncated.
    ## No ad-hoc functions or GLM-specific Hobs code, evaluating derivatives of the Hexp weights
    #
    ##, The mandatory function for LLM-capable Hobs are the D<n>logLDmu<n> functions,
    ## used to compute $Md2logcLdeta2 and $Md3logcLdeta3 in .muetafn(), and then generic code:
    ## dlW_deta <- muetablob$Md3logcLdeta3/muetablob$Md2logcLdeta2 
    ## if (calcCoef1) coef1 <- dlW_deta/muetablob$Md2logcLdeta2 
    dlW_Hexp__detafun <- d_dlWdmu_detafun <- coef1fun <- NULL
  } else {
    # ad-hoc functions or GLM-specific Hobs code, evaluating derivatives of the Hexp weights
    dlW_Hexp__detafun <- switch( # where  W is that of Hexp
      linktemp,
      "log"= function(mu) shape/(mu+shape),
      "identity"= function(mu) -(2*mu + shape)/(mu^2 + mu*shape),
      "sqrt"= function(mu) -2 * sqrt(mu)/(mu + shape)
    )
    d_dlWdmu_detafun <- switch( # see Hessian_weights.nb # where  W is that of Hexp
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
    #DvarDmu <- function(mu)1+ 2*mu/shape
  } 
  ########################
  D2muDeta2 <- .D2muDeta2(linktemp)
  D3muDeta3 <- .D3muDeta3(linktemp)
  ## all closures defined here have parent.env the environment(spaMM_Gamma) ie <environment: namespace:spaMM>
  ## changes the parent.env of all these functions (aic, dev.resids, simfun, validmu, variance): 
  # parent.env(environment(aic)) <- environment(stats::binomial) ## parent = <environment: namespace:stats>
  structure(list(family = structure("negbin",
                                    withArgs=quote(paste0("negbin2(shape=",signif(shape,4),")"))), 
                 link = linktemp, linkfun = linkfun, 
                 linkinv = linkinv, variance = variance, dev.resids = dev.resids, sat_logL=sat_logL, logl=logl, dev.resids=dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun, 
                 DlogLDmu = DlogLDmu, D2logLDmu2 = D2logLDmu2, D3logLDmu3 = D3logLDmu3, D2muDeta2 = D2muDeta2, D3muDeta3 = D3muDeta3,
                 # functions not to be used by the LLF_only methods:
                 Dtheta.Dmu=Dtheta.Dmu, D2theta.Dmu2=D2theta.Dmu2,
                 dlW_Hexp__detafun=dlW_Hexp__detafun, coef1fun=coef1fun, d_dlWdmu_detafun=d_dlWdmu_detafun,
                 #
                 flags=list(obs=TRUE, 
                            exp= FALSE, # so that negbin2's obsInfo is always be TRUE
                            ## Currently negbin2 does not provide Hexp fits (setting exp=TRUE => obsInfo to FALSE gives result distinct from negbin())
                            ## despite all the necessary GLM family components being present. 
                            ## This must be because other spaMM procedures do not handle negbin2 appropriately for that case.
                            ## This does not matter at user level, because THIS negbin2 has no reason to belong to the API.
                            ## The API negbin() family has all capacities including truncated Hobs. 
                            canonicalLink=FALSE, GLMalgo= ! LLF_only),
                 zero_truncated=(trunc==0L)), 
            class = c("family","LLF"))
}
