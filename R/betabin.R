betabin <- function (prec = stop("betabin's 'prec' must be specified"), link = "logit") {
  trunc <- FALSE
  .BinomialDen <- NULL # hack so that non need for special handling of BinomialDen in llm.fit()
  resid.model <- list2env(list(off=0)) # env so that when we assign to it we don't create a new instance of the family object
  mc <- match.call()
  if (inherits(shch <- substitute(prec),"character") ||
      (inherits(shch,"name") && inherits(prec, "function")) # "name" is for e.g. beta_resp(logit)
      # (but testing only "name" would catch e.g. negbin(prec=prec) )
  ) { 
    if (inherits(shch,"character")) shch <- paste0('"',shch,'"')
    errmess <- paste0('It looks like betabin(',shch,') was called, which absurdly means betabin(prec=',shch,
                      ').\n  Use named argument: betabin(link=',shch,') instead.')
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
  
  # not needed for LLM fit algo but needed for residVar computation  
  variance <- function(mu, new_fampar=NULL, sizes) {
    
    if (missing(sizes)) return(mu*(1-mu)) #residVar computation in particular
    
    if ( ! is.null(new_fampar)) prec <- new_fampar
    mu *(1- mu)* sizes*(sizes+prec)/(1+prec) 
  }
  validmu <- function(mu) all(mu > 0 & mu < 1)
  
  # in the following functions I must distingusih the wt= prior weights for the prec model, and the BinomialDen; by contrast wt is used for BinomialDen in the binomial code
  logl <- function(y, n=.BinomialDen, mu, wt) { # mu is muFREQS is alpha/(alpha+beta) = alpha/prec so mu*prec=alpha
    y <- drop(y)
    prec <- prec*c(wt) 
    lchoose(n,y)  +  lbeta(y+mu*prec, n-y+(1-mu)*prec)  -  lbeta(mu*prec, (1-mu)*prec)
  }
  # phi not to be used: see .add_Md_logcLdeta_terms() 
  DlogLDmu <- function(muFREQS, y, BinomialDen=.BinomialDen, wt, prec_it=NULL) {
    # dlogL/dmu
    y <- drop(y)
    if ( ! is.null(prec_it)) {
      prec <- prec_it*wt
    } else prec <- prec*c(wt) 
    # Binomial: (y-muCOUNT)/(muFREQS*(1-muFREQS)) 
    muCOUNT <- muFREQS * BinomialDen
    #cat(crayon::yellow("dL"))
    #str((y-muCOUNT)/(muFREQS*(1-muFREQS)))
    ## \[Phi] (-PolyGamma[0, \[Mu] \[Phi]] + PolyGamma[0, \[Phi] - \[Mu] \[Phi]] - PolyGamma[0, n - y + \[Phi] - \[Mu] \[Phi]] + PolyGamma[0, y + \[Mu] \[Phi]])
    #str(prec*(digamma(y+prec * muFREQS) - digamma(BinomialDen-y+prec * (1-muFREQS)) - digamma(prec * muFREQS) + digamma(prec * (1-muFREQS))))
    prec*(digamma(y+prec * muFREQS) - digamma(BinomialDen-y+prec * (1-muFREQS)) - digamma(prec * muFREQS) + digamma(prec * (1-muFREQS)))
  }
  D2logLDmu2 <- function(muFREQS, y, BinomialDen=.BinomialDen, wt) {
    y <- drop(y)
    prec <- prec*c(wt) 
    # Binomial:  (-BinomialDen + y)/(1 - muFREQS)^2  -  y/muFREQS^2 
    #cat(crayon::red("d2L"))
    #str((-BinomialDen + y)/(1 - muFREQS)^2  -  y/muFREQS^2 )
    #str(prec^2 *(trigamma(y+prec * muFREQS) + trigamma(BinomialDen-y+prec * (1-muFREQS)) -trigamma(prec*muFREQS)-trigamma(prec * (1-muFREQS))))
    prec^2 *(trigamma(y+prec * muFREQS) + trigamma(BinomialDen-y+prec * (1-muFREQS)) -trigamma(prec*muFREQS)-trigamma(prec * (1-muFREQS)))
  }
  D3logLDmu3 <- function(muFREQS, y, BinomialDen=.BinomialDen, wt) {
    y <- drop(y)
    prec <- prec*c(wt) 
    prec^3 *(psigamma(y+prec * muFREQS, deriv=2) - psigamma(BinomialDen-y+prec * (1-muFREQS), deriv=2) - psigamma(prec * muFREQS, deriv=2) + psigamma(prec * (1-muFREQS), deriv=2))
  }
  
  get_sat_mu <- function(y, wt, BinomialDen=.BinomialDen, prec_it=NULL) {
    if (y==0L) return(1e-12)
    if (y==BinomialDen) return(1-1e-12)
    freq <- y/BinomialDen
    if (freq>0.5) {
      interval <- c(0.5, freq)
    } else interval <- c(freq, 0.500000001) # ... handling y=1/2 exactly 
    uniroot(DlogLDmu, interval=interval, y = y, wt=wt, BinomialDen=BinomialDen, prec_it=prec_it)$root
  }
  
  sat_logL <- function(y, BinomialDen=.BinomialDen, wt) { # cf call through dev.resids: wt is typically rep(1,length(fv)) here
    nr <- length(y)
    muFREQS <- numeric(nr)
    if (length(prec)>1L) {
      for (it in seq_len(nr)) muFREQS[it] <- get_sat_mu(y[it], BinomialDen=BinomialDen[it], wt=wt[it], prec_it=prec[it])
    } else for (it in seq_len(nr)) muFREQS[it] <- get_sat_mu(y[it], BinomialDen=BinomialDen[it], wt=wt[it])
    logl(y=y,n=BinomialDen, mu=muFREQS, wt=wt)
  } 
  
  dev.resids <- function(y,mu,BinomialDen=.BinomialDen, wt) { # input (y, mu=fv, BinomialDen=BinomialDen, ,wt="argh"), cf .dev_resids() extractor; the mu are muFREQS
    2*(sat_logL(y, BinomialDen=BinomialDen, wt=wt)-logl(y,n=BinomialDen, mu=mu,wt=wt)) } #
  
  aic <- function(y, n=.BinomialDen, mu, wt, ...) { # 'wt' for the prior weights on the precision parameter
    - 2 * sum(logl(y, n, mu, wt))
  }
  # binomial()$initialize:
  initialize <- expression({
    if (NCOL(y) == 1L) {
      if (is.factor(y)) 
        y <- y != levels(y)[1L]
      .BinomialDen <- rep.int(1L, nobs)
      y[weights == 0] <- 0L
      if (any(y < 0L | y > 1L)) stop("y values must be 0 <= y <= 1")
      mustart <- (weights * y + 0.5)/(weights + 1)
    } else if (NCOL(y) == 2L) {
      if (any(abs(y - round(y)) > 0.001)) 
        warning("non-integer counts in a betabinomial model!", domain = NA)
      .BinomialDen <- y[, 2L] + (y <- y[, 1L]) # different handling of y relative do stats::binomial
      if (any(n0 <- .BinomialDen == 0L)) y[n0] <- 0
      mustart <- (y + 0.5)/(.BinomialDen + 1)
    } else stop("for the 'betabin' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures", 
                domain = NA)
    environment(family$aic)$.BinomialDen <- .BinomialDen
  })
  simulate <- function(object, nsim) { # This function is not used by spaMM. See .r_resid_var() instead.
                                     # Calling simulate on <glm() fit using a spaMM family> will call them.
                                     # Cf example in test.LLM.   
    ftd <- fitted(object)
    # prec <- .get_family_par(family=object$family) ## fails as the spaMM namespace is not accessible form here... and not needed anyway. 
    precW <- prec*object$prior.weights
    ntot <- nsim * length(ftd)
    rmu <- rbeta(n=ntot, shape1=ftd*precW,shape2=(1-ftd)*precW)
    rbinom(n=ntot, size=object$BinomialDen, prob=rmu)
  }
  ## all closures defined here have parent.env the environment(beta_resp) ie <environment: namespace:spaMM>
  ## changes the parent.env of all these functions (aic, dev.resids, simulate, validmu, variance): 
  D2muDeta2 <- .D2muDeta2(linktemp)
  D3muDeta3 <- .D3muDeta3(linktemp)
  parent.env(environment(aic)) <- environment(stats::binomial) ## parent = <environment: namespace:stats>; 
  # before the change this is <environment: namespace:spaMM>, necessary to access .D2muDeta2...
  structure(list(family = structure("betabin",
                                    withArgs=quote({
                                      if ( ! is.null(resid.model$beta)) {
                                        paste0("betabin(prec=",paste0(signif(prec[1:min(3L,length(prec))],4),collapse=" "),"...)")
                                      } else paste0("betabin(prec=",signif(prec,4),")")
                                    })), 
                 link = linktemp, linkfun = stats$linkfun, 
                 linkinv = stats$linkinv, variance = variance,  sat_logL=sat_logL, logl=logl, dev.resids=dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simulate, 
                 DlogLDmu = DlogLDmu, D2logLDmu2 = D2logLDmu2, D3logLDmu3 = D3logLDmu3, 
                 D2muDeta2 = D2muDeta2, D3muDeta3 = D3muDeta3, resid.model=resid.model, 
                 flags=list(obs=TRUE, exp=FALSE, canonicalLink=FALSE, LLgeneric=TRUE)), 
            class = c("LLF","family"))
}
