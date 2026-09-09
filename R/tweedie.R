# jmax_sc is intended to control the number of terms in summations
# There is j_max = arma::pow(yg, 2.0 - r) / (phig * (2.0 - r)) ...
# => reduce phi to increase jmax, and use Tweedie rescaling identity (12.3) in DunnS18.
.get_dtweedie <-  local({
  latest_jmax_sc <- -1L
  dtweedie <- NULL
  function(jmax_sc=1L, ...) {
    if (jmax_sc != latest_jmax_sc) {
      dtweedie <<- NULL
      latest_jmax_sc <<- jmax_sc
    }
    
    if (is.null(dtweedie)) {
      .dtweedie <- get0("dtweedie", asNamespace("tweedieDistr"), inherits=FALSE) 
      if (is.null(.dtweedie)) {
        warn <- paste0( "If the 'tweedieDistr' package were installed,\n", 
                        "fitting tweedie models might be faster.")
        .dtweedie <- get0("dtweedie", asNamespace("tweedie"), inherits=FALSE) 
        if (is.null(.dtweedie)) stop("No backend package is available to fit tweedie models using spaMM.")
        # Case using .dtweedie = tweedie::dtweedie:
        dtweedie <<- function(y, mu, phi, p, log) { 
          if (jmax_sc==1L) {
            res <- .dtweedie(y, mu=mu, phi=phi, power = p) 
            if (log) res <- log(res)
          } else {
            res <- numeric(length(mu))
            isy0 <- y==0
            if (any(isy0)) {
              if (length(phi)>1L) {
                phi0 <- phi[isy0]
                phin0 <- phi[ ! isy0]
              } else phi0 <- phin0 <- phi
              if (log) {
                res[isy0] <- log(.dtweedie(0, mu=mu[isy0], phi=phi0, power = p)) 
              } else  res[isy0] <- .dtweedie(0, mu=mu[isy0], phi=phi0, power = p) 
            } else phin0 <- phi
            sc <- jmax_sc^(1/(p-2))
            sc_d <- .dtweedie(sc*y[ ! isy0], mu=sc*mu[ ! isy0], phi=phin0/jmax_sc, power = p) 
            if (log) {
              res[! isy0] <- log(sc_d) +log(sc)
            } else res[! isy0] <- sc_d*sc
          }
          res
        }
      } else { # Case using .dtweedie = tweedieDistr::dtweedie:
        dtweedie <<- function(y, mu, phi, p, log) { 
          if (jmax_sc==1L) {
            res <- .dtweedie(x=y, mean = mu, dispersion = phi, power = p, log = log)
          } else {
            res <- numeric(length(mu))
            isy0 <- y==0
            if (any(isy0)) {
              if (length(phi)>1L) {
                phi0 <- phi[isy0]
                phin0 <- phi[ ! isy0]
              } else phi0 <- phin0 <- phi
              res[isy0] <- .dtweedie(x=0, mean = mu[isy0], dispersion = phi0, power = p, log = log)
            } else phin0 <- phi
            sc <- jmax_sc^(1/(p-2))
            sc_d <- .dtweedie(x=sc*y[ ! isy0], mean = sc*mu[ ! isy0],
                              dispersion = phin0/jmax_sc, power = p, log = log)
            if (log) {
              res[! isy0] <- sc_d +log(sc) 
            } else res[! isy0] <- sc_d*sc
          }
          res
        } 
      }
    }
    dtweedie
  }
})

.get_ptweedie <-  local({
  ptweedie <- NULL
  function() {
    if (is.null(ptweedie)) {
      .ptweedie <- get0("ptweedie", asNamespace("tweedieDistr"), inherits=FALSE) 
      if (is.null(.ptweedie)) {
        .ptweedie <- get0("ptweedie", asNamespace("tweedie"), inherits=FALSE) 
        if (is.null(.ptweedie)) stop("No backend package is available to fit tweedie models using spaMM.")
        ptweedie <<- function(q, mu, phi, p, lower.tail=TRUE, log.p=FALSE) {
          res <- .ptweedie(q, mu=mu, phi=phi, power = p, verbose = FALSE) 
          if ( ! lower.tail) res <- 1-res
          if (log.p) res <- log(res)
          res
        }
      } else {
        ptweedie <<- function(q, mu, phi, p, lower.tail=TRUE, log.p=FALSE) {
          .ptweedie(q, mean = mu, dispersion = phi, power = p, lower.tail = lower.tail, 
                    log.p = log.p)
        } 
      }
    }
    ptweedie
  }
})

# A better version would have dtweedie_dlogfdphi in C++:
# It would be faster, would remove dependence on tweedie::, and could be protected against NaNs
.get_dlogLdphi_dldphi <-  function(dtweedie) { # this dtweedie should be the family member version
                                               # with arguments y, mu, phi, p, log
  tweedie_dlogfdphi <- get0("dtweedie_dlogfdphi", # 'f' for a fn that is not 'L' !!!!
                            asNamespace("tweedie"), inherits=FALSE) 
  dlogLdphi <- function(y, mu, p, phi) { # assuming 1<p<2
    if (y==0) return(mu^(2-p)/(phi^2 * (2-p)))
    
    k <- phi^(1/(p - 2)) 
    cond <- 0 < k & k < 1 # When 1<p<2, k is >1 if phi<1
    k_is_lo <- which(cond)  # so this is phi> 1
    
    res <- numeric(length(mu))
    if (any(k_is_lo)) {
      klo <- k[k_is_lo]
      ylo <-  y[k_is_lo]
      mulo <- mu[k_is_lo]
      kylo <- klo*ylo
      f <- dtweedie(y = kylo, p = p, mu = klo * mulo, phi = 1, log=FALSE)
      d <- tweedie_dlogfdphi(y = kylo, p = p, 
                             mu = klo * mulo, phi = 1)
      top <- d * f
      res[k_is_lo] <- top/f * klo^(2 - p)
    }
    
    k_is_hi <- which( ! cond) # phi < 1
    if (any (k_is_hi)) {
      res[k_is_hi] <- tweedie_dlogfdphi(y = y[k_is_hi], power = p, 
                                        mu = mu[k_is_hi], phi = phi[k_is_hi])
    }
    # if (is.nan(res)) browser()
    res
  }
  dlogLdphi
}

# This requires only one of the two packages 
# and is rel fast if dtweedie uses C++ code of tweedieDistr
.get_dlogLdphi_numDeriv <- function(dtweedie, # this dtweedie should be the family member version,
                                    # with arguments y, mu, phi, p, log
                                    method="Richardson", method.args=list(eps=1e-4),
                                    ...) {   
  force(dtweedie)
  force(method)
  force(method.args)
  function(y, mu, p, phi) {
    if (y==0) return(mu^(2-p)/(phi^2 * (2-p)))
    
    objfn <- function(disp) dtweedie(y=y, mu=mu, phi=disp, p=p, log=TRUE)
    if (phi<1.01*method.args$eps) {side <- 1} else side <- NA
    numDeriv::grad(objfn, x=phi, method=method, method.args = method.args, side=side)
  } 
}



# $dlogfdphi easily returns NaN, e.g family$dlogfdphi(y=7.988925, mu=7.834178, power=1.08, phi=0.04418144)
# where it evaluates as the ratio of two infinite values.
# dtweedie_dldphi acts as a wrapper for -2*sum(dtweedie_dlogfdphi()),
#   and does not vectorize over phi.

# code can be checked by comparing cas with p=1.99 to Gamma(log) case 
.calc_levphi_corr_tweedie <- function(y, mu, phiscaled, family, dlogLs) {
  p <- environment(family$aic)$"p"
  nobs <- length(y)
  qcorr <- numeric(nobs)
  isy0 <- y==0
  mu_0y <- mu[isy0]
  phisc_0y <- phiscaled[isy0]
  qcorr[isy0] <- 1 # (the code for !isy0 also gives this result)
  
  for (ii in which(! isy0)) {
    yi <- y[ii]
    mui <- mu[ii]
    phisci <- phiscaled[ii]
    unwei_devres <- family$dev.resids(y=yi,mu=mui, wt=1)
    dlogLdphi <- dlogLs[ii]
    if (is.nan(dlogLdphi) || is.infinite(dlogLdphi)) {
      delta <- 1.0e-5
      a1 <- family$dtweedie(p = p, 
                              phi = phisci, 
                              mu = mui, 
                              y = yi,
                            log=FALSE)
      a2 <- family$dtweedie(p = p, 
                              phi = phisci+delta, 
                              mu = mui, 
                              y = yi,
                              log=FALSE)
      if (a1==0 && a2==0) {
        dlogLdphi <- 0 # quick patch
      } else if (a2==0) {
        a2 <- family$dtweedie(p = p, 
                                phi = phisci-delta, 
                                mu = mui, 
                                y = yi,
                                log=FALSE)
        if (a2==0) {
          dlogLdphi <- 0 # quick patch
        } else dlogLdphi <- (log(a1) - log(a2) ) / delta
      } else dlogLdphi <- (log(a2) - log(a1) ) / delta
    }
    # if (is.infinite(dlogLdphi)) browser()
    qcorr[ii] <- 1 + (2*phisci^2* dlogLdphi - unwei_devres)/phisci
  }
  qcorr
}


.tidy_Tw_index <- function(index) {
  if ( ! inherits(index, c("character","numeric"))) index <- eval(index, parent.frame()) 
  if (inherits(index,"character")) { # eg tweedie("Gamma")
    m <- match(index, c("gaussian", "poisson", "Gamma", "gamma", 
                        "inverse.gaussian"))
    if (is.na(m)) 
      stop("Tw_index should be a number")
    else {
      index <- c(0, 1, 2, 2, 3)[m]
    }
  }
  index     
}

.deparse_link <- function(link) {
  linktemp <- substitute(link, parent.frame()) # if link was char LHS is char ; else deparse will create a char from a language object 
  if ( ! inherits(linktemp,c("character","numeric"))) {
    linktemp <- deparse(linktemp)
    okLinks <- c("identity", "log", "inverse")
    if (linktemp %in% okLinks) {
      link <- linktemp
    } else link <- eval(link, parent.frame())    
  }
  link
}

.tidy_Tw_link <- function(link) {
  if (link=="fit") {
  } else {
    if (inherits(link,"character")) {
      m <- match(link, c("identity", "log", "inverse"))
      if (is.na(m)) 
        stop("Tw_link should be a number")
      else link <- c(1, 0, -1)[m]
    }
  } 
  link
}

# The tweedieDistr version is pure R.
.rtweedie <- function(n, mu = 1, phi = 1, p = 1.5) {
  lambda <- (mu^(2 - p)) / (phi * (2 - p))
  shape1 <- (2 - p) / (p - 1)
  scale <- phi * (p - 1) * (mu^(p - 1))
  Pdraws <- rpois(n, lambda)
  rgamma(n, shape=Pdraws * shape1, scale=scale)
}


# Smyth02 discusses divergence of phi estimates but see my comments on it:
# maybe not relevant as even the direct solution for a phi intercept 
# may need correction for HLfit inner iterations to converge. 
# Instead, here, algo using info from dlogLdphi  
# but NOT assuming that they are the correct gradients of
# a one dimensional fn: they depend on all other fitted parameters.
#
# Example of case (*): case where we moved left (dx<0),
# and both dlogLs are <<0 (but *not* more negative right:  
# combination dlogLdphi>prev_dlogLdphi and prev_dlogLdphi < 0)
# both x values appear too high but the gradients cannot 
# be used to guess a solution. We only use dy/dx.
#
# In case (**) we avoid phi jumps to extreme values 
# (=> fit be trapped in extreme phi, lambda local minimum).
#
# Note that ..calcPHI has 'iter' argument so it might eventually be used.
.fix_div_phi <- function(next_phi_est, prev_PHIblob, dev_res_info) {
  # keep (x,y) before y is corrected:
  input_phi_est <- prev_PHIblob$next_phi_est
  dlogLdphi <- dev_res_info$dlogLdphi
  next_fix_div_blob <- list(x=input_phi_est,
                            y=next_phi_est, dlogLdphi=dlogLdphi)
  prev_fix_div_blob <- prev_PHIblob$fix_div_blob
  if ( ! is.null(ante_x <- prev_fix_div_blob$x)) {
    ante_y <- prev_fix_div_blob$y
    dx <- input_phi_est - ante_x
    dy <- next_phi_est - ante_y # both uncorrected fitted values from Gamma GLM
    
    prev_dlogLdphi <- prev_fix_div_blob$dlogLdphi
    if (dlogLdphi*prev_dlogLdphi < 0) {
      # Solution presumably inbetween, we use the gradient values to guess the solution
      DdlogL <- dlogLdphi - prev_fix_div_blob$dlogLdphi
      next_phi_est <- input_phi_est - dx*dlogLdphi/DdlogL
    } else if (dlogLdphi>prev_dlogLdphi && prev_dlogLdphi > 0) {
      # do not move left if both grads are positive and more positive left (but can move right)
      if (all(next_phi_est < input_phi_est)) next_phi_est <- input_phi_est
    } else if (dlogLdphi<prev_dlogLdphi && prev_dlogLdphi < 0) {
      # symmetric case, do not move right if...
      if (all(next_phi_est > input_phi_est)) next_phi_est <- input_phi_est
    } else if ((m <- mean(dy/dx)) < - 0.9999999) { 
      # See comments above (*)
      locfac <- 1+(next_phi_est-input_phi_est)/(dx-dy) 
      next_phi_est <- ante_x+locfac*dx  # =  ante_y+locfac*dy  
    }
  } else { # otherwise first iters, so above corrections not feasible. See (**)
    prev_phi_est <- dev_res_info$phi_est 
    next_phi_est <- pmin(pmax(prev_phi_est/1.1, next_phi_est),prev_phi_est*1.1)
  }
  list(fix_div_blob=next_fix_div_blob, next_phi_est=next_phi_est)
}

.get_fix_div_phi <- function(fix_div_phi=NULL, ...) {
  if (is.null(fix_div_phi)) fix_div_phi <- .fix_div_phi
  fix_div_phi
}

# Handled \dots args: 
#  method & method.args for numDeriv::grad() in .get_dlogLdphi_numDeriv();
#  jmax_sc for .get_dtweedie() to provide it in the environment of $dtweedie();
#  fix_div_phi for .get_fix_div_phi() 
#     (allows expert user to provide alternative to .fix_div_phi)
tweedie <- function (index, link = "log", numderiv=TRUE, ...) {
  resid.model <- list2env(list(off=0)) # for outer phiGLM 
  
  if (p_missing <- missing(index)) {
    delayedAssign("p", stop("tweedie's 'index' must be specified"))
  } else p <-  .tidy_Tw_index(index)

  link <- .deparse_link(link)
  Tw_link <-  .tidy_Tw_link(link)
  if (q_missing <- (Tw_link=="fit")) {
    # This case is recognized in two places in other spaMM function by if (get("q_missing", ...))
    # (rather than trying to test q which may be missing, then no longer missing...)
  } else q <- Tw_link
  remove("Tw_link")
  
  if ( ! p_missing) {    
    if (p == 0) {
      validmu <- function(mu) TRUE
    } else if (p>0) {
      validmu <- function(mu) all(mu >= 0)
    } else  validmu <- function(mu) all(mu > 0)
  } else {
    validmu <- function(mu) {
      if (p==0) {
        TRUE
      } else if (p>0) {
        all(mu >= 0)
      } else all(mu > 0)
    }
  }

  if ( ! q_missing) {    
    if (q == 0) {
      linkfun <- log
      linkinv <- mu.eta <- D2muDeta2 <- D3muDeta3 <- .safe_exp
    } else {
      linkfun <- function(mu) mu^q
      linkinv <- function(eta) eta^(1/q)
      mu.eta <- function(eta) (1/q) * eta^(1/q - 1)
      D2muDeta2 <- function(eta) ((1-q)*eta^(1/q-2))/(q^2)
      D3muDeta3 <- function(eta) ((1-q)*(1-2*q)*eta^(1/q-3))/(q^3)
    }
  } else {
    linkfun <- function(mu) {
      if (q==0) {log(mu)} else mu^q
    }
    linkinv <- function(eta) {
      if (q==0) {.safe_exp(eta)} else eta^(1/q)
    }
    mu.eta <- function(eta) {
      if (q==0) {.safe_exp(eta)} else (1/q) * eta^(1/q - 1)
    }
    D2muDeta2 <- function(eta) { # .D2muDeta2("power")
      if (q==0)  {
        .safe_exp(eta)
      } else ((1-q)*eta^(1/q-2))/(q^2)
    } 
    D3muDeta3 <- function(eta) { # .D3muDeta3("power")
      if (q==0)  {
        .safe_exp(eta)
      } else ((1-q)*(1-2*q)*eta^(1/q-3))/(q^3)
    }
  }
  
  # same relationship between valideta and q as between validmu and p
  if ( ! q_missing) {    
    if (q == 0) {
      valideta <- function(mu) TRUE
    } else if (q>0) {
      valideta <- function(mu) all(mu >= 0)
    } else  valideta <- function(mu) all(mu > 0)
  } else {
    valideta <- function(mu) {
      if (q==0) {
        TRUE
      } else if (q>0) {
        all(mu >= 0)
      } else all(mu > 0)
    }
  }
  
  variance <- function(mu) mu^p
  
  dev.resids <- function(y, mu, wt) { # this is notably used as ..calc.PHI()'s dev.res argument.
    if (p == 1) { # dev.resids from poisson()
      r <- mu * wt
      p <- which(y > 0)
      r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
      2 * r
    } else if (p == 2) { # dev.resids from Gamma()
      -2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
    } else if (p == 0) { # gaussian
      wt * ((y - mu)^2)
    } else { # 1 < p < 2 
      if (FALSE) { # version from statmod::tweedie
        y1 <- y + 0.1 * (y == 0) 
        theta1 <- (y1^(1 - p) - mu^(1 - p))/(1 - p)
        kappa <- (y^(2 - p) - mu^(2 - p))/(2 - p)
        2 * wt * (y * theta1 - kappa)
      } else { 
        fac <- numeric(length(mu))
        theta <- (y^(1 - p) - mu^(1 - p))/(1 - p)
        kappa <- (y^(2 - p) - mu^(2 - p))/(2 - p)
        y0sing <- y==0 | is.infinite(theta) # y^(1 - p) singular near 0
        # mu~y => all terms ->0 may become ~ -1e-16 =>bug. Fix without testing:
        fac[ ! y0sing] <- (y * theta - kappa)[ ! y0sing] + 1e-12
        # and assuming that infinite values are those for y close to zero:
        fac[y0sing] <- mu[y0sing]^(2 - p) /(2 - p) # value for y=0
        2 * wt * fac
      }
    }
  }
  
  initialize <- expression({
    n <- rep(1, nobs)
    mustart <- y + 0.1 * (y == 0)
  })
  
  DlogLDmu <- function(mu,y,wt, phi) drop(wt*(y-mu)*mu^(-p))/phi # {dlogLdTh=(y - mu)/phi} * {dThdmu=mu^{-p}}
  D2logLDmu2 <- function(mu,y,wt, phi) {
    mufac <- - p*y*mu^(-p-1)     - (1-p)*mu^(-p) 
    drop( wt*mufac /phi )
  }
  D3logLDmu3 <- function(mu,y,wt, phi) {
    mufac <-  p*(p+1)*y*mu^(-p-2)     + p*(1-p)*mu^(-p-1) 
    drop( wt*mufac /phi )
  }
  
  aic <- function(y, n, mu, wt, dev) { # formal arguments assumed by glm.fit
    n <- sum(wt)
    disp <- dev/n
    logls <- dtweedie(y, mu=mu, phi=disp, p = p, log=TRUE) # vector
    - 2 * sum(logls)
  }

  if ( ! (p_missing || q_missing)) {
    canonicalLink <- abs(diff(c(q,1 - p)))<1e-8
  } else canonicalLink <-  FALSE 

  ptweedie <- .get_ptweedie()
  dtweedie <- .get_dtweedie(...)
  
  .prettify_index <- function() {
    if (p_missing || is.numeric(index)) { 
      paste0("index p=",signif(p,4))
    } else index
  }
  
  .prettify_link <- function() {
    if (q_missing || is.numeric(link) || link=="fit") { 
      paste0("link=",signif(q,4))
    } else paste0("link=",link)
  }
  
  res <- list(family = structure("tweedie",
                          withArgs=quote(paste0("tweedie(",.prettify_index(),
                                                ", ",.prettify_link(),")"))), 
       variance = variance, dev.resids = dev.resids, 
       aic = aic, link ="power", linkfun = linkfun, linkinv = linkinv, 
       mu.eta = mu.eta, initialize = initialize, validmu = validmu, 
       valideta = valideta,
       DlogLDmu=DlogLDmu, D2logLDmu2=D2logLDmu2, D3logLDmu3=D3logLDmu3,  # D3 affects the result even for fixed phi, lambda
       D2muDeta2=D2muDeta2, D3muDeta3=D3muDeta3, # D3 affects the result even for fixed phi, lambda
       dtweedie=dtweedie, ptweedie=ptweedie,
       flags=list(obs=TRUE, # as all spaMM families, allows full Laplace
                  exp=TRUE, # It is a GLM family so 'Hexp' approx is defined.
                  LLgeneric=TRUE, # all spaMM families have it TRUE by default
                  canonicalLink=canonicalLink)
  )
  
  has_tweedie_pack <- requireNamespace("tweedie", quietly = TRUE) 
  if (numderiv) { # default
    res$dlogLdphi <- .get_dlogLdphi_numDeriv(dtweedie, ...)
    if (has_tweedie_pack) { # include the fn for devel purposes
      res$dlogLdphi_dldphi <- .get_dlogLdphi_dldphi(dtweedie) 
    } else res$dlogLdphi_dldphi <- "tweedie package not available for dlogLdphi_dldphi."
  } else if (has_tweedie_pack) {
    res$dlogLdphi <- .get_dlogLdphi_dldphi(dtweedie)
    res$dlogLdphi_numDeriv <- .get_dlogLdphi_numDeriv(dtweedie, ...) # include the fn for devel purposes
  } else { # numderiv=FALSE but tweedie pack not available
    warning("numderiv=FALSE overriden because tweedie package is not available.")
    res$dlogLdphi <- .get_dlogLdphi_numDeriv(dtweedie, ...)
    res$dlogLdphi_dldphi <- "tweedie package not available for dlogLdphi_dldphi."
  }
  
  res$fix_div_phi <- .get_fix_div_phi(...)

  class(res) <- "family"
  res
}
