.CMP_use_asympto <- function(pow_lam_nu, nu, lambda) {
  ## Gaunt et al suggest lambda>1.5 and pow_lam_nu>1.5
  ## this ensures that nu*pow_lam_nu>1.1 but the expansion is clearly better for nu*pow_lam_nu >> 1 
  return(.spaMM.data$options$CMP_asympto_cond(pow_lam_nu, nu, lambda))
}

.COMP_maxn <- function(lambda,nu) { ## this uses the simplest asympto approx to determine the discrete sum bound 
  # => This is so often called as to become quite time consuming... I rewrote it notably to avoid .COMP_Pr_moments() calls 
  # but without changing the results except at (*)
  pow_lam_nu  <- lambda^(1/nu)
  safe_maxn <- .spaMM.data$options$COMP_maxn
  nu_pow_lam_nu <- nu*pow_lam_nu
  c1 <- (nu^2-1)/24
  c1_nu_pow2 <- c1/(nu_pow_lam_nu^2)
  c1_app <- pow_lam_nu*( 1 -(nu-1)/(2*nu_pow_lam_nu) - c1_nu_pow2- c1_nu_pow2/nu_pow_lam_nu) # asympto approx in .COMP_Pr_moments()
  if (nu<1) { # <=> c1_app > pow_lam_nu
    # most pbatic case is nu->0, lambda->1^- where E -> lambda/(1-lambda) ->Inf while the c1_app approx <1 or <<1
    app_Vn <- pow_lam_nu/nu
    if (c1_app < safe_maxn) {
      res <- 4+c1_app+6*sqrt(app_Vn)  ## real rather than integer, to allow continuity correction 
      # VAGUELY inspired from 6*sqrt(app_Vn) to keep all but 1e-6 of the distrib (Feller p.245). 
      # But in fact this uses Gaussian approx and fails for e.g. y=0, mu=0.09568245 => visible diff between dpois and .dCOMP
      # => add 4. 
      # 1-ppois(4+1+6*sqrt(1),lambda = 1)=8.316108e-10 ; 1-ppois(4+10000+6*sqrt(10000),lambda = 10000)=1.067214e-09
      #   For same <small-value ~ 1e-N >= 1-pgeom(maxn,prob = 1-<lambda>) one needs maxn=-N/log10(lambda) 
      # here N ~ 9 
      # 1-pgeom(-9/log10(0.01),prob = 1-0.01)= 1e-10, 1-pgeom(-9/log10(0.99),prob = 1-0.99)=9.994734e-10
      # and with log(lambda)~lambda-1 we need N*sqrt(app_Vn)~ N/[log(10)*(1-lambda)] To ensure continuity ( & with log10() argument always <1):
      # N ~ 
      if (nu-1< -1e-12 && lambda<1) { # fatally, testing more simply nu<1 failed to prevent a ~ (-1e-16)/log10(1) calculation...  
        res <- -((1-nu)^2)*9/log10(nu+(1-nu)*lambda)+(1-(1-nu)^2)*res  ## real, not integer, to allow continuity correction 
        # => in nu->1 the correction is zero for all 1>lambda>0, in nu->0 res is -9/log10(lambda)
      }
    } else {  # the initial 'c1_app' may be way off, it is Inf if lambda -> 0, nu small, and second moment is NaN
      ## Now using a quick Chebyshev-like approx to recover "at least" ~0.999999 of the distribution: 
      # res <- app_En + 1e3*sqrt(app_Vn) # so that tail beyond upper has proba ~ var/guess^2 =1e-6 # but .COMP_maxn() should produce appropriate guesses
      ## I've had some sort of continuity and derivability concern here and the actual code meant
      # res <- 4+pow_lam_nu+(6+994*(1-1/sqrt(1+(c1_app/safe_maxn-1)^2)))*sqrt(app_Vn)  ## real, to allow continuity correction
      ## which is perhaps not continuous in the way I first thought about it. Continuity with the c1_app < safe_maxn case requires
      res <- 4+c1_app+(6+994*(1-1/sqrt(1+(c1_app/safe_maxn-1)^2)))*sqrt(app_Vn) # (*) ## real, to allow continuity correction
      # which remains > safe_maxn so the returned value may be structure(safe_maxn, upper=res) 
      # maintaining some info in the attribute, but no form of derivability.
    }  
  } else { # mu >=1  => app_Vn will be >= pow_lam_nu/nu
    app_En <- pow_lam_nu
    if (c1_app < safe_maxn) {
      app_Vn <- (pow_lam_nu/nu)*( 1 + c1_nu_pow2 + 2*c1_nu_pow2/nu_pow_lam_nu) # asympto approx in .COMP_Pr_moments()
      res <- 4+app_En+6*sqrt(app_Vn)  ## real rather than integer, to allow continuity correction 
    } else {  # the initial res may be way off, that first 'res' is Inf if lambda -> 0, nu small, and second moment is NaN
      ## Same idea as in the nu<1 case, but now one may have c1_app < safe_maxn <- aap
      app_Vn <- pow_lam_nu/nu
      res <- 4+c1_app+(6+994*(1-1/sqrt(1+(c1_app/safe_maxn-1)^2)))*sqrt(app_Vn)  ## again a real to allow continuity correction 
    }  
  }
  
  if (nu-1< -1e-12 && lambda<1) { # fatally, testing more simply nu<1 failed to prevent a ~ (-1e-16)/log10(1) calculation...  
    res <- -((1-nu)^2)*9/log10(nu+(1-nu)*lambda)+(1-(1-nu)^2)*res  ## real, not integer, to allow continuity correction 
    # => in nu->1 the correction is zero for all 1>lambda>0, in nu->0 res is -9/log10(lambda)
  }
  # this may still be far off if initial app_En, app_Vn were far too great, which occurs 
  #  as the asympto approx is not valid (and diverges) for low lambda and low nu. So we add another simple correction: 
  if (lambda<0.98) res <- min(res, -9/log10(lambda)) # always <= than the value for the geometric
  
  if (res>safe_maxn) {
    attr(safe_maxn,"upper") <- res
    res <- safe_maxn
    if ( ! identical(.spaMM.data$options$COMP_maxn_warned,TRUE)) {
      warning(paste0("maxn truncated to ",res," for (lambda,nu)= (",lambda,",",nu,") and possibly other values."))
      .spaMM.data$options$COMP_maxn_warned <- TRUE
      ### This is first called through 
      # .get_inits_by_xLM(processed, reset = quote(family$family %in% c("COMPoisson", "negbin2")))
      # -> .calc_inits_by_xLM(processed)
      # -> .tryCatch_W_E(glm.fit ...)
      ### so the first warning is suppressed, and no further ones will be emitted UNLESS... 
      # unless we reset .spaMM.data$options$COMP_maxn_warned <- FALSE after the .tryCatch_W_E()
      # and as .get_inits_by_xLM() is called repeatedly for COMPoisson, we do so only if a warning was emitted by the glm
      # which means that COMP_maxn_warned was FALSE in entry to .get_inits_by_xLM()
      # which means that COMP_maxn_warned is not set to FALSE by .get_inits_by_xLM() when it was TRUE in entry to it.
      # if (.spaMM.data$options$need_memoise_warning) {
      #   warning("If the memoise package had been installed (prior to loading spaMM), faster computation could be possible.")
      #   .spaMM.data$options$need_memoise_warning <- FALSE
      # }
    }
  }
  res
}

# debug version including R version replaced by Rcpp version
.debug_COMP_Z_integrand <- function(z,lambda,eta=log(lambda),nu,moment=0,logScaleFac=0) {
  if (missing(lambda)) {
    Cres <- .COMP_Z_integrand(z=z,eta=eta,nu=nu,moment=moment,logScaleFac=logScaleFac)
  } else Cres <- .COMP_Z_integrand(z,lambda,eta,nu,moment,logScaleFac)
  logint <- moment*log(z)+ z*eta - nu*lfactorial(z)
  res <- exp(logint-logScaleFac)
  res <- pmin(.Machine$double.xmax, res)
  res
}

## pure documentation for Rcpp version:
.pure_R_COMP_Z_asympto <- function(nu, pow_lam_nu) {
  logScaleFac <- nu*pow_lam_nu[[1L]] ## drop any name ! otherwise return value has wrong names
  # using Gaunt et al.:
  c1 <- (nu^2-1)/24
  c2 <- (nu^2-1)*(nu^2+23)/1152
  # c3 <- (nu^2-1)*(5*nu^4-298*nu^2+11237)/414720
  # c4 <- (nu^2-1)*(5*nu^6-1887*nu^4-241041*nu^2+2482411)/39813120
  # c5 <- (nu^2-1)*(7*nu^8-7420*nu^6+1451274*nu^4-220083004*nu^2+1363929895)/6688604160
  # c6 <- (nu^2-1)*(35*nu^10 - 78295*nu^8 + 76299326*nu^6 + 25171388146*nu^4 
  #                 - 915974552561*nu^2 + 4175309343349)/4815794995200
  # c7 <- (nu^2-1)*(5*nu^12- 20190*nu^10 + 45700491*nu^8 - 19956117988*nu^6
  #                 +7134232164555*nu^4-142838662997982*nu^2 + 525035501918789)/115579079884800
  # scaled <- (1+c1/logScaleFac+c2/logScaleFac^2+c3/logScaleFac^3+c4/logScaleFac^4+
  #              c5/logScaleFac^5+c6/logScaleFac^6+c7/logScaleFac^7)
  scaled <- (1+c1/logScaleFac+c2/logScaleFac^2)
  scaled <- (scaled/((2*pi*pow_lam_nu)^((nu-1)/2)*sqrt(nu)))[[1L]] ## drop any name !
  resu <- c(logScaleFac=logScaleFac, scaled=scaled)
  resu
} 

.get_quadinf <- local(
  {
    success <- quadinf <- NULL
    function() {
      if (is.null(success)) {
        tryres <- suppressWarnings(do.call("require",list(package="pracma", quietly = TRUE)))  ## package not declared in DESCRIPTION
        if (success <<- ! inherits(tryres,"try-error")) {
          quadinf <<- get("quadinf",envir = asNamespace("pracma"))
        } else {
          message(paste0("If the 'pracma' package were available,\n", "more accurate evaluation of an integral would be possible."))
          # quadinf remains NULL
        }
      }
      quadinf # NULL or the pracma::quadinf function
    }
  }
)

if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") {
  .COMP_local_check <- function(test) {
    if(test) browser("devel check for COMP ('__F I X M E___' tag: .COMP_Z_moment() includes check only conditionally on _LOCAL_TESTS_)")
    # from version 4.1.26, 2022/12/25 upwards
  }
} else {
  .COMP_local_check <- function(test) NULL # promise ignored...
}

## N O T the moments of the proba distr, Rather the 'num' for a moment obtained by .COMP_Z_ratio(num,denum)
# we call it repeatedly for nu=1 (and variable lambda) but this may be needed to get the continuity correction in nu=1
.COMP_Z_moment <- function(eta,nu,lambda=exp(eta), moment, 
                           pow_lam_nu = exp(eta/nu),
                           maxn=.COMP_maxn(lambda,nu),
                           use_asympto=.CMP_use_asympto(pow_lam_nu,nu,lambda)) {
  if (use_asympto) { # for moment=0
    if (moment) {
      stop("Execution should not reach this point.")  
    } else { 
      # resu <- .pure_R_COMP_Z_asympto(nu, pow_lam_nu)
      resu <- .Rcpp_COMP_Z_asympto(nu, pow_lam_nu)
    }
  } else {
    resu <- .Rcpp_COMP_Z(moment=moment,nu=nu,lambda=lambda,maxn=maxn)
    Z_mode <- ceiling(pow_lam_nu) ## practically locates the mode 
    lfac_Z_mode <- lfactorial(Z_mode)
    if (is.infinite(lfac_Z_mode)) { ## Should not occur as the asymptotic approx for the moments of the *PDF* should have been used 
      stop(paste("Practically infinite sum for COMPoisson's nu=",nu,". The asymptotic approx for the moments of the *PDF* should have been used."))
    } 
    if ( ! is.null(upper <- attr(maxn,"upper"))) { # If the NSum upper bound maxn was truncated by the safety control => 
      ## Add approximation for tail beyond summation from zero to maxn, using numerical integration.
      ## But integrating directly from maxn to Inf may not work, as integrate() may then miss the mode of the integrand. 
      ## So will be done in 2 stpes.
      logScaleFac <- (Z_mode*eta - nu*lfac_Z_mode)[[1L]] ## drop any name !
      ## (1) Optionally add integral from maxn to Z_mode, if Z_mode is higher 
      if (maxn<Z_mode) {  
        scaled <- max(0,integrate(.COMP_Z_integrand,lower=maxn,upper=Z_mode,eta=eta,nu=nu,moment=moment,
                                           logScaleFac=logScaleFac)$value)
        lower <- Z_mode # maxn < Z_mode < upper   (upper>maxn must be ensured by .COMP_maxn())
      } else {
        scaled <- 0
        lower <- maxn # Z_mode < maxn < upper
      }
      # (2) Always add tail beyond lower := the highest of maxn and Z_mode
      if (lower>1e7) { # has been >1e10
        
        # checking that the integrand is not more that numerical imprecision at upper bound:
        .COMP_local_check(test= {
          crit <- ((upper)^moment) *
            .COMP_Z_integrand(upper, eta = eta, lambda = NULL, nu = nu, moment = moment, logScaleFac = logScaleFac)
          crit > 1e-12
        }) 
        # E.g., compare the value of integral in (spaMM:::.COMP_Z_n(lambda=exp(7.321), nu=0.25)) verss (spaMM:::.COMP_Z_n(lambda=exp(7.32), nu=0.25))
        
        if (is.null(quadinf <- .get_quadinf())) {
          # nintegrate behaves more poorly that what can be handled by max(0,.)
          ## With upper=Inf, integrate once failed, with a message suggesting a diverging integrand, which it was not. 
          scaled <- max(0,integrate(.COMP_Z_integrand,lower=lower,upper=upper, # has been # , upper=Inf # => cf above comment
                                    eta=eta,nu=nu,moment=moment,
                                    logScaleFac=logScaleFac)$value)
        } else {
          scaled <- quadinf(f=.COMP_Z_integrand, xa = lower, xb = upper,    # has been # , xb=Inf #
                            eta = eta, nu = nu, moment = moment, logScaleFac = logScaleFac)$Q #Bug fixed in v4.1.25 here
        }
        # quadinf_res <- .do_call_wrap("quadinf", 
        #                         pack="pracma",
        #                         info_mess=paste0("If the 'pracma' package were available,\n",
        #                                          "more accurate evaluation of an integral would be possible.")
        # )$Q
      }
      scaled <- pmin(.Machine$double.xmax,scaled)
      resu <- .COMP_Z_sum(resu, c(logScaleFac=logScaleFac,scaled=scaled))
    }
  }
  return(resu)
}

##  the moments of the proba distr

.COMP_Pr_moments <- function(pow_lam_nu=lambda^(1/nu), lambda, nu, moments,
                             use_asympto=.CMP_use_asympto(pow_lam_nu,nu, lambda)) {
  whichM <- which(moments)
  resu <- rep(NA,length(whichM))
  names(resu) <- whichM
  if (use_asympto) {
    # using Gaunt et al.: ...3.22...
    c1 <- (nu^2-1)/24
    mu1 <- pow_lam_nu*( 1 -(nu-1)/(2*nu*pow_lam_nu) - c1/((nu*pow_lam_nu)^2)- c1/((nu*pow_lam_nu)^3))
    if (moments[1L]) resu["1"] <- mu1
    if (moments[2L] || moments[3L] || moments[4L]) {
      Var <- (pow_lam_nu/nu)*( 1 + c1/((nu*pow_lam_nu)^2) + 2*c1/((nu*pow_lam_nu)^3)) # >0 when nu>1
      if (moments[2L]) resu["2"] <- mu1^2 + Var # supp to mu1^2 when nu>1
    }
    if (moments[3L] || moments[4L]) {
      k3 <- (pow_lam_nu/nu^2)*( 1 - c1/((nu*pow_lam_nu)^2) -4*c1/((nu*pow_lam_nu)^3))
      resu["3"] <- mu1^3 + 3*mu1*Var +k3
    }
    if (moments[4L]) {
      k4 <- (pow_lam_nu/nu^3)*( 1 + c1/((nu*pow_lam_nu)^2) +8*c1/((nu*pow_lam_nu)^3))
      resu["4"] <- mu1^4 + 6*Var*mu1^2 + 3*Var^2 + 4*k3*mu1 + k4
    }
  } else {
    denum <- .COMP_Z(lambda=lambda,nu=nu) # -> calling .COMP_Pr_moments( use_asympto=TRUE) ... 
    denum_corr <- .COMP_Z(lambda=lambda,nu=1) # for continuity correction in nu=1
    if (moments[1L]) {
      num <- .COMP_Z_n(lambda=lambda,nu=nu)
      resu["1"] <- .COMP_Z_ratio(num,denum)
      # continuity correction wrt poisson: corr =lambda + error of the COMP_Z... functions
      corr <- .COMP_Z_ratio(.COMP_Z_n(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
      resu["1"] <- resu["1"]+(lambda-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
    }
    if (moments[2L]) {
      num <- .COMP_Z_n2(lambda=lambda,nu=nu)
      resu["2"] <- .COMP_Z_ratio(num,denum)
      # cotinuity correction wrt poisson: 
      corr <- .COMP_Z_ratio(.COMP_Z_n2(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
      resu["2"] <- resu["2"]+(lambda*(1+lambda)-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
    }
    if (moments[3L]) {
      num <- .COMP_Z_n3(lambda=lambda,nu=nu)
      resu["3"] <- .COMP_Z_ratio(num,denum)
      # cotinuity correction wrt poisson: 
      corr <- .COMP_Z_ratio(.COMP_Z_n3(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
      resu["3"] <- resu["3"]+(lambda*(1+lambda*(3+lambda))-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
    }
    if (moments[4L]) {
      num <- .COMP_Z_n4(lambda=lambda,nu=nu)
      resu["4"] <- .COMP_Z_ratio(num,denum)
      # cotinuity correction wrt poisson: 
      corr <- .COMP_Z_ratio(.COMP_Z_n4(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
      resu["4"] <- resu["4"]+(lambda*(1+lambda*(7+lambda*(6+lambda)))-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
    }
  }
  return(resu)
}


.COMP_Z <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  if (lambda==0) return(c(logScaleFac=0,scaled=1))  
  ## ELSE
  if (missing(eta)) eta <- log(lambda)
  resu <- .COMP_Z_moment(moment=0,eta=eta,nu=nu,lambda=lambda,maxn=maxn)
  return(resu)
  # # R version of the summation term:
  # if (nu==0) {
  #   return(c(logScaleFac=0,scaled=as.numeric(1/(1-lambda))))
  # } 
  # floorn <- floor(maxn)
  # epsn <- maxn - floorn
  # facs <- c(1,lambda/seq(floorn+1L)^nu)
  # cumprodfacs <- cumprod(facs)
  # cumprodfacs[floorn+2L] <- cumprodfacs[floorn+2L]*epsn
  # scaled <- sum(cumprodfacs)
  # if ( ! is.finite(scaled)) {
  #   logfacs <- log(facs)
  #   cumsumlogfacs <- cumsum(logfacs)
  #   refi <- which.max(cumsumlogfacs)
  #   logScaleFac <- cumsumlogfacs[refi]
  #   scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  # } else logScaleFac <- 0
  # return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z = exp(logScaleFac)*scaled
}

.COMP_Z_n <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  if (missing(eta)) eta <- log(lambda)
  resu <- .COMP_Z_moment(moment=1,eta=eta,nu=nu,lambda=lambda,maxn=maxn, use_asympto=FALSE)
  return(resu)
  # # R version:
  # if (nu==0) {
  #   return(c(logScaleFac=0,scaled=as.numeric(lambda/(1-lambda)^2)))
  # } 
  # floorn <- floor(maxn)
  # epsn <- maxn - floorn
  # seqn <- seq(2L,floorn+1L)
  # facs <- c(lambda,(seqn/(seqn-1L))*lambda/(seqn^nu))
  # cumprodfacs <- cumprod(facs)
  # cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  # scaled <- sum(cumprodfacs)
  # if ( ! is.finite(scaled)) {
  #   logfacs <- log(facs)
  #   cumsumlogfacs <- cumsum(logfacs)
  #   refi <- which.max(cumsumlogfacs)
  #   logScaleFac <- cumsumlogfacs[refi]
  #   scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  # } else logScaleFac <- 0
  # return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n = exp(logScaleFac)*scaled
}

.COMP_Z_n2 <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  if (missing(eta)) eta <- log(lambda)
  resu <- .COMP_Z_moment(moment=2,eta=eta,nu=nu,lambda=lambda,maxn=maxn, use_asympto=FALSE)
  return(resu)
  # # R version:
  # if (nu==0) {
  #   return(c(logScaleFac=0,scaled=as.numeric(lambda*(1+lambda)/(1-lambda)^3)))
  # } 
  # floorn <- floor(maxn)
  # epsn <- maxn - floorn
  # seqn <- seq(2L,floorn+1L)
  # facs <- c(lambda,((seqn/(seqn-1L))^2)*lambda/(seqn^nu))
  # cumprodfacs <- cumprod(facs)
  # cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  # scaled <- sum(cumprodfacs)
  # if ( ! is.finite(scaled)) {
  #   logfacs <- log(facs)
  #   cumsumlogfacs <- cumsum(logfacs)
  #   refi <- which.max(cumsumlogfacs)
  #   logScaleFac <- cumsumlogfacs[refi]
  #   scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  # } else logScaleFac <- 0
  # return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n2 = exp(logScaleFac)*scaled
}

.COMP_Z_n3 <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  if (missing(eta)) eta <- log(lambda)
  resu <- .COMP_Z_moment(moment=3,eta=eta,nu=nu,lambda=lambda,maxn=maxn, use_asympto=FALSE)
  return(resu)
  # # R version:
  # if (nu==0) {
  #   return(c(logScaleFac=0,scaled=as.numeric(lambda*(1+4*lambda+lambda^2)/(1-lambda)^4)))
  # } 
  # floorn <- floor(maxn)
  # epsn <- maxn - floorn
  # seqn <- seq(2L,floorn+1L)
  # facs <- c(lambda,((seqn/(seqn-1L))^3)*lambda/(seqn^nu))
  # cumprodfacs <- cumprod(facs)
  # cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  # scaled <- sum(cumprodfacs)
  # if ( ! is.finite(scaled)) {
  #   logfacs <- log(facs)
  #   cumsumlogfacs <- cumsum(logfacs)
  #   refi <- which.max(cumsumlogfacs)
  #   logScaleFac <- cumsumlogfacs[refi]
  #   scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  # } else logScaleFac <- 0
  # return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n3 = exp(logScaleFac)*scaled
}

.COMP_Z_n4 <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  if (missing(eta)) eta <- log(lambda)
  resu <- .COMP_Z_moment(moment=4,eta=eta,nu=nu,lambda=lambda,maxn=maxn, use_asympto=FALSE)
  return(resu)
}



.COMP_Z_lfacn <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,lambda/(seqn^nu)) ## n from 1 to floorn+1
  cumprodfacs <- cumprod(facs) ## sequence of lambda^n / n!^nu 
  cumsumlogn <- cumsum(c(0,log(seqn))) ## sequence of log(n!) for n from 1 to floorn+1
  cumprodfacs <- cumprodfacs*cumsumlogn ## sequence of log(n!) lambda^n / n!^nu 
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  scaled <- sum(cumprodfacs)
  if ( ! is.finite(scaled)) {
    logfacs <- log(facs)
    cumsumlogfacs <- cumsum(logfacs)
    refi <- which.max(cumsumlogfacs)
    logScaleFac <- cumsumlogfacs[refi]
    scaled <- sum(exp(cumsumlogfacs-logScaleFac))
  } else logScaleFac <- 0
  return(c(logScaleFac=logScaleFac,scaled=scaled)) ## Z_n = exp(logScaleFac)*scaled
}


.COMP_Z_ratio <- function(Z1,Z2,log=FALSE) {
  logratio <- Z1[["logScaleFac"]]-Z2[["logScaleFac"]]+log(Z1[["scaled"]]/Z2[["scaled"]]) # more overflow-robust than directly using exp()
  if (log) {
    return(logratio)
  } else return(exp(logratio))
} 

.COMP_Z_sum <- function(Z1,Z2) {
  z12sc <- Z1[["logScaleFac"]] - Z2[["logScaleFac"]]
  #keep largest logScaleFac
  if (z12sc>0) {
    logScaleFac <- Z1[["logScaleFac"]]
    Z1sc <- Z1[["scaled"]]
    Z2sc <- Z2[["scaled"]] / exp(z12sc)
  } else {
    logScaleFac <- Z2[["logScaleFac"]]
    Z1sc <- Z1[["scaled"]] * exp(z12sc)
    Z2sc <- Z2[["scaled"]]
  }
  scaled <- pmin(.Machine$double.xmax, Z1sc+Z2sc)
  return(c(logScaleFac=logScaleFac,scaled=scaled))
} 

.COMP_dnu_objfn <- function(nu,y,eta,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)) { ## not currently used
  summand <- numeric(length(y))
  for (i in seq_len(length(y))) {
    yi <- y[i]
    lambdai <- lambda[i]
    comp_z_lfacn <- .COMP_Z_lfacn(nu=nu,lambda=lambdai,maxn=maxn)
    comp_z <- .COMP_Z(nu=nu,lambda=lambdai,maxn=maxn)
    summand[i] <- lfactorial(yi)-.COMP_Z_ratio(comp_z_lfacn,comp_z,log=TRUE)
  }
  return(sum(summand))
}

# for .r_resid_var():
.dCOMP <- function(x, mu, family_env,
                   nu=family_env$nu,
                   lambda= family_env$mu2lambda(mu), # COMPoisson(nu=nu)$linkfun(mu,log=FALSE),
                   log = FALSE, maxn=.COMP_maxn(lambda,nu)) {
  compz <- .COMP_Z(lambda=lambda,nu=nu,maxn=maxn)
  logd <- x * log(lambda) - nu* lfactorial(x) - compz[["logScaleFac"]] -log(compz[["scaled"]])
  if (log) { 
    return(logd)
  } else return(exp(logd))
}

.COMP_simulate <- function(lambda,nu,nsim=1L) { ## for scalar lambda
  maxn <- .COMP_maxn(lambda=lambda,nu=nu)
  Z <- .COMP_Z(lambda=lambda,nu=nu)
  floorn <- floor(maxn)
  epsn <- maxn - floorn
  seqn <- seq(2L,floorn+1L)
  facs <- c(lambda,(seqn/(seqn-1L))*lambda/(seqn^nu))
  cumprodfacs <- cumprod(facs)
  cumprodfacs[floorn+1L] <- cumprodfacs[floorn+1L]*epsn
  sample(floorn+1L,size=nsim,prob=cumprodfacs/sum(cumprodfacs), replace=TRUE)
}

.CMP_lambda2mui <- function(lambda,nu) { # mu(lambda)
  if (lambda==0) {
    return(1e-8)  
  } else {
    mu <- .COMP_Pr_moments(lambda=lambda, nu=nu, moments=c(TRUE, FALSE, FALSE, FALSE))
    return(mu) #(max(mu,1e-8)) would avoid mu=0 which fails validmu(mu) (as for poisson family)
    # but it's better to use a consistent code to sanitize the input lambda than to fix the return value.
  }
}

.CMP_lambda2mu <- function(lambda) { # mu(lambda)  
  mus <- attr(lambda,"mu") 
  if (is.null(mus)) {
    if (nu==0) {
      mus <- lambda/(1-lambda) 
    } else {
      nu <- parent.env(environment())$nu
      # mus <- sapply(lambda, .CMP_lambda2mui,nu=nu) 
      mus <- numeric(length(lambda))
      for (it in seq_along(lambda)) mus[it] <- .CMP_lambda2mui(lambda[it],nu=nu)
    }
    dim(mus) <- dim(lambda) ## may be NULL
  }
  attributes(lambda) <- NULL
  return(structure(mus,lambda=lambda))
}

# different input
.CMP_loglambda_linkinv <- function(eta,lambda=exp(eta), muetaenv=NULL) { # \equiv .CMP_lambda2mu(lambda) but with nu properly passed...
  mus <- attr(lambda,"mu") 
  if (is.null(mus)) {
    if (is.null(muetaenv)) {
      if (nu==0) {
        mus <- lambda/(1-lambda) 
      } else {
        nu <- parent.env(environment())$nu # for canonical link .CMP_loglambda_linkinv has been copied to $linkinv with nu in its parent env.
        #mus <- sapply(lambda, .CMP_lambda2mui,nu=nu)
        mus <- numeric(length(lambda))
        for (it in seq_along(lambda)) mus[it] <- .CMP_lambda2mui(lambda[it],nu=nu)
      }
      dim(mus) <- dim(lambda) ## may be NULL
    } else {
      mus <- muetaenv$EX
    } 
  }
  attributes(lambda) <- NULL
  return(structure(mus,lambda=lambda))
}


# called by the next two functions
..CMP_mu2lambda <- function(mu,nu, CMP_linkfun_objfn) { 
  if (nu==1) {
    lambda <- mu ## pb du code general est qu'objfn n'est alors que l'erreur numérique de linkinv()
  } else if (mu==Inf) {
    warning(paste("Approximating lambda as 1-1e-8 in COMPoisson(nu=",nu,") for mu = Inf"))
    lambda <- 1 - 1e-8 ## 1 -> mu.eta = NaN. mu.eta being Inf would not be OK bc Inf GLMweights -> 1/w.resid=0 -> singular Sig or d2hdv2
  } else if (mu==0) { # occurs when computing attr(y,"CMP")
    lambda <- 0
  } else {
    app_lambda <-  max(c(0,(mu+max(0,(nu-1)/(2*nu)))))^(nu)
    lambdamin <- max(c(.Machine$double.eps, app_lambda-1)) ## (the ultimate constraint is that log(lambda)> -Inf)
    lambdamax <- max(c(.Machine$double.eps, (app_lambda+1)*c(1.01))) ## not perfect...
    # last one for low lambda,nu values
    fupper <- CMP_linkfun_objfn(lambdamax, mu=mu)
    while (fupper<0) {
      lambdamax <- lambdamax*2
      fupper <- CMP_linkfun_objfn(lambdamax, mu=mu)
      if (is.nan(fupper)) break
    }
    ## normally >0 
    if (is.nan(fupper)) {
      if ( ! identical(spaMM.getOption("COMP_geom_approx_warned"),TRUE)) {
        warning(paste("Geometric approximation tried in COMPoisson(nu=",nu,") for mu=",mu," and possibly for other nu,mu values"))
        .spaMM.data$options$COMP_geom_approx_warned <- TRUE
      }
      lambda <- mu/(1+mu) ## but it should be checked hence the warning. Ideally for large nu ?
    } else {
      flower <- CMP_linkfun_objfn(lambdamin, mu=mu) ## normally <0 
      while ((flower)>0) {
        lambdamin <- lambdamin/2
        flower <- CMP_linkfun_objfn(lambdamin, mu=mu)
        if (is.nan(flower)) break
      }
      interval <- c(lambdamin,lambdamax)
      # linkinv(0)=0 => objfn(0)= -mu
      lambda <- uniroot(CMP_linkfun_objfn,interval=interval,f.lower = flower,f.upper = fupper, mu=mu)$root
    }
  }
  return(lambda)
}

# different arguments
.CMP_mu2lambda <- function(mu) {## scalar or vector
  lambdas <- attr(mu,"lambda")
  if ( is.null(lambdas)) {
    nu <- parent.env(environment())$nu
    if (nu==0) {
      lambdas <- mu/(1+mu) 
    } else {
      uniqmu <- unique(mu)
      uniqlambda <- numeric(length(uniqmu))
      for (it in seq_along(uniqmu)) uniqlambda[it] <- ..CMP_mu2lambda(uniqmu[it],nu=nu, CMP_linkfun_objfn=parent.env(environment())$CMP_linkfun_objfn)
      lambdas <- uniqlambda[match(mu, uniqmu)]
    }
  } 
  attributes(mu) <- NULL ## avoids 'mise en abime'
  return(structure(lambdas,mu=mu))
}


# different return value
# Defines COMPoisson(link="loglambda")$linkfun 
.CMP_loglambda_linkfun <- function(mu) {## scalar or vector
  lambdas <- attr(mu,"lambda")
  if ( is.null(lambdas)) {
    nu <- parent.env(environment())$nu
    if (nu==0) {
      lambdas <- mu/(1+mu) 
    } else {
      uniqmu <- unique(mu)
      uniqlambda <- numeric(length(uniqmu))
      for (it in seq_along(uniqmu)) uniqlambda[it] <- ..CMP_mu2lambda(uniqmu[it],nu=nu, CMP_linkfun_objfn=parent.env(environment())$CMP_linkfun_objfn)
      lambdas <- uniqlambda[match(mu, uniqmu)]
    }
  } 
  attributes(mu) <- NULL ## avoids 'mise en abime'
  return(structure(log(lambdas),mu=mu)) ## eta (ie standard linkfun value) for the Shmueli version
}

# quickly becoming obsolete
.CMP_dev.resid <- function(yi,lambdai,nu, CMP_linkfun_objfn) {
  Z2 <- .COMP_Z(lambda=lambdai,nu=nu)
  if (yi==0) { # lambda = 0,  Z1 = 1
    dev <- 2*(Z2[["logScaleFac"]]+log(Z2[["scaled"]]))
  } else {
    # eval lambda for Z1(lambda...)
    if (nu==0) {
      lambda <- yi/(1+yi)  
    } else {
      app_lambda <-  max(c(0,(yi+max(0,(nu-1)/(2*nu)))))^(nu)
      lambdamin <- max(c(.Machine$double.eps, app_lambda-1))
      while (CMP_linkfun_objfn(lambdamin, mu=yi)>0) {
        lambdamin <- lambdamin/2 ## greedy 
        if (lambdamin<1e-8) stop("lambdamin<1e-8 in COMPoisson$dev.resids()")
      }
      lambdamax <- max(c(.Machine$double.eps, (app_lambda+1)*c(1.01))) 
      while (CMP_linkfun_objfn(lambdamax, mu=yi)<0) {
        lambdamax <- lambdamax*10 ## greedy but resolves errors for low nu 
        if (lambdamax>1e100) stop("lambdamax>1e100 in COMPoisson$dev.resids()")
      }
      interval <- c(lambdamin,lambdamax)
      lambda <- uniroot(CMP_linkfun_objfn,interval=interval, mu=yi)$root
    }
    Z1 <- .COMP_Z(lambda=lambda,nu=nu)
    #
    dev <- 2*(yi*log(lambda/lambdai)-(Z1[["logScaleFac"]]-Z2[["logScaleFac"]]+log(Z1[["scaled"]]/Z2[["scaled"]])))
  }
  dev
}

.CMP_attr_y <- function(y, family) { # the 'family'  can be the family working envir
  Z1_y <- vector("list", length(y))
  lambdas_y <- family$mu2lambda(mu=y) # function of y and nu, with mu=y ("saturated model") 
  nu <- environment(family$aic)$nu
  for(i in seq_along(y)) if (y[i]>0) Z1_y[[i]] <- .COMP_Z(lambda=lambdas_y[[i]],nu=nu)
  list(lambdas_y=lambdas_y,Z1_y=Z1_y)
}

.CMP_dev_resids <- function(y, mu, 
                            wt, # currently ignored, see comment below 
                            muetaenv=NULL){
  # must accept, among others, vector y and scalar mu.
  familyenv <- parent.env(environment()) # parent envir of COMPoisson()$dev.resids
  if (is.null(muetaenv)) {
    lambdas <- familyenv$mu2lambda(mu=mu) 
  } else lambdas <- muetaenv$lambdas
  nu <- familyenv$nu
  n <- length(y)
  #
  if (is.null(y_CMP <- attr(y,"CMP"))) { # In spaMM_glm.fit, the attribute is precomputed before the main loop.
    # But calculation is mostly avoided in internal spaMM procedures HLfit_body() and glm.nodev.fit().
    # In HLfit_body, $dev.resids seems to be called only by .calcPHIs() for families without a fixed phi.
    # glm.nodev.fit avoids it.
    # So a spaMM fit still requests it only optionally post xLM fit in .get_inits_by_xLM() is called, 
    # notably for the initial lambda value in outer-optimization.
    # (however, there may still be an inefficiency for mv fits)
    y_CMP <- .CMP_attr_y(y, family=familyenv)  # adds mu2lambda(mu=y) and the resulting COMP_Z
  }
  lambdas_y <- y_CMP$lambdas_y
  Z1_y <- y_CMP$Z1_y
  #
  if (length(mu)==1L) lambdas <- rep(lambdas,n)
  devs <- numeric(n)
  #CMP_linkfun_objfn <- parent.env(environment())$CMP_linkfun_objfn
  #for(i in seq(n)) { devs[i] <- .CMP_dev.resid(y[i],lambdai=lambdas[i],nu=nu, CMP_linkfun_objfn=CMP_linkfun_objfn) }
  for(i in seq(n)) {
    if (is.null(muetaenv)) {
      Z2 <- .COMP_Z(lambda=lambdas[i],nu=nu)
    } else Z2 <- muetaenv$denum_Z[[i]]
    if (y[i]==0) { #lambda = 0,  Z1 = 1
      devs[i] <- 2*(Z2[["logScaleFac"]]+log(Z2[["scaled"]]))
    } else {
      Z1_yi <- Z1_y[[i]]
      devs[i] <- 2*(y[i]*log(lambdas_y[[i]]/lambdas[i])-(Z1_yi[["logScaleFac"]]-Z2[["logScaleFac"]]+log(Z1_yi[["scaled"]]/Z2[["scaled"]]))) 
    }
  }
  # devs <- devs*wt # seems as inappropriate as for other families where the dispersion parameter is not a scaling factor for the logl 
  devs[devs==Inf] <- .Machine$double.xmax/n ## so that total deviance may be finite 
  return(devs) 
}

# Related to the 'denominator' of .CMP_dev_resids(), but as function of the canonical parameter theta.
.CMP_clik_fn <- function(theta,y, 
                         nu, # presumably a vector of 1
                         COMP_nu, 
                         muetaenv) { ## theta = log(lambda) 
  logLs <- numeric(length(y))
  for (i in seq(length(y))) {
    if (is.null(muetaenv)) {
      comp_z <- .COMP_Z(lambda=exp(theta[i]),nu=COMP_nu) 
    } else comp_z <- muetaenv$denum_Z[[i]]
    logLs[i] <- nu[i]*(theta[i]*y[i]-comp_z[[1]]-log(comp_z[[2]])) - COMP_nu * lfactorial(y[i])
  }
  logLs[theta== -Inf & y==0] <- 1
  logLs
}

.CMP_aic <- function(y, n, mu, wt, dev) {
  family_env <- parent.env(environment())
  aici <- numeric(length(y))
  for (i in seq_len(length(y))) { aici[i] <- .dCOMP(y[i], mu[i], family_env=family_env, log = TRUE) }
  -2 * sum(aici * wt)
}

.CMP_simfun <- function(object,nsim) {
  wts <- object$prior.weights
  if (any(wts != 1)) 
    warning("ignoring prior weights")
  if (inherits(object,"glm")) { # if COMPoisson() used in a glm() call simulate.lm -> family$simulate -> here
    lambdas <- attr(object$fitted.values,"lambda")
    if (is.null(lambdas)) lambdas <- sapply(mu, family$mu2lambda)
  } else { # HLfit object: 
    mu <- object$muetablob$mu
    lambdas <- attr(mu,"lambda") # exp(object$eta)
    if (is.null(lambdas)) lambdas <- sapply(mu, family$mu2lambda)
  }
  nu <- parent.env(environment())$nu
  resu <- sapply(lambdas,.COMP_simulate,nu=nu,nsim=nsim) # vector if nsim=1, nsim-row matrix otherwise
  if (nsim>1) resu <- t(resu)
  return(resu)  # vector if nsim=1, nsim-col matrix otherwise
}

.CMP_mu.eta <- function (eta,lambda=exp(eta), muetaenv=NULL) {
  nu <- parent.env(environment())$nu ## not necess for R to find it, but to avoid a NOTE from CHECK
  if (nu==0) {
    resu <- lambda/(1-lambda)^2
  } else {
    if (is.null(muetaenv)) {
      resu <- numeric(length(lambda))
      for (it in seq_len(length(lambda))) {
        lambdai <- lambda[[it]]
        if (lambdai==0) {
          resu[[it]] <- .Machine$double.eps
        } else {
          moments <- .COMP_Pr_moments(lambda=lambdai, nu=nu, moments=c(TRUE, TRUE, FALSE, FALSE))
          res <- moments["2"] - moments["1"]^2
          resu[[it]] <- pmax(res, .Machine$double.eps)
        }
      }
    } else {
      resu <- muetaenv$EX2 - muetaenv$EX^2
      resu <- pmax(resu, .Machine$double.eps)
    }
  }
  resu
}

.CMP_variance <- function(mu, muetaenv=NULL) { # same as .CMP_mu.eta but with different argument
  nu <- parent.env(environment())$nu
  if (nu==0) {
    return(mu*(1+mu)) 
  } else {
    if (is.null(muetaenv)) {
      lambdas <- parent.env(environment())$mu2lambda(mu)
      len <- length(lambdas)
      resu <- numeric(len)
      for (it in seq_len(len)) {
        lambdai <- lambdas[[it]]
        if (lambdai==0) {
          resu[[it]] <- 1e-8
        } else {
          moments <- .COMP_Pr_moments(lambda=lambdai,nu=nu,moments=c(TRUE, TRUE, FALSE, FALSE))
          resu[[it]] <- max(moments["2"] -moments["1"]^2,1e-8) ## pmax otherwise for low mu, Vmu=0, -> ... -> w.resid=0
        }
      } 
    } else {
      resu <- pmax(muetaenv$EX2 -muetaenv$EX^2,1e-8) 
    }
    ## for geom V(mu) = mu(1+mu) but this cannot be a lower bound for resu.
    return(resu)
  }
}

.CMP_thetaMuDerivs_2 <- function(muetaenv) { # computation of derivatives of inverse of .CMP_calc_dlW_deta_locfn...
  
  # computation of derivatives of inverse of .CMP_calc_dlW_deta_locfn...
  if (muetaenv$nu==1) {
    mu <- muetaenv$mu
    return(list(Dtheta.Dmu=1/mu, D2theta.Dmu2= - 1/mu^2))
  }
  # ELSE
  len <- muetaenv$len
  lambdas <- muetaenv$lambdas
  EX <- muetaenv$EX
  EX2 <- muetaenv$EX2
  EX3 <- muetaenv$EX3
  Dtheta.Dmu <- D2theta.Dmu2 <- numeric(len)
  for (i in seq_len(len)) {
    lambdai <- lambdas[[i]]
    if (lambdai==0) {
      Dtheta.Dmu[i] <- 1e8 # 1/lambda
      D2theta.Dmu2[i] <- - 1e16 # -1/lambda^2
    } else {
      EXi <- EX[[i]]
      V <- EX2[[i]] - EXi^2 # ...that's the family $variance()...
      Dtheta.Dmu[i] <- 1/V
      # D2theta.Dmu2[i] <- (EX2[[i]] * (1 + EXi) - EX3[[i]] - V*(1-2*EXi) - EXi^2)/V^3 # Same as:
      D2theta.Dmu2[i] <- (3*V*EXi - EX3[[i]] + EXi^3)  /V^3 # computation to COMP_Z.nb
    }
  }
  
  return(list(Dtheta.Dmu=Dtheta.Dmu,D2theta.Dmu2=D2theta.Dmu2))
}

.CMP_thetaMuDerivs_3 <- function(muetaenv) { # computation of derivatives of inverse of .CMP_calc_dlW_deta_locfn... 
  # code redundancy with previous fn seems acceptable (OK from efficiency viewpoint: functions not called on the same muetaenv)
  
  # computation of derivatives of inverse of .CMP_calc_dlW_deta_locfn...
  if (muetaenv$nu==1) {
    mu <- muetaenv$mu
    return(list(Dtheta.Dmu=1/mu, D2theta.Dmu2= - 1/mu^2, D3theta.Dmu3= 2/mu^3)) # for theta= log(mu)
  }
  # ELSE
  len <- muetaenv$len
  lambdas <- muetaenv$lambdas
  EX <- muetaenv$EX
  EX2 <- muetaenv$EX2
  EX3 <- muetaenv$EX3
  EX4 <- muetaenv$EX4
  Dtheta.Dmu <- D2theta.Dmu2 <- D3theta.Dmu3 <- numeric(len)
  for (i in seq_len(len)) {
    lambdai <- lambdas[[i]]
    if (lambdai==0) {
      Dtheta.Dmu[i] <- 1e8 # 1/lambda
      D2theta.Dmu2[i] <- - 1e16 # -1/lambda^2
      D3theta.Dmu3[i] <- - 2e24 # 2/lambda^3
    } else {
      EXi <- EX[[i]]
      V <- EX2[[i]] - EXi^2 # ...that's the family $variance()...
      Dtheta.Dmu[i] <- 1/V
      # D2theta.Dmu2[i] <- (EX2[[i]] * (1 + EXi) - EX3[[i]] - V*(1-2*EXi) - EXi^2)/V^3 # Same as:
      D2theta.Dmu2[i] <- (3*V*EXi - EX3[[i]] + EXi^3)  /V^3 # computation in COMP_Z.nb
      D3theta.Dmu3[i] <- 3*V*D2theta.Dmu2[i]^2 +
        (3*V*EX2[[i]]- EX4[[i]]+EX3[[i]]*EXi-3*EXi*D2theta.Dmu2[i]*V^3)/V^4
    }
  }
  
  return(list(Dtheta.Dmu=Dtheta.Dmu, D2theta.Dmu2=D2theta.Dmu2, D3theta.Dmu3=D3theta.Dmu3))
}



COMPoisson <- function(nu = stop("COMPoisson's 'nu' must be specified"), 
                       link = "loglambda" # eta <-> mu link, not the eta <-> lambda log link
) { 
  .spaMM.data$options$COMP_maxn_warned <- FALSE # much better here than in .preprocess(); works with glm()
  .spaMM.data$options$COMP_geom_approx_warned <- FALSE
  resid.model <- list2env(list(off=0)) # env so that when we assign to it we don't create a new instance of the family object
  if (inherits(nuch <- substitute(nu),"character") ||
      (inherits(nuch,"name") && inherits(nu, "function")) # "name" is for e.g. COMPoisson(log)
      # (but testing only "name" would catch e.g. COMPoisson(nu=nu) )
     ) { 
    if (inherits(nuch,"character")) nuch <- paste0('"',nuch,'"')
    errmess <- paste0('It looks like COMPoisson(',nuch,') was called, which absurdly means COMPoisson(nu=',nuch,
                      ').\n  Use named argument: COMPoisson(link=',nuch,') instead.')
    stop(errmess)
  }
  # When 'nu' is recognized as as call to some function ! = stop(), we eval it so it is no longer recognized as a call by .calc_optim_args()
  if (inherits(nuch,"call") && deparse(nuch[[1]])!="stop") nu <- eval(nuch) 
  
  linktemp <- substitute(link) # if link was char LHS is char ; else deparse will create a char from a language object 
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("loglambda", "log", "identity", "sqrt")
  if (inherits(link, "link-glm")) { # a make.link() object was provided
    stats <- link
    if ( ! is.null(stats$name)) linktemp <- stats$name
  } else if (linktemp=="loglambda") {
    linkfun <- .CMP_loglambda_linkfun # need to associate nu to this one (not to more classical linkfun's)
    linkinv <- .CMP_loglambda_linkinv
    mu.eta <- .CMP_mu.eta
    environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment()  ## containing nu
  } else if (linktemp %in% okLinks) {
    stats <- make.link(linktemp)
    linkfun <- stats$linkfun
    linkinv <- stats$linkinv
    mu.eta <- stats$mu.eta
  } # at this point 'stats' is always the result of make.link and 'linktemp' is always char. 
  if ( ! linktemp %in% okLinks) stop(gettextf("link \"%s\" not available for COMPoisson family; available links are %s", 
                                              linktemp, paste(sQuote(okLinks), collapse = ", ")),           domain = NA)
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0) ## from poisson()
  valideta <- function(eta) TRUE ## from poisson()
  mu2lambda <- .CMP_mu2lambda # copy to pass nu, used in in .CMP_dev_resids2()
  lambda2mu <- .CMP_lambda2mu # ## copied locally to pass nu to .CMP_lambda2mui() through the lambda2mu parent envir
  CMP_linkfun_objfn <- function(lambda, mu) {lambda2mu(lambda=lambda) -mu} # objective fn of ~.CMP_mu2lambda # called by .CMP_dev.resid() to find lambda=lambda(mu)
  # CMP_linkfun_objfn_log <- function(loglambda, logmu) {log(lambda2mu(lambda=exp(loglambda))) -logmu} # objective fn of ~.CMP_mu2lambda # called by .CMP_dev.resid() to find lambda=lambda(mu)
  variance <- .CMP_variance
  dev.resids <- .CMP_dev_resids
  aic <- .CMP_aic
  initialize <- expression({
    if (any(y < 0L)) stop("negative values not allowed for the 'COMPoisson' family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  DlogLDmu <- function(mu, y, thetaMuDerivs) { drop((y-mu)*thetaMuDerivs$Dtheta.Dmu) } # inefficiency as optional muetaenv not available # also (y-mu)*.CMP_thetaMuDerivs(mu, family,muetaenv, )$Dtheta.Dmu
  D2logLDmu2 <- function(mu, y, thetaMuDerivs) { drop( -thetaMuDerivs$Dtheta.Dmu +(y-mu)*thetaMuDerivs$D2theta.Dmu2) }
  D3logLDmu3 <- function(mu, y, thetaMuDerivs) { drop( -2*thetaMuDerivs$D2theta.Dmu2 + (y-mu)*thetaMuDerivs$D3theta.Dmu3) }
  D2muDeta2 <- .D2muDeta2(linktemp)
  D3muDeta3 <- .D3muDeta3(linktemp)
  simulate <- .CMP_simfun 

  environment(dev.resids) <- environment(aic) <- environment(simulate) <- environment(mu2lambda) <- 
    environment(variance) <- environment(lambda2mu) <- environment() ## containing nu
  
  ## Change the parent.env of all functions that have the same envir as aic(): 
  parent.env(environment(aic)) <- environment(.dCOMP) ## gives access to spaMM:::.dCOMP and other .COMP_ fns
  structure(
    list(
      family = structure("COMPoisson",
                         withArgs=quote(paste0("COMPoisson(nu=",signif(nu,4),")"))), 
      link = linktemp, linkfun = linkfun,           mu2lambda=mu2lambda,
      linkinv = linkinv, variance = variance, dev.resids = dev.resids, 
      aic = aic, mu.eta = mu.eta, initialize = initialize, 
      validmu = validmu, valideta = valideta, simulate = simulate, 
      DlogLDmu = DlogLDmu, D2logLDmu2 = D2logLDmu2, D3logLDmu3 = D3logLDmu3, 
      D2muDeta2 = D2muDeta2, D3muDeta3 = D3muDeta3, resid.model=resid.model, 
      flags =list(obs=TRUE, exp=TRUE, canonicalLink=(linktemp=="loglambda"), LLgeneric=TRUE) 
    ), class = c("LLF","family")  
  )
}

# A lot of small function so that their variables are local
# Otherwise, they are those of the muetaenve, and nested loops interfere with each other.
# Which is also why there are distinct EX_it, EX2_it, EX3_it in the delayedAssign's

.CMP_asympto_EX <- function(pow_lam_nu, nu, c1) {
  pow_lam_nu*( 1 -(nu-1)/(2*nu*pow_lam_nu) - c1/((nu*pow_lam_nu)^2)- c1/((nu*pow_lam_nu)^3))
}

.CMP_asympto_Var <- function(pow_lam_nu, nu, c1) {
  (pow_lam_nu/nu)*( 1 + c1/((nu*pow_lam_nu)^2) + 2*c1/((nu*pow_lam_nu)^3))
}

.CMP_asympto_k3 <- function(pow_lam_nu, nu, c1) {
  (pow_lam_nu/nu^2)*( 1 - c1/((nu*pow_lam_nu)^2) -4*c1/((nu*pow_lam_nu)^3))
}

.CMP_asympto_k4 <- function(pow_lam_nu, nu, c1) {
  (pow_lam_nu/nu^3)*( 1 + c1/((nu*pow_lam_nu)^2) +8*c1/((nu*pow_lam_nu)^3))
}


.CMP_series_EX <- function(lambda, nu, denum_Z, denum_corr) {
  num <- .COMP_Z_n(lambda=lambda,nu=nu)
  uncorr <- .COMP_Z_ratio(num,denum_Z)
  # continuity correction wrt poisson: corr =lambda + error of the COMP_Z... functions
  corr <- .COMP_Z_ratio(.COMP_Z_n(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
  uncorr+(lambda-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
}

.CMP_series_EX2 <- function(lambda, nu, denum_Z, denum_corr) {
  # cat(crayon::red(EX2_it)," ")
  num <- .COMP_Z_n2(lambda=lambda,nu=nu)
  uncorr <- .COMP_Z_ratio(num,denum_Z)
  # cotinuity correction wrt poisson: 
  corr <- .COMP_Z_ratio(.COMP_Z_n2(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
  uncorr+(lambda*(1+lambda)-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
}

.CMP_series_EX3 <- function(lambda, nu, denum_Z, denum_corr) {
  num <- .COMP_Z_n3(lambda=lambda,nu=nu)
  uncorr <- .COMP_Z_ratio(num,denum_Z)
  # cotinuity correction wrt poisson: 
  corr <- .COMP_Z_ratio(.COMP_Z_n3(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
  uncorr + (lambda*(1+lambda*(3+lambda))-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
}

.CMP_series_EX4 <- function(lambda, nu, denum_Z, denum_corr) {
  num <- .COMP_Z_n4(lambda=lambda,nu=nu)
  uncorr <- .COMP_Z_ratio(num,denum_Z)
  # cotinuity correction wrt poisson: 
  corr <- .COMP_Z_ratio(.COMP_Z_n4(lambda=lambda,nu=1), denum_corr) ## poisson value by general approx
  uncorr + (lambda*(1+lambda*(7+lambda*(6+lambda)))-corr) ## approx_any_nu+(exact_poi-approx_poi): exact in nu=1
}

.CMP_muetaenv <- function(family, pw, eta) {
  # cat(crayon::bgRed("NEW muetaenv"))
  EX <- uniqEX <- c1 <- denum_Z <- uniqdenum_Z <-lambdas <- uniqlambdas <- uniqpow_lams_nu <- this <- use_asympto <- 
    uniquse_asympto <- Vmu <- dmudeta <- mu <- sane_eta <- etamatch <- uniq_asympto_var <- uniq_asympto_k3 <- NULL 
  nu <- environment(family$aic)$nu
  uniqeta <- unique(eta)
  uniqlen <- length(uniqeta)
  len <- length(eta)
  muetaenv <- list2env(list(pw=eval(pw), 
                            family=family, 
                            sane_eta=eta, 
                            uniqeta=uniqeta,
                            uniqlen=uniqlen,
                            len=len,
                            etamatch=match(eta,uniqeta),
                            nu=nu,
                            c1=(nu^2-1)/24,
                            denum= vector("list", len),
                            denum_corr= vector("list", len), # each element being a vector with elements (logScaleFac, scaled)
                            uniqdenum= vector("list", uniqlen),
                            uniqdenum_corr= vector("list", uniqlen), # each element being a vector with elements (logScaleFac, scaled)
                            Var= numeric(len)
  ), parent=environment(.muetafn))
  muetaenv$this <- muetaenv
  delayedAssign("pow_lams_nu", {
    # cat(crayon::bgRed("pow_lams_nu"))
    lambdas^(1/nu)}, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("uniqpow_lams_nu", {
    # cat(crayon::bgRed("pow_lams_nu"))
    uniqlambdas^(1/nu)}
    , assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("uniquse_asympto", {
    uniquse_asympto <- logical(uniqlen)
    for (it_ua in seq_len(uniqlen)) uniquse_asympto[[it_ua]] <- .CMP_use_asympto(uniqpow_lams_nu[[it_ua]],nu, uniqlambdas[[it_ua]])
    uniquse_asympto
  }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("use_asympto", { uniquse_asympto[etamatch] }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("uniqdenum_Z", {
    # cat(crayon::bgRed("denum_Z"))
    for (den_it in seq_len(uniqlen)) {
      if (uniquse_asympto[[den_it]]) {
        # that should not be used in this case
        uniqdenum[[den_it]] <- .Rcpp_COMP_Z_asympto(nu, uniqpow_lams_nu[[den_it]])
        uniqdenum_corr[[den_it]] <- "denum_corr sought despite use_asympto"
      } else {
        uniqdenum[[den_it]] <- .COMP_Z(lambda=uniqlambdas[[den_it]],nu=nu)
        uniqdenum_corr[[den_it]] <- .COMP_Z(lambda=uniqlambdas[[den_it]],nu=1)
      }
    }
    uniqdenum
  }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("denum_Z", { uniqdenum_Z[etamatch] }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("uniqEX", {
    # cat(crayon::bgRed("EX"))
    uniqEX <- numeric(uniqlen)
    for (EX_it in seq_len(uniqlen)) {
      if (uniqlambdas[[EX_it]]==0) {
        uniqEX[[EX_it]] <- 0
      } else {
        if (uniquse_asympto[[EX_it]]) {
          # using Gaunt et al.:
          uniqEX[[EX_it]] <- .CMP_asympto_EX(pow_lam_nu=uniqpow_lams_nu[[EX_it]], nu=nu, c1=c1)
        } else {
          uniqEX[[EX_it]] <-  .CMP_series_EX(lambda=uniqlambdas[[EX_it]], nu=nu, denum_Z=uniqdenum_Z[[EX_it]], 
                                             denum_corr=uniqdenum_corr[[EX_it]])
        }
      }
    } 
    uniqEX
  }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("EX", { uniqEX[etamatch] }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("uniq_asympto_var", {
    uniq_asympto_var <- numeric(uniqlen)
    for (var_it in seq_len(uniqlen)) {
      if (uniqlambdas[[var_it]]>0 && uniquse_asympto[[var_it]]) {      
        uniq_asympto_var[[var_it]] <- .CMP_asympto_Var(pow_lam_nu=uniqpow_lams_nu[[var_it]], nu=nu, c1=c1) 
      }
    } 
    uniq_asympto_var
  }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("uniq_asympto_k3", {
    uniq_asympto_k3 <- numeric(uniqlen)
    for (k3_it in seq_len(uniqlen)) {
      if (uniqlambdas[[k3_it]]>0 && uniquse_asympto[[k3_it]]) {      
        uniq_asympto_k3[[k3_it]] <- .CMP_asympto_k3(pow_lam_nu=uniqpow_lams_nu[[k3_it]], nu=nu, c1=c1)
      }
    } 
    uniq_asympto_k3
  }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("EX2", {
    # cat(crayon::bgRed("EX2"))
    uniqEX2 <- numeric(uniqlen)
    for (EX2_it in seq_len(uniqlen)) {
      if (uniqlambdas[[EX2_it]]==0) {
        uniqEX2[[EX2_it]] <- 0
      } else {
        if (uniquse_asympto[[EX2_it]]) {      
          vari <- uniq_asympto_var[[EX2_it]]
          EXi <- uniqEX[[EX2_it]]
          uniqEX2[[EX2_it]] <- EXi*EXi + vari # supp to mu1^2 when nu>1
        } else {
          uniqEX2[[EX2_it]] <- .CMP_series_EX2(lambda=uniqlambdas[[EX2_it]], nu=nu, denum_Z=uniqdenum_Z[[EX2_it]], 
                                               denum_corr=uniqdenum_corr[[EX2_it]])
        }
      }
    } 
    EX2 <- uniqEX2[etamatch]
    EX2
  }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("EX3", {
    uniqEX3 <- numeric(uniqlen)
    for (EX3_it in seq_len(uniqlen)) {
      if (uniqlambdas[[EX3_it]]==0) {
        uniqEX3[[EX3_it]] <- 0
      } else {
        if (uniquse_asympto[[EX3_it]]) {    
          k3 <-  uniq_asympto_k3[[EX3_it]]
          vari <- uniq_asympto_var[[EX3_it]]
          EXi <- uniqEX[[EX3_it]]
          uniqEX3[[EX3_it]] <- EXi*EXi*EXi + 3*EXi*vari + k3
        } else {
          uniqEX3[[EX3_it]] <- .CMP_series_EX3(lambda=uniqlambdas[[EX3_it]], nu=nu, denum_Z=uniqdenum_Z[[EX3_it]], 
                                               denum_corr=uniqdenum_corr[[EX3_it]])
        }
      }
    } 
    EX3 <- uniqEX3[etamatch]
    EX3
  }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("EX4", {
    uniqEX4 <- numeric(len)
    for (EX4_it in seq_len(uniqlen)) {
      if (uniqlambdas[[EX4_it]]==0) {
        uniqEX4[[EX4_it]] <- 0
      } else {
        if (uniquse_asympto[[EX4_it]]) {      
          k4 <-  .CMP_asympto_k4(pow_lam_nu=uniqpow_lams_nu[[EX4_it]], nu=nu, c1=c1)
          k3 <-  uniq_asympto_k3[[EX4_it]]
          vari <- uniq_asympto_var[[EX4_it]]
          EXi <- uniqEX[[EX4_it]]
          EXi_sq <- EXi*EXi
          uniqEX4[[EX4_it]] <- EXi_sq*EXi_sq + 6*vari*EXi_sq + 3*vari^2 + 4*k3*EXi + k4
        } else {
          uniqEX4[[EX4_it]] <- .CMP_series_EX4(lambda=uniqlambdas[[EX4_it]], nu=nu, denum_Z=uniqdenum_Z[[EX4_it]],
                                               denum_corr=uniqdenum_corr[[EX4_it]])
        }
      }
    } 
    EX4 <- uniqEX4[etamatch]
    EX4
  }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("Vmu", {family$variance(mu, muetaenv=this)}, assign.env = muetaenv, eval.env = muetaenv)
  
  delayedAssign("GLMweights", {
    GLMweights <- pw * dmudeta^2 /Vmu
    attr(GLMweights, "unique") <- FALSE ## might actually be true sometimes
    attr(GLMweights, "is_unit") <- FALSE
    GLMweights
  }, assign.env = muetaenv, eval.env = muetaenv)
  
  if (family$link=="loglambda") {
    muetaenv$uniqlambdas <- exp(uniqeta)
    muetaenv$lambdas <- muetaenv$uniqlambdas[muetaenv$etamatch]
    muetaenv$mu <- family$linkinv(eta, muetaenv=muetaenv) # calls muetaenv$EX hence must be after its definition
    delayedAssign("dmudeta", {family$mu.eta(sane_eta, muetaenv=this)}, assign.env = muetaenv, eval.env = muetaenv)
  } else {
    muetaenv$mu <- family$linkinv(eta)
    delayedAssign("lambdas", { family$mu2lambda(mu) }, assign.env = muetaenv, eval.env = muetaenv)
    delayedAssign("uniqlambdas", {unique(lambdas)}, assign.env = muetaenv, eval.env = muetaenv) #TEMPO
    delayedAssign("dmudeta", {family$mu.eta(sane_eta)}, assign.env = muetaenv, eval.env = muetaenv)
  }
  
  delayedAssign("thetaMuDerivs_2", { # which are derivatives of loglambda needed for fixed-effect models
    .CMP_thetaMuDerivs_2(muetaenv=this)
  }, assign.env = muetaenv, eval.env = muetaenv)
  delayedAssign("thetaMuDerivs_3", { # which are derivatives of loglambda needed for Mixed-effect models
    .CMP_thetaMuDerivs_3(muetaenv=this)
  }, assign.env = muetaenv, eval.env = muetaenv)
  return(muetaenv)
}

