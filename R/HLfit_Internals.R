
## declared "arglist" to print a clean summary instead of very long list
## the print(summary()) is required if the "arglist" is a member (eg as.list(<call>));
## summary alone would print nothing
print.arglist <-function(x,...) {print(summary(x,...))}
##

.sweepZ1Wwrapper <- function(ZZ,WW) { ## for *m*atrix input
  if (nrow(ZZ)!=length(WW)) {
    stop("From .sweepZ1Wwrapper(): nrow(ZZ)!=length(WW) ") ## fatal error for eigen code...
  } else if (ncol(ZZ)==0L) {
    return(ZZ)
  } else sweepZ1W(ZZ,WW) ##Rcpp
}

## the following fns try to keep the input class in output, but are called with dense matrices (except irst tested case).
# les Matrix::(t)crossprod  paraissent (parfois au moins) remarquablement inefficaces !!
# idem pour Diagonal()
.ZWZtwrapper <- function(ZAL,w) { ## in .AUGI0_ZX_sparsePrecision; in calc_asDmLR_invV_from_fitobject...
  if (inherits(ZAL,"Matrix")) {
    if (inherits(ZAL,"ddiMatrix")) {
      if ( ZAL@diag=="U") {
        return(Diagonal(x=w))
      } else {
        ZAL@x <- ZAL@x * w ## OK for diagonal Z
        return(ZAL)
      }
    } else  {
      return(Matrix::tcrossprod( Matrix_times_Dvec(ZAL,w) ,ZAL)) ## (Z W) %*% Zt
    } 
  } else return(ZWZt(ZAL,w))
}

.ZtWZwrapper <- function(ZAL,w) { ## used in seval contexts
  if (ncol(ZAL)==0L) {
    stop(".ZtWZwrapper called with ncol(ZAL)=0") ## temporary devel code since all calls are in principle protected 
  } else if (inherits(ZAL,"Matrix")) {
    if (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U") {
      return(Diagonal(x=w))
    } else  {
      DZAL <- ZAL
      DZAL@x <- DZAL@x * w[DZAL@i+1L] ## W Z
      return(Matrix::crossprod(ZAL, DZAL)) ## t(Z) %*% (W Z)
    }
  } else return(ZtWZ(ZAL,w))
}

Sigwrapper <- function(ZAL,wa,wb,ZALtZAL=NULL) { 
  if (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U") {
    Sig <- diag( wa + wb )
  } else if (! is.null(ZALtZAL)) { ## FR->FR currently always FALSE bc
    # only for constant w.ranef (LMM with single ranef (CAR?)) this would be useful
    Sig <- ZALtZAL*wa[1] ## FR->FR valid only for constant wa; this is why it is not currently used
    nc <- ncol(Sig)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    Sig[diagPos] <- Sig[diagPos] + wb 
  } else if (inherits(ZAL,"Matrix")) {
    Sig <- Matrix::tcrossprod(ZAL %*% Diagonal(x=wa),ZAL)
    nc <- ncol(Sig)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    Sig[diagPos] <- Sig[diagPos] + wb 
  } else { ## valid for matrix but not for Matrix
    Sig <- Rcpp_Sig(ZAL,wa,wb)
  }
  return(Sig)
}


calcRanefPars <- function(corrEstList=NULL, ## potentially a list...
                          lev_lambda,
                          ranefEstargs,
                          lambda.Fix,
                          rand.families,
                          lcrandfamfam,
                          psi_M,
                          verbose,
                          control,
                          iter ## ajustement gamma(identity...)
) {
  ## Build pseudo response for lambda GLM/HGLM
  glm_lambda <- NULL
  next_LMatrices <- NULL
  LMatricesList <- ranefEstargs$prev_LMatrices
  if ( ! is.null(LMatricesList) && ! is.list(LMatricesList)) LMatricesList <- list(dummyid=LMatricesList)
  if ( ! is.null(corrEstList) && ! is.list(corrEstList)) corrEstList <- list(corr_est=corrEstList)
  next_corrEstList <- list()
  ranefs <- attr(ranefEstargs$ZAlist,"ranefs")
  nrand <- length(ranefEstargs$ZAlist) ## Cf notes du 3/6/2015
  ## Male/Female has one apparent (.|.) but ZA is duplicated in ZAlist, with ranefs evaluating to "( 1 | Female )" "( 1 | Male )" 
  ## so that nrand = length(ranefs) = length(ranefEstargs$ZAlist)
  done <- rep(FALSE,nrand)
  u_h <- ranefEstargs$u_h
  cum_n_u_h <- ranefEstargs$cum_n_u_h
  resp_lambda <- matrix(0,cum_n_u_h[nrand+1L],1L)
  #########################
  isRandomSlope <- (attr(ranefEstargs$ZAlist,"Xi_cols") > 1L)
  if (any(isRandomSlope)) {
    ## handling correlation in random slope models # slmt pr gaussian ranefs, verif dans preprocess
    ranefEstargs$lcrandfamfam <- lcrandfamfam
    LMatricesBlob <- do.call(.spaMM.data$options$covEstmethod, ranefEstargs)            
    next_LMatrices <- LMatricesBlob$next_LMatrices ## updated list of matrices where only random-slope elements have been updated
    next_lambda_est <- LMatricesBlob$next_lambda_est ## a full-length vector with values only in the appropriate u ranges 
    ## only for testing convergence: 
    next_corrEstList$cov_est_vec <- LMatricesBlob$optr_par
    for (it in which(isRandomSlope)) {
      u.range <- (cum_n_u_h[it]+1L):cum_n_u_h[it+1L]
      resp_lambda[u.range] <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=1) ## must give d1 in table p 989 de LeeN01
    }
    done[ isRandomSlope ] <- TRUE
  } else { ## only 'declarations' for all further code
    next_lambda_est <- numeric(length(u_h)) ## next_LMatrices remains an empty list()
    next_corrEstList$cov_est_vec <- NULL
  }
  ### next the other LMatrix models
  for (it in seq_along(LMatricesList) ) { ## this loop will ignore ranefs not affected by any lmatrix
    lmatrix <- LMatricesList[[it]]
    ## find ZAlist elements affected by LMatrix element
    affected <- which(ranefs %in% attr(lmatrix,"ranefs") & ! done)
    ## then for each L matrix we need to select the relevant blocks of random effects
    if (length(affected)>1L) {
      stop("code needed for length(affected)>1")
    } else if (length(affected)==1L) {
      if ( ! is.null(corrEstList$corr_est) && attr(lmatrix,"corr.model") %in% c("adjacency","ar1")) { 
        ## the following conforms to an interface where 
        ##  next_lambda_est is a vector (heterosc) of what should go in the sigma_aug matrix
        ## consequently, the LMatrix for this model should be decomp$u, constant 
        next_LMatrices[it] <- LMatricesList[it]
        adj_symSVD <- ranefEstargs$processed$AUGI0_ZX$envir$adj_symSVD ## may be NULL
        if (is.null(adj_symSVD)) {
          stop("is.null(adj_symSVD)")
          adj_symSVD <- attr(lmatrix,attr(lmatrix,"type"))  ## older conception
        }
        adjd <- adj_symSVD$adjd
        locdf <- data.frame(adjd=adjd) ## $adjd, not $d which is (1/(1-rho * $adjd)): adj, not corr
        u.range <- (cum_n_u_h[affected]+1L):cum_n_u_h[affected+1L]
        locdf$resp <- resp_lambda[u.range] <- u_h[u.range]^2
        ## here CAR allows REML contrary to the SEM CAR, hence leverages
        glm_lambda <- calc_CARdispGammaGLM(data=locdf, lambda.Fix=lambda.Fix[affected], lev=lev_lambda[u.range],control=control)
        attr(glm_lambda,"whichrand") <- affected
        next_lambda_est[u.range] <- fitted(glm_lambda) ## prediction of heteroscedastic variances
        coeffs <- coefficients(glm_lambda)
        if (is.na(lambda.Fix[affected])) {  ## pour l'instant tjrs vrai mais évoluera
          next_corrEstList$corr_est <- list(rho = - coeffs[["adjd"]]/ coeffs[["(Intercept)"]]) ## FR->FR different corr_est s'écrasent les uns les autres
        } else {
          next_corrEstList$corr_est <- list(rho = - coeffs[1]*lambda.Fix)  
        }
        done[affected] <- TRUE
      } ## Matern remains undone
    } ## else (length(affected)=0), for previously processed random slope models
  }        
  ### next the (no L matrices) or (Matern model and other fixed L matrix cases)
  for (it in which( ! done )) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    if (is.na(unique.lambda <- lambda.Fix[it])) {
      resp_lambda[u.range] <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=1) ## must give d1 in table p 989 de LeeN01
      unique.lambda <- sum(resp_lambda[u.range])/sum(1-lev_lambda[u.range]) ## NOT in linkscale 
      unique.lambda <- max(unique.lambda,1e-8) # FR->FR still corrected
      unique.lambda <- min(unique.lambda,.spaMM.data$options$maxLambda)  
      if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { ## Gamma(identity)
        unique.lambda <- pmin(unique.lambda,1-1/(2^(iter+1)))  ## impose lambda<1 dans ce cas 
      }
    } 
    next_lambda_est[u.range] <- rep(unique.lambda,length(u.range))
  }
  if (verbose["trace"]) { print(paste("lambda(s)=",paste(signif(unique(next_lambda_est),4),collapse=" ")),quote=F)}
  return(list(next_LMatrices=next_LMatrices,
              resp_lambda=resp_lambda, ## for final glm...
              next_corrEstList=next_corrEstList,
              next_lambda_est=next_lambda_est, ## heterosc
              glm_lambda=glm_lambda ## potentially to be replaced by a list of glms later
              ))
}

process_resglm_list <- function(resglm_lambdaS, ## les 2 autres args for handling errors
                                nrand) {
  lambda_seS <- as.list(rep(NA,nrand)) # return value will be a list of length nrand
  coefficients_lambdaS <- as.list(rep(NA,nrand)) ## idem
  linkS <- list() # list of same length as resglm_lambdaS
  linkinvS <- list() # list of same length as resglm_lambdaS 
  warnmesses <- list() 
  for (glmit in seq_len(length(resglm_lambdaS))) {
    glm_lambda <- resglm_lambdaS[[glmit]]
    # next line for SEM
    coeff_lambdas <- coefficients(glm_lambda)
    coeffs_substitute <- glm_lambda$coeffs_substitute ## convoluted but avoids permutations by using the names:
    ## code for SEM 09/2015:
    ##The coefficients have different names whether a precomputed design matrix was used in ~X-1 or the data=<data.frame> syntax
    if (class(glm_lambda$data)=="environment") {
      locform <- glm_lambda$formula
      if (deparse(locform[[length(locform)]]) == "X - 1") {
        names(coeff_lambdas) <- colnames(glm_lambda$model$X)
      } else stop(paste("code missing for names(coeff_lambdas) with formula:",deparse(glm_lambda$formula)))
    } ## else names are already OK (and model$X does not exist)
    ## FR->FR mais je n'ai encore aucun controle sur l'identite des noms de params
    if ( ! is.null(coeffs_substitute)) coeff_lambdas <- coeffs_substitute[names(coeff_lambdas)]
    coefficients_lambdaS[[attr(glm_lambda,"whichrand")]] <- coeff_lambdas
    #
    p_lambda <- length(coeff_lambdas) ## only for the local resglm
    lambda_seS[[attr(glm_lambda,"whichrand")]] <- summary(glm_lambda,dispersion=1)$coefficients[(p_lambda+1L):(2L*p_lambda)]
    linkS[[glmit]] <- glm_lambda$family$link
    linkinvS[[glmit]] <- glm_lambda$family$linkinv
    rf <- attr(resglm_lambdaS,"rand.families")
    for (it in seq_len(length(rf))) { ## iteration over ranefs only for the local resglm
      if (tolower(rf[[it]]$family)=="gamma" && rf[[it]]$link=="identity" && coeff_lambdas[it]>0) {
        message("lambda failed to converge to a value < 1 for gamma(identity) random effects.")
        message("This suggests that the gamma(identity) model with lambda < 1 is misspecified, ")
        message("and Laplace approximations are unreliable for lambda > 1. ")
      }
    }
    warnmesses[[glmit]] <- glm_lambda$warnmess
  }
  return( list(lambda_seS =lambda_seS, #list
             coefficients_lambdaS = coefficients_lambdaS, #list
             linkS = linkS, # list of same length as resglm_lambdaS
             linkinvS = linkinvS, # list of same length as resglm_lambdaS 
             warnmesses = warnmesses  ))
}

calcPHI <- function(oriFormula, ## with offset
                    dev.res,family,data,
                    lev_phi,
                    phimodel,verbose,method="glm",control.phi=list(),
                    control) {
  if (phimodel!="phiHGLM") {  
    if ( identical(deparse(oriFormula),"~1")) { ## one case where we can easily avoid an explicit call to a glm (but one will be used to compute SEs later) 
      next_phi_est <- sum(dev.res)/sum(1-lev_phi) ## NOT in linkscale
      beta_phi <- c("(Intercept)"=family$linkfun(next_phi_est)) ## linkscale value
      glm_phi <- NULL
    } else { ## which means that calcPHI cannot be used for final estim phi
      glm_phi <- calc_dispGammaGLM(formula=oriFormula, dev.res=dev.res,
                                   data=data,lev=lev_phi, family=family,
                                   etastart=control.phi$etastart, ## private, availabble, but NULL so far
                                   control=control)
      beta_phi <- coefficients(glm_phi)
      next_phi_est <- fitted(glm_phi)
      if (family$link!="log" && any(next_phi_est<=0)) { stop("Gamma GLM for dispersion yields negative phi estimates.") }
    }
    if (verbose["trace"]) {print(paste("phi_est=",signif(next_phi_est,4)),quote=F)}
  } else { ## random effect(s) in predictor for phi
    stop("This function should not be called when a residual model with random effects is fitted.")
  } 
  return(list(next_phi_est=next_phi_est,  #low phi values are handled in calc_APHLs...
              glm_phi=glm_phi,
              beta_phi=beta_phi ## used at least to initiate final GLM in "~1" case
  )) ## 
}


`calc.w.resid` <- function(GLMweights,phi_est) { ## One should not correct this phi_est argument by prior.weights (checked)
  phi_est[phi_est<1e-12] <- 1e-11 ## 2014/09/04 local correction, cf comment in calc_APHLS...
  structure(as.vector(GLMweights/phi_est),unique= (attr(GLMweights,"unique") && length(phi_est)==1L))
}


## spaMM_Gamma() fixes Gamma()$dev.resids(1e10+2,1e10,1) is < 0
# dev.resids() must be >0 for coputation deviance_residual in fitting Gamma GLMM, and alos for $aic() computation.
spaMM_Gamma <- function (link = "inverse") {
  mc <- match.call()
  linktemp <- substitute(link) ## does not evaluate
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp) ## converts to char the unevaluated expression
  okLinks <- c("inverse", "log", "identity")
  if (linktemp %in% okLinks) 
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link) ## evals expression converted to char (with  $name in particular); but returns link=linktemp, not link=stats$name
    # problem is that the families fns return link=linktemp, which seems weird: better is   
    linktemp <- stats$name ## line not in Gamma() [and different in binomial()], which absence prevents programming with link argument...  
  } else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name)) 
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for gamma family; available links are %s", 
                    linktemp, paste(sQuote(okLinks), collapse = ", ")), 
           domain = NA)
    }
  }
  variance <- function(mu) mu^2
  validmu <- function(mu) all(mu > 0)
  dev.resids <- function(y, mu, wt) {
    if (any(mu < 0)) return(Inf) ## 2015/04/27; maybe not useful
    ## otherwise see deviance.gamma function locally defined in statmod::glmgam.fit
    dev_res <- -2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
    dev_res[dev_res<.Machine$double.eps] <- .Machine$double.eps ##  ## FR: added this
    return(dev_res)
  }
  aic <- function(y, n, mu, wt, dev) {
    n <- sum(wt)
    disp <- dev/n
    -2 * sum(dgamma(y, 1/disp, scale = mu * disp, log = TRUE) * 
               wt) + 2
  }
  initialize <- expression({
    if (any(y <= 0)) stop("non-positive values not allowed for the gamma family") 
    n <- rep.int(1, nobs)
    mustart <- y
  })
  simfun <- function(object, nsim) {
    wts <- object$prior.weights
    if (any(wts != 1)) 
      message("using weights as shape parameters")
    ftd <- fitted(object)
    shape <- MASS::gamma.shape(object)$alpha * wts 
    resu <- rgamma(nsim * length(ftd), shape = shape, rate = shape/ftd)
    if (nsim>1L) resu <- matrix(resu,ncol=nsim)
    resu
  }
  # linkinv <- function (eta) pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax) 
  # : permet des plantages severes dans glm.fit ensuite (R CMD check en detecte) 
  ## all closures defined here have parent.env the environment(spaMM_Gamma) ie <environment: namespace:spaMM>
  ## changes the parent.env of all these functions (aic, dev.resids, simfun, validmu, variance): 
  # as.list(environment(aic)) ## this has an unexplained effet on saveSize!
  parent.env(environment(aic)) <- environment(stats::Gamma) ## parent = <environment: namespace:stats>
  ## That _does_ reduce the size of the fitted objects using spaMM_Gamma (eg in a phi.object)
  ## That does not eliminate an hidden environment shared among member functions 
  #    _after_ compiling the package:  
  ##  spaMM:::.saveSize(attr(attr(spaMM_Gamma()$aic,"srcref"),"srcfile")) grows
  ## compared to the non-compiled version
  ## It _is_ an environment : ls(attr(attr(spaMM_Gamma()$aic,"srcref"),"srcfile")) lists it.  
  ## But its size is not explained by its contents...
  structure(list(family =  structure("Gamma",patch="spaMM_Gamma"), 
                 link = linktemp, linkfun = stats$linkfun, 
                 linkinv = stats$linkinv, 
                 variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun), 
            class = "family")
}


`selectLoglfn` <- function(family) {
  famfam <- tolower(family$family)
  # return value of each fn must be a vector if y is a vector
  switch(famfam,
         gaussian = function(theta,y,nu) {nu*(theta*y-(theta^2)/2)- ((y^2)*nu+log(2*pi/nu))/2}, 
         poisson = function(theta,y,nu) { ## theta = log(mu)
           res <- nu*(theta*y-attr(theta,"mu"))   -  lfactorial(y)
           res[theta== -Inf & y==0] <- 1
           res
         },
         binomial = function(theta,freqs,sizes,nu) {nu*sizes*(freqs*theta-log(1+exp(theta))) +lchoose(sizes,round(sizes*freqs))},
         # gamma = function(theta,y,nu) {nu*(y*theta+log(-theta))+nu*(log(nu*y))-lgamma(nu)-log(y)} ## mean mu=-1/th, **** var = mu^2 / vu ****
         # same by using ad hoc C library...
         gamma = function(theta,y,nu) {
           dgamma(y, shape=nu , scale = attr(theta,"mu")/nu, log = TRUE) 
         },
         compoisson = function(theta,y,nu) { ## theta = log(lambda) 
           COMP_nu <- environment(family$aic)$nu
           logLs <- sapply(seq(length(y)), function(i) {
             comp_z <- .COMP_Z(lambda=exp(theta[i]),nu=COMP_nu)
             res <- nu[i]*(theta[i]*y[i]-comp_z[[1]]-log(comp_z[[2]])) - COMP_nu * lfactorial(y[i])
             res
           })
           logLs[theta== -Inf & y==0] <- 1
           logLs
         },
         negbin = function(theta,y,nu) { ## theta is the canonical param, -log(1+shape/mu)
           NB_shape <- environment(family$aic)$shape
           ( nu*(theta*y-NB_shape*log(1+attr(theta,"mu")/NB_shape))
             + lgamma(y+NB_shape) - lgamma(y + 1) - lgamma(NB_shape))
           # doit etre -negbin()$aic(y=y,n,mu=NB_shape/(exp(-theta)-1),wt=1,dev)/2
         }
  )
}


`.theta.mu.canonical` <- function(mu,family) { 
  ## the (fixed) canonical link between theta and mu, not the family link between eta and mu 
  if (inherits(family,"family")) {
    famfam <- family$family
  } else famfam <- family
  switch(tolower(famfam),
         gaussian = mu ,
         poisson = structure(log(mu),mu=mu) ,
         binomial = make.link("logit")$linkfun(mu),  # correct syntax, does not use 'non-public API' such as .Call to access code from dlls from the base packages...
         ## if this does no work, use 
         #                 { 
         #                    theta <- logit(mu)
         #                    theta[theta>27.6310] <- 27.6310 ## mu>1-1e-12
         #                    theta[theta < -27.6310] <- -27.6310 
         #                    theta
         #                 },
         gamma = structure(-1/mu,mu=mu), ## "-": McC&N p. 290
         compoisson = {
           if (is.null(lambda <- attr(mu,"lambda"))) {
             lambda <-  family$linkfun(mu,log=FALSE)
           }  
           structure(log(lambda),mu=mu)
         },
         negbin = structure(-log(1+(environment(family$aic)$shape)/mu),mu=mu) ## keep mu b/c useful for loglfn.fix
  )
} ## returns values for given mu

`theta.mu.conjugate` <-function(mu,family) { 
  ## theta(u) in LeeN01... this is more pedagogy than efficient code
  switch(tolower(family),
         gaussian = .theta.mu.canonical(mu,"gaussian") , ## mu 
         gamma = .theta.mu.canonical(mu,"poisson"), ## log(mu)
         beta = .theta.mu.canonical(mu,"binomial"), ## improved logit(mu)      
         "inverse.gamma" = .theta.mu.canonical(mu,"gamma") ## -1/mu
  )
} ## returns values for given mu

## another approach is
#theta.mu.call <- switch(tolower(family),
#      poisson =  call("log",quote(mu)) ,
#   )
## followed by eval(theta.mu.call)

### there is a nice syntax
#eta.mu.expr <- parse(text=paste(family$link,"(mu)",sep="")) # "log" -> expression(log(mu))
#etadmu<-D(theta.mu.expr,"mu") ## a call; eval(dthetadmu) then suffices if mu has a value
### except that D(identity... does not work)

.thetaMuDerivs <-function(mu,BinomialDen,family) { ## used for non-canonical links
  familyfam <- family$family
  if (familyfam=="binomial") muFREQS <- mu/BinomialDen
  if (familyfam=="negbin") NB_shape <- environment(family$aic)$shape
  ## these definitions depend only on the canonical link
  Dtheta.Dmu <- switch(tolower(familyfam),
                       gaussian = rep(1,length(mu)) ,
                       poisson = 1/mu ,
                       binomial = 1/(muFREQS*(1-muFREQS)),
                       gamma = 1/mu^2,
                       negbin = 1/(mu*(1+mu/NB_shape))
                       # COMPoisson has no implemented non-canonical link
  ) ## values for given mu
  if (familyfam=="binomial") Dtheta.Dmu <- Dtheta.Dmu/BinomialDen
  D2theta.Dmu2 <- switch(tolower(familyfam),
                         gaussian = rep(0,length(mu)) ,
                         poisson = -1/mu^2 ,
                         binomial = -(1-2*muFREQS)/(muFREQS*(1-muFREQS))^2,
                         gamma = -2/mu^3,
                         negbin = -(1+2*mu/NB_shape)/(mu*(1+mu/NB_shape))^2
                         # COMPoisson has no implemented non-canonical link
  ) ## values for given mu
  if (familyfam=="binomial") D2theta.Dmu2 <- D2theta.Dmu2/(BinomialDen^2)
  return(list(Dtheta.Dmu=Dtheta.Dmu,D2theta.Dmu2=D2theta.Dmu2))
}

muetafn <- function(eta,BinomialDen,processed) { ## note outer var BinomialDen 
  family <- processed$family
  ## patches for borderline eta's
  if (family$link =="log") {
    eta[eta>30] <-30 ## 100 -> mu = 2.688117e+43 ; 30 -> 1.068647e+13
  } else if (family$family == "COMPoisson" && family$link =="loglambda") {
    etamax <- 30*environment(family$aic)$nu
    eta[eta>etamax] <- etamax ## using log(mu) ~ eta/nu for large nu
  } else if (family$link=="inverse" && family$family=="Gamma") {
    etamax <- sqrt(.Machine$double.eps)
    eta[eta>etamax] <- etamax ## both eta and mu must be >0
  }
  mu <- family$linkinv(eta) ## linkinv(eta) is FREQS for binomial, COUNTS for poisson...
  if (family$link %in% c("logit","probit","cloglog","cauchit")) {
    mu[mu > (1-1e-12)] <- (1-1e-12)
    mu[mu < (1e-12)] <- (1e-12)
  }
  dmudeta <- family$mu.eta(eta) ## aberrant at hoc code for cloglog 'elsewhere'...
  Vmu <- family$variance(mu) 
  if (family$family=="binomial") {
    Vmu <- Vmu * BinomialDen 
    mu <- mu * BinomialDen
    dmudeta <- dmudeta * BinomialDen
  } 
  if (processed$LMMbool) {
    GLMweights <- eval(processed$prior.weights) ## with attr(.,"unique")
    attr(GLMweights,"unique") <- attr(processed$prior.weights,"unique") ## might actually be true sometimes
  } else  if (family$family=="Gamma" && family$link=="log") {
    GLMweights <- eval(processed$prior.weights)  ## with attr(.,"unique")
    attr(GLMweights,"unique") <- attr(processed$prior.weights,"unique") ## might actually be true sometimes
  } else {
    GLMweights <- eval(processed$prior.weights) * dmudeta^2 /Vmu ## must be O(n) in binomial cases
    attr(GLMweights,"unique") <- FALSE ## might actually be true sometimes
  }
  return(list(mu=mu,dmudeta=dmudeta,GLMweights=GLMweights))
} ## end local def muetafn

.updateWranef <- function(rand.family,lambda,u_h,v_h) {
  dudv <- rand.family$mu.eta(v_h) ## general cf InvGamma with log link rand.family$mu.eta(v_h) = exp(v) =u is du/d(log u)   
  ## compute w.ranef := - d^2 log dens(v)/dv^2 := 1/Sigma^2_v (= 1/lambda for LMM). See Appendix 3 of LeeN01 + my notes
  ## computed either directly or as (dudv/V_M)*(dudv/lambda)
  ## compute dlogWran_dv_h := d log w.ranef/dv
  if (rand.family$family=="gaussian") {
    if (rand.family$link=="identity") {
      V_M <- rand.family$variance(u_h) ##rep(1,length(u_h)) ## GLMMs in general
      dlogWran_dv_h <- rep(0L,length(u_h))
    }
  } else if (rand.family$family=="Gamma") { 
    if (rand.family$link=="log") {
      V_M <- u_h ## V(u), canonical conjugate Gamma as in canonical Poisson Gamma HGLM
      dlogWran_dv_h <- rep(1L,length(u_h))
    } else if (rand.family$link=="identity") { ## gamma(identity)
      w.ranef <- as.numeric((1-lambda)/(lambda * u_h^2)) ## vanishes for lambda=1 and negative above... (in which case the Laplace approx is bad anyway)
      dlogWran_dv_h <- -2/as.numeric(u_h)
      return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))  ###### return here !
    } 
  } else if (rand.family$family=="inverse.Gamma") { ## for Gamma HGLM 
    ## the canonical form gives the density of theta(u)
    if (rand.family$link=="log") {
      w.ranef <- as.numeric(1/(u_h * lambda)) ## W1/lambda, W1 computation shown in appendix 3 of LeeN01; also in Noh and Lee's code.
      dlogWran_dv_h <- rep(-1L,length(u_h)) ## v=log u, dlogW/dv= dlog(1/u)/dv=-1
      return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))  ###### return here !
    } else if (rand.family$link=="-1/mu") {
      ## D[\[Nu] (\[Theta][u] - (-Log[-\[Theta][u]])), {\[Theta][u], 2}]
      V_M <- rand.family$variance(u_h) ## u_h^2 ## V(u), canonical conjugate HGLM 
      dlogWran_dv_h <- 2 * u_h ## no independent check 
    }
  } else if (rand.family$family=="Beta") {
    if (rand.family$link=="logit") {
      V_M <- rand.family$variance(u_h) ##  u_h*(1-u_h) ## canonical conjugate HGLM
      dlogWran_dv_h <- 1 - 2 * u_h ## D[Log[u (1 - u)] /. u -> 1/(1 + E^-v), v] /. v -> Log[u/(1 - u)] ; no independent check
    }
  }
  ## dudv/V_M may be 1 as both diverge: 
  w.ranef <- as.numeric((dudv/V_M)*(dudv/lambda)) ## semble valide quand v=g(u) = th(u): not previous return()
  return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))
}

.updateW_ranefS <- function(cum_n_u_h,rand.families,lambda,u_h,v_h) {
  nrand <- length(rand.families)
  blob <- lapply(seq(nrand), function(it) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    .updateWranef(rand.family=rand.families[[it]],lambda[u.range],u_h[u.range],v_h[u.range])
  })
  w.ranef <- unlist(lapply(blob,function(b) {b$w.ranef}))
  w.ranef[w.ranef>1e10] <- 1e10 ## patch useful to avoid singular d2hdv2 in PLoG model
  dlogWran_dv_h <- unlist(lapply(blob,function(b) {b$dlogWran_dv_h}))
  dvdu <- unlist(lapply(blob,function(b) {b$dvdu}))
  return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=dvdu))
}

d2mudeta2fn <- function(link,mu=NULL,eta=NULL,BinomialDen=NULL) { ## d2 MuCOUNTS d etaFREQS^2
  switch(link,
         identity = 0,
         log = mu, 
         inverse = 2 * mu^3 , ## canonical for Gamma()
         ## next three make sense for Binomial response data
         logit = {muFREQS <- mu/BinomialDen;
                  d2muFREQS <- muFREQS*(1-muFREQS)*(1-2*muFREQS);
                  d2muFREQS * BinomialDen
         },
         probit = -eta*dnorm(eta) * BinomialDen,
         cloglog = exp(eta-exp(eta))*(1-exp(eta)) * BinomialDen, ## D[1 - E^-E^\[Eta], {\[Eta], 2}]
         cauchit = -2 *eta/(pi * (1+eta^2)^2),
         stop(paste("unhandled link'",link,"'in d2mudeta2fn()"))
  )
} 

#derivatives of GLM weights wrt eta 
.calc_dlW_deta_locfn <- function(i,lambdas,mu,COMP_nu) {
  lambdai <- lambdas[[i]]
  mui <- mu[[i]]
  compz <- .COMP_Z(lambda=lambdai,nu=COMP_nu)
  compzn <- .COMP_Z_n(lambda=lambdai,nu=COMP_nu)
  compzn2 <- .COMP_Z_n2(lambda=lambdai,nu=COMP_nu)
  compzn3 <- .COMP_Z_n3(lambda=lambdai,nu=COMP_nu)
  rn3 <- .COMP_Z_ratio(compzn3,compz)
  rn2 <- .COMP_Z_ratio(compzn2,compz)
  dmu.dlogLambda <- rn2 - mui^2 # =family$mu.eta() without repeating some computations
  d2mu.dlogLambda2 <- rn3-mui*rn2-2*mui*dmu.dlogLambda
  #         if (is.nan(dmu.dlogLambda) ) { 
  #           if (COMP_nu<0.05 && lambdai<1) dmu.dlogLambda <- lambdai/(1-lambdai)^2
  #         } ## geometric approx
  #         if (is.nan(dmu.dlogLambda) ) { 
  #           if (COMP_nu<0.05 && lambdai<1) {dmu.dlogLambda <- lambdai*(1+lambdai)/(1-lambdai)^3
  #         } ## geometric approx
  return(c(dmudeta=dmu.dlogLambda, d2mudeta2=d2mu.dlogLambda2))
}

.calc_dlW_deta <- function(dmudeta,family,mu,eta,BinomialDen,canonicalLink,calcCoef1=FALSE) {
  coef1 <- NULL
  ## We first handle the canonical link cases, where comput. of coef1 depends only on the link  
  ## here w=dmudeta; d1=dwdmu dmudeta /w^2 = dlogwdeta/w = (d2mu/deta2)/(dmu/deta) /w =
  ##      (d2mu/deta2)/(dmu/deta)^2 = (d(dmudeta)/dmu)/dmudeta where d(dmudeta)/dmu is the numerator as detailed:
  if (canonicalLink) {
    #dlW_deta <- d2mudeta2 / dmudeta or :
    if (family$family=="gaussian") {
      if (calcCoef1) coef1 <- rep(0L,length(mu))
      dlW_deta <- rep(0L,length(mu))
    } else if (family$family=="poisson") {
      ## numerator is D[D[E^\[Eta], \[Eta]] /. {E^\[Eta] -> \[Mu]}, \[Mu]] =1 
      if (calcCoef1) coef1 <- 1/dmudeta
      dlW_deta <- rep(1L,length(mu))
    } else if (family$family=="binomial") {
      ## numerator is D[D[1/(1 + E^-\[Eta]), \[Eta]] /. {E^-\[Eta]->(1-\[Mu])/\[Mu]} ,\[Mu]]=1-2 mu 
      if (calcCoef1) coef1 <-(1-2*mu/BinomialDen)/dmudeta  
      dlW_deta <-(1-2*mu/BinomialDen)  
    } else if (family$family=="Gamma") { ## link= "inverse" !
      ## numerator is D[D[-1/\[Eta], \[Eta]] /. {\[Eta] -> -1/\[Mu]}, \[Mu]] =2 mu 
      if (calcCoef1) coef1 <- 2*mu /dmudeta
      dlW_deta <- 2*mu
    } else if (family$family=="COMPoisson") { 
      COMP_nu <- environment(family$aic)$nu
      lambdas <- exp(eta) ## pmin(exp(eta),.Machine$double.xmax) ##  FR->FR lambdas missing as mu attribute here ?
      blob <- sapply(seq(length(lambdas)), .calc_dlW_deta_locfn,lambdas=lambdas,mu=mu,COMP_nu=COMP_nu)
      dlW_deta <- blob["d2mudeta2",] / blob["dmudeta",]
      if (family$family=="COMPoisson") dlW_deta[is.nan(dlW_deta)] <- 0 ## quick patch for cases that should have low lik
      if (calcCoef1) {
        coef1 <- dlW_deta / blob["dmudeta",]
        if (family$family=="COMPoisson") coef1[is.nan(coef1)] <- 0 ## idem
      }
    } 
  } else if (family$family=="binomial" && family$link=="probit") { ## ad hoc non canonical case 
    muFREQS <- mu/BinomialDen
    dlW_deta <- -2*eta - dnorm(eta)*(1-2*muFREQS)/(muFREQS*(1-muFREQS))
    if (calcCoef1) {
      coef1 <- dlW_deta *(muFREQS*(1-muFREQS))/ (BinomialDen * dnorm(eta)^2) 
      coef1[coef1>1e100] <- 1e100
      coef1[coef1< -1e100] <- -1e100
    }
  } else if (family$family=="Gamma" && family$link=="log") { ## ad hoc non canonical case 
    if (calcCoef1) coef1 <- rep(0L,length(mu))
    dlW_deta <- rep(0L,length(mu)) ## because they both involve dW.resid/dmu= 0
  } else {
    ## we need to update more functions of mu...
    tmblob <- .thetaMuDerivs(mu,BinomialDen,family)
    Dtheta.Dmu <- tmblob$Dtheta.Dmu # calcul co fn de muFREQS puis / BinomialDen
    D2theta.Dmu2 <- tmblob$D2theta.Dmu2 # calcul co fn de muFREQS puis / BinomialDen ^2
    d2mudeta2 <- d2mudeta2fn(link=family$link,mu=mu,eta=eta,BinomialDen=BinomialDen)
    ## ... to compute this:
    D2theta.Deta2_Dtheta.Deta <- (D2theta.Dmu2 * dmudeta^2 + Dtheta.Dmu * d2mudeta2)/(Dtheta.Dmu * dmudeta)
    dlW_deta <- d2mudeta2 / dmudeta + D2theta.Deta2_Dtheta.Deta
    if (calcCoef1) coef1 <- dlW_deta / (Dtheta.Dmu * dmudeta^2) ## note that coef2 is indep of the BinomialDen, but coef1 depends on it 
  }
  return(list(dlW_deta=dlW_deta,coef1=coef1)) ## dlW_deta equiv coef2
}

`safesolve.qr.matrix` <- function(qr.a,B,silent=TRUE,stop.on.error=TRUE) { ## solve.qr with fall-back; qr.a should be a qr object, B a matrix
  ## there was a 'Matrix' subcode prior to 10/03/2013
  res <- try(solve.qr(qr.a,B),silent=silent)
  if (inherits(res,"try-error")) { ##FR->FR sb systematique qd phi -> 0; slow step.
    # pivI <- sort.list(qr.a$pivot)  ## inverse perm such as pivI[$pivot]=$pivot[pivI]= identity
    #solveA <- try(solve(qr.R(qr.a)[,pivI])%*%t(qr.Q(qr.a)),silent=silent)  
    ## solveA is solve(<original 'a' matrix>) using the QR decomp... but this may work when solve.qr fails !
    solveA <- try(backsolve(qr.R(qr.a),t(qr.Q(qr.a))[qr.a$pivot,]))
    if (inherits(solveA,"try-error")) {
      if (stop.on.error) {
        mess <- pastefrom("inherits(solveA,'try-error').",prefix="(!) From ")
        message(mess)
        stop("More code is needed to handle this... I exit.") ## perhaps recover A by qr.X and solve(A) ?
      } else return(solveA) ## passes control to calling function
    } else res <- solveA %*% B  
  }
  return(res)
}

`safesolve.qr.vector` <- function(qr.a,b,silent=TRUE,stop.on.error=TRUE) { ## solve.qr with fall-back; qr.a should be a qr object, b must be a vector
  if (class(qr.a)=="sparseQR") { ## pas de 'essai' en var locale !
    ## there was a 'Matrix' subcode prior to 10/03/2013; another try on 11/2013
    res <- qr.coef(qr.a,b)
  } else {
    res <- try(solve.qr(qr.a,b),silent=silent)
    if (inherits(res,"try-error")) {   ## then some weird code, but...
      ## we try to solve(<original 'a' matrix>) using the QR decomp... this may work when solve.qr fails !
      ## The following is equivalent to solveA <- try(solve(qr.R(qr.a)[,pivI])%*%t(qr.Q(qr.a)),silent=silent) then return solveA %*% b 
      ## but uses only backward mutliplication of vectors and transpose of vectors  
      res <- t(crossprod(b, (qr.Q(qr.a)))) ## not yet res
      #pivI <- sort.list(qr.a$pivot) ## inverse perm such as pivI[$pivot]=$pivot[pivI]= identity
      #res <- try(solve(qr.R(qr.a)[, pivI]) %*% res, silent=silent)
      res <- try(backsolve(qr.R(qr.a),res[qr.a$pivot]))
      if (inherits(res,"try-error")) {
        if (stop.on.error) {
          mess <- pastefrom("inherits(res,'try-error').",prefix="(!) From ")
          message(mess)
          stop("More code is needed to handle this... I exit.") ## perhaps recover A by qr.X and solve(A) ?
        } else return(res) ## passes control to calling function
      } 
    }  
  }
  return(res)
}

solveWrap.vector <- function(A,b,...) { ## chol versions not used...
  if (inherits(A,"diagonalMatrix")) return(b/diag(A)) 
  if (inherits(A,"RcppChol")) { ## no pivoting; A$L is .L.ower tri
    return(forwardsolve(A$L,forwardsolve(A$L,b),transpose=TRUE)) ## inv(LLt).b=inv(Lt).(invL.b)
  }
  if (inherits(A,"cholL_LLt")) { ## as results from t(base::chol) with default pivot=FALSE; A is .L.ower tri 
    return(forwardsolve(A,forwardsolve(A,b),transpose=TRUE)) 
  }
  if (inherits(A,"Rcpp_sparseQR")) { ## PIVOTING
    if (length(b)==0L) {return(A$Q_ap)} ## appropriate dimensions
    dim(b) <- c(1,length(b)) ## conversion to "1-rwo matrix" without copy contrary to t(b)
    b <- b %*% (A$Q_ap)
    dim(b) <- c(length(b),1) ## transposition without copy    
    solved <- try(solve(A$R_ap,b)[A$pivI,],silent=TRUE)
    return(solved) ## gives control to calling function 
  } 
  ## next line should become obsolete ?
  if (inherits(A,"sparseQR")) return(Matrix::solve(A,b)) ## as produced by Matrix::qr; return value not documented, but sparse storage is used 
  ## all other cases   
  safesolve.qr.vector(A,b,...)
}

solveWrap.matrix <- function(A,B,...) { ## chol versions not used...
  if (inherits(A,"diagonalMatrix")) return(B/diag(A)) ## works if A is the matrix, not its diagonal...
  if (inherits(A,"RcppChol")) { ## no pivoting; A$L is .L.ower tri
    return(forwardsolve(A$L,forwardsolve(A$L,B),transpose=TRUE)) 
  }
  ## *!* FR->FR confusions ahead
  # class is never "cholL_LLt"; (identical(attr(A,"type"),"cholL_LLt"))  would be a more meaningful test
  # But the code would be wrong anyway, whe A is an LMatrix, not a list with an $L element.
  # => Distinction type de decomp / classe de matrice à garder
  if (inherits(A,"cholL_LLt")) { ## as results from t(base::chol) with default pivot=FALSE; A is .L.ower tri 
    return(forwardsolve(A,forwardsolve(A$L,B),transpose=TRUE)) 
  }
  if (inherits(A,"Rcpp_sparseQR")) { ## PIVOTING
    ## Q and R need not be sparse (even if stored as sparse matrices), can still be sparse in simple aplications
    if (is.identity(B,matrixcheck=FALSE)) {
      solved <- try(solve(A$R_ap,t(A$Q_ap))[A$pivI,],silent=TRUE)
#    } else if (inherits(B,"sparseMatrix")) { 
#      solved <- try(solve(as(A$R_ap,"dtCMatrix"),t(A$Q_ap) %*% B,sparse=TRUE)[A$pivI,],silent=TRUE)
    } else solved <- try(backsolve(A$R_ap,crossprod(A$Q_ap, B))[A$pivI,],silent=TRUE)    
    return(solved) ## gives control to calling function 
  }
  if (inherits(A,"sparseQR")) { ## as produced by Matrix::qr; return value not documented (!), but sparse storage is used
    return(suppressMessages(solve(A,B,sparse=inherits(B,"sparseMatrix")))) 
  }
  ## all other cases
  safesolve.qr.matrix(A,B,...)
}

QRwrap <- function(mat, ## now M or m
                   useEigen=TRUE ## TRUE seems important for large square matrices such as d2hdv2
                   ## FALSE is faster for wAugX for linear solving
                   ) {
  if (inherits(mat,"diagonalMatrix")) { 
    return(mat)
  } else if (inherits(mat,"Matrix") && ncol(mat)<=nrow(mat)) { ## ncol may be > nrow in get_predVar -> ...
    # ... -> .get_logdispObject -> QRwrap(ZAL [cf example cov1 <- get_predVar(fitobject,newdata=moregroups,...) ]
    QR <- Matrix::qr(mat) ##
  } else {
    QR <- qr(mat)
  }
  return(QR)
} 

sym_eigen <- function(X) {
  if (inherits(X,"sparseMatrix")) {
    X <- as.matrix(X) ## dumb, but this is what RSpectra:::eigs_real_sym() does when full eigen is required. 
  }
  if (is.integer(X)) X <- 1.0*X
  return(.selfAdjointSolverCpp(X))
}

.Cholwrap <- function(mat) {
  if (inherits(mat,"diagonalMatrix")) { 
    chol <- diag(ncol(mat))
    class(chol) <- c("RcppChol",class(chol)) ## FR->FR alternatively return a diagonalMatrix ?  
  } else if (inherits(mat,"Matrix")) {
    chol <- t(Matrix::chol(mat))
  } else if (.spaMM.data$options$USEEIGEN) {
    chol <- RcppChol(mat) ##
    if ( chol$Status==1L) { 
      return(chol$L) 
    } else stop("chol$Status !=1L") ## best used in combination with try()
  } else chol <- t(chol(mat)) ## !! this is the R matrix; pivot=FALSE by default
} 

LogAbsDetWrap <- function(mat,logfac=0,provide.qr=FALSE) { ## M or m
  if (ncol(mat)==0) return(0) ## GLM fitted by ML: d2hdbv2 is 0 X 0 matrix 
  # un piege est que mat/(2*pi) conserve les attributes de mat (telle qu'une décomp QR de mat...)
  # il nefaut  dont pas demander LogAbsDetWrap(mat/(2*pi))
  if ( ! is.null(envir <- attr(mat,"envir"))) {
    qrmat <- .get_qr(mat,provide=provide.qr) ##  LogAbsDetWrap' own provide.qr=FALSE
  } else qrmat <- NULL
  if ( ! is.null(qrmat)) { ## 
    if (inherits(qrmat,"qr")) {
      lad <- sum(log(abs(diag(qr.R(qrmat)))))
    } else if (inherits(qrmat,"sparseQR")) { ## if d2hdv2 is a Matrix
      lad <- sum(log(abs(diag(qrmat@R)))) # _@_ ... Matrix::qr.R() serait surement plus recommande 
    } else if (inherits(qrmat,"diagonalMatrix")) { 
      lad <- sum(log(abs(diag(qrmat))))  
    }
  } else if (inherits(mat,"Matrix")) {
    lad <- Matrix::determinant(mat)$modulus[1]
  } else if (.spaMM.data$options$USEEIGEN) {
    lad <- LogAbsDetCpp(mat)
  } else lad <- determinant(mat)$modulus[1]
  # pb general est cert eigenvalues peuvent -> +inf et d'autres -inf auquel cas logabsdet peut être innocuous mais pas estimaable précisément   
  if (is.nan(lad) || is.infinite(lad)){## because of determinant of nearly singular matrix
    zut <- abs(eigen(mat,only.values = TRUE)$values) 
    zut[zut<1e-12] <- 1e-12
    lad <- sum(log(zut)) 
  }
  lad <- lad + nrow(mat)*logfac
  return(lad)
}

singularSigmaMessagesStop <- function(lambda_est,phi_est,corrPars) {
  message("the augmented 'Sigma' matrix appears singular.")
  maxLambda <- max(lambda_est)
  if (min(phi_est)/maxLambda < .Machine$double.eps) {
    cat("This may occur because of small phi/lambda")
    largeLambdaMessages()
    cat(paste("max(lambda estimates)=",maxLambda))
  }
  if (length(corrPars)>0) {
    cat("; correlation parameters=")
    cat(paste(names(corrPars),"=",corrPars))
  }
  stop()
}

.tcrossprod <-  function(x,y=NULL) {
  if (is.null(x)) return(NULL) ## allows lapply(,.tcrossprod) on a listof (matrix or NULL)
  if (inherits(x,"Matrix") || inherits(y,"Matrix")) {
    if (is.null(y)) {
      return(Matrix::tcrossprod(x))
    } else return(Matrix::tcrossprod(x,y))
  } else {
    resu <- tcrossprodCpp(x,y)
    if (is.null(y)) {
      colnames(resu) <- rownames(resu) <- rownames(x)
    } else {
      rownames(resu) <- rownames(x)
      colnames(resu) <- rownames(y)
    }
    return(resu)
  }
}

.crossprod <- function(x,y=NULL) {
  if (is.null(x)) return(NULL) ## allows lapply(,.tcrossprod) on a listof (matrix or NULL)
  if (inherits(x,"Matrix") || inherits(y,"Matrix")) {
    if (is.null(y)) {
      return(Matrix::crossprod(x))
    } else return(suppressWarnings(Matrix::crossprod(x,y))) ## suppressWarnings for case (Matrix,matrix)
  } else {
    resu <- crossprodCpp(x,y)
    if (is.null(y)) {
      colnames(resu) <- rownames(resu) <- colnames(x)
    } else {
      rownames(resu) <- colnames(x)
      colnames(resu) <- colnames(y)
    }
    return(resu)
  }
}

.get_beta_w_cov <- function(res) {
  if (is.null(beta_w_cov <- res$envir$beta_w_cov)) { 
    beta_cov <- .get_beta_cov(res) ## beta_v_cov needed
    beta_w_cov <- attr(beta_cov,"beta_v_cov")
    invL <- .calc_invL(res) ## correlation matrix of ranefs is solve((t(invL)%*%(invL)))
    # invL is currently a single matrix for allranefs. de facto a block matrix when several ranefs
    if ( ! is.null(invL)) {
      pforpv <- ncol(beta_cov)
      v.range <- pforpv+seq(ncol(invL))
      beta_w_cov[v.range,] <- .crossprod(invL, beta_w_cov[v.range,])
      beta_w_cov[,v.range] <- beta_w_cov[,v.range] %*% invL # invL' %*% . %*% invL on the v.range,v.range block
      attr(beta_w_cov,"min_eigen") <- min(eigen(beta_w_cov,only.values = TRUE)$values)
    }
    res$envir$beta_w_cov <- beta_w_cov
  } 
  return(beta_w_cov)
}

get_ZALMatrix <- function(object,as_matrix) {
  # if (object$spaMM.version<"1.11.57") {
  #   stop("This version of get_ZALMatrix() works only on fit objects produced by spaMM v1.11.57 or later.")
  # }
  if (length(ZAlist <- object$ZAlist)>0L) { ## ou tester if (object$models[["eta"]]=="etaGLM")
    if (is.null(object$envir$ZALMatrix)) {
      if (object$spaMM.version<"1.11.57") {
        strucList <- list(dummyid=attr(object$predictor,"LMatrix")) ## back compat
      } else strucList <- object$strucList
      ZALlist <- .compute_ZAXlist(XMatrix=strucList, ZAlist=ZAlist)
      object$envir$ZALMatrix <- .post_process_ZALlist(ZALlist,as_matrix=FALSE) 
    }
    if (as_matrix) {
      object$envir$ZALmatrix <- as.matrix(object$envir$ZALMatrix)
      return(object$envir$ZALmatrix)
    } else return(object$envir$ZALMatrix)
  } else return(NULL) 
}

.get_beta_cov <- function(res) {
  if (is.null(res$envir$beta_cov)) { 
    ## This call gets args (except res) from the envir of locfn def'd below:
    ZAL <- get_ZALMatrix(res,as_matrix=.eval_as_mat_arg(res)) ## should later simplify as =(res$QRmethod=="dense")) FIXME 
    nrd <- length(res$w.ranef)
    pforpv <- ncol(res$X.pv)
    if (inherits(ZAL,"Matrix")) {
      AUGI0_ZX <- list(I=suppressWarnings(as(Diagonal(n=nrd),"CsparseMatrix")),
                       ZeroBlock=Matrix(0,nrow=nrd,ncol=pforpv),X.pv=res$X.pv)
    } else {
      AUGI0_ZX <- list(I=diag(nrow=nrd),ZeroBlock=matrix(0,nrow=nrd,ncol=pforpv),X.pv=res$X.pv)
    }
    res$envir$beta_cov <- .calc_beta_cov(AUGI0_ZX=AUGI0_ZX,ZAL=ZAL,ww=c(res$w.resid,res$w.ranef))
  } 
  return(res$envir$beta_cov)
}

.eval_as_mat_arg <- function(object) { 
  (
    object$HL[1L]=="SEM" || ## SEM code does not yet handle sparse as it uses a dense Sig matrix
    ! identical(object$QRmethod,"sparse") ## => conversion to matrix if object$QRmethod is NULL
  )
}

.get_invColdoldList <- function(res,regul.threshold=1e-7) {
  if (is.null(invColdoldList <- res$envir$invColdoldList)) { 
    ## returns a list of inv(Corr) from the LMatrix
    if (res$spaMM.version<"1.11.57") {
      strucList <- list(dummyid=attr(res$predictor,"LMatrix")) ## back compat
    } else strucList <- res$strucList
    if ( ! is.null(strucList)) {
      cum_n_u_h <- attr(res$lambda,"cum_n_u_h")
      vec_n_u_h <- diff(attr(res$lambda,"cum_n_u_h")) 
      invColdoldList <- lapply(vec_n_u_h,Diagonal)
      ranefs <- attr(res$ZAlist,"ranefs") 
      for (Lit in seq_len(length(strucList))) {
        lmatrix <- strucList[[Lit]]
        if (!is.null(lmatrix)) {
          affecteds <- which(ranefs %in% attr(lmatrix,"ranefs"))
          ## end of designL.from.corr implies either type is cholL_LLt or we have decomp $u and $d 
          type <-  attr(lmatrix,"type")
          invCoo <- NULL
          if (type == "cholL_LLt")  {
            Rmatrix <- t(lmatrix)
          } else { ## Rcpp's symSVD, or R's eigen() => LDL (also possible bad use of R's svd, not normally used)
            condnum <- kappa(lmatrix,norm="1")
            if (condnum<1/regul.threshold) {
              decomp <- attr(lmatrix,attr(lmatrix,"type")) ## of corr matrix !
              invCoo <-  try(ZWZt(decomp$u,1/decomp$d),silent=TRUE)
              if (inherits(invCoo,"try-error")) invCoo <- NULL
            }
            if (is.null(invCoo)) Rmatrix <- qr.R(qr(t(lmatrix))) 
          }
          if (is.null(invCoo)){ ## see comments on chol2inv in .calc_invL()
            singular <- which(abs(diag(Rmatrix))<regul.threshold) 
            if (length(singular)>0L) {
              if (spaMM.getOption("wRegularization")) warning("regularization required.")
              nc <- ncol(Rmatrix)
              diagPos <- seq.int(1L,nc^2,nc+1L)[singular]
              Rmatrix[diagPos] <- sign(Rmatrix[diagPos])* regul.threshold
            }
            invCoo <- chol2inv(Rmatrix) ## 
          }
          for (aff in affecteds) invColdoldList[[aff]] <- invCoo
        } 
      }
      res$envir$invColdoldList <- invColdoldList
    } else return(NULL)
  } 
  return(invColdoldList)
}

.get_info_crits <- function(object) {
  if (is.null(info_crits <- object$envir$info_crits)) { 
    pforpv <- object$dfs[["pforpv"]]
    p_phi <- object$dfs[["p_phi"]]
    p_lambda <- object$dfs[["p_lambda"]]
    APHLs <- object$APHLs
    w.resid <- object$w.resid
    predictor <- object$predictor
    info_crits <- list()
    if  ( ! is.null(resid_fit <- object$resid_fit)) { ## indicates a phiHGLM: we need to get info from it
      # input p_phi (above) is typically set to NA, and will be ignored
      resid_fit <- object$resid_fit
      info_crits_phi <- .get_info_crits(resid_fit)
      phi_pd <- length(resid_fit$y)-info_crits_phi$GoFdf
      p_phi <- phi_pd+sum(resid_fit$dfs) ## all df's absorbed by the phi model
    }
    names_est_ranefPars <- (unlist(.get_methods_disp(object)))  
    p_GLM_family <- length(intersect(names_est_ranefPars,c("NB_shape","NU_COMP")))
    p_phi <- p_phi+p_GLM_family ## effectively part of the model for residual error structure
    # poisson-Gamma and negbin should have similar similar mAIC => NB_shape as one df or lambda as one df   
    forAIC <- APHLs
    if (object$models[[1]]=="etaHGLM") {
      if (object$HL[1]=="SEM") {
        forAIC <- list(p_v=APHLs$logLapp,p_bv=APHLs$logLapp,clik=APHLs$clik)
      } 
      if (object$spaMM.version<"1.11.57") {
        strucList <- list(dummyid=attr(object$predictor,"LMatrix")) ## back compat
      } else strucList <- object$strucList
      ZALlist <- .compute_ZAXlist(XMatrix=strucList, ZAlist=object$ZAlist)
      ZAL <- .post_process_ZALlist(ZALlist,as_matrix=.eval_as_mat_arg(object)) 
      d2hdv2 <- calcD2hDv2(ZAL,w.resid,object$w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
      # if standard ML: there is an REMLformula ~ 0 (or with ranefs ?); processed$X.Re is 0-col matrix
      # if standard REML: REMLformula is NULL: $X.Re is X.pv, processed$X.Re is NULL
      # non standard REML: other REMLformula: $X.Re and processed$X.Re identical, and may take essentially any value
      # if (identical(attr(object$REMLformula,"isML"),TRUE)) {
      #   Md2hdbv2 <- - d2hdv2 
      #   Md2clikdbv2 <-  as.matrix(.ZtWZwrapper(ZAL,w.resid))
      # } else {
      #   ## REML standard || REML non standard
      #   X.Re <- object$distinctX.Re ## null if not distinct from X.pv
      #   if (is.null(X.Re)) X.Re <- object$X.pv ## standard REML
      #   ## diff de d2hdbv2 slmt dans dernier bloc (-> computation pd)
      #   hessnondiag <- .crossprod(ZAL, sweep(X.Re, MARGIN = 1, w.resid, `*`))
      #   Md2hdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
      #                                cbind2(hessnondiag, - d2hdv2))) 
      #   Md2clikdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
      #                                   cbind2(hessnondiag, .ZtWZwrapper(ZAL,w.resid))))            
      # }
      X.pv <- object$X.pv
      if ( ncol(X.pv)>0 ) { ## the projection matrix for the response always includes X even for REML!
        hessnondiag <- .crossprod(ZAL, sweep(X.pv, MARGIN = 1, w.resid, `*`))
        Md2hdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.pv,w.resid), t(hessnondiag)),
                                     cbind2(hessnondiag, - d2hdv2))) 
        Md2clikdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.pv,w.resid), t(hessnondiag)),
                                        cbind2(hessnondiag, .ZtWZwrapper(ZAL,w.resid))))            
      } else {
        Md2hdbv2 <- - d2hdv2 
        Md2clikdbv2 <-  as.matrix(.ZtWZwrapper(ZAL,w.resid))
      }
      p_corrPars <- length(intersect(names_est_ranefPars,names(object$corrPars)))
      info_crits$mAIC <- -2*forAIC$p_v + 2 *(pforpv+p_lambda+p_corrPars+p_phi)
      info_crits$dAIC <- -2*forAIC$p_bv + 2 * (p_lambda+p_phi+p_corrPars) ## HaLM07 (eq 10) focussed for dispersion params
      #                                                                             including the rho param of an AR model
      eigvals <- eigen(Md2hdbv2/(2*pi),only.values = TRUE)$values
      eigvals[eigvals<1e-12] <- 1e-12
      if (min(eigvals)>1e-11) {
        qr.Md2hdbv2 <- QRwrap(Md2hdbv2)
        ## dans un LMM avec estimation ML, pd = sum(lev_phi), mais pas de simplif plus generale 
        pd <- sum(diag(solveWrap.matrix(qr.Md2hdbv2,Md2clikdbv2,stop.on.error=FALSE)))
        if (inherits(pd,"try-error")) {
          warning("Computation of cAIC/GoF df's failed because the 'd2hdbv2' matrix appears singular")
          pd <- NA
        }
      } else pd <- Inf
      info_crits$GoFdf <- length(object$y) - pd ## <- nobs minus # df absorbed in inference of ranefs
      ## eqs 4,7 in HaLM07
      info_crits$cAIC <- -2*forAIC$clik + 2*(pd+p_phi) ## no p_lambda !
      # print(c(pd,p_phi))
      # but Yu and Yau then suggest caicc <- -2*clik + ... where ... involves d2h/db d[disp params] and d2h/d[disp params]2
    } else { ## fixed effect model
      info_crits$mAIC <- -2*forAIC$p_v+2*(pforpv+p_phi) 
    }
    object$envir$info_crits <- info_crits
  } 
  return(info_crits)
}

calcD2hDv2 <- function(ZAL,w.resid,w.ranef) { ## FR->FR review usagesof this function
  ## Si Gamma(identity) avec lambda >1 et w.ranef approche de de -1e6, et si on dit phi <- 1e-06, w.resid = 1e6 
  #    d2hdv2 peut etre une diag matrix with zome 0 elements => logabsdet=log(0)
  if (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U") {
    d2hdv2 <- Diagonal(x= - w.resid - w.ranef)
  } else if (attr(w.resid,"unique")) {
    crossprodZAL <- .crossprod(ZAL)
    d2hdv2 <- - w.resid[1L] * crossprodZAL
    nc <- ncol(d2hdv2)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    d2hdv2[diagPos] <- d2hdv2[diagPos] - w.ranef 
  } else if (inherits(ZAL,"Matrix")) {
    d2hdv2 <- Matrix::crossprod(x=ZAL,y= Diagonal(x= - w.resid) %*% ZAL)    
    nc <- ncol(d2hdv2)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    d2hdv2[diagPos] <- d2hdv2[diagPos] - w.ranef 
  } else d2hdv2 <- Rcpp_d2hdv2(ZAL,w.resid,w.ranef) 
  #cat("\n new d2hdv2")
  structure(d2hdv2,envir=list2env(list(tag="d2hdv2",callcount=0L),parent=environment(HLfit_body)))
}


.get_logdispObject <- function(res) {
  if (is.null(res$envir$logdispObject) && res$models[["eta"]]=="etaHGLM" ) { 
    asDmLR_invV <- .calc_asDmLR_invV_from_fitobject(res)
    dvdloglamMat <- res$envir$dvdloglamMat
    dvdloglamMat_needed <- ( is.null(dvdloglamMat) && 
                               # length(attr(res$ZAlist,"namesTerms"))==1L && ## not only one ranef but not random slope ## but TRUE for CAR
                               # attr(res$ZAlist,"namesTerms")=="(Intercept)" && ## also TRUE for CAR
                               all(unlist(attr(res$ZAlist,"namesTerms"))=="(Intercept)") && ## (1|.) or CAR or Matern
                               any(res$lambda.object$type!="fixed") ) ## some lambda params were estimated
    dvdlogphiMat <- res$envir$dvdlogphiMat
    dvdlogphiMat_needed <- (is.null(dvdlogphiMat) && 
                              res$models[["phi"]]=="phiScal") ## cf comment in calc_logdisp_cov
    if (dvdloglamMat_needed || dvdlogphiMat_needed) {
      ZAL <- get_ZALMatrix(res,as_matrix=.eval_as_mat_arg(res)) ## should later simplify as =(res$QRmethod=="dense")) FIXME       
      d2hdv2 <- calcD2hDv2(ZAL,res$w.resid,res$w.ranef) 
    }
    if (dvdloglamMat_needed) { 
      cum_n_u_h <- attr(res$ranef,"cum_n_u_h")
      psi_M <- rep(attr(res$rand.families,"unique.psi_M"),diff(cum_n_u_h))
      dlogfthdth <- (psi_M - res$ranef)/res$lambda.object$lambda_est ## the d log density of th(u)
      dvdloglamMat <- .calc_dvdloglamMat_new(dlogfthdth=dlogfthdth,
                                             cum_n_u_h=cum_n_u_h,
                                             lcrandfamfam=attr(res$rand.families,"lcrandfamfam"),
                                             rand.families=res$rand.families,
                                             u_h=res$ranef,d2hdv2=d2hdv2,stop.on.error=TRUE)
    }
    if (dvdlogphiMat_needed) {
      muetablob <- res$muetablob
      dh0deta <- ( res$w.resid *(res$y-muetablob$mu)/muetablob$dmudeta ) ## (soit Bin -> phi fixe=1, soit BinomialDen=1)
      dvdlogphiMat  <- .calc_dvdlogphiMat_new(dh0deta=dh0deta, ZAL=ZAL,
                                              d2hdv2=d2hdv2, stop.on.error=TRUE)
    }
    ## This call gets args (except res) from the envir of locfn def'd below:
    res$envir$logdispObject <- calc_logdisp_cov(res, dvdloglamMat=dvdloglamMat, 
                                       dvdlogphiMat=dvdlogphiMat, asDmLR_invV=asDmLR_invV)
  } 
  return(res$envir$logdispObject)
}

.calc_Sig_from_fitobject <- function(object) { ## seems not to be used !?
  predictor <- object$predictor
  if (object$spaMM.version<"1.11.57") {
    strucList <- list(dummyid=attr(object$predictor,"LMatrix")) ## back compat
  } else strucList <- object$strucList
  ZALlist <- .compute_ZAXlist(XMatrix=strucList, ZAlist=object$ZAlist)
  ZAL <- .post_process_ZALlist(ZALlist,as_matrix=.eval_as_mat_arg(object))  
  if (object$models[[1]]=="etaHGLM") {
    if (nrow(ZAL)>3000L) { 
      Sig <- NA
    } else Sig <- Sigwrapper(ZAL,1/object$w.ranef,1/object$w.resid)
  } else Sig <- Diagonal(x=1/object$w.resid)
  return(Sig)
}

.calc_asDmLR_invV_from_fitobject <- function(object) { ## used by .get_logdispObject
  # use representation of V as ZLDL'Z' +D, then use QR repres to control dimensions of matrices
  predictor <- object$predictor
  ZAlist <- object$ZAlist
  ## we need matrices which have the dimensions of Q and R from a QR factorization
  # when there is a single ranef, wecan use ZA and L;
  # otherwise, we need to construct a Q and a R from a global ZAL
  
  if (object$spaMM.version<"1.11.57") {
    strucList <- list(dummyid=attr(object$predictor,"LMatrix")) ## back compat
  } else strucList <- object$strucList
  ZALlist <- .compute_ZAXlist(XMatrix=strucList, ZAlist=ZAlist)
  ZAL <- .post_process_ZALlist(ZALlist,as_matrix=.eval_as_mat_arg(object)) 
  #ZAL <- as.matrix(ZAL)
  qrZAL <- QRwrap(ZAL,useEigen=FALSE)
  if(inherits(qrZAL,"sparseQR")) {
    ZA <- qr.Q(qrZAL)
    Rmat <- qrR(qrZAL)
  } else {
    ZA <- qr.Q(qrZAL)
    Rmat <- qr.R(qrZAL)[,sort.list(qrZAL$pivot)]
  }
  
  invd <- object$w.resid
  ZtinvDZ <- .ZWZtwrapper(t(as.matrix(ZA)), invd) 
  if (is.null(Rmat)) { # no LMatrix
    invRWRt <- Diagonal(x=object$w.ranef)
  } else {
    RWRt <- .ZWZtwrapper(Rmat, 1/object$w.ranef)
    invRWRt <- try(solve(RWRt),silent=TRUE)
    if (inherits(invRWRt,"try-error") || anyNA(invRWRt)) {
      #singularSigmaMessagesStop(lambda_est=lambda,phi_est=object$phi,corrPars=object$corrPars)
      #warning("Generalized inverse used ")
      invRWRt <- MASS::ginv(RWRt) ## FR->FR quick patch at least
    }
  }
  ## central term of Sherman-M-W formula in "(Not) inverting V" section of doc:
  inv2 <- suppressMessages(invRWRt+ZtinvDZ) ## suppress signature message
  invinv <- try(solve(inv2),silent=TRUE)
  if (inherits(invinv,"try-error") || anyNA(invinv)) {
    #singularSigmaMessagesStop(lambda_est=lambda,phi_est=object$phi,corrPars=object$corrPars)
    invinv <- MASS::ginv(inv2) ## FR->FR quick patch at least
  }
  QpinvD <- suppressMessages(sweep(t(ZA), 2L, invd,`*`)) ## suppress message("method with signature...") [found by debug(message)] 
  ## avoid formation of a large nxn matrix:
  return(list(r_x_n=invinv %*% QpinvD, n_x_r=t(QpinvD), invD=invd)) ## invSig = invD- t(QpinvD) %*% invinv %*% QpinvD = invD- n_x_r %*% r_x_n
}

.get_glm_phi <- function(fitobject) {
  if (is.null(fitobject$envir$glm_phi)) { 
    glm_phi_args <- c(fitobject$phi.object$glm_phi_args, 
                      list(formula=attr(fitobject$resid.predictor,"oriFormula"),
                           lev=fitobject$lev_phi, data=fitobject$data,  
                           family= fitobject$resid.family)
    )
    fitobject$envir$glm_phi <-  do.call("calc_dispGammaGLM", glm_phi_args)
  } 
  return(fitobject$envir$glm_phi)
}



calc_wAugX <- function(augX,sqrt.ww) {
  #cat("\n New_AugX ")
  if (inherits(augX,"Matrix")) { 
    wAugX <- augX 
    wAugX@x <- wAugX@x * sqrt.ww[wAugX@i+1L]
  } else wAugX <- .sweepZ1Wwrapper(augX,sqrt.ww) # sweep(TT,MARGIN=1,sqrt.ww,`*`) # rWW %*%TT
  return(  structure(wAugX,envir=list2env(list(tag="wAugX",callcount=0L),parent=environment(HLfit_body))) )
}

.calc_TT <- function(AUGI0_ZX,ZAL) {
  if (inherits(ZAL,"Matrix")) {
    TT <- suppressMessages(cbind2(
      rbind2(AUGI0_ZX$X.pv,AUGI0_ZX$ZeroBlock), 
      rbind2(ZAL, AUGI0_ZX$I)
    )) 
  } else {
    TT <- cbind(
      rbind(AUGI0_ZX$X.pv,AUGI0_ZX$ZeroBlock), 
      rbind(ZAL, AUGI0_ZX$I)
    ) 
  }
  return(TT) ## aug design matrix 
}


.calc_beta_cov <- function(wAugX=NULL, 
                          AUGI0_ZX, 
                          ZAL, ww) {
  if (is.null(wAugX)) {
    if (is.null(ZAL)) {
      wAugX <- calc_wAugX(augX=AUGI0_ZX$X.pv,sqrt.ww=sqrt(ww))
    } else {
      TT <- .calc_TT(AUGI0_ZX=AUGI0_ZX,ZAL)
      wAugX <- calc_wAugX(augX=TT,sqrt.ww=sqrt(ww))
    }
  }
  if (inherits(wAugX,"Matrix")) {
    mMatrix_method <- .spaMM.data$options$Matrix_method
  } else mMatrix_method <- .spaMM.data$options$matrix_method
  # hack to recycle sXaug code; 
  wAugX <- do.call(mMatrix_method,list(Xaug=wAugX, weight_X=rep(1,nrow(AUGI0_ZX$X.pv)), 
                                       w.ranef=rep(1,ncol(AUGI0_ZX$I)), ## we need at least its length for get_from Matrix methods
                                       H_global_scale=1))
  beta_v_cov <- get_from_MME(wAugX,"beta_v_cov_from_wAugX")
  pforpv <- ncol(AUGI0_ZX$X.pv)
  beta_cov <- beta_v_cov[seq_len(pforpv),seq_len(pforpv),drop=FALSE]
  colnames(beta_cov) <- rownames(beta_cov) <- colnames(AUGI0_ZX$X.pv)
  attr(beta_cov,"beta_v_cov") <- beta_v_cov
  return(beta_cov)
}

# not doc'ed (no mention of augmented model in doc)
get_LSmatrix <- function(object,augmented=FALSE) { 
  ## gets inv(tX_a invSig_a X_a).tX_a invSig_a that gives hat(beta,v_h)
  ww <- c(object$w.resid, object$w.ranef)
  sqrt.ww <- sqrt(ww)
  pforpv <- ncol(object$X.pv)
  nrd <- length(object$w.ranef)
  nobs <- nrow(object$X.pv)
  ZAL <- get_ZALMatrix(object)  
  augX <- cbind2(
    rbind2(object$X.pv, matrix(0,nrow=nrd,ncol=pforpv)), 
    rbind2(ZAL, diag(nrow=nrd))
  ) ## template with ZAL block to be filled later
  wAugX <- calc_wAugX(augX=augX,sqrt.ww=sqrt.ww)
  # next line ~ .get_beta_cov(object), but assuming it's not useful to keep the result in memory
  beta_cov <- .calc_beta_cov(wAugX=wAugX,AUGI0_ZX=object$AUGI0_ZX) ## beta_v_cov needed
  beta_v_cov <- attr(beta_cov,"beta_v_cov")
  augXWXXW <- beta_v_cov %*% crossprod(wAugX, diag(x=sqrt.ww))
  if (augmented) {
    return(augXWXXW)
  } else {
    return(augXWXXW[seq_len(pforpv),seq_len(nobs)])
  }
}


## returns a list !!
## input XMatrix is either a single LMatrix whcih is assumed to be the spatial one, or a list of matrices 
.compute_ZAXlist <- function(XMatrix,ZAlist) {
  ## ZAL is nobs * (# levels ranef) and ZA too
  ## XMatrix is (# levels ranef) * (# levels ranef) [! or more generally a list of matrices!]
  ## the levels of the ranef must match each other in multiplied matrices
  ## the only way to check this is to have the levels as rownames and colnames and to check these
  if (is.null(ZAlist)) return(list())
  ## ELSE
  ZAX <- ZAlist
  if ( ! is.null(XMatrix) && length(ZAlist)>0 ) {
    if (inherits(XMatrix,"blockDiag")) {
      stop(".compute_ZAXlist code should be revised to handle blockDiag objects")
    } ## ELSE
    if ( ! inherits(XMatrix,"list")) XMatrix <- list(dummyid=XMatrix)
    LMlen <- length(XMatrix)
    for (ii in seq_len(LMlen)) {
      xmatrix <- XMatrix[[ii]]
      ## find ZAlist elements affected by LMatrix element
      affecteds <- which(attr(ZAlist,"ranefs") %in% attr(xmatrix,"ranefs"))
      for (it in affecteds) {
        ZA <- ZAlist[[it]]
        if (is.identity(ZA)) {
          ZAX[[it]] <- xmatrix          
        } else {
          locnc <- ncol(ZA)
          locnr <- nrow(xmatrix)
          if ( locnc %% locnr !=0) {
            mess <- paste("The number of levels of the grouping variable in random term ", attr(ZAlist,"ranefs")[it],sep="")
            mess <- paste(mess,"\n  is not a multiple of the dimension of the correlation matrix.") ## by distMatrix checking in corrHLfit or no.info check somewhere...
            stop(mess)
          }         
          nblocks <- locnc %/% locnr 
          if (nblocks>1) {
            locZA <- ZA
            for (bt in 1:nblocks) 
              locZA[,locnr*(bt-1)+(1:locnr)] <- locZA[,locnr*(bt-1)+(1:locnr)] %*% xmatrix[] ## [] to handle ff_matrix
            ZAX[[it]] <- locZA
          } else {
            ### it's difficult to make checks on names at this step because the colnames of ZA are not controlled.
            ## ZAlist inherits anything from the .spMMFactorList call which input does not include info about rownames of data
            ## (If this was solved then:
            ##  LMatrix inherits its names from those of uniqueGeo. 
            ##  rownames(xmatrix) are the names of first occurrences of unique geographic locations, 
            ## )
            #             if ( ! all(attr(ZA,"colnames")==rownames(xmatrix))) {
            #               stop("The colnames of the design matrix Z in eta=...+Zv should be the rownames of the design matrix L  in v=Lu")
            #             }
            ###
            ## With a proxy::dist or 'crossdist' matrix, it is likely that ZA was = I and we don't reach this code;
            ## However, exceptions can occur: cf Infusion with CIpoint = MLE => replicate in points where MSEs are to be estimated
            ## Then the xmatrix must have been converted from proxy style to a matrix.
            # In random slope models, xmatrix can be a Matrix
            zax <- ZA %*% xmatrix[]
            attr(zax,"corr.model") <- attr(xmatrix,"corr.model")
            if (identical(attr(zax,"corr.model"),"random-coef")) colnames(zax) <- colnames(ZA)
            #if (named) colnames(ZAX[[it]) <- colnames(ZA) ## names needed for match_old_new_levels()
            ZAX[[it]] <- zax
          }
        }
        attr(ZAX[[it]],"userLfixed") <- attr(xmatrix,"userLfixed") ## TRUE or NULL
      }
    }
  }
  attr(ZAX,"userLfixeds") <- unlist(lapply(ZAX,function(mat) { 
    att <- attr(mat,"userLfixed") ## TRUE or NULL   
    if (is.null(att)) att <- FALSE
    att
  })) ## vector of TRUE or FALSE
  return(ZAX)
}


# even though the Z's were sparse postmultplication by LMatrix leads some of the ZAL's to dgeMatrix (dense)
.post_process_ZALlist <- function(ZALlist, as_matrix ) {
  nrand <- length(ZALlist)
  if ( as_matrix ) {
    ZALlist <- lapply(ZALlist,as.matrix) 
    ZAL <- do.call(cbind,ZALlist)
  } else {
    for (it in seq_len(length(ZALlist))) if (inherits(ZALlist[[it]],"dgeMatrix")) ZALlist[[it]] <- as(ZALlist[[it]],"dgCMatrix")
    ## but leave diagonal matrix types unchanged 
    ZAL <- do.call(cbind,ZALlist)
  } 
  ## desactiv'e 2015/02/15
  #     if (nrand==1L && ( ! inherits(ZAL,"Matrix") ) && 
  #           ## detect a nested random effect: 
  #           (! is.null(attr(predictor,"%in%"))) && attr(predictor,"%in%") && ncol(ZAL)==nrow(ZAL)) {
  #       ## findblocks should become useless...
  #       ## test of the attribute is a heuristic way of detecting when using the block structure will lead to faster analysis
  #       partition <- findblocks(ZAL) 
  #       if ( length(partition)>1 ) {
  #         partition <- cumsum(c(0,partition))
  #         attr(ZAL,"partition") <- partition
  #       }
  #     }
  return(ZAL)
}


## cette fonction marche que si on a fixed effect + un terme aleatoire....
eval.corrEst.args <- function(family,rand.families,predictor,data,X.Re,
                              REMLformula,ranFix,
                              term=NULL,
                              Optimizer) {
  ## ici on veut une procedure iterative sur les params de covariance
  #  HLCor.args$processed <- processed ## FR->FR dangerous in early development
  corrEst.args <- list(family=family,rand.family=rand.families) ## but rand.families must only involve a single spatial effect 
  loc.oriform <- attr(predictor,"oriFormula")
  loc.lhs <- paste(loc.oriform)[[2]]
  ## build formula, by default with only spatial effects
  if (is.null(term)) term <- findSpatial(loc.oriform)
  corrEst.form <-  as.formula(paste(loc.lhs," ~ ",paste(term)))
  corrEst.args$data <- data ## FR->FR way to use preprocess ???                    
  # if standard ML: there is an REMLformula ~ 0; ____processed$X.Re is 0-col matrix, not NULL____
  # if standard REML: REMLformula is NULL: processed$X.Re is NULL
  # non standard REML: other REMLformula: processed$X.Re may take essentially any value
  if (is.null(X.Re) ) { ## processed$X.Re should be NULL => actual X.Re=X.pv => standard REML 
    corrEst.args$REMLformula <- predictor ## standard REML 
  } else corrEst.args$REMLformula <- REMLformula ## _ML_ _or_ non-standard REML
  if (NCOL(X.Re)>0) { ## some REML correction (ie not ML)
    corrEst.args$objective <- "p_bv" ## standard or non-standard REML
  } else corrEst.args$objective <- "p_v" ## ML
  corrEst.args$ranFix <- ranFix ## maybe not very useful
  corrEst.args$control.corrHLfit$Optimizer<- Optimizer ## (may be NULL => L-BFGS-B) 
  corrEst.args$control.corrHLfit$optim$control$maxit <- 1 
  corrEst.args$control.corrHLfit$optimize$tol <- 1e10 
  corrEst.args$control.corrHLfit$maxcorners <- 0 
  return(list(corrEst.args=corrEst.args,corrEst.form=corrEst.form))
}




corr.notEQL.lambda <- function(nrand,cum_n_u_h,lambda_est,lcrandfamfam) {  
  ## d h/ d !log! lambda correction (nul for gaussian ranef)
  ## ! correction for not using the deviance residuals as approx for the distribution of the random effects. It's not specifically ReML !
  ## this is a trick for still using deviances residuals in the Gamma GLM
  notEQL <- unlist(lapply(seq(nrand), function(it) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    loclambda <- lambda_est[u.range]
    blob <- switch(lcrandfamfam[it], 
                   gaussian=rep(0,length(u.range)),
                   gamma=1+2*(log(loclambda)+digamma(1/loclambda))/loclambda,## cf notes on p. 89 of the book
                   "inverse.gamma"=1+2*(log(loclambda)-loclambda+digamma(1+(1/loclambda)) )/loclambda, ## appears to be the same as for the gamma case [digamma(1+x)=digamma(x)+1/x]... 
                   beta=1-2*(digamma(1/loclambda)/loclambda)+2*(digamma(1/(2*loclambda))/loclambda)+log(4)/loclambda
    ) ## as in HGLMMM PoissonEst.r
    blob    
  }))
  return(notEQL)
}

initialize_v_h <- function(psi_M,etaFix,init.HLfit,cum_n_u_h,rand.families,port_env) {
  v_h <- etaFix$v_h
  if (is.null(v_h) ) v_h <- port_env$port_fit_values$v_h
  if (is.null(v_h) ) v_h <- init.HLfit$v_h
  if (is.null(v_h) ) {
    v_h <- unlist(lapply(seq(length(rand.families)), function(it){
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      rand.families[[it]]$linkfun(psi_M[u.range]) ## v as link(mean(u)) 
    }))
  }
  return(v_h)
}

## u_h and v_h in box constraints fro m unconstrainted v_h
.u_h_v_h_from_v_h <- function(v_h,rand.families,cum_n_u_h,lcrandfamfam,lower.v_h,upper.v_h) {
  if(!is.null(lower.v_h)) {v_h[v_h<lower.v_h] <- lower.v_h}
  if(!is.null(upper.v_h)) {v_h[v_h>upper.v_h] <- upper.v_h}
  nrand <- length(rand.families)
  u_list <- lapply(seq(nrand), function(it){
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    uh <- rand.families[[it]]$linkinv(v_h[u.range])
    if (any(is.infinite(uh))) {
      mess <- pastefrom("infinite 'u_h'.",prefix="(!) From ") 
      warning(mess) ## and wait for problems to happen...
    } 
    return(uh) 
  })
  u_h <- unlist(u_list)
  ## if there were box constr, v_h may have been modified, we put it in return value
  if ( ! (is.null(lower.v_h) && is.null(upper.v_h))) attr(u_h,"v_h") <- v_h
  return(u_h)
}

## Aggregate info on corrpars, inner-estimated and inner-fixed.
## $CorrPars is only for info in messages() and return value, 
.get_CorrEst_and_RanFix <- function(ranFix, ## has "fix", "outer", and also "var" values ! code corrently assumes "var" <=> corr_est
                             corr_est ## correlation parameters estimated within HLfit_body
                             ) {
  CorrEst_and_RanFix <- ranFix
  corrNames_in_corr_est <- names(corr_est) # rho,nu,  pas trRho, trNu 
  CorrEst_and_RanFix[corrNames_in_corr_est] <- corr_est
  if (is.null(alltypelist <- attr(ranFix,"type"))) {
    alltypelist <- list()
    alltypelist[names(ranFix)] <- "fix"
  }
  alltypelist[corrNames_in_corr_est] <- "var" ## eg rho from HLCor(adjacency) ## maybe not useful as alredy set to "var"
  attr(CorrEst_and_RanFix,"type") <- alltypelist
  #
  corrNames_in_both <- intersect(names(CorrEst_and_RanFix),c("nu","rho","Nugget","ARphi"))
  corrPars <- CorrEst_and_RanFix[corrNames_in_both] ## as for functions in corrMM.LRT that always look in phi, lambda, rather than .Fix  
  attr(corrPars,"type") <- attr(CorrEst_and_RanFix,"type")[corrNames_in_both]
  return(list(corrPars=corrPars, ## correlation parameters with the appropriate types "fix" or "var"
              CorrEst_and_RanFix=CorrEst_and_RanFix) ## correlation params + whatever else was in ranFix
         ) 
}

.make_lambda_object <- function(nrand, lambda_models, cum_n_u_h, lambda_est, init.lambda,                        
                               process_resglm_blob, coefficients_lambdaS, ZAlist, next_LMatrices) {
  namesTerms <- attr(ZAlist,"namesTerms") ## a list, which names correspond to the grouping variable, and elements are the names of the coefficients fitted
  namesnames <- unlist(lapply(names(namesTerms),function(st) {
    if (nchar(st)>10) st <- paste(substr(st,0,9),".",sep="")
    st
  }))
  names(namesTerms) <- make.unique(namesnames,sep=".") ## makes group identifiers unique (names of coeffs are unchanged); using base::make.unique
  print_lambda <- lapply(seq(nrand), function(it) {
    plam <- process_resglm_blob$print_lambdas[[it]]
    if (anyNA(plam)) { ## those for which no glm was available, such as fixed lambdas...
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      plam <- unique(lambda_est[u.range])
    }
    structure(plam, names=names(namesTerms)[it])
  })
  attr(print_lambda,"cum_n_u_h") <- cum_n_u_h
  lambda.object <- list(lambda_est = lambda_est,  ## full vector for simulate() calc_logdisp_cov()
                        lambda=print_lambda)  ## nrand-elements list, in output -> used by simulate, useful for init another fit, may substitute to the coefficients_lambdaS when the latter have not bee computed from  glm, etc.
  lambda.object$type <- type <- attr(init.lambda,"type")
  if (any(type=="inner")) { ## modifies default namesTerms
    for (it in seq_len(length(coefficients_lambdaS))) { ## detect exceptions
      if (type[it]=="inner") {
        coefficients <- names(coefficients_lambdaS[[it]]) 
        if ("adjd" %in% coefficients) namesTerms[[it]] <- coefficients
      }
    }
    lambda.object <- c(lambda.object,
                       list(coefficients_lambdaS=coefficients_lambdaS,  
                            rand_to_glm_map=process_resglm_blob$rand_to_glm_map,
                            lambda_se=unlist(process_resglm_blob$lambda_seS),
                            linkS = process_resglm_blob$linkS,
                            linkinvS = process_resglm_blob$linkinvS ) 
    )
    attr(lambda.object,"warning") <- unlist(process_resglm_blob$warnmesses) ## may be NULL
  } 
  lambda.object$namesTerms <-  namesTerms
  return(lambda.object)
}

.timerraw <- function(time1) {
  return(round(as.numeric(difftime(Sys.time(), time1, units = "secs")), 1))
}

.make_strucList <- function(ZAlist, next_LMatrices, coefficients_lambdaS,LMatrix) {
  strucList <- replicate(length(ZAlist), list(NULL)) 
  if (is.list(next_LMatrices)) { ## excludes primitive case of a single Matern matrix
    for (it in seq_len(length(ZAlist))) {
      lmatrix <- next_LMatrices[[it]]
      if ( ! is.null(lmatrix)) {
        corr.model <- attr(lmatrix, "corr.model") 
        if (is.null(corr.model)) warning('attr(next_LMatrix, "corr.model") is NULL')
        if (corr.model=="random-coef") {
          longLmatrix <- .makelong(attr(lmatrix,"Lcompact"),longsize=ncol(ZAlist[[it]]))
          ## extra attributes are hidden from str() because lmatrix is an S4 object...
          ## ranefs itself has attribute type="(.|.)"
          ## attr(lmatrix,"Lcompact") doit etre ce que .calc_latentL()$u retrouve
          lmatrix <- do.call(structure,c(list(.Data=longLmatrix),
                                         attributes(lmatrix)[c("Lcompact","par","ranefs","corr.model")]))
          # keep in mind that str(S4...) does not show extra attributes
          attr(lmatrix,"cov.mat") <- ZWZt(attr(lmatrix,"Lcompact"),
                                          exp(coefficients_lambdaS[[which(attr(ZAlist,"Xi_cols")>1L)]])) 
        } 
        strucList[[it]] <- lmatrix  
      } 
    }
  }
  if (!is.null(LMatrix)) {
    ## typically with attributes type="cholL_LLt", corr.model="Matern", ranefs
    ## find ZAlist elements affected by LMatrix element
    affecteds <- which(attr(ZAlist,"ranefs") %in% attr(LMatrix,"ranefs"))
    for (it in affecteds) { strucList[[it]] <- LMatrix }
  }
  return(strucList)
}
