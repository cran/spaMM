
## declared "arglist" to print a clean summary instead of very long list
## the print(summary()) is required if the "arglist" is a member (eg as.list(<call>));
## summary alone would print nothing
print.arglist <-function(x,...) {print(summary(x,...))}
##

calc_TT <- function(X,ZAL) {
  nfix <- ncol(X)
  nrnd <- NCOL(ZAL)
  nobs <- nrow(X)
  if (inherits(ZAL,"Matrix")) {
    OO1 <- Matrix(0L,ncol=nfix,nrow=nrnd)
    X001 <- rbind2(as(X,"sparseMatrix"),OO1) 
    ZALI <- rbind2(ZAL,Diagonal(n=nrnd)) ## this function iscalled only from calc beta cov and attrZALI has been removed
    TT <- cbind2(X001,ZALI)  ## aug design matrix 
  } else {
    OO1 <- matrix(0,ncol=nfix,nrow=nrnd)
    X001 <- rbind2(X,OO1) 
    ZALI <- rbind2(as.matrix(ZAL), diag(nrow=nrnd))
    TT <- cbind2(X001,ZALI)  
  }
  return(TT) ## aug design matrix 
}

sweepZ1Wwrapper <- function(ZZ,WW) { ## for *m*atrix input
  if (nrow(ZZ)!=length(WW)) {
    stop("From sweepZ1Wwrapper(): nrow(ZZ)!=length(WW) ") ## fatal error for eigen code...
  } else if (ncol(ZZ)==0L) {
    return(ZZ)
  } else sweepZ1W(ZZ,WW) ##Rcpp
}

## the following fns try to keep the input class in output, but are called with dense matrices (except irst tested case).
# les Matrix::(t)crossprod  paraissent (parfois au moins) remarquablement inefficaces !!
# idem pour Diagonal()
ZWZtwrapper <- function(ZAL,w) { ## USED ONLY BY Sigwrapper ! => only in calc_Sig_from_fitobject and betaFirst
  if (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U") {
    return(diag(w))
  } else if (inherits(ZAL,"sparseMatrix")) {
    return(as(ZWZt(as.matrix(ZAL),w),"sparseMatrix"))
  } else if (inherits(ZAL,"Matrix")) {
    return(as(ZWZt(as.matrix(ZAL),w),"Matrix"))
  } else return(ZWZt(ZAL,w))
}

ZtWZwrapper <- function(ZAL,w) { ## used in seval contexts
  if (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U") {
    return(diag(w))
  } else if (ncol(ZAL)==0L) {
    return(matrix(NA,ncol=0L,nrow=0L))    
  } else if (inherits(ZAL,"sparseMatrix")) {
    return(as(ZtWZ(as.matrix(ZAL),w),"sparseMatrix"))
  } else if (inherits(ZAL,"Matrix")) {
    return(as(ZtWZ(as.matrix(ZAL),w),"Matrix"))
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
  } else { ## valid for matrix but not for Matrix
    Sig <- Rcpp_Sig(ZAL,wa,wb)
  }
  return(Sig)
}

calcD2hDv2 <- function(ZAL,w.resid,w.ranef) {
  ## Si Gamma(identity) avec lambda >1 et w.ranef approche de de -1e6, et si on dit phi <- 1e-06, w.resid = 1e6 
  #    d2hdv2 peut etre une diag matrix with zome 0 elements => logabsdet=log(0)
  if (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U") {
    d2hdv2 <- Diagonal( - w.resid - w.ranef)
  } else if (attr(w.resid,"unique")) {
    crossprodZAL <- crossprod(ZAL)
    d2hdv2 <- - w.resid[1L] * crossprodZAL
    nc <- ncol(d2hdv2)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    d2hdv2[diagPos] <- d2hdv2[diagPos] - w.ranef 
  } else d2hdv2 <- Rcpp_d2hdv2(ZAL,w.resid,w.ranef) 
  #cat("\n new d2hdv2")
  structure(d2hdv2,get_qr=create_get_qr(tag="d2hdv2"))
}


## uses either X.Re or wAugX
## interesting signed.wAugX concept of version 1.5.3 had not been extended to the leverage computations 
`calc_hatval`  <- function( X.Re,  sqrt.ww, wAugX, levQ=NULL) {
  if ( ! is.null(X.Re) ) { ## not standard REML
    distinct.X.ReML <- attr(X.Re,"distinct.X.ReML")
    if ( distinct.X.ReML[1L] ) wAugX <- wAugX[,-attr(X.Re,"unrestricting_cols")] ## test TRUE for standard ML 
    ## [, -integer(0)] would empty the matrix...
    if ( distinct.X.ReML[2L] ) { ## test FALSE for standard ML, TRUE only for some non-standard REML cases
      extra_vars <- attr(X.Re,"extra_vars")
      suppl_cols <- matrix(0,ncol=length(extra_vars),nrow=nrow(wAugX))
      suppl_cols[1:nrow(X.Re),] <- X.Re[,extra_vars]
      suppl_cols <- sweepZ1Wwrapper(suppl_cols,sqrt.ww)## not obtimized...
      wAugX <- cbind(suppl_cols,wAugX)
    }
    hatval <- leverages(wAugX)
  } else { ## basic REML, leverages from the same matrix used for estimation of betaV (even simply V)
    if (is.null(levQ)) {
      if (FALSE)  {
        ## wAugX updated not only by change in lambda, phi but also GLM weights -> leverage comput difficult to optimize  
        ## the following could be useful if the GLM wights are unity, phiScal et lamScal...
        # mais il manque bouts pour pour produire <u> et <d> | unperturbed RpR = u d u'                
        #                hatval <- LevPerturbedQCpp(perturbedwAugX=wAugX,pforREML=ncol(X.Re),
        #                                           RpRu = <u>,RpRd=<d>,lamOverLam0=lambda/lambda0,phiOverPhi0=phi/phi0)
      } else {
        qr_wAugX <- attr(wAugX,"get_qr")(wAugX,provide=FALSE) ## provide=TRUE tested
        if (is.null(qr_wAugX)) {
          hatval <- leverages(wAugX)
        } else {
          #if (inherits(qr_wAugX,"Rcpp_QR")) {
          #  hatval <- rowSums(qr_wAugX$Q^2)
          #} else 
          if (inherits(qr_wAugX,"sparseQR")) {
            Qy <- Matrix::qr.qy(qr_wAugX, diag(1, nrow = nrow(wAugX), ncol = qr_wAugX$rank))
            hatval <- rowSums(Qy*Qy)
          } else hatval <- stats::hat(qr_wAugX) ## = rowSums(qr.qy(qr_wAugX, diag(1,nrow=nrow(wAugX),ncol=qr_wAugX$rank))^2) 
        }
      }
    } else {hatval <- rowSums(levQ^2)} ## levQ important for intervalStep's local wAugX   
  }
  return(hatval)
}

## direct if no QR repres is available
## FR->FR comparer à leverages() et poss fusionner cette fn et `calc_hatval`
`calc.Pdiag` <- function(ZAL,ww,wAugZALI) {
  nrI <- ncol(ZAL) ## ZALI est rbind(ZAL, I=diag(ncol(ZAL)) donc le nrow(I)=ncol(ZAL)
  if (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U") {
    wZAL <- ww[1L:nrI]    
    wI <- ww[nrI+(1L:nrI)]
    denom <- wZAL+wI
    Pdiag <- c(wZAL/denom,wI/denom) ## sum(ith element, (nrI+i)th element)= 1 for all i
  } else {
    ## next two lines will be slow for very large matrices but the leverages() function using RcppEigen is even slower
  #   partition <- attr(ZAL,"partition")
  #   if ( !is.null(partition) ) { ## use block structure in ZAL;  
  #     Pdiag <- rep(0,nrow(wAugZALI))
  #     abyss <- sapply(seq_len(length(partition)-1),function(v) {
  #       sequ <- (partition[v]+1):partition[v+1] 
  #       colonnes <- wAugZALI[,sequ,drop=F]
  #       qrcolonnes <- qr(colonnes)
  #       levs <- rowSums(qr.qy(qrcolonnes, diag(1, nrow = nrow(colonnes), ncol = ncol(colonnes)))[c(sequ,sequ+nrI),,drop=FALSE]^2)
  #       Pdiag[c(sequ,sequ+nrI)] <<- levs ## local fn proceeds by this side-effect
  #     })
  #   } else { ## straightforward but does not use block structure
     qrwAugZALI <- qr(wAugZALI) 
     Pdiag <- rowSums(qr.qy(qrwAugZALI, diag(1, nrow = nrow(wAugZALI), ncol = ncol(wAugZALI)))^2)
  #   }      
  }
  return(Pdiag)
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
    next_corrEstList$cov12_est <- LMatricesBlob$latest.unique.cov
    for (it in which(isRandomSlope)) {
      u.range <- (cum_n_u_h[it]+1L):cum_n_u_h[it+1L]
      resp_lambda[u.range] <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=1) ## must give d1 in table p 989 de LeeN01
    }
    done[ isRandomSlope ] <- TRUE
  } else { ## only 'declarations' for all further code
    next_lambda_est <- numeric(length(u_h)) ## next_LMatrices remains an empty list()
    next_corrEstList$cov12_est <- NULL
  }
  ### next the other LMatrix models
  for (st in names(LMatricesList) ) { ## this loop will ignore ranefs not affected by any lmatrix
    lmatrix <- LMatricesList[[st]]
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
        next_LMatrices[st] <- LMatricesList[st]
        adjd <- attr(lmatrix,attr(lmatrix,"type"))$adjd
        locdf <- data.frame(adjd=adjd) ## suppose que le type est sym SVD; $d est (1/(1-rho * $adjd)) (ie. corr, not L)
        u.range <- (cum_n_u_h[affected]+1L):cum_n_u_h[affected+1L]
        locdf$resp <- resp_lambda[u.range] <- u_h[u.range]^2
        ## here CAR allows REML conctrary to the SEM CAR, hence leverages
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
      unique.lambda <- pmax(unique.lambda,1e-8) # FR->FR still corrected
      unique.lambda <- pmin(unique.lambda,.spaMM.data$options$maxLambda)  
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
    }
    if (verbose["trace"]) {print(paste("phi_est=",signif(next_phi_est,4)),quote=F)}
  } else { ## random effect(s) in predictor for phi
    stop("This function should not be called when a residual model with random effects is fitted.")
  } 
  return(list(next_phi_est=next_phi_est,  #low phi values are handled in calc.p_v
              glm_phi=glm_phi,
              beta_phi=beta_phi ## used at least to initiate final GLM in "~1" case
  )) ## 
}


`calc.w.resid` <- function(GLMweights,phi_est) { ## One should not correct this phi_est argument by prior.weights (checked)
  phi_est[phi_est<1e-12] <- 1e-11 ## 2014/09/04 local correction, cf comment in calc.p_v
  structure(as.vector(GLMweights/phi_est),unique= (attr(GLMweights,"unique") && length(phi_est)==1L))
}


## spaMM_Gamma() fixes Gamma()$dev.resids(1e10+2,1e10,1) is < 0
# dev.resids() must be >0 for coputation deviance_residual in fitting Gamma GLMM, and alos for $aic() computation.
spaMM_Gamma <- function (link = "inverse") {
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
    pmax(-2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu),.Machine$double.eps)  ## FR: added the pmax
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
    rgamma(nsim * length(ftd), shape = shape, rate = shape/ftd)
  }
  # linkinv <- function (eta) pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax) 
  # : permet des plantages severes dans glm.fit ensuite (R CMD check en detecte) 
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
             comp_z <- COMP_Z(lambda=exp(theta[i]),nu=COMP_nu)
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


`theta.mu.canonical` <- function(mu,family) { 
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
         gaussian = theta.mu.canonical(mu,"gaussian") , ## mu 
         gamma = theta.mu.canonical(mu,"poisson"), ## log(mu)
         beta = theta.mu.canonical(mu,"binomial"), ## improved logit(mu)      
         "inverse.gamma" = theta.mu.canonical(mu,"gamma") ## -1/mu
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

thetaMuDerivs <-function(mu,BinomialDen,family) { ## used for non-canonical links
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
    eta <- pmin(30,eta) ## 100 -> mu = 2.688117e+43 ; 30 -> 1.068647e+13
  } else if (family$family == "COMPoisson" && family$link =="loglambda") {
    eta <- pmin(30*environment(family$aic)$nu,eta) ## using log(mu) ~ eta/nu for large nu
  } else if (family$link=="inverse" && family$family=="Gamma") {
    eta <- pmin(eta,sqrt(.Machine$double.eps)) ## both eta and mu must be >0
  }
  mu <- family$linkinv(eta) ## linkinv(eta) is FREQS for binomial, COUNTS for poisson...
  if (family$link %in% c("logit","probit","cloglog")) {
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

updateWranef <- function(rand.family,lambda,u_h,v_h) {
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

updateW_ranefS <- function(cum_n_u_h,rand.families,lambda,u_h,v_h) {
  nrand <- length(rand.families)
  blob <- lapply(seq(nrand), function(it) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    updateWranef(rand.family=rand.families[[it]],lambda[u.range],u_h[u.range],v_h[u.range])
  })
  w.ranef <- unlist(lapply(blob,function(b) {b$w.ranef}))
  w.ranef <- pmin(w.ranef,1e10) ## patch useful to avoid singular d2hdv2 in PLoG model
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
         cloglog = exp(eta-exp(eta))*(1-exp(eta)) * BinomialDen ## D[1 - E^-E^\[Eta], {\[Eta], 2}]
  )
} 

#derivatives of GLM weights wrt eta 
calc.dlW_deta <- function(dmudeta,family,mu,eta,BinomialDen,canonicalLink,calcCoef1=FALSE) {
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
      blob <- sapply(seq(length(lambdas)), function(i) {
        lambdai <- lambdas[i]
        mui <- mu[i]
        compz <- COMP_Z(lambda=lambdai,nu=COMP_nu)
        compzn <- COMP_Z_n(lambda=lambdai,nu=COMP_nu)
        compzn2 <- COMP_Z_n2(lambda=lambdai,nu=COMP_nu)
        compzn3 <- COMP_Z_n3(lambda=lambdai,nu=COMP_nu)
        rn3 <- COMP_Z_ratio(compzn3,compz)
        rn2 <- COMP_Z_ratio(compzn2,compz)
        dmu.dlogLambda <- rn2 - mui^2 # =family$mu.eta() without repeating some computations
        d2mu.dlogLambda2 <- rn3-mui*rn2-2*mui*dmu.dlogLambda
        #         if (is.nan(dmu.dlogLambda) ) { 
        #           if (COMP_nu<0.05 && lambdai<1) dmu.dlogLambda <- lambdai/(1-lambdai)^2
        #         } ## geometric approx
        #         if (is.nan(dmu.dlogLambda) ) { 
        #           if (COMP_nu<0.05 && lambdai<1) {dmu.dlogLambda <- lambdai*(1+lambdai)/(1-lambdai)^3
        #         } ## geometric approx
        return(c(dmudeta=dmu.dlogLambda, d2mudeta2=d2mu.dlogLambda2))
      })
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
    tmblob <- thetaMuDerivs(mu,BinomialDen,family)
    Dtheta.Dmu <- tmblob$Dtheta.Dmu # calcul co fn de muFREQS puis / BinomialDen
    D2theta.Dmu2 <- tmblob$D2theta.Dmu2 # calcul co fn de muFREQS puis / BinomialDen ^2
    d2mudeta2 <- d2mudeta2fn(link=family$link,mu=mu,eta=eta,BinomialDen=BinomialDen)
    ## ... to compute this:
    D2theta.Deta2_Dtheta.Deta <- (D2theta.Dmu2 * dmudeta^2 + Dtheta.Dmu * d2mudeta2)/(Dtheta.Dmu * dmudeta)
    dlW_deta <- d2mudeta2 / dmudeta + D2theta.Deta2_Dtheta.Deta
    if (calcCoef1) coef1 <- dlW_deta / (Dtheta.Dmu * dmudeta^2) ## note that coef2 is indep of the BinomialDen, but coef1 depends on it 
  }
  return(list(dlW_deta=as.numeric(dlW_deta),coef1=coef1)) ## dlW_deta equiv coef2
}

calc.dvdloglamMat <- function(dlogfthdth,cum_n_u_h,lcrandfamfam,rand.families,u_h,d2hdv2,stop.on.error) {
  neg.d2f_dv_dloglam <- unlist(lapply(seq(length(lcrandfamfam)), function(it) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    ## First the cases where g(u) differs from theta(u) : cf oklink dans preprocess pour detection des cas
    ## same computation as canonical case, except that first we consider dlogfthdv=dlogfthdth * [dth/dv]
    if (lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") { 
      return(dlogfthdth[u.range] / u_h[u.range])  ## [dth/dv=1/u] for th(u)=-1/u, v=log(u)
    } else if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { 
      return(dlogfthdth[u.range] / u_h[u.range]) ## gamma(identity)  ## [dth/dv=1/u] for th(u)=log(u), v=u
    } else { ## v=g(u) = th(u) : random effect model is canonical conjugate
      return(dlogfthdth[u.range]) ## (neg => -) (-)(psi_M-u)/lambda^2    *    lambda.... 
    } 
  }))
  qr_d2hdv2 <- attr(d2hdv2,"get_qr")(d2hdv2,provide=FALSE)
  if (is.null(qr_d2hdv2)) {
    dvdloglamMat <- try(solve(d2hdv2,diag( as.vector(neg.d2f_dv_dloglam) ))) ## FR->FR new 08/2016... ! K2needed but dlW_deta_or_v!=0L
  } else dvdloglamMat <- solveWrap.matrix(qr_d2hdv2, diag( as.vector(neg.d2f_dv_dloglam) ), stop.on.error=stop.on.error) # rXr
  if (inherits(dvdloglamMat,"try-error")) {
    mess <- pastefrom("problem in dvdloglamMat computation.",prefix="(!) From ")
    warning(mess)
    dvdloglamMat <- sweep(ginv(d2hdv2),MARGIN=2,as.vector(neg.d2f_dv_dloglam),`*`) ## ginv(d2hdv2) %*% diag( as.vector(neg.d2f_dv_dloglam))
  }
  return(dvdloglamMat)
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
      res <- t(t(b) %*% (qr.Q(qr.a))) ## not yet res
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
  if (inherits(A,"Rcpp_QR")) { ## no pivoting, $R is upper triangular
    if (length(b)==0L) {return(A$Q)} ## appropriate dimensions wwith zero cols
    dim(b) <- c(1,length(b)) ## conversion to "1-row matrix" without copy contrary to t(b)
    b <- b %*% (A$Q)
    dim(b) <- c(length(b),1) ## transposition without copy    
    solved <- try(backsolve(A$R,b),silent=TRUE)
    return(solved) ## gives control to calling function 
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
  if (inherits(A,"Rcpp_QR")) { ## no pivoting, $R is upper triangular
    if (is.identity(B,matrixcheck=FALSE)) {
      solved <- try(backsolve(A$R,t(A$Q)),silent=TRUE)
#    } else if (inherits(B,"sparseMatrix")) { ## backsolve has no Matrix method and returns a matrix =>
      # use solve on upper triangular sparse Matrix [http://r.789695.n4.nabble.com/backsolve-chol-Matrix-and-SparseM-td4712746.html]
#      solved <- try(solve(as(A$R,"dtCMatrix"),t(A$Q) %*%B,sparse=TRUE),silent=TRUE)
    } else solved <- try(backsolve(A$R,t(A$Q) %*%B),silent=TRUE)
    return(solved) ## gives control to calling function 
  }
  if (inherits(A,"Rcpp_sparseQR")) { ## PIVOTING
    ## Q and R need not be sparse (even if stored as sparse matrices), can still be sparse in simple aplications
    if (is.identity(B,matrixcheck=FALSE)) {
      solved <- try(solve(A$R_ap,t(A$Q_ap))[A$pivI,],silent=TRUE)
#    } else if (inherits(B,"sparseMatrix")) { 
#      solved <- try(solve(as(A$R_ap,"dtCMatrix"),t(A$Q_ap) %*% B,sparse=TRUE)[A$pivI,],silent=TRUE)
    } else solved <- try(backsolve(A$R_ap,t(A$Q_ap) %*% B)[A$pivI,],silent=TRUE)    
    return(solved) ## gives control to calling function 
  }
  ## past version of calc_logdisp_cov -> calc_invS.X(Sig... -> here; suppress signature message;
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
    QR <- list(Q=diag(ncol(mat)),R=mat)
    class(QR) <- c("Rcpp_QR",class(QR)) ## means it must match the Rcpp_QR output  ## FR->FR alternatiely return adiagonalMatrix ?
  } else if (inherits(mat,"Matrix")) {
    QR <- Matrix::qr(mat) ##
    # I have tried  Eigen::SparseMatrix in the past but that failed on win-builder at that time
  #} else if (useEigen) {
  #  QR <- Rcpp_QR(mat) ## the main benefit is that the solve method for Rcpp-QR objects is numerically more stable 
  } else {
    QR <- qr(mat)
  }
  return(QR)
} 

## Chol with fall-back on QR. But this appears slower than QRwrap!
Cholwrap <- function(mat) {
  stop("Cholwrap is buggy") ## FR->FR bc sometimes returns R::chol, sometimes RcppChol's $L=t(R::chol)
  if (inherits(mat,"diagonalMatrix")) { 
    chol <- list(L=diag(ncol(mat)))
    class(chol) <- c("RcppChol",class(chol)) ## FR->FR alternatively return a diagonalMatrix ?  
  } else if (.spaMM.data$options$USEEIGEN) {
    if (inherits(mat,"Matrix")) {
      chol <- Matrix::chol(mat)
    } else {
      chol <- RcppChol(mat) ##
      if ( ! chol$Status==1L) return(QRwrap(mat)) 
    }
  } else {
    chol <- chol(mat) ## !! this is the R matrix; pivot=FALSE by default
  }
  return(chol)
} 

LogAbsDetWrap <- function(mat,logfac=0,provide.qr=FALSE) { ## M or m
  if (ncol(mat)==0) return(0) ## GLM fitted by ML: d2hdbv2 is 0 X 0 matrix 
  # un piege est que mat/(2*pi) conserve les attributes de mat (telle qu'une décomp QR de mat...)
  # il nefaut  dont pas demander LogAbsDetWrap(mat/(2*pi))
  if ( ! is.null(get_qr <- attr(mat,"get_qr"))) {
    qrmat <- get_qr(mat,provide=provide.qr) ##  LogAbsDetWrap' own provide.qr=FALSE
  } else qrmat <- NULL
  if ( ! is.null(qrmat)) { ## 
    if (inherits(qrmat,"qr")) {
      lad <- sum(log(abs(diag(qr.R(qrmat)))))
    } else if (inherits(qrmat,"sparseQR")) { ## if d2hdv2 is a Matrix (if QRmethod == "Matrix::qr") 
      lad <- sum(log(abs(diag(qrmat@R)))) # _@_ ... Matrix::qr.R() serait surement plus recommande 
    } else if (inherits(qrmat,"Rcpp_QR")) { ## true Rcpp output or imitated one
      lad <- sum(log(abs(diag(qrmat$R)))) # _$_ 
    }
  } else if (inherits(mat,"Matrix")) {
    lad <- Matrix::determinant(mat)$modulus[1]
  } else if (.spaMM.data$options$USEEIGEN) {
    lad <- LogAbsDetCpp(mat)
  } else lad <- determinant(mat)$modulus[1]
  # pb general est cert eigenvalues peuvent -> +inf et d'autres -inf auquel cas logabsdet peut être innocuous mais pas estimaable précisément   
  if (is.nan(lad) || is.infinite(lad)){## because of determinant of nearly singular matrix
    zut <- abs(eigen(mat,only.values = T)$values) 
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

# returns t( Sig^{-1}.X )
## this serves mainly for the SEs of beta and disp param (but also in betaFirst)
## One needs Xt.InvS.X, a small matrix but InvS is the slow step.  
calc_tXinvS <- function(Sig,X.pv,stop.on.error) { ## slow... only if betaFirst
  if (ncol(X.pv)==0L) return(matrix(nrow=0,ncol=0))
  invS.X <- calc_invS.X(Sig,X.pv,stop.on.error = stop.on.error) ## invSig %*% X.pv
  if (inherits(invS.X,"try-error")) {
    return(invS.X) ## returns a try-error
  } else tXinvS <- t(invS.X) ## we have to transpose either this one or X.pv, which are of the same size 
  return(tXinvS)
}

`calc_invS.X` <- function(Sig,X.pv,stop.on.error) { ## slow... only used by calc_tXinvS itself only if betaFirst
  if (ncol(X.pv)==0L) return(X.pv)
  if (isDiagonal(Sig)) { ## matrix or Matrix
    if (inherits(Sig,"ddiMatrix") ) {
      invS.X <- solve(Sig,X.pv) ## efficient solve
    } else {
      invS.X <-solve(Diagonal(x=diag(Sig)),X.pv) ## quick patch... but ideally Sig would not be a matrix here
    }
  } else {
    qr.Sig <- QRwrap(Sig) ## Cholwrap tested
    invS.X <- solveWrap.matrix(qr.Sig,X.pv,stop.on.error=stop.on.error) ## invSig %*% X.pv
  }
  return(invS.X) ## may be a try-error. handling code outside as it requires additional info
}


calc_dvdlogphiMat <- function(dh0deta=dh0deta,ZAL=ZAL,d2hdv2=d2hdv2,stop.on.error) {
  ## cf calcul dhdv, but here we want to keep each d/d phi_i distinct hence not sum over observations i 
  ## code corrected here 12/2013; this is dh0dv = neg.d2h0_dv_dlogphi (eta always linear in v :-) and w.resid always propto 1/phi)
  neg.d2h0_dv_dlogphi <- sweep(t(ZAL),MARGIN=2,as.vector(dh0deta),`*`) ## dh0dv <- t(ZAL) %*% diag(as.vector(dh0deta)) ## nXr each ith column is a vector of derivatives wrt v_k
  qr_d2hdv2 <- attr(d2hdv2,"get_qr")(d2hdv2)
  dvdlogphiMat <- solveWrap.matrix(qr_d2hdv2, neg.d2h0_dv_dlogphi , stop.on.error=stop.on.error)  # rXn       
  if (inherits(dvdlogphiMat,"try-error")) {
    mess <- pastefrom("problem in dvdlogphiMat computation.",prefix="(!) From ")
    stop(mess) ## ou warning + ginv  comme dans calc.dvdloglamMat
  }
  return(dvdlogphiMat)
}

tcrossprodWrap <- function(X) {
  if (inherits(X,"Matrix")) {
    return(tcrossprod(X))
  } else return(tcrossprodCpp(X))
}

create_get_w_h_coeffs <- function() {
  ## creates the function in a local envir different from the closure of HLfit
  # so that this envir can be manipulated easily
  w_h_coeffs <- NULL## so that no build/check note on the <<- 
  locfn <- function(res) { 
    ## with res$get_w_h_coeffs <- create_get_w_h_coeffs(), 'res' will be the argument of object$get_w_h_coeffs(object)
    if (is.null(w_h_coeffs)) { 
      ## This call gets args (except res) from the envir of locfn def'd below:
      w_h_coeffs <<- calc_invL_coeffs(res,res$v_h)
    } 
    return(w_h_coeffs)
  } ## this function changes environment(<res object>$get_info_crits)$info_crits
  environment(locfn) <- list2env(list(w_h_coeffs=NULL),
                                 parent=environment(calc_invL_coeffs))
  return(locfn)
}

create_get_beta_w_cov <- function() {
  ## creates the function in a local envir different from the closure of HLfit
  # so that this envir can be manipulated easily
  beta_w_cov <- NULL## so that no build/check note on the <<- 
  locfn <- function(res) {
    if (is.null(beta_w_cov)) { 
      ## This call gets args (except res) from the envir of locfn def'd below:
      beta_w_cov <<- calc_beta_w_cov(res)
    } 
    return(beta_w_cov)
  } ## this function changes environment(<res object>$get_info_crits)$info_crits
  environment(locfn) <- list2env(list(beta_w_cov=NULL),
                                 parent=environment(calc_beta_w_cov))
  return(locfn)
}

create_get_ZALMatrix <- function() {
  ## creates the function in a local envir different from the closure of HLfit
  # so that this envir can be manipulated easily
  ZALMatrix <- NULL## so that no build/check note on the <<- 
  locfn <- function(res) {
    if (is.null(ZALMatrix)) { 
      ## This call gets args (except res) from the envir of locfn def'd below:
      ZALMatrix <<- calc_ZALMatrix(LMatrix=attr(res$predictor,"LMatrix"),ZAlist=res$ZAlist)
    } 
    return(ZALMatrix)
  } ## this function changes environment(<res object>$get_info_crits)$info_crits
  environment(locfn) <- list2env(list(ZALMatrix=NULL),
                                 parent=environment(calc_ZALMatrix))
  return(locfn)
}

calc_ZALMatrix <- function(LMatrix,ZAlist) {
  if (length(ZAlist)>0L) {
    ZALlist <- computeZAXlist(XMatrix=LMatrix,ZAlist=ZAlist)
    ZAL <- post.process.ZALlist(ZALlist,as_matrix=TRUE) ## TRUE is provisional
    return(ZAL)
  } else return(NULL)
}



create_get_beta_cov <- function() {
  ## creates the function in a local envir different from the closure of HLfit
  # so that this envir can be manipulated easily
  beta_cov <- NULL## so that no build/check note on the <<- 
  locfn <- function(res) {
    if (is.null(beta_cov)) { 
      ## This call gets args (except res) from the envir of locfn def'd below:
      ZAL <- res$get_ZALMatrix(res)
      beta_cov <<- calc_beta_cov(X.pv=res$X.pv,ZAL=ZAL,ww=c(res$w.resid,res$w.ranef))
    } 
    return(beta_cov)
  } ## this function changes environment(<res object>$get_info_crits)$info_crits
  environment(locfn) <- list2env(list(beta_cov=NULL),
                                 parent=environment(calc_beta_cov))
  return(locfn)
}


create_get_invColdoldList <- function() {
  ## creates the function in a local envir different from the closure of HLfit
  # so that this envir can be manipulated easily
  invColdoldList <- NULL## so that no build/check note on the <<- 
  locfn <- function(res) {
    if (is.null(invColdoldList)) { 
      ## This call gets args (except res) from the envir of locfn def'd below:
      invColdoldList <<- calc_invColdoldList(res)
    } 
    return(invColdoldList)
  } ## this function changes environment(<res object>$get_info_crits)$info_crits
  environment(locfn) <- list2env(list(invColdoldList=NULL),
                                 parent=environment(calc_invColdoldList))
  return(locfn)
}

create_get_info_crits <- function(pforpv, p_lambda, p_phi) {
  ## creates the function in a local envir different from the closure of HLfit
  # so that this envir can be manipulated easily
  info_crits <- NULL## so that no build/check note on the <<- 
  locfn <- function(res) {
    if (is.null(info_crits)) { 
      ## This call gets args (except res) from the envir of locfn def'd below:
      info_crits <<- calc_info_crits(res, pforpv=pforpv, p_lambda=p_lambda, p_phi=p_phi)
    } 
    return(info_crits)
  } ## this function changes environment(<res object>$get_info_crits)$info_crits
  environment(locfn) <- list2env(list(info_crits=NULL,
                                      pforpv=pforpv, p_lambda=p_lambda, p_phi=p_phi),
                                 parent=environment(calc_info_crits))
  return(locfn)
}

create_get_logdispObject <- function(dvdloglamMat, ## may be NULL, in which case it may be recreated 
                                     dvdlogphiMat, ## may be NULL, in which case it may be recreated
                                     muetablob,
                                     stop.on.error) {
  ## creates the function in a local envir different from the closure of HLfit
  # so that this envir can be manipulated easily
  logdispObject <- NULL## so that no build/check note on the <<- 
  locfn <- function(res) {
    if (is.null(logdispObject) && res$models[["eta"]]=="etaHGLM" ) { 
      asDmLR_invV <- calc_asDmLR_invV_from_fitobject(res)
      dvdloglamMat_needed <- ( is.null(dvdloglamMat) && 
        length(attr(res$ZAlist,"namesTerms"))==1L && ## not only one ranef but not random slope ## but TRUE for CAR
        attr(res$ZAlist,"namesTerms")=="(Intercept)" && ## also TRUE for CAR
        length(res$lambda.object$coefficients_lambda)>0L ) ## some lambda params were estimated
      dvdlogphiMat_needed <- (is.null(dvdlogphiMat) && 
                                res$models[["phi"]]=="phiScal") ## cf comment in calc_logdisp_cov
      if (dvdloglamMat_needed || dvdlogphiMat_needed) {
        ZAL <- res$get_ZALMatrix(res)
        ddi_or_matrix_ZAL <- post_process_ZAL(ZAL,attr(res$w.resid,"unique"))
        d2hdv2 <- calcD2hDv2(ddi_or_matrix_ZAL,res$w.resid,res$w.ranef) 
      }
      if (dvdloglamMat_needed) { 
        cum_n_u_h <- attr(res$ranef,"cum_n_u_h")
        psi_M <- rep(attr(res$rand.families,"unique.psi_M"),attr(cum_n_u_h,"vec_n_u_h"))
        dlogfthdth <- (psi_M - res$ranef)/res$lambda.object$lambda_est ## the d log density of th(u)
        dvdloglamMat <- calc.dvdloglamMat(dlogfthdth=dlogfthdth,
                                          cum_n_u_h=cum_n_u_h,
                                          lcrandfamfam=attr(res$rand.families,"lcrandfamfam"),
                                          rand.families=res$rand.families,
                                          u_h=res$ranef,d2hdv2=d2hdv2,stop.on.error=stop.on.error)
      }
      if (dvdlogphiMat_needed) {
        dh0deta <- ( res$w.resid *(res$y-muetablob$mu)/muetablob$dmudeta ) ## (soit Bin -> phi fixe=1, soit BinomialDen=1)
        dvdlogphiMat  <- calc_dvdlogphiMat(dh0deta=dh0deta, ZAL=ZAL,
                                           d2hdv2=d2hdv2, stop.on.error=stop.on.error)
      }
      ## This call gets args (except res) from the envir of locfn def'd below:
      logdispObject <<- calc_logdisp_cov(res, dvdloglamMat=dvdloglamMat, 
                                         dvdlogphiMat=dvdlogphiMat, asDmLR_invV=asDmLR_invV,
                                         stop.on.error=stop.on.error)
    } 
    return(logdispObject)
  } ## Inner locfn() changes members of environment(outer fn <res object>$get_logdispObject)
  ## INITIALIZATION of envir of outer fn function created by create_get_logdispObject(): 
  environment(locfn) <- list2env(list(logdispObject=NULL, dvdloglamMat=dvdloglamMat, 
                                      dvdlogphiMat=dvdlogphiMat, 
                                      muetablob=muetablob, stop.on.error=stop.on.error),
                                 parent=environment(calc_logdisp_cov))
  return(locfn)
}

calc_Sig_from_fitobject <- function(object) {
  predictor <- object$predictor
  ZALlist <- computeZAXlist(XMatrix=attr(predictor,"LMatrix"),
                            ZAlist=object$ZAlist)
  ZAL <- post.process.ZALlist(ZALlist,as_matrix=TRUE) ## TRUE is provisional 
  if (object$models[[1]]=="etaHGLM") {
    if (nrow(ZAL)>3000L) { 
      Sig <- NA
    } else Sig <- Sigwrapper(ZAL,1/object$w.ranef,1/object$w.resid)
  } else Sig <- Diagonal(x=1/object$w.resid)
  return(Sig)
}

calc_asDmLR_invV_from_fitobject <- function(object) {
  predictor <- object$predictor
  ZAlist <- object$ZAlist
  if (length(ZAlist)==1L) {
    ZA <- ZAlist[[1]]
    Rmat <- attr(predictor,"LMatrix") 
  } else {
    ZALlist <- computeZAXlist(XMatrix=attr(predictor,"LMatrix"),
                              ZAlist=ZAlist)
    ZAL <- post.process.ZALlist(ZALlist,as_matrix=TRUE) ## TRUE is provisional 
    qrZAL <- QRwrap(ZAL,useEigen=FALSE)
    ZA <- qr.Q(qrZAL)
    Rmat <- qr.R(qrZAL)
  }
  invd <- object$w.resid
  ZtinvDZ <- ZWZtwrapper(t(as.matrix(ZA)), invd) 
  if (is.null(Rmat)) { # no LMatrix
    invRWRt <- Diagonal(x=object$w.ranef)
  } else {
    RWRt <- ZWZtwrapper(Rmat, 1/object$w.ranef)
    invRWRt <- try(solve(RWRt),silent=TRUE)
    if (inherits(invRWRt,"try-error") || anyNA(invRWRt)) {
      #singularSigmaMessagesStop(lambda_est=lambda,phi_est=object$phi,corrPars=object$corrPars)
      #warning("Generalized inverse used ")
      invRWRt <- ginv(RWRt) ## FR->FR quick patch at least
    }
  }
  inv2 <- suppressMessages(invRWRt+ZtinvDZ) ## suppress signature message
  invinv <- try(solve(inv2),silent=TRUE)
  if (inherits(invinv,"try-error") || anyNA(invinv)) {
    #singularSigmaMessagesStop(lambda_est=lambda,phi_est=object$phi,corrPars=object$corrPars)
    invinv <- ginv(inv2) ## FR->FR quick patch at least
  }
  QpinvD <- suppressMessages(sweep(t(ZA), 2L, invd,`*`)) ## suppress message("method with signature...") [found by debug(message)] 
  ## avoid formation of a large nxn matrix:
  return(list(r_x_n=invinv %*% QpinvD, n_x_r=t(QpinvD), invD=invd)) ## invSig = invD- t(QpinvD) %*% invinv %*% QpinvD = invD- n_x_r %*% r_x_n
}

create_get_glm_phi <- function(fitobject) {
  ## creates the function in a local envir different from the closure of HLfit
  # so that this envir can be manipulated easily
  glm_phi <- NULL## so that no build/check note on the <<- 
  locfn <- function(fitobject) {
    if (is.null(glm_phi)) { 
      glm_phi_args <- c(fitobject$phi.object$glm_phi_args, 
                        list(formula=attr(fitobject$resid.predictor,"oriFormula"),
                             lev=fitobject$lev_phi, data=fitobject$data,  
                             family= fitobject$resid.family)
                        )
      glm_phi <<-  do.call("calc_dispGammaGLM", glm_phi_args)
    } 
    return(glm_phi)
  } ## this function changes environment(<res object>$get_glm_phi)$glm_phi
  environment(locfn) <- list2env(list(glm_phi=NULL),
                                 parent=environment(calc_info_crits))
  return(locfn)
}


calc_info_crits <- function(object, pforpv, p_lambda, p_phi) {
  APHLs <- object$APHLs
  w.resid <- object$w.resid
  predictor <- object$predictor
  info_crits <- list()
  if  ( ! is.null(resid_fit <- object$resid_fit)) { ## indicates a phiHGLM: we need to get info from it
    # input p_phi (above) is typically set to NA, and will be ignored
    resid_fit <- object$resid_fit
    info_crits_phi <- resid_fit$get_info_crits(resid_fit)
    phi_pd <- length(resid_fit$y)-info_crits_phi$GoFdf
    phi_other_pS <- unlist(as.list(environment(resid_fit$get_info_crits))[c("pforpv","p_lambda","p_phi")])
    p_phi <- phi_pd+sum(phi_other_pS) ## all df's absorbed by the phi model
  }
  if (object$models[[1]]=="etaHGLM") {
    if (object$HL[1]=="SEM") {
      forAIC <- list(p_v=APHLs$logLapp,p_bv=APHLs$logLapp)
    } else forAIC <- APHLs
    ZALlist <- computeZAXlist(XMatrix=attr(predictor,"LMatrix"),ZAlist=object$ZAlist)
    ZAL <- post.process.ZALlist(ZALlist,as_matrix=TRUE) ##TRUE is provisional
    ddi_or_matrix_ZAL <- post_process_ZAL(ZAL,attr(w.resid,"unique"))
    d2hdv2 <- calcD2hDv2(ddi_or_matrix_ZAL,w.resid,object$w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    X.Re <- object$distinctX.Re
    if (is.null(X.Re)) X.Re <- object$X.pv
    if ( ncol(X.Re)>0 ) {  ## diff de d2hdbv2 slmt dans dernier bloc (-> computation pd)
      ## REML standard || REML non standard
      hessnondiag <- suppressWarnings(crossprod(ZAL, sweep(X.Re, MARGIN = 1, w.resid, `*`))) ## Matrix or matrix depending on ZAL; suppressWarnings for case crossprod(Matrix,matrix)
      Md2hdbv2 <- as.matrix(rbind2(cbind2(ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
                                  cbind2(hessnondiag, - d2hdv2))) 
      Md2clikdbv2 <- as.matrix(rbind2(cbind2(ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
                                     cbind2(hessnondiag, ZtWZwrapper(ZAL,w.resid))))            
    } else {
      Md2hdbv2 <- - d2hdv2 
      Md2clikdbv2 <-  as.matrix(ZtWZwrapper(ZAL,w.resid))
    }
    ## see readable account of aic in HaLM07
    info_crits$mAIC <- -2*forAIC$p_v + 2 *(pforpv+p_lambda+p_phi)
    info_crits$dAIC <- -2*forAIC$p_bv + 2 * (p_lambda+p_phi) ## HaLM (10) focussed for dispersion params
    if (object$HL[1]=="SEM") {
      message("Effective df and cAIC not (yet) computed for SEM.")
      # Md2hdbv2and $clik missing
    } else {
      # pd for GoFdf and cAIC
      eigvals <- eigen(Md2hdbv2/(2*pi),only.values = TRUE)$values
      eigvals <- pmax(eigvals,1e-12)
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
      info_crits$cAIC <- -2*APHLs$clik + 2*(pd+p_lambda) ## FR-FR note that HLfit does not determine which correlation parameters are estimated, but they enter through pd?
      # but Yu and Yau then suggest caicc <- -2*clik + ... where ... involves d2h/db d[disp params] and d2h/d[disp params]2
    }
  } else info_crits$mAIC <- -2*APHLs$p_v+2*pforpv
  return(info_crits)
}


calc_logdisp_cov <- function(object, dvdloglamMat=NULL, dvdlogphiMat=NULL, asDmLR_invV=NULL, stop.on.error) { 
  
  # evaluate tr(A %*% B)= sum(A*B) where A and B are large matrices but each of the form l %*% r for 'narrow' l
  # this function avoids the formation of the large matrices, using a form of commutation of trace arguments.
  traceAB <- function(lA,rA,lB,rB) {
    ll <- t(lA) %*% lB
    rr <- rA %*% t(rB)
    return(sum(ll*rr))
  }
  # and same concept for trace( D %*% B)
  traceDB <- function(dD,lB,rB) { sum( sweep(lB * t(rB),1L,dD,`*`) )}
  # lA <- matrix(runif(6),ncol=2)
  # rA <- matrix(runif(6),ncol=3)
  # lB <- matrix(runif(6),ncol=2)
  # rB <- matrix(runif(6),ncol=3)
  # traceAB(lA,rA,lB,rB)
  # dD <- runif(3)
  # sum(diag(diag(dD) %*% lB %*% rB))
  # traceDB(dD,lB,rB)

  lambda.object <- object$lambda.object
  LMatrix <- attr(object$predictor,"LMatrix")
  adjacency <- identical(attr(LMatrix,"corr.model"),"adjacency")
  ar1 <- identical(attr(LMatrix,"corr.model"),"ar1")
  ## debut de code de detection parmi +s ranefs
  #whichadj <- attr(attr(ZAlist,"ranefs"),"type")=="adjacency"  
  #notscalar <- unlist(lapply(coefficients_lambdaS,length)==2L)
  dwdlogphi <- dwdloglam <- NULL ## always a cbind at the end of calc_logdisp_cov
  dispnames <- c()
  problems <- list() ## Its elements' names are tested in calcPredVar, and the strings are 'development info'
  if (length(lambda.object$coefficients_lambda)>0L) { ## not lambda.Fix
    ## in the first cases, a dvdloglamMat may exist (eg salamander, male and female effects)
    if (adjacency || ar1) {
      ## cf summary.HLfit for extracting info from object
      locit <- 1L
      it <- 1L
      if ("adjd" %in% object$lambda.object$namesTerms[[it]]) {
        rho <- - with(lambda.object,coefficients_lambda[locit+1L]/coefficients_lambda[locit])
        lambda <- with(lambda.object,linkinvS[[rand_to_glm_map[it]]](coefficients_lambda[locit]))
      } else {
        lambda <- with(lambda.object,linkinvS[[rand_to_glm_map[it]]](coefficients_lambda[locit]))        
        rho <- object$corrPars$rho
      }
      adjd <- attr(LMatrix,"symSVD")$adjd
      denom <- 1-rho*adjd
      problems$warnadj <- warning("lambda dispVar component not implemented for adjacency model.") 
      # il me manque dwdrho (et meme dwdloglam pour ce modele ?) donc on inactive les lignes suivantes:
      #       if (is.null(lambda.object$lambda.fix)) dispnames <- c(dispnames,"loglambda")
      #       corrFixNames <- names(unlist(object$corrPars[which(attr(corrPars,"type")=="fix")]))
      #       if (! ("rho" %in% corrFixNames) ) dispnames <- c(dispnames,"rho")
    } else if (length(lambda.object$coefficients_lambda)>1L) { 
      problems$warnmult <- warning("lambda dispVar component not implemented for model with several lambda parameters.")
    } else if (is.null(dvdloglamMat)) {
      problems$stopmiss <- warning("is.null(dvdloglamMat) in a case where it should be available.") 
    }else {
      dvdloglam <- rowSums(dvdloglamMat) ## assuming each lambda_i = lambda
      dwdloglam <- calc_invL_coeffs(object,dvdloglam)
      if ( all(lambda.object$type=="inner")) dispnames <- c(dispnames,"loglambda")
    }
  } else lambda <- NULL # no ranefs
  
  phimodel <- object$models[["phi"]]
  if (phimodel=="phiScal") { ## semble impliquer pas outer phi.Fix... => no need to test object$phi.object$phi_outer,"type")
    phi_est <- object$phi ## no need to get object$phi.object$phi_outer
    if (length(phi_est)!=1L) problems$stopphi <- warning("phimodel=\"phiScal\" but length(phi_est)!=1L.")
    if ( ! is.null(dvdlogphiMat)) {
      dvdlogphi <- rowSums(dvdlogphiMat) ## using each phi_i = phi
      dwdlogphi <- calc_invL_coeffs(object,dvdlogphi)
      dispnames <- c(dispnames,"logphi")
    } else if (object$models[["eta"]]=="etaHGLM") stop("phimodel=='phiScal' but is.null(dvdlogphiMat)")
  } else {  ## else phimodel="", e.g. binomial
    # if binomial or poisson, phimodel=""; warning for other phimodels
    if (phimodel!="") problems$structphi <- warning("phi dispVar component not yet available for phi model != ~1.")
  }
  ## compute info matrix:
  if ((nrc <- length(dispnames))==0L) {
    return(list(problems=problems))
  } else {
    # cf my documentation, based on McCullochSN08 6.62 and 6.74
    # lambda and phi factors enter in dV/dlog(.), computed instead of dV/d(.) to match dwdlog(.) vectors.
    ZAL <- object$get_ZALMatrix(object) # attr(object$predictor,"ZALMatrix")
    logdispInfo <- matrix(NA,nrow=nrc,ncol=nrc)
    colnames(logdispInfo) <- rownames(logdispInfo) <- dispnames
    if ("loglambda" %in% dispnames) {
      if (adjacency || ar1) { 
        #  lambda already evaluated 
        ZALd <- ZAL %id*id% Diagonal(x=sqrt(1/denom))
        lhs_invV.dVdlam <- asDmLR_invV$n_x_r %*% (asDmLR_invV$r_x_n %*% ZALd)  
        lhs_invV.dVdlam <- sweep( ZALd,1L,asDmLR_invV$invD,`*`) - lhs_invV.dVdlam
        rhs_invV.dVdlam <- t(ZALd)
      } else { ## standard lamScal model
        lambda <- exp(lambda.object$coefficients_lambda)
        lhs_invV.dVdlam <- asDmLR_invV$n_x_r %*% (asDmLR_invV$r_x_n %*% ZAL)  
        lhs_invV.dVdlam <- sweep( ZAL,1L,asDmLR_invV$invD,`*`) - lhs_invV.dVdlam
        rhs_invV.dVdlam <- t(ZAL)
      }
      # lambda^2 *sum(invV.dVdlam* t(invV.dVdlam)) :
      logdispInfo["loglambda","loglambda"] <- lambda^2 *traceAB(lhs_invV.dVdlam,rhs_invV.dVdlam,
                                                                t(rhs_invV.dVdlam),t(lhs_invV.dVdlam)) 
    }
    if ("rho" %in% dispnames) {
      # no use of sqrt because adjd can be negative
      #invV.dVdrho <- (invV %id*id% ZAL) %*% ( Diagonal(x=lambda*adjd/(denom^2)) %id*id% t(ZAL))
      lhs_invV.dVdrho <- asDmLR_invV$n_x_r %*% (asDmLR_invV$r_x_n %*% ZAL)  
      lhs_invV.dVdrho <- sweep( ZAL,1L,asDmLR_invV$invD,`*`) - lhs_invV.dVdrho
      rhs_invV.dVdrho <- ( Diagonal(x=lambda*adjd/(denom^2)) %id*id% t(ZAL))
      #logdispInfo["rho","rho"] <- sum(invV.dVdrho*t(invV.dVdrho))
      logdispInfo["rho","rho"] <- traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho))
      if ("loglambda" %in% dispnames) {
        logdispInfo["loglambda","rho"] <- 
          logdispInfo["rho","loglambda"] <- lambda* ( 
            traceDB(asDmLR_invV,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho)) -
            traceAB(lhs_invV.dVdlam,rhs_invV.dVdlam,t(rhs_invV.dVdrho),t(lhs_invV.dVdrho)) 
          ) #lambda * sum(invV.dVdlam*t(invV.dVdrho))
      }
      } 
    ## if (! is.null(dwdlogphi)) { ## currently length(phi)==1L && ! is.null(dvdlogphiMat)
    if ("logphi" %in% dispnames) { ## more transparent, but error if mismatch of conditions
      ## next lines assume that  the design matrix for the residual error is I
      logdispInfo["logphi","logphi"] <- phi_est^2 * (
        sum(asDmLR_invV$invD^2) -2*traceDB(asDmLR_invV$invD,asDmLR_invV$n_x_r, asDmLR_invV$r_x_n) + 
          traceAB(asDmLR_invV$n_x_r, asDmLR_invV$r_x_n, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)
        ) # phi_est^2 * sum(invV^2)
      if ("loglambda" %in% dispnames) {
        logdispInfo["loglambda","logphi"] <- 
          logdispInfo["logphi","loglambda"] <- lambda * phi_est * (
            traceDB(asDmLR_invV$invD, lhs_invV.dVdlam, rhs_invV.dVdlam) -
            traceAB(lhs_invV.dVdlam, rhs_invV.dVdlam, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)
          ) # lambda * phi_est * sum(invV.dVdlam * invV)
      }
      if ("rho" %in% dispnames) {
        logdispInfo["rho","logphi"] <- 
          logdispInfo["logphi","rho"] <- phi_est * traceAB(lhs_invV.dVdrho,rhs_invV.dVdrho, asDmLR_invV$n_x_r, asDmLR_invV$r_x_n)  
        # phi_est * sum(invV.dVdrho * invV)  
      }
    } 
    logdispInfo <- logdispInfo/2
    logdisp_cov <- try(solve(logdispInfo),silent=TRUE)
    if (inherits(logdisp_cov,"try-error")) logdisp_cov <- ginv(logdispInfo) ## quick patch for uninteresting case
    dwdlogdisp <- cbind(dwdloglam,dwdlogphi) ## typically nobs * 2
    return(list(dwdlogdisp=dwdlogdisp,logdisp_cov=logdisp_cov,problems=problems)) ## more compact than storing ww %*% logdisp_cov %*% t(ww) which is nobs*nobs 
  }
}

calc_wAugX <- function(augX,sqrt.ww) {
  #cat("\n New_AugX ")
  if (inherits(augX,"Matrix")) { 
    wAugX <- sweepZ1Wwrapper(as.matrix(augX),sqrt.ww) 
    # wAugX <- Diagonal(x=sqrt.ww) %*% augX ## apparemment sweep vriament pas efficace sur Matrix
  } else wAugX <- sweepZ1Wwrapper(augX,sqrt.ww) # sweep(TT,MARGIN=1,sqrt.ww,`*`) # rWW %*%TT
  return(  structure(wAugX,get_qr=create_get_qr(tag="wAugX")) )
}


calc_beta_cov <- function(wAugX=NULL, 
                          X.pv, ## for pforpv if for nothing else (rethink ?)
                          ZAL, ww) {
  if (is.null(wAugX)) {
    if (is.null(ZAL)) {
      wAugX <- calc_wAugX(augX=X.pv,sqrt.ww=sqrt(ww))
    } else {
      TT <- calc_TT(X.pv,ZAL)
      wAugX <- calc_wAugX(augX=TT,sqrt.ww=sqrt(ww))
    }
  }
  qr_wAugX <- attr(wAugX,"get_qr")(wAugX)
  if (inherits(qr_wAugX,"sparseQR")) {
    qrR <- suppressWarnings(as.matrix(Matrix::qr.R(qr_wAugX))) ## suppress qrR warning. We need a triangular matrix
    beta_v_cov <- chol2inv(qrR)
    perm <- qr_wAugX@q + 1L
    beta_v_cov[perm, perm] <- beta_v_cov
  } else if (! is.null(qrR <- qr_wAugX$R)) { ## non-uniform output of QRwrap bc of useeigen
    beta_v_cov <- chol2inv(qrR) 
  } else beta_v_cov <- chol2inv(qr.R(qr_wAugX))
  pforpv <- ncol(X.pv)
  beta_cov <- beta_v_cov[seq_len(pforpv),seq_len(pforpv),drop=FALSE]
  colnames(beta_cov) <- rownames(beta_cov) <- colnames(X.pv)
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
  ZAL <- calc_ZALMatrix(LMatrix=attr(object$predictor,"LMatrix"),ZAlist=object$ZAlist)  
  augX <- cbind2(
    rbind2(object$X.pv, matrix(0,nrow=nrd,ncol=pforpv)), 
    rbind2(ZAL, diag(nrow=nrd))
  ) ## template with ZAL block to be filled later
  wAugX <- calc_wAugX(augX=augX,sqrt.ww=sqrt.ww)
  beta_cov <- calc_beta_cov(wAugX=wAugX,X.pv=object$X.pv) ## X.pv for its ncol only...
  beta_v_cov <- attr(beta_cov,"beta_v_cov")
  augXWXXW <- beta_v_cov %*% t(wAugX) %*% diag(x=sqrt.ww)
  if (augmented) {
    return(augXWXXW)
  } else {
    return(augXWXXW[seq_len(pforpv),seq_len(nobs)])
  }
}

## for SEM

rntneg <- function(n,mu,sigma2)
{
  # produce n samples from the
  # specified rigth-truncated to 0 gaussian
  # distribution
  pn <- runif(n)*pnorm(0,mu,sqrt(sigma2))
  pn[pn==0] <-  .Machine$double.eps ## because if mu is large -> qnorm(0) is -Inf which later cause NaN's
  pn[pn==1] <-  1-.Machine$double.eps 
  qnorm(pn,mu,sqrt(sigma2))
  ## alternatively use ... qn[pn==0] <- ... 
}

################################################################################

rntpos <- function(n,mu,sigma2)
{
  # produce n samples from the
  # specified left-truncated to 0 gaussian
  # distribution
  -rntneg(n,-mu,sigma2)
}

## returns a list !!
## input XMatrix is either a single LMatrix whcih is assumed to be the spatial one, or a list of matrices 
computeZAXlist <- function(XMatrix,ZAlist) {
  ## ZAL is nobs * (# levels ranef) and ZA too
  ## XMatrix is (# levels ranef) * (# levels ranef) [! or more generally a list of matrices!]
  ## the levels of the ranef must match each other in multiplied matrices
  ## the only way to check this is to have the levels as rownames and colnames and to check these
  if (is.null(ZAlist)) return(list())
  ## ELSE
  namesTerms <- attr(ZAlist,"namesTerms")
  ZAX <- ZAlist
  if ( ! is.null(XMatrix) && length(ZAlist)>0 ) {
    if (inherits(XMatrix,"blockDiag")) {
      stop("computeZAXlist code should be revised to handle blockDiag objects")
    } ## ELSE
    if ( ! inherits(XMatrix,"list")) XMatrix <- list(dummyid=XMatrix)
    LMlen <- length(XMatrix)
    for (ii in seq_len(LMlen)) {
      lmatrix <- XMatrix[[ii]]
      ## find ZAlist elements affected by LMatrix element
      affecteds <- which(attr(ZAlist,"ranefs") %in% attr(lmatrix,"ranefs"))
      for (it in affecteds) {
        ZA <- ZAlist[[it]]
        if (is.identity(ZA)) {
          ZAX[[it]] <- lmatrix          
        } else {
          locnc <- ncol(ZA)
          locnr <- nrow(lmatrix)
          if ( locnc %% locnr !=0) {
            mess <- paste("The number of levels of the grouping variable in random term (...|",names(namesTerms)[it],")",sep="")
            mess <- paste(mess,"\n  is not the dimension of the correlation matrix.") ## by distMatrix checking in corrHLfit or no.info check somewhere...
            stop(paste(mess," I exit."))
          }         
          nblocks <- locnc %/% locnr 
          if (nblocks>1) {
            locZA <- ZA
            for (bt in 1:nblocks) 
              locZA[,locnr*(bt-1)+(1:locnr)] <- locZA[,locnr*(bt-1)+(1:locnr)] %*% lmatrix[] ## [] to handle ff_matrix
            ZAX[[it]] <- locZA
          } else {
            ### it's difficult to make checks on names at this step because the colnames of ZA are not controlled.
            ## ZAlist inherits anything from the spMMFactorList call which input does not include info about rownames of data
            ## (If this was solved then:
            ##  LMatrix inherits its names from those of uniqueGeo. 
            ##  rownames(lmatrix) are the names of first occurrences of unique geographic locations, 
            ## )
            #             if ( ! all(attr(ZA,"colnames")==rownames(lmatrix))) {
            #               stop("The colnames of the design matrix Z in eta=...+Zv should be the rownames of the design matrix L  in v=Lu")
            #             }
            ###
            ## With a proxy::dist or 'crossdist' matrix, it is likely that ZA was = I and we don't reach this code;
            ## However, exceptions can occur: cf Infusion with CIpoint = MLE => replicate in points where MSEs are to be estimated
            ## Then the lmatrix must have been converted from proxy style to a matrix
            ZAX[[it]] <- ZA %*% lmatrix[]
            # rajout 02/2016: but in random slope models, lmatrix can be a Matrix
          }
        }
        attr(ZAX[[it]],"userLfixed") <- attr(lmatrix,"userLfixed") ## TRUE or NULL
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
post.process.ZALlist <- function(ZALlist, as_matrix ) {
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

post_process_ZAL <- function(ZAL,attr_wresid_unique) {
  if (inherits(ZAL,"Matrix") && ! (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U")) {
    ddi_or_matrix_ZAL <- as.matrix(ZAL)
  } else ddi_or_matrix_ZAL <- ZAL
  return(ddi_or_matrix_ZAL)
}

intervalStep <- function(old_betaV,wAugX,wAugz,currentlik,intervalInfo,
                         ranFix=ranFix,corr_est=corr_est, ## both onnly for message()
                         likfn,QRmethod) {
  parmcol <- attr(intervalInfo$parm,"col")
  #print((control.HLfit$intervalInfo$fitlik-currentlik)/(control.HLfit$intervalInfo$MLparm-old_betaV[parmcol]))
  ## voir code avant 18/10/2014 pour une implem rustique de VenzonM pour debugage  
  ## somewhat more robust algo (FR->FR: still improvable ?), updates according to a quadratic form of lik near max
  ## then target.dX = (current.dX)*sqrt(target.dY/current.dY) where dX,dY are relative to the ML x and y 
  ## A nice thing of this conception is that if the target lik cannot be approached, 
  ##   the inferred x converges to the ML x => this x won't be recognized as a CI bound (important for locoptim) 
  currentDx <- (old_betaV[parmcol]-intervalInfo$MLparm)
  targetDy <- (intervalInfo$fitlik-intervalInfo$targetlik)
  currentDy <- (intervalInfo$fitlik-currentlik)
  betaV <- rep(NA,length(old_betaV))
  if (currentDy <0) { 
    locmess <- paste("A higher",likfn,"was found than for the original fit.",
                     "\nThis suggests the original fit did not fully maximize",likfn,"\n or numerical accuracy issues.")
    message(locmess)
    dispCorrPars <- get_DispCorrPars(ranFix,corr_est)
    if (length(dispCorrPars)>0) message(paste("Current dispersion and correlation parameters are ",
                                          paste(names(dispCorrPars),"=",signif(unlist(dispCorrPars),6),collapse=", ")))
    message("Current likelihood is =",currentlik)                    
    message("lik of the fit=",intervalInfo$fitlik)    
    betaV[parmcol] <- old_betaV[parmcol]
  } else {
    betaV <- rep(NA,length(old_betaV))
    Dx <- currentDx*sqrt(targetDy/currentDy)
    ## pb is if Dx=0 , Dx'=0... and Dx=0 can occur while p_v is still far from the target, because other params have not converged.
    ## patch:
    if (currentDy<targetDy) { ## we are close to the ML: we extrapolate a bit more confidently
      min_abs_Dx <- intervalInfo$asympto_abs_Dparm/1000
    } else min_abs_Dx <- 1e-6 ## far from ML: more cautious move our of Dx=0
    Dx <- sign(currentDx)*max(abs(Dx),min_abs_Dx)
    betaV[parmcol] <- intervalInfo$MLparm + Dx 
  }
  locwAugX <- wAugX[,-parmcol,drop=FALSE]
  locwAugz <- as.matrix(wAugz-wAugX[,parmcol]*betaV[parmcol])
  if (inherits(wAugX,"Matrix")) {
    if (QRmethod == "Matrix::qr") {
      locqr_wAugX <- Matrix::qr(locwAugX)
      betaV[-parmcol] <- as.matrix(Matrix::qr.coef(locqr_wAugX,locwAugz))  
      betaVQ <- list(Q=as.matrix(Matrix::qr.Q(locqr_wAugX))) ## un peu nouille, mais unifie l'interface
    } else {stop("Unknown (processed) QRmethod")}
  } else { ## wAugX is matrix not Matrix (lmwithQ_denseZAL), with useEigen
    locqr_wAugX <- QRwrap(locwAugX,useEigen=FALSE) ## already matrix
    betaV[-parmcol] <- solveWrap.vector(locqr_wAugX,locwAugz)
    betaVQ <- list(Q=qr.Q(locqr_wAugX)) ## nouille aussi
  }
  return(list(levQ=betaVQ$Q,betaV=betaV)) # levQ ispresumably always dense
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
  corrEst.args$control.corrHLfit$nlminb$control$iter.max <- 2 ## 1 => convergence trop lente
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

initialize.v_h <- function(psi_M,etaFix,init.HLfit,cum_n_u_h,rand.families) {
  v_h <- etaFix$v_h
  if (is.null(v_h) ) v_h <- init.HLfit$v_h
  if (is.null(v_h) ) {
    v_h <- unlist(lapply(seq(length(rand.families)), function(it){
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      rand.families[[it]]$linkfun(psi_M[u.range]) ## v as link(mean(u)) 
    }))
  }
}

## u_h and v_h in box constraints fro m unconstrainted v_h
u_h_v_h_from_v_h <- function(v_h,rand.families,cum_n_u_h,lcrandfamfam,lower.v_h,upper.v_h) {
  if(!is.null(lower.v_h)) {v_h <-pmax(v_h,lower.v_h)}
  if(!is.null(upper.v_h)) {v_h <-pmin(v_h,upper.v_h)}
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


selfAdjointWrapper <- function(mat) {
  if (is.integer(mat)) mat <- 1.0*mat
  selfAdjointSolverCpp(mat)
}

pMVN <- function (L, limits, ismax, rand_seq=NULL,nrep=1000) {
  if (length(ismax) != length(limits)) {
    stop("lengths of arguments 'limits' and 'ismax' are not identical in pMVN()")
  } else dim <- length(ismax)
  if (nrow(L)!=dim) stop("nrow(L) != dim in pMVN()")
  if (is.null(rand_seq)) {
    rand_seq <- matrix(nrow=0L,ncol=0L)
  } else if (ncol(rand_seq)!=dim) stop("ncol(rand_seq) != nrow(L) in pMVN()")
  blob <- Rcpp_pMVN(L, limits, ismax, rand_seq=rand_seq, nrep=nrep)
  return(blob) ## FR->FR NaN seInt's are produced... not necess harmful (try arasample with n=40 locations)
}

provide_init_lambdas <- function(processed,stillNAs,
                               init_lambda_by_glm) {
  nrand <-  length(processed$ZAlist)
  preproFix <- processed$lambda.Fix
  family <- processed$family
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  resu <- rep(NA,nrand)
  
  ### (1) Generate a single value for all ranefs ## FR->FR this could  be simplified to generate missing values only
  init.lambda <- init_lambda_by_glm
  ## les tests et examples sont mieux sans la $variance, mais ce n'est pas attendu theor pour large N
  ## with variance: cf CoullA00 p. 78 top for motivation, et verif par simul:
  ## the excess rel var (given by their corr_rr' for r=r') is of the order of (their rho_rr=)lambda
  ## hence excess var (overdisp) is lambda Npq hence lambda ~ overdisp/Npq ~ overdisp/(resglm$prior.weights*family$variance(fv)) ?
  ## pas convainquants en temps:
  #init.lambda <- max(0.01,sum(((resid(resglm,type="pearson")/resglm$prior.weights)^2)/resglm$family$variance(fv)-1)/resglm$df.residual )
  #init.lambda <- max(0.01,sum((resid(resglm,type="pearson")/resglm$prior.weights)^2-resglm$family$variance(fv))/resglm$df.residual )
  ## pas loin du temps de réference ?:
  #init.lambda <- sum(resid(resglm,type="pearson")^2/(resglm$prior.weights*family$variance(fv)))/resglm$df.residual
  if (family$link=="log") { ## test of family, not rand.family...
    init.lambda <- log(1.00001+init.lambda) ## max(0.0001,log(init.lambda))
  } else init.lambda <- init.lambda/5 ## assume that most of the variance is residual
  init.lambda <- init.lambda/nrand        
  ###
  #
  ### (2) generate missing values from this unique value, each adapted to a family and link
  # allows for different rand.family
  valuesforNAs <- unlist(lapply(stillNAs, function(it) {
    if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity" && init.lambda==1) {
      adhoc <- 0.9999
    } else if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="log") {
      objfn <- function(lambda) {psigamma(1/lambda,1)-init.lambda}
      adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
    } else if(lcrandfamfam[it]=="beta" && rand.families[[it]]$link=="logit") {
      #ad hoc approximation which should be quite sufficient; otherwise hypergeometric fns.
      objfn <- function(lambda) {8* lambda^2+3.2898*lambda/(1+lambda)-init.lambda}
      adhoc <- uniroot(objfn,interval=c(2.5e-6,1e8))$root
    } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") {
      ## this tries to controle var(v)<-init.lambda of v=log(u) by defining the right lambda for var(u)=lambda/(1-lambda)
      ## log (X~inverse.G) = - log (~Gamma), ie
      ##  log (X~inverse.G(shape=1+1/init.lambda,scale=1/init.lambda))~ - log (rgamma(shape=1+1/init.lambda,scale=1/init.lambda)
      ## where + log (~Gamma) has known mean (=...) and variance=trigamma(shape)
      ## (pi^2)/6=1.644934... is upper bound for trigamma(x->1) ie for lambda ->infty.
      ## however, there is something wrong with setting very large initial lambda values. Hence
      # if (init.lambda > 1.64491 ) { ## trigamma(1+1/100000) 
      #   adhoc <- 100000 
      if (init.lambda > 0.3949341 ) { ## trigamma(3) [for lambda=0.5]
        adhoc <- 0.5 
      } else { ## loer init.lambda -> lower lambda
        objfn <- function(lambda) {psigamma(1+1/lambda,1)-init.lambda}
        adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
      }
      ## but the mean of v  is function of lambda and I should rather try to control the second moment of v
    } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="-1/mu") {
      adhoc <- (sqrt(1+4*init.lambda)-1)/2 # simple exact solution
    } else adhoc <- init.lambda
    adhoc
  }))
  resu[stillNAs] <- valuesforNAs
  return(resu)
} 


## dispCorrPars is only for info in messages() and return value, 
get_DispCorrPars <- function(ranFix, corr_est) {
  corrNames_in_ranFix <- intersect(names(ranFix),c("nu","rho","Nugget","ARphi"))
  dispCorrPars <- ranFix[corrNames_in_ranFix] ## as for functions in corrMM.LRT that always look in phi, lambda, rather than .Fix. 
  corrNames_in_init_HLfit <- names(corr_est) # rho,nu,  pas trRho, trNu 
  dispCorrPars[corrNames_in_init_HLfit] <- NA ## will be filled at the end of the fit
  typelist <- list()
  typelist[corrNames_in_ranFix] <- "fix"
  if (!is.null(rFtype <- attr(ranFix,"type"))) { 
    corrNames_in_ranFix_type <- intersect(corrNames_in_ranFix,names(rFtype))
    typelist[corrNames_in_ranFix_type] <- rFtype[corrNames_in_ranFix_type]
  }
  typelist[corrNames_in_init_HLfit] <- "var" 
  attr(dispCorrPars,"type") <- typelist
  dispCorrPars[corrNames_in_init_HLfit] <- corr_est #[corrNames_in_init_HLfit]
  return(dispCorrPars)
}

make_lambda_object <- function(nrand, lambda_models, cum_n_u_h, lambda_est, 
                               init.lambda, ## for attr(.,"type")
                               coefficients_lambdaS, 
                               process_resglm_blob, ZAlist, next_LMatrices) {
  namesTerms <- attr(ZAlist,"namesTerms") ## a list, which names correspond to the grouping variable, and elements are the names of the coefficients fitted
  namesnames <- unlist(lapply(names(namesTerms),function(st) {
    if (nchar(st)>10) st <- paste(substr(st,0,9),".",sep="")
    st
  }))
  names(namesTerms) <- make.unique(namesnames,sep=".") ## makes group identifiers unique (names of coeffs are unchanged); using base::make.unique
  print_lambda <- lapply(seq(nrand), function(it) {
    if (lambda_models[it]=="lamScal") {
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      structure(unique(lambda_est[u.range]),names=names(namesTerms)[it]) ## using list names, not element names   
    } else {
      print_lambda <- lambda_est ## pseudocode,not clear what would be best
      message("not lamScal: no namescreated for res$lambda")
    }
  })
  if (all(lambda_models=="lamScal")) print_lambda <- unlist(print_lambda) ## format for easy display... 
  ## but also in output -> used by simulate, useful for init another fit etc.
  attr(print_lambda,"cum_n_u_h") <- cum_n_u_h
  lambda.object <- list(lambda_est = lambda_est,  ## full vector for simulate() calc_logdisp_cov()
                        lambda=print_lambda)  ## nrand-vector, duplicate info
  lambda.object$type <- type <- attr(init.lambda,"type")
  if (any(type=="inner")) { ## modifies default namesTerms
    for (it in seq_len(length(coefficients_lambdaS))) { ## detect exceptions
      if (type[it]=="inner") {
        coefficients <- names(coefficients_lambdaS[[it]]) 
        if ("adjd" %in% coefficients) namesTerms[[it]] <- coefficients
      }
    }
    lambda.object <- c(lambda.object,
                       list(coefficients_lambda=unlist(coefficients_lambdaS),  
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

