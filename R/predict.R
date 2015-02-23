`%id*id%` <- function(A,B) {
  if (is.identity(A)) return(B)
  if (is.identity(B)) return(A)
  return(A %*% B)
}

#t.identityMatrix <- function(x) {x} ## (t is a generic fn)  ## added benefit of keeping attributes... 

calcPredVar <- function(phi,lambda,Coldnew,Cnewnew,oldX.pv,X.pv,ZA,newZA,beta_cov,Sig) {
  if (is.list(newZA)) {
    nrand <- length(newZA) ## or of any other of the lists of matrices
    sumLambdaZCZtfn <- function(lambda_short,leftZA,rightZA,Clist) {
      locfn <- function(it) {
        ## still assuming a single lambda for each u.range (otherwise needs true cov mat instead of Clist[[it]]):
        terme <- lambda_short[it] * leftZA[[it]] %id*id% Clist[[it]] %id*id% t(rightZA[[it]]) ## %id*id% further keeps col names of C if tZA==I
        as.matrix(terme) ## loses names but they are not useful here 
      }
      varMatList <- lapply(seq_len(nrand), locfn)
      Reduce("+", varMatList) ## matrix sum, https://stat.ethz.ch/pipermail/r-help/2008-June/163817.html 
    } ## keeps first explict names -> no clean pattern at this point
    lamZCZt_nn <- sumLambdaZCZtfn(lambda_short=lambda,
                                  leftZA=newZA,rightZA=newZA,
                                  Clist=Cnewnew
    )
    lamZCZt_on <- sumLambdaZCZtfn(lambda_short=lambda,
                                  leftZA=ZA,rightZA=newZA,
                                  Clist=Coldnew
    )  
    invSig.lamZCZt_on <- solve(Sig,lamZCZt_on) ## obs*obs.obs*new = obs*new
  } else {
    tnewZA <- t(newZA) ## to keep the attribute...
    # Cnn -c'.inv(Coo).c + (x-X'.inv(V).c)'.inv(X'inv(V)X). (x-X'.inv(V).c)
    lamZCZt_on <- lambda * ZA %id*id% Coldnew %id*id% tnewZA ## obs*new
    invSig.lamZCZt_on <- solve(Sig,lamZCZt_on) ## obs*obs.obs*new = obs*new
    if (is.identity(newZA)) {
      lamZCZt_nn <- lambda * Cnewnew
    } else lamZCZt_nn <- lambda * newZA %*% Cnewnew %*% tnewZA
  }
  ## for known beta:
  ## several ranefs -> different lambda, not handled
  # term for known beta:
  predVar <- lamZCZt_nn- t(lamZCZt_on) %*% invSig.lamZCZt_on ## new*new- new*obs.obs*new = new*new
  ## correction for estimated beta:
  if (! is.null(beta_cov)) { ## FR->FR test post 1.4.1 CRAN
    K <- t(X.pv)-t(oldX.pv) %*% invSig.lamZCZt_on ## p*new -(p*obs).(obs*new) = p*new
    predVar <- predVar + t(K) %*% beta_cov %*% K 
  } ## for lambda=0 predVar reduces to X.pv %*% beta_cov %*% t(X.pv) hence it's still the variance of etz 
  return(as.matrix(predVar)) ## need to make sure it's a matrix not a Matrix, for later usage
}

calcResidVar <- function(object,newdata=NULL) {
  phi.object <- object$phi.object
  if ( ! is.null(phi.object$phi.Fix)) {
    if (length(phi.object$phi.Fix)==1L) {
      residVar <- rep(phi.object$phi.Fix,nrow(newdata))
    } else if (is.null(newdata)) {
      residVar <- rep(phi.object$phi.Fix)          
    } else {
      stop("Unable to compute 'residVar' given 'newdata' argument and given vector of phi values in fit object.")
    }
  } else {
    ## creates pseudo HLfit object
    phipred <- object$resid.predictor
    if (!inherits(phipred,"predictor")) phipred <- Predictor(phipred) ## but the offset is not always set...
    pseudoObject <- list(fixef=object$phi.object$beta_phi, predictor=phipred,family=object$resid.family)
    residVar <- predict.HLfit(pseudoObject, newdata=newdata,
                              residVar=FALSE) ## avoid recursive call        
    attributes(residVar) <- NULL
  }
  residVar
}  

calcNewCorrs <- function(object,locdata,predVar,spatial.model) {
  olduniqueGeo <- attr(object,"info.uniqueGeo")
  geonames <- colnames(olduniqueGeo)
  newgeo <- locdata[,geonames,drop=FALSE] ## not unique
  newuniqueGeo <- unique(newgeo) ## It is essential that it keeps the same order as spMMfactorlist -> ULI -> unique. 
  ### rho only used to compute scaled distances
  rho <- getPar(object$ranFix,"rho")
  if( !is.null(rho.mapping <- object$control.dist$rho.mapping)) rho <- rho[rho.mapping] # has been missing
  ## rows from newuniqueGeo, cols from olduniqueGeo:
  msd.arglist <- list(uniqueGeo=newuniqueGeo,uniqueGeo2=olduniqueGeo,rho=rho)
  ### rho not further used 
  if( ! is.null(dist.method <- object$control.dist$dist.method)) msd.arglist$dist.method <- dist.method
  uuCnewold <- do.call(make.scaled.dist,msd.arglist)
  uuCnewold <- matrix(uuCnewold,ncol=ncol(uuCnewold),nrow=nrow(uuCnewold)) ## ultimately allows products with Matrix ## '*cross*dist' has few methods, not even as.matrix
  if (predVar)  {
    msd.arglist$uniqueGeo2 <- NULL
    Cnewnew <- do.call(make.scaled.dist,msd.arglist) ## dist matrix and then call to correl fn
    Cnewnew <- as.matrix(Cnewnew) 
  }
  if (spatial.model[[1]]=="AR1") {
    args <-object$ranFix[which(names(object$ranFix) %in% c("ARphi"))]
    uuCnewold <- args$ARphi^uuCnewold  
    if (predVar) Cnewnew <- args$ARphi^Cnewnew  
  } else {
    args <-object$ranFix[which(names(object$ranFix) %in% c("nu","Nugget"))] ## so that rho=1 in Matern.corr
    uuCnewold <- do.call(Matern.corr,args=c(args,list(d=uuCnewold)))  
    if (predVar) Cnewnew <- do.call(Matern.corr,args=c(args,list(d=Cnewnew)))  
  }
  resu <- list(uuCnewold=uuCnewold)
  if (predVar) resu <- c(resu,list(Cnewnew=Cnewnew))
  return(resu)
}

makenewname <- function(base,varnames) { ## FR->FR post CRAN 1.4.1
  varnames <- varnames[which(substring(varnames,1,nchar(base))==base)] 
  allremainders <- substring(varnames,nchar(base)+1) 
  allremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  if (length(allremainders) == 0) {
    fvname <- base
  } else {
    num <- max(allremainders)+1
    fvname <-paste ( base , num , sep="") 
  }
  fvname
}

## (1) for surface prediction: (developed in InferentialSimulation/InferentialSimulation.R)
## (2) But also for generation of fixed effects in simulation of nested-effects models
predict.HLfit <- function(object,newdata=newX,newX = NULL,coeffs=NULL,re.form = NULL,
                          variances=list(fixef=FALSE,linPred=FALSE,resid=FALSE,sum=FALSE,cov=FALSE),
                          predVar=variances$linPred,residVar=variances$resid,
                          binding = FALSE, 
                          ...) { ## but not new Y
  if ( ! is.null(variances$ranef)) {
    message("'variances$ranef' is obsolete: I interpret this as 'variances$linPred'.")
    variances$linPred <- variances$ranef
    variances$ranef <- NULL
  }
  if (is.null(variances$sum)) variances$sum <- FALSE ## transient value; so that the following lines do not leave the lhs as NULL
  if (is.null(variances$linPred)) variances$linPred <- variances$sum
  if (is.null(variances$resid)) variances$resid <- variances$sum
  if (is.null(variances$fixef)) variances$fixef <- FALSE ## 12/2014 post 1.4.4
  if (is.null(variances$cov)) variances$cov <- FALSE
  if ( ! is.null(residVar) && residVar) variances$resid <- TRUE
  if (predVar == "Cov") {
    variances$linPred <- TRUE; variances$cov <- TRUE
  } else if ( ! is.null(predVar) && predVar) variances$linPred <- TRUE
  variances$sum <- (variances$linPred && variances$resid)
  ##
  locform <- attr(object$predictor,"oriFormula")
  #if (length(locform)==2L) locform <- as.formula(paste('"LHS"',paste(locform,collapse=" "))) ## ## syntax no longer requires LHS
  spatial.terms <- findSpatial(locform) 
  spatial.model <- spatial.terms[[1]] 
  if (length(locform)==3L) locform <- locform[-2] ## R CMD CHECK if  this is removed...
  allvars <- all.vars(locform) 
  if (is.vector(newdata)) { ## ## less well controlled case, but useful for maximization
    binding <- binding ## :-) binding must be evaluated before newdata is modified
    newdata <- data.frame(matrix(newdata,nrow=1))
    if (length(allvars)==ncol(newdata)) {
      names(newdata) <- allvars
    } else {
      stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",paste(allvars,collapse=" ")))
    }
  } 
  ## it is important that newdata remains NULL if it was so initially because it is tested below. Hence copy in locdata
  if (is.null(newdata)) {
    locdata <- object$data[,allvars,drop=FALSE]
  } else {
    if( is.matrix(newdata) ) newdata <- as.data.frame(newdata)  
    # so that matrix 'newdata' arguments can be used as in some other predict methods.
    locdata <- newdata ## FR->FR [,allvars,drop=FALSE] ?
    if (any(is.na(locdata))) {
      stop("NA's in required variables from 'newdata'. Prediction not possible.")
    }
  }
  #   if (is.null(newdata)) {
  #     locdata <- object$data
  #   } else {
  #     locdata <- newdata
  #   }
  ## possible change of random effect terms
  if (noReForm(re.form)) locform <- nobarsMM(locform)
  ## preparation for fixed effects
  allFrames <- HLframes(formula=locform,data=locdata)
  if ( ( ! is.null(newdata) ) || variances$linPred || noReForm(re.form)) {
    etaFix <- newetaFix(object,allFrames)  ## newfixed effects      
  }  
  if (noReForm(re.form)) {
    nrand <- 0
  } else {
    nrand <- length(attr(object$ZAlist,"ranefs"))
    spatialOne <- which(attr(object$ZAlist,"ranefs") == spatial.model) ## strictly a spatial one, not other correlated ones      
    oldLMatrix <- attr(object$predictor,"LMatrix") ## may be NULL
  }
  if (variances$linPred) { ## then we need X.pv; if newdata present then we also compute etaFix 
    X.pv <- allFrames$X  
  } else X.pv <- NULL
  if ( ! is.null(spatial.model)) {
    if (! is.null(newdata) && spatial.model[[1]]=="adjacency") {
      stop("Prediction in newdata not implemented or not possible in the 'adjacency' model")
    } ## FR->FR would be possible for new non-spatial predictor values in the original locations... il faudrait un test sur les elements de la distance matrix
  } 
  if ( ! is.null(spatial.model) &&
         ( ! is.null(newdata) || variances$linPred) ) { ## all cases where we need distance matrices
    blob <- calcNewCorrs(object=object,locdata=locdata,predVar=variances$linPred,spatial.model=spatial.model)
    uuCnewold <- blob$uuCnewold
    attr(uuCnewold,"ranefs") <- attr(oldLMatrix,"ranefs")  
    if (variances$linPred) Cnewnew <- blob$Cnewnew
  } 
  ## (1) computes fv (2) compute predVar
  ##### fv
  if (noReForm(re.form)) {
    fv <- object$family$linkinv(etaFix) 
  } else if (is.null(newdata)) {
    fv <- object$fv ## same length including replicates
    newZAlist <- object$ZAlist ## useful if predVar
  } else { ## 
    if ( nrand==0 ) {
      eta <- etaFix
      newZAlist <- NULL
    } else {
      FL <- spMMFactorList(locform, allFrames$mf, 0L, drop=TRUE) 
      newZAlist <- FL$Design
      #       ## previously to 1.3.11 a block of code matched the new matrices to the old v_h
      #       ## In 1.3.11 we match a new v_h to the new matrices (see augm_v_h_coeffs below)
      attr(newZAlist,"ranefs") <- attr(object$ZAlist,"ranefs") 
      if ( is.null(spatial.model)) {
        v_h_coeffs <- object$v_h ## v_h for non spatial and will later contain coeffs for spatial 
        ZALlist <- compute.ZALlist(LMatrix=NULL,ZAlist=newZAlist,Groupings=FL$Groupings) 
      } else {
        #### precomputation of coeffs
        ## on the gaussian scale, L.v_ori ~ lam C (lam C + phi I)^{-1}y 
        ## new random autocorr term ~ lam c (lam C + phi I)^{-1}y = c C^{-1} L_ori.v_ori = c [t(L_ori)]^{-1} v_ori
        ## [t(L_ori)]^{-1} v_ori can be computed once for all predictions => 'coeffs' argument
        if (is.null(coeffs)) { ## if not provided as argument
          v_h_coeffs <- predictionCoeffs(object) ## changes the coefficients in the right u_range
        } else v_h_coeffs <- coeffs
        ZALlist <- compute.ZALlist(CMatrix=uuCnewold,ZAlist=newZAlist,Groupings=FL$Groupings) ## ZAL's for ZA's and L's (typically some ZA's are unaffected)
      }  
      for (it in seq_len(length(ZALlist))) { 
        if ( ! is.matrix(ZALlist[[it]])) { ## this is not in compute ZALlist because the latter does not runs other all list elements
          cnames <- colnames(ZALlist[[it]]) ## used a few lines belwo to match levels of ranef
          rnames <- rownames(ZALlist[[it]]) ## used a few lines belwo to match levels of ranef
          ZALlist[[it]] <- as.matrix(ZALlist[[it]]) ## pfff  lose names
          colnames(ZALlist[[it]]) <- cnames
          rownames(ZALlist[[it]]) <- rnames
        }
      }
      if (nrand>1) {
        ZAL <- do.call(cbind,ZALlist)
      } else {
        ZAL <- ZALlist[[1]]
      }
      old_cum_n_u_h <- cumsum(c(0,attr(object$lambda,"n_u_h")))
      lcrandfamfam <- attr(object$`rand.families`,"lcrandfamfam")
      augm_v_h_coeffs <- lapply(seq_len(length(ZALlist)),function(it) {
        oldu.range <- (old_cum_n_u_h[it]+1L):(old_cum_n_u_h[it+1L])
        if (it == spatialOne) {
          return(v_h_coeffs[oldu.range])          
        } else {
          oldlevels <- colnames(object$ZAlist[[it]])
          newlevels <- colnames(ZALlist[[it]])
          interlevels <- intersect(oldlevels,newlevels)
          oldv <- v_h_coeffs[oldu.range]
          names(oldv) <- oldlevels
          psi_M <- switch(lcrandfamfam[it], 
                           gaussian = 0,
                           gamma = 1, 
                           beta = 1/2, 
                           "inverse.gamma" = 1
          )
          vpsi_M <- object$rand.families[[it]]$linkfun(psi_M)
          newv <- rep(vpsi_M,length(newlevels))
          names(newv) <- newlevels
          newv[interlevels] <- oldv[interlevels] 
          return(newv)
        }
      })
      augm_v_h_coeffs <- unlist(augm_v_h_coeffs)
      eta <- etaFix + ZAL %*% augm_v_h_coeffs ## (length(eta)) col vector from coeffs = length(eta) row vector...
    }
    # done with eta
    fv <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  }
  resu <- fv ## suitable for objective function of optim() etc ## matrix ! maybe more suitable than data frame as objective function
  if ( ! is.logical(binding) ) { ## expecting a string
    binding <- makenewname(base=binding,varnames=colnames(locdata)) ## FR->FR 09/11/2014 = posterior to CRAN 1.4.1 
    resu <- data.frame(resu)
    colnames(resu) <- binding
    resu <- cbind(locdata,resu) ## becomes a data frame !
    attr(resu,"fittedName") <- binding
  } else { ## expecting binding= FALSE
    if (ncol(locdata)>0)  attr(resu,"frame") <- locdata 
  }
  ##### (2) predVar
  if(variances$linPred) {
    ## list for Cnewnew, which enters in  newZA %*% Cnewnew %*% tnewZA, hence should not represent newZA itself 
    if (nrand>0) newnewClist <- lapply(seq_len(length(newZAlist)),function(it) {
      if (length(spatialOne)>0 && it == spatialOne) { ## avoids running the next algo which is slow on large matrices
        Cnn <- Cnewnew ## already computed for point prediction
      } else { 
        Cnn <- Diagonal(ncol(newZAlist[[it]])) ## diag(ncol(newZAlist[[it]]))
      }
      return(Cnn)
    })
    ## list for Coldnew, which enters in ZA %id*id% Coldnew %id*id% tnewZA 
    if (nrand>0L) oldnewClist <- lapply(seq_len(length(object$ZAlist)),function(it) {
      if (length(spatialOne)>0L && it == spatialOne) { ## avoids running the next algo which is slow on large matrices
        oldnewC <- t(uuCnewold)
      } else {
        oldlevels <- colnames(object$ZAlist[[it]])
        newlevels <- colnames(newZAlist[[it]])
        if (identical(oldlevels,newlevels)) {
          oldnewC <- Diagonal(length(oldlevels)) ## replaces old identityMatrix
        } else {
          oldornew <- unique(c(oldlevels,newlevels))
          oldnewC <- diag(length(oldornew))
          colnames(oldnewC) <- rownames(oldnewC) <- oldornew
          oldnewC <- oldnewC[oldlevels,newlevels]
        }
      }
      return(oldnewC)
    })
    if (nrand>0L) {
      oldX.pv <- object$`X.pv` 
      ZAL <- object$ZALMatrix 
      Sig <- ZWZtwrapper(ZAL, 1/object$w.ranef) + diag(1/object$w.resid)  ## same as in HLfit
      if (nrand==1L) {
        loclist <- list(phi=object$phi, lambda=object$lambda, Coldnew=oldnewClist[[1]],
                        Cnewnew=newnewClist[[1]], oldX.pv=oldX.pv, X.pv=X.pv, 
                        ZA=object$ZAlist[[1]], newZA=newZAlist[[1]],
                        beta_cov=object$beta_cov, Sig=Sig)
        ## but the code for nrand >1 should work for nrand==1:
      } else {
        loclist <- list(phi=object$phi, lambda=object$lambda, Coldnew=oldnewClist,
                        Cnewnew=newnewClist, oldX.pv=oldX.pv, X.pv=X.pv, 
                        ZA=object$ZAlist, newZA=newZAlist,
                        beta_cov=object$beta_cov, Sig=Sig)
      }
      predVarMat <- do.call(calcPredVar,loclist) ## matrix, not Matrix (partly required by call to Reduce in calcPredVar, and assumed below)
      rownames(predVarMat) <- colnames(predVarMat) <- rownames(locdata)
    } else predVarMat <- matrix(0,nrow=nrow(locdata),ncol=nrow(locdata))
    if (variances$cov) {
      sumVar <- predVarMat
    } else sumVar <- diag(predVarMat)
    attr(resu,"predVar") <- sumVar ## vector or matrix
  } else sumVar <- 0  
  if (variances$resid) {
    if (object$family$family %in% c("poisson","binomial")) {
      attr(resu,"residVar") <- object$family$variance(fv)
    } else attr(resu,"residVar") <- calcResidVar(object,newdata=locdata) 
    if (inherits(sumVar,"matrix")) {
      diag(sumVar) <- diag(sumVar) + attr(resu,"residVar")
    } else sumVar <- sumVar + attr(resu,"residVar")
  }
  if ( variances$fixef && ! is.null(object$beta_cov)) {
    fixefcov <- X.pv %*% object$beta_cov %*% t(X.pv)
    if (variances$cov) {
      attr(resu,"fixefVar") <- fixefcov 
    } else attr(resu,"fixefVar") <- diag(fixefcov)
    ## sumVar <- sumVar + attr(resu,"fixefVar") ## there is already such a term in predVar
  }
  if (variances$sum) attr(resu,"sumVar") <- sumVar
  return(resu)
}
