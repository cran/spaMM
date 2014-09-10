predict.HLfit <-
function(object,newX = NULL,coeffs=NULL,predVar=FALSE,re.form = NULL,
                          binding = if(is.vector(newX)) {FALSE} else {"fitted"},...) { ## but not new Y
  if (predVar == "Cov") {
    predVar <- TRUE; predCov=TRUE
  } else predCov <- FALSE
  predictor <- object$predictor
  spatial.terms <- findSpatial(predictor)
  spatial.model <- spatial.terms[[1]] 
  if ( ! is.null(spatial.model) ) {
      olduniqueGeo <- attr(object,"info.uniqueGeo")
      geonames <- colnames(olduniqueGeo)
  } else geonames <- NULL
  allvars <- all.vars(attr(predictor,"oriFormula")[[3]]) 
  if (is.vector(newX)) { ## ## less well controlled case, but useful for maximization
    newX <- data.frame(matrix(newX,nrow=1))
    if (length(allvars)==ncol(newX)) {
      names(newX) <- allvars
    } else {
      stop(paste("(!) newX has incorrect length. It should match the following variables:\n",paste(allvars,collapse=" ")))
    }
  } 
  ## it is important that newX remains NULL if it was so initially because it is tested below
  if (is.null(newX)) {
    locdata <- object$data
  } else locdata <- newX
  locform <- attr(predictor,"oriFormula")
  if (noReForm(re.form)) locform <- nobarsMM(locform)
  allFrames <- HLframes(formula=locform[-2],data=locdata)
  if ( ( ! is.null(newX) ) || predVar || noReForm(re.form)) {
    etaFix <- newetaFix(object,allFrames)       
    oldLMatrix <- attr(predictor,"LMatrix") ## may be NULL
    if (noReForm(re.form)) {
      nrand <- 0
    } else {
      nrand <- length(attr(object$ZAlist,"ranefs"))
      spatialOne <- which(attr(object$ZAlist,"ranefs") == spatial.model) ## strictly a spatial one, not other correlated ones
    }
    if ( ! is.null(attr(predictor,"AMatrix"))) {
      stop("No appropriate code for prediction with an 'AMatrix'.")
    }
    if ( ! is.null(attr(predictor,"ZALMatrix"))) {
      stop("No appropriate code for prediction with a 'ZALMatrix'.")
    }
  }  
  if (predVar) { ## then we need X.pv; if newX present then we also compute etaFix 
    if (is.null(oldLMatrix)) {stop("no predVar code for models without spatial effects")} 
    X.pv <- allFrames$X  
  } else X.pv <- NULL
  if (is.null(newX)) {
    locdata <- object$data[,allvars,drop=FALSE]
  } else {
    locdata <- newX
  }
  if ( ! is.null(spatial.model)) {
    if (! is.null(newX) && spatial.model[[1]]=="adjacency") {
      stop("Prediction in newX not implemented or not possible in the 'adjacency' model")
    } ## FR->FR would be possible for new non-spatial predictor values in the original locations... il faudrait un test sur les elements de la distance matrix
  } 
  if ( ( ! is.null(newX)  && ! is.null(spatial.model)) 
      || predVar) { ## all cases where we need distance matrices
    newgeo <- locdata[,geonames,drop=FALSE] ## not unique
    newuniqueGeo <- unique(newgeo) ## It is essential that it keeps the same order as spMMfactorlist -> ULI -> unique. 
    rho <- object$ranFix$rho
    if (length(rho)==1) {
      oldgeoScal <- olduniqueGeo * rho 
      newgeoScal <- newuniqueGeo * rho 
    } else if (ncol(olduniqueGeo)==length(rho)) {
      oldgeoScal <-t(t(olduniqueGeo) * rho) ## valid for vectorial rho...
      newgeoScal <-t(t(newuniqueGeo) * rho) ## valid for vectorial rho...
    } else {
      mess <- pastefrom("invalid length(rho).",prefix="(!) From ")
      stop(mess)
    }
    ## rows from newgeo cols from olduniqueGeo:
    ## very slow: Cnew <-apply(as.matrix(oldgeoScal),1,function(v) {apply(newgeoScal,1,function(vv){ sqrt(sum((v-vv)^2))}) })
    ## faster but assumesvector: Cnew <- sqrt(colSums((apply(oldgeoScal,1,function(v){v-as.numeric(newgeoScal)}))^2))
    uuCnewold <- proxy::dist(newgeoScal,oldgeoScal) ## unique, unique
    if (predCov) {
      Cnewnew <- proxy::dist(newgeoScal,newgeoScal) ## dist matrix and then call to correl fn
    } else if (predVar) Cnewnew <-diag(rep(1,nrow(newgeoScal))) ## directly "corr" of u_h with itself at zero distance
    if (spatial.model[[1]]=="AR1") {
      args <-object$ranFix[which(names(object$ranFix) %in% c("ARphi"))]
      uuCnewold <- args$ARphi^uuCnewold  
      if (predCov) Cnewnew <- args$ARphi^Cnewnew  
    } else {
      args <-object$ranFix[which(names(object$ranFix) %in% c("nu","Nugget"))]
      uuCnewold <- do.call(Matern.corr,args=c(args,list(d=uuCnewold)))  
      if (predCov) Cnewnew <- do.call(Matern.corr,args=c(args,list(d=Cnewnew)))  
    }
    attr(uuCnewold,"ranefs") <- attr(oldLMatrix,"ranefs")  ## includes replicate locations !
  } 
  ## (1) computes fv (2) compute predVar
  ##### fv
  if (noReForm(re.form)) {
    fv <- object$family$linkinv(etaFix) 
  } else if (is.null(newX)) {
    fv <- object$fv ## same length including replicates
    newZAlist <- object$ZAlist ## useful if predVar
  } else { ## 
    if ( nrand==0 ) {
      eta <- etaFix 
    } else {
      FL <- spMMFactorList(locform, allFrames$mf, 0L, 0L) 
      tempZAlist <- FL$Design 
      ## next we convert these to design matrices matching the predicted v_h 
      newZAlist <- lapply(seq_len(length(tempZAlist)),function(it) {
        newlevels <- colnames(tempZAlist[[it]])
        if (it == spatialOne) { ## avoids running the next algo which is slow on large matrices
          newZA <- tempZAlist[[it]]
        } else {
          oldlevels <- colnames(object$ZAlist[[it]])
          template <- data.frame(matrix(rep(0,ncol(object$ZAlist[[it]])),nrow=1))
          newZA <- apply(tempZAlist[[it]],1, function(ro) {
            nonzero <- which(ro!=0)## only one 
            whichcol <- which(oldlevels==newlevels[nonzero]) ## identifies column in old ZAlist <=> level of v_h
            template[whichcol] <- ro[nonzero] ## local modif of template
            template
          })
          newZA <- do.call(rbind,newZA) 
          attr(newZA,"identityMatrix") <- is.identity(newZA)
          colnames(newZA) <- oldlevels ## 06/2014 not used...
        }
        as.matrix(newZA)
      })
      attr(newZAlist,"ranefs") <- attr(object$ZAlist,"ranefs") 
      if ( is.null(spatial.model)) {
        v_h_coeffs <- object$v_h ## v_h for non spatial and will later contain coeffs for spatial 
        ZALlist <- compute.ZALlist(LMatrix=NULL,ZAlist=newZAlist,Groupings=FL$Groupings)
      } else {
        #### precomputation of coeffs
        ## on the gaussian scale, L.v_ori ~ lam C (lam C + phi I)^{-1}y 
        ## new random autocorr term ~ lam c (lam C + phi I)^{-1}y = c C^{-1} L_ori.v_ori = c [t(L_ori)]^{-1} v_ori
        ## [t(L_ori)]^{-1} v_ori can be computed once for all predictions => 'coeffs' argument
        v_h_coeffs <- predictionCoeffs(object)
        ZALlist <- compute.ZALlist(LMatrix=uuCnewold,ZAlist=newZAlist,Groupings=FL$Groupings)
      }  
      if (nrand>1) {
        ZAL <- do.call(cbind,ZALlist)
      } else ZAL <- ZALlist[[1]]    
      eta <- etaFix + ZAL %*% v_h_coeffs ## (length(eta)) col vector from coeffs = length(eta) row vector...
    }
    # done with eta
    fv <- object$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  }
  resu <- fv ## suitable for objective function of optim() etc ## matrix ! maybe more suitable than data frame as objective function
  if ( ! is.logical(binding) ) { ## expecting a string
    colnames(resu) <- binding
    resu <- cbind(locdata,resu) ## becomes a data frame !    
  } else { ## expecting binding= FALSE
    if (ncol(locdata)>0)  attr(resu,"frame") <- locdata 
  }
  ##### (2) predVar
  if(predVar) {
    if (is.null(spatial.model)) {
      stop("No code for predVar without spatial random effects.")
    }
    if (nrand==1) {
      ZA <- object$ZAlist[[1]]
    } else {
      stop("No code for predVar with non-spatial random effects.")
    }
    oldX.pv <- object$X ## the X.pv of the fit of the original data ## FR->FR confusing names: X for X.pv (desgin fixed) et newX for new geo coord+pred vars
    ZAL <- ZA %id*id% oldLMatrix
    Sig <- ZWZt(ZAL,1/object$w.ranef) + diag(1/object$w.resid)  ## same as in HLfit
    loclist <- list(phi=object$phi,lambda=object$lambda,Coldnew=t(uuCnewold),
                    Cnewnew=Cnewnew,oldX.pv=oldX.pv,X.pv=X.pv,ZA=ZA,newZA=newZAlist[[spatialOne]],
                    beta_cov=object$beta_cov,Sig=Sig)
    predVarMat <- do.call(calcPredVar,loclist)
    if (predCov) {
      attr(resu,"predVar") <- predVarMat
    } else attr(resu,"predVar") <- diag(predVarMat) ## only the diag is correct
  } ## else attr(resu,"predVar") will be NULL  
  return(resu)
}
