predict.HLfit <-
function(object,newX = NULL,coeffs=NULL,predVar=FALSE,re.form = NULL,
                          binding = if(is.vector(newX)) {FALSE} else {"fitted"},...) { ## but not new Y
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
  if ( ( ! is.null(newX) ) || predVar || noReForm(re.form)) {
    if (is.null(newX)) {
      newetaX <- newetaFix(object,object$data)       
    } else {
      newetaX <- newetaFix(object,newX)       
    }
    oldLMatrix <- attr(predictor,"LMatrix") ## may be NULL
    if (noReForm(re.form)) {
      nrand <- 0
    } else {
      nrand <- length(attr(object$ZAlist,"ranefs"))
    }
    if ( ! is.null(attr(predictor,"AMatrix"))) {
      stop("No appropriate code for prediction with an 'AMatrix'.")
    }
    if ( ! is.null(attr(predictor,"ZALMatrix"))) {
      stop("No appropriate code for prediction with an 'ZALMatrix'.")
    }
  }  
  if (predVar) { ## then we need X.pv; if newX present then we also compute etaFix 
    if (is.null(oldLMatrix)) {stop("no predVar code for models without spatial effects")} 
    X.pv <- newetaX$newMeanFrames$X  
  } else X.pv <- NULL
  if (is.null(newX)) {
    locdata <- object$data[,allvars,drop=FALSE]
  } else {
    locdata <- newX
    ## new a new design matrix for random effects, which columns do not match the predicted v_h
    allFrames <- HLframes(formula=attr(predictor,"oriFormula")[-2],data=locdata)
    FL <- spMMFactorList(predictor, allFrames$mf, 0L, 0L) 
    tempZAlist <- FL$Design 
    ## next we convert these to design matrices matching the predicted v_h 
    spatialOnes <- which(attr(object$ZAlist,"ranefs") %in% attr(oldLMatrix,"ranefs"))
    newZAlist <- lapply(seq_len(length(tempZAlist)),function(it) {
      newlevels <- colnames(tempZAlist[[it]])
      if (it %in% spatialOnes) { ## quick patch which makes newZA an identity matrix
        #oldlevels <- newlevels 
        #template <- data.frame(matrix(rep(0,nrow(locdata)),nrow=1)) ## Cnewold
        newZA <- diag(nrow(locdata))
        ## note that running the next algo on 41*41 locations is slow...
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
        colnames(newZA) <- oldlevels ## 06/2014 not used...
      }
      as.matrix(newZA)
    })
    attr(newZAlist,"ranefs") <- attr(object$ZAlist,"ranefs") 
  }
  if ( ! is.null(spatial.model)) {
    if (! is.null(newX) && spatial.model[[1]]=="adjacency") {
      stop("Prediction in newX not implemented or not possible in the 'adjacency' model")
    } ## FR->FR would be possible for new non-spatial predictor values in the original locations... il faudrait un test sur les elements de la distance matrix
  } 
  if ( ( ! is.null(newX)  && ! is.null(spatial.model)) 
      || predVar) { ## all cases where we need distance matrices
    newgeo <- locdata[,geonames,drop=FALSE] ## not unique       
    rho <- object$ranFix$rho
    if (length(rho)==1) {
      oldgeoScal <- olduniqueGeo * rho 
      newgeoScal <- newgeo * rho 
    } else if (ncol(olduniqueGeo)==length(rho)) {
      oldgeoScal <-t(t(olduniqueGeo) * rho) ## valid for vectorial rho...
      newgeoScal <-t(t(newgeo) * rho) ## valid for vectorial rho...
    } else {
      mess <- pastefrom("invalid length(rho).",prefix="(!) From ")
      stop(mess)
    }
    ## rows from newgeo cols from olduniqueGeo:
    ## very slow: Cnew <-apply(as.matrix(oldgeoScal),1,function(v) {apply(newgeoScal,1,function(vv){ sqrt(sum((v-vv)^2))}) })
    ## faster but assumesvector: Cnew <- sqrt(colSums((apply(oldgeoScal,1,function(v){v-as.numeric(newgeoScal)}))^2))
    Cnewold <- proxy::dist(newgeoScal,oldgeoScal) ## 
    if (predVar) Cnewnew <- proxy::dist(newgeoScal,newgeoScal)
    if (spatial.model[[1]]=="AR1") {
      args <-object$ranFix[which(names(object$ranFix) %in% c("ARphi"))]
      Cnewold <- args$ARphi^Cnewold  
      if (predVar) Cnewnew <- args$ARphi^Cnewnew  
    } else {
      args <-object$ranFix[which(names(object$ranFix) %in% c("nu","Nugget"))]
      Cnewold <- do.call(Matern.corr,args=c(args,list(d=Cnewold)))  
      if (predVar) Cnewnew <- do.call(Matern.corr,args=c(args,list(d=Cnewnew)))  
    }
    attr(Cnewold,"ranefs") <- attr(oldLMatrix,"ranefs")  ## includes replicate locations !
  } 
  ## (1) computes fv (2) compute predVar
  ##### fv
  if (is.null(newX)) {
    fv <- object$fv ## same length including replicates
  } else { ## 
    if ( nrand==0 ) {
      eta <- newetaX$etaFix 
    } else {
      etaFix <- newetaX$etaFix ## 
      v_h_coeffs <- object$v_h ## v_h for non spatial and will later contain coeffs for spatial 
      if ( is.null(spatial.model)) {
        ZALlist <- compute.ZALlist(LMatrix=NULL,ZAlist=newZAlist,Groupings=FL$Groupings)
      } else {
        ZALlist <- compute.ZALlist(LMatrix=Cnewold,ZAlist=newZAlist,Groupings=FL$Groupings)
        #### precomputation of coeffs
        ## on the gaussian scale, L.v_ori ~ lam C (lam C + phi I)^{-1}y 
        ## new random autocorr term ~ lam c (lam C + phi I)^{-1}y = c C^{-1} L_ori.v_ori = c [t(L_ori)]^{-1} v_ori
        ## FR->FR [t(L_ori)]^{-1} v_ori could be computed once for all predictions
        vec_n_u_h <- attr(object$lambda,"n_u_h")
        cum_n_u_h <- cumsum(c(0,vec_n_u_h))
        for (spit in spatialOnes) {
          u.range <- (cum_n_u_h[spit]+1L):(cum_n_u_h[spit+1L])
          if (is.null(coeffs)) coeffs <- solve(t(oldLMatrix),v_h_coeffs[u.range])   
          v_h_coeffs[u.range] <- coeffs 
        }
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
    Corr <- tcrossprodCpp(oldLMatrix)
    if (is.null(spatial.model)) {
      stop("No code for predVar without spatial random effects.")
    }
    if (nrand==1) {
      ZA <- object$ZAlist[[1]]
    } else {
      stop("No code for predVar with non-spatial random effects.")
    }
    oldX.pv <- object$X ## the X.pv of the fit of the original data ## FR->FR confusing names: X for X.pv (desgin fixed) et newX for new geo coord+pred vars
    predVar <- calcPredVar(object$phi,object$lambda,Corr,t(Cnewold),Cnewnew,oldX.pv,X.pv,ZA)
    attr(resu,"predVar") <- predVar
  } ## else attr(resu,"predVar") will be NULL  
  return(resu)
}
