predict.HL <-
function(object,newX = NULL,coeffs=NULL,...) { ## but not new Y
  hlfit <- object$hlfit
  spatial.model <- findSpatial(hlfit$predictor)[[1]] 
  geonames <- object$HLCor.info$coordinates
  if (is.null(geonames)) { ## reconstructs the info
    bars <- spatial.model[[2]] 
    coordinates <- deparse(bars[[3]]) ## "x + y"
    coordinates <-  strsplit(coordinates," ")[[1]]
    geonames <- coordinates[coordinates != "+"]    
  }
  oldgeo <- hlfit$data[,geonames,drop=FALSE]
  if (is.null(newX)) { ## trivial but convenient option
    return(cbind(oldgeo,fv=hlfit$fv))
  } ## ELSE
  ########### newX not NULL (won't work with adjacency matrix...) #############
  if (spatial.model[[1]]=="adjacency") {
    stop("Prediction in new spatial positions not possible in the autoregressive model")
  }
  newgeo <- newX[,geonames,drop=FALSE]
  newMeanFrames <- HLframes(formula=attr(hlfit$predictor,"oriFormula"),data=newX) ## may need to reconstruct offset using formula term
  rho <- hlfit$ranFix$rho
  if (length(rho)==1) {
    oldgeoScal <- oldgeo * rho 
    newgeoScal <- newgeo * rho 
  } else if (ncol(oldgeo)==length(rho)) {
    oldgeoScal <-t(t(oldgeo) * rho) ## valid for vectorial rho...
    newgeoScal <-t(t(newgeo) * rho) ## valid for vectorial rho...
  } else {
    mess <- pastefrom("invalid length(rho).",prefix="(!) From ")
    stop(mess)
  }
  ## rows from newgeo cols from oldgeo:
  Cnew <-apply(oldgeoScal,1,function(v) {apply(newgeoScal,1,function(vv){ sqrt(sum((v-vv)^2))}) })
  args <-hlfit$ranFix[which(names(hlfit$ranFix) %in% c("nu","Nugget"))]
  Cnew <- do.call(Matern.corr,args=c(args,list(d=Cnew)))  
  ## reconstructs linear predictor 
  X.pv <- newMeanFrames$X  
  if (ncol(X.pv)>0) {eta <- X.pv %*% hlfit$fixef} else {eta <- rep(0,nrow(newX))}
  ## newX -> offset must be recomputed
  off <- model.offset(newMeanFrames$mf) ## Predictor has checked that there was either an offset term XOR an $offset 
  if ( is.null(off) ) { ## no offset oriFormula term
    off <- attr(hlfit$predictor,"offset") ## a PROCESSED predictor always has a non-NULL offset term 
    if (length(unique(off))>1) { ## that means there was a non trivial offset in the original call 
      message("Prediction in new locations with an offset from original positions is suspect.")
    }
  }
  eta <- eta + off 
  ## on the gaussian scale, L.v_ori ~ lam C (lam C + phi I)^{-1}y 
  ## new random autocorr term ~ lam c (lam C + phi I)^{-1}y = c C^{-1} L_ori.v_ori = c [t(L_ori)]^{-1} v_ori
  ## [t(L_ori)]^{-1} v_ori could be computed once for all predictions
  ## newX -> designL must be recomputed
  LMatrix <- designL.from.Corr(Cnew)
  if (is.null(coeffs)) coeffs <- solve(t(LMatrix),hlfit$v_h)   
  eta <- eta + Cnew %*% coeffs
  fv <- hlfit$family$linkinv(eta) ## ! freqs for binomial, counts for poisson
  return(cbind(newgeo,fv))
}
