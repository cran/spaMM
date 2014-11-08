## was CKrigcoefs in Rmigraine
# xy should have attribute fittedPars
# smoothnessRange is a new argument, no default to catch future problems  in assembling code => smoothnessRange=c(minSmoothness,maxSmoothness)
# nuniquerows was given by selectFn

ordinaryKrigingCoefs <- function(xy,nuniquerows,covfnparamA=NA,lambdaA=NA,hglmLambdaA=NA,hglmPhiA=NA,optimise=F,method="GCV",smoothnessRange){
  nrowxy<-nrow(xy)
  ncolxy<-ncol(xy)
  c<-rep(0,nuniquerows)
  d<-rep(0,nuniquerows)
  D<-rep(0,nuniquerows)
  u<-rep(0,nuniquerows)
  #if (covFamily=="Matern") {fneval<-0;} else fneval<-1
  fneval <- 1
  Tmatrix <- 1
  GCV<-0 ## input value controls Cross validation method; 0 for GCV, 1 for match_fs2hat_pure_error
  if (method =="HGLM") hglmP_bv<-0 ## some initialization is required...
  if (optimise) {
    covfnparam<-c(rep(0,length(attr(xy,"fittedPars"))),smoothnessRange[1]) ## 0=>dummy; initial values are set in x(ii) in CKrigcoefs' Krig.cpp. Last value is min value for smoothness
    lambda<-0 ## just a dummy value, not tested in C code
    if (method =="HGLM") {hglmLambda<-0;hglmPhi<-0} ## just dummy values, not tested in C code
    #optBool<-0
  } else {
    if (any(is.na(covfnparamA))) {
      mess <- pastefrom("NA in given covariance function parameters (including smoothness).",prefix="(!) From ")
      message(covfnparamA)
      stop(mess)
    }
    covfnparam<-covfnparamA
    if (method =="HGLM") {
      if (is.na(hglmLambdaA)) {
        mess <- pastefrom("given hglmLambdaA parameter is NA.",prefix="(!) From ")
        stop(mess)
      } else {hglmLambda<-hglmLambdaA;hglmPhi<-hglmPhiA} ## then returns NAs 
    } else {
      if (is.na(lambdaA)) {
        mess <- pastefrom("given lambda parameter is NA.",prefix="(!) From ")
        stop(mess)
      } else lambda<-lambdaA ## then returns NAs 
    }
    #optBool<-1
  }
  CKrigidx <- -1
  ## attention quand on rajoute un element aux arguments de CKrigcoefs il faut changer le return(...) aussi !!!
  if (method =="HGLM") {
    blob<-.C("CHGLMsmoothing",
             as.double(t(xy)), 
             as.integer(nrowxy),
             as.integer(ncolxy),
             as.integer(nuniquerows),
             as.double(smoothnessRange[2]), 
             as.double(c),
             as.double(d),
             #                  as.double(D),
             #                  as.double(u),
             as.double(hglmLambda),
             as.double(hglmPhi),
             as.double(hglmP_bv),
             as.double(covfnparam),  ## minSmoothness is passed as the last value of covfnparam.
             ## minSmoothness is then constrained to [1.001,4]
             ## init Smoothness is then determined as x(ii)=minSmoothness+2.*maxrange[ii-1]/3 ~ 3by default.
             as.integer(fneval),
             as.integer(optimise),
             as.integer(CKrigidx)
    ) 
    return(list(
      txy=blob[[1]],
      nrowxy=blob[[2]], ## should be input value
      ncolxy=blob[[3]], ## should be input value
      uniquerows=blob[[4]],
      maxSmoothness=blob[[5]],
      c=blob[[6]],
      d=blob[[7]],
      ##        D=blob[[8]],
      ##         u=blob[[9]],
      hglmLambda=blob[[8]],
      hglmPhi=blob[[9]],
      hglmP_bv=blob[[10]],
      covfnparam=blob[[11]],
      fneval=blob[[12]],
      optimise=blob[[13]],
      CKrigidx=blob[[14]]
    )) # is a list
  } else {
    blob<-.C("CKrigcoefs",
             as.double(t(xy)), 
             as.integer(nrowxy),
             as.integer(ncolxy),
             as.integer(nuniquerows),
             as.double(smoothnessRange[2]), 
             as.double(c),
             as.double(d),
             as.double(D),
             as.double(u),
             as.double(lambda),
             as.double(GCV),
             as.double(covfnparam),  ## minSmoothness is passed as the last value of covfnparam.
             ## minSmoothness is then constrained to [1.001,4]
             ## init Smoothness is then determined as x(ii)=minSmoothness+2.*maxrange[ii-1]/3 ~ 3by default.
             as.integer(fneval),
             as.integer(optimise),
             as.integer(CKrigidx)
    ) 
    return(list(
      txy=blob[[1]],
      nrowxy=blob[[2]], ## should be input value
      ncolxy=blob[[3]], ## should be input value
      uniquerows=blob[[4]],
      maxSmoothness=blob[[5]],
      c=blob[[6]],
      d=blob[[7]],
      D=blob[[8]],
      u=blob[[9]],
      lambda=blob[[10]],
      GCV=blob[[11]],
      covfnparam=blob[[12]],
      fneval=blob[[13]],
      optimise=blob[[14]],
      CKrigidx=blob[[15]]
    )) # is a list
  }
} ## end def ordinaryKrigingCoefs
