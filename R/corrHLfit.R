dispFn <- function(x,xref=1) {log(x+xref)}
dispInv <- function(x,xref=1) {exp(x)-xref}
### transfos useful (test Vn phiFix...)
## for rho/trRho: rho=0 exactly is meaningful in adjacency model. rhoFn and rhoInv should be used only for other models.
rhoFn <- function(x) {log(x/(.spaMM.data$options$RHOMAX-x))} ## rho should be constrained to <RHOMAX and trRho should diverge as rho approaches RHOMAX
rhoInv <- function(trRho) {.spaMM.data$options$RHOMAX*exp(trRho)/(1+exp(trRho))} 
nuFn <- function(nu,rho) {log(nu/(.spaMM.data$options$NUMAX-nu))} ## nu should be constrained to <NUMAX and trNu should diverge as nu approaches NUMAX
nuInv <- function(trNu,trRho) {.spaMM.data$options$NUMAX*exp(trNu)/(1+exp(trNu))}
## FR->FR rho/sqrt(nu) scaling => should be fixed nu and transformed rho to handle vectorial rho
## thus nu fns should be indep of nu and rho fns should be functions of nu ! :
if (FALSE) {
 futurerhoFn <- function(rho,nu) {
  rhosc <- rho/sqrt(nu)
  log(rhosc/(.spaMM.data$options$RHOMAX-rhosc))
 }
 futurerhoInv <- function(trRhosc,trNu) {
  nu <- .spaMM.data$options$NUMAX*exp(trNu)/(1+exp(trNu))
  rhosc <- .spaMM.data$options$RHOMAX*exp(trRhosc)/(1+exp(trRhosc))
  rhosc*sqrt(nu)
 }
}


makeLowerUpper <- function(canon.init, ## cf calls: ~ in user scale, must be a full list of relevant params
                           lower,upper, ## ~in transformed scale
                           user.lower=list(),user.upper=list(),
                           corr.model="Matern",nbUnique,ranFix=list(),
                           lowerbound=list(),upperbound=list(),
                           control.dist=list(),
                           optim.scale) {
  ## init.optim not further used...
  if (corr.model %in% c("SAR_WWt","adjacency","ar1")) { ## adjacency model
    lower$rho <- user.lower$rho ## no transfo for adjacency model
    if (is.null(lower$rho)) lower$rho <- lowerbound$rho
    upper$rho <- user.upper$rho ## no transfo again
    if (is.null(upper$rho)) upper$rho <- upperbound$rho
  } else {
    if (corr.model=="AR1") {
      if ( ! is.null(canon.init$ARphi)) {
        ARphi <- user.lower$ARphi
        if (is.null(ARphi)) ARphi <- -0.999999
        lower$ARphi <- ARphi
        ARphi <- user.upper$ARphi
        if (is.null(ARphi)) ARphi <- 0.999999
        upper$ARphi <- ARphi
      }    
    } else { ## then Matern model....
      if (! is.null(canon.init$rho)) {
        rho <- user.lower$rho
        if (is.null(rho)) rho <- canon.init$rho/150
        if (optim.scale=="transformed") {
          lower$trRho <- rhoFn(rho)
        } else lower$rho <- rho
        rho <- user.upper$rho
        if (is.null(rho)) {
          if (inherits(nbUnique,"list")) nbUnique <- mean(unlist(nbUnique))
          rho <- canon.init$rho*2*nbUnique ## The following was a bit too low for experiments with nu=0.5 : 1/(maxrange/(2*nbUnique)) ## nb => unique rows !
          ## *modify* upper rho so that it does not exceed RHOMAX => /($RHOMAX+...)
          if (optim.scale=="transformed") rho <- 2*rho * .spaMM.data$options$RHOMAX/(.spaMM.data$options$RHOMAX+2*rho)
        }
        if (optim.scale=="transformed") {
          upper$trRho <- rhoFn(rho) 
        } else upper$rho <- rho
        rhoForNu <- canon.init$rho
      } else rhoForNu <- getPar(ranFix,"rho")
      if (! is.null(canon.init$nu)) {
        nu <- user.lower$nu
        if (is.null(nu)) nu <- canon.init$nu/100
        if (optim.scale=="transformed") {
          lower$trNu <- nuFn(nu,rhoForNu)
          #print(c(rhoForNu,nu,lower$trNu))
        } else lower$nu <-nu
        nu <- user.upper$nu
        if (is.null(nu)) {
          if ( ! is.null(dm <- control.dist$`dist.method`) && dm %in% c("Geodesic","Earth")) {
            nu <- 0.5
          } else {
            ## constructs upper nu from NUMAX => /(1+...)
            ## nu should not diverge otherwise it will diverge in Bessel_lnKnu, whatever the transformation used
            nu <- .spaMM.data$options$NUMAX * canon.init$nu/(1+canon.init$nu) 
          }
        }
        if (optim.scale=="transformed") {
          upper$trNu <- nuFn(nu,rhoForNu)
        } else upper$nu <- nu
        #print(c(rhoForNu,nu,upper$trNu))
      }
    } 
    ##### common to the different models except adjacency (because there are several places where NUgget+adjacency is not handled)
    if ( ! is.null(canon.init$Nugget)) {
      lower$Nugget <- 0
      upper$Nugget <- 0.999999
    }
  }
  if (! is.null(canon.init$phi)) {
    phi <- user.lower$phi
    if (is.null(phi)) phi <- canon.init$phi/1000
    lower$trPhi <- dispFn(phi)
    phi <- user.upper$phi
    if (is.null(phi)) phi <- canon.init$phi*1000
    ## if phi is badly initialized then it gets a default which may cause hard to catch problems in the bootstrap...
    upper$trPhi <- dispFn(phi)
  }
  if (! is.null(canon.init$lambda)) {
    lambda <- user.lower$lambda
    if (is.null(lambda)) lambda <- canon.init$lambda/1000
    lower$trLambda <- dispFn(lambda)
    lambda <- user.upper$lambda
    if (is.null(lambda)) lambda <- canon.init$lambda*1000
    upper$trLambda <- dispFn(lambda)
  }
  return(list(lower=lower,upper=upper))
}

checkDistMatrix <- function(distMatrix,data,coordinates) {
  classDistm <- class(distMatrix)
  if ( ! classDistm %in% c("matrix","dist")) {
    message(paste("(!) 'distMatrix' argument appears to be a '",classDistm,"',",sep=""))
    stop("not a 'matrix' or 'dist'. Check the input. I exit.")
  }
  ## chol() fails on distances matrices with repeated locations (which are pos SD)... but chol() not used by default
  ## the following code assumes that distMatrix deals only with unique locations, and checks this
  ## HENCE ******* distMatrix must refer to unique values of a grouping variable *********
  usernames <- rownames(distMatrix)
  checknames <- all(sapply(usernames,function(v) {v %in% rownames(data)})) ## 
  if (!checknames) {
    warning("The rownames of 'distMatrix' are not rownames of the 'data'. Further checking of 'distMatrix' is not possible.")
    nbUnique <- NA
  } else {
    uniqueGeo <- calcUniqueGeo(data=data[usernames,coordinates,drop=FALSE]) ## check that this corresponds to unique locations
    nbUnique <- nrow(uniqueGeo)
    if (nbUnique != nrow(distMatrix)) {
      stop("The dimension of 'distMatrix' does not match the number of levels of the grouping variable")
    } else { ## check order
      redondGeo <- data[,coordinates,drop=F]
      designRU <- apply(redondGeo,1,function(v) {which(apply(v==t(uniqueGeo),2,all))}) ## has no names
      ## eg 1 1 2 2 3 2 3 4 is valid for 8 obs, 4 unique locations
      designRU <- unique(as.vector(designRU)) ## should then be 1 2 3 4
      ## but if distMatrix in reverse order, the first row of redondGeo would match the 4th of uniqueGeo and then the following test is FALSE:
      if ( ! all (designRU==seq_len(length(designRU))) ) {
        stop("The rows of 'distMatrix' are not ordered as rows of the 'data'.")
      }
    } 
  }
  nbUnique ## if stop() did not occur
}

makeCheckGeoMatrices <- function(data,distMatrix=NULL,uniqueGeo=NULL,coordinates,dist.method="Euclidean") {
  isListData <- inherits(data,"list")
  if (is.null(distMatrix)) { 
    if ( is.null(uniqueGeo) ) { ## then construct it from the data ## this should be the routine case
      if (isListData) {
        uniqueGeo <- lapply(data,function(dd) {calcUniqueGeo(data=dd[,coordinates,drop=FALSE])})
        nbUnique <- lapply(uniqueGeo,nrow) 
      } else {
        uniqueGeo <- calcUniqueGeo(data=data[,coordinates,drop=FALSE])
        nbUnique <- nrow(uniqueGeo) 
      }
    } 
    ## (2): we need distMatrix *here* in all cases for the check
    if (isListData) {
      distMatrix <- lapply(uniqueGeo,proxy::dist,method=dist.method)
    } else distMatrix <- proxy::dist(uniqueGeo,method=dist.method)
  } else { ## there is a distMatrix, this is what will be used by HLCor
    if (isListData) {
      nbUnique <- lapply(seq_len(length(data)),function(dd) {checkDistMatrix(distMatrix,dd,coordinates)})
    } else nbUnique <- checkDistMatrix(distMatrix,data,coordinates)
    ## stops if problems, otherwise checkDistMatrix has no useful return value
  }
  return(list(nbUnique=nbUnique,uniqueGeo=uniqueGeo,distMatrix=distMatrix))
}




alternating <-function(init.optim,LowUp,anyOptim.args,maxIter,ranPars,HLCor.args,trace,Optimizer="L-BFGS-B",optimizers.args,maxcorners) {
  nam <- names(init.optim)
  if (any(c("trPhi","trLambda") %in% nam )) {
    mess <- pastefrom("Dispersion parameters non allowed in 'init.corrHLfit' with alternating algorithm.",prefix="(!) From ")
    stop(mess)
  }
  initcorr <- init.optim[nam %in% c("trRho","trNu","Nugget","ARphi")]
  HLfitLowUp <- LowUp
  HLfitLowUp$lower[c("trRho","trNu","Nugget","ARphi")] <- NULL
  HLfitLowUp$upper[c("trRho","trNu","Nugget","ARphi")] <- NULL
  corrLowUp <- LowUp
  corrLowUp$lower[c("trPhi","trLambda")] <- NULL
  corrLowUp$upper[c("trPhi","trLambda")] <- NULL
  anycorrOptim.args <- anyOptim.args
  iter <- 0
  conv <- 1
  currentLik <- -Inf
  while (iter < maxIter && conv > 1e-5 ) { ## if alternating: alternate HLCor and locoptim
    ranPars[names(initcorr)] <- initcorr
    attr(ranPars,"type")[names(initcorr)] <- "var"
    HLCor.args$ranPars <- ranPars
    oldLik <- currentLik
    if (is.character(trace$file)) {
      if(.spaMM.data$options$TRACE.UNLINK) unlink("HLCor.args.*.RData")
      zut <- paste(unlist(initcorr),collapse="")  
      save(HLCor.args,file=paste("HLCor.args.",zut,".RData",sep="")) ## for replicating the problem
    }
    givencorr <- do.call("HLCor",HLCor.args) ## optim disp and beta given corr param
    currentLik <- givencorr$APHLs$p_v ## iterations maximize p_v
    conv <- currentLik-oldLik
    anycorrOptim.args$ranPars$lambda <- givencorr$lambda
    anycorrOptim.args$ranPars$phi <- givencorr$phi
    anycorrOptim.args$etaFix <- list(beta=givencorr$fixef) ## LeeN01sm imply that eta should be fixed too 
    loclist <- list(initcorr,corrLowUp,anyObjfnCall.args=anycorrOptim.args,trace,Optimizer=Optimizer,
                    optimizers.args=optimizers.args,maxcorners=maxcorners,maximize=TRUE) 
    initcorr <- do.call("locoptim",loclist) 
    iter <- iter+1
  }
  optPars <- c(initcorr,givencorr$lambda,givencorr$phi)
  optPars
}

scalefnfn <- function(st) {switch(st,
                                  trRho="rhoFn",
                                  trNu="nuFn",
                                  trLambda="dispFn",
                                  "identity" ## a base:: function
)}
invscalefnfn <- function(st) {switch(st,
                                     trRho="rhoInv",
                                     trNu="nuInv",
                                     trLambda="dispInv",
                                     "identity" ## a base:: function
)}
validRangefn <- function(st) {switch(st,
                                     trRho=c(-Inf,.spaMM.data$options$RHOMAX),
                                     trNu=c(-Inf,.spaMM.data$options$NUMAX),
                                     c(-Inf,Inf)
)}
labelfn <- function(st) {switch(st,
                                trRho=expression(rho),
                                trNu=expression(nu),
                                trLambda=expression(lambda),
                                st
)}


# plots Krigobj$data$logLobj i.e. *not* the smoothing
SEMdiagnosticPlot2D <- function(Krigobj,smoothingOK,titlemain,titlesub, nextpoints, info, optrPar) {
  axesNames <- names(optrPar) 
  ticks_x <- do.call(invscalefnfn(axesNames[1]),list(Krigobj$data[,axesNames[1]]))
  tmp <- makeTicks(ticks_x, scalefn=scalefnfn(axesNames[1]),axis=1,logticks=TRUE,validRange=validRangefn(axesNames[1]))
  xlabs <- tmp$labels;xat <- tmp$at; #xph <- tmp$phantomat
  ticks_y <- do.call(invscalefnfn(axesNames[2]),list(Krigobj$data[,axesNames[2]]))
  tmp <- makeTicks(ticks_y, scalefn=scalefnfn(axesNames[2]),axis=2,logticks=TRUE,validRange=validRangefn(axesNames[2]))
  ylabs <- tmp$labels;yat <- tmp$at; #yph <- tmp$phantomat
  spaMMplot2D(Krigobj$data[,axesNames[1]],Krigobj$data[,axesNames[2]],Krigobj$data$logLobj,
              plot.title={
                title(main=titlemain,line=2,xlab=labelfn(axesNames[1]),ylab=labelfn(axesNames[2])) ## default cex.main=1.2, line ~1.7
                title(main=titlesub,line=0.8,cex.main=1.1) ## 
              },
              plot.axes={
                axis(1, at=xat, labels=xlabs)
                axis(2, at=yat, labels=ylabs)
              },
              decorations= {
                if ( ! is.null(nextpoints)) points(nextpoints,pch=15,cex=0.4)
                if (smoothingOK) apply(info$simplicesTable,1,function(v){
                  polygon(info$vertices[v,])
                });
                points(matrix(optrPar,nrow=1),pch="+",col="#008080")  ## col= colortools::complementary(spaMM.colors()[64])
              },
              map.asp=1/2 ## -(1-sqrt(5))/2 not wide enough
  )      
}

SEMdiagnosticPlot <- function(Krigobj,titlemain, optr) {
  nPredictors <- length(optr$par)
  intsqrt <- as.integer(sqrt(nPredictors))
  oldpar <- par(no.readonly=TRUE,mfrow=c(ceiling(nPredictors/intsqrt), intsqrt))
  for (st in names(optr$par)) {
    invscalefn <- invscalefnfn(st)
    ticks_x <- do.call(invscalefn,list(Krigobj$data[,st]))
    tmp <- makeTicks(ticks_x,scalefn=scalefnfn(st),axis=1,logticks=TRUE,validRange=validRangefn(st))
    xlabs <- tmp$labels;xat <- tmp$at
    plot(Krigobj$data[,st],Krigobj$data$logLobj,
           axes=FALSE,frame=TRUE,xlab=labelfn(st),ylab="log(L)"
    )
    axis(1, at=xat, labels=xlabs)
    axis(2)
    points(optr$par[st],optr$value,pch=19,cex=2,col="red")  
  }
  title( titlemain, outer = TRUE )
  par(oldpar)                              
  invisible(NULL)
}



## wrapper for optimization of HLCor.obj OR (iterateSEMSmooth -> HLCor directly)
corrHLfit <- function(formula,data, ## matches minimal call of HLfit
                      init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      trace=list(file=NULL,append=TRUE),
                      objective="p_bv", ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                      control.dist=list(),
                      control.corrHLfit=list(), ## alternating, optim.scale, Optimizer, optimizer.args, maxIter, maxcorners, 
                      processed=NULL, ## added 2014/02 for programming purposes
                      family=gaussian(),
                      ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) { 
  mc <- match.call() 
  init.optim <- init.corrHLfit ## usages of init.corrHLfit$rho, etc. must be tracked through init.optim too  
  alternating <- control.corrHLfit$alternating 
  if (is.null(alternating)) alternating <- FALSE
  optim.scale <- control.corrHLfit$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  Optimizer <- control.corrHLfit$Optimizer ## either "nlminb" or one of the methods of optim()
  if (is.null(Optimizer)) Optimizer="L-BFGS-B" ## default for locoptim but still requested for "Optimizer=Optimizer"
  optimizers.args <- control.corrHLfit[c("nlminb","optim","optimize")] 
  maxcorners <- control.corrHLfit$maxcorners
  ## 2015/12/24 : setting maxcorners to 0L till makes a difference for
  # test-Nugget -> different p_v / p_bv compromise (!)
  # nullfit of blackcap (test-fixedLRT) is better with the corners.
  if (is.null(maxcorners)) maxcorners <- 2^11 
  maxIter <- control.corrHLfit$maxIter 
  if (is.null(maxIter)) maxIter<- 10000
  if ( ! (objective %in% c("p_v","p_bv"))) {
    mess <- pastefrom("invalid value of the 'objective' argument.",prefix="(!) From ")
    stop(mess)
  }
  Fixrho <- getPar(ranFix,"rho")
  if ( (! is.null(Fixrho)) && (! is.null(init.corrHLfit$rho)) ) {
    stop("(!) 'rho' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  } else rho.size <- max(length(Fixrho),length(init.corrHLfit$rho))
  if ( (! is.null(getPar(ranFix,"nu"))) && (! is.null(init.corrHLfit$nu)) ) {
    stop("(!) 'nu' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  }
  if ( (! is.null(getPar(ranFix,"ARphi"))) && (! is.null(init.corrHLfit$ARphi)) ) {
    stop("(!) 'ARphi' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  }
  if ( (! is.null(getPar(ranFix,"Nugget"))) && (! is.null(init.corrHLfit$Nugget)) ) {
    stop("(!) 'Nugget' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  }
  if (length(trace)>0 && ! all(names(trace) %in% c("file","append"))) {
    mess <- pastefrom("'trace' elements should be named and = 'file','append'.",prefix="(!) From ")
    message(mess)
  } 
  if (is.character(trace$file) && ( ! trace$append) ) {
    try(unlink(trace$file))
  }  ## the file is written in by HLCor()                   
  dotlist <-list(...) 
  ## Preventing obsolete options
  if (!is.null(dotlist$rho.mapping)) {
    stop("'rho.mapping' is obsolete. Use control.dist$rho.mapping instead.")
  }
  if (!is.null(dotlist$ranPars)) {
    stop("incorrect 'ranPars' argument in corrHLfit call. Use ranFix (ranPars is for HLCor only)")
  }
  if (!is.null(dotlist$LamFix)) {
    stop("obsolete LamFix argument in corrHLfit call")
  }
  if (!is.null(dotlist$PhiFix)) {
    stop("obsolete PhiFix argument in corrHLfit call")
  }  
  predictor <- formula
  if (! inherits(predictor,"predictor")) predictor <- Predictor(formula) 
  spatial.terms <- findSpatial(predictor)
  spatial.model <- spatial.terms[[1]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1]]) 
  } else {
    # corr.model <- "Matern" ## up to v 1.0; for the defective syntax of the scripts for the Ecography paper
    stop("spatial correlation model not specified in 'formula': was valid in version 1.0 but not later.")
  } ## FR->FR ce test n'a peut être plsu de sens
  if (!is.null(HLmethod <- dotlist$HLmethod)) {
      if (HLmethod=="SEM") dotlist$`try.chol` <- FALSE
  } ## FR->FR ??? why not handled by preprocess ? actuellement longue construction de HLCor.args puis copie dans anyHLCor_obj_args apres preprocess
  HLCor.formals <- names(formals(HLCor))
  HLfit.formals <- names(formals(HLfit))
  designL.formals <- names(formals(designL.from.Corr))
  makescaled.formals <- names(formals(make_scaled_dist))
  HLnames <- (c(HLCor.formals,HLfit.formals,designL.formals,makescaled.formals))  ## cf parallel code in HLCor.obj
  argcheck <- names(dotlist)[which(! names(dotlist) %in% HLnames)]
  if (length(argcheck)>0) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in corrHLfit call."))
  HLCor.args <- dotlist[intersect(names(dotlist),HLnames)] ## those specifically for HLCor and those for HLfit and designL.from.Corr
  typelist <- list() ## on veut une list pour pouvoir supp des elements par <- NULL
  typelist[names(ranFix)] <- "fix"
  attr(ranFix,"type") <- typelist
  ranPars <- ranFix ## ranPars argument of HLCor contains both fixed and estimated parameters:
  HLCor.args$trace <- trace$file
  HLCor.args$formula <- predictor
  HLCor.args$corr.model <- corr.model
  verbose <- dotlist$verbose 
  if (is.null(verbose)) verbose <- logical(0)
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["SEM"])) verbose["SEM"] <- FALSE
  if (is.na(verbose["warn"])) verbose["warn"] <- FALSE ## important!
  if (is.na(verbose["summary"])) {
    verbose["corrHLfitSummary"] <- TRUE
  } else verbose["corrHLfitSummary"] <- verbose["summary"]
  verbose["summary"] <- FALSE ## this is for HLCor
  HLCor.args$verbose <- verbose
  ## one more modif of HLCor.args : distMatrix or uniqueGeo below
  #  rhoobj.formals <- names(formals(rho.obj)) ## list with default values !
  #  rhoobj.args <- dotlist[intersect(names(dotlist),rhoobj.formals)]
  ############### (almost) always check geo info ###################
  lowerbound <- list()
  upperbound <- list()
  if (corr.model  %in% c("SAR_WWt","adjacency","ar1")) {
    ## adjMatrix should become weightMatrix
    if ( is.null(HLCor.args$adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
    if (isSymmetric(HLCor.args$adjMatrix)) {
      decomp <- selfAdjointWrapper(HLCor.args$adjMatrix)
      attr(HLCor.args$adjMatrix,"symSVD") <- decomp
    } else {
      if (corr.model  %in% c("SAR_WWt")) {
        decomp <- eigen(HLCor.args$adjMatrix,symmetric=FALSE) ## FR->FR not RcppEigen-optimized
        attr(HLCor.args$adjMatrix,"UDU.") <- list(u=decomp$vectors,d=decomp$values,u.=solve(decomp$vectors))
      } else stop("'adjMatrix' is not symmetric") ## => invalid cov mat for MVN
    }
    rhorange <- sort(1/range(decomp$d)) ## keeping in mind that the bounds can be <>0
    if(verbose["SEM"]) cat(paste("Feasible rho range: ",paste(signif(rhorange,6),collapse=" -- "),"\n"))
    if (is.null(getPar(ranFix,"rho")) && (! is.numeric(init.HLfit$rho))) {
      lowerbound$rho <- rhorange[1]
      upperbound$rho <- rhorange[2]
    } ## else the bounds should remain empty otherwise makeLowerUpper willuse them to compute lower$rho etc...     
  } else {
    if ( is.null(spatial.model)) {
      stop("An obsolete syntax for the adjacency model appears to be used.")
      ## coordinates <- c("x","y") ## backward compatibility... old syntax with (1|pos) and default values of the coordinates argument
    } else {
      if ( inherits(data,"list")) {
        dataForCheck <- data[[1]]
      } else dataForCheck <- data
      coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(dataForCheck))
    }
  } ## 
  #########################################################
  family <- checkRespFam(family)
  famfam <- family$family
  if ( ! is.null(famfam) && famfam=="multi") {
    if ( ! inherits(data,"list")) {
      familyargs <- family
      familyargs$family <- NULL
      familyargs$binfamily <- NULL
      data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
    } ## then always a list:
  }
  if ( inherits(data,"list")) {
    data <- lapply(data,function(dt) {
      validdata <- validData(formula=formula,resid.formula=dotlist$resid.formula,data=dt) ## will remove rows with NA's in required variables
      dt[rownames(validdata),,drop=FALSE] ## ## before Predictor is called and an LMatrix is added, etc.
    })
  } else {
    validdata <- validData(formula=formula,resid.formula=dotlist$resid.formula,data=data) ## will remove rows with NA's in required variables
    if (!inherits(data,"environment")) {
      data <- data[rownames(validdata),,drop=FALSE] ## 
    } else data <- validdata
  }
  HLCor.args$family <- family
  HLCor.args$data <- data
  ## fills init.optim with all necessary values. There must be values for all parameters that are to be optimized 
  init <- list() ## will keep the initial values in untransformed scale
  if ( corr.model  %in% c("SAR_WWt","adjacency","ar1") ) { 
    if (is.null(getPar(ranFix,"rho")) && (! is.numeric(init.HLfit$rho))) {
      init$rho <- init.optim$rho 
      if (is.null(init$rho)) {
        lr <- lower$rho
        if (is.null(lr)) lr <- lowerbound$rho 
        ur <- upper$rho
        if (is.null(ur)) ur <- upperbound$rho 
        init$rho <- (lr+ur)/2
      }
      init.optim$rho <- init$rho
    }
    if ( (! is.null(init.HLfit$rho)) && (! is.numeric(init.HLfit$rho))) { ## init.HLfit$rho is NA
      init.HLfit$rho <- init.optim$rho
      init.optim$rho <- lowerbound$rho <- upperbound$rho <- NULL ## makeLowerUpper will use non-NULL lowerbound...
    }
    nbUnique <- NULL
  } else { ## if (corr.model %in% c("Matern","corMatern","AR1"))
    locarglist<- list(data=data,distMatrix=dotlist$distMatrix,
                      uniqueGeo=HLCor.args$uniqueGeo,coordinates=coordinates)
    if(!is.null(dist.method <- control.dist$dist.method)) locarglist$dist.method <- dist.method
    geoMats <- do.call(makeCheckGeoMatrices,locarglist)
    nbUnique <- geoMats$nbUnique  
    distMatrix <- geoMats$distMatrix   
    uniqueGeo <- geoMats$uniqueGeo   
    rho.mapping <- control.dist$rho.mapping ## may be NULL
    if (rho.size<2) { ## can be 0 if no explicit rho in the input  
      HLCor.args$distMatrix <- distMatrix   
    } else {
      HLCor.args$uniqueGeo <- uniqueGeo 
    }
    if (corr.model=="AR1") {
      if (is.null(getPar(ranFix,"ARphi")) && (! is.numeric(init.HLfit$ARphi))) { 
        init$ARphi <- init.optim$ARphi 
        if (is.null(init$ARphi)) init$ARphi <- 0. 
        if (! is.null(init.HLfit$ARphi)) {
          init.HLfit$ARphi <- init$ARphi 
        } else {
          init.optim$ARphi <- init$ARphi
        }
      } 
      if ( (! is.null(init.HLfit$ARphi)) && (! is.numeric(init.HLfit$ARphi))) {
        init.HLfit$ARphi <- init.optim$ARphi
        init.optim$ARphi <- NULL
      }      
    } else { ## there is a scale (rho) parameter
      if (rho.size==0L) {
        if ( ! is.null(rho.mapping)) {rho.size <- length(unique(rho.mapping))}
      }
      if (rho.size<2L) { ## can be 0 if no explicit rho in the input  
        if (inherits(distMatrix,"list")) {
          maxrange <- max(unlist(lapply(distMatrix,function(dd) max(c(-Inf,dd))))) ## les Inf to handle dist(0)...
                     -min(unlist(lapply(distMatrix,function(dd) min(c( Inf,dd)))))
        } else maxrange <- max(distMatrix)-min(distMatrix)        
      } else { ## rho.size >1
        if (inherits(uniqueGeo,"list")) {
          if (is.null(rho.mapping)) {
            if (ncol(uniqueGeo[[1]])==rho.size) { ## ncol>1, rho of length ncol
              rho.mapping <- seq_len(rho.size)           
            } else stop("'control.dist$rho.mapping' missing with no obvious default from the other arguments.")
          }
          maxrange <- lapply(unique(rho.mapping), function(idx) {
            ranges <- matrix(unlist(lapply(uniqueGeo,function(uu){
              if (nrow(uu)>1) {
                range(proxy::dist(uu[,rho.mapping==idx],method=dist.method))
              } else c(Inf,-Inf) ## encore des Inf to handle dist(0)...
              })),ncol=2)
            max(ranges[,2])-min(ranges[,1]) 
          })
        } else { ## single data set
          if (is.null(rho.mapping)) {
            if (ncol(uniqueGeo)==rho.size) {
              rho.mapping <- seq_len(rho.size)           
            } else stop("'rho.mapping' missing with no obvious default from the other arguments.")
          }
          maxrange <- lapply(unique(rho.mapping), function(idx) {
            rng <- range(proxy::dist(uniqueGeo[,rho.mapping==idx],method=dist.method))
            rng[2]-rng[1] 
          })
        }  
        maxrange <- unlist(maxrange)
        HLCor.args$`control.dist` <- control.dist
      }
      if (is.null(getPar(ranFix,"rho")) && (! is.numeric(init.HLfit$rho))) {
        init$rho <- init.optim$rho 
        if (is.null(init$rho)) {
          init$rho <- 30/(2*maxrange)
        } else if (any( narho <- is.na(init$rho))) init$rho[narho] <- 30/(2*maxrange[narho]) ## 05/2015 allows init.corrHLfit=list(rho=rep(NA,...
        if (! is.null(init.HLfit$rho)) { ## ! is.null, but ! is.numeric as tested above; but now it becomes numeric
          init.HLfit$rho <- init$rho ## avant transformation
        } else {
          if (optim.scale=="transformed") {
            init.optim$trRho <- rhoFn(init$rho) ## we're in Matern model here
            init.optim$rho <- NULL
          } else init.optim$rho <- init$rho
        }
      } 
      if ( (! is.null(init.HLfit$rho)) && (! is.numeric(init.HLfit$rho))) {
        init.HLfit$rho <- init.optim$rho
        init.optim$rho <- NULL
      }
      if (is.null(getPar(ranFix,"nu")) && (! is.numeric(init.HLfit$nu))) { 
        init$nu <- init.optim$nu 
        if (is.null(init$nu)) {
          if ( ! is.null(dm <- control.dist$`dist.method`) && dm %in% c("Geodesic","Earth")) {
            init$nu <- 0.25  
          } else init$nu <- 0.5 
        }
        if (! is.null(init.HLfit$nu)) {
          init.HLfit$nu <- init$nu ## avant transformation
        } else {
          if (optim.scale=="transformed") {
            Fixrho <- getPar(ranFix,"rho")
            if (is.null(Fixrho)) { 
              init.optim$trNu <- nuFn(init$nu,init$rho) 
            } else init.optim$trNu <- nuFn(init$nu,Fixrho)
            init.optim$nu <- NULL
          } else init.optim$nu <- init$nu
        }
      } 
      if ( (! is.null(init.HLfit$nu)) && (! is.numeric(init.HLfit$nu))) {
        init.HLfit$nu <- init.optim$nu
        init.optim$nu <- NULL
      }
    }
    if (is.null(getPar(ranFix,"Nugget"))) { init$Nugget <- init.optim$Nugget }  ## this may be null, but in this case we leave it so and the Nugget keeps its default value through all computations
  }
  if (is.null(getPar(ranFix,"lambda"))) { ## no ranFix$lambda: process init.optim
    init$lambda <- init.optim$lambda 
    if (!is.null(init$lambda)) {
      if (init$lambda<1e-4) init$lambda <- 1e-4
      init.optim$trLambda <- dispFn(init$lambda) 
      init.optim$lambda <- NULL
    }
  } else { ## ranFix$lambda present, do NOT put it in init.optim
    if (!is.null(init.optim$lambda)) stop("(!) Arguments 'ranFix$lambda' and 'init.corrHLfit$lambda' conflict with each other.")  
  } 
  if (is.null(getPar(ranFix,"phi"))) {
    init$phi <- init.optim$phi 
    if (!is.null(init$phi)) {
      if (init$phi<1e-4) init$phi <- 1e-4
      init.optim$trPhi <- dispFn(init$phi)
      init.optim$phi <- NULL
    }
  } else {
    if (!is.null(init.optim$phi)) stop("(!) Arguments 'ranFix$phi' and 'init.corrHLfit$phi' conflict with each other.")  
  } 
  ################
  ## done with init.optim
  # { ## ad hoc maximization of p_bv VIA HLCor.obj
  ## by maxim over (corrpars,phi,lambda, (beta)_[corrpars,phi,lambda])
  ##     if trPhi,trLambda are in the init.optims
  ## or by maxim over (corrpars,(beta,phi,lambda)_corrpars)
  ##     otherwise.
  ################ construct intervals for this maximization
  ## construct default upper and lower values ; on transformed scale by default
  user.lower <- lower; user.upper <- upper ## keep user input 
  if ("lambda" %in% c(names(user.lower),names(user.lower)) 
      && is.null(init$lambda)) {
    stop("'lambda' in 'lower' or 'upper' has no effect if absent from 'init.corrHLfit'.")
  }
  ################
  LUarglist <- list(canon.init=init,
                    lower=init.optim, upper=init.optim, ## initially with right transformed variables but wrong values
                    user.lower=user.lower,user.upper=user.upper,
                    corr.model=corr.model,nbUnique=nbUnique,
                    ranFix=ranFix,control.dist=control.dist,
                    optim.scale=optim.scale)
  if (corr.model %in% c("SAR_WWt","adjacency","ar1") ) {
    LUarglist$lowerbound <- lowerbound
    LUarglist$upperbound <- upperbound
  }
  LowUp <- do.call("makeLowerUpper",LUarglist)
  ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  varNames <- names(init.HLfit) ## hence those that will be variable within HLfit
  ranPars[varNames] <- init.HLfit[varNames] ## FR->FR duplicat (?) qui montre qu'un attribute serait mieux
  attr(ranPars,"type")[varNames] <- "var"  
  ################
  ## common stuff to both optim and optimize
  HLfit.formal.args <- formals(HLfit) ## makes sure about default values 
  nondefault.HLfit.formal.args <- HLCor.args[which(names(HLCor.args) %in% HLfit.formals)]
  HLfit.formal.args[names(nondefault.HLfit.formal.args)] <- nondefault.HLfit.formal.args ## full HLfit args
  if (is.null(processed)) {
    preprocess.formal.args <- HLfit.formal.args[which(HLfit.formals %in% names(formals(preprocess)))] ## only those for preprocess
    preprocess.formal.args$rand.families <- HLfit.formal.args$rand.family ## because preprocess expects $rand.families 
    preprocess.formal.args$predictor <- HLfit.formal.args$formula ## because preprocess stll expects $predictor 
    preprocess.formal.args$resid.predictor <- HLfit.formal.args$resid.formula ## because preprocess stll expects $predictor 
    if (family$family %in% c("poisson","binomial","multi")) { ## multi added 2015/12/12
      phi.Fix <- 1 
    } else {
      phi.Fix <- getPar(ranFix,"phi")
      if (any(phi.Fix==0)) {
        mess <- pastefrom("phi cannot be fixed to 0.",prefix="(!) From ")
        stop(mess)
      }
    } 
    preprocess.formal.args$phi.Fix <- phi.Fix
    preprocess.formal.args$prior.weights <- dotlist$prior.weights
    processed <- do.call("preprocess",preprocess.formal.args)
    if ( ! is.null(attr(processed,"multiple"))) {
      pnames <- names(processed[[1]])
    } else pnames <- names(processed)
    for (st in pnames) HLCor.args[st] <- NULL ## this keeps the data in HLCor.args as they are not returned in processed, but predictor is suppressed
    HLCor.args$HLmethod <- NULL ## because the processed<...>$HL element must be used 
  } else preprocess.formal.args <- as.list(getProcessed(processed,"callargs"))[-1] ## we'll use them in SEM code
  HLCor.args$processed <- processed
  HLCor.args$processed <- setProcessed(HLCor.args$processed,"SEMargs$SEMseed",
                                              value="NULL") ## removing default SEMseed
  HLCor.args$ranPars <- ranPars
  HLCor.args$control.dist <- control.dist
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  anyHLCor_obj_args$skeleton <- init.optim ## logscale, only used by HLCor.obj
  attr(anyHLCor_obj_args$skeleton,"type") <- list() 
  attr(anyHLCor_obj_args$skeleton,"type")[names(init.optim)] <- "fix" 
  anyHLCor_obj_args$`HLCor.obj.value` <- objective ## p_v when fixedLRT-> corrHLfit
  anyHLCor_obj_args$traceFileName <- trace$file
  ## optim/optimize specific code
  initvec <- unlist(init.optim)
  ####    tmpName <- generateName("HLtmp") ## tmpName is a string such as "HLtmp0"
  #    anyHLCor_obj_args$init.HLfit <- tmpName 
  ####    assign(tmpName,list(),pos=".GlobalEnv") ## sets HLtmp0 (or a similarly named variable) at the global level
  processedHL1 <- getProcessed(processed,"HL[1]") ## there's also HLmethod in processed<[[]]>$callargs
  if (!is.null(processedHL1) && processedHL1=="SEM" && length(lower)>0) {   
    ## : bc SEMseed OK to control individual SEMs but not  series of SEM 
    optr <- iterateSEMSmooth(processed=processed, anyHLCor_obj_args=anyHLCor_obj_args, 
                             LowUp=LowUp,init.corrHLfit=init.corrHLfit, 
                             preprocess.formal.args=preprocess.formal.args, 
                             control.corrHLfit=control.corrHLfit)
    optPars <- as.list(optr$par)
    if (!is.null(optPars)) attr(optPars,"method") <-"optimthroughSmooth"
  } else if (alternating) { ## renewed coding of the iterative algo (only p_v); not documented => not checked for a long time
    optPars <- alternating(init.optim=init.optim,LowUp=LowUp,
                           anyHLCor_obj_args=anyHLCor_obj_args,maxIter=maxIter,
                           ranPars=ranPars,HLCor.args=HLCor.args,trace=trace,Optimizer=Optimizer,
                           optimizers.args=optimizers.args,maxcorners=maxcorners)
    if (!is.null(optPars)) attr(optPars,"method") <-"alternating"
  } else { ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null 
    loclist<-list(init.optim=init.optim,LowUp=LowUp,anyObjfnCall.args=anyHLCor_obj_args,trace=trace,Optimizer=Optimizer,
                  optimizers.args=optimizers.args,maxcorners=maxcorners,maximize=TRUE) 
    optPars <- do.call("locoptim",loclist)
    if (!is.null(optPars)) attr(optPars,"method") <- "locoptim"
  }
  ranPars[names(optPars)] <- optPars ## avoids overwriting fixed ranPars 
  attr(ranPars,"type")[names(optPars)] <- "outer" ##  
  HLCor.args$ranPars <- ranPars
  verbose["warn"] <- TRUE ## important!
  HLCor.args$verbose <- verbose ##
  hlcor <- do.call("HLCor",HLCor.args) ## recomputation post optimization (or only computation, if length(lower)=0)
  if ( is.null(HLCor.args$adjMatrix) && is.null(attr(hlcor,"info.uniqueGeo")) ) { ## typically if DistMatrix was passed to HLCor...
    attr(hlcor,"info.uniqueGeo") <- uniqueGeo ## uniqueGeo should have been computed in all relevant cases where test is true (tricky)
  }
  attr(hlcor,"objective") <- anyHLCor_obj_args$`HLCor.obj.value` 
  attr(hlcor,"corrHLfitcall") <- mc
  ## 
  attr(hlcor,"optimInfo") <- list(optim.pars=optPars, 
                                  init.optim=init.optim,
                                  lower=lower,upper=upper,
                                  user.lower=user.lower,user.upper=user.upper)
  if ( ! is.null(locoptr <- attr(optPars,"optr")) && locoptr$convergence>0L) {
    hlcor$warnings$optimMessage <- paste("optim() message: ",locoptr$message," (convergence=",locoptr$convergence,")",sep="")
  }
  if ( ( ! is.null(optPars)) && attr(optPars,"method")== "optimthroughSmooth") {
    # provide logL estimate from the smoothing, to be used rather than the hlcor logL :
    logLapp <- optr$value
    attr(logLapp,"method") <- "  logL (smoothed)" 
    if (FALSE) { ## essai correction biais
      message("Trying to correct bias of SEM procedure...")
      p_lambda <- NROW(hlcor$lambda.object$coefficients_lambda)
      np <- NCOL(hlcor$`X.pv`) + p_lambda
      lmframe <- replicate(4*((np+1)*(np+2)/2), {## ((np+1)*(np+2)/2) is number of params to be fitted
        rephlcor <- do.call("HLCor",HLCor.args)
        if (p_lambda==0L) {
          c(fixef(rephlcor), logLik(rephlcor)[1])
        } else if (p_lambda==1L) {
          c(fixef(rephlcor), lambda=rephlcor$lambda, logLik(rephlcor)[1])
        } else stop("code missing for p_lambda>0 for bias correction of SEM.")
        ## FR->FR + problem if one of the predict vaiables in named logLapp
      } )
      lmframe <-as.data.frame(t(lmframe))
      colnames(lmframe)[colnames(lmframe)=="(Intercept)"] <- "fixef.intercept" ## because "(Intercept)" does not work below
      quadfit <- eval(parse(text=paste("lm(logLapp ~ polym(",
                                       paste(colnames(lmframe)[1:np],collapse=","),
                                       ", degree=2,raw=TRUE),data=lmframe)"))) ## predict does not work if raw=FALSE
      lmlower <- apply(lmframe[,1:np,drop=FALSE],2,min)
      lmupper <- apply(lmframe[,1:np,drop=FALSE],2,max)
      parnames <- names(lmlower)
      ## R issue: [.data.frame loses col name for single col.
      ## => lmframe[which.max(lmframe$logLapp),1:np] is 1-row data.frame is np>1; is *unnamed* numeric if np=1 unless drop=FALSE; 
      ##    where we would wish it to be named numeric in both cases
      initpar <- unlist(lmframe[which.max(lmframe$logLapp),1:np,drop=FALSE]) 
      optpolym <- optim(par=initpar, 
                        fn=function(v) {
                          v <- data.frame(matrix(v,nrow=1)) ## matrix() loses names... ## data.frame uses names as col or row names dependeing on length(v)...
                          colnames(v) <- parnames ## hence give names as a last step; again drop=F will be needed for 1-col 
                          predict(quadfit, newdata=v[c(1,1),,drop=FALSE])[1] ## [c(1,1),] is an awful patch for predict.poly
                        }, 
                        lower=lmlower,upper=lmupper,method="L-BFGS-B",
                        control=list(fnscale=-1,parscale=lmupper-lmlower)) 
      ## computation of dfs as in compare.model.structures
      attr(logLapp,"optpolym") <- optpolym$value
      attr(logLapp,"lmframe") <- lmframe
      attr(logLapp,"replicates") <- replicate(20,logLik(do.call("HLCor",HLCor.args)))
      attr(logLapp,"stddevres") <- ((predict(optr$Krigobj)-optr$Krigobj$data$logLobj)/(1-optr$Krigobj$lev_phi))
      attr(logLapp,"df_replicates") <- np
    }
    hlcor$APHLs$logLapp <- logLapp
  }
  if (is.character(trace$file)) {
    ## crude display of variable names in the trace file
    traceNames <- paste("# ",paste(names(hlcor$APHLs),collapse=" "))
    traceNames <- paste(traceNames,"lambda",sep=" ")
    if ( ! is.null(hlcor$phi)) traceNames <- paste(traceNames,"phi",sep=" ")
    traceNames <- paste(traceNames,paste(names(anyHLCor_obj_args$skeleton),collapse=" "),sep=" ")
    traceNames <- paste(traceNames," and optim parameters in canonical scale  ",sep=" ")
    write(traceNames,file=trace$file,append=T)   
  }
  ####  rm(list=c(tmpName),pos=".GlobalEnv") ## removes HLtmp0 at the global level
  return(hlcor) ## it's the call which says it was returned by corrHLfit
}

    
