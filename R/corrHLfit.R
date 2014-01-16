corrHLfit <-
function(formula,data, ## matches minimal call of HLfit
                      init.corrHLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      trace=list(file=NULL,append=T),
                      objective="p_bv",
                      control.corrHLfit=list(), ## optim.scale, optim.method, optim.args, maxIter
                      ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) { 
  init.optim <- init.corrHLfit
  # verbose <- control.corrHLfit$verbose 
  alternating <- control.corrHLfit$alternating 
  if (is.null(alternating)) alternating <- FALSE
  optim.scale <- control.corrHLfit$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  optim.method <- control.corrHLfit$optim.method 
  if (is.null(optim.method)) optim.method="L-BFGS-B" ## either "nlminb" or one of the methods of optim()
  optim.args <- control.corrHLfit$optim.args 
  if (is.null(optim.args)) optim.args=list() ##works for optim() only
  maxIter <- control.corrHLfit$maxIter 
  if (is.null(maxIter)) maxIter<- 10000
  if ( ! (objective %in% c("p_v","p_bv"))) {
    mess <- pastefrom("invalid value of the 'objective' argument.",prefix="(!) From ")
    stop(mess)
  }
  if ( (! is.null(ranFix$rho)) && (! is.null(init.corrHLfit$rho)) ) {
    stop("(!) 'rho' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  } else rho.size <- max(length(ranFix$rho),length(init.corrHLfit$rho))
  if ( (! is.null(ranFix$nu)) && (! is.null(init.corrHLfit$nu)) ) {
    stop("(!) 'nu' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  }
  if ( (! is.null(ranFix$Nugget)) && (! is.null(init.corrHLfit$Nugget)) ) {
    stop("(!) 'Nugget' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  }
  if (length(trace)>0 && ! all(names(trace) %in% c("file","append"))) {
    mess <- pastefrom("'trace' elements should be named and = 'file','append'.",prefix="(!) From ")
    message(mess)
  } 
  if (is.character(trace$file) && ( ! trace$append) ) {
    try(unlink(trace$file))
  }  ## the file is written in by HLCor()                   
  datanames <- names(data)
  ## FR->FR ? ici essayer une autre syntaxe avec match.call(expand.dots=T)
  dotlist <-list(...) 
  ## Preventing obsolete options
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
  if (! "predictor" %in% class(predictor)) predictor <- Predictor(formula) 
  spatial.model <- findSpatial(predictor)[[1]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1]]) 
  } else {
    # corr.model <- "Matern" ## up to v 1.0; for the defective syntax of the scripts for the Ecography paper
    stop("spatial correlation model not specified in 'formula': was valid in version 1.0 but not later.")
  }
  HLCor.formals <- names(formals(HLCor))
  HLfit.formals <- names(formals(HLfit))
  designL.formals <- names(formals(designL.from.Corr))
  HLnames <- (c(HLCor.formals,HLfit.formals,designL.formals))  ## cf parallel code in HLCor.obj
  argcheck <- names(dotlist)[which(! names(dotlist) %in% HLnames)]
  if (length(argcheck)>0) warning(paste("suspect arguments(s) ",paste(argcheck, collapse=",")," in corrHLfit call."))
  HLCor.args <- dotlist[intersect(names(dotlist),HLnames)] ## those specifically for HLCor and those for HLfit and designL.from.Corr
  HLCor.args$ranPars <- ranFix
  HLCor.args$ranFix <- NULL ## make sure it's no longer there
  HLCor.args$data <- data
  HLCor.args$trace <- trace$file
  HLCor.args$formula <- predictor
  HLCor.args$corr.model <- corr.model
  verbose <- dotlist$verbose 
  ## we want the default to be FFF in HLCOr, except for the final fit 
  ## a coherent user could use TFF to display all HLCor results, or worse...
  lv <- length(verbose)
  if (lv<3) verbose<- c(verbose,c(FALSE,FALSE,FALSE)[(lv+1):3]) ## ensures that the middle one is F, whch is not default for HLfit
  if (lv==3) verbose <- c(TRUE,verbose) ## encodes a corrHLfit default for distinct HLCor calls...
  HLCor.args$verbose <- verbose[-1] ## default $verbose should here be F F F from verbose = T F F F
  ## one more modif of HLCor.args : distMatrix or uniqueGeo below
  #  rhoobj.formals <- names(formals(rho.obj)) ## list with default values !
  #  rhoobj.args <- dotlist[intersect(names(dotlist),rhoobj.formals)]
  ############### (almost) always check geo info ###################
  if (corr.model=="adjacency") {
    if ( is.null(HLCor.args$adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
  } else {
    if ( ! is.null(spatial.model)) {
      bars <- spatial.model[[2]] 
      coordinates <- deparse(bars[[3]]) ## "x + y"
      coordinates <-  strsplit(coordinates," ")[[1]]
      coordinates <- coordinates[coordinates != "+"]
    } else {
      stop("An obsolete syntax for the adjacency model appears to be used.")
      ## coordinates <- c("x","y") ## backward compatibility... old syntax with (1|pos) and default values of the coordinates argument
    }
    coordcols <- which(datanames %in% coordinates)
    if ( length(coordcols) != length(coordinates) ) {
      stop("variables 'coordinates' not all found in the 'data'")
    }
  } ## 
  #########################################################
  ## fills init.optim with all necessary values. There must be values for all parameters that are to be optimized 
  init <- list() ## will keep the initial values in untransformed scale
  if ( is.null(HLCor.args$adjMatrix) ) { ## then MATERN MODEL #  if (corr.model %in% c("Matern","corMatern"))
    distMatrix <- dotlist$distMatrix
    if (is.null(dotlist$distMatrix)) { ## if there is no distMatrix then 
      ## (1) we need to construct auxiliary matrices for HLCor/HLfit
      ## (2) we need distMatrix *here* to determine good search values for rho
      ## (1):
      if ( is.null(HLCor.args$uniqueGeo) ) { ## then construct it from the data ## this should be the routine case
        uniqueGeo <- unique(data[,coordinates,drop=F])
      }
      nbUnique <- nrow(uniqueGeo) 
      ## (2): we need distMatrix *here* in all cases
      ##distMatrix <- as.matrix(dist(uniqueGeo))
      distMatrix <- dist(uniqueGeo)
      if (rho.size<2) { ## can be 0 if no explicit rho in the input  
        HLCor.args$distMatrix <- distMatrix   
      } else {
        HLCor.args$uniqueGeo <- uniqueGeo 
      }
    } else { ## there is a distMatrix, this is what will be used by HLCor
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
      } else {
        uniqueGeo <- unique(data[usernames,coordinates,drop=F]) ## check that this corresponds to unique locations
        nbUnique <- nrow(uniqueGeo)
        if (nbUnique != nrow(distMatrix)) {
          stop("The dimension of 'distMatrix' does not match the number of levels of the grouping variable")
        } else { ## check order
          redondGeo <- data[,coordinates,drop=F]
          designRU <- apply(redondGeo,1,function(v) {which(apply(v==t(uniqueGeo),2,all))}) ## has no names
          ## eg 1 1 2 2 3 2 3 4 is valid for 8 obs, 4 unique locations
          designRU <- unique(as.vector(designRU)) ## should then be 1 2 3 4
          ## but if distMatrix in reverse order, the first row of redondGeo would match the 4th of uniqueGeo and then the following test is FALSE:
          if ( ! all (designRU==seq(length(designRU))) ) {
            stop("The rows of 'distMatrix' are not ordered as rows of the 'data'.")
          }
        } 
      }
      HLCor.args$distMatrix <- distMatrix   ## passed the checks
    }
    maxrange<-max(distMatrix)-min(distMatrix)
    if (is.null(ranFix$rho)) {
      init$rho <- init.optim$rho 
      if (is.null(init$rho)) init$rho <- 30/(2*maxrange) 
      if (optim.scale=="transformed") {
        init.optim$trRho <- rhoFn(init$rho) ## we're in Matern model here
        init.optim$rho <- NULL
      } else init.optim$rho <- init$rho
    } 
    if (is.null(ranFix$nu)) { 
      init$nu <- init.optim$nu 
      if (is.null(init$nu)) init$nu <- 0.5 
      if (optim.scale=="transformed") {
        if (is.null(ranFix$rho)) { 
          init.optim$trNu <- nuFn(init$nu,init$rho) 
        } else init.optim$trNu <- nuFn(init$nu,ranFix$rho)
        init.optim$nu <- NULL
      } else init.optim$nu <- init$nu
    } 
    if (is.null(ranFix$Nugget)) { init$Nugget <- init.optim$Nugget }  ## this may be null, but in this case we leave it so and the Nugget keeps its default value through all computations
  } else { ## NEIGHBOR MODEL there is a explicit adjMatrix provided but then users must provide bounds for rho for non-euclidian models...
    if (is.null(ranFix$rho)) {
      if (is.null(lower$rho)) {
        mess <- pastefrom("lower$rho required.",prefix="(!) From ")
        stop(mess)
      }
      if (is.null(upper$rho)) {
        mess <- pastefrom("upper$rho required.",prefix="(!) From ")
        stop(mess)
      }
      init$rho <- init.optim$rho 
      if (is.null(init$rho)) init$rho <- (lower$rho+upper$rho)/2 
      init.optim$rho <- init$rho
    }
  }
  if (is.null(ranFix$lambda)) { ## no ranFix$lambda: process init.optim
    init$lambda <- init.optim$lambda 
    if (!is.null(init$lambda)) {
      if (init$lambda<1e-4) init$lambda <- 1e-4
      init.optim$trLambda <- dispFn(init$lambda) 
      init.optim$lambda <- NULL
    }
  } else { ## ranFix$lambda present, do NOT put it in init.optim
    if (!is.null(init.optim$lambda)) stop("(!) Arguments 'ranFix$lambda' and 'init.corrHLfit$lambda' conflict with each other.")  
  } 
  if (is.null(ranFix$phi)) {
    init$phi <- init.optim$phi 
    if (!is.null(init$phi)) {
      if (init$phi<1e-4) init$phi <- 1e-4
      init.optim$trPhi <- dispFn(init$phi)
      init.optim$phi <- NULL
    }
  } else {
    if (!is.null(init.optim$phi)) stop("(!) Arguments 'ranFix$phi' and 'init.corrHLfit$phi' conflict with each other.")  
  } 
  ## done with init.optim
  # { ## ad hoc maximization of p_bv VIA HLCor.obj
  ## by maxim over (corrpars,phi,lambda, (beta)_[corrpars,phi,lambda])
  ##     if trPhi,trLambda are in the init.optims
  ## or by maxim over (corrpars,(beta,phi,lambda)_corrpars)
  ##     otherwise.
  ################ construct intervals for this maximization
  ## construct default upper and lower values ; on transformed scale by default
  user.lower <- lower; user.upper <- upper ## keep user input 
  lower <- init.optim; upper <- init.optim ## initialization with the right variables but wrong values
  if (is.null(HLCor.args$adjMatrix) ) { ## then Matern model....
    if (! is.null(init$rho)) {
      rho <- user.lower$rho
      if (is.null(rho)) rho <- init$rho/150
      if (optim.scale=="transformed") {
        lower$trRho <- rhoFn(rho)
      } else lower$rho <- rho
      rho <- user.upper$rho
      if (is.null(rho)) {
        rho <- init$rho*2*nbUnique ## The following was a bit too low for experiments with nu=0.5 : 1/(maxrange/(2*nbUnique)) ## nb => unique rows !
        if (optim.scale=="transformed") rho <- rho*.SpaMM$RHOMAX/(1+rho) ## so that it does not exceed RHOMAX
      }
      if (optim.scale=="transformed") {
        upper$trRho <- rhoFn(rho) 
      } else upper$rho <- rho
      rhoForNu <- init$rho
    } else rhoForNu <-ranFix$rho
    if (! is.null(init$nu)) {
      nu <- user.lower$nu
      if (is.null(nu)) nu <- init$nu/100
      if (optim.scale=="transformed") {
        lower$trNu <- nuFn(nu,rhoForNu)
        #print(c(rhoForNu,nu,lower$trNu))
      } else lower$nu <-nu
      nu <- user.upper$nu
      if (is.null(nu)) nu <- init$nu*.SpaMM$NUMAX/(1+init$nu) ## nu should not diverge otherwise it will diverge in Bessel_lnKnu, whatever the transformation used
      if (optim.scale=="transformed") {
        upper$trNu <- nuFn(nu,rhoForNu)
      } else upper$nu <- nu
      #print(c(rhoForNu,nu,upper$trNu))
    }
  } else { ## adjacency model
    ## no default value, user values are required 
    lower$rho <- user.lower$rho ## no transfo for adjacency model
    upper$rho <- user.upper$rho ## idem
  }
  if ( ! is.null(init$Nugget)) {
    lower$Nugget <- 0
    upper$Nugget <- 0.999999
  }
  if (! is.null(init$phi)) {
    phi <- user.lower$phi
    if (is.null(phi)) phi <- init$phi/1000
    lower$trPhi <- dispFn(phi)
    phi <- user.upper$phi
    if (is.null(phi)) phi <- init$phi*1000
    ## if phi is badly initialized then it gets a default which may cause hard to catch problems in the bootstrap...
    upper$trPhi <- dispFn(phi)
  }
  if (! is.null(init$lambda)) {
    lambda <- user.lower$lambda
    if (is.null(lambda)) lambda <- init$lambda/1000
    lower$trLambda <- dispFn(lambda)
    lambda <- user.upper$lambda
    if (is.null(lambda)) lambda <- init$lambda*1000
    upper$trLambda <- dispFn(lambda)
  }
  LowUp <- list(lower=lower,upper=upper) ## inherits names from init.optim, must be logscale as init.optim is by construction
  ################
  ########## common stuff to both optim and optimize
  HLCor.args$ranPars <- ranFix 
  ####
  HLfit.formal.args <- formals(HLfit) ## makes sure about default values 
  nondefault.HLfit.formal.args <- HLCor.args[which(names(HLCor.args) %in% HLfit.formals)]
  HLfit.formal.args[names(nondefault.HLfit.formal.args)] <- nondefault.HLfit.formal.args ## full HLfit args
  preprocess.formal.args <- HLfit.formal.args[which(HLfit.formals %in% names(formals(preprocess)))] ## only those for preprocess
  preprocess.formal.args$predictor <- HLfit.formal.args$formula ## because preprocess stll expects $predictor 
  preprocess.formal.args$resid.predictor <- HLfit.formal.args$resid.formula ## because preprocess stll expects $predictor 
  processed <- do.call(preprocess,preprocess.formal.args)
  for (st in names(processed)) HLCor.args[st] <- NULL ## this keeps the data as they are not returned in processed, but predictor is suppressed
  HLCor.args$processed <- processed
  anyOptim.args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  anyOptim.args$skeleton <- init.optim ## logscale, only used by HLCor.obj
  anyOptim.args$f <- HLCor.obj ##### that is the function optimized
  anyOptim.args$HLCor.obj.value <- objective ## moved here 25/11/2012                    
  ## optim/optimize specific code
  initvec <- unlist(init.optim)
  ####    tmpName <- generateName("HLtmp") ## tmpName is a string such as "HLtmp0"
  #    anyOptim.args$init.HLfit <- tmpName 
  ####    assign(tmpName,list(),pos=".GlobalEnv") ## sets HLtmp0 (or a similarly named variable) at the global level 
  if (alternating) { ## renewed coding of the iterative algo (only p_v); 
    nam <- names(init.optim)
    if (any(c("trPhi","trLambda") %in% nam )) {
      mess <- pastefrom("Dispersion parameters non allowed in 'init.corrHLfit' with alternating algorithm.",prefix="(!) From ")
      stop(mess)
    }
    initHLfit <- init.optim[! nam %in% c("trRho","trNu","Nugget")] ## vector, not list
    initcorr <- init.optim[nam %in% c("trRho","trNu","Nugget")]
    HLfitLowUp <- LowUp
    HLfitLowUp$lower[c("trRho","trNu","Nugget")] <- NULL
    HLfitLowUp$upper[c("trRho","trNu","Nugget")] <- NULL
    corrLowUp <- LowUp
    corrLowUp$lower[c("trPhi","trLambda")] <- NULL
    corrLowUp$upper[c("trPhi","trLambda")] <- NULL
    anycorrOptim.args <- anyOptim.args
    iter <- 0
    conv <- 1
    currentLik <- -Inf
    while (iter < maxIter && conv > 1e-5 ) { ## alternate HLCor and locoptim
      HLCor.args$ranPars[names(initcorr)] <- initcorr
      oldLik <- currentLik
      if (is.character(trace$file)) {
        if(.SpaMM$TRACE.UNLINK) unlink("HLCor.args.*.RData")
        zut <- paste(unlist(initcorr),collapse="")  
        save(HLCor.args,file=paste("HLCor.args.",zut,".RData",sep="")) ## for replicating the problem
      }
      givencorr <- do.call(HLCor,HLCor.args)$hlfit ## optim disp and beta given corr param
      currentLik <- givencorr$APHLs$p_v ## iterations maximize p_v
      conv <- currentLik-oldLik
      anycorrOptim.args$ranPars$lambda <- givencorr$lambda
      anycorrOptim.args$ranPars$phi <- givencorr$phi
      #### anycorrOptim.args$etaFix <- list(beta=givencorr$fixef,v_h=givencorr$v_h) ## that's what LeeN01sm say, but this does not work
      anycorrOptim.args$etaFix <- list(beta=givencorr$fixef) 
      initcorr <- locoptim(initcorr,corrLowUp,anycorrOptim.args,trace,optim.method,optim.args) 
      iter <- iter+1
    }
    optPars <- c(initcorr,givencorr$lambda,givencorr$phi)
  } else {
    optPars <- locoptim(init.optim,LowUp,anyOptim.args,trace,optim.method,optim.args) 
  }
  HLCor.args$ranPars[names(optPars)] <- optPars
  HLCor.args$verbose[1] <- verbose[1] ## default $verbose should here be T F F from verbose = T F F F
  hlcor <- do.call(HLCor,HLCor.args) ## recomputation post optimization
  hlfit <- hlcor$hlfit
  resu <- list(hlfit=hlfit,HLCor.info=hlcor$HLCor.info,objective=objective,ranFixNames=names(ranFix))
  if (is.character(trace$file)) {
    ## crude display of variable names in the trace file
    traceNames <- paste(paste(names(hlfit$APHLs),collapse=" "))
    traceNames <- paste(traceNames,"lambda",sep=" ")
    if ( ! is.null(hlfit$phi)) traceNames <- paste(traceNames,"phi",sep=" ")
    traceNames <- paste(traceNames,paste(names(anyOptim.args$skeleton),collapse=" "),sep=" ")
    write(traceNames,file=trace$file,append=T)   
  }
  #    ranPars <- HLCor.args$ranPars
  #    class(ranPars) <- c("ranPars",class(ranPars)) ## ready for a print.ranPars method 
  #  }
  ####  rm(list=c(tmpName),pos=".GlobalEnv") ## removes HLtmp0 at the global level
  class(resu) <- c("corrHLfit",class(hlfit)) 
  return(resu)
}
