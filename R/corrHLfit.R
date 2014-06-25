corrHLfit <-
function(formula,data, ## matches minimal call of HLfit
                      init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      trace=list(file=NULL,append=T),
                      objective="p_bv", ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                      rho.mapping,
                      control.corrHLfit=list(), ## optim.scale, Optimizer, optimizer.args, maxIter, corners
                      processed=NULL, ## added 2014/02 for programming purposes
                      family=gaussian(),
                      ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) { 
  init.optim <- init.corrHLfit
  alternating <- control.corrHLfit$alternating 
  if (is.null(alternating)) alternating <- FALSE
  optim.scale <- control.corrHLfit$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  Optimizer <- control.corrHLfit$Optimizer 
  if (is.null(Optimizer)) Optimizer="L-BFGS-B" ## either "nlminb" or one of the methods of optim()
  optimizers.args <- control.corrHLfit[c("nlminb","optim","optimize")] 
  corners <- control.corrHLfit$corners 
  if (is.null(corners)) corners <- TRUE 
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
  if ( (! is.null(ranFix$ARphi)) && (! is.null(init.corrHLfit$ARphi)) ) {
    stop("(!) 'ARphi' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
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
  spatial.terms <- findSpatial(predictor)
  spatial.model <- spatial.terms[[1]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1]]) 
  } else {
    # corr.model <- "Matern" ## up to v 1.0; for the defective syntax of the scripts for the Ecography paper
    stop("spatial correlation model not specified in 'formula': was valid in version 1.0 but not later.")
  }
  if (!is.null(HLmethod <- dotlist$HLmethod)) {
      if (HLmethod=="SEM") dotlist$`try.chol` <- FALSE
  } ## FR->FR ??? why not handled by preprocess ? actuellement longue construction de HLCor.args puis copie dans anyOptim.args apres preprocess
  HLCor.formals <- names(formals(HLCor))
  HLfit.formals <- names(formals(HLfit))
  designL.formals <- names(formals(designL.from.Corr))
  makescaled.formals <- names(formals(make.scaled.dist))
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
  if (corr.model=="adjacency") {
    if ( is.null(HLCor.args$adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
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
      validrows <- validRows(formula=formula,resid.formula=dotlist$resid.formula,data=dt) ## will remove rows with NA's in required variables
      dt[validrows,,drop=FALSE] ## ## before Predictor is called and an LMatrix is added, etc.
    })
  } else {
    validrows <- validRows(formula=formula,resid.formula=dotlist$resid.formula,data=data) ## will remove rows with NA's in required variables
    data <- data[validrows,,drop=FALSE] ## 
  }
  HLCor.args$family <- family
  HLCor.args$data <- data
  ## fills init.optim with all necessary values. There must be values for all parameters that are to be optimized 
  init <- list() ## will keep the initial values in untransformed scale
  if ( ! is.null(HLCor.args$adjMatrix) ) { ## NEIGHBOR MODEL there is a explicit adjMatrix provided but then users must provide bounds for rho for non-euclidian models...
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
  } else { ## if (corr.model %in% c("Matern","corMatern","AR1"))
    geoMats <- makeCheckGeoMatrices(data,dotlist$distMatrix,HLCor.args$uniqueGeo,coordinates)
    nbUnique <- geoMats$nbUnique  
    distMatrix <- geoMats$distMatrix   
    uniqueGeo <- geoMats$uniqueGeo   
    if (rho.size<2) { ## can be 0 if no explicit rho in the input  
      HLCor.args$distMatrix <- distMatrix   
    } else {
      HLCor.args$uniqueGeo <- uniqueGeo 
    }
    if (corr.model=="AR1") {
      if (is.null(ranFix$ARphi) && (! is.numeric(init.HLfit$ARphi))) { 
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
      if (rho.size==0) {
        if ( ! missing(rho.mapping)) {rho.size <- length(unique(rho.mapping))}
      }
      if (rho.size<2) { ## can be 0 if no explicit rho in the input  
        if (inherits(distMatrix,"list")) {
          maxrange <- max(unlist(lapply(distMatrix,function(dd) max(c(-Inf,dd))))) ## les Inf to handle dist(0)...
                     -min(unlist(lapply(distMatrix,function(dd) min(c( Inf,dd)))))
        } else maxrange <- max(distMatrix)-min(distMatrix)        
      } else { ## rho.size >1
        if (inherits(uniqueGeo,"list")) {
          if (missing(rho.mapping)) {
            if (ncol(uniqueGeo[[1]])==rho.size) { ## ncol>1, rho of length ncol
              rho.mapping <- seq_len(rho.size)           
            } else stop("'rho.mapping' missing with no obvious default from the other arguments.")
          }
          maxrange <- lapply(unique(rho.mapping), function(idx) {
            ranges <- matrix(unlist(lapply(uniqueGeo,function(uu){
              if (nrow(uu)>1) {
                range(proxy::dist(uu[,rho.mapping==idx]))
              } else c(Inf,-Inf) ## encore des Inf to handle dist(0)...
              })),ncol=2)
            max(ranges[,2])-min(ranges[,1]) 
          })
        } else { ## single data set
          if (missing(rho.mapping)) {
            if (ncol(uniqueGeo)==rho.size) {
              rho.mapping <- seq_len(rho.size)           
            } else stop("'rho.mapping' missing with no obvious default from the other arguments.")
          }
          maxrange <- lapply(unique(rho.mapping), function(idx) {
            rng <- range(proxy::dist(uniqueGeo[,rho.mapping==idx]))
            rng[2]-rng[1] 
          })
        }  
        maxrange <- unlist(maxrange)
        HLCor.args$`rho.mapping` <- rho.mapping
      }
      if (is.null(ranFix$rho) && (! is.numeric(init.HLfit$rho))) {
        init$rho <- init.optim$rho 
        if (is.null(init$rho)) init$rho <- 30/(2*maxrange) 
        if (! is.null(init.HLfit$rho)) {
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
      if (is.null(ranFix$nu) && (! is.numeric(init.HLfit$nu))) { 
        init$nu <- init.optim$nu 
        if (is.null(init$nu)) init$nu <- 0.5 
        if (! is.null(init.HLfit$nu)) {
          init.HLfit$nu <- init$nu ## avant transformation
        } else {
          if (optim.scale=="transformed") {
            if (is.null(ranFix$rho)) { 
              init.optim$trNu <- nuFn(init$nu,init$rho) 
            } else init.optim$trNu <- nuFn(init$nu,ranFix$rho)
            init.optim$nu <- NULL
          } else init.optim$nu <- init$nu
        }
      } 
      if ( (! is.null(init.HLfit$nu)) && (! is.numeric(init.HLfit$nu))) {
        init.HLfit$nu <- init.optim$nu
        init.optim$nu <- NULL
      }
    }
    if (is.null(ranFix$Nugget)) { init$Nugget <- init.optim$Nugget }  ## this may be null, but in this case we leave it so and the Nugget keeps its default value through all computations
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
  lower <- init.optim; upper <- init.optim ## initialization with the right variables but wrong values
  if (corr.model=="adjacency") { ## adjacency model
    ## no default value, user values are required 
    lower$rho <- user.lower$rho ## no transfo for adjacency model
    upper$rho <- user.upper$rho ## idem
  } else {
    if (corr.model=="AR1") {
      if ( ! is.null(init$ARphi)) {
        ARphi <- user.lower$ARphi
        if (is.null(ARphi)) ARphi <- -0.999999
        lower$ARphi <- ARphi
        ARphi <- user.upper$ARphi
        if (is.null(ARphi)) ARphi <- 0.999999
        upper$ARphi <- ARphi
      }    
    } else { ## then Matern model....
      if (! is.null(init$rho)) {
        rho <- user.lower$rho
        if (is.null(rho)) rho <- init$rho/150
        if (optim.scale=="transformed") {
          lower$trRho <- rhoFn(rho)
        } else lower$rho <- rho
        rho <- user.upper$rho
        if (is.null(rho)) {
          if (inherits(nbUnique,"list")) nbUnique <- mean(unlist(nbUnique))
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
    } 
    ##### common to the different models except adjacency (because there are several places where NUgget+adjacency is not handled)
    if ( ! is.null(init$Nugget)) {
      lower$Nugget <- 0
      upper$Nugget <- 0.999999
    }
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
  ################
  LowUp <- list(lower=lower,upper=upper) ## inherits names from init.optim, must be logscale as init.optim is by construction
  varNames <- names(init.HLfit)
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
    if (family$family %in% c("poisson","binomial")) {
      phi.Fix <- 1 
    } else {
      phi.Fix <- ranFix$phi
      if (any(phi.Fix==0)) {
        mess <- pastefrom("phi cannot be fixed to 0.",prefix="(!) From ")
        stop(mess)
      }
    } 
    preprocess.formal.args$phi.Fix <- phi.Fix 
    if (inherits(data,"list")) {
      processed <- lapply(data,function(dd) {
        locargs <- preprocess.formal.args
        locargs$data <- dd
        if ( ! is.null(famfam) && famfam=="multi") {
          locargs$family <- family$binfamily ## only local. multinomial family will be passed to HLCor.obj 
        }
        do.call(preprocess,locargs)}
      )
      for (st in names(processed[[1]])) HLCor.args[st] <- NULL ## this keeps the data as they are not returned in processed, but predictor is suppressed
    } else {
      processed <- do.call(preprocess,preprocess.formal.args)
      for (st in names(processed)) HLCor.args[st] <- NULL ## this keeps the data as they are not returned in processed, but predictor is suppressed
    }
  }
  HLCor.args$processed <- processed
  HLCor.args$ranPars <- ranPars
  anyOptim.args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  anyOptim.args$skeleton <- init.optim ## logscale, only used by HLCor.obj
  attr(anyOptim.args$skeleton,"type") <- list() 
  attr(anyOptim.args$skeleton,"type")[names(init.optim)] <- "fix" 
##################  anyOptim.args$fn <- HLCor.obj ##### that is the function optimized 
  anyOptim.args$HLCor.obj.value <- objective ## moved here 25/11/2012                    
  ## optim/optimize specific code
  initvec <- unlist(init.optim)
  ####    tmpName <- generateName("HLtmp") ## tmpName is a string such as "HLtmp0"
  #    anyOptim.args$init.HLfit <- tmpName 
  ####    assign(tmpName,list(),pos=".GlobalEnv") ## sets HLtmp0 (or a similarly named variable) at the global level
  processedHL1 <- processed$HL[1] ## don't assume an HLmethod argument here... although that should be valid
  if (!is.null(processedHL1) && processedHL1=="SEM" && length(lower)>0) {## FR->FR il faudra integrer length lower = 0 dans la function pour homogénéiser le code.
    pargrid <- sampleGridFromLowUp(LowUp,n=init.corrHLfit$nSmoothed) ## n may be NULL
    ## using PQL to find a good starting region
    PQLarglist <- list(pargrid=pargrid,anyHLCor.args=anyOptim.args)
    ## PQLarglist$anyHLCor.args$HLmethod <- "PQL/L" ## does not work since the HLmethod info has already been processed
    ###### 
    PQLarglist$anyHLCor.args$processed$HL <- c(0,0,1)
    # this will perform an RE PQL fit; it is difficult to have more control wihtout reprocessing the arguments, which implies that 
    # the call to corrHLfit ditd not involve preprocessed arguments...
    ###### 
    PQLoptr <- do.call(locoptimthroughSmooth,PQLarglist)
    Krigobj <- PQLoptr$Krigobj
    predVar <- as.numeric(attr(predict(Krigobj,newX=PQLoptr$par,predVar=TRUE),"predVar"))
    oldpredVar <- 0
    ## new sampling for SEM guided by the PQL results
    blocksize <- 30
    nextpoints <- sampleNextPoints(n=blocksize,optr=PQLoptr,minPtNbr=3,expand=1,D.resp=sqrt(predVar)/2) ## random sample
    info <- attr(nextpoints,"info") ## only used if diagnostic plot but removed by the rbind
    nearbypts <- sampleNearby(nextpoints,n=6,stepsizes=(unlist(upper)-unlist(lower))/100)      
    nextpoints <- rbind(nextpoints,nearbypts)
    ## diagnostic plot for previous and next computation
    if (interactive() && length(lower)==2) {
      zut <- signif(unlist(Krigobj$corrPars$rho),4)
      titlemain <- bquote(paste("PQL initialization: ",rho[f(rho)],"=",.(zut[1]),", ",rho[f(nu)],"=",.(zut[2]),"; predVar=",.(signif(predVar,4))))
      mapMM(Krigobj,plot.title=title(main=titlemain),
            decorations= {
              points(nextpoints,pch=15,cex=0.4);
              apply(info$simplicesTable,1,function(v){
                polygon(info$vertices[v,])
              });
              points(matrix(PQLoptr$par,nrow=1),pch="+",col="red");
            }
      ) 
    }
    ## now the SEM computations
    pargrid <- rbind(pargrid,nextpoints)
    arglist <- list(pargrid=pargrid,anyHLCor.args=anyOptim.args)
    arglist[["control.smooth"]]$nrepl <- 20 ## number of points for which replicate estimates of likelihood are computed (modified later)
    optr <- do.call(locoptimthroughSmooth,arglist)    
    smoothingOK <- FALSE
    control.smooth <- list()
    it <- 1
    continue <- TRUE
    while ( continue ) {
      oldpredVar <- predVar
      Krigobj <- optr$Krigobj
      predVar <- as.numeric(attr(predict(Krigobj,newX=optr$par,predVar=TRUE),"predVar"))
      if (interactive() ) {cat(it," ");print(paste("(",paste(signif(optr$value,4),collapse=","),")",signif(optr$value,4),
                                                   "+/-",signif(sqrt(predVar),4),
                                                   "; n_points=",nrow(Krigobj$data),
                                                   "; lambda=",signif(Krigobj$lambda,4),sep=""))} 
      prevPtls <- optr$forSmooth
      smoothRho <- unlist(Krigobj$corrPars$rho)
      ## tests whether some correlation structure has been detected and adjust smoothing controls accordingly:
      if ( predVar < oldpredVar && predVar < Krigobj$lambda /80 ){ # nrow(Krigobj$data) ) {
        smoothingOK <- TRUE 
        trysize <- blocksize <- 1
        expand <- 1
        control.smooth$ranFix <- Krigobj$corrPars
        control.smooth$nrepl <- 1 
      } else {
        smoothingOK <- FALSE
        trysize <- 90
        blocksize <- 30
        expand <- min(it,2)
        control.smooth$ranFix <- Krigobj$corrPars["nu"] 
        control.smooth$nrepl <- ceiling(20/it - 0.0001)
        ## get rid of some possibly aberrant points that prevent good smoothing 
        prevPtls <- prevPtls[order(prevPtls$logL)[-c(1:2)],] ## FR->FR but aberrant points may not be the lowest... 
      }
      trypoints <- sampleNextPoints(n=trysize,optr=optr,expand=expand,D.resp=sqrt(predVar)/2) ## random sample
      info <- attr(trypoints,"info") ## only used if diagnostic plot but removed by the rbind
      ###### selection of points
      obspred <- predict(Krigobj,predVar=TRUE)
      obsSE <- diag(attr(obspred,"predVar"))
      obsSE[obsSE<0] <- 0
      obsSE <- sqrt(obsSE)
      Qmax <- max(predict(Krigobj)$fitted+1.96 * obsSE)
      #
      trypred <- predict(Krigobj,trypoints,predVar=TRUE)
      trySE <- diag(attr(trypred,"predVar"))
      trySE[trySE<0] <- 0
      trySE <- sqrt(trySE)
      tryQ <- trypred$fitted + 1.96*trySE
      #
      expectedImprovement <- trySE*dnorm((Qmax-tryQ)/trySE)+(tryQ-Qmax)*pnorm((tryQ-Qmax)/trySE)
      trypoints <- cbind(trypoints,EI=expectedImprovement)
      trypoints <- trypoints[order(trypoints[,"EI"],decreasing=TRUE)[seq_len(blocksize)],,drop=FALSE]
      nextpoints <- trypoints[which(trypoints[,"EI"]>0),names(lower),drop=FALSE] ## maybe no point...
      ## 
      if(nrow(nextpoints)==0) nextpoints <- rbind(nextpoints,optr$par)
      if( ! smoothingOK) { ## need close pairs to estimate better the smoothing parameters
        ## FR->FR problem the following points may extrapolate... particularly for small smoothRho
        nearbypts <- sampleNearby(nextpoints,n=min(nrow(nextpoints),6),stepsizes=(unlist(upper)-unlist(lower))/(100*smoothRho))     
        nextpoints <- rbind(nextpoints,nearbypts)
      }
      ## and a bit of extrapolation
      #       if (it>1) {
      #         cS <- connectedSets(info$simplicesTable)
      #         outerpoints <- lapply(cS, function(v){
      #           v <- intersect(v,info$innerVertexIndices) ## only the really good points in the set
      #           pts <- info$vertices[v,,drop=FALSE]
      #           if (nrow(pts)>length(lower)+1) { ## more vertices than a simplex => can be redundant
      #             return(pts[unique(as.vector(convhulln(info$vertices[v,],"Pp"))),])
      #           } else return(pts) ## extrapolhull will handle special cases
      #         })
      #         extrap <- lapply(outerpoints,extrapolhull)
      #         extrap <- do.call(rbind,extrap)
      #         nextpoints <- rbind(nextpoints,extrap)
      #       }
      ##
      if (interactive() && length(lower)==2 && 
            ! smoothingOK) {
        zut <- signif(smoothRho,4)
        titlemain <- bquote(paste(.(DEPARSE(processed$predictor)),", iter=",.(it)))
        if (nchar(eval(titlemain))>50) {
          titlemain <- bquote(paste(.(DEPARSE(nobarsNooffset(processed$predictor))),"+..., iter=",.(it)))
        }
        if (nchar(eval(titlemain))>50) {
          titlemain <- bquote(paste(.(substr(aschar,0,35)),"+... [length(",beta,")=",.(ncol(processed$X.pv)),"], iter=",.(it)))
        }
        titlesub <- bquote(paste(rho[f(rho)],"=",.(zut[1]),", ",rho[f(nu)],"=",.(zut[2]),
                                  "; max=",.(signif(optr$value,4)),"; predVar=",.(signif(predVar,4))))
        mapMM(Krigobj,plot.title={
          ## inconsistent behaviour: one can title(main=<bquote stuff> ...) 
          ## but not title(main=<eval(bquote ...) stuff>) => error : cannot evaluate f in f(rho), f(nu)
          title(main=titlemain,line=2) ## default cex.main=1.2, line ~1.7
          title(main=titlesub,line=0.8,cex.main=1.1) ## 
        }, 
              decorations= {
                points(nextpoints,pch=15,cex=0.4);
                apply(info$simplicesTable,1,function(v){
                  polygon(info$vertices[v,])
                });
                points(matrix(optr$par,nrow=1),pch="+",col="red");
              }
        ) 
      }
      #browser()
      arglist <- list(pargrid=nextpoints,control.smooth=control.smooth,anyHLCor.args=anyOptim.args,prevPtls=prevPtls)
      optr <- do.call(locoptimthroughSmooth,arglist)     
      it <- it+1
      ## terminates if either of these two considtions are reached *...* :
      if (predVar < 0.02) continue <- FALSE
      if (it > 50) continue <- FALSE
      ## ... UNLESS one of these conditions are true
      if (it < 10) continue <- TRUE
      if (predVar> oldpredVar) continue <- TRUE
    } ## end 'while' loop
    Krigobj <- optr$Krigobj
    predVar <- as.numeric(attr(predict(Krigobj,newX=optr$par,predVar=TRUE),"predVar"))
    if (interactive() && length(lower)==2) {
      zut <- signif(unlist(Krigobj$corrPars$rho),4)
      titlemain <- bquote(paste(rho[f(rho)],"=",.(zut[1]),", ",rho[f(nu)],"=",.(zut[2]),
                                "; max=",.(signif(optr$value,4)),"; predVar=",.(signif(predVar,4))))
      mapMM(Krigobj,plot.title=title(main=titlemain),
            decorations= {
              points(matrix(optr$par,nrow=1),pch="+");
            }
      ) 
    }
    optPars <- as.list(optr$par)
    attr(optPars,"method") <-"locoptimthroughSmooth"
  } else if (alternating) { ## renewed coding of the iterative algo (only p_v); not documented => not checked for a long time
    optPars <- alternating(init.optim=init.optim,LowUp=LowUp,anyOptim.args=anyOptim.args,maxIter=maxIter,
                           ranPars=ranPars,HLCor.args=HLCor.args,trace=trace,Optimizer=Optimizer,
                           optimizers.args=optimizers.args,corners=corners)
    if (!is.null(optPars)) attr(optPars,"method") <-"alternating"
  } else {
    optPars <- locoptim(init.optim,LowUp,anyOptim.args=anyOptim.args,trace,Optimizer=Optimizer,optimizers.args=optimizers.args,corners=corners) 
    if (!is.null(optPars)) attr(optPars,"method") <-"locoptim"
  }
  ranPars[names(optPars)] <- optPars
  attr(ranPars,"type")[names(optPars)] <- "fix" ## fix dans HLfit !
  HLCor.args$ranPars <- ranPars
  verbose["warn"] <- TRUE ## important!
  HLCor.args$verbose <- verbose ##
  hlcor <- do.call(HLCor,HLCor.args) ## recomputation post optimization
  if ( is.null(HLCor.args$adjMatrix) && is.null(attr(hlcor,"info.uniqueGeo")) ) { ## typically if DistMatrix was passed to HLCor...
    attr(hlcor,"info.uniqueGeo") <- uniqueGeo ## uniqueGeo should have been computed in all relevant cases where this is NULL (tricky)
  }
  attr(hlcor,"objective") <- objective
  attr(hlcor,"ranFixNames") <- names(ranFix)
  if ( ( ! is.null(optPars)) && attr(optPars,"method")== "locoptimthroughSmooth") {
    attr(optr$value,"predVar") <- predVar
    APHLs$estlogL <- optr$value
  }
  if (is.character(trace$file)) {
    ## crude display of variable names in the trace file
    traceNames <- paste(paste(names(hlcor$APHLs),collapse=" "))
    traceNames <- paste(traceNames,"lambda",sep=" ")
    if ( ! is.null(hlcor$phi)) traceNames <- paste(traceNames,"phi",sep=" ")
    traceNames <- paste(traceNames,paste(names(anyOptim.args$skeleton),collapse=" "),sep=" ")
    write(traceNames,file=trace$file,append=T)   
  }
  #    ranPars <- HLCor.args$ranPars
  #    class(ranPars) <- c("ranPars",class(ranPars)) ## ready for a print.ranPars method 
  #  }
  ####  rm(list=c(tmpName),pos=".GlobalEnv") ## removes HLtmp0 at the global level
  return(hlcor)
}
