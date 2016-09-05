fitme <- function(formula,data, ## matches minimal call of HLfit
                  family=gaussian(),
                  init=list(),
                  fixed=list(), ## replaces ranFix
                  lower=list(),upper=list(),
                  resid.model=~1,
                  init.HLfit=list(),
                  control=list(), ## optim.scale (private), nloptr, refit
                  control.dist=list(),
                  method="ML", 
                  HLmethod=method, ## LRT fns assume HLmethod when they are called and when calling
                  processed=NULL, 
                  ... 
) {
  oricall <- mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  
  ## Preventing confusions
  if (!is.null(mc$ranPars)) {
    stop("incorrect 'ranPars' argument in fitme() call. Use 'fixed' (ranPars is for HLCor() only)")
  }
  if (!is.null(mc$ranFix)) {
    stop("incorrect 'ranFix' argument in fitme() call. Use 'fixed' (ranFix is for HLfit() and corrHLfit() only)")
  }
  if ( ! (is.list(lower) && is.list(upper))) {
    wrongclass <- setdiff(unique(c(class(lower),class(upper))),"list")
    stop(paste("'lower' and 'upper' must be of class list, not",paste(wrongclass,collapse=" or ")))
  } ## as.list() would flatten rho vectors
  #
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(fitme)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)>0) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in fitme call."))
  # 
  if (is.null(processed)) {
    family <- checkRespFam(family)
    FHF <- formals(HLfit) ## makes sure about default values 
    names_FHF <- names(FHF)
    if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
    names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
    FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
    preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(preprocess)))] 
    preprocess.formal.args$family <- family ## already checked 
    preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
    preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
    preprocess.formal.args$ranFix <- fixed ## because preprocess expects ranFix
    preprocess.formal.args$HLmethod <- HLmethod ## forces evaluation
    #
    famfam <- family$family
    if ( identical(family$family,"multi")) {
      ## then data are reformatted as a list. Both HLCor and HLfit can analyse such lists for given corrPars and return the joint likelihood
      ## By contrast HLCor should not fit different corrPars to each data, so it does not lapply("corrHLfit",...)
      ## Rather, it calls preprocess which will construct a list of processed objects, to be used conjointly with the data list.
      ## But then we should not attempt to modify an element of 'pocessed' as if it was a single processed object
      ## We must use setProcessed / getProcessed to access list elements.
      if ( ! inherits(data,"list")) {
        familyargs <- family
        familyargs$family <- NULL
        familyargs$binfamily <- NULL
        ## we need the data list in the corrHLfit envir for the call to makeCheckGeoMatrices
        preprocess.formal.args$data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }     
    }
    mc$processed <- do.call(preprocess,preprocess.formal.args,envir=environment(formula))
    mc$verbose <- reformat_verbose(eval(mc$verbose),For="corrHLfit")
    ## removing all elements that are matched in processed:
    mc$data <- NULL
    mc$family <- NULL
    mc$formula <- NULL
    mc$prior.weights <- NULL
    mc$HLmethod <- NULL ## processed$HL  
    mc$rand.family <- NULL ## processed$rand.families  
    mc$control.glm <- NULL ## processed$control.glm  
    mc$resid.formula <- NULL ## mc$resid.model  
    mc$REMLformula <- NULL ## processed$REMLformula
  }  
  
  mc[[1L]] <- quote(spaMM::fitme_body) 
  hlcor <- eval(mc,parent.frame()) 
  attr(hlcor,"fitmecall") <- oricall ## this says the hlcor was returned by corrHLfit (which should be fitme?)
  return(hlcor)
}


## next steps are COMPoisson optim, detection des modeles lambda, gestion des lambda multiples
#  et meilleurs default phi, lambda 
fitme_body <- function(processed,
                       init=list(),
                       init.HLfit=list(),
                       fixed=list(), ## NULL not valid (shouldbe handled in preprocess?)
                       lower=list(),upper=list(),
                       control.dist=list(),
                       control=list(), ## optim.scale, Optimizer, optimizer.args, maxIter, precision, refit
                       verbose, ## provided by fitme
                       ... ## cf dotnames processing below 
) {
  
  dotlist <- list(...) ## forces evaluations, which makes programming easier...
  data <- getProcessed(processed,"data") ## gets a data list if processed is a list
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
  ## fill HLCor.args
  good_dotnames <- intersect(names(dotlist),HLnames) ## those specifically for the called fns as def'd by HLnames
  if (length(good_dotnames)>0L) {
    HLCor.args <- dotlist[good_dotnames]
  } else HLCor.args <- list() 
  typelist <- list() ## on veut une list pour pouvoir supp des elements par <- NULL
  typelist[names(fixed)] <- "fix"
  attr(fixed,"type") <- typelist
  ## replace some HLCor.args members  
  if ( ! is.null(attr(processed,"multiple"))) {
    pnames <- names(processed[[1]])
  } else pnames <- names(processed)
  for (st in pnames) HLCor.args[st] <- NULL 
  # 'processed' may be modified below, then will be copied in HLCor.args (and then removed from this envir for safety)
  HLCor.args$verbose <- verbose ## fitme_body fn argument
  
  RHOMAX <- NULL
  NUMAX <- 50
  
  init.optim <- init ## usages of init$rho, etc. must be tracked through init.optim too  
  optim.scale <- control$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  #
  spatial.terms <- findSpatial(getProcessed(processed,"predictor",from=1L))
  spatial.model <- spatial.terms[[1L]] 
  corr.model <- paste("",as.character(spatial.model[[1L]]),sep="") ## so that LHS= "" if spatial.model[[1L]] is NULL
  #
  ############### (almost) always check geo info ###################
  if (corr.model  %in% c("Matern","ar1","AR1")) {
    if ( inherits(data,"list")) {
      dataForCheck <- data[[1]]
    } else dataForCheck <- data
    coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(dataForCheck))
  } ##
  #
  ## distMatrix or uniqueGeo potentially added to HLCor.args:
  if (corr.model %in% "ar1") {
    # we need all intermediate levels even those not represented in the data
    levelrange <- range(as.integer(data[,coordinates]))
    tmax <- levelrange[2]-levelrange[1]+1L ## time zero IS the first level of 'coordinate' 
    control.dist$ar1_tmaxs <- tmax
    if (useAdjacencyApprox <-FALSE) {
      adj <- Diagonal(tmax,0L)
      diag(adj[-1,]) <- 1L
      diag(adj[,-1]) <- 1L
      HLCor.args$adjMatrix <- as.matrix(adj) ## suitable for *approximation* of AR1 by adjacency model
    }
  } 
  #
  # modify HLCor.args and <>bounds
  if (corr.model  %in% c("SAR_WWt","adjacency")) { #!# "ar1"
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
  }
  #
  Fixrho <- getPar(fixed,"rho")
  if ( (! is.null(Fixrho)) && (! is.null(init$rho)) ) {
    stop("(!) 'rho' given as element of both 'fixed' and 'init'. Check call.")    
  } else rho.size <- max(length(Fixrho),length(init$rho))
  if ( (! is.null(getPar(fixed,"nu"))) && (! is.null(init$nu)) ) {
    stop("(!) 'nu' given as element of both 'fixed' and 'init'. Check call.")    
  }
  if ( (! is.null(getPar(fixed,"ARphi"))) && (! is.null(init$ARphi)) ) {
    stop("(!) 'ARphi' given as element of both 'fixed' and 'init'. Check call.")    
  }
  if ( (! is.null(getPar(fixed,"Nugget"))) && (! is.null(init$Nugget)) ) {
    stop("(!) 'Nugget' given as element of both 'fixed' and 'init'. Check call.")    
  }
  #
  if (corr.model %in% c("Matern","AR1")) {
    rho_mapping <- control.dist$rho.mapping ## may be NULL
    if (is.null(rho_mapping) ) { ## && length(coordinates)>1L ?
      if (length(coordinates)==rho.size) { ## rho.size comes from explicit rho from user
        rho_mapping <- seq_len(rho.size)           
        names(rho_mapping) <- coordinates
        control.dist$rho.mapping <- rho_mapping
      } else if (length(rho.size)>1L) stop("'rho.mapping' missing with no obvious default from the other arguments.")
    } ## then (for given corr.model's) there is rho_mapping
    locarglist <- list(data=data,distMatrix=dotlist$distMatrix,
                      uniqueGeo=HLCor.args$uniqueGeo, ## typically NULL
                      coordinates=coordinates)
    if(!is.null(dist.method <- control.dist$dist.method)) locarglist$dist.method <- dist.method
    geoMats <- do.call(makeCheckGeoMatrices,locarglist) ## typically computes uniqueGeo
    nbUnique <- geoMats$nbUnique  
    if (inherits(nbUnique,"list")) nbUnique <- mean(unlist(nbUnique))
    distMatrix <- geoMats$distMatrix   
    uniqueGeo <- geoMats$uniqueGeo   
    if (rho.size<2) { ## can be 0 if no explicit rho in the input  
      HLCor.args$distMatrix <- distMatrix   ## determined locally
    } else {
      HLCor.args$uniqueGeo <- uniqueGeo   ## determined locally
    }
  } else nbUnique <- NULL
  #
  phiform <- getProcessed(processed,"resid.predictor",from=1L)
  ##### init.optim$phi/lambda will affect calc_inits -> calc_inits_dispPars.
  # outer estim seems useful when we can suppress all inner estim (thus the hatval calculations). 
  # ./. Therefore, we need to identify all cases where phi is fixed, 
  # ./. or can be outer optimized jointly with lambda:  
  phimodel <- getProcessed(processed,"models",from=1L)[["phi"]]
  if (phimodel=="phiGLM") {message("'fitme' not optimized for models with structured dispersion.")} ## FR->FR test this later"
  phi.Fix <- getProcessed(processed,"phi.Fix",from=1L)
  if ( is.null(phi.Fix) && is.null(init.optim$phi)) {
    get_init_phi <- getProcessed(processed,"get_init_phi",from=1L)
  }   ## v1.8.23: $get_init_phi is non-NULL for "phiScal" but not "phiGLM"
  ## All cases where phi is fixed, or can be outer optimized jointly with lambda:  
  if ( ! is.null(phi.Fix) ||
       ! is.null(init.optim$phi) || ## but its dubious that phiGLM is handled 
       ! is.null(get_init_phi) ) {
    # (1) set phi
    if (is.null(phi.Fix)) {
      if (is.null(init.optim$phi)) {
        if (is.call(processed$prior.weights)) {
          init.optim$phi <- get_init_phi(processed,weights=NULL)
        } else init.optim$phi <- get_init_phi(processed,weights=processed$prior.weights)
      }  
      ## FR->FR fitme may fail obscurely if get_init_phi(processed) fails silently
      if (is.null(init.optim$phi)) { ## get_init_phi returned NULL if ./. 
        # ./. no replicate obs is available, or
        # ./. ULI too large ?? Grey area here FR->FR
        init.optim$phi <- 0.1 ## assuming a single phi coefficient (as the enclosing code)
      }  
    }
    # (2) set lambda (handling incomplete ranFix$lambda vectors)
    if (is.null(init.optim$lambda)) {
      lFix <- getProcessed(processed,"lambda.Fix",1L)
      if (anyNA(lFix)) {
        nranefs <- length(getProcessed(processed,"ZAlist",1L))
        # Low initial values  of lambda should be used to avoid spuriously high logLik from Laplace for high lambda  
        ## .../... unless Laplace is exact...
        init.optim$lambda <- rep(0.1/nranefs,length(which(is.na(lFix)))) 
      }
    }
  }
  family <- getProcessed(processed,"family",from=1L)
  if (family$family=="COMPoisson") {
    checknu <- try(environment(family$aic)$nu,silent=TRUE)
    if (inherits(checknu,"try-error")) init.optim$COMP_nu <- 1 ## FR->FR what if not incuded in the range ?
  } else if (family$family=="negbin") {
    checktheta <- try(environment(family$aic)$shape,silent=TRUE)
    if (inherits(checktheta,"try-error")) init.optim$NB_shape <- 1 ## FR->FR idem ...
  }
  arglist <- list(init.optim=init.optim,init.HLfit=init.HLfit,ranFix=fixed,
                  corr.model=corr.model,optim.scale=optim.scale,
                  control.dist=control.dist,For="fitme")
  if (corr.model %in% c("Matern")) {
    maxrange <- calc_maxrange(rho.size,distMatrix,uniqueGeo,rho_mapping,dist.method) 
    RHOMAX <- 1.1*30*nbUnique/maxrange ## matches init$rho in calc_inits()
    arglist <-c(arglist,list(maxrange=maxrange,RHOMAX=RHOMAX,NUMAX=NUMAX))
  }
  if (corr.model %in% c("SAR_WWt","adjacency") #!# "ar1" rhorange suppressed 
      &&  is.null(getPar(fixed,"rho"))
      && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho$rho NULL or NA) 
  ) arglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
  inits <- do.call("calc_inits",arglist)
  init <- inits$`init` ## keeps all init values, all in untransformed scale
  init.optim <- inits$`init.optim` ## subset of all estimands, as name implies, and in transformed scale
  init.HLfit <- inits$`init.HLfit` ## subset as name implies 
  #
  user.lower <- lower; user.upper <- upper ## keep user input 
  if ("lambda" %in% c(names(user.lower),names(user.lower)) 
      && is.null(init$lambda)) {
    stop("'lambda' in 'lower' or 'upper' has no effect if absent from 'init'.")
  }
  ################
  LUarglist <- list(canon.init=init,
                    init.optim=init.optim, ## initially with right transformed variables but wrong values
                    user.lower=user.lower,user.upper=user.upper,
                    corr.model=corr.model,nbUnique=nbUnique,
                    ranFix=fixed,control.dist=control.dist,
                    optim.scale=optim.scale, RHOMAX=RHOMAX,NUMAX=NUMAX)
  if (corr.model %in% c("SAR_WWt","adjacency") #!# "ar1" rhorange suppressed 
      &&  is.null(getPar(fixed,"rho"))
      && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho$rho NULL or NA) 
  ) { 
    LUarglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
  }
  LowUp <- do.call("makeLowerUpper",LUarglist)
  ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  #
  processedHL1 <- getProcessed(processed,"HL[1]",from=1L) ## there's also HLmethod in processed<[[]]>$callargs
  if (!is.null(processedHL1) && processedHL1=="SEM" && length(lower)>0) {
    optimMethod <- "iterateSEMSmooth"
    # processed <- setProcessed(processed,"SEMargs$SEMseed",value="NULL") ## removing default SEMseed
    ## : bc SEMseed OK to control individual SEMs but not  series of SEM 
    if (is.null(getProcessed(processed,"SEMargs$control_pmvnorm$maxpts",from=1L))) {
      if (length(LowUp$lower>0L)) {
        processed <- setProcessed(processed,"SEMargs$control_pmvnorm$maxpts",value="quote(250L*nobs)") 
      } ## else default visible in SEMbetalambda
    }
  } else optimMethod <- "nloptr"
  needHLCor_specific_args <- (length(setdiff(names(lower),c("trPhi","trLambda","COMP_nu","NB_shape")))>0L
                              # 'empty' corr.model is "" and the match() is then NA
                              || ! is.na(match(corr.model,c("Matern","adjacency","AR1","ar1","corrMatrix")))  
                              #|| length(intersect(c("distMatrix","uniqueGeo","adjMatrix","corrMatrix"),HLCor.args))>0L
                              ) 
  if (needHLCor_specific_args) {
    HLcallfn.obj <- "HLCor.obj"
    HLcallfn <- "HLCor"
    HLCor.args$corr.model <- corr.model 
    HLCor.args$control.dist <- control.dist ## modified locally
    # Subtlety is in HLCor_body: ranPars argument of HLCor contains both fixed and estimated parameters:
    # HL.info$init.HLfit[varNames] <- ranPars[varNames] 
    # ranPars is used to carry both part of the init_HLfit info, and other info.  
    # BUT this is not true if we directly call HLfit_body
    RanFixOrVar <- fixed 
    varNames <- names(init.HLfit) ## hence those that will be variable within HLfit
    varNames <- setdiff(varNames,c("fixef","v_h"))
    RanFixOrVar[varNames] <- init.HLfit[varNames] ## FR->FR duplicat (?) qui montre qu'un attribute serait mieux
    attr(RanFixOrVar,"type")[varNames] <- "var"  
    HLCor.args$ranPars <- RanFixOrVar  
  } else {
    HLcallfn.obj <- "HLfit.obj"
    HLcallfn <- "HLfit"
    HLCor.args$ranFix <- fixed  
  }
  HLCor.args$init.HLfit <- init.HLfit[intersect(names(init.HLfit),c("fixef","v_h"))] ## 
  HLCor.args$processed <- processed ## for the <...>.obj and <...>_body functions  
  processed <- "'processed' erased after copy in 'HLCor.args' to make sure it is not modified later"
  ## 
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  if (needHLCor_specific_args) {
    anyHLCor_obj_args$skeleton <- structure(init.optim,RHOMAX=RHOMAX,NUMAX=NUMAX) ## logscale, only used by HLCor.obj
  } else {
    anyHLCor_obj_args$skeleton <- init.optim
  }
  attr(anyHLCor_obj_args$skeleton,"type") <- list() ## declares a list of typeS of elemnts of skeleton
  attr(anyHLCor_obj_args$skeleton,"type")[names(init.optim)] <- "fix" # fixed within the HLCor call 
  anyHLCor_obj_args$processed <- setProcessed(anyHLCor_obj_args$processed,
                                              "return_only",
                                              value = paste("\"",anyHLCor_obj_args$processed$objective,"APHLs\"",sep=""))
  initvec <- unlist(init.optim)
  ####    tmpName <- generateName("HLtmp") ## tmpName is a string such as "HLtmp0"
  #    anyHLCor_obj_args$init.HLfit <- tmpName 
  ####    assign(tmpName,list(),pos=".GlobalEnv") ## sets HLtmp0 (or a similarly named variable) at the global level
  if (length(initvec)>0L) {
    if (identical(anyHLCor_obj_args$verbose["getCall"][[1L]],TRUE)) {
      refitPars <- init.optim
    } else {
      if (optimMethod=="iterateSEMSmooth") {
        stop("reimplement iterateSEMSmooth() in fitme() later")
        
        ## and then we need to put back the code for logL from smoothing at the end of this function
        
      } else { ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null 
        objfn_nloptr <- function(x,anyHLCor_obj_args) { ## all functions should have the same args.
          arglist <- c(list(ranefParsVec=x),anyHLCor_obj_args)
          return( - do.call(HLcallfn.obj,arglist))
        }
        nloptr_controls <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1.0e-4,maxeval=-1,print_level=0) ## DEFAULT
        nloptr_controls[names(control$nloptr)] <- control$nloptr ## Overwrite defaults with any element of $nloptr
        optr <- nloptr(x0=initvec,eval_f=objfn_nloptr,lb=unlist(lower),ub=unlist(upper),
                       opts=nloptr_controls,anyHLCor_obj_args=anyHLCor_obj_args)
        optPars <- relist(optr$solution,init.optim)
        refitPars <- optPars <- structure(optPars,method="nloptr",optr=optr) 
        ## full optr is big. Contains the call and an environment... 
      }
    }
    RanFixOrVar <- fixed 
    if (( identical(control$refit,TRUE) || identical(control$refit$phi,TRUE)) && ! is.null(refitPars$trPhi)) {
      HLCor.args$init.HLfit$phi <- dispInv(refitPars$trPhi)
      RanFixOrVar$trPhi <- refitPars$trPhi <- NULL
    }
    if (( identical(control$refit,TRUE) || identical(control$refit$lambda,TRUE)) && ! is.null(refitPars$trLambda)) {
      HLCor.args$init.HLfit$lambda <- dispInv(refitPars$trLambda)
      RanFixOrVar$trLambda <- refitPars$trLambda <- NULL
    }
    RanFixOrVar[names(refitPars)] <- refitPars ## avoids overwriting fixed ran pars 
    attr(RanFixOrVar,"type")[names(refitPars)] <- "outer" ##  
    if (needHLCor_specific_args) {
      attr(RanFixOrVar,"RHOMAX") <- RHOMAX
      attr(RanFixOrVar,"NUMAX") <- NUMAX
      HLCor.args$ranPars <- RanFixOrVar ## variable locally
    } else HLCor.args$ranFix <- RanFixOrVar ## variable locally  
  }
  #
  verbose["warn"] <- TRUE ## important!
  HLCor.args$verbose <- verbose ## modified locally
  hlcor <- do.call(HLcallfn,HLCor.args) ## recomputation post optimization, or only computation if length(initvec)=0
  if (identical(HLCor.args$verbose["getCall"][[1L]],TRUE)) {
    return(hlcor) ## HLCorcall
  }
  #
  if (needHLCor_specific_args) {
    if ( is.null(HLCor.args$adjMatrix) && is.null(HLCor.args$corrMatrix)) { ## Matern, AR1/ar1; 
      ## info.uniqueGeo,  used by predict -> calcNewCorrs and mapMM,  and getDistMat for Matern,
      ##   may already have been set by HLCor_body (if uniqueGeo was one of its explict args)
      if (is.null(attr(hlcor,"info.uniqueGeo"))) attr(hlcor,"info.uniqueGeo") <- HLCor.args$uniqueGeo ## Matern with rho.size>1
      if (is.null(attr(hlcor,"info.uniqueGeo"))) attr(hlcor,"info.uniqueGeo") <- uniqueGeo ## Matern with rho.size=1
      ## still NULL in AR1/ar1 model
    }
  }
  #
  if (length(initvec)>0L) {
    attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, 
                                    objective=anyHLCor_obj_args$processed$objective) ## processed was erased for safety
  }
  return(hlcor)
}