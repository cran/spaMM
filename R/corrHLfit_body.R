corrHLfit_body <- function(processed,
                           init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, Optimizer, <optimizer controls>
                      nb_cores=NULL,
                      ... ## cf dotnames processing below
) { 

  

  #########################################################
  #
  dotlist <- list(...) ## forces evaluations, which makes programming easier...
  data <- .getProcessed(processed,"data") ## gets a data list if processed is a list
  verbose <-  .getProcessed(processed,"verbose",from=1L)
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
  ## fill HLCor.args
  good_dotnames <- intersect(names(dotlist),HLnames) ## those specifically for the called fns as def'd by HLnames
  if (length(good_dotnames)>0L) {
    HLCor.args <- dotlist[good_dotnames]
  } else HLCor.args <- list() 
  typelist <- list() ## on veut une list pour pouvoir supp des elements par <- NULL
  typelist[names(ranFix)] <- "fix"
  attr(ranFix,"type") <- typelist
  ## replace some HLCor.args members  
  if ( is.list(processed)) { ## list of environments
    pnames <- names(processed[[1]])
  } else pnames <- names(processed)
  for (st in pnames) HLCor.args[st] <- NULL 
  # 'processed' may be modified below, then will be copied in HLCor.args (and then removed from this envir for safety) => NO NEED to copy any processed element in HLCor.args now.

  RHOMAX <- NULL
  NUMAX <- 50
  
  init.optim <- init.corrHLfit ## usages of init.corrHLfit$rho, etc. must be tracked through init.optim too  
  optim.scale <- control.corrHLfit$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  # test-Nugget -> different p_v / p_bv compromise (!)
  spatial.terms <- .findSpatial(.getProcessed(processed,"predictor",from=1L))
  spatial.model <- spatial.terms[[1]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1]]) 
  } else {
    # corr.model <- "Matern" ## up to v 1.0; for the defective syntax of the scripts for the Ecography paper
    stop("spatial correlation model not specified in 'formula': was valid in version 1.0 but not later.")
  } 
  HLCor.args$corr.model <- corr.model ## determined locally
  if (corr.model== "corrMatrix") { ## could cover all specifications of a constant 'cov' Matrix without correlation parameter
    corrMatrix <- HLCor.args$corrMatrix
    if (is.null(corrMatrix)) corrMatrix <- .get_corr_prec_from_covStruct(HLCor.args$covStruct)
    .check_corrMatrix(corrMatrix)
    sparse_precision <- inherits(corrMatrix,"precision")
    if ( ! sparse_precision) sparse_precision <- .determine_sparse_precision(processed, corr.model, init.HLfit) ## checks spaMM option
  } else sparse_precision <- .determine_sparse_precision(processed, corr.model, init.HLfit) ## checks spaMM option
  ignored_as_no_need_for_local_copy_processed <- .reset_processed_sparse_precision(processed, sparse_precision)
  #
  ############### (almost) always check geo info ###################
  if (corr.model  %in% c("Matern","AR1")) {
    coordinates <- .get_coordinates(spatial.model=spatial.model, data=data)
    coord_within <- .extract_check_coords_within(spatial.model=spatial.model) 
  }
  if (corr.model %in% "AR1") { ## FIXME only necess for sparse precision ?
    # we need all intermediate levels even those not represented in the data
    levelrange <- range(as.integer(data[,coord_within]))
    tmax <- levelrange[2]-levelrange[1]+1L ## time zero IS the first level of 'coordinate' 
    control.dist$AR1_tmax <- tmax
  }
  #
  # modify HLCor.args and <>bounds
  cov_matrices__from_covStruct <- .get_cov_matrices__from_covStruct(HLCor.args$covStruct, corr.model, HLCor.args)
  if ( ! is.null(cov_matrices__from_covStruct)) for (st in names(cov_matrices__from_covStruct)) HLCor.args[[st]] <- cov_matrices__from_covStruct[[st]]
  HLCor.args$covStruct <- NULL # so that is it not re-processed by HLCor_body later 
  
  if (corr.model  %in% c("SAR_WWt","adjacency")) { 
    decomp <- .provide_AR_factorization(HLCor.args, sparse_precision, corr.model)
    if (corr.model=="SAR_WWt") attr(HLCor.args$adjMatrix,"UDU.") <- decomp
    if (corr.model=="adjacency") attr(HLCor.args$adjMatrix,"symSVD") <- decomp 
    rhorange <- sort(1/range(decomp$d)) ## keeping in mind that the bounds can be <>0
    if(.getProcessed(processed,"verbose",from=1L)["SEM"]) cat(paste("Feasible rho range: ",paste(signif(rhorange,6),collapse=" -- "),"\n"))
    ## added 11/2016:
    if ( ! is.null(lower$rho)) rhorange[1L] <- max(rhorange[1L],lower$rho)
    if ( ! is.null(upper$rho)) rhorange[2L] <- min(rhorange[2L],upper$rho)
  }
  #
  rho.size <- .check_conflict_init_fixed(ranFix,init.corrHLfit,
                             "given as element of both 'ranFix' and 'init.corrHLfit'. Check call.") ## NOTfixme: difference bw fitme et corrHLfit
  #
  ## distMatrix or uniqueGeo potentially added to HLCor.args:
  if (corr.model %in% c("Matern","AR1")) {
    coords_nesting <- setdiff(coordinates,coord_within)
    if (length(coords_nesting)) {
      if (.getProcessed(processed,"sparsePrecisionBOOL",1L)) {
        message("sparse precision method not available for nested AR1 effects.")
        .setProcessed(processed,"sparsePrecisionBOOL","FALSE")
      }
    }
    locarglist <- list(data=data,distMatrix=dotlist$distMatrix,
                       uniqueGeo=HLCor.args$uniqueGeo, ## typically NULL
                       coords_nesting=coords_nesting, coordinates=coordinates, dist.method = control.dist$dist.method)
    geoMats <- do.call(".makeCheckGeoMatrices",locarglist) ## computes matrices which are NULL on input
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
  family <- .getProcessed(processed,"family",from=1L)
  if (family$family=="COMPoisson") {
    checknu <- try(environment(family$aic)$nu,silent=TRUE)
    if (inherits(checknu,"try-error")) init.optim$COMP_nu <- 1
  } else if (family$family=="negbin") {
    checktheta <- try(environment(family$aic)$shape,silent=TRUE)
    if (inherits(checktheta,"try-error")) init.optim$NB_shape <- 1
  }
  #
  calc_inits_arglist <- list(init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix,
                  corr.model=corr.model,optim.scale=optim.scale,
                  control.dist=control.dist,For="corrHLfit") ## NOTfixme difference fitme/corrHLfit
  if (corr.model %in% c("Matern")) {
    control.dist$rho.mapping <- .provide_rho_mapping(control.dist, coordinates, rho.size)
    maxrange <- .calc_maxrange(rho.size,distMatrix,uniqueGeo,control.dist$rho.mapping,control.dist$dist.method) 
    RHOMAX <- 1.1*30*nbUnique/maxrange ## matches init$rho in calc_inits()
    calc_inits_arglist <-c(calc_inits_arglist,list(maxrange=maxrange,RHOMAX=RHOMAX,NUMAX=NUMAX))
  }
  if (corr.model %in% c("SAR_WWt","adjacency")
      &&  is.null(.getPar(ranFix,"rho")) ## NOTfixme difference fitme/corrHLfit
      && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho NULL or NA) 
  ) calc_inits_arglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
  inits <- do.call(".calc_inits",calc_inits_arglist)
  #
  init <- inits$`init` ## keeps all init values, all in untransformed scale
  init.optim <- inits$`init.optim` ## subset of all optim estimands, as name implies, and in transformed scale
  init.HLfit <- inits$`init.HLfit` ## subset as name implies 
  #
  ## maximization VIA HLCor.obj
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
  LUarglist <- list(canon.init=init, ## canonical scale
                    init.optim=init.optim, ## transformed scale, used only to initialize lower and upper 
                    user.lower=user.lower,user.upper=user.upper,
                    corr.model=corr.model,nbUnique=nbUnique,
                    ranFix=ranFix,control.dist=control.dist,
                    optim.scale=optim.scale, RHOMAX=RHOMAX,NUMAX=NUMAX)
  if (corr.model %in% c("SAR_WWt","adjacency") 
      &&  is.null(.getPar(ranFix,"rho"))
      && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho$rho NULL or NA) 
      ) { 
    LUarglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
  }
  LowUp <- do.call(".makeLowerUpper",LUarglist)
  ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  #
  ranPars <- ranFix ## ranPars argument of HLCor contains both fixed and estimated parameters:
  varNames <- names(init.HLfit) ## hence those that will be variable within HLfit
  varNames <- setdiff(varNames,c("fixef","v_h"))
  ranPars[varNames] <- init.HLfit[varNames] ## FR->FR duplicat (?) qui montre qu'un attribute serait mieux
  attr(ranPars,"type")[varNames] <- "var"  
  if (spaMM.getOption("wDEVEL2")) {
    parlist <- .merge_parlist(,new=ranPars,types="fix")
    parlist <- .merge_parlist(parlist,new=init.HLfit,types="var")
    attr(ranPars,"parlist") <- parlist
  } 
  HLCor.args$ranPars <- ranPars  ## variable locally
  #
  HLCor.args$control.dist <- control.dist ## modified locally
  processedHL1 <- .getProcessed(processed,"HL[1]",from=1L) ## there's also HLmethod in processed<[[]]>$callargs
  if (!is.null(processedHL1) && processedHL1=="SEM" && length(lower)>0) {
    optimMethod <- "iterateSEMSmooth"
    if (is.null(.getProcessed(processed,"SEMargs$control_pmvnorm$maxpts",from=1L))) {
      if (length(LowUp$lower)>0L) {
        processed <- .setProcessed(processed,"SEMargs$control_pmvnorm$maxpts",value="quote(250L*nobs)") 
      } ## else default visible in SEMbetalambda
    }
  } else optimMethod <- ".locoptim"
  HLCor.args$processed <- processed
  ## 
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  anyHLCor_obj_args$skeleton <- structure(init.optim,RHOMAX=RHOMAX,NUMAX=NUMAX) ## logscale, only used by HLCor.obj
  attr(anyHLCor_obj_args$skeleton,"type") <- list() ## declares a list of typeS of elemnts of skeleton
  attr(anyHLCor_obj_args$skeleton,"type")[names(init.optim)] <- "fix" # fixed with the HLCor call 
  # if processed is an envir or a list of envirs, the following is not local to anyHLCor_obj_args$processed but change processed globally
  anyHLCor_obj_args$processed <- .setProcessed(anyHLCor_obj_args$processed,"return_only",
                                              value=paste("\"",HLCor.args$processed$objective,"APHLs\"",sep="")) ## laborious for strings
  initvec <- unlist(init.optim)
  ### .do_TRACE(verbose["TRACE"])
  if (optimMethod=="iterateSEMSmooth") {   
    MAX <- list(trRho=RHOMAX, trNu=NUMAX) ##  MAX is used in the SEMdiagnosticplot...;
    ## its names should match the colnames of the data in Krigobj = the  parameters of the likelihood surface. Current code maybe not general.
    loclist <- list(anyHLCor_obj_args=anyHLCor_obj_args,  ## contains $processed
                    LowUp=LowUp,init.corrHLfit=init.corrHLfit, 
                    #preprocess.formal.args=preprocess.formal.args, 
                    control.corrHLfit=control.corrHLfit,
                    verbose=verbose[["iterateSEM"]],
                    nb_cores=nb_cores,
                    MAX=MAX)
    #optr <- do.call(probitgem::iterateSEMSmooth,loclist) ## will pass CHECK when CRAN knows probitgem
    optr <- eval(as.call(c(quote(iterateSEMSmooth),loclist))) ## if probitgem unknown
    optPars <- relist(optr$par,init.optim)
    if (!is.null(optPars)) attr(optPars,"method") <-"optimthroughSmooth"
  } else { ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null
    if (identical(processed$verbose["getCall"][[1L]],TRUE)) {
      optPars <- init.optim
    } else {
      loclist<-list(init.optim=init.optim,LowUp=LowUp,anyObjfnCall.args=anyHLCor_obj_args,
                    control=control.corrHLfit,maximize=TRUE) 
      optPars <- do.call(".locoptim",loclist)
    }
  }
  ranPars[names(optPars)] <- optPars ## avoids overwriting fixed ranPars; and keep the list class of ranPars...
  attr(ranPars,"type")[names(optPars)] <- "outer" ##  
  attr(ranPars,"RHOMAX") <- RHOMAX
  attr(ranPars,"NUMAX") <- NUMAX
  if (spaMM.getOption("wDEVEL2")) {
    parlist <- .merge_parlist(,new=ranPars,types="fix")
    parlist <- .merge_parlist(parlist,new=optPars,types="outer") ## consistent with earlier code
    # so as to combine outer estimation and SEs of inner estimation
    attr(ranPars,"parlist") <- parlist
  }
  HLCor.args$ranPars <- ranPars ## variable locally
  # if processed is an envir, the following is not local to anyHLCor_obj_args$processed but change processed globally
  HLCor.args$processed <- .setProcessed(HLCor.args$processed,"return_only", value="NULL") 
  HLCor.args$processed <- .setProcessed(HLCor.args$processed,"verbose['warn']", value="TRUE") ## important!
  hlcor <- do.call("HLCor",HLCor.args) ## recomputation post optimization (or only computation, if length(lower)=0)
  if (is.call(hlcor)) { return(hlcor[]) } ## HLCorcall
  if ( is.null(HLCor.args$adjMatrix)  && is.null(HLCor.args$corrMatrix)) { ## Matern, AR1; 
    ## info.uniqueGeo,  used by predict -> calcNewCorrs and mapMM, and getDistMat for Matern,
    ##   may already have been set by HLCor_body (if uniqueGeo was one of its explict args)
    if (is.null(attr(hlcor,"info.uniqueGeo"))) attr(hlcor,"info.uniqueGeo") <- HLCor.args$uniqueGeo ## Matern with rho.size>1
    if (is.null(attr(hlcor,"info.uniqueGeo"))) attr(hlcor,"info.uniqueGeo") <- uniqueGeo ## Matern with rho.size=1
    ## still NULL in AR1 model
  }
  attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, objective=HLCor.args$processed$objective)
  if ( ! is.null(optPars)) {
    locoptr <- attr(optPars,"optr")
    if (attr(optPars,"method")=="nloptr") {
      if (locoptr$status<0L) hlcor$warnings$optimMessage <- paste("nloptr() message: ",
                                                                  locoptr$message," (status=",locoptr$status,")",sep="")
    } else if ( attr(optPars,"method")=="optim" ) {
      if (locoptr$convergence>0L) hlcor$warnings$optimMessage <- paste("optim() message: ",locoptr$message,
                                           " (convergence=",locoptr$convergence,")",sep="")
    } else if ( attr(optPars,"method")== "optimthroughSmooth") {
      # provide logL estimate from the smoothing, to be used rather than the hlcor logL :
      logLapp <- optr$value
      attr(logLapp,"method") <- "  logL (smoothed)" 
      hlcor$APHLs$logLapp <- logLapp
    }
  }
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  return(hlcor) ## it's the call which says it was returned by corrHLfit
}

    
