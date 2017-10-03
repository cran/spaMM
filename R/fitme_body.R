.objfn_nloptr <- function(x, anyHLCor_obj_args, HLcallfn.obj) { ## all functions should have the same args.
  anyHLCor_obj_args$ranefParsVec <- x 
  return( - do.call(HLcallfn.obj, anyHLCor_obj_args))
}

## next steps are COMPoisson optim, detection des modeles lambda, gestion des lambda multiples
#  et meilleurs default phi, lambda 
fitme_body <- function(processed,
                       init=list(),
                       init.HLfit=list(),
                       fixed=list(), ## NULL not valid (shouldbe handled in preprocess?)
                       lower=list(),upper=list(),
                       control.dist=list(),
                       control=list(), ## Optimizer, <optimizer controls>, precision
                       ... ## cf dotnames processing below 
) {
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
  typelist[names(fixed)] <- "fix"
  attr(fixed,"type") <- typelist
  ## replace some HLCor.args members  
  if (  is.list(processed) ) {
    pnames <- names(processed[[1]])
  } else pnames <- names(processed)
  for (st in pnames) HLCor.args[st] <- NULL 
  # 'processed' may be modified below, then will be copied in HLCor.args (and then removed from this envir for safety)

  RHOMAX <- NULL
  NUMAX <- 50
  
  init.optim <- init ## usages of init$rho, etc. must be tracked through init.optim too  
  optim.scale <- control$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  #
  spatial.terms <- .findSpatial(.getProcessed(processed,"predictor",from=1L))
  spatial.model <- spatial.terms[[1L]] 
  corr.model <- paste("",as.character(spatial.model[[1L]]),sep="") ## so that LHS= "" if spatial.model[[1L]] is NULL
  if (corr.model== "corrMatrix") { ## could cover all specifications of a constant 'cov' Matrix without correlation parameter
    corrMatrix <- HLCor.args$corrMatrix
    if (is.null(corrMatrix)) corrMatrix <- .get_corr_prec_from_covStruct(HLCor.args$covStruct)
    .check_corrMatrix(corrMatrix)
    sparse_precision <- inherits(corrMatrix,"precision")
    if ( ! sparse_precision) sparse_precision <- .determine_sparse_precision(processed, corr.model, init.HLfit) ## checks spaMM option
  } else { ## all other models: "", "adjacency", "AR1"...
    sparse_precision <- .determine_sparse_precision(processed, corr.model, init.HLfit)
  }
  ignored_as_no_need_for_local_copy_processed <- .reset_processed_sparse_precision(processed, sparse_precision)
  #
  ############### (almost) always check geo info ###################
  #
  if (corr.model  %in% c("Matern","AR1")) {
    coordinates <- .get_coordinates(spatial.model=spatial.model, data=data)
    coord_within <- .extract_check_coords_within(spatial.model=spatial.model) 
  }
  #
  if (corr.model %in% "AR1") { ## 
    # we need all intermediate levels even those not represented in the data
    levelrange <- range(as.integer(data[,coordinates]))
    tmax <- levelrange[2]-levelrange[1]+1L ## time zero IS the first level of 'coordinate' 
    control.dist$AR1_tmax <- tmax  ## fixme but .spMMFactorList_locfn does something similar
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
    if(verbose["SEM"]) cat(paste("Feasible rho range: ",paste(signif(rhorange,6),collapse=" -- "),"\n"))
    ## added 11/2016:
    if ( ! is.null(lower$rho)) rhorange[1L] <- max(rhorange[1L],lower$rho)
    if ( ! is.null(upper$rho)) rhorange[2L] <- min(rhorange[2L],upper$rho)  
  }
  #
  rho.size <- .check_conflict_init_fixed(fixed,init,
                             "given as element of both 'fixed' and 'init'. Check call.")
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
    geoMats <- do.call(".makeCheckGeoMatrices",locarglist) ## typically computes uniqueGeo
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
  #phiform <- .getProcessed(processed,"resid.predictor",from=1L)
  ##### init.optim$phi/lambda will affect calc_inits -> calc_inits_dispPars.
  # outer estim seems useful when we can suppress all inner estim (thus the hatval calculations). 
  # ./. Therefore, we need to identify all cases where phi is fixed, 
  # ./. or can be outer optimized jointly with lambda:  
  phimodel <- .getProcessed(processed,"models",from=1L)[["phi"]]
  #if (phimodel=="phiGLM") {message("'fitme' not optimized for models with structured dispersion.")} ## FR->FR test this later"
  phi.Fix <- .getProcessed(processed,"phi.Fix",from=1L)
  ## All cases where phi is fixed, or can be outer optimized jointly with lambda:  
  if ( ! is.null(phi.Fix) ||
       ! is.null(init.optim$phi) || ## but its dubious that phiGLM is handled 
       phimodel == "phiScal" ) {
    # (1) set phi
    if (is.null(phi.Fix)) {
      if (is.null(init.optim$phi)) {
        if (is.call(processed$prior.weights)) {
          init.optim$phi <- .get_init_phi(processed,weights=NULL)
        } else init.optim$phi <- .get_init_phi(processed,weights=processed$prior.weights)
      }  
      ## FR->FR fitme may fail obscurely if .get_init_phi(processed) fails silently
      if (is.null(init.optim$phi)) { ## .get_init_phi returned NULL if ./. 
        # ./. no replicate obs is available, or
        # ./. ULI too large ?? Grey area here FR->FR
        init.optim$phi <- 0.1 ## assuming a single phi coefficient (as the enclosing code)
      }  
    }
    # (2) set lambda (handling incomplete ranFix$lambda vectors)
    if (is.null(init.optim$lambda)) {
      lFix <- .getProcessed(processed,"lambda.Fix",1L)
      if (anyNA(lFix)) {
        nrand <- length(.getProcessed(processed,"ZAlist",1L))
        if (nrand>0L && attr(.getProcessed(processed,"ZAlist"),"anyRandomSlope")) 
          stop("fitme() does not yet fit random-slope models.")
        # Low initial values  of lambda should be used to avoid spuriously high logLik from Laplace for high lambda  
        ## .../... unless Laplace is exact...
        init.optim$lambda <- rep(0.1/nrand,length(which(is.na(lFix)))) ## complement to lambda.Fix
      }
    }
  }
  family <- .getProcessed(processed,"family",from=1L)
  if (family$family=="COMPoisson") {
    checknu <- try(environment(family$aic)$nu,silent=TRUE)
    if (inherits(checknu,"try-error")) init.optim$COMP_nu <- 1 ## FR->FR what if not incuded in the range ?
  } else if (family$family=="negbin") {
    checktheta <- try(environment(family$aic)$shape,silent=TRUE)
    if (inherits(checktheta,"try-error")) init.optim$NB_shape <- 1 ## FR->FR idem ...
  }
  # a few 'global' vars here:
  calc_inits_arglist <- list(init.optim=init.optim,init.HLfit=init.HLfit,ranFix=fixed,
                  corr.model=corr.model,optim.scale=optim.scale,
                  control.dist=control.dist,For="fitme")
  if (corr.model %in% c("Matern")) {
    control.dist$rho.mapping <- .provide_rho_mapping(control.dist, coordinates, rho.size)
    maxrange <- .calc_maxrange(rho.size,distMatrix,uniqueGeo,control.dist$rho.mapping,control.dist$dist.method) 
    RHOMAX <- 1.1*30*nbUnique/maxrange ## matches init$rho in calc_inits()
    calc_inits_arglist <-c(calc_inits_arglist,list(maxrange=maxrange,RHOMAX=RHOMAX,NUMAX=NUMAX))
  }
  if (corr.model %in% c("SAR_WWt","adjacency") 
      &&  is.null(.getPar(fixed,"rho"))
      && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho$rho NULL or NA) 
  ) calc_inits_arglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
  inits <- do.call(".calc_inits",calc_inits_arglist)
  #
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
  if (corr.model %in% c("SAR_WWt","adjacency")
      &&  is.null(.getPar(fixed,"rho"))
      && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho$rho NULL or NA) 
  ) { 
    LUarglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
  }
  LowUp <- do.call(".makeLowerUpper",LUarglist)
  ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  #
  processedHL1 <- .getProcessed(processed,"HL[1]",from=1L) ## there's also HLmethod in processed<[[]]>$callargs
  if (!is.null(processedHL1) && processedHL1=="SEM" && length(lower)>0) {
    optimMethod <- "iterateSEMSmooth"
    if (is.null(.getProcessed(processed,"SEMargs$control_pmvnorm$maxpts",from=1L))) {
      if (length(LowUp$lower)>0L) {
        processed <- .setProcessed(processed,"SEMargs$control_pmvnorm$maxpts",value="quote(250L*nobs)") 
      } ## else default visible in SEMbetalambda
    }
  } else optimMethod <- "nloptr"
  needHLCor_specific_args <- (length(setdiff(names(lower),c("trPhi","trLambda","COMP_nu","NB_shape")))>0L
                              # 'empty' corr.model is "" and the match() is then NA
                              || ! is.na(match(corr.model,c("Matern","adjacency","AR1","corrMatrix")))  
                              #|| length(intersect(c("distMatrix","uniqueGeo","adjMatrix","corrMatrix"),HLCor.args))>0L
                              )
  if (spaMM.getOption("wDEVEL2")) {
    parlist <- .merge_parlist(,new=fixed,types="fix")
    parlist <- .merge_parlist(parlist,new=init.HLfit,types="var")
    attr(fixed,"parlist") <- parlist
  } 
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
  anyHLCor_obj_args$processed <- .setProcessed(anyHLCor_obj_args$processed,
                                              "return_only",
                                              value = paste("\"",anyHLCor_obj_args$processed$objective,"APHLs\"",sep=""))
  initvec <- unlist(init.optim)
  ### .do_TRACE(verbose["TRACE"])
  if (length(initvec)>0L) {
    if (identical(processed$verbose["getCall"][[1L]],TRUE)) {
      refitPars <- init.optim
    } else {
      refit_ <- control$refit ## may be a boolean or a list of booleans 
      if (optimMethod=="iterateSEMSmooth") {
        if (is.null(refit_)) refit_ <- FALSE
        stop("reimplement iterateSEMSmooth() in fitme() later")
        
        ## and then we need to put back the code for logL from smoothing at the end of this function
        
        ## beware of SEM case with length(lower)=0
        
      } else { ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null 
        nloptr_controls <- spaMM.getOption("nloptr")
        nloptr_controls[names(control$nloptr)] <- control$nloptr ## Overwrite defaults with any element of $nloptr
        lowerb <- unlist(lower)
        upperb <- unlist(upper)
        optr <- nloptr::nloptr(x0=initvec,eval_f=.objfn_nloptr,lb=lowerb,ub=upperb,
                       opts=nloptr_controls, anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj)
        while (optr$status==5L) { ## optr$status=5 => termination bc maxeval has been reached 
          # met status=4: nloptr message in normal termination due toxtol_rel, but is this true ?
          prevlik <- optr$objective
          reinit <- pmax(lowerb,pmin(upperb,optr$solution))
          optr <- nloptr::nloptr(x0=reinit,eval_f=.objfn_nloptr,lb=lowerb,ub=upperb,
                                 opts=nloptr_controls, anyHLCor_obj_args=anyHLCor_obj_args, HLcallfn.obj=HLcallfn.obj)
          loc_ftol <- max(1e-8, optr$options$ftol_abs)
          if (- optr$objective < - prevlik+loc_ftol) break ## no progress in <= maxeval iterations
        }
        optPars <- relist(optr$solution,init.optim)
        if (is.null(refit_)) {
          refit_ <- (nloptr_controls$xtol_rel > (5e-6 + 1e-8))
          ## or perhaps: 
          # canon <- .canonizeRanPars(optPars,corr.model="",checkComplete = FALSE) ## fixme should not need the checkComplete when model=""?
          # canon <- unlist(canon$ranPars[c("phi","lambda")])
          # if (! is.null(canon) && any(canon< 5e-5)) refit_ <- TRUE
        } 
        ## full optr is big. We take out the two items that contribute much to saveSize:
        optr$eval_f <- NULL
        optr$nloptr_environment <- NULL
        refitPars <- optPars <- structure(optPars,method="nloptr",optr=optr) 
      }
    }
    RanFixOrVar <- fixed 
    # refit_ is list if so provided by user, else typically boolean
    if ( is.list(refit_)) {refit_phi <- refit_$phi} else refit_phi <- refit_ ## result may be NULL in list case
    if (identical(refit_phi,TRUE) && ! is.null(refitPars$trPhi)) {
      HLCor.args$init.HLfit$phi <- .dispInv(refitPars$trPhi)
      #if (any(HLCor.args$init.HLfit$phi<nloptr_controls$xtol_rel)) .setProcessed(processed,"spaMM_tol$Xtol_abs","1e-7") ## FIXME ? not set back after refit
      RanFixOrVar$trPhi <- refitPars$trPhi <- NULL
    }
    if ( is.list(refit_)) {refit_lambda <- refit_$lambda} else refit_lambda <- refit_
    if (identical(refit_lambda,TRUE) && ! is.null(refitPars$trLambda)) {
      HLCor.args$init.HLfit$lambda <- .dispInv(refitPars$trLambda)
      #if (any(HLCor.args$init.HLfit$lambda<nloptr_controls$xtol_rel)) .setProcessed(processed,"spaMM_tol$Xtol_abs","1e-7") ## FIXME ? not set back after refit
      RanFixOrVar$trLambda <- refitPars$trLambda <- NULL
    }
    RanFixOrVar[names(refitPars)] <- refitPars ## avoids overwriting fixed ran pars 
    attr(RanFixOrVar,"type")[names(refitPars)] <- "outer" ##  
    if (spaMM.getOption("wDEVEL2")) {
      parlist <- .merge_parlist(,new=RanFixOrVar,types="fix")
      parlist <- .merge_parlist(parlist,new=refitPars,types="outer") ## consistent with earlier code
      # so as to combine outer estimation and SEs of inner estimation
      attr(RanFixOrVar,"parlist") <- parlist
    }
    if (needHLCor_specific_args) {
      attr(RanFixOrVar,"RHOMAX") <- RHOMAX
      attr(RanFixOrVar,"NUMAX") <- NUMAX
      HLCor.args$ranPars <- RanFixOrVar ## variable locally
    } else HLCor.args$ranFix <- RanFixOrVar ## variable locally  
  }
  #
  # if processed is an envir, the following is not local to anyHLCor_obj_args$processed but change processed globally
  HLCor.args$processed <- .setProcessed(HLCor.args$processed,"return_only", value="NULL") 
  HLCor.args$processed <- .setProcessed(HLCor.args$processed,"verbose['warn']", value="TRUE") ## important!
  hlcor <- do.call(HLcallfn,HLCor.args) ## recomputation post optimization, or only computation if length(initvec)=0
  if (is.call(hlcor)) {
    ## then do.call(HLcallfn,HLCor.args) has retuned the call, not the fit. 
    ## see def of get_HLCorcall() for further explanation
    return(hlcor) ## HLCorcall
  }
  #
  if (needHLCor_specific_args) {
    if ( is.null(HLCor.args$adjMatrix) && is.null(HLCor.args$corrMatrix)) { ## Matern, AR1; 
      ## info.uniqueGeo,  used by predict -> calcNewCorrs and mapMM,  and getDistMat for Matern,
      ##   may already have been set by HLCor_body (if uniqueGeo was one of its explict args)
      if (is.null(attr(hlcor,"info.uniqueGeo"))) attr(hlcor,"info.uniqueGeo") <- HLCor.args$uniqueGeo ## Matern with rho.size>1
      if (is.null(attr(hlcor,"info.uniqueGeo"))) attr(hlcor,"info.uniqueGeo") <- uniqueGeo ## Matern with rho.size=1
      ## still NULL in AR1 model
    }
  }
  #
  if (length(initvec)>0L) {
    attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, 
                                    objective=anyHLCor_obj_args$processed$objective) ## processed was erased for safety
  }
  ## substantial effect on object size! :
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  ##
  return(hlcor)
}