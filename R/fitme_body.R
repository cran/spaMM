## next steps are COMPoisson optim, detection des modeles lambda, gestion des lambda multiples
#  et meilleurs default phi, lambda 
fitme_body <- function(processed,
                       init=list(),
                       init.HLfit=list(),
                       fixed=list(), ## NULL not valid (shouldbe handled in preprocess?)
                       lower=list(),upper=list(),
                       control.dist=list(),
                       control=list(), ## optimizer, <optimizer controls>, precision
                       ... ## cf dotnames processing below 
) {
  dotlist <- list(...) ## forces evaluations, which makes programming easier...
  if (is.list(processed)) {
    proc1 <- processed[[1L]]
  } else proc1 <- processed
  verbose <-  proc1$verbose
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
  ## fill HLCor.args
  good_dotnames <- intersect(names(dotlist),HLnames) ## those specifically for the called fns as def'd by HLnames
  if (length(good_dotnames)) {
    HLCor.args <- dotlist[good_dotnames]
  } else HLCor.args <- list() 
  typelist <- list() ## on veut une list pour pouvoir supp des elements par <- NULL
  # for (st in names(fixed)) {
  #   typ <- rep("fix",length(fixed[[st]]))
  #   if (is.list(fixed[[st]])) {
  #     typ[sapply(fixed[[st]],is.null)] <- "var"
  #   } else typ[is.na(fixed[[st]])] <- "var"
  #   typelist[[st]] <- typ
  # }
  ## Il semble difficile de garder une info sur chaque element de lambda ou ranCoefs, particulierement parce que 
  #  les elements NULL de ranCoefs poseraient probleme pour relist(). Il faut plutôt utiliser les noms.
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
  
  init.optim <- user_init_optim <- init ## user_init_optim only for a check in new_locoptim, true initial value init.optim is modified below
  optim.scale <- control[["optim.scale"]] 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  sparse_precision <- processed$sparsePrecisionBOOL
  #
  corr_info <- proc1$corr_info 
  #
  # modify HLCor.args and <>bounds;   ## distMatrix or uniqueGeo potentially added to HLCor.args:
  corr_types <- corr_info$corr_types
  spatial_terms <- attr(proc1$ZAlist,'exp_spatial_terms')
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if (! is.na(corr_type)) {
      if (corr_type  %in% c("SAR_WWt","adjacency")) { 
        # We need the $d t odetermine rhorange.
        # given the decomp is computed, we try to make it available for other computations.
        # But others computs shouls not assume it is available.
        # Further, if is.list(processed) is is made available only in proc1...
        # HLCor_body can also compute "symSVD" attribute and add in each proc environment
        # If we put it in .preprocess (1) it will be computed in each envir, 
        #                             (2) we should condition on For=corrHLfit/fitme to avoid some useless computations
        ##                            (3) It cannot be in .assign_cov_matrices__from_covStruct bc this bit requires sparse_precision value
        decomp <- .provide_AR_factorization(corr_info$adjMatrices[[it]], proc1$sparsePrecisionBOOL, corr_type)
        if (corr_type=="SAR_WWt") attr(corr_info$adjMatrices[[it]],"UDU.") <- decomp
        if (corr_type=="adjacency") attr(corr_info$adjMatrices[[it]],"symSVD") <- decomp 
        rhorange <- sort(1/range(decomp$d)) ## keeping in mind that the bounds can be <>0
        if(verbose["SEM"]) cat(paste("Feasible rho range: ",paste(signif(rhorange,6),collapse=" -- "),"\n"))
        ## added 11/2016:
        if ( ! is.null(lower$rho)) rhorange[1L] <- max(rhorange[1L],lower$rho)
        if ( ! is.null(upper$rho)) rhorange[2L] <- min(rhorange[2L],upper$rho)
      }
    }
  }
  #
  ##### init.optim$phi/lambda will affect calc_inits -> calc_inits_dispPars.
  # outer estim seems useful when we can suppress all inner estim (thus the hatval calculations). 
  # ./. Therefore, we need to identify all cases where phi is fixed, 
  # ./. or can be outer optimized jointly with lambda:  
  phimodel <- proc1$models[['phi']]
  #if (phimodel=="phiGLM") {message("'fitme' not optimized for models with structured dispersion.")} ## FR->FR test this later"
  ranCoefs_blob <- proc1$ranCoefs_blob
  ## trying to guess all cases where optimization is useful. But FIXME: create all init and decide afterwardsS
  if ( (is_MixedM <- ( ! is.null(ranCoefs_blob) )) && (
    (var_ranCoefs <- (ranCoefs_blob$isRandomSlope & ! ranCoefs_blob$is_set)) ||
    phimodel == "phiScal" || # lambda + phi
    length(var_ranCoefs)>1L || # +s lambda
    length(corr_types[ ! is.na(corr_types)]) # lambda + corr pars
  )
  ) { ## prefer outer optimization if var_ranCoefs or not resid.model
    # (1) set (or not) outer optimization for phi: 
    phi.Fix <- proc1$phi.Fix
    ## All cases where phi is fixed, or can be outer optimized jointly with lambda:  
    not_inner_phi <- (! is.null(phi.Fix) ||
      ! is.null(init.optim$phi) || ## but its dubious that phiGLM is handled 
      phimodel == "phiScal" )
    if (not_inner_phi) {
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
    }
    # (2) set outer optimization for lambda and ranCoefs (handling incomplete ranFix$lambda vectors)
    ## tests show it is very inefficient to use outer optim on lambda (at least) when phi mus be inner optimized
    if (is_MixedM) { 
      lFix <- proc1$lambda.Fix ## only for 'simple" ranefs with Xi_cols=1
      nranterms <- sum(var_ranCoefs | is.na(lFix))
      which_simplelambda <- which(is.na(lFix) & ! ranCoefs_blob$isRandomSlope) ## exclude random slope whether set or not
      if (any(var_ranCoefs) && is.null(init.optim$ranCoefs)) {
        which_var_ranCoefs <- which(var_ranCoefs)
        Xi_cols <- attr(proc1$ZAlist,'Xi_cols')
        ranCoefs <- vector("list",length(which_var_ranCoefs))
        it <- 0L
        for (rt in which_var_ranCoefs) {
          it <- it+1L
          Xi_ncol <- Xi_cols[rt]
          rc <- rep(0,Xi_ncol*(Xi_ncol+1L)/2L)
          rc[cumsum(seq(Xi_ncol))] <- 0.1/(Xi_ncol*nranterms)
          ranCoefs[[it]] <- rc
        }
        names(ranCoefs) <- which_var_ranCoefs
        init.optim$ranCoefs <- ranCoefs
      }
      if (is.null(init.optim$lambda) && not_inner_phi ) {
        init.optim$lambda <- rep(0.1/nranterms,length(which_simplelambda))
        names(init.optim$lambda) <- which_simplelambda
      } ## complement to lambda.Fix
      # Low initial values  of lambda should be used to avoid spuriously high logLik from Laplace for high lambda  
      ## .../... unless Laplace is exact... (fixme)
    }
  }
  family <- proc1$family
  if (family$family=="COMPoisson") {
    checknu <- try(environment(family$aic)$nu,silent=TRUE)
    if (inherits(checknu,"try-error") && is.null(init.optim$COMP_nu)) init.optim$COMP_nu <- 1 ## FR->FR FIXME what if not incuded in the range ?
  } else if (family$family == "negbin") {
    checktheta <- try(environment(family$aic)$shape,silent=TRUE)
    if (inherits(checktheta,"try-error") && is.null(init.optim$NB_shape)) init.optim$NB_shape <- 1 ## FR->FR FIXME idem ...
  }
  user.lower <- lower; user.upper <- upper ## keep user input 
  # a few 'global' vars here:
  calc_inits_arglist <- list(init.optim=init.optim,init.HLfit=init.HLfit,ranFix=fixed,
                             corr_types=corr_types,optim.scale=optim.scale,
                  user.lower=user.lower,user.upper=user.upper,
                  control.dist=control.dist,For="fitme")
  ## adjustments of calc_inits_arglist:
  range_info_blob <- NULL
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if (! is.na(corr_type)) {
      if (corr_type == "Matern") {
        rho.size <- .check_conflict_init_fixed(fixed,init,
                                               "given as element of both 'fixed' and 'init'. Check call.")
        range_info_blob <- .calc_range_info(rho.size, processed, it, control.dist) 
        RHOMAX <- 1.1*30*range_info_blob$nbUnique/range_info_blob$maxrange## matches init$rho in calc_inits() ## not yet spaMM 3.0 
        control.dist$rho.mapping <- range_info_blob$rho_mapping
        ## not yet spaMM 3.0 (bug expected here if several Matern):
        calc_inits_arglist <-c(calc_inits_arglist,list(maxrange=range_info_blob$maxrange,RHOMAX=RHOMAX,NUMAX=NUMAX)) 
      }
      ## not yet spaMM 3.0 (bug expected here if several adj ):
      if (corr_type %in% c("SAR_WWt","adjacency") 
          &&  is.null(.getPar(fixed,"rho"))
          && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho$rho NULL or NA) 
      ) calc_inits_arglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
    }
  }
  inits <- do.call(".calc_inits",calc_inits_arglist)
  #
  init <- inits$`init` ## keeps all init values, all in untransformed scale
  init.optim <- inits$`init.optim` ## subset of all estimands, as name implies, and in transformed scale
  init.HLfit <- inits$`init.HLfit` ## subset as name implies 
  #
  if ("lambda" %in% c(names(user.lower),names(user.lower)) 
      && is.null(init$lambda)) {
    stop("'lambda' in 'lower' or 'upper' has no effect if absent from 'init'.")
  }
  ################
  LUarglist <- list(canon.init=init,
                    init.optim=init.optim,
                    user.lower=user.lower,user.upper=user.upper,
                    corr_types=corr_types,nbUnique=range_info_blob$nbUnique,
                    ranFix=fixed,control.dist=control.dist,
                    optim.scale=optim.scale, RHOMAX=RHOMAX,NUMAX=NUMAX,
                    rhorange=calc_inits_arglist$rhorange)
  LowUp <- do.call(".makeLowerUpper",LUarglist)
  ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  #
  processedHL1 <- proc1$HL[1] ## there's also HLmethod in processed<[[]]>$callargs
  needHLCor_specific_args <- (length(setdiff(names(lower),c("trPhi","trLambda","COMP_nu","trNB_shape","trRanCoefs")))>0L
                              || length(intersect(corr_types,c("Matern","adjacency","AR1","corrMatrix")))  
                              )
  if (spaMM.getOption("wDEVEL2")) {
    parlist <- .merge_parlist(,new=fixed,types="fix")
    parlist <- .merge_parlist(parlist,new=init.HLfit,types="var")
    attr(fixed,"parlist") <- parlist
  } 
  if (needHLCor_specific_args) {
    HLcallfn.obj <- "HLCor.obj" 
    HLcallfn <- "HLCor"
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
  .assignWrapper(anyHLCor_obj_args$processed,
                   paste("return_only <- \"",anyHLCor_obj_args$processed$objective,"APHLs\"",sep=""))
  initvec <- unlist(init.optim)
  if (length(initvec)) {
    if (identical(verbose["getCall"][[1L]],TRUE)) { ## toget an optim call with its initial value. Then HLcallfn is called and its call returned.
      refitPars <- init.optim
      refit_info <- FALSE
    } else {
      use_SEM <- (!is.null(processedHL1) && processedHL1=="SEM")
      if (use_SEM) {
        if (is.null(proc1$SEMargs$control_pmvnorm$maxpts)) {
          if (length(LowUp$lower)) {
            .assignWrapper(processed,"SEMargs$control_pmvnorm$maxpts <- quote(250L*nobs)") 
          } ## else default visible in SEMbetalambda
        }
        stop("reimplement iterateSEMSmooth() in fitme() later")
        ## and then we need to put back the code for logL from smoothing at the end of this function
        ## beware of SEM case with length(lower)=0
        refit_info <- control[["refit"]] ## may be a NULL/ NA/ boolean or a list of booleans 
        if (is.null(refit_info) || is.na(refit_info)) refit_info <- FALSE ## alternatives are TRUE or an explicit list, but no longer NULL
      } else {
        refitPars <- optPars <- .new_locoptim(init.optim, 
                                              LowUp, 
                                              control, objfn_locoptim=.objfn_locoptim, 
                                              HLcallfn.obj=HLcallfn.obj, anyHLCor_obj_args=anyHLCor_obj_args, 
                                              user_init_optim=user_init_optim)
        refit_info <- attr(optPars,"refit_info") ## 'derives' from control[["refit"]]
      }
    }
    RanFixOrVar <- fixed 
    # refit_info is list if so provided by user, else typically boolean. An input NA should have been converted to something else (not documented).
    if ( is.list(refit_info)) {refit_phi <- refit_info$phi} else refit_phi <- refit_info ## result may be NULL in list case
    if (identical(refit_phi,TRUE) && ! is.null(refitPars$trPhi)) {
      HLCor.args$init.HLfit$phi <- .dispInv(refitPars$trPhi)
      RanFixOrVar$trPhi <- refitPars$trPhi <- NULL
    }
    if ( is.list(refit_info)) {refit_lambda <- refit_info$lambda} else refit_lambda <- refit_info
    if (identical(refit_lambda,TRUE) && ! is.null(refitPars$trLambda)) {
      HLCor.args$init.HLfit$lambda <- .dispInv(refitPars$trLambda)
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
  .assignWrapper(HLCor.args$processed,"return_only <- NULL") 
  .assignWrapper(HLCor.args$processed,"verbose['warn'] <- TRUE") ## important!
  hlcor <- do.call(HLcallfn,HLCor.args) ## recomputation post optimization, or only computation if length(initvec)=0
  if (is.call(hlcor)) {
    ## then do.call(HLcallfn,HLCor.args) has retuned the call, not the fit. 
    ## see def of get_HLCorcall() for further explanation
    return(hlcor) ## HLCorcall
  }
  # hlcor may have received attr(.,"info.uniqueGeo") from HLCor_body.
  if (length(initvec)) {
    attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, 
                                    objective=anyHLCor_obj_args$processed$objective) ## processed was erased for safety
  }
  ## substantial effect on object size! :
  lsv <- c("lsv",ls())
  rm(list=setdiff(lsv,"hlcor")) 
  ##
  return(hlcor)
}