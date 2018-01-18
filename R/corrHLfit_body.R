corrHLfit_body <- function(processed, ## possibly a list of environments
                           init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, optimizer, <optimizer controls>
                      nb_cores=NULL,
                      ... ## cf dotnames processing below
) { 

  

  #########################################################
  #
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
  sparse_precision <- processed$sparsePrecisionBOOL
  # test-Nugget -> different p_v / p_bv compromise (!)
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
        if(verbose["SEM"])  cat(paste("Feasible rho range: ",paste(signif(rhorange,6),collapse=" -- "),"\n"))
        ## added 11/2016:
        if ( ! is.null(lower$rho)) rhorange[1L] <- max(rhorange[1L],lower$rho)
        if ( ! is.null(upper$rho)) rhorange[2L] <- min(rhorange[2L],upper$rho)
      }
    }
  }
  family <- proc1$family
  if (family$family=="COMPoisson") {
    checknu <- try(environment(family$aic)$nu,silent=TRUE)
    if (inherits(checknu,"try-error") && is.null(init.optim$COMP_nu)) init.optim$COMP_nu <- 1
  } else if (family$family == "negbin") {
    checktheta <- try(environment(family$aic)$shape,silent=TRUE)
    if (inherits(checktheta,"try-error") && is.null(init.optim$NB_shape)) init.optim$NB_shape <- 1
  }
  #
  user.lower <- lower; user.upper <- upper ## keep user input 
  calc_inits_arglist <- list(init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix,
                             user.lower=user.lower,user.upper=user.upper,
                             corr_types=corr_types,optim.scale=optim.scale,
                  
                  control.dist=control.dist,For="corrHLfit") ## NOTfixme difference fitme/corrHLfit
  ## adjustments of calc_inits_arglist:
  range_info_blob <- NULL
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if (! is.na(corr_type)) {
      if (corr_type %in% c("Matern")) {
        rho.size <- .check_conflict_init_fixed(ranFix,init.corrHLfit, ## NOTfixme: difference bw fitme et corrHLfit 
                                               "given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")
        range_info_blob <- .calc_range_info(rho.size, processed, it, control.dist) 
        control.dist$rho.mapping <- range_info_blob$rho_mapping
        RHOMAX <- 1.1*30*range_info_blob$nbUnique/range_info_blob$maxrange ## matches init$rho in calc_inits() ## not yet spaMM 3.0 
        ## not yet spaMM 3.0 (bug expected here if several Matern):
        calc_inits_arglist <-c(calc_inits_arglist,list(maxrange=range_info_blob$maxrange,RHOMAX=RHOMAX,NUMAX=NUMAX)) 
      }
      ## not yet spaMM 3.0 (bug expected here if several adj ):
      if (corr_type %in% c("SAR_WWt","adjacency")
          &&  is.null(.getPar(ranFix,"rho")) ## NOTfixme difference fitme/corrHLfit
          && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho NULL or NA) 
      ) calc_inits_arglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
    }
  }
  inits <- do.call(".calc_inits",calc_inits_arglist)
  #
  init <- inits$`init` ## keeps all init values, all in untransformed scale
  #  code correct mais opaque pour l'utilisateur: init.optim va être ignore en cas de 1D optimization... (OK in fitme)
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
    ## corrHLfit-specific logic: if lambda absent from init.corrHLfit it is not outer optimized and then 'as the message says'.
  }
  ################
  LUarglist <- list(canon.init=init, ## canonical scale
                    init.optim=init.optim, ## transformed scale, used only to initialize lower and upper 
                    user.lower=user.lower,user.upper=user.upper,
                    corr_types=corr_types,nbUnique=range_info_blob$nbUnique,
                    ranFix=ranFix,control.dist=control.dist,
                    optim.scale=optim.scale, RHOMAX=RHOMAX,NUMAX=NUMAX,
                    rhorange=calc_inits_arglist$rhorange)
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
  processedHL1 <- proc1$HL[1] ## there's also HLmethod in processed<[[]]>$callargs
  if (!is.null(processedHL1) && processedHL1=="SEM" && length(lower)>0) {
    optimMethod <- "iterateSEMSmooth"
    if (is.null(proc1$SEMargs$control_pmvnorm$maxpts)) {
      if (length(LowUp$lower)) {
        .assignWrapper(processed,"SEMargs$control_pmvnorm$maxpts <- quote(250L*nobs)") 
      } ## else default visible in SEMbetalambda
    }
  } else optimMethod <- ".new_locoptim"
  HLCor.args$processed <- processed
  ## 
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  anyHLCor_obj_args$skeleton <- structure(init.optim,RHOMAX=RHOMAX,NUMAX=NUMAX) ## logscale, only used by HLCor.obj
  attr(anyHLCor_obj_args$skeleton,"type") <- list() ## declares a list of typeS of elemnts of skeleton
  attr(anyHLCor_obj_args$skeleton,"type")[names(init.optim)] <- "fix" # fixed with the HLCor call 
  .assignWrapper(anyHLCor_obj_args$processed,
                   paste("return_only <- \"",anyHLCor_obj_args$processed$objective,"APHLs\"",sep=""))
  initvec <- unlist(init.optim)
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
    if (identical(verbose["getCall"][[1L]],TRUE)) {
      optPars <- init.optim
    } else {
      optPars <- .new_locoptim(init.optim, LowUp=LowUp,anyHLCor_obj_args=anyHLCor_obj_args,
                               objfn_locoptim=.objfn_locoptim,
                               control=control.corrHLfit)
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
  .assignWrapper(HLCor.args$processed,"return_only <- NULL") 
  .assignWrapper(HLCor.args$processed,"verbose['warn'] <- TRUE") ## important!
  hlcor <- do.call("HLCor",HLCor.args) ## recomputation post optimization (or only computation, if length(lower)=0)
  #
  if (is.call(hlcor)) { return(hlcor[]) } ## HLCorcall
  #
  # hlcor should have received attr(.,"info.uniqueGeo") from HLCor_body.
  attr(hlcor,"optimInfo") <- list(LUarglist=LUarglist, optim.pars=optPars, objective=HLCor.args$processed$objective)
  if ( ! is.null(optPars)) {
    locoptr <- attr(optPars,"optr")
    if (attr(optPars,"method")=="nloptr") {
      if (locoptr$status<0L) hlcor$warnings$optimMessage <- paste("nloptr() message: ",
                                                                  locoptr$message," (status=",locoptr$status,")",sep="")
    } else if ( attr(optPars,"method")=="optim" ) {
      if (locoptr$convergence) hlcor$warnings$optimMessage <- paste("optim() message: ",locoptr$message,
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

    
