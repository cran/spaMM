## the following must match the 'unique' method in .ULI as explained there
.calcUniqueGeo <- function(data) { # unique() on a character representation of the data
  redondGeo <- apply(data,1,paste,collapse=" ") ## creates character string
  dfforunique <- cbind(data,redondGeo) ## associates rownames of data to redondGeo
  uniqueGeo <- unique(dfforunique[,ncol(dfforunique),drop=FALSE]) ## keeps rownames of first instances
  uniqueGeo <- data[rownames(uniqueGeo),,drop=FALSE] ## uses rownames, 'unique' numeric values based on character representations 
  return(uniqueGeo)
}


.extract_check_coords <- function(spatial_term,datanames, check=TRUE) {
  if ( ! is.null(spatial_term)) {
    bars <- spatial_term[[2]] 
    coordinates <- .DEPARSE(bars[[3]]) ## "x + y"
    coordinates <-  strsplit(coordinates," ")[[1]]
    coordinates <- setdiff(coordinates,c("+","%in%",":","/","")) ## "" for hidden linebreaks (?)
  } else {
    stop("Spatial term is NULL.")
  }
  if (check) {
    coordcols <- which(datanames %in% coordinates)
    if ( length(coordcols) != length(coordinates) ) {
      stop("variables 'coordinates' not all found in the 'data'")
    }
  }
  return(coordinates) ## should be ordered as bars[[3]] (important for predict)
}

.extract_check_coords_within <- function(spatial_term) {
  bars <- spatial_term[[2]] 
  coordinates <- .DEPARSE(bars[[3]]) ## "x + y"
  coordinates <-  strsplit(coordinates," ")[[1]]
  if (length(grep("/",coordinates))) {
    stop(paste("'/' not yet handled in",spatial_term))
  } else if (length(grep_in <- grep("%in%|:",coordinates))>1L) {
    stop(paste("multiple nesting not yet handled in",spatial_term))
  } else if (length(grep_in <- grep("%in%|:",coordinates))) {
    coordinates <- coordinates[1L:(min(grep_in)-1L)]
  }
  coordinates <- setdiff(coordinates,c("+","")) ## "" for hidden linebreaks (?)
  return(coordinates) ## should be ordered as bars[[3]] (important for predict)
}



.preprocess_covStruct <- function(covStruct) {
  if ( ! inherits(covStruct,"list")) stop("covStruct must inherit from class 'list'.")
  types <- attr(covStruct,'types') ## 1st way of specifying types
  if (is.null(types)) {
    types <- names(covStruct) ## 2nd way of specifying types
  } else names(covStruct) <- types ## repeated names possible
  known_types <- c("adjMatrix","corrMatrix","precision","SAR_WWt","distMatrix", "IMRF") 
  checktypes <- setdiff(types,c(known_types,"", paste(seq_along(covStruct)))) ## "" for unhandled ranefs
  if (length(checktypes)) stop(paste("Unhandled name(s)/type(s)",
                                     paste0("'",checktypes,"'",collapse=", "),"in 'covStruct'."))
  for (lit in seq_along(covStruct)) {
    if (types[[lit]]=="precision") {
      if (is.list(covStruct[[lit]])) { # always the case ?
        covStruct[[lit]]$matrix <- forceSymmetric(covStruct[[lit]]$matrix)
      } else covStruct[[lit]] <- forceSymmetric(covStruct[[lit]])
    } 
  }
  if ( ! is.null(AMatrices <- attr(covStruct,'AMatrices')) && ! is.list(AMatrices)) {stop("attr(covStruct,'AMatrices') must be a list.")} 
  return(covStruct)
}

.check_corrMatrix <- function(corrMatrix) {
  if (is.list(corrMatrix)) {
    dim_corrMatrix <- dim(corrMatrix[[1]])
  } else dim_corrMatrix <- dim(corrMatrix)
  if (dim_corrMatrix[1L]!=dim_corrMatrix[2L])  stop("corrMatrix is not square") 
}


.calc_AR1_t_chol_Q_block <- function(n_u, ARphi) {
  # denom <- sqrt(1-ARphi^2)
  # t_chol_Q <- .symDiagonal(x=c(rep(1/denom,n_u-1L),1)) 
  # if (n_u>1L) diag(t_chol_Q[,-1,drop=FALSE]) <- -ARphi/denom 
  # return(t_chol_Q) # equivalent to nlme's AR1_fact() in corStruct.c
  seqn <- seq(n_u)
  denom <- sqrt(1-ARphi^2)
  list(i=c(seqn,seqn[-n_u]),
       j=c(seqn,seqn[-1L]),
       x=c(rep(1/denom, n_u-1L),1,rep(-ARphi/denom, n_u-1L))
  )
}


HLCor_body <- function(processed, ## single environment
                  ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  control.dist=list(), # possibly distinct from processed info bc corrHLfit_body/fitme_body may have modified it => 
                                       # this one is used, not the one from processed
                  ...) { ## dots for HLfit
  dotlist <- list(...)
  spatial_terms <- attr(processed$ZAlist,"exp_spatial_terms")
  corr_types <- processed$corr_info$corr_types
  ## convert back ranPars to canonical scale:
  ranPars <- .reformat_corrPars(ranPars,corr_families=processed$corr_info$corr_families) 
  ranPars <- .canonizeRanPars(ranPars=ranPars, corr_info=processed$corr_info, rC_transf=.spaMM.data$options$rC_transf) ## with init.HLfit as attribute # expands hyper params.
  ########################################################################################################################
  # * assigns geo_envir <- .get_geo_info(...)
  # * modifies processed$AUGI0_ZX$envir by .init_AUGI0_ZX_envir_spprec_info(...) 
  # * computes processed$AUGI0_ZX$envir$LMatrices except for ranCoefs (the latter being filled in HLfit_body)
  .assign_geoinfo_and_LMatrices_but_ranCoefs(processed, corr_types, spatial_terms, ranPars, control.dist, 
                                argsfordesignL=dotlist[intersect(names(dotlist),names(formals(mat_sqrt)))] )
  if (any(vec_normIMRF <- processed$AUGI0_ZX$vec_normIMRF ) ) {
    processed$ZAlist <- .normalize_IMRF(processed$ZAlist, 
                                        vec_normIMRF=vec_normIMRF,
                                        #ranges=processed$hyper_info$ranges,
                                        strucList=processed$AUGI0_ZX$envir$LMatrices)
    .assign_ZAfix(processed)
  }
  ########################################################################################################################
  ###
  if ( (! is.null(processed$return_only)) && processed$augZXy_cond) {
    hlfit <- do.call(.spaMM.data$options$augZXy_fitfn,list(processed=processed, ranFix=ranPars))
    if (FALSE) { ## this test interfered with the results (fitme3, fitme6 tests), presumably bc of the inner attribute (now added here, ignored during the test):
      processed$augZXy_cond <- structure(FALSE, inner=attr(processed$augZXy_cond,"inner"))
      HLFormals <- names(formals(HLfit)) 
      good_dotnames <- intersect(names(dotlist),HLFormals)
      if (length(good_dotnames)) {
        HL.info <- dotlist[good_dotnames]  ## including init.HLfit: possibly modified from processed$init_HLfit by <corrfitme>_body 
      } else HL.info <- list()
      ## all printing in HLfit is suppressed by default
      HL.info$processed <- processed
      HL.info$init.HLfit <- .modify_list(HL.info$init.HLfit, attr(ranPars,"init.HLfit")) 
      locrp <- ranPars
      attr(locrp,"init.HLfit") <- NULL
      locrp$phi <- hlfit$APHLs$phi_est
      locrp$lambda <- locrp$lambda * hlfit$APHLs$phi_est
      HL.info$ranFix <- locrp
      vanilla <- do.call("HLfit",HL.info) 
      processed$augZXy_cond <- structure(TRUE, inner=attr(processed$augZXy_cond,"inner"))
      #print(vanilla$APHLs$p_bv-hlfit$APHLs$p_bv)
      #if ((vanilla$APHLs$p_bv-hlfit$APHLs$p_bv)>1e-2) browser() # to catch in part. when phi-profiled p_bv is lower than non profiled !
      # hmmm aug_ZXy's p_bv is suspicious close to vanilla's p_v
      # debug(.calc_APHLs_by_augZXy_or_sXaug)
    }
  } else {
    HLFormals <- names(formals(HLfit)) 
    good_dotnames <- intersect(names(dotlist),HLFormals)
    if (length(good_dotnames)) {
      HL.info <- dotlist[good_dotnames]  ## including init.HLfit: possibly modified from processed$init_HLfit by <corrfitme>_body 
    } else HL.info <- list()
    ## all printing in HLfit is suppressed by default
    HL.info$processed <- processed
    HL.info$init.HLfit <- .modify_list(HL.info$init.HLfit, attr(ranPars,"init.HLfit")) 
    attr(ranPars,"init.HLfit") <- NULL
    HL.info$ranFix <- ranPars
    hlfit <- do.call("HLfit",HL.info) 
  }  
  ## Here there was debug code that saved HL.info in case of error; before 1.8.5
  if ( ! is.null(processed$return_only)) {
    return(hlfit)    ########################   R E T U R N   a list with $APHLs, with class "list"
  } # else class(hlfit) <- c("HLCor", class(hlfit)) seems totally obsolete; object has class HLfit.
  #### Infos in the final fit object: 
  hlfit$spatial_terms <- spatial_terms
  info_uniqueGeo <- msd_arglist <- list()
  is_uniqueGeo_needed <- ( (! is.na(corr_types)) & (corr_types=="Matern" | corr_types=="Cauchy" | corr_types=="AR1"
                           | corr_types=="IMRF")) ## IMRF: mapMM expects it.
  for (rd in which(is_uniqueGeo_needed)) {
    char_rd <- as.character(rd)
    geo_envir <- .get_geo_info(processed, which_ranef=rd, which="uniqueGeo", 
                               dist_method_rd=control.dist[[char_rd]]$dist.method) 
    info_uniqueGeo[[char_rd]] <- geo_envir$uniqueGeo 
  }
  attr(hlfit,"info.uniqueGeo") <- info_uniqueGeo
  #attr(hlfit,"control_dists") <- control.dist
  #
  hlfit$call <- "$call removed by HLCor_body. Use getCall() to extract the call from the object." ## instead of the $call with evaluated arguments
  return(hlfit) ## 
}


## wrapper for HLCor, suitable input and output for optimization
`HLCor.obj` <- function(ranefParsVec,skeleton,objective=processed$objective,processed,...) { ## name of first arg MUST differ from names in dotlist...
  mc <- match.call(expand.dots=TRUE) ## (1) expand.dots added 11/04/2014 for the multinomial... eval 
  
  if (is.null(processed)) { stop("Call to HLCor.obj() without a 'processed' argument is invalid") }

  if (  is.list(processed) )  { ## "multiple" processed list 
    ## RUN THIS LOOP and return
    fitlist <- lapply(seq_len(length(processed)), function(proc_it){
      locmc <- mc
      locmc[[1L]] <- get("HLCor.obj", asNamespace("spaMM"), inherits=FALSE) # as.name("HLCor.obj") ## replaces "f" !
      locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
      locmc$processed <- processed[[proc_it]] 
      eval(locmc)
    }) ## a pure list of HLCor objects
    resu <- sum(unlist(fitlist))
    return(resu)
  } else { ## there is one processed for a single data set 
    family <- processed$family
    #data <- processed$data
  }
  
  HLCor.formals <- names(formals(HLCor))
  names_formals_HLfit <- names(formals(HLfit))
  designL.formals <- names(formals(mat_sqrt))
  makescaled.formals <- names(formals(make_scaled_dist))
  HLnames <- (c(HLCor.formals,names_formals_HLfit,designL.formals,makescaled.formals))  ## cf parallel code in corrHLfit
  HLCor.call <- mc[c(1,which(names(mc) %in% HLnames))] ## keep the call structure
  ranefParsList <- relist(ranefParsVec,skeleton)
  print_phiHGLM_info <- ( ! is.null(processed$residProcessed) && processed$verbose["phifit"]) 
  if (print_phiHGLM_info) {
    # set a 'prefix' for the line to be printed for each iteration of the phi fit when outer optimization is used for the main response. 
    # In that case a *distinct line* of the form HLCor for <outer opt pars>: phi fit's iter=<say up to 6>, .phi[1]=... 
    # is written for each call of the outer objfn (=> multi-line output).
    # Currently there is no such 'prefix' for mv (_F I X M E_)
    # That would require checking processed$residProcesseds (with -'s') and some further effort.
    urP <- unlist(.canonizeRanPars(ranefParsList, corr_info=processed$corr_info,checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf))
    processed$port_env$prefix <- paste0("HLCor for ", paste(signif(urP,6), collapse=" "), ": ")
  } 
  ranPars <- .modify_list(HLCor.call$ranPars, ranefParsList)
  rpType <- .modify_list(attr(HLCor.call$ranPars,"type"),attr(skeleton,"type"))
  moreargs <- attr(skeleton,"moreargs")
  ranPars <- .expand_hyper(ranPars, processed$hyper_info,moreargs=moreargs) ## input ranPars contains both unconstrained ranPars and $hyper
  # => failing to expand leads to unconstrained optimization
  # removed 'ranPars$resid' code here [ v3.5.52
  HLCor.call$ranPars <- structure(ranPars, ## adds given values of the optimized variables 
                                  type=rpType, ## adds "fix"'s... somewhat confusing 
                                  moreargs=moreargs )
  # ranPars may have $trLambda (from notlambda) for what is optimized,
  #              and $lambda (from ranPars$lambda) for what was fixed in the whole outer fit  
  HLCor.call[[1L]] <- get("HLCor", asNamespace("spaMM"), inherits=FALSE) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
  #HLCor.call[[1L]] <- quote(spaMM::HLCor)
  hlfit <- eval(HLCor.call) ## returns fit or call 
  #
  if (is.call(hlfit)) {return(hlfit)} ## HLCorcall
  #
  aphls <- hlfit$APHLs
  resu <- aphls[[objective]]
  if (print_phiHGLM_info) cat(paste0(objective,"=",resu)) # verbose["phifit"]
  if (objective=="cAIC") resu <- - resu ## for minimization of cAIC (private & experimental)
  if (processed$augZXy_cond && resu>processed$augZXy_env$objective) {
    processed$augZXy_env$objective <- resu
    processed$augZXy_env$phi_est <- aphls[["phi_est"]]
  }
  return(resu) #
}


