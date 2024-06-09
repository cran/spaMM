## the following must match the 'unique' method in .ULI as explained there
.calcUniqueGeo <- function(data) { # unique() on a character representation of the data
  redondGeo <- .pasteCols(t(data), collapse=" ") # apply(data,1,paste,collapse=" ") ## creates character string
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
  types <- attr(covStruct,"types") ## 1st way of specifying types; immediately converted to names of list 
                                   ## ** so the attribute can be neglected in further code.**
  if (is.null(types)) {
    types <- names(covStruct) ## 2nd way of specifying types
  } else names(covStruct) <- types ## repeated names possible
  known_types <- c("adjMatrix","corrMatrix",
                   "precision", # => not one of the built_in ranef types
                   "SAR_WWt","distMatrix", "IMRF","corrFamily") 
  checktypes <- setdiff(types,c(known_types,"", paste(seq_along(covStruct)))) ## "" for unhandled ranefs
  if (length(checktypes)) stop(paste("Unhandled name(s)/type(s)",
                                     paste0("'",checktypes,"'",collapse=", "),"in 'covStruct'."))
  for (lit in which(types=="precision")) {
    if (is.list(covStruct[[lit]])) { # NOT always the case. In the tests, fitme(cases ~ I(prop.ag/10) + corrMatrix(1 | gridcode) .. , covStruct = list(precision = precmat) 
      covStruct[[lit]]$matrix <- forceSymmetric(covStruct[[lit]]$matrix)
    } else covStruct[[lit]] <- forceSymmetric(covStruct[[lit]])
  }
  if ( ! is.null(AMatrices <- attr(covStruct,'AMatrices')) && ! is.list(AMatrices)) {stop("attr(covStruct,'AMatrices') must be a list.")} 
  return(covStruct)
}

.check_corrMatrix <- function(corrMatrix, element) {
  if (is.list(corrMatrix)) {
    dim_corrMatrix <- dim(corrMatrix[[element]])
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

# Creating an ad hoc list may not be memory efficient but is quite convenient for prediction.
# The .get_old_info_uniqueGeo() extractor should be used to access this info.
.reformat_info_uniqueGeo <- function(processed, control.dist, corr_types=processed$corr_info$corr_types) {
  info_uniqueGeo <- msd_arglist <- list() # list with char_rd-indexed elements
  
  # for old builtin types:
  is_uniqueGeo_needed <- ( (! is.na(corr_types)) & 
                             (corr_types=="Matern" | corr_types=="Cauchy" | 
                                corr_types=="AR1" | # probably redundant with next condition
                                processed$corr_info$levels_types=="time_series" | # includes ARp
                                corr_types=="IMRF" ## IMRF: mapMM expects it.
                             )
  ) 
  for (rd in which(is_uniqueGeo_needed)) {
    char_rd <- as.character(rd)
    geo_envir <- .get_geo_info(processed, which_ranef=rd, which="uniqueGeo", 
                               dist_method_rd=control.dist[[char_rd]]$dist.method) 
    info_uniqueGeo[[char_rd]] <- geo_envir$uniqueGeo 
  }
  
  # for newer corr families:
  corr_families <- processed$corr_info$corr_families
  for (rd in which(corr_types=="corrFamily")) {
    char_rd <- as.character(rd)
    if ( ! inherits(corr_families[[rd]]$calc_corr_from_dist, "stub")) {
      geo_envir <- .get_geo_info(processed, which_ranef=rd, which="uniqueGeo", 
                                 dist_method_rd=control.dist[[char_rd]]$dist.method) 
      info_uniqueGeo[[char_rd]] <- geo_envir$uniqueGeo 
    }
  }
  info_uniqueGeo
}

HLCor_body <- function(processed, ## single environment
                  #ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  control.dist=list(), # possibly distinct from processed info bc corrHLfit_body/fitme_body may have modified it => 
                                       # this one is used, not the one from processed
                  fixed=NULL,
                  ...) { ## dots for HLfit
  dotlist <- list(...)
  spatial_terms <- attr(processed$ZAlist,"exp_spatial_terms")
  corr_types <- processed$corr_info$corr_types
  ## convert back ranPars to canonical scale:
  fixed <- .reformat_corrPars(fixed,corr_families=processed$corr_info$corr_families) 
  fixed <- .canonizeRanPars(ranPars=fixed, corr_info=processed$corr_info, rC_transf=.spaMM.data$options$rC_transf) ## with init.HLfit as attribute # expands hyper params.
  ########################################################################################################################
  # * assigns geo_envir <- .get_geo_info(...)
  # * modifies processed$AUGI0_ZX$envir by .init_AUGI0_ZX_envir_spprec_info(...) 
  # * computes processed$AUGI0_ZX$envir$LMatrices except for ranCoefs (the latter being filled in HLfit_body)
  .assign_geoinfo_and_LMatrices_but_ranCoefs(processed, corr_types, spatial_terms, ranPars=fixed, control.dist, 
                                argsfordesignL=dotlist[intersect(names(dotlist),names(formals(mat_sqrt)))] )
  if (any(vec_normIMRF <- processed$AUGI0_ZX$vec_normIMRF ) ) {
    processed$ZAlist <- .normalize_IMRF(processed, 
                                        vec_normIMRF=vec_normIMRF,
                                        #ranges=processed$hyper_info$ranges,
                                        strucList=processed$AUGI0_ZX$envir$LMatrices)
    .assign_ZAfix(processed)
  }
  ########################################################################################################################
  ###
  if ( (! is.null(processed$return_only)) && processed$augZXy_cond) {
    hlfit <- do.call(.spaMM.data$options$augZXy_fitfn,list(processed=processed, fixed=fixed))
    if (FALSE) { # check consistency of aug_ZXy and non-aug_ZXy procedures
      # This code works: reused 01/2022 to check corrFamily model. 
      ## this test has previously interfered with the results (fitme3, fitme6 tests), presumably bc of the inner attribute (now added here, ignored during the test):
      processed$augZXy_cond <- structure(FALSE, inner=attr(processed$augZXy_cond,"inner"))
      HLFormals <- names(formals(HLfit)) 
      good_dotnames <- intersect(names(dotlist),HLFormals)
      if (length(good_dotnames)) {
        HL.info <- dotlist[good_dotnames]  ## including init.HLfit: possibly modified from processed$init_HLfit by <corrfitme>_body 
      } else HL.info <- list()
      ## all printing in HLfit is suppressed by default
      HL.info$processed <- processed
      HL.info$init.HLfit <- .modify_list(HL.info$init.HLfit, attr(fixed,"init.HLfit")) 
      locrp <- fixed
      attr(locrp,"init.HLfit") <- NULL
      locrp$phi <- hlfit$APHLs$phi_est # Imporant: if we don't fix that, 'vanilla' and 'hlfit' will be inconsistent. This is expected: 
      # "vanilla" estimates phi for lambda fixed, while "hlfit" estimate a (phi+lambda) parameter for fixed phi/lambda ratio (a given lambda being a lambda factor)
      locrp$lambda <- locrp$lambda * hlfit$APHLs$phi_est
      HL.info$ranFix <- locrp
      oldr <- processed$return_only
      processed$return_only <- FALSE # not essential to catch problems, but should help to diagnose problems
      vanilla <- do.call("HLfit",HL.info) 
      #print(vanilla$APHLs$p_bv-hlfit$APHLs$p_bv)
      if (is.null(hlfit$APHLs$p_bv)) { # for REML fits, check p_bv
        if (abs(vanilla$APHLs$p_v-hlfit$APHLs$p_v)>1e-4) browser("Problem in check of aug_ZXy procedure") # browser() only in devel code
      } else if ((vanilla$APHLs$p_bv-hlfit$APHLs$p_bv)>1e-4) browser("Problem in check of aug_ZXy procedure") # browser() only in devel code
      #if (abs(vanilla$APHLs$p_bv-hlfit$APHLs$p_bv)>1e-4) browser() # 
      # hmmm aug_ZXy's p_bv is suspicious close to vanilla's p_v
      # debug(.calc_APHLs_by_augZXy_or_sXaug)
      
      #  clean debugging mess before leaving
      processed$return_only <- oldr
      processed$augZXy_cond <- structure(TRUE, inner=attr(processed$augZXy_cond,"inner"))
    }
  } else {
    HLFormals <- names(formals(HLfit)) 
    good_dotnames <- intersect(names(dotlist),HLFormals)
    if (length(good_dotnames)) {
      HL.info <- dotlist[good_dotnames]  ## including init.HLfit: possibly modified from processed$init_HLfit by <corrfitme>_body 
    } else HL.info <- list()
    ## all printing in HLfit is suppressed by default
    HL.info$processed <- processed
    HL.info$init.HLfit <- .modify_list(HL.info$init.HLfit, attr(fixed,"init.HLfit")) 
    attr(fixed,"init.HLfit") <- NULL
    HL.info$fixed <- fixed
    HL.info$ranFix <- NULL
    hlfit <- do.call("HLfit",HL.info) 
  }  
  ## Here there was debug code that saved HL.info in case of error; before 1.8.5
  if ( ! is.null(processed$return_only)) {
    return(hlfit)    ########################   R E T U R N   a list with $APHLs, with class "list"
  } # else class(hlfit) <- c("HLCor", class(hlfit)) seems totally obsolete; object has class HLfit.
  #### Infos in the final fit object: 
  hlfit$spatial_terms <- spatial_terms
  hlfit$ranef_info$info_oldUniqueGeo <- .reformat_info_uniqueGeo(processed, control.dist, corr_types) # list; previously stored as attr(hlfit,"info.uniqueGeo")
  #
  hlfit$call <- "$call removed by HLCor_body. Use getCall() to extract the call from the object." ## instead of the $call with evaluated arguments
  return(hlfit) ## 
}

# Called by HLCor.obj():
.merge_fixed <- function(fixed, ranefParsList, skeleton, processed, HLCor.call) {
  fixed <- .modify_list(fixed, ranefParsList) # merges variable and fixed params, 
  # but result may be messy (lambda+trLambda... corrFamily parameters in wrong order...) it will be
  # HLfit|HLCor_body -> .canonizeRanPars 's task to put this in order.
  rpType <- .modify_list(attr(fixed,"type"),attr(skeleton,"type"))
  moreargs <- attr(skeleton,"moreargs")
  fixed <- .expand_hyper(fixed, processed$hyper_info,moreargs=moreargs) ## input ranPars contains both unconstrained ranPars and $hyper
  # => failing to expand leads to unconstrained optimization
  # removed 'ranPars$resid' code here [ v3.5.52
  attr(fixed,"type") <- rpType
  attr(fixed,"moreargs") <- moreargs
  fixed
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
  HLCor.call <- mc[c(1L,which(names(mc) %in% HLnames))] ## keep the call structure
  ranefParsList <- relist(ranefParsVec,skeleton) # converts back a vector of variable parameters to a structured list of variable parameters
  print_phiHGLM_info <- ( ! is.null(processed$residProcessed) && processed$verbose["phifit"]) 
  if (print_phiHGLM_info) {
    # set a 'prefix' for the line to be printed for each iteration of the phi fit when outer optimization is used for the mean response. 
    # In that case a *distinct line* of the form HLCor for <outer opt pars>: phi fit's iter=<say up to 6>, .phi[1]=... 
    # is written for each call of the outer objfn (=> multi-line output).
    # Currently there is no such 'prefix' for mv (_F I X M E_)
    # That would require checking processed$residProcesseds (with -'s') and some further effort.
    urP <- unlist(.canonizeRanPars(ranefParsList, corr_info=processed$corr_info,checkComplete=FALSE, 
                                   rC_transf=.spaMM.data$options$rC_transf))
    processed$port_env$prefix <- paste0("HLCor for ", paste(signif(urP,6), collapse=" "), ": ")
  } 
  if ( ! is.null(processed$X_off_fn)) { # beta outer-optimisation
    if ( ! is.null(trBeta <- ranefParsList$trBeta)) { # outer beta
      ranefParsList$trBeta <- NULL
      HLCor.call$etaFix$beta <- .spaMM.data$options$.betaInv(trBeta)
    } else if ( ! is.null(beta <- ranefParsList$beta)) { # outer beta
      ranefParsList$beta <- NULL
      HLCor.call$etaFix$beta <- beta
    }
  }
  HLCor.call$fixed <- .merge_fixed(HLCor.call$fixed, ranefParsList, skeleton, processed, HLCor.call)
  # 'fixed' may have $trLambda (from notlambda) for what is optimized,
  #              and $lambda (from ranPars$lambda) for what was fixed in the whole outer fit  
  HLCor.call[[1L]] <- processed$HLCor
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


