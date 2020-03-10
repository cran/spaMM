## the following must match the'unique' method in .ULI as explained there
.calcUniqueGeo <- function(data) {
  redondGeo <- apply(data,1,paste,collapse=" ") ## creates character string
  dfforunique <- cbind(data,redondGeo) ## associates rownames of data to redondGeo
  uniqueGeo <- unique(dfforunique[,ncol(dfforunique),drop=FALSE]) ## keeps rownames of first instances
  uniqueGeo <- data[rownames(uniqueGeo),,drop=FALSE] ## uses rownames, 'unique' numeric values based on character representations 
  return(uniqueGeo)
}


.extract_check_coords <- function(spatial_term,datanames) {
  if ( ! is.null(spatial_term)) {
    bars <- spatial_term[[2]] 
    coordinates <- .DEPARSE(bars[[3]]) ## "x + y"
    coordinates <-  strsplit(coordinates," ")[[1]]
    coordinates <- setdiff(coordinates,c("+","%in%",":","/","")) ## "" for hidden linebreaks (?)
  } else {
    stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
    ## very old code handling old syntax with (1|pos) and default values of the coordinates argument
    coordinates <- c("x","y") ## back compat
  }
  coordcols <- which(datanames %in% coordinates)
  if ( length(coordcols) != length(coordinates) ) {
    stop("variables 'coordinates' not all found in the 'data'")
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
  return(coordinates) ## should be ordered as bars[[3]] (important for predict)
}



.preprocess_covStruct <- function(covStruct) {
  if ( ! inherits(covStruct,"list")) stop("covStruct must inherit from class 'list'.")
  types <- attr(covStruct,'types') ## 1st way of specifying types
  if (is.null(types)) types <- names(covStruct) ## 2nd way of specifying types
  known_types <- c("adjMatrix","corrMatrix","precision","SAR_WWt","distMatrix", "IMRF") 
  checktypes <- setdiff(types,c(known_types,"")) ## "" for unhandled ranefs
  if (length(checktypes)) stop(paste("Unhandled name(s)/type(s)",
                                     paste0("'",checktypes,"'",collapse=", "),"in 'covStruct'."))
  resu <- vector("list",length(covStruct)) ## list with sublists(?); compatible with spaMM3.0 extended syntax
  for (lit in seq_along(covStruct)) {
    if (types[[lit]]=="precision") {
      resu[[lit]] <- forceSymmetric(covStruct[[lit]])
    } else resu[[lit]] <- covStruct[[lit]]
  }
  names(resu) <- types ## repeated names possible
  attr(resu,'AMatrices') <- AMatrices <- attr(covStruct,'AMatrices')
  if ( ! is.null(AMatrices) && ! is.list(AMatrices)) {stop("attr(covStruct,'AMatrices') must be a list.")} 
  return(resu)
}

.check_corrMatrix <- function(corrMatrix) {
  if (is.list(corrMatrix)) {
    dim_corrMatrix <- dim(corrMatrix[[1]])
  } else dim_corrMatrix <- dim(corrMatrix)
  if (dim_corrMatrix[1L]!=dim_corrMatrix[2L])  stop("corrMatrix is not square") 
}


.check_subset_corrMatrix <- function(corrMatrix,ZA) {
  ZAnames <- colnames(ZA) ## set by .calc_Zlist() or .calc_ZAlist(), with two cases for corrMatrix 
  if (is.null(ZAnames)) {
    stop("NULL colnames in (a block of) the design matrix for random effects. Some mishandling of 'AMatrices'?")
  }
  if (inherits(corrMatrix,"dist")) {
    corrnames <- labels(corrMatrix) ## unclear
  } else if (inherits(corrMatrix,c("matrix","Matrix"))) {
    corrnames <- rownames(corrMatrix)
  } else if ( inherits(corrMatrix,"precision")) {
    corrnames <- rownames(corrMatrix[["matrix"]])
  } else stop("Unhandled class of corrMatrix object.")
  if (is.null(corrnames)) {
    mess <- paste("(!) corrMatrix without labels or row names: the grouping levels, in order",
                  paste0(ZAnames[1L:min(5L,length(ZAnames))], collapse=" "),if(length(ZAnames)>5L){"...,"} else{","},
                  "\n are matched in this order to rows and columns of corrMatrix, without further check.",
                  "\n This may cause later visible errors (notably, wrongly dimensioned matrices)",
                  "\n or even silent errors. See help(\"corrMatrix\") for a safer syntax.")
    warning(mess)
  } else if (is.null(colnames(corrMatrix))) {
    if (inherits(corrMatrix, c("matrix", "Matrix"))) {
      colnames(corrMatrix) <- corrnames
    }
    else if (inherits(corrMatrix, "precision")) {
      colnames(corrMatrix[["matrix"]]) <- corrnames
    }
  }
  if ( length(setdiff(ZAnames,corrnames)) ==0L ) { ## i.e. all ZAnames in corrnames
    ## : should be the case when generator = "as.factor"
    if ( inherits(corrMatrix,"precision")) { ## reordering only 
      if (any(corrnames!=ZAnames)) {
        cov_info_mat <- corrMatrix[ZAnames,ZAnames] 
        if ( morelevels <- length(setdiff(corrnames,ZAnames))) {
          message(paste("Note: precision matrix has", morelevels, "more levels than there are in the data."))
        }
      } else cov_info_mat <- corrMatrix
    } else if ( length(setdiff(corrnames,ZAnames)) || any(corrnames!=ZAnames) ) { # reordering and subsetting
      if (inherits(corrMatrix,"dist")) {
        cov_info_mat <- (proxy::as.matrix(corrMatrix,diag=1)[ZAnames,ZAnames]) ## IF diag missing in input corrMatrix THEN assume a correlation matrix
        ## it's not useful to convert back to dist (either uglily by as.dist(); or package 'seriation' has (permute.dist-> C code)
      } else cov_info_mat <- corrMatrix[ZAnames,ZAnames]  
    } else cov_info_mat <- corrMatrix ## orders already match
  } else {
    ## : expected when generator = ".ULI"
    if ( ! is.null(corrnames)) {
      message("Incompletely checked case: corrMatrix may be invalid, or with complex grouping term\n that spaMM is not able to match to the names of corrMatrix.")
      message("First grouping levels will be matched\n  to first rows of corrMatrix, without further check. \n See help(\"corrMatrix\") for a safer syntax.")
      if ( length(corrnames)!=length(ZAnames)){ 
        message("First grouping levels could not be matched to first rows of corrMatrix, because of inconsistent dimensions.")
        stop("The dimension of corrMatrix does not match the number of levels of the grouping variable.")
        #message("Summary of grouping levels that do not appear in the corrMatrix:")
        #str(checknames)
      }
    }
    cov_info_mat <- corrMatrix ## no clear reordering
  }
  return(cov_info_mat)
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
  ranPars <- .post_process_parlist(ranPars,corr_families=processed$corr_info$corr_families) 
  ranPars <- .canonizeRanPars(ranPars=ranPars, corr_info=processed$corr_info, rC_transf=.spaMM.data$options$rC_transf) ## with init.HLfit as attribute # expands hyper params.
  ########################################################################################################################
  # * assigns geo_envir <- .get_geo_info(...)
  # * modifies processed$AUGI0_ZX$envir by .init_precision_info(...) 
  # * computes processed$AUGI0_ZX$envir$LMatrices except for ranCoefs (the latter being filled in HLfit_body)
  .assign_geoinfo_and_LMatrices_but_ranCoefs(processed, corr_types, spatial_terms, ranPars, control.dist, 
                                argsfordesignL=dotlist[intersect(names(dotlist),names(formals(mat_sqrt)))] )
  if (any(vec_normIMRF <- processed$AUGI0_ZX$vec_normIMRF ) ) {
    processed$ZAlist <- .normalize_IMRF(processed$ZAlist, 
                                        vec_normIMRF=vec_normIMRF,
                                        #ranges=processed$hyper_info$ranges,
                                        strucList=processed$AUGI0_ZX$envir$LMatrices)
    ZAfix <- .ad_hoc_cbind(processed$ZAlist, as_matrix=FALSE)   
    if (processed$sparsePrecisionBOOL) {
      if ( ! inherits(ZAfix,"sparseMatrix")) ZAfix <- as(ZAfix,"dgCMatrix") # .Rcpp_as_dgCMatrix(ZAfix) ## 
      processed$AUGI0_ZX$is_unitary_ZAfix <- FALSE
    }
    processed$AUGI0_ZX$ZAfix <- ZAfix
  }
  ########################################################################################################################
  ###
  if ( (! is.null(processed$return_only)) && processed$augZXy_cond) {
    hlfit <- do.call(.spaMM.data$options$augZXy_fitfn,list(processed=processed, ranFix=ranPars))
    if (FALSE) { ## this test interferes with the results (fitme3, fitme6 tests)
      processed$augZXy_cond <- FALSE
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
      processed$augZXy_cond <- TRUE
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
  } else class(hlfit) <- c("HLCor", class(hlfit))
  #### Infos in the final fit object: 
  hlfit$spatial_terms <- spatial_terms
  info_uniqueGeo <- msd_arglist <- list()
  is_uniqueGeo_needed <- ( (! is.na(corr_types)) & (corr_types=="Matern" | corr_types=="Cauchy" | corr_types=="AR1"
                           | corr_types=="IMRF")) ## IMRF: mapMM expects it.
  for (rd in which(is_uniqueGeo_needed)) {
    char_rd <- as.character(rd)
    geo_envir <- .get_geo_info(processed, which_ranef=rd, which="uniqueGeo", 
                               dist.method=control.dist[[char_rd]]$dist.method) 
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
      locmc[[1L]] <- get("HLCor.obj", asNamespace("spaMM")) # as.name("HLCor.obj") ## replaces "f" !
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
    urP <- unlist(.canonizeRanPars(ranefParsList, corr_info=processed$corr_info,checkComplete=FALSE, rC_transf=.spaMM.data$options$rC_transf))
    processed$port_env$prefix <- paste0("HLCor for ", paste(signif(urP,6), collapse=" "), ": ")
  } 
  ranPars <- .modify_list(HLCor.call$ranPars, ranefParsList)
  rpType <- .modify_list(attr(HLCor.call$ranPars,"type"),attr(skeleton,"type"))
  moreargs <- attr(skeleton,"moreargs")
  ranPars <- .expand_hyper(ranPars, processed$hyper_info,moreargs=moreargs) ## input ranPars contains both unconstrained ranPars and $hyper
  # => failing to expand leads to unconstrained optimization
  if ( ! is.null(ranPars$resid) ) { ## mmm FIXME never operational ?
    resid_ranPars <- structure(ranPars$resid, ## but not sure that the attributes are necessary...
                               type=rpType$resid, 
                               moreargs=moreargs$resid)
    # canonize bc fitme_body(,fixed=.) does not handle transformed parameters
    processed$residProcessed$envir$ranPars <- .canonizeRanPars(ranPars=resid_ranPars,
                                corr_info=processed$residProcessed$corr_info,
                                checkComplete = FALSE, rC_transf=.spaMM.data$options$rC_transf) 
    ranPars$resid <- NULL
    rpType$resid <- NULL
    moreargs$resid <- NULL
  }
  HLCor.call$ranPars <- structure(ranPars, ## adds given values of the optimized variables 
                                  type=rpType, ## adds "fix"'s... somewhat confusing 
                                  moreargs=moreargs )
  # ranPars may have $trLambda (from notlambda) for what is optimized,
  #              and $lambda (from ranPars$lambda) for what was fixed in the whole outer fit  
  HLCor.call[[1L]] <- get("HLCor", asNamespace("spaMM")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
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


