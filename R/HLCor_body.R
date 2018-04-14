## the following must match the'unique' method is ULI as explained there
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
  known_types <- c("adjMatrix","corrMatrix","precision","SAR_WWt","distMatrix") 
  checktypes <- setdiff(types,c(known_types,"")) ## "" for unhandled ranefs
  if (length(checktypes)) stop(paste("Unhandled name(s)/type(s)",paste("'",checktypes,"'",sep="",collapse=", "),"in 'covStruct'."))
  resu <- vector("list",length(covStruct)) ## list with sublists(?); compatible with spaMM3.0 extended syntax
  for (lit in seq_along(covStruct)) {
    if (types[[lit]]=="precision") {
      resu[[lit]] <- forceSymmetric(covStruct[[lit]])
    } else resu[[lit]] <- covStruct[[lit]]
  }
  names(resu) <- types ## repeated names possible
  return(resu)
}

.check_corrMatrix <- function(corrMatrix) {
  if (is.list(corrMatrix)) {
    dim_corrMatrix <- dim(corrMatrix[[1]])
  } else dim_corrMatrix <- dim(corrMatrix)
  if (dim_corrMatrix[1L]!=dim_corrMatrix[2L])  stop("corrMatrix is not square") 
}


.check_subset_corrMatrix <- function(corrMatrix,ZA) {
  ## ELSE check descriptors of square matrix:
  if (inherits(corrMatrix,"dist")) {
    corrnames <- labels(corrMatrix)
  } else if (inherits(corrMatrix,c("matrix","Matrix"))) {
    corrnames <- rownames(corrMatrix)
  } else if ( inherits(corrMatrix,"precision")) {
    corrnames <- rownames(corrMatrix[["matrix"]])
  } else stop("Unhandled class of corrMatrix object.")
  if (is.null(corrnames)) {
    message("corrMatrix without labels: first grouping levels are matched\n  to first rows of corrMatrix, without further check.\n This may cause later errors (notably, wrongly dimensioned matrices) \n See help(\"corrMatrix\") for a safer syntax.")
  }
  ZAnames <- colnames(ZA) ## set by .spMMFactorList(), with two cases for corrMatrix 
  if (is.null(ZAnames)) {
    stop("NULL colnames in (a block of) the design matrix for random effects. Some mishandling of 'AMatrix'?")
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
        cov_info_mat <- (as.matrix(corrMatrix)[ZAnames,ZAnames]) 
        ## it's not useful to convert back to dist (either uglily by as.dist(); or package 'seriation' has (permute.dist-> C code)
        diag(cov_info_mat) <- 1L ## IF diag missing in input corrMatrix THEN assume a correlation matrix
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
  t_chol_Q <- Diagonal(x=c(rep(1,n_u-1L),sqrt(1-ARphi^2))) 
  if (n_u>1L) diag(t_chol_Q[,-1,drop=FALSE]) <- -ARphi 
  t_chol_Q <- t_chol_Q/sqrt(1-ARphi^2) # equivalent to nlme's AR1_fact() in corStruct.c
  return(t_chol_Q)
}


HLCor_body <- function(processed, ## single environment
                  ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  control.dist=list(), # info NOT from processed bc modified by corrHLfit_body or fitme_body
                  ...) { ## dots for HLfit
  dotlist <- list(...)
  spatial_terms <- attr(processed$ZAlist,"exp_spatial_terms")
  corr_types <- processed$corr_info$corr_types
  ## convert back ranPars to canonical scale:
  ranPars <- .post_process_parlist(ranPars,corr_types=corr_types) 
  ranPars <- .canonizeRanPars(ranPars=ranPars,corr_types=corr_types) ## with init.HLfit as attribute
  ########################################################################################################################
  # * assigns geo_envir <- .get_geo_info(...)
  # * modifies processed$AUGI0_ZX$envir by .init_precision_info(...) 
  # * computes processed$AUGI0_ZX$envir$LMatrices except for ranCoefs (the latter being filled in HLfit_body)
  .assign_geoinfo_and_LMatrices_but_ranCoefs(processed, corr_types, spatial_terms, ranPars, control.dist, 
                                argsfordesignL=dotlist[intersect(names(dotlist),names(formals(designL.from.Corr)))] )
  ########################################################################################################################
  ###
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
  ## Here there was debug code that saved HL.info in case of error; before 1.8.5
  if ( ! is.null(processed$return_only)) {
    return(hlfit)    ########################   R E T U R N   a list with $APHLs, with class "list"
  } else class(hlfit) <- c(class(hlfit),"HLCor")
  #### Infos in the final fit object: 
  hlfit$spatial_terms <- spatial_terms
  info_uniqueGeo <- msd_arglist <- list()
  is_uniqueGeo_needed <- ( (! is.na(corr_types)) & (corr_types=="Matern" | corr_types=="Cauchy" | corr_types=="AR1"))
  for (rd in which(is_uniqueGeo_needed)) {
    char_rd <- as.character(rd)
    geo_envir <- .get_geo_info(processed, which_ranef=rd, which="", 
                               dist.method=control.dist[[char_rd]]$dist.method) 
    info_uniqueGeo[[char_rd]] <- geo_envir$uniqueGeo 
  }
  attr(hlfit,"info.uniqueGeo") <- info_uniqueGeo
  #attr(hlfit,"control_dists") <- control.dist
  #
  hlfit$call <- "$call removed by HLCor. Use getCall() to extract the call from the object." ## instead of the $call with evaluated arguments
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
      locmc[[1L]] <- as.name("HLCor.obj") ## replaces "f" !
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
  designL.formals <- names(formals(designL.from.Corr))
  makescaled.formals <- names(formals(make_scaled_dist))
  HLnames <- (c(HLCor.formals,names_formals_HLfit,designL.formals,makescaled.formals))  ## cf parallel code in corrHLfit
  HLCor.call <- mc[c(1,which(names(mc) %in% HLnames))] ## keep the call structure
  HLCor.call$ranPars <- structure(.modify_list(HLCor.call$ranPars, relist(ranefParsVec,skeleton)), ## adds given values of the optimized variables 
                                  type=.modify_list(attr(HLCor.call$ranPars,"type"),attr(skeleton,"type")), ## adds "fix"'s... somewhat confusing 
                                  moreargs=attr(skeleton,"moreargs") )
  # ranPars may have $trLambda (from notlambda) for what is optimized,
  #              and $lambda (from ranPars$lambda) for what was fixed in the whole outer fit  
  HLCor.call[[1L]] <- quote(spaMM::HLCor)
  hlfit <- eval(HLCor.call) ## retruns fit or call 
  #
  if (is.call(hlfit)) {return(hlfit)} ## HLCorcall
  #
  aphls <- hlfit$APHLs
  resu <- aphls[[objective]]
  if (objective=="cAIC") resu <- - resu ## for minimization of cAIC (private & experimental)
  return(resu) #
}


