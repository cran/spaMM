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
  if (length(grep("/",coordinates))>0) {
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
  checktypes <- setdiff(types,known_types)
  if (length(checktypes)) stop(paste("Unhandled name(s)/type(s)",paste("'",checktypes,"'",sep="",collapse=", "),"in 'covStruct'."))
  resu <- list() ## list with sublists; next few lines handle multiple values of each type.
  for (st in known_types) if (any(match_ <- types==st)) {
    if (st=="precision") {
      resu[[st]] <- forceSymmetric(covStruct[[which(match_)]])
    } else resu[[st]] <- covStruct[[which(match_)]]
  }
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
    if (  (is_superset <- (length(setdiff(corrnames,ZAnames))>0)) || 
          (is_not_in_order_Zcols <- any(corrnames!=ZAnames)) ) { ## ...but superset, or not same order
      if ( inherits(corrMatrix,"precision")) {
        if ( is_superset ) {
          # could be corrected though subsampling of the _correlation_ factor matrix, but:
          stop("Precision matrix involves levels of the grouping variable absent from the data. 
               spaMM asks the user to correct this rather than itself to correct it automatically.")
        } else if (is_not_in_order_Zcols) cov_info_mat <- corrMatrix[ZAnames,ZAnames] ## reordering
        ## do nothing, because it is incorrect to subset a precision matrix. return is a list...
      } else if (inherits(corrMatrix,"dist")) {
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
                  control.dist=list(),
                  ...) { ## dots for HLfit
  dotlist <- list(...)
  spatial_terms <- attr(processed$ZAlist,"exp_spatial_terms")
  corr_types <- processed$corr_info$corr_types
  ## convert back ranPars to canonical scale:
  rpblob <- .canonizeRanPars(ranPars=ranPars,corr_types=corr_types) ## also provides some init.HLfit 
  ranPars <- rpblob$ranPars
  ########################################################################################################################
  # * assigns geo_envir <- .get_geo_info(...)
  # * modifies processed$AUGI0_ZX$envir by .init_precision_info(...) 
  # * computes processed$AUGI0_ZX$envir$LMatrices except for ranCoefs (the latter being filled in HLfit_body)
  .assign_geoinfo_and_LMatrices_but_ranCoefs(processed, corr_types, spatial_terms, ranPars, control.dist, rpblob$trueCorrpars, 
                                argsfordesignL=dotlist[intersect(names(dotlist),names(formals(designL.from.Corr)))] )
  ########################################################################################################################
  ###
  HLFormals <- names(formals(HLfit)) 
  good_dotnames <- intersect(names(dotlist),HLFormals)
  if (length(good_dotnames)) {
    HL.info <- dotlist[good_dotnames]
  } else HL.info <- list()
  ## all printing in HLfit is suppressed by default
  HL.info$processed <- processed
  ## convert ranPars to ranFix + init.HLfit
  ## allows log and not log:
  varNames <- names(which(attr(ranPars,"type")=="var"))
  HL.info$init.HLfit[varNames] <- ranPars[varNames] ## inherits values from corrHLfit(...,init.HLfit(...))... or thorugh fitme !
  fixNames <- setdiff(names(ranPars),varNames) ## "fix" or "outer"
  if (!is.null(fixNames)) { ## could be NULL for corrMatrix case
    ranFix <- ranPars[fixNames] ## 11/2014 as there is no other source for ranFix
    typelist <- list() 
    typelist[fixNames] <- "fix" ## corrected by next lines
    if ( ! is.null(rPtype <- attr(ranPars,"type"))) { ## it may not exist in case of direct call of HLCor, 
      ## else it has elements "fix" or "outer"
      typelist[names(rPtype)] <- rPtype ## "outer" may override "fix"
    }
    attr(ranFix,"type") <- typelist 
    if (spaMM.getOption("wDEVEL2")) {  
      ## rPparlist NULLin direct call of HLCor:
      if ( is.null(rPparlist <- attr(ranPars,"parlist"))) {
        parlist <- .merge_parlist(NULL,new=ranFix,types="fix") ## will use attr(ranFix,"types") --> "fix"
        parlist <- .merge_parlist(parlist,new=HL.info$init.HLfit,types="var")
        attr(ranFix,"parlist") <- parlist ## pallist has fix+var contrary to ranFix[]
        #str(ranFix)
      } else {
        ## ranPars may have a preesisting parlist but eg with trLambda instead of lambda
        if ( ! identical(rpType <- attr(ranPars,"type"), lapply(rpTypes <- attr(rPparlist,"types"),tolower))) {
          #stop(" ! identical(attr(ranPars,\"type\"),attr(attr(ranPars,\"parlist\"),\"types\"))")
          # can be different bc rpType handles named vectors (rho,lambda) differently...
          str(rpType)
          str(rpTypes)
        }
      }
    }
    HL.info$ranFix <- ranFix
  }
  hlfit <- do.call("HLfit",HL.info) 
  ## Here there was debug code that saved HL.info in case of error; before 1.8.5
  if ( ! is.null(processed$return_only)) {
    return(hlfit)    ########################   R E T U R N   a list with $APHLs, with class "list"
  } else class(hlfit) <- c(class(hlfit),"HLCor")
  #### Infos in the final fit object: 
  hlfit$control.dist <- control.dist
  hlfit$spatial_terms <- spatial_terms
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if ( (! is.na(corr_type)) && (corr_type=="Matern" || corr_type=="AR1")) {
      geo_envir <- .get_geo_info(processed, which_ranef=it, which="", 
                                 dist.method=control.dist$dist.method) 
      attr(hlfit,"info.uniqueGeo") <- geo_envir$uniqueGeo ## spaMM 3.0 will require a list of matrices
      if (corr_type=="Matern") {
        msd.arglist <- attr(processed$AUGI0_ZX$envir$LMatrices[[it]],"msd.arglist") 
        ## we try to remove the big matrix if it can be reconstructed
        if ( ! is.null(dM <- msd.arglist$distMatrix) ) { 
          if ( ! is.null(distcall <- attr(dM,"call"))) {
            msd.arglist$distcall <- distcall ## save the call, eg language proxy::dist(x = uniqueGeo, method = dist.method)
            msd.arglist$distMatrix <- NULL ## removes the big matrix
          }
        }
        attr(hlfit,"dist_info") <- msd.arglist ## spaMM 3.0 soudl store a list of msd.arglists
      }
    }
  }
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
  forGiven <- relist(ranefParsVec,skeleton) ## given values of the optimized variables 
  ## ... relist keeps the RHOMAX... attributes from the skeleton, but the partial copy into ranPars does not.
  if (spaMM.getOption("wDEVEL2")) {
    parlist <- attr(HLCor.call$ranPars,"parlist")
    parlist <- .merge_parlist(parlist,new=forGiven,types="fix")## consistently with previous code
    attr(parlist,"RHOMAX") <- attr(skeleton,"RHOMAX")
    attr(parlist,"NUMAX") <- attr(skeleton,"NUMAX")
    attr(HLCor.call$ranPars,"parlist") <- parlist
  }
  notlambda <- setdiff(names(forGiven),"lambda")
  HLCor.call$ranPars$lambda[names(forGiven$lambda)] <- forGiven$lambda
  HLCor.call$ranPars[notlambda] <- forGiven[notlambda] ## do not wipe out other fixed, non optimized variables
  # ranPars may have $trLambda (from notlambda) for what is optimized,
  #              and $lambda (from ranPars$lambda) for what was fixed in the whole outer fit  
  attr(HLCor.call$ranPars,"RHOMAX") <- attr(skeleton,"RHOMAX")
  attr(HLCor.call$ranPars,"NUMAX") <- attr(skeleton,"NUMAX")
  types <- attr(skeleton,"type")
  attr(HLCor.call$ranPars,"type")[names(types)] <- types
  HLCor.call[[1L]] <- quote(spaMM::HLCor)
  hlfit <- eval(HLCor.call) ## retruns fit or call 
  #
  if (is.call(hlfit)) {return(hlfit)} ## HLCorcall
  #
  aphls <- hlfit$APHLs
  resu <- aphls[[objective]]
  return(resu) #
}


