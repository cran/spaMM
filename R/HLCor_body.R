## the following must match the'unique' method is ULI as explained there
.calcUniqueGeo <- function(data) {
  redondGeo <- apply(data,1,paste,collapse=" ") ## creates character string
  dfforunique <- cbind(data,redondGeo) ## associates rownames of data to redondGeo
  uniqueGeo <- unique(dfforunique[,ncol(dfforunique),drop=FALSE]) ## keeps rownames of first instances
  uniqueGeo <- data[rownames(uniqueGeo),,drop=FALSE] ## uses rownames, 'unique' numeric values based on character representations 
  return(uniqueGeo)
}


.extract_check_coords <- function(spatial.model,datanames) {
  if ( ! is.null(spatial.model)) {
    bars <- spatial.model[[2]] 
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

.extract_check_coords_within <- function(spatial.model) {
  bars <- spatial.model[[2]] 
  coordinates <- .DEPARSE(bars[[3]]) ## "x + y"
  coordinates <-  strsplit(coordinates," ")[[1]]
  if (length(grep("/",coordinates))>0) {
    stop(paste("'/' not yet handled in",spatial.model))
  } else if (length(grep_in <- grep("%in%|:",coordinates))>1L) {
    stop(paste("multiple nesting not yet handled in",spatial.model))
  } else if (length(grep_in <- grep("%in%|:",coordinates))>0L) {
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

.get_corr_prec_from_covStruct <- function(covStruct) {
  corrMatrix <- covStruct[["corrMatrix"]]
  if (is.null(corrMatrix)) {
    matrix_ <- covStruct[["precision"]]
    if (is.null(matrix_)) stop("missing covariance structure for corrMatrix model")
    corrMatrix <- structure(list(matrix=matrix_),class=c("list","precision")) ## so that further code can detect its a precision matrix
  }
  return(corrMatrix)
}

.check_corrMatrix <- function(corrMatrix) {
  if (is.list(corrMatrix)) {
    dim_corrMatrix <- dim(corrMatrix[[1]])
  } else dim_corrMatrix <- dim(corrMatrix)
  if (dim_corrMatrix[1L]!=dim_corrMatrix[2L])  stop("corrMatrix is not square") 
}


.check_subset_corrMatrix <- function(corrMatrix,ZAlist) {
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
  whichranef <- which(attr(attr(ZAlist,"ranefs"),"type")=="corrMatrix")
  ZAnames <- colnames(ZAlist[[whichranef]]) ## set by .spMMFactorList(), with two cases for corrMatrix 
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


HLCor_body <- function(processed,
                  ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  distMatrix,uniqueGeo=NULL,adjMatrix,corrMatrix=NULL, covStruct=NULL,
                  control.dist=list(),
                  ...) { 
  dotlist <- list(...)
  ################# 
  data <- processed$data
  verbose <- processed$verbose
  predictor <- processed$predictor 
  spatial.terms <- .findSpatial(predictor)
  spatial.model <- spatial.terms[[1L]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1L]]) 
  } else {
    if ( ! is.null(corrMatrix)) {
      mess <- pastefrom("corrMatrix argument despite no corrMatrix term in formula:",prefix="(!) From ")
      message(mess)
      stop("This syntax is obsolete; add a corrMatrix(...) term in the formula.")
    } ## ELSE more generic message: 
    stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
  }
  ## convert back ranPars to canonical scale:
  rpblob <- .canonizeRanPars(ranPars=ranPars,corr.model=corr.model) ## also provides some init.HLfit 
  ranPars <- rpblob$ranPars
  trueCorrpars <- rpblob$trueCorrpars
  rho <- ranPars$rho
  #
  coordinates <- NULL
  test.in <- FALSE
  if (!is.null(covStruct)) covStruct <- .preprocess_covStruct(covStruct)
  ### ensure LMatrix in predictor: 
  ## if it is currently absent, first provide corr matrix or its symSVD, from which Lunique will be computed using designL.from.Corr
  if (is.null(Lunique <- attr(predictor,"LMatrix"))) { ## should exist only if the user provided it explicitly. Should not be used in programming.
    symSVD <- NULL
    adj_rho_is_inner_estimated <- (corr.model=="adjacency"
#                               && ! is.null(attr(ranPars,"type")) ## through corrHLfit (or fitme !?) or some direct HLCor call
                               && identical(attr(ranPars,"type")$rho,"var") ## can occur in direct call of HLCor 
    )    
    if (corr.model %in% c("adjacency")) {
      if ( missing(adjMatrix) ) {
        adjMatrix <- covStruct[["adjMatrix"]]
        if (is.null(adjMatrix)) stop("missing 'adjMatrix' for adjacency model")
      }
      ## no nugget in the adjacency model... ## (use additional ranef instead)      
      symSVD <- attr(adjMatrix,"symSVD")
      if (is.null(symSVD) && adj_rho_is_inner_estimated) { ## can occur in direct call of HLCor 
        if (isSymmetric(adjMatrix)) {
          symSVD <- sym_eigen(adjMatrix)
          attr(adjMatrix,"symSVD") <- symSVD
        }             
      }
      if (is.null(symSVD)) {
        cov_info_mat <- as.matrix(solve(diag(nrow(adjMatrix))-rho*(adjMatrix))) 
      } else {
        symSVD$adjd <- symSVD$d
        # FR->FR remove $d from symSVD?
      }
    }  else if (corr.model %in% c("SAR_WWt")) { 
      adjMatrix <- covStruct[["SAR_WWt"]]
      if (is.null(adjMatrix)) stop("missing covariance structure for SAR_WWt model")
      UDU. <- attr(adjMatrix,"UDU.")
      if (is.null(UDU.)) {
        cov_info_mat <- as.matrix(solve(diag(nrow(adjMatrix))-rho*(adjMatrix))) 
      } else {
        cov_info_mat <- UDU.$u %*% sweep(UDU.$u.,MARGIN=1,1/(1-rho*UDU.$d),`*`) 
      }
      cov_info_mat <- .tcrossprodCpp(cov_info_mat,NULL)
    }  else if (corr.model=="AR1" && ! identical(processed$sparsePrecisionBOOL,TRUE)) {
      if (is.null(distMatrix)) { ## if not precomputed in corrHLfit_body/fitme_body
        coordinates <- .get_coordinates(spatial.model=spatial.model, data=data)
        coord_within <- .extract_check_coords_within(spatial.model=spatial.model) 
        coords_nesting <- setdiff(coordinates,coord_within)
        if (length(coords_nesting)) {
          if (.getProcessed(processed,"sparsePrecisionBOOL",1L)) {
            message("sparse precision method not available for nested AR1 effects.")
            .setProcessed(processed,"sparsePrecisionBOOL","FALSE")
          }
        }
        locarglist <- list(data=data,distMatrix=dotlist$distMatrix, uniqueGeo=uniqueGeo, 
                           coords_nesting=coords_nesting, coordinates=coordinates, dist.method = control.dist$dist.method)
        geoMats <- do.call(".makeCheckGeoMatrices",locarglist) ## computes matrices which are NULL on input
        distMatrix <- geoMats$distMatrix   
        uniqueGeo <- geoMats$uniqueGeo   
      }
      cov_info_mat <- trueCorrpars$ARphi^distMatrix  
      cov_info_mat[distMatrix==Inf] <- 0 ## should not be necess, but is.
    } else  if (corr.model %in% c("Matern")) {
      txt <- paste(spatial.model[[2]][[3]]) ## the RHS of the ( . | . ) 
      if (length(grep("%in%",txt))>0) {
        stop("(!) Matern( . | <coord> %in% <grp>) is not yet handled.")
        test.in <- TRUE ## should be useful when this case will be handled
      } 
      ## in a typical call from corrHLfit the following test should be FALSE because uniqueGeo and maybe distMatrix should have been precomputed
      if ((length(rho)>1 || missing(distMatrix)) && is.null(uniqueGeo)) { ## all cases where we need uniqueGeo
        coordinates <- .extract_check_coords(spatial.model=spatial.model,datanames=names(data))
        uniqueGeo <- .calcUniqueGeo(data=data[,coordinates,drop=FALSE]) ## keeps the names of first instances of the coordinates in data
      } 
      ## then compute scaled distances from unscaled info, for HLfit call
      msd.arglist <- list(rho = rho)
      msd.arglist$`dist.method` <- control.dist$`dist.method` ## may be NULL
      if (length(rho)>1L) {
        msd.arglist <- c(msd.arglist,list(uniqueGeo=uniqueGeo))
        msd.arglist$`rho.mapping` <- control.dist$`rho.mapping` ## may be NULL
      } else {
        if ( missing(distMatrix)) { 
          dist.arglist <- list(x=uniqueGeo)
          dist.arglist$method <- control.dist$dist.method ## may be NULL
          distMatrix <- do.call(proxy::dist,dist.arglist)
        } 
        msd.arglist <- c(msd.arglist,list(distMatrix=distMatrix))
      }
      cov_info_mat <- do.call("make_scaled_dist",msd.arglist)
      ## at this point if a single location, dist_mat should be dist(0) and make_scaled_dist was modified to that effect
      if ( nrow(cov_info_mat)>1 ) { ## >1 locations
        norho <- trueCorrpars 
        norho$rho <- NULL ## because the MaternCorr input will be an already scaled distance 'cov_info_mat'
        cov_info_mat <- do.call(MaternCorr,args=c(norho,list(d=cov_info_mat)))        
      } 
    }
    if (is.null(processed$sparsePrecisionBOOL)) { ## if next operations not already performed by fitme / corrHLfit
      if (corr.model== "corrMatrix") { ## could cover all specifications of a constant 'cov' Matrix without correlation parameter
        if (is.null(corrMatrix)) corrMatrix <- .get_corr_prec_from_covStruct(covStruct)
        .check_corrMatrix(corrMatrix)
        sparse_precision <- inherits(corrMatrix,"precision")
        if ( ! sparse_precision) sparse_precision <- .determine_sparse_precision(processed, corr.model, dotlist$init.HLfit) ## checks spaMM option
        ignored_as_no_need_for_local_copy_processed <- .reset_processed_sparse_precision(processed, inherits(corrMatrix,"precision"))
      } else sparse_precision <- .determine_sparse_precision(processed, corr.model, dotlist$init.HLfit) ## checks spaMM option
    }
    if (corr.model== "corrMatrix") { ## could cover all specifications of a constant 'cov' Matrix without correlation parameter
      cov_info_mat <- .check_subset_corrMatrix(corrMatrix,ZAlist=processed$ZAlist) ## correlation or precision...
    } 
    if (verbose["trace"] && length(trueCorrpars)>0) print(unlist(trueCorrpars))
    ## call designL.from.Corr if Lunique not available
    ## Lunique can be available here either bc the user provided it explictly (will not occur in public usage)
    ##   or if we add an explicit calculation above
    sparse_Qmat <- NULL
    if (is.null(Lunique)) { ##
      ## The following code implies that sparse_precision is not alway enforced and $method not always set
      if (corr.model=="corrMatrix" && identical(processed$sparsePrecisionBOOL,TRUE)) {
        processed$AUGI0_ZX$envir$method <- "def_AUGI0_ZX_sparsePrecision"
        sparse_Qmat <- Matrix::drop0(cov_info_mat[["matrix"]])
      } else if (corr.model=="adjacency" 
          && ! adj_rho_is_inner_estimated 
          && identical(processed$sparsePrecisionBOOL,TRUE)) {
        ## Implies call from fitme_body with outer rho estim.
        processed$AUGI0_ZX$envir$method <- "def_AUGI0_ZX_sparsePrecision"
        sparse_Qmat <- - rho * Matrix::drop0(adjMatrix)
        diag(sparse_Qmat) <- diag(sparse_Qmat)+1
        ################################L_Q <- .designL.from.Qmat(Qmat) ## solve(tcrossprod(LMatrix)) = Qmat or
        ## Cholesky gives proper LL' (think LDL')  while chol() gives L'L...
      } else if (corr.model=="AR1" && identical(processed$sparsePrecisionBOOL,TRUE)) { # F I X M E nesting
        processed$AUGI0_ZX$envir$method <- "def_AUGI0_ZX_sparsePrecision"
        ARphi <- ranPars$ARphi
        # equivalent to nlme's AR1_fact() in corStruct.c
        t_chol_Q <- Diagonal(x=c(rep(1,control.dist$AR1_tmax-1L),sqrt(1-ARphi^2))) 
        diag(t_chol_Q[,-1]) <- -ARphi 
        t_chol_Q <- t_chol_Q/sqrt(1-ARphi^2)
        ## bc we nned the CHMfactor per se, we do not directly provide chol_Q. This is not elegant.
        sparse_Qmat <- crossprod(t_chol_Q) 
        ## as(Matrix::Cholesky(sparse_Qmat,perm=FALSE,LDL=FALSE),"sparseMatrix") gives back Linv...
        ## older code used: 
        #Lunique <- as.matrix(solve(tLinv/sqrt(1-ARphi^2))) ## corrmat is tcrossprod(Lunique): we keep tcrossprod but L' is tri.sup ! 
        #attr(Lunique,"type") <- "cholR_RRt" ## not equivalent to base::chol() which whould produce cholR_tRR 
      } else {
        # this is HLCor, hence there must be a dense correlation structure
        ## densePrecision cases with outer rho estimation, and some residual inner rho estimation cases
        if ( ! is.null(symSVD)) {
          symSVD$d <- 1/(1-rho*symSVD$adjd) ## from adjMatrix to correlation matrix
          # outer optim -> LMatrix recomputed from this for each rho  
          Lunique <- try(designL.from.Corr(symSVD=symSVD)) ## usin $d not $adjd
          # Lunique is num, with attr(*, "type")= chr "symsvd" ....
          if ( adj_rho_is_inner_estimated ) { # contrived way of construction Lunique with the correct attributes.
            Lunique[] <- attr(Lunique,"symsvd")$u   ## "[] <- " keeps attributes... except for Matrix...
          }
        } else { ## cov_info_mat must exist
          argsfordesignL <- dotlist[intersect(names(dotlist),names(formals(designL.from.Corr)))] 
          if (processed$HL[1L]=="SEM") argsfordesignL$try.chol <- FALSE
          if (inherits(cov_info_mat,"dist")) {
            cov_info_mat <- as.matrix(cov_info_mat)
            diag(cov_info_mat) <- 1L ## IF diag missing in input corrMatrix THEN assume a correlation matrix
          } ## else full matrix may be a COV matrix with non-unit diag
          Lunique <- try(do.call(designL.from.Corr,c(list(m=cov_info_mat),argsfordesignL)))
        }
        if (inherits(Lunique,"try-error")) { 
          print("correlation parameters were:",quote=FALSE) ## makes sense if designL.from.Corr already issued some warning
          print(unlist(trueCorrpars))    
          stop()
        }
      }
      if (identical(processed$sparsePrecisionBOOL,TRUE)) {
        if (is.null(processed$AUGI0_ZX$envir$precisionFactorList)) .init_precision_info(processed) ## modifies processed$AUGI0_ZX$envir  
        if ( ! is.null(sparse_Qmat)) {
          #processed$AUGI0_ZX$envir$Qmat <- sparse_Qmat ## fixme do we need this copy ?
          Q_CHMfactor <- Matrix::Cholesky(sparse_Qmat,LDL=FALSE,perm=FALSE) ## called for each corrPars
        }
        if (identical(processed$AUGI0_ZX$envir$method, "def_AUGI0_ZX_sparsePrecision")) {
          ## fixme this block of code only to identify the good ranef...
          types <- processed$AUGI0_ZX$envir$types
          for (it in seq_len(length(types))) {
            if (types[it] %in% c("adjacency","AR1") || (types[it]=="corrMatrix" && inherits(cov_info_mat,"precision"))) { 
              processed$AUGI0_ZX$envir$precisionFactorList[[it]] <- list(Qmat=sparse_Qmat,
                                                                         chol_Q=as(Q_CHMfactor, "sparseMatrix"))
            } # else precisionFactorList[[it]] <- Diagonal(n=ncol(processed$ZAlist[[it]])) ## already pre-filled
          }
          ## fitme sparse_precision has an incomplete symSVD=> corr mattrix not computed, 
          ##    and try(designL.from.Corr(symSVD=symSVD)) fails. Instead use code always valid:
          Lunique <- solve(Q_CHMfactor,system="Lt") #  L_Q^{-\top}=LMatrix_correlation 
          attr(Lunique, "type") <- "from_Q_CHMfactor"
          attr(Lunique,"Q_CHMfactor") <- Q_CHMfactor ## FIXME store directly the Q_CHMfactor rather than Lunique ?
          # Lunique is a dtCMatrix; it is for single correlated effect, and still used in HLfit_body; 
          ## whether it is used or not in MME_method (sXaug_...), a lot of other code still expects it
        }
      }
    }
    attr(predictor,"%in%") <- test.in
    attr(Lunique,"corr.model") <- corr.model
    attr(Lunique,"ranefs") <- unlist(lapply(spatial.terms,.DEPARSE)) ## essentiel pour la construction de ZAL!
  }
  processed$AUGI0_ZX$envir$LMatrix <- Lunique ## either the constant from the predictor or something just defined
  processed$AUGI0_ZX$envir$adj_symSVD <- symSVD ## may be NULL
  processed$predictor <- predictor
  ###
  HLFormals <- names(formals(HLfit))
  good_dotnames <- intersect(names(dotlist),HLFormals)
  if (length(good_dotnames)>0L) {
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
          utils::str(rpType)
          utils::str(rpTypes)
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
  #
  hlfit$control.dist <- control.dist
  attr(hlfit,"info.uniqueGeo") <- uniqueGeo ## whether Matern or not (eg AR1)
  if (corr.model %in% c("Matern")) {
    ## we try to remove the big matrix if it can be reconstructed
    if ( ! is.null(dM <- msd.arglist$distMatrix) ) { 
      if ( ! is.null(distcall <- attr(dM,"call"))) {
        msd.arglist$distcall <- distcall ## save the call, eg language proxy::dist(x = uniqueGeo, method = dist.method)
        msd.arglist$distMatrix <- NULL ## removes the big matrix
      }
    }
    attr(hlfit,"dist_info") <- msd.arglist
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
    fitlist <- lapply(seq_len(length(processed)), function(it){
      locmc <- mc
      locmc[[1L]] <- as.name("HLCor.obj") ## replaces "f" !
      locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
      locmc$processed <- processed[[it]] ## The data are in processed !
      locmc$distMatrix <- mc$distMatrix[[it]] ## but the matrices are not HLfit args hence not in processed ! 
      locmc$uniqueGeo <- mc$uniqueGeo[[it]]
      eval(locmc)
    }) ## a pure list of HLCor objects
    resu <- sum(unlist(fitlist))
    return(resu)
  } else { ## there is one processed for a single data set 
    family <- processed$family
    data <- processed$data
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


