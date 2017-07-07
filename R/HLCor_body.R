## the following must match the'unique' method is ULI as explained there
calcUniqueGeo <- function(data) {
  redondGeo <- apply(data,1,paste,collapse=" ") ## creates character string
  dfforunique <- cbind(data,redondGeo) ## associates rownames of data to redondGeo
  uniqueGeo <- unique(dfforunique[,ncol(dfforunique),drop=FALSE]) ## keeps rownames of first instances
  uniqueGeo <- data[rownames(uniqueGeo),,drop=FALSE] ## uses rownames, 'unique' numeric values based on character representations 
  return(uniqueGeo)
}


`extract.check.coords` <- function(spatial.model,datanames) {
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

HLCor_body <- function(processed,
                  ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  distMatrix,uniqueGeo=NULL,adjMatrix,corrMatrix,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),control.dist=list(),
                  ...) { 
  dotlist <- list(...)
  ################# 
  data <- processed$data
  predictor <- processed$predictor 
  spatial.terms <- findSpatial(predictor)
  spatial.model <- spatial.terms[[1L]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1L]]) 
  } else {
    if ( ! missing(corrMatrix)) {
      mess <- pastefrom("corrMatrix argument despite no corrMatrix term in formula:",prefix="(!) From ")
      message(mess)
      stop("This syntax is obsolete; add a corrMatrix(...) term in the formula.")
    } ## ELSE more generic message: 
    stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
  }
  ## convert back ranPars to canonical scale:
  rpblob <- canonizeRanPars(ranPars=ranPars,corr.model=corr.model) ## also provides some init.HLfit 
  ranPars <- rpblob$ranPars
  trueCorrpars <- rpblob$trueCorrpars
  rho <- ranPars$rho
  #
  coordinates <- NULL
  test.in <- FALSE
  ### ensure LMatrix in predictor: 
  ## if it is currently absent, first provide corr matrix or its symSVD, from which Lunique will be computed using designL.from.Corr
  if (is.null(Lunique <- attr(predictor,"LMatrix"))) { ## Lunique already exists if it is constant 
    symSVD <- NULL
    if (corr.model %in% c("ar1")) {
      ## inv(corrmat) has a simple repres as Linv.tLinv <=> corrmat as R.tR
      ARphi <- ranPars$ARphi
      # equivalent to nlme's AR1_fact() in corStruct.c
      tLinv <- Diagonal(x=c(rep(1,control.dist$ar1_tmax-1L),sqrt(1-ARphi^2))) 
      diag(tLinv[,-1]) <- -ARphi 
      # efficient solve on sparse matrix:
      Lunique <- as.matrix(solve(tLinv/sqrt(1-ARphi^2))) ## corrmat is tcrossprod(Lunique): we keep tcrossprod but L' is tri.sup ! 
      attr(Lunique,"type") <- "cholR_RRt" ## not equivalent to base::chol() which whould produce cholR_tRR 
    }
    adj_rho_is_inner_estimated <- (corr.model=="adjacency"
#                               && ! is.null(attr(ranPars,"type")) ## through corrHLfit (or fitme !?) or some direct HLCor call
                               && identical(attr(ranPars,"type")$rho,"var") ## can occur in direct call of HLCor 
    )    
    if (corr.model %in% c("adjacency")) { ## "ar1" != "AR1" is a tempo name for a futur generic model  #!# "ar1"
      if ( missing(adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
      ## no nugget in the adjacency model... ## (use additional ranef instead)      
      symSVD <- attr(adjMatrix,"symSVD")
      if (is.null(symSVD) && adj_rho_is_inner_estimated) { ## can occur in direct call of HLCor 
        if (isSymmetric(adjMatrix)) {
          symSVD <- sym_eigen(adjMatrix)
          attr(adjMatrix,"symSVD") <- symSVD
        }             
      }
      if (is.null(symSVD)) {
        corrm <- as.matrix(solve(diag(nrow(adjMatrix))-rho*(adjMatrix))) ## dense for designL.from.corr()
      } else {
        symSVD$adjd <- symSVD$d
        # FR->FR remove $d from symSVD?
      }
    }  else if (corr.model %in% c("SAR_WWt")) { ## "ar1" != "AR1" is a tempo name for a futur generic model  
      if ( missing(adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
      UDU. <- attr(adjMatrix,"UDU.")
      if (is.null(UDU.)) {
        corrm <- as.matrix(solve(diag(nrow(adjMatrix))-rho*(adjMatrix))) ## dense for designL.from.corr()
      } else {
        corrm <- UDU.$u %*% sweep(UDU.$u.,MARGIN=1,1/(1-rho*UDU.$d),`*`) 
      }
      corrm <- tcrossprodCpp(corrm,NULL)
    }  else if (corr.model=="AR1") {
      coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(data))
      uniqueGeo <- calcUniqueGeo(data=data[,coordinates,drop=FALSE])
      txt <- paste(spatial.model[[2]][[3]]) ## the RHS of the ( . | . ) 
      if (length(grep("%in%",txt))>0) {
        stop("HLCor code should be allowed again to handle blockDiag objects")
        #scaled.dist <- as.blockDiag.bar(spatial.model[[2]],formula,data=uniqueGeo)
        #test.in <- TRUE
      } else scaled.dist <- proxy::dist(uniqueGeo)
      corrm <- trueCorrpars$ARphi^scaled.dist
    } else  if (corr.model %in% c("Matern")) {
      txt <- paste(spatial.model[[2]][[3]]) ## the RHS of the ( . | . ) 
      if (length(grep("%in%",txt))>0) {
        stop("(!) Matern( . | <coord> %in% <grp>) is not yet handled.")
        test.in <- TRUE ## should be useful when this case will be handled
      } 
      ## in a typical call from corrHLfit the following test should be FALSE because uniqueGeo and maybe distMatrix should have been precomputed
      if ((length(rho)>1 || missing(distMatrix)) && is.null(uniqueGeo)) { ## all cases where we need uniqueGeo
        coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(data))
        uniqueGeo <- calcUniqueGeo(data=data[,coordinates,drop=FALSE]) ## keeps the names of first instances of the coordinates in data
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
      corrm <- do.call("make_scaled_dist",msd.arglist)
      ## at this point is a single location, corrm should be dist(0) and make_scaled_dist was modified to that effect
      if ( nrow(corrm)>1 ) { ## >1 locations
        norho <- trueCorrpars 
        norho$rho <- NULL ## because the MaternCorr input will be an already scaled distance 'corrm'
        corrm <- do.call(MaternCorr,args=c(norho,list(corrm)))        
      } 
    } else if (corr.model== "corrMatrix") {
      if (missing(corrMatrix)) {
        mess <- pastefrom("missing(corrMatrix) argument despite corrMatrix term in formula.",prefix="(!) From ")
        stop(mess)
      } else {
        dim_corrMatrix <- dim(corrMatrix)
        if (dim_corrMatrix[1L]!=dim_corrMatrix[2L])  stop("corrMatrix is not square") 
      }
      ## ELSE square matrix:
      if (inherits(corrMatrix,"TDT.")) {
        ## stub for Triangular LDLt factorization being provided (Henderson's)
      } else if (inherits(corrMatrix,"dist")) {
          corrnames <- labels(corrMatrix)
        if (is.null(corrnames)) {
          message("corrMatrix without labels: first grouping levels are matched\n  to first rows of corrMatrix, without further check.\n This may cause later errors (notably, wrongly dimensioned matrices) \n See help(\"corrMatrix\") for a safer syntax.")
        }
      } else if (inherits(corrMatrix,"matrix") || inherits(corrMatrix,"Matrix")) {
        corrnames <- rownames(corrMatrix)
        if (is.null(corrnames)) {
          message("corrMatrix without row names: first grouping levels are matched\n  to first rows of corrMatrix, without further check.\n This may cause later errors (notably, wrongly dimensioned matrices) \n See help(\"corrMatrix\") for a safer syntax.")
        }
      } else stop(paste("(!) 'corrMatrix' is neither a 'matrix', 'Matrix' or 'dist' object."))
      whichranef <- which(attr(attr(processed$ZAlist,"ranefs"),"type")=="corrMatrix")
      ZAnames <- colnames(processed$ZAlist[[whichranef]]) ## set by .spMMFactorList(), with two cases for corrMatrix 
      generator <- attr(processed$ZAlist[[whichranef]],"generator")
      if ( length(setdiff(ZAnames,corrnames)) ==0L ) { ## i.e. all ZAnames in corrnames
        ## : should be the case when generator = "as.factor"
        if ( length(corrnames)>length(ZAnames) || any(corrnames!=ZAnames) ) { ## ...but superset, or not same order
          if (inherits(corrMatrix,"dist")) {
            corrm <- as.dist(as.matrix(corrMatrix)[ZAnames,ZAnames]) ## fairly ugly. package 'seriation' has (permute.dist-> C code)
          } else corrm <- corrMatrix[ZAnames,ZAnames]  
        } else corrm <- corrMatrix ## orders already match
      } else {
        ## : expected when generator = "ULI"
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
        corrm <- corrMatrix ## no clear reordering
      }
      # sparse_precision <- spaMM.getOption("sparse_precision")
      # if ( sparse_precision) { ## FIXME : test aussi GLMMbool...
      #   processed <- .reset_processed_sparse_precision(processed)
      #   Qmat <- solve(corrm) ## F I X M E solve <=> devel hack  ## and the resulting fit is wrong
      #   Q_CHMfactor <- Matrix::Cholesky(Matrix::drop0(Qmat),LDL=FALSE,perm=FALSE) 
      #   processed$AUGI0_ZX$envir$Qmat <- Qmat 
      #   processed$AUGI0_ZX$envir$Q_CHMfactor <- Q_CHMfactor ## *existence of this is criterion for use of sparsePrecision algo*
      # }
    } 
    if (verbose["trace"] && length(trueCorrpars)>0) print(unlist(trueCorrpars))
    ## call designL.from.Corr if Lunique not available
    if (is.null(Lunique)) { ## test FR 11/2013 ## modif 2015/04. Noter un calcul de Lunique ci dessus
      if (corr.model=="adjacency" 
          && ! adj_rho_is_inner_estimated 
          && identical(processed$sparse_precision,TRUE)) {
        ## Implies call from fitme_body with outer rho estim.
        Qmat <- - rho * adjMatrix
        diag(Qmat) <- diag(Qmat)+1
        ################################L_Q <- .designL.from.Qmat(Qmat) ## solve(tcrossprod(LMatrix)) = Qmat or
        ## Cholesky gives proper LL' (think LDL')  while chol() gives L'L...
        Q_CHMfactor <- Matrix::Cholesky(Matrix::drop0(Qmat),LDL=FALSE,perm=FALSE) ## FIXME couteux (appele repetitivement)
        processed$AUGI0_ZX$envir$Qmat <- Qmat 
        processed$AUGI0_ZX$envir$Q_CHMfactor <- Q_CHMfactor ## *existence of this is criterion for use of sparsePrecision algo*
        ## fitme sparse_precision has an incomplete symSVD=> corrm not computed, and try(designL.from.Corr(symSVD=symSVD)) fails. Instead:
        Lunique  <- solve(Q_CHMfactor,system="Lt") #  L_Q^{-\top}=LMatrix_correlation
        # Lunique is a dtCMatrix...
      } else {
        ## densePrecision cases with outer rho estimation, and some residula inner rho estimation cases
        if ( ! is.null(symSVD)) {
          symSVD$d <- 1/(1-rho*symSVD$adjd) ## from adjMatrix to correlation matrix
          # outer optim -> LMatrix recomputed from this for each rho  
          Lunique <- try(designL.from.Corr(symSVD=symSVD)) ## usin $d not $adjd
          # Lunique is num, with attr(*, "type")= chr "symsvd" ....
          if ( adj_rho_is_inner_estimated ) { # contrived way of construction Lunique with the correct attributes.
            Lunique[] <- attr(Lunique,"symsvd")$u   ## "[] <- " keeps attributes... except for Matrix...
          }
        } else { ## corrm must exist
          argsfordesignL <- dotlist[intersect(names(dotlist),names(formals(designL.from.Corr)))] 
          if (processed$HL[1L]=="SEM") argsfordesignL$try.chol <- FALSE
          if (inherits(corrm,"dist")) {
            corrm <- as.matrix(corrm)
            diag(corrm) <- 1L ## IF diag missing in input corrMatrix THEN assume a correlation matrix
          } ## else full matrix may be a COV matrix with non-unit diag
          
          Lunique <- try(do.call(designL.from.Corr,c(list(m=corrm),argsfordesignL)))
        }
        if (inherits(Lunique,"try-error")) { 
          print("correlation parameters were:",quote=FALSE) ## makes sense if designL.from.Corr already issued some warning
          print(unlist(trueCorrpars))    
          stop()
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
  HL.info$verbose <- verbose #[intersect(names(verbose),c("warn","trace","summary","SEM"))] 
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
        #attr(ranFix,"types") <- typelist ## used for:
        #parlist <- merge_parlist(,new=ranFix,types=NULL) ## will use attr(ranFix,"types") --> "fix"
        parlist <- merge_parlist(,new=ranFix,types="fix") ## will use attr(ranFix,"types") --> "fix"
        parlist <- merge_parlist(parlist,new=HL.info$init.HLfit,types="var")
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
  hlfit$sparse_precision <- processed$sparse_precision
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

  multiple <- attr(processed,"multiple")
  if ( ( ! is.null(multiple)) && multiple)  { ## "multiple" processed list 
    ## RUN THIS LOOP and return
    fitlist <- lapply(seq_len(length(processed)),function(it){
      locmc <- mc
      locmc[[1L]] <- as.name("HLCor.obj") ## replaces "f" !
      locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
      locmc$processed <- processed[[it]] ## The data are in processed !
      locmc$distMatrix <- mc$distMatrix[[it]] ## but the matrices are not HLfit args hence not in processed ! 
      locmc$uniqueGeo <- mc$uniqueGeo[[it]]
      eval(locmc)
    }) ## a pure list of HLCor objects
    resu <- sum(unlist(fitlist))
    if(mc$verbose["objective"]) {
      unrelist <- unlist(relist(ranefParsVec,skeleton)) ## handles elements of lemgth>1
      cat(" parent HLCor: ",objective,"=",resu,
          " given: ",paste(names(forGiven),"=",signif(unrelist,6),sep="",collapse=", "),
          "\n",sep="")
      # cat(paste(names(forGiven),"=",signif(unrelist,6),sep="",collapse=", "),
      #                               ", ",objective,"=",resu,"\n",sep="")
    }
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
    parlist <- merge_parlist(parlist,new=forGiven,types="fix")## consistently with previous code
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
  if (identical(HLCor.call$verbose["getCall"][[1L]],TRUE)) {return(hlfit)} ## HLCorcall
  #
  aphls <- hlfit$APHLs
  resu <- aphls[[objective]]
  if(mc$verbose["objective"]) {
    unrelist <- unlist(forGiven) ## handles elements of lemgth>1
    cat(" parent HLCor: ",objective,"=",resu,
        " given: ",paste(names(forGiven),"=",signif(unrelist,6),sep="",collapse=", "),
        "\n",sep="")
    # cat(paste(names(forGiven),"=",signif(unrelist,6),sep="",collapse=", "),
    #     ", ",objective,"=",resu,"\n",sep="")
  }
  return(resu) #
}


