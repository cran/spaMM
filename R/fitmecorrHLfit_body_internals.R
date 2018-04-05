.post_process_parlist <- function(parlist,corr_types) {
  if (is.null(parlist$corrPars)) { ## must be exclusive of other rho, nu...
    corrPars <- list()
    for (rd in seq_along(corr_types)) {
      corr_type <- corr_types[rd]
      if (! is.na(corr_type)) {
        char_rd <- as.character(rd) ## avoid creating NULL list elements 
        if (corr_type == "Matern") {
          corrPars[[char_rd]]$rho <- parlist$rho ## does not create NULL element if rhs is NULL
          corrPars[[char_rd]]$nu <- parlist$nu
          corrPars[[char_rd]]$Nugget <- parlist$Nugget
          parlist$rho <- parlist$nu <- parlist$Nugget <- NULL
        } else if (corr_type =="Cauchy") {
          corrPars[[char_rd]]$rho <- parlist$rho ## does not create NULL element if rhs is NULL
          corrPars[[char_rd]]$shape <- parlist$shape
          corrPars[[char_rd]]$longdep <- parlist$longdep
          corrPars[[char_rd]]$Nugget <- parlist$Nugget
          parlist$rho <- parlist$shape <- parlist$longdep <- parlist$Nugget <- NULL
        } else if (corr_type=="AR1") {
          corrPars[[char_rd]]$ARphi <- parlist$ARphi
          parlist$ARphi <- NULL
        } else if (corr_type=="adjacency") {
          corrPars[[char_rd]]$rho <- parlist$rho
          parlist$rho <- NULL
        }
      }
    }
    parlist$corrPars <- corrPars
  }
  return(parlist)
}

.reformat_init_lambda_with_NAs <- function(init_lambda, nrand, default=NA) {
  if ( ! is.null(init_lambda)) {
    if (length(init_lambda)==nrand) {
      ## According to help(ranPars), the lambda's should be indexed "1",... and presumably init.optim to
      # but they may also be taken from a previous fit and will have complex names (eg test-dhglm) 
      names(init_lambda) <- seq_len(nrand) ## erases complex names
    } else { 
      lambda <- structure(rep(default,nrand),names=seq_len(nrand)) ## default=NA implies that fitme will outer optimize them
      lambda[names(init_lambda)] <- init_lambda ## the lambda's should be indexed "1",...
      if (length(lambda)>nrand) stop("length(lambda)>nrand: case not handled")
      init_lambda <- lambda
    }
  } else init_lambda <- structure(rep(default,nrand),names=seq_len(nrand))
  return(init_lambda)
}

.post_process_fixed <- function(fixed,corr_types) {
  if (is.null(fixed)) return(NULL)
  fixed <- .post_process_parlist(fixed,corr_types)
  attr(fixed,"type") <- relist(rep("fix",length(unlist(fixed))),fixed) ## on veut une list pour pouvoir supp des elements par <- NULL
  return(fixed)
}

.calc_moreargs <- function(processed, # possibly a list of environments -> .calc_range_info -> scans then to compute a mean(nbUnique) 
                           corr_types, fixed, init.optim, control_dist, NUMAX=50, LDMAX=50,
                           init.HLfit, corr_info, verbose, lower, upper) {
  moreargs <- list()
  for (rd in seq_along(corr_types)) {
    corr_type <- corr_types[rd]
    if (! is.na(corr_type)) {
      char_rd <- as.character(rd)
      if (corr_type %in% c("Matern","Cauchy")) {
        if (is.null(fixed$corrPars)) { ## before spaMM3.0
          rho_size <- max(length(fixed$rho),length(init.optim$rho))
        } else { ## spaMM3.0
          rho_size <- max(length(fixed$corrPars[[char_rd]][["rho"]]),length(init.optim$corrPars[[char_rd]][["rho"]]))
        }
        range_info_blob <- .calc_range_info(rho_size, processed, rd, control_dist[[char_rd]]) 
        RHOMAX <- 1.1*30*range_info_blob$nbUnique/range_info_blob$maxrange## matches init$rho in calc_inits() 
        control_dist[[char_rd]]$rho.mapping <- range_info_blob$rho_mapping 
        moreargs[[char_rd]] <- list(maxrange=range_info_blob$maxrange, RHOMAX=RHOMAX,
                                    NUMAX=NUMAX, LDMAX=LDMAX,## not variable...
                                    nbUnique=range_info_blob$nbUnique, control.dist=control_dist[[char_rd]]## needed for LUarglist, not calc_inits
                                    ) 
      }
      if (is.null(fixed$corrPars)) { ## before spaMM3.0
        adj_test <- (corr_type %in% c("SAR_WWt","adjacency") 
                     &&  is.null(.getPar(fixed,"rho")) 
                     && (! is.numeric(init.HLfit$rho))) ## init.HLfit$rho NULL or NA) 
      } else { ## spaMM3.0
        adj_test <- (corr_type %in% c("SAR_WWt","adjacency") 
                     &&  is.null(fixed[[char_rd]]$rho) 
                     && (! is.numeric(init.HLfit[[char_rd]]$rho))) ## init.HLfit[[]]$rho NULL or NA) 
      }
      if (adj_test) {
        if (corr_type=="SAR_WWt") decomp <- attr(corr_info$adjMatrices[[rd]],"UDU.")
        if (corr_type=="adjacency") decomp <- attr(corr_info$adjMatrices[[rd]],"symSVD")
        rhorange <- sort(1/range(decomp$d)) ## keeping in mind that the bounds can be <>0
        if(verbose["SEM"])  cat(paste("Feasible rho range: ",paste(signif(rhorange,6),collapse=" -- "),"\n"))
        if (is.null(lower$corrPars)) { ## TRUE before spaMM3.0
          rhorange[1L] <- max(rhorange[1L],lower$rho)
          rhorange[2L] <- min(rhorange[2L],upper$rho)
        } else { ## spaMM3.0
          rhorange[1L] <- max(rhorange[1L],lower$corrPars[[char_rd]][["rho"]])
          rhorange[2L] <- min(rhorange[2L],upper$corrPars[[char_rd]][["rho"]])
        }
        moreargs[[char_rd]] <- list(rhorange=rhorange) # for spaMM3.0
      }
    }
  }
  return(moreargs)
}


.get_coordinates <- function(spatial_term, data) {
  if ( is.null(spatial_term)) {
    stop("An obsolete syntax for the adjacency model appears to be used.")
    ## coordinates <- c("x","y") ## backward compatibility... old syntax with (1|pos) and default values of the coordinates argument
  } else {
    if ( inherits(data,"list")) {
      dataForCheck <- data[[1]]
    } else dataForCheck <- data
    coordinates <- .extract_check_coords(spatial_term=spatial_term,datanames=names(dataForCheck))
  }
}

.get_dist_nested_or_not <- function(spatial_term, data, distMatrix, uniqueGeo, dist.method, as_matrix=FALSE,
                                    needed=c()) {
  if (is.na(needed["distMatrix"])) needed["distMatrix"] <- FALSE 
  if (is.na(needed["uniqueGeo"])) needed["uniqueGeo"] <- FALSE 
  if (is.na(needed["nbUnique"])) needed["nbUnique"] <- FALSE 
  coordinates <- .get_coordinates(spatial_term=spatial_term, data=data)
  geoMats <- list(coordinates=coordinates)
  if (is.null(distMatrix)) { 
    coord_within <- .extract_check_coords_within(spatial_term=spatial_term) 
    coords_nesting <- setdiff(coordinates,coord_within)
    coord_within <- setdiff(coordinates, coords_nesting)
    ## needed in most cases for evaluation of further elements:
    if (is.null(uniqueGeo) ) { ## then construct it from the data ## this should be the routine case, except for AR1
      uniqueGeo <- .calcUniqueGeo(data=data[,coord_within,drop=FALSE])
    } 
    geoMats$nbUnique <- nrow(uniqueGeo) ## only serves to control RHOMAX and should not be computed from the final uniqueGeo in case of nesting
    if (length(coords_nesting) && any(needed[c("distMatrix","uniqueGeo")]) ) {
      ## should now (>2.3.9) work on data.frames, nesting factor not numeric.
      e_uniqueGeo <- .calcUniqueGeo(data=data[,coordinates,drop=FALSE])
      ## The rows_bynesting <- by(e_uniqueGeo ,e_uniqueGeo[,coords_nesting],rownames) line in.expand_GeoMatrices()
      #  works only if the cols are not factors ! (an as.matrix() ?). Unless we have a fix for this,
      #  e_uniqueGeo classes should not be factor; and as.integer(<factor>) can produce NA's hence as.character()
      # same operation was performed to generate uniqueGeo in .calc_AR1_sparse_Q_ranges()
      for (nam in coords_nesting) {if (is.factor(fac <- e_uniqueGeo[[nam]])) e_uniqueGeo[[nam]] <- as.character(levels(fac))[fac]}
      eGM <- .expand_GeoMatrices(w_uniqueGeo=uniqueGeo, e_uniqueGeo=e_uniqueGeo, 
                                 coords_nesting=coords_nesting, coord_within=coord_within,dist.method=dist.method)
      distMatrix <- eGM$distMatrix 
      uniqueGeo <- e_uniqueGeo
    } else if (needed["distMatrix"]) {
      notnumeric <- ! unlist(lapply(uniqueGeo,is.numeric))
      if (any(notnumeric)) stop(paste(paste(names(which(notnumeric)),collapse=" "),
                                      "are not numeric, hence not suitable for computation of distance matrix."))
      distMatrix <- proxy::dist(uniqueGeo,method=dist.method)
      if (as_matrix) distMatrix <- as.matrix(distMatrix) ## useful to select rows or cols in predict() 
    }
    if (needed["distMatrix"]) geoMats$distMatrix <- distMatrix ## not always computed
    geoMats$uniqueGeo <- uniqueGeo ## always computed (needed for distMatrix)
  } else { ## there is a distMatrix, this is what will be used by HLCor
    if (needed["nbUnique"]) geoMats$nbUnique <- .checkDistMatrix(distMatrix,data,coordinates)
  }
  return(geoMats)
}

.provide_AR_factorization <- function(adjMatrix, sparse_precision, corr.model) {
  if (corr.model  %in% c("SAR_WWt")) {
    decomp <- eigen(adjMatrix,symmetric=FALSE) ## FIXME not RcppEigen-optimized
    return(list(u=decomp$vectors,d=decomp$values,u.=solve(decomp$vectors)))
  }
  # ELSE
  
  if (corr.model  %in% c("adjacency")) {
    ## eigenvalues needed in all cases for the bounds. Full decomp not always needed
    if (isSymmetric(adjMatrix)) { ## FALSE if row and col names differ !
      if ( sparse_precision) { 
        decomp <- list(d=eigen(adjMatrix,only.values = TRUE)$values) ## only eigenvalues
      } else {
        decomp <- sym_eigen(adjMatrix)
      }
      return(decomp)
    } else stop("'adjMatrix' is not symmetric") ## => invalid cov mat for MVN
  }

}

.check_conflict_init_fixed <- function(fixed, init, errstring) {
  fixed_cP <- fixed$corrPars
  init_cP <- init$corrPars
  for (char_rd in unique(names(fixed_cP),names(init_cP))) {
    fixed_rd <- fixed_cP[[char_rd]]
    init_rd <- init_cP[[char_rd]]
    for (st in c("nu","ARphi","Nugget","rho")) {
      if ( (! is.null(fixed_rd[[st]])) && (! is.null(init_rd[[st]])) ) {
        stop(paste("(!) '",st,"'" ,errstring,sep=""))    
      }
    }
  }
}



.provide_rho_mapping <- function(control.dist, coordinates, rho.size) {
  rho_mapping <- control.dist$rho.mapping ## may be NULL
  if (is.null(rho_mapping) ) { ## && length(coordinates)>1L ?
    if (length(coordinates)==rho.size) { ## rho.size comes from explicit rho from user
      rho_mapping <- seq_len(rho.size)           
      names(rho_mapping) <- coordinates
    } else if (length(rho.size)>1L) stop("'rho.mapping' missing with no obvious default from the other arguments.")
  } ## then (for given corr.model's) there is rho_mapping
  return(rho_mapping)
}

.calc_range_info <- function(rho.size, processed, it, control_dist_rd) {
  if (rho.size<2) { ## can be 0 if no explicit rho in the input  
    if (is.list(processed)) {
      dist.method <- control_dist_rd$dist.method 
      locdist <- vector("list",length(processed))
      maxs <- numeric(length(processed))
      mins <- numeric(length(processed))
      nbUnique <- 0L
      for (lit in seq_along(processed)) {
        geo_envir <- .get_geo_info(processed[[lit]], which_ranef=it, which=c("distMatrix","uniqueGeo","nbUnique"), 
                                   dist.method=dist.method) ## this is all for ranef [[it]]:
        maxs[lit] <- max(c(-Inf,geo_envir$distMatrix)) ## les Inf to handle dist(0)..
        mins[lit] <- min(c(Inf,geo_envir$distMatrix)) ## list over proc envirs not over ranefs!
        nbUnique <- nbUnique + geo_envir$nbUnique
      }
      maxrange <- max(maxs)-min(mins)
      nbUnique <- nbUnique/length(processed)
    } else {
      geo_envir <- .get_geo_info(processed, which_ranef=it, which=c("distMatrix","uniqueGeo","nbUnique"), 
                                 dist.method=control_dist_rd$dist.method) ## this is all for ranef [[it]]:
      ## => assuming, if (is.list(processed)), the same matrices accross data for the given ranef term.
      locdist <- geo_envir$distMatrix   ## specific to the it'th ranef
      maxrange <- max(geo_envir$distMatrix)-min(geo_envir$distMatrix)
      nbUnique <- geo_envir$nbUnique
    }
    return(list(maxrange=maxrange,nbUnique=nbUnique))
  } else { 
    if (is.list(processed)) { ## rho.size >1 FIXME remove the local functions...
      dist.method <- control_dist_rd$dist.method 
      uniqueGeo <- vector("list",length(processed))
      nbUnique <- 0L
      for (lit in seq_along(processed)) {
        geo_envir <- .get_geo_info(processed[[lit]], which_ranef=it, which=c("distMatrix","uniqueGeo","nbUnique"), 
                                   dist.method=dist.method) ## this is all for ranef [[it]]:
        uniqueGeo[[lit]] <- geo_envir$uniqueGeo   
        nbUnique <- nbUnique + geo_envir$nbUnique
      }
      rho_mapping <- .provide_rho_mapping(control_dist_rd, geo_envir$coordinates, rho.size) ## using the last geo_envir
      maxrange <- lapply(unique(rho_mapping), function(idx) {
        ranges <- matrix(unlist(lapply(uniqueGeo, function(uu){
          if (nrow(uu)>1) {
            range(proxy::dist(uu[,rho_mapping==idx],method=dist.method))
          } else c(Inf,-Inf) ## encore des Inf to handle dist(0)...
        })),ncol=2)
        max(ranges[,2])-min(ranges[,1]) 
      })
      nbUnique <- nbUnique/length(processed)
    } else { 
      geo_envir <- .get_geo_info(processed, which_ranef=it, which=c("distMatrix","uniqueGeo","nbUnique"), 
                                 dist.method=control_dist_rd$dist.method) ## this is all for ranef [[it]]:
      ## => assuming, if (is.list(processed)), the same matrices accross data for the given ranef term.
      rho_mapping <- .provide_rho_mapping(control_dist_rd, geo_envir$coordinates, rho.size)
      u_rho_mapping <- unique(rho_mapping)
      maxrange <- numeric(length(u_rho_mapping))
      for (uit in u_rho_mapping) {
        idx <- u_rho_mapping[uit]
        maxrange[uit] <- diff(range(proxy::dist(geo_envir$uniqueGeo[,rho_mapping==idx],method=control_dist_rd$dist.method)))
      }
      nbUnique <- geo_envir$nbUnique
    }
    maxrange <- unlist(maxrange)
    return(list(maxrange=maxrange,nbUnique=nbUnique, rho_mapping=rho_mapping))
  }
}

.do_TRACE <- function(level) {
  if (level) {
    suppressMessages(trace(HLfit_body,print=FALSE, 
                           tracer=quote(try({
                             ranPars <- .canonizeRanPars(ranFix,corr_types=NULL,checkComplete = FALSE)
                             ntC <- names(ranPars)
                             for (lit in seq_along(ranPars)) {
                               urP <- unlist(ranPars[[lit]]) ## ranPars$corrPars can be list() in which case urP is NULL 
                               if (!is.null(urP)) cat(ntC[lit],"=",paste(signif(urP,6),collapse=" ")," ")
                             }
                           })),
                           exit=quote({
                             aphl <- unlist(res$APHLs[c("p_bv","p_v")])[1L] ## unlist drops a NULL p_bv
                             print(paste(names(aphl),"= ",.prettify_num(aphl,nsmall=4),sep=""),quote=FALSE)
                           }))) 
    suppressMessages (trace(.solve_IRLS_as_ZX,print=FALSE,tracer=quote(cat(">"))))
    suppressMessages (trace(.solve_IRLS_as_spprec,print=FALSE,tracer=quote(cat(">"))))
    suppressMessages (trace(spaMM.getOption("matrix_method"),print=FALSE,tracer=quote(cat("."))))
    suppressMessages(trace(spaMM.getOption("Matrix_method"),print=FALSE,tracer=quote(cat("."))))
    suppressMessages(trace("def_AUGI0_ZX_sparsePrecision",print=FALSE,tracer=quote(cat("."))))
    if (level>3L) {
      fn <- paste("get_from_MME",strsplit(spaMM.getOption("matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_EigenDense_QRP_Chol_scaled
      suppressMessages(trace(fn,print=FALSE,tracer=quote(cat(which))))
      fn <- paste("get_from_MME",strsplit(spaMM.getOption("Matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_Matrix_QRP_CHM_scaled
      suppressMessages(trace("fn",print=FALSE,tracer=quote(cat(which))))
      suppressMessages(trace("get_from_MME.AUGI0_ZX_sparsePrecision",print=FALSE,tracer=quote(cat(which))))
    }
  } else {
    suppressMessages(untrace(spaMM::HLfit_body))
    suppressMessages(try(untrace(.solve_IRLS_as_ZX),silent=TRUE)) ## untracing untraced internal functions fails
    suppressMessages(try(untrace(.solve_IRLS_as_spprec),silent=TRUE))
    suppressMessages(untrace(spaMM.getOption("matrix_method")))
    suppressMessages(untrace(spaMM.getOption("Matrix_method")))
    suppressMessages(untrace("def_AUGI0_ZX_sparsePrecision"))
    fn <- paste("get_from_MME",strsplit(spaMM.getOption("matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_EigenDense_QRP_Chol_scaled
    suppressMessages(untrace(fn))
    fn <- paste("get_from_MME",strsplit(spaMM.getOption("Matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_Matrix_QRP_CHM_scaled
    suppressMessages(untrace(fn))
    suppressMessages(untrace("get_from_MME.AUGI0_ZX_sparsePrecision"))
  } 
} 

## creates precisionFactorList if it doesn't exist, and partially fills it with Diagonal()'s
## updates precisionFactorList with LMatrices info
.init_precision_info <- function(processed, LMatrices=NULL) {
  envir <- processed$AUGI0_ZX$envir
  precisionFactorList <- envir$precisionFactorList
  if (is.null(precisionFactorList)) {
    nranef <- length(envir$finertypes)
    precisionFactorList <- vector("list",nranef) ## will contain diagonal matrices/info for non-trivial (diagonal) precision matrices 
    for (it in seq_len(nranef)) {
      if ( envir$finertypes[it] %in% c("adjacency","corrMatrix","AR1") ) {
        ## terms of these types must have been dealt with by ad hoc code for each type
      } else if ( envir$finertypes[it] =="(.|.)" ) { # EXcludes "ranCoefs"
        nc <- ncol(processed$ZAlist[[it]])
        precisionFactorList[[it]] <- list(Qmat=Diagonal(n=nc), ## used independently of chol_Q_list, see precisionBlocks
                                          chol_Q=new("dtCMatrix",i= 0:(nc-1L), p=0:(nc), Dim=c(nc,nc),x=rep(1,nc)) )
        # All chol_Q's must be dtCMatrix so that bdiag() gives a dtCMatrix
      } else if ( envir$finertypes[it] =="ranCoefs" ) { 
        ## leave precisionFactorList[[it]] NULL
      } else stop(paste("sparse precision was selected, but",envir$finertypes[it],"terms are not yet handled by sparse precision code."))
      # Matern fails with Error in solve(Q_CHMfactor, system = "Lt") : object 'Q_CHMfactor' not found 
    }
  }
  if ( ! is.null(LMatrices)) {
    is_given_by <- attr(LMatrices,"is_given_by")
    for (it in which(is_given_by=="ranCoefs")) {
      latentL_blob <- attr(LMatrices[[it]],"latentL_blob")
      compactchol_Q <- latentL_blob$compactchol_Q
      chol_Q <- .makelong(compactchol_Q,longsize=ncol((LMatrices[[it]])) ) 
      if ( ! inherits(chol_Q,"dtCMatrix")) stop("chol_Q should be a 'dtCMatrix'.") ## FIXME check lower tri too
      #
      precisionFactorList[[it]] <- list(chol_Q=chol_Q,
                                  precmat=.makelong(latentL_blob$compactprecmat,longsize=ncol((LMatrices[[it]])) ))
    }
  }
  envir$precisionFactorList <- precisionFactorList
  ## environment modified, no return value
}

