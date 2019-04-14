.post_process_parlist <- function(parlist,corr_types,corr_families=NULL) {
  if (is.null(parlist$corrPars)) { ## must be exclusive of other rho, nu...
    corrPars <- list()
    if ( ! is.null(corr_families)) {
      for (rd in seq_along(corr_families)) {
        corr_type <- corr_families[[rd]]$corr_family
        if (! is.null(corr_type)) {
          char_rd <- as.character(rd) ## avoidrho <- parlist$shape <- parlist$longdep <- parlist$Nugget creating NULL list elements 
          parnames <- corr_families[[rd]]$names_for_post_process_parlist
          for (st in intersect(parnames,names(parlist))) corrPars[[char_rd]][st] <- parlist[st] ## creates sublist; does not create NULL element if rhs is NULL
          parlist[parnames] <- NULL
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

.post_process_fixed <- function(fixed,corr_families, hyper_info) {
  if (is.null(fixed)) return(NULL)
  fixed <- .post_process_parlist(fixed,corr_families=corr_families) ## standardize use of corrPars in the parlist
  # I write "F I X" as a TAG for this modif type attribute:
  #attr(fixed,"type") <- relist(rep("fix",length(unlist(fixed))),fixed) ## on veut une list pour pouvoir supp des elements par <- NULL
  # if ( ! is.null(fixed$hyper)) {
  #   for (char_rd in as.character(unlist(hyper_info$ranges))) {
  #     if (is.null(fixed$corrPars[[char_rd]])) fixed$corrPars[[char_rd]] <- list() ## seems no longer necessary for .expand_hyper()
  #   }
  # }
  return(fixed)
}

.calc_moreargs <- function(processed, # possibly a list of environments -> .calc_range_info -> scans then to compute a mean(nbUnique) 
                           corr_types, fixed, init.optim, control_dist, NUMAX=50, LDMAX=50, 
                           KAPPAMAX=100.000001, # so that users can set it to 100...
                           init.HLfit, corr_info, verbose, lower, upper) {
  moreargs <- list()
  for (rd in seq_along(corr_types)) {
    corr_type <- corr_types[rd]
    if (! is.na(corr_type)) {
      char_rd <- as.character(rd)
      adj_test <- (corr_type %in% c("SAR_WWt","adjacency") 
                   &&  is.null(fixed[[char_rd]]$rho) 
                   && (! is.numeric(init.HLfit[[char_rd]]$rho))) ## init.HLfit[[]]$rho NULL or NA) 
      if (adj_test) {
        if (corr_type=="SAR_WWt") decomp <- attr(corr_info$adjMatrices[[rd]],"UDU.")
        if (corr_type=="adjacency") decomp <- attr(corr_info$adjMatrices[[rd]],"symSVD")
      }
      moreargs[[char_rd]] <- corr_info$corr_families[[rd]]$calc_moreargs(fixed=fixed, char_rd=char_rd, init.optim=init.optim, 
                                                                         processed=processed, rd=rd, control_dist=control_dist, 
                                                                         NUMAX=NUMAX, LDMAX=LDMAX, KAPPAMAX=KAPPAMAX,
                                                                         # quite distinct arguments for adjacency:
                                                                         decomp=decomp, verbose=verbose, lower=lower, upper=upper,
                                                                         # IMRF...:
                                                                         IMRF_pars=attr(attr(attr(processed$ZAlist,"exp_spatial_terms")[[rd]],"type"),"pars")
      )
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
                                    needed=c(), geo_envir) {
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
      # activelevels are data-ordered levels whose ranef values affect the likelihood
      # .calcUniqueGeo must produce data-ordered values.
      if (!is.null(geo_envir$activelevels)) uniqueGeo <- uniqueGeo[geo_envir$activelevels,,drop=FALSE]
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
        decomp <- eigen(adjMatrix, symmetric=TRUE)
        svdnames <- names(decomp)
        svdnames[svdnames=="values"] <- "d"
        svdnames[svdnames=="vectors"] <- "u"
        names(decomp) <- svdnames
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
        stop(paste0("(!) '",st,"'" ,errstring))    
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
  return(rho_mapping) ## if length(rho)>1 and input rho.mapping was NULL, output rho.mapping is named 1 2 3... 
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
    return(list(maxrange=maxrange,nbUnique=nbUnique, rho_mapping=rho_mapping)) ## possibly modified (more explicit) rho_mapping
  }
}

.do_TRACE <- function(processed) { ## no need for an 'unTRACE' call at the end of each fit since the next fit will cope with everything.
  ## trouble when called from another package while not attached (bboptim example)
  # THe syntax spaMM::HLfit_body, where=spaMM::fitme does not stop() in that case, 
  #     but HLfit_body is not effectively traced when using spaMM directly attached (standard library(spaMM))
  # The syntax spaMM::HLfit_body without where=also does not trace when using spaMM directly attached
  level <- processed$verbose["TRACE"]
  if ("package:spaMM" %in% search()) {
    if (processed$augZXy_cond) {
      traced_fn <- quote(.HLfit_body_augZXy)
    } else traced_fn <- quote(HLfit_body)
    if (level) {
      # F I X M E ? provide alternative tracer that only prints "." or so.
      suppressMessages(trace(traced_fn, where=asNamespace("spaMM"), print=FALSE, 
                             tracer=quote(try({
                               ranPars <- .canonizeRanPars(ranFix,corr_info=NULL,checkComplete = FALSE)
                               #
                               if ( ! is.null(ranPars$hyper)) {
                                 ranges <- processed$hyper_info$ranges
                                 for (char_hyper_it in names(ranPars$hyper)) {
                                   rd_range <- ranges[[char_hyper_it]]
                                   char_rd_range <- as.character(rd_range)
                                   first_char_rd <- char_rd_range[1L]
                                   ranPars$hyper[[char_hyper_it]] <- list(hy_kap=ranPars$corrPars[[first_char_rd]]$kappa,
                                                                  hy_lam=sum(ranPars$lambda[char_rd_range]))
                                   for (char_rd in char_rd_range) {
                                     ranPars$corrPars[[char_rd]]$kappa <- NULL
                                   }
                                   ranPars$lambda <- ranPars$lambda[setdiff(names(ranPars$lambda),char_rd_range)]
                                 }
                                  if ( ! length(ranPars$lambda)) ranPars$lambda <- NULL
                               }
                               #
                               if (length(ranPars$phi)>1) {
                                 cat("(phi fixed) ")
                                 ranPars[["phi"]] <- NULL
                               }
                               ntC <- names(ranPars)
                               for (lit in seq_along(ranPars)) {
                                 urP <- unlist(ranPars[[lit]]) ## ranPars$corrPars can be list() in which case urP is NULL 
                                 if (!is.null(urP)) cat(ntC[lit],"=",paste(signif(urP,6),collapse=" ")," ")
                               }
                             })),
                             exit=quote({
                               aphl <- unlist(res$APHLs[c("p_bv","p_v","logLapp")])[1L] ## unlist drops NULL values
                               if (is.null(aphl)) {
                                 print("(objective not found)",quote=FALSE)
                               } else print(paste0(names(aphl),"= ",.prettify_num(aphl,nsmall=4)),quote=FALSE)
                             }))) 
      if (processed$sparsePrecisionBOOL) {
        suppressMessages(trace(.solve_IRLS_as_spprec, where=asNamespace("spaMM"),print=FALSE,tracer=quote(cat(">"))))
      } else suppressMessages(trace(.solve_IRLS_as_ZX, where=asNamespace("spaMM"), print=FALSE,tracer=quote(cat(">"))))
      suppressMessages(trace(spaMM.getOption("matrix_method"),print=FALSE,tracer=quote(cat("."))))
      suppressMessages(trace(spaMM.getOption("Matrix_method"),print=FALSE,tracer=quote(cat("."))))
      suppressMessages(trace(spaMM.getOption("spprec_method"),print=FALSE,tracer=quote(cat("."))))
      if (level>3L) {
        fn <- paste("get_from_MME",strsplit(spaMM.getOption("matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_EigenDense_QRP_Chol_scaled
        suppressMessages(trace(fn,print=FALSE,tracer=quote(cat(which))))
        fn <- paste("get_from_MME",strsplit(spaMM.getOption("Matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_Matrix_QRP_CHM_scaled
        suppressMessages(trace(fn,print=FALSE,tracer=quote(cat(which))))
        fn <- paste("get_from_MME",strsplit(spaMM.getOption("spprec_method"),"def_")[[1L]][2],sep=".") #
        suppressMessages(trace(fn,print=FALSE,tracer=quote(cat(which))))
      }
    } else {
      suppressMessages(try(untrace(HLfit_body, where=asNamespace("spaMM")), silent=TRUE))      
      suppressMessages(try(untrace(.HLfit_body_augZXy, where=asNamespace("spaMM")), silent=TRUE))      
      if (processed$sparsePrecisionBOOL) {
        suppressMessages(try(untrace(.solve_IRLS_as_spprec, where=asNamespace("spaMM")), silent=TRUE))
      } else suppressMessages(try(untrace(.solve_IRLS_as_ZX, where=asNamespace("spaMM")), silent=TRUE))
      suppressMessages(untrace(spaMM.getOption("matrix_method")))
      suppressMessages(untrace(spaMM.getOption("Matrix_method")))
      suppressMessages(untrace(spaMM.getOption("spprec_method")))
      fn <- paste("get_from_MME",strsplit(spaMM.getOption("matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_EigenDense_QRP_Chol_scaled
      suppressMessages(untrace(fn))
      fn <- paste("get_from_MME",strsplit(spaMM.getOption("Matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_Matrix_QRP_CHM_scaled
      suppressMessages(untrace(fn))
      fn <- paste("get_from_MME",strsplit(spaMM.getOption("spprec_method"),"def_")[[1L]][2],sep=".") 
      suppressMessages(untrace(fn))
    } 
  } else if (level) {warning("The 'spaMM' package must be *attached* for verbose(TRACE=...) tracing to operate",
                             immediate.=TRUE)}
} 

## creates precisionFactorList if it doesn't exist, and partially fills it with Diagonal()'s
## updates precisionFactorList with LMatrices info
.init_precision_info <- function(processed, LMatrices=NULL) {
  envir <- processed$AUGI0_ZX$envir
  precisionFactorList <- envir$precisionFactorList
  # I N I T
  if (is.null(precisionFactorList)) {
    nranef <- length(envir$finertypes)
    precisionFactorList <- vector("list",nranef) ## will contain diagonal matrices/info for non-trivial (diagonal) precision matrices 
    for (it in seq_len(nranef)) {
      if ( envir$finertypes[it] %in% c("adjacency","corrMatrix","AR1", "IMRF") ) {
        ## terms of these types must be dealt with by ad hoc code for each type elsewhere
      } else if ( envir$finertypes[it] =="(.|.)" ) { # EXcludes "ranCoefs"
        nc <- ncol(processed$ZAlist[[it]])
        precisionFactorList[[it]] <- list(Qmat=.symDiagonal(n=nc), ## used independently of chol_Q_list, see precisionBlocks
                                          chol_Q=new("dtCMatrix",i= 0:(nc-1L), p=0:(nc), Dim=c(nc,nc),x=rep(1,nc)) )
        # All chol_Q's must be dtCMatrix so that bdiag() gives a dtCMatrix
      } else if ( envir$finertypes[it] %in% c("ranCoefs", "Matern", "Cauchy") ) { 
        ## leave precisionFactorList[[it]] NULL
      } else stop(paste("sparse-precision methods were requested, but",envir$finertypes[it],"terms are not yet handled by sparse precision code."))
    }
  }
  # F I L L
  if ( ! is.null(LMatrices)) {
    is_given_by <- attr(LMatrices,"is_given_by")
    for (rt in which(is_given_by %in% c("ranCoefs"))) {
      latentL_blob <- attr(LMatrices[[rt]],"latentL_blob")
      compactchol_Q <- latentL_blob$compactchol_Q
      chol_Q <- .makelong(compactchol_Q,longsize=ncol((LMatrices[[rt]])), 
                          template=processed$ranCoefs_blob$longLv_templates[[rt]] ) 
      if ( ! inherits(chol_Q,"dtCMatrix")) stop("chol_Q should be a 'dtCMatrix'.") ## FIXME check lower tri too
      #
      precisionFactorList[[rt]] <- list(chol_Q=chol_Q,
                                        precmat=.makelong(latentL_blob$compactprecmat,longsize=ncol((LMatrices[[rt]])), 
                                                          template=processed$ranCoefs_blob$longLv_templates[[rt]] ))
    }
    for (rt in which(is_given_by=="inner_ranCoefs")) {
      # : need to initialize because spprec IRLS requires .init_precision_info to have filled all matrices
      nc <- ncol(processed$ranCoefs_blob$longLv_templates[[rt]] )
      precisionFactorList[[rt]] <- list(precmat=.symDiagonal(n=nc), 
                                        chol_Q=new("dtCMatrix",i= 0:(nc-1L), p=0:(nc), Dim=c(nc,nc),x=rep(1,nc)) )
    }
  }
  envir$precisionFactorList <- precisionFactorList
  ## environment modified, no return value
}

.update_precision_info <- function(processed, LMatrices, which.) {
  envir <- processed$AUGI0_ZX$envir
  precisionFactorList <- envir$precisionFactorList
  is_given_by <- attr(LMatrices,"is_given_by")
  for (rt in which(is_given_by %in% which.)) {
    if (is_given_by[rt] =="inner_ranCoefs") { 
      latentL_blob <- attr(LMatrices[[rt]],"latentL_blob")
      compactchol_Q <- latentL_blob$compactchol_Q
      chol_Q <- .makelong(compactchol_Q,longsize=ncol((LMatrices[[rt]])), 
                          template=processed$ranCoefs_blob$longLv_templates[[rt]] ) 
      if ( ! inherits(chol_Q,"dtCMatrix")) stop("chol_Q should be a 'dtCMatrix'.") ## FIXME check lower tri too
      #
      precisionFactorList[[rt]] <- list(chol_Q=chol_Q,
                                        precmat=.makelong(latentL_blob$compactprecmat,longsize=ncol((LMatrices[[rt]])), 
                                                          template=processed$ranCoefs_blob$longLv_templates[[rt]] ))
    } else stop("'which.' value not handled ")
  }
  envir$precisionFactorList <- precisionFactorList
  ## environment modified, no return value
}

.init_optim_phi <- function(phi.Fix, phimodel, processed, init.optim, nrand) {
  if (is.null(phi.Fix)) {
    if (phimodel == "phiScal") {
      if (NROW(processed$y)>200L || 
          processed$cum_n_u_h[length(processed$cum_n_u_h)]>200 || # leverage computation also costly in that case (FIXME single criterion ?) 
          processed$family$family=="Gamma") { # FIXME: processed ? (more occurrences below)
        ## default method is outer but can be reversed if init is NaN (and to test NaN we need to test NULL)
        not_inner_phi <- is.null(init.optim$phi) || ! is.nan(init.optim$phi) ## outer if NULL, NA or numeric
      } else {
        ## default is inner but can be reversed if numeric or NA (hence neither NULL nor NaN)
        not_inner_phi <- ! (is.null(init.optim$phi) || is.nan(init.optim$phi))
      }
  } else not_inner_phi <- FALSE
    if (not_inner_phi) {
      ## .get_init_phi takes more account of ranefs than .get_inits_by_glm(.)$phi_est, but it may return NULL.
      #  It seems (?) marginally better hence we try it first and check the result.
      if (is.null(init.optim$phi)) init.optim$phi <- .get_init_phi(processed,weights=eval(processed$prior.weights)) 
      ## (fixme): fitme may fail obscurely if .get_init_phi(processed) fails silently
      if (is.null(init.optim$phi)) { ## .get_init_phi returned NULL if ./. 
        # ./. no replicate obs is available, or
        # ./. ULI too large ?? Grey area here (fixme)
        init.optim$phi <- .get_inits_by_glm(processed)$phi_est/(nrand+1L) ## at least one initial value should represent high guessed variance
        # if init.optim$phi too low (as in min(.,2)) then fitme(Reaction ~ Days + AR1(1|Days) + (Days|Subject), data = sleepstudy) is poor
      }  
    }
  } else not_inner_phi <- TRUE ## not outer as well. Hence not_inner_phi != outer phi...
  return(list(not_inner_phi=not_inner_phi, init.optim=init.optim))
}

.init_optim_lambda_ranCoefs <- function(proc1, not_inner_phi, init.optim, nrand, ranCoefs_blob, var_ranCoefs) {
  lFix <- proc1$lambda.Fix ## only for 'simple" ranefs with Xi_cols=1
  if (proc1$augZXy_cond || 
      not_inner_phi) { ## Tests show it is very inefficient to use outer optim on lambda (at least) when phi must be inner optimized
    optim_lambda_with_NAs <- .reformat_init_lambda_with_NAs(init.optim$lambda, nrand=nrand, default=NA)
    ## handling fitme call for resid fit with meanfit-optimized parameters (if input is NULL, output is all NA):
    optim_resid_lambda_with_NAs <- .reformat_init_lambda_with_NAs(proc1$envir$ranPars$lambda, nrand=nrand, default=NA)
    which_NA_simplelambda <- which(is.na(lFix) & 
                                     is.na(optim_resid_lambda_with_NAs) & 
                                     (is.na(optim_lambda_with_NAs) & ! is.nan(optim_lambda_with_NAs)) & ## explicit NaN's will be inner-optimized
                                     ! ranCoefs_blob$isRandomSlope) ## exclude random slope whether set or not
    if (length(which_NA_simplelambda)) { 
      init_lambda <- .eval_init_lambda_guess(proc1, stillNAs=which_NA_simplelambda, For="optim")
      optim_lambda_with_NAs[which_NA_simplelambda] <- init_lambda[which_NA_simplelambda]
      init.optim$lambda <- optim_lambda_with_NAs[ ! is.na(optim_lambda_with_NAs)] ## NaN now rmoved if still there (cf is.na(c(1,NA,NaN))) BUT
      # ... it's dubious that we have augZXy_cond || not_inner_phi if we requested inner estimation of lambda by a NaN                                                  
    }
  } ## else use inner optimization  for simple lambdas if inner_phi is necessary
  init.optim$lambda <- init.optim$lambda[ ! is.nan(init.optim$lambda)] ## removes users's explicit NaN, which effect is documented in help(fitme)
  #
  if (any(var_ranCoefs)) {
    # Not super clear why I considered nranterms (user level ??) instead of nrand. FIXME.
    nranterms <- sum(var_ranCoefs | is.na(lFix)) ## var_ranCoefs has FALSE elements for non-ranCoefs (hence it is full-length)
    guess_from_glm_lambda <- .get_inits_by_glm(proc1)$lambda * (3L*nranterms)/((nranterms+1L)) # +1 for residual
    fam_corrected_guess <- .calc_fam_corrected_guess(guess=guess_from_glm_lambda, For="optim", processed=proc1) ## divides by nrand...
    for (rt in which(var_ranCoefs)) {
      char_rt <- as.character(rt)
      if (is.null(init.optim$ranCoefs[[char_rt]])) {
        Xi_cols <- attr(proc1$ZAlist,'Xi_cols')
        Xi_ncol <- Xi_cols[rt]
        rc <- rep(0,Xi_ncol*(Xi_ncol+1L)/2L)
        #
        if (FALSE) { # guess to find good initial value, but no useful impact...
          # gaussian at least: fits the model y~Zv and uses the corr of the v's...
          ZAlist <- proc1$ZAlist
          ZA <- .compute_ZAL(NULL,proc1$ZAlist[rt], as_matrix=FALSE)
          coef <-.lmwith_sparse_QRp(ZA,1.0*(as.numeric(proc1$y)*1.0-proc1$off),returntQ = FALSE,returnR = TRUE)$coef 
          cor <- cov2cor(cov(matrix(coef,ncol=Xi_ncol)))
          rc <- cor[lower.tri(cor,diag = TRUE)]
        }
        #
        lampos <- rev(length(rc) -cumsum(seq(Xi_ncol))+1L)  ## NOT cumsum(seq(Xi_cols))
        rc[lampos] <- fam_corrected_guess/(Xi_ncol)
        init.optim$ranCoefs[[char_rt]] <- rc ## see help(ranCoefs)
      }
    }
  }
  return(init.optim)
}


  .more_init_optim <- function(proc1, corr_types, init.optim) {
  phimodel <- proc1$models[['phi']]
  #if (phimodel=="phiGLM") {message("'fitme' not optimized for models with structured dispersion.")} ## FR->FR test this later"
  ranCoefs_blob <- proc1$ranCoefs_blob
  ## trying to guess all cases where optimization is useful. But FIXME: create all init and decide afterwardsS
  if ( (is_MixedM <- ( ! is.null(ranCoefs_blob) )) && (
    any(var_ranCoefs <- (ranCoefs_blob$isRandomSlope & ! ranCoefs_blob$is_set)) ||
    phimodel == "phiScal" || # lambda + phi
    var(proc1$y)<1e-3 || # mixed model with low response variance (including case where phi is fixed (phimodel="") )
    (phimodel == "phiHGLM" && identical(spaMM.getOption("outer_optim_resid"),TRUE)) || 
    length(var_ranCoefs)>1L || # +s lambda
    length(corr_types[ ! is.na(corr_types)]) # lambda + corr pars
  )
  ) { ## All cases where phi is fixed, or can be outer optimized jointly with lambda:  
    ## prefer outer optimization if var_ranCoefs or not resid.model
    nrand <- length(proc1$ZAlist)
    # (1) set (or not) outer optimization for phi: 
    phi.Fix <- proc1$phi.Fix
    ## Then determine a single phi value
    init_optim_phi_blob <- .init_optim_phi(phi.Fix, phimodel, proc1, init.optim, nrand)
    not_inner_phi <- init_optim_phi_blob$not_inner_phi
    init.optim <- init_optim_phi_blob$init.optim
    # (2) set outer optimization for lambda and ranCoefs (handling incomplete ranFix$lambda vectors)
    if (is_MixedM) {
      init.optim <- .init_optim_lambda_ranCoefs(proc1, not_inner_phi, init.optim, nrand, ranCoefs_blob, var_ranCoefs)
    }
  }
  return(init.optim)
}

.calc_optim_args <- function(proc1, processed, init, fixed, lower, upper, verbose, optim.scale, For) {
  corr_info <- proc1$corr_info 
  # modify HLCor.args and <>bounds;   ## distMatrix or uniqueGeo potentially added to HLCor.args:
  corr_types <- corr_info$corr_types
  fixed <- .post_process_fixed(fixed, corr_families=corr_info$corr_families, processed$hyper_info)
  if (processed$models[["phi"]]!="phiScal" && ! is.null(init$phi)) {
    warning("initial value for 'phi' is ignored when the is a non-default resid.model") # i.e. anything but Intercept model
    init$phi <- NULL # otherwise the residModel would be ignored!
  }
  init.optim <- .post_process_parlist(init, corr_families=corr_info$corr_families)
  init.HLfit <- proc1$init_HLfit #to be modified below ## dhglm uses fitme_body not (fitme-> .preprocess) => dhglm code modifies processed$init_HLfit
  init <- NaN ## make clear it's not to be used 
  #
  ##### init.optim$phi/lambda will affect calc_inits -> calc_inits_dispPars.
  # outer estim seems useful when we can suppress all inner estim (thus the hatval calculations). 
  # ./. Therefore, we need to identify all cases where phi is fixed, 
  # ./. or can be outer optimized jointly with lambda:  
  # Provide ranCoefs inits by calling .init_optim_lambda_ranCoefs:
  if (For=="fitme" && proc1$HL[1]!="SEM") { ## for SEM, it's better to let SEMbetalambda find reasonable estimates
    init.optim <- .more_init_optim(proc1=proc1, corr_types=corr_types, init.optim=init.optim) # This decides for outer/inner optimisations
    if (proc1$augZXy_cond) init.optim$phi <- NULL
  }
  #
  family <- proc1$family
  if (family$family=="COMPoisson") {
    checknu <- suppressWarnings(try(environment(family$aic)$nu,silent=TRUE))
    if (inherits(checknu,"try-error") && is.null(init.optim$COMP_nu)) init.optim$COMP_nu <- 1 ## FR->FR FIXME what if not included in the range ?
  } else if (family$family == "negbin") {
    checktheta <- suppressWarnings(try(environment(family$aic)$shape,silent=TRUE))
    if (inherits(checktheta,"try-error") && is.null(init.optim$NB_shape)) init.optim$NB_shape <- 1 ## FR->FR FIXME idem ...
  }
  user.lower <- .post_process_parlist(lower,corr_families=corr_info$corr_families)
  user.upper <- .post_process_parlist(upper,corr_families=corr_info$corr_families) ## keep user input 
  .check_conflict_init_fixed(fixed,init.optim, "given as element of both 'fixed' and 'init'. Check call.")
  .check_conflict_init_fixed(init.HLfit,init.optim, "given as element of both 'init.HLfit' and 'init'. Check call.") ## has quite poor effect on fits
  moreargs <- .calc_moreargs(processed=processed, # possibly a list of environments -> .calc_range_info -> scans then to compute a mean(nbUnique) 
                             corr_types=corr_types, fixed=fixed, init.optim=init.optim, control_dist=proc1$control_dist, 
                             init.HLfit=init.HLfit, corr_info=corr_info, verbose=verbose, lower=lower, upper=upper)
  fixed <- .expand_hyper(fixed, hyper_info=processed$hyper_info, moreargs=moreargs)
  inits <- .calc_inits(init.optim=init.optim, init.HLfit=init.HLfit,
                       ranFix=fixed,  corr_info=corr_info,
                       moreargs=moreargs,
                       user.lower=user.lower, user.upper=user.upper,
                       optim.scale=optim.scale, 
                       For=For, hyper_info=processed$hyper_info
  )
  if ("lambda" %in% c(names(user.lower),names(user.lower)) 
      && is.null(inits$init$lambda)) {
    stop("'lambda' in 'lower' or 'upper' has no effect if absent from 'init'.")
  }
  ################
  LUarglist <- list(canon.init=inits$`init`, 
                    init.optim=inits$init.optim,
                    user.lower=user.lower,user.upper=user.upper,
                    corr_types=corr_types,
                    ranFix=fixed, # inits$ranFix, # Any change in $ranFix would be ignored 
                    optim.scale=optim.scale, 
                    moreargs=moreargs) ## list needed as part of attr(,"optimInfo")
  ################
  LowUp <- do.call(".makeLowerUpper",LUarglist)
  ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  return(list(inits=inits, moreargs=moreargs, 
              fixed=fixed, corr_types=corr_types, LUarglist=LUarglist,LowUp=LowUp))
} 

