.reformat_lambda <- function(user_lFix, nrand, namesTerms=NULL, full_lambda) {
  seq_nrand <- seq_len(nrand)
  if (full_lambda) { # 
    template <- rep(NA,nrand)
    names(template) <- seq_nrand # (character)
  } else if (is.null(user_lFix)) { # NULL input -> NULL output
    return(NULL)
  } else {
    template <- user_lFix # reformats but not full-length; do not introduce (nor remove) NA/NaN's
    if (is.null(names(template))) names(template) <- seq_along(template)
  }
  if ( ! is.null(user_lFix)) {
    user_names <- names(user_lFix)
    if (length(unique(user_names))!=length(user_names)) stop("Repeated names in names of user-specified lambda. Check your input.")
    # 'user'_names are not necessarily 1,2... (fitted values have names from rhs of ranef terms) 
    matchnames <- paste(seq_nrand)
    extranames <- setdiff(user_names, matchnames)
    # if names are not in "", "2", "3"... try to match to other names 
    # it's easy to specify repeated names inadvertently (cf tests with lambda=c(fitmale$lambda,fitfem$lambda));
    # .cosy_names_lambdas(names(namesTerms)) disambiguates repeated names, 
    # 
    if (othernames_tried <- (length(extranames) &&  ! is.null(namesTerms))) { 
      ## a list, which names correspond to the grouping variable, and elements are the names of the coefficients fitted
      matchnames <- .cosy_names_lambdas(names(namesTerms)) # because these are possible names in some uses; namesTerms must be in the 'canonical' order
      extranames <- setdiff(user_names, matchnames)
    }
    if (length(extranames)) { ## no set of names match
      if (length(user_lFix)!=nrand) {
        stop("'fixed lambda' vector cannot be matched to random effects (different lengths & non-standard names).") ## needs to include NA's
      } else template[] <- user_lFix
    } else if (othernames_tried) {
      pos_in_lFix <- seq_nrand[match(user_names,matchnames)]
      if ( ! full_lambda) names(template) <- sort(pos_in_lFix) ## do not modify existing names
      template[as.character(pos_in_lFix)] <- user_lFix 
    } else if (length(user_names)) { # names were 2,3... not assuming completeness nor order
      template[user_names] <- user_lFix
    } else template[] <- user_lFix ## keep lFix names # user did not use any names
  }
  return(template)
}

.reformat_phi <- function(user_pFix, n_models, full_phi) { # similar to .reformat_lambda() but given it formats to list, it's OK only for mv phi
  seq_nmodels <- seq_len(n_models)
  if (full_phi) { # 
    template <- vector(list(NULL),n_models)
    names(template) <- seq_nmodels # (character)
  } else if (is.null(user_pFix)) { # NULL input -> NULL output
    return(NULL)
  } else {
    template <- user_pFix # reformats but not full-length; do not introduce (nor remove) NA/NaN's
    if (is.null(names(template))) names(template) <- seq_along(template)
  }
  if ( ! is.null(user_pFix)) {
    user_names <- names(user_pFix)
    if (length(unique(user_names))!=length(user_names)) stop("Repeated names in names of user-specified phi. Check your input.")
    # 'user'_names are not necessarily 1,2... (fitted values have names from rhs of ranef terms) 
    matchnames <- paste(seq_nmodels)
    extranames <- setdiff(user_names, matchnames)
    if (length(extranames)) { ## no set of names match
      if (length(user_pFix)!=n_models) {
        stop("'fixed phi' vector cannot be matched to submodels (different lengths & non-standard names).") ## needs to include NA's
      } else template[] <- user_pFix
    } else if (length(user_names)) { # names were 2,3... not assuming completeness nor order
      template[user_names] <- user_pFix
    } else template[] <- user_pFix ## keep lFix names # user did not use any names
  }
  return(template)
}


.reformat_corrPars <- function(parlist,corr_families=NULL) {
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

# .post_process_fixed <- function(fixed,corr_families, hyper_info) { 
#   if (is.null(fixed)) return(NULL)
#   fixed <- .reformat_corrPars(fixed,corr_families=corr_families) ## standardize use of corrPars in the parlist
#   # I write "F I X" as a TAG for this modif type attribute:
#   #attr(fixed,"type") <- relist(rep("fix",length(unlist(fixed))),fixed) ## on veut une list pour pouvoir supp des elements par <- NULL
#   # if ( ! is.null(fixed$hyper)) {
#   #   for (char_rd in as.character(unlist(hyper_info$ranges))) {
#   #     if (is.null(fixed$corrPars[[char_rd]])) fixed$corrPars[[char_rd]] <- list() ## seems no longer necessary for .expand_hyper()
#   #   }
#   # }
#   return(fixed)
# }

.calc_moreargs <- function(processed, # possibly a list of environments -> .calc_range_info -> scans then to compute a mean(nbUnique) 
                           corr_types, fixed, init.optim, control_dist, NUMAX=50, LDMAX=50, 
                           KAPPAMAX=100.000001, # so that users can set it to 100...
                           init.HLfit, corr_info, verbose, lower, upper) {
  moreargs <- list() # structure(vector("list", length=length(corr_types)), names=seq_along(corr_types))
  for (rd in seq_along(corr_types)) {
    corr_type <- corr_types[rd]
    if (! is.na(corr_type)) {
      char_rd <- as.character(rd)
      has_adj <- (corr_type %in% c("SAR_WWt","adjacency") 
                   #&&  is.null(fixed$corrPars[[char_rd]]$rho) 
                   && (! is.numeric(init.HLfit[[char_rd]]$rho))) ## init.HLfit[[]]$rho NULL or NA) 
      if (has_adj) {
        if (corr_type=="SAR_WWt") decomp <- attr(corr_info$adjMatrices[[rd]],"UDU.")
        if (corr_type=="adjacency") decomp <- attr(corr_info$adjMatrices[[rd]],"symSVD")
        # decomp may be missing: cf comments in adjacency()$calc_moreargs
      }
      moreargs[[char_rd]] <- corr_info$corr_families[[rd]]$calc_moreargs(
        fixed=fixed, char_rd=char_rd, init.optim=init.optim, 
        processed=processed, rd=rd, control_dist=control_dist, 
        NUMAX=NUMAX, LDMAX=LDMAX, KAPPAMAX=KAPPAMAX,
        # quite distinct arguments for adjacency:
        decomp=decomp, verbose=verbose, lower=lower, upper=upper,
        # IMRF...:
        IMRF_pars=attr(attr(attr(processed$ZAlist,"exp_spatial_terms")[[rd]],"type"),"pars")
      )
      if (has_adj && ! is.null(rho <- fixed$corrPars[[char_rd]]$rho) ) {
        rhorange <- moreargs[[char_rd]]$rhorange
        if ((dif <- rho-rhorange[1L])<=0) stop(paste0("User-given rho too low by ",signif(dif,6)))
        if ((dif <- rho-rhorange[2L])>=0) stop(paste0("User-given rho too high by ",signif(dif,6)))
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

#.match_coords_in_uniqueGeo <- function(coords,uniqueGeo) {any(apply(uniqueGeo,1L,identical,y=coords))}

.checkDistMatrix <- function(distMatrix,data,coordinates) {
  if (inherits(distMatrix,"dist")) {
    userdistnames <- labels(distMatrix)
  } else if (inherits(distMatrix,"matrix")) {
    userdistnames <- rownames(distMatrix)
  } else message(paste("(!) 'distMatrix' is neither a 'matrix' or 'dist' object. Check the input. I exit."))
  ## chol() fails on distances matrices with repeated locations
  ## the following code assumes that distMatrix deals only with unique locations, and checks this
  ## HENCE ******* distMatrix must refer to unique values of a grouping variable *********
  #
  # .checkDistMatrix() is called after a number of operations on the input distMatrix,
  # including a possible reduction of the .distMatrix for a Matern(LHS|.) term. 
  # In that case, left_of_bar variables in Matern(X|.) are not available here. 
  # Therefore, I cannot distinguish here whether :
  # a distance is missing from the distmat after reduction because the position is not retained according to X (which may be correct), 
  # or if the position should really be there.
  # This leads to a perhaps weaker check than ideal : (_F I X M E_ ?)
  checknames <- setdiff(userdistnames, rownames(data)) 
  if (length(checknames)) {
    warning("The rownames of 'distMatrix' are not contained in those of the 'data'. Further checking of 'distMatrix' is not possible.", immediate=TRUE)
    # nbUnique <- NA # => init rho is NA =>bug (and one does not even see the warning)    
    nbUnique<- nrow(unique(distMatrix))
  } else { # further checking on uniqueness of coordinates (ie on their 'levels') and their order
    dist_ordered_uniqueGeo <- .calcUniqueGeo(data=data[userdistnames, coordinates, drop=FALSE]) ## check that this corresponds to unique locations
    nbUnique <- nrow(dist_ordered_uniqueGeo)
    if (nbUnique != nrow(distMatrix)) {
      stop(paste0("The dimension of 'distMatrix' (",nrow(distMatrix),
                  " rows) does not match the number of levels (",nbUnique,") of the grouping variable"))
    } else { ## check order
      dataordered_namesindist <- intersect(rownames(data),userdistnames)
      if (! all(userdistnames==dataordered_namesindist)) {
        stop("The order of appearance of locations in 'distMatrix' is not the same as that in the 'data'.")
      }
    } 
  }
  nbUnique ## if stop() did not occur
}


.get_dist_nested_or_not <- function(spatial_term, data, distMatrix, uniqueGeo, dist.method, as_matrix=FALSE,
                                    needed=c(), geo_envir) {
  if (is.na(needed["distMatrix"])) needed["distMatrix"] <- FALSE 
  if (is.na(needed["uniqueGeo"])) needed["uniqueGeo"] <- FALSE 
  if (is.na(needed["nbUnique"])) needed["nbUnique"] <- FALSE 
  if (is.na(needed["notSameGrp"])) needed["notSameGrp"] <- FALSE 
  coordinates <- .get_coordinates(spatial_term=spatial_term, data=data)
  geoMats <- list(coordinates=coordinates)
  if (is.null(distMatrix) || needed["notSameGrp"]) { # In principle we avoid computing anything new if 
    #                     distMatrix (implying uniqueGeo too )has been computed
    #                     but for grouped effects, defining the optim range may have used distMatrix in a case where we now need notSameGrp.
    coord_within <- .extract_check_coords_within(spatial_term=spatial_term) 
    coords_nesting <- setdiff(coordinates,coord_within)
    coord_within <- setdiff(coordinates, coords_nesting)
    ## needed in most cases for evaluation of further elements:
    if (is.null(uniqueGeo) ) { ## then construct it from the data ## this should be the routine case, except for AR1
      uniqueGeo <- .calcUniqueGeo(data=data[,coord_within,drop=FALSE])
      # activelevels are data-ordered levels whose ranef values affect the likelihood
      # .calcUniqueGeo must produce data-ordered values.
      if (!is.null(geo_envir$activelevels)) uniqueGeo <- uniqueGeo[geo_envir$activelevels,,drop=FALSE]
      # Setting 'raw_levels' as row names to allow levels check in .calc_ZAlist(): 
      x <- t(uniqueGeo)
      pastestring <- paste("list(",paste("x","[",seq_len(nrow(x)),",]",sep="",collapse=","),")",sep="")
      raw_levels <- do.call(paste,c(eval(parse(text = pastestring)),sep=":"))
      attr(uniqueGeo,"names_ori") <- rownames(uniqueGeo) 
      rownames(uniqueGeo) <- raw_levels 
    } 
    geoMats$nbUnique <- nrow(uniqueGeo) ## only serves to control RHOMAX and should not be computed from the final uniqueGeo in case of nesting
    if (length(coords_nesting) && any(needed[c("distMatrix","uniqueGeo","notSameGrp")]) ) {
      ## should now (>2.3.9) work on data.frames, nesting factor not numeric.
      e_uniqueGeo <- .calcUniqueGeo(data=data[,coordinates,drop=FALSE])
      ## The rows_bynesting <- by(e_uniqueGeo ,e_uniqueGeo[,coords_nesting],rownames) line in.expand_GeoMatrices()
      #  works only if the cols are not factors ! (an as.matrix() ?). Unless we have a fix for this,
      #  e_uniqueGeo classes should not be factor; and as.integer(<factor>) can produce NA's hence as.character()
      # same operation was performed to generate uniqueGeo in .calc_AR1_sparse_Q_ranges()
      for (nam in coords_nesting) {if (is.factor(fac <- e_uniqueGeo[[nam]])) e_uniqueGeo[[nam]] <- as.character(levels(fac))[fac]}
      if (needed["notSameGrp"]) geoMats$notSameGrp <- .expand_grouping(w_uniqueGeo=uniqueGeo, e_uniqueGeo=e_uniqueGeo, 
                                                            coords_nesting=coords_nesting, coord_within=coord_within)
      if (any(needed[c("distMatrix","uniqueGeo")])) {
        eGM <- .expand_GeoMatrices(w_uniqueGeo=uniqueGeo, e_uniqueGeo=e_uniqueGeo, 
                                   coords_nesting=coords_nesting, coord_within=coord_within,dist.method=dist.method)
        distMatrix <- eGM$distMatrix
      }
      uniqueGeo <- e_uniqueGeo
    } else if (needed["distMatrix"]) {
      notnumeric <- ! unlist(lapply(uniqueGeo,is.numeric))
      if (any(notnumeric)) stop(paste0("Variables(s) '",paste(names(which(notnumeric)),collapse="' '"),"'",
                                      " are not numeric, hence not suitable for computation of distance matrix."))
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


.provide_AR_factorization_info <- local({
  RSpectra_warned <- FALSE
  function(adjMatrix, sparse_precision, corr.model) {
    if (corr.model  %in% c("SAR_WWt")) {
      decomp <- eigen(adjMatrix,symmetric=FALSE) ## could be symmetric=TRUE if adjMatrix is dsC as in adjacency case.
      return(list(u=decomp$vectors,d=decomp$values,u.=solve(decomp$vectors)))
    }
    # ELSE
    
    if (corr.model  %in% c("adjacency")) { # adjMatrix is dsC from .preprocess -> ... -> .sym_checked(adjMatrix)
      ## extreme eigenvalues needed in all cases for the bounds. Full decomp not always needed
      if ( sparse_precision) { 
        if (requireNamespace("RSpectra",quietly=TRUE)) { #https://scicomp.stackexchange.com/questions/26786/eigen-max-and-minimum-eigenvalues-of-a-sparse-matrix
          eigrange <- RSpectra::eigs(adjMatrix, k=2, which="BE", opts=list(retvec=FALSE))$values
          decomp <- list(eigrange=eigrange) # only the extreme eigenvalues
        } else {
          if ( ! RSpectra_warned) { #if ( ! identical(spaMM.getOption("RSpectra_warned"),TRUE)) {
            message("If the 'RSpectra' package were installed, an eigenvalue computation could be faster.")
            RSpectra_warned <<- TRUE # .spaMM.data$options$RSpectra_warned <- TRUE
          }
          eigvals <- eigen(adjMatrix, symmetric=TRUE, only.values = TRUE)$values # first converts to dense matrix, so quite inefficient.
          decomp <- list(eigrange=range(eigvals)) 
        }
      } else {
        decomp <- eigen(adjMatrix, symmetric=TRUE)
        svdnames <- names(decomp)
        svdnames[svdnames=="values"] <- "d"
        svdnames[svdnames=="vectors"] <- "u"
        names(decomp) <- svdnames
        decomp$adjd <- decomp$d
        decomp$eigrange=range(decomp$adjd)
      }
      return(decomp)
    }
  }
})

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
      # handling infinite values used in nested geostatistical models
      maxrange <- max(geo_envir$distMatrix[! is.infinite(geo_envir$distMatrix)])-min(geo_envir$distMatrix)
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

.TRACE_fn <- function(ranFix, processed) {
  ranPars <- .canonizeRanPars(ranFix,corr_info=NULL,checkComplete = FALSE, rC_transf=.spaMM.data$options$rC_transf)
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
    if (is.null(processed$vec_nobs)) { # Is this to catch the case where a full phi vector is provided ? 
      cat("(phi fixed) ")
      ranPars[["phi"]] <- NULL
    } else { # mv case
      if (is.list(phi <- ranPars[["phi"]])) { # HYPOTHETICAL mv case
        for (mv_it_st in names(phi)) {
          if (length(phi[[mv_it_st]])>1) {
            cat(paste0("(phi_",mv_it_st, "fixed) "))
            ranPars[["phi"]][mv_it_st] <- list(NULL)
          }
        }
      } else { # but the simple mv case is that of a phi vector
        
      }
    }
  }
  ntC <- names(ranPars)
  for (lit in seq_along(ranPars)) {
    urP <- unlist(ranPars[[lit]]) ## ranPars$corrPars can be list() in which case urP is NULL 
    if (!is.null(urP)) cat(ntC[lit],"=",paste(signif(urP,6),collapse=" ")," ")
  }
}


## creates precisionFactorList if it doesn't exist, and partially fills it with Diagonal()'s + corrMatrix info
.init_AUGI0_ZX_envir_spprec_info <- function(processed) {
  AUGI0_ZX_envir <- processed$AUGI0_ZX$envir
  if (is.null(precisionFactorList <- AUGI0_ZX_envir$precisionFactorList)) {
    nranef <- length(AUGI0_ZX_envir$finertypes)
    updateable <- rep(FALSE,nranef)
    precisionFactorList <- latent_d_list <- vector("list",nranef) ## will contain diagonal matrices/info for non-trivial (diagonal) precision matrices
    cum_n_u_h <- processed$cum_n_u_h
    for (rd in seq_len(nranef)) {
      finertype <- AUGI0_ZX_envir$finertypes[rd]
      if ( finertype %in% c("adjacency","AR1") ) {
        ## terms of these types must be dealt with by ad hoc code for each type elsewhere
      } else if ( finertype=="corrMatrix") {
        cov_info_mat <- processed$corr_info$cov_info_mats[[rd]]
        if (inherits(cov_info_mat,"precision")) {
          sparse_Qmat <- drop0(cov_info_mat[["matrix"]])
        } else stop(' ! inherits(cov_info_mat,"precision")')
        ####
        # Compared to more general case, we don't need to store a template! This is always ONE-TIME CODE.
        #  I once had a problem with perm_Q=TRUE in test-predVar-Matern-corrMatrix -> predict(f2,.) so make sure to check this when changing the code
        if (is.null(perm_Q <- .spaMM.data$options$perm_Q)) { # assuming dsC...
          perm_Q <- (length(sparse_Qmat@x)/ncol(sparse_Qmat)^2)<0.8 # quick exclusion of precision matrices that are fully dense so that no permutation would be useful. 
        }
        Q_CHMfactor <- Cholesky(sparse_Qmat,LDL=FALSE,perm=perm_Q) 
        if (perm_Q) {
          .ZA_update(rd, Q_CHMfactor, processed, 
                     Amat=attr(processed$ZAlist,"AMatrices")[[as.character(rd)]])
          # One-time ZA updating, but also here specifically for corrMatrix, one-time construction of sparse_Qmat
          permuted_Q <- attr(attr(processed$ZAlist,"AMatrices")[[as.character(rd)]],"permuted_Q") 
          # $Qmat <- sparse_Qmat will be used together with ZA independently from the CHM to construct the Gmat
          # If we use permuted Chol, then we must permute sparse_Qmat, by 
          if (identical(permuted_Q,TRUE)) sparse_Qmat <- tcrossprod(as(Q_CHMfactor,"sparseMatrix")) 
          # precisionFactorList[[rd]]$template <- Q_CHMfactor # No updating ever needed
        }
        precisionFactorList[[rd]] <- list(Qmat=sparse_Qmat, # should by *dsC* 
                                          chol_Q=as(Q_CHMfactor, "sparseMatrix")) # Linv
        ####
        .assign_Lunique_from_Q_CHM(processed=processed, rd=rd, Q_CHMfactor=Q_CHMfactor, corr_type="corrMatrix", 
                                   spatial_term=attr(processed$ZAlist,"exp_spatial_terms")[[rd]], type="from_Q_CHMfactor") ## assigns to AUGI0_ZX$envir
      } else if ( finertype=="IMRF" ) {
        ## terms of these types must be dealt with by ad hoc code for each type elsewhere, but we can set
        updateable[rd] <- TRUE
      } else if ( finertype =="(.|.)" ) { # EXcludes "ranCoefs"
        nc <- ncol(processed$ZAlist[[rd]])
        precisionFactorList[[rd]] <- list(Qmat=.symDiagonal(n=nc), ## used independently of chol_Q_list, see precisionBlocks
                                          chol_Q=new("dtCMatrix",i= 0:(nc-1L), p=0:(nc), Dim=c(nc,nc),x=rep(1,nc)) )
        # All chol_Q's must be dtCMatrix so that bdiag() gives a dtCMatrix
        updateable[rd] <- TRUE
      } else if ( finertype %in% c("Matern", "Cauchy") ) { 
        ## leave precisionFactorList[[rd]] NULL
        updateable[rd] <- TRUE
      } else if ( finertype != "ranCoefs") stop(paste("sparse-precision methods were requested, but",
                                                      finertype,"terms are not yet handled by sparse precision code."))
    }
    diff_n_u_h <- diff(cum_n_u_h)
    for (rd in seq_along(diff_n_u_h)) latent_d_list[[rd]] <- rep(1,diff_n_u_h[rd])
    AUGI0_ZX_envir$precisionFactorList <- precisionFactorList
    AUGI0_ZX_envir$latent_d_list <- latent_d_list
    AUGI0_ZX_envir$updateable <- updateable
  }
  ## environment modified, no return value
}

## updates precisionFactorList with LMatrices info
.update_AUGI0_ZX_envir_ranCoefs_info <- function(processed, LMatrices=NULL) {
  if ( ! is.null(LMatrices)) {
    AUGI0_ZX_envir <- processed$AUGI0_ZX$envir
    precisionFactorList <- AUGI0_ZX_envir$precisionFactorList
    latent_d_list <- AUGI0_ZX_envir$latent_d_list
    updateable <- AUGI0_ZX_envir$updateable
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
      d_rt <- latentL_blob[["d"]]
      nc <- length(d_rt) 
      latent_d_list[[rt]] <- d_rt[gl(nc, length(latent_d_list[[rt]])/nc)]
      updateable[rt] <- (length(which(compactchol_Q@x!=0)) == nc*(nc+1L)/2L)
    }
    for (rt in which(is_given_by=="inner_ranCoefs")) {
      # : need to initialize because spprec IRLS requires .update_AUGI0_ZX_envir_ranCoefs_info() to have filled all matrices
      nc <- ncol(processed$ranCoefs_blob$longLv_templates[[rt]] )
      precisionFactorList[[rt]] <- list(precmat=.symDiagonal(n=nc), 
                                        chol_Q=new("dtCMatrix",i= 0:(nc-1L), p=0:(nc), Dim=c(nc,nc),x=rep(1,nc)) )
      updateable[rt] <- TRUE
    }
    AUGI0_ZX_envir$precisionFactorList <- precisionFactorList
    AUGI0_ZX_envir$latent_d_list <- latent_d_list
    AUGI0_ZX_envir$updateable <- updateable
  }
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
  # envir$updateable[rt] should be TRUE and remain so for those terms
  envir$precisionFactorList <- precisionFactorList
  ## environment modified, no return value
}

.init_optim_phi <- function(phimodel, processed, init.optim, nrand, reasons_for_outer) {
  if (phimodel == "phiScal") {
    if (reasons_for_outer) {
      ## default method is outer but can be reversed if init is NaN (and to test NaN we need to test NULL)
      not_inner_phi <- is.null(init.optim$phi) || ! is.nan(init.optim$phi) ## outer if NULL, NA or numeric
    } else {
      ## default is inner but can be reversed if numeric or NA (hence neither NULL nor NaN)
      not_inner_phi <- ! (is.null(init.optim$phi) || is.nan(init.optim$phi))
    }
  } else not_inner_phi <- FALSE ## complex phi model, we weed inner optim
  if (not_inner_phi) {
    if (is.null(init.optim$phi)) { 
      init.optim$phi <- .get_inits_by_glm(processed)$phi_est/(nrand+1L) ## at least one initial value should represent high guessed variance
      # if init.optim$phi too low (as in min(.,2)) then fitme(Reaction ~ Days + AR1(1|Days) + (Days|Subject), data = sleepstudy) is poor
    }  
  } else {
    init.optim$phi <- init.optim$phi[ ! is.nan(init.optim$phi)]
    if ( ! length(init.optim$phi)) init.optim["phi"] <- NULL
  }
  return(list(not_inner_phi=not_inner_phi, init.optim=init.optim))
}

.eval_init_lambda_guess <- function(processed, stillNAs, ZAL=NULL, cum_n_u_h,For) {
  nrand <-  length(processed$ZAlist)
  if (is.null(processed$main_terms_info$Y)) { ## for resid model
    if (For=="optim") { 
      guess_from_glm_lambda <- 0.1 ## (FIXME ad hoc) Cannot let it NA as this will be ignored by further code, namely
      #                               optim_lambda_with_NAs[which_NA_simplelambda] <- init_lambda[which_NA_simplelambda]
      #                               init.optim$lambda <- optim_lambda_with_NAs[!is.na(optim_lambda_with_NAs)]
    } else guess_from_glm_lambda <- NA
  } else {
    inits_by_glm <- .get_inits_by_glm(processed) 
    if (For=="optim") { 
      guess_from_glm_lambda <- inits_by_glm$lambda*(3L*nrand)/((nrand+1L)) # +1 for residual
    } else guess_from_glm_lambda <- inits_by_glm$lambda*(3L*nrand+2L)/((nrand+1L)) # +1 for residual  ## old code was *5/(nr+1) ## f i x m e super ad hoc
  }
  init_lambda <- rep(NA,nrand)
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  corr_types <- processed$corr_info$corr_types
  for (rd in stillNAs) { ## fam_corrected_guess for each ranef in stillNAs
    if ( ! is.null(processed$families)) {
      which_mv <- attr(processed$ZAlist[[rd]],"which_mv")
      link_ <- unlist(lapply(processed$families[which_mv],`[[`,i="link")) 
    } else link_ <- processed$family$link
    if (is.null(ZAL)) {
      ZA <- processed$ZAlist[[rd]]
    } else if (inherits(ZAL,"ZAXlist")) {
      ZA <- ZAL@LIST[[rd]]
      if (inherits(ZA,"ZA_QCHM")) ZA <- ZA$ZA
    } else {
      u.range <- (cum_n_u_h[rd]+1L):cum_n_u_h[rd+1L]
      ZA <- ZAL[,u.range,drop=FALSE]
    }
    denom <- colSums(ZA*ZA) # colSums(ZA*ZA) are diag crossprod(ZA) (memo: length=ncol(ZA)). 
    denom <- denom[denom!=0] ## so that same result in dense and sparse if ZA has empty cols in sparse
    ZA_corrected_guess <- guess_from_glm_lambda/sqrt(mean(denom)) 
    #if (corr_types[it]=="AR1") ZA_corrected_guess <- log(1.00001+ZA_corrected_guess) ## ad hoc fix but a transformation for ARphi could be better FIXME
    fam_corrected_guess <- .calc_fam_corrected_guess(guess=ZA_corrected_guess, link_=link_, For=For, processed=processed, nrand=nrand)
    init_lambda[rd] <- .preprocess_valuesforNAs(rd, lcrandfamfam=lcrandfamfam, 
                                                rand.families=rand.families, init.lambda=fam_corrected_guess)
  }
  if (For != "optim") {   ## If called by HLfit: the present pmax() matters.
    init_lambda[stillNAs] <- pmax(init_lambda[stillNAs],1e-4) 
  } ## ELSE If called by fitme: calc_init_dispPars runs pmax on the whole vector
  return(init_lambda)
}



.init_optim_lambda_ranCoefs <- function(proc1, not_inner_phi, init.optim, nrand, ranCoefs_blob, var_ranCoefs) {
  lFix <- proc1$lambda.Fix ## only for 'simple" ranefs with Xi_cols=1
  if (proc1$augZXy_cond || anyNA(init.optim$lambda) ||
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
        #rc <- rep(0.0001,Xi_ncol*(Xi_ncol+1L)/2L)
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
        #rc[lampos] <- 1e-3*fam_corrected_guess/(Xi_ncol)
        #rc[lampos[1]] <- fam_corrected_guess/(Xi_ncol)
        init.optim$ranCoefs[[char_rt]] <- rc ## see help(ranCoefs)
      }
    }
  }
  return(init.optim)
}


.more_init_optim <- function(proc1, processed, corr_types, init.optim, phi_by_augZXy) {
  ## trying to guess all cases where optimization is useful. But FIXME: create all init and decide afterwardsS
  phimodel1 <- proc1$models[['phi']]
  allPhiScalorFix <- all(processed$models[['phi']] == "phiScal" | processed$models[['phi']] == "")
  ranCoefs_blob <- processed$ranCoefs_blob
  is_MixedM <- ( ! is.null(ranCoefs_blob) )
  if (is_MixedM) {
    var_ranCoefs <- (ranCoefs_blob$isRandomSlope & ! ranCoefs_blob$is_set) # vector !
    has_corr_pars <- length(corr_types[ ! is.na(corr_types)])
  } else var_ranCoefs <- has_corr_pars <- FALSE
  sufficient_reasons_for_outer <- (anyNA(init.optim$lambda) || # first one meaning that the user explictly set a NA init lambda
                                 any(var_ranCoefs) || 
                                 has_corr_pars ||
                                 ( has_family_par <- (( ! is.null(init.optim$COMP_nu)) || ( ! is.null(init.optim$NB_shape))) ) || # lambda + family pars # motivated by 'ahzut' example in private test-COMPoisson-difficult.R
                                 ( has_estim_families_par <- identical(attr(processed$families,"has_estim_families_par"), TRUE))) ## identical bc absent in multi() case
  if (
    sufficient_reasons_for_outer ||  
    var(proc1$y)<1e-3 || # mixed model with low response variance (including case where phi is fixed (phimodel="") )
    (phimodel1 == "phiHGLM" && identical(spaMM.getOption("outer_optim_resid"),TRUE)) ||
    allPhiScalorFix || # lambda + phi
    { # reach here if several (lambda || fixed ranCoefs) and no var_ranCoefs; 
      lFix <- processed$lambda.Fix
      if (! is.null(ranCoefs_blob)) lFix <- lFix[ ! ranCoefs_blob$isRandomSlope] # Fixed ranCoefs count as a NA in $lambda.Fix so we remove them to count variable lambda's
      length(which(is.na(lFix)))>1L # so that two lambdas is handled as 1 lambda+phi
    } 
  ) { ## All cases where outer optimization MAY be better for dispersion parameters (or forced outer optim)
    nrand1 <- length(proc1$ZAlist)
    reasons_for_outer <- sufficient_reasons_for_outer
    if ( ! reasons_for_outer) {
      reasons_for_outer <- outer_spares_costly_dleve_comput <- ( 
        # the dleve computations for HL[2L] with costly .calc_dvdlog...Mat_new() occur in inner for HL[2L] irrespective of vecdisneeded[3]
        # and in particular even for canonical link GLMM. The dleve terms should be nonzero when GLM weights !=1 i.e vecdisneeded[2] (or [1])
        (nrand1 > 0L && processed$cum_n_u_h[length(processed$cum_n_u_h)] > 1000L) &&
          ( calc_dvdlogdisp_needed_for_inner_ML <-  (processed$vecdisneeded[2] && processed$HL[2L]) ) 
      )
    }
    if ( ! reasons_for_outer) {
      ## we look whether the hatval_Z are needed or not in .calc_sscaled_new(), hence in fitme() with outer estim of dispersion params
      hatval_Z_is_costly <- ( NROW(processed$y)>1000L || 
                                (nrand1>0L && processed$cum_n_u_h[length(processed$cum_n_u_h)]>200L))
      reasons_for_outer <- outer_spares_costly_hatval_Z <- ( ## "but inner requires it"
        hatval_Z_is_costly && allPhiScalorFix && (
          hatval_Z_needed_for_inner_ML <-  (! NCOL(processed$X.Re)) ## X.Re is a 0-col matrix
        ) && ( 
          hatval_Z_not_needed_for_sscaled <- ! (processed$HL[1L] && any(processed$vecdisneeded)) 
        )  
      )
    }
    # so that reasons_for_outer = ( sufficient_reasons_for_outer || outer_spares_costly_dleve_comput || outer_spares_costly_hatval_Z)
    if (is.null(proc1$phi.Fix) && ! phi_by_augZXy ) { # Set (or not) outer optimization for phi: 
      init_optim_phi_blob <- .init_optim_phi(phimodel1, proc1, init.optim, nrand1, reasons_for_outer)
      not_inner_phi <- init_optim_phi_blob$not_inner_phi # which may indicate outer phi or no phi estimation
      init.optim <- init_optim_phi_blob$init.optim
    } else {
      not_inner_phi <- reasons_for_outer
    }
    if (length(proc1$ZAlist)>0L) {
      # Set outer optimization for lambda and ranCoefs (handling incomplete ranFix$lambda vectors) through call to .init_optim_lambda_ranCoefs()
      init.optim <- .init_optim_lambda_ranCoefs(proc1, not_inner_phi, init.optim, nrand1, proc1$ranCoefs_blob, var_ranCoefs)
    } 
    if (length(init.optim$lambda) > nrand1) {
      mess <- paste0("Initial value for lambda: (",
                     paste(#names(init.optim$lambda),"=",
                       init.optim$lambda, collapse=","),
                     ") has more elements than there are random effects.\n   The fit will likely fail.")
      warning(mess, immediate.=TRUE)
    }
  } ## else no addition to the inits, everything will be inner optimized
  return(init.optim)
}

# called by fitme_body/corrHLfit_body/(fitmv_body on individual models) hence no multivariate input except possibly through 'processed'
.calc_optim_args <- function(proc_it, 
                             processed, # multi() and mv
                             init, fixed, lower, upper, verbose, optim.scale, For) {
  corr_info <- proc_it$corr_info 
  # modify HLCor.args and <>bounds;   ## distMatrix or uniqueGeo potentially added to HLCor.args:
  corr_types <- corr_info$corr_types
  if ( ! is.null(fixed)) fixed <- .reformat_corrPars(fixed, corr_families=corr_info$corr_families)
  if ( ! is.null(init$phi) ) {
    if (proc_it$models[["phi"]]=="") {
      warning("initial value for 'phi' is ignored when there is no phi parameter (e.g. poisson or binomial families)") # i.e. anything but Intercept model
      init$phi <- NULL # erase the un-usable value otherwise the residModel would be ignored!
    } else if (proc_it$models[["phi"]]!="phiScal") {
      warning("initial value for 'phi' is ignored when there is a non-default resid.model") # i.e. anything but Intercept model
      init$phi <- NULL # erase the un-usable value otherwise the residModel would be ignored!
    }
  }
  init.optim <- .reformat_corrPars(init, corr_families=corr_info$corr_families)
  init.HLfit <- proc_it$init_HLfit #to be modified below ## dhglm uses fitme_body not (fitme-> .preprocess) => dhglm code modifies processed$init_HLfit
  init <- NaN ## make clear it's not to be used
  # used by .more_init_optim:
  family <- proc_it$family
  if (family$family=="COMPoisson") {
    if (inherits(substitute(nu, env=environment(family$aic)),"call")) {
      if (is.null(init.optim$COMP_nu)) init.optim$COMP_nu <- 1 # template: .calc_inits will modify it according to lower, upper 
    } else {
      if ( ! is.null(init.optim$COMP_nu)) {
        warning("initial value is ignored when 'COMP_nu' is fixed.") # i.e. anything but Intercept model
        init.optim$COMP_nu <- NULL
      } # and this should have the effect that user lower and upper values should be ignored too.
    }  
  } else if (family$family == "negbin") {
    if (inherits(substitute(shape, env=environment(family$aic)),"call")) {
      if (is.null(init.optim$NB_shape)) init.optim$NB_shape <- 1 # idem
    } else {
      if ( ! is.null(init.optim$NB_shape)) {
        warning("initial value is ignored when 'NB_shape' is fixed.") # i.e. anything but Intercept model
        init.optim$NB_shape <- NULL
      } # and this should have the effect that user lower and upper values should be ignored too.
    }  
  }
  #
  ##### init.optim$phi/lambda will affect calc_inits -> calc_inits_dispPars.
  # outer estim seems useful when we can suppress all inner estim (thus the hatval calculations). 
  # ./. Therefore, we need to identify all cases where phi is fixed, 
  # ./. or can be outer optimized jointly with lambda:  
  # Provide ranCoefs inits by calling .init_optim_lambda_ranCoefs:
  user.lower <- .reformat_corrPars(lower,corr_families=corr_info$corr_families)
  user.upper <- .reformat_corrPars(upper,corr_families=corr_info$corr_families) ## keep user input 
  if (For %in% c("fitme","fitmv") && proc_it$HL[1]!="SEM") { ## for SEM, it's better to let SEMbetalambda find reasonable estimates
    augZXy_cond <- proc_it$augZXy_cond # processed$augZXy_cond would fail in multi() case
    # We may use augZXy to estimate lambda and phi (more exactly, a scaling factor common to them), or lambda only
    # In .preprocess_augZXy(), we decided not to use it to estimate lambda when there is an init phi. 
    # But we could have decided otherwise (still use augZXy to estimate lambda in that case)
    # Now we check if there is an upper|lower phi, in which case we cannot use it to estimate phi
    # But we face the same question of using augZXy to estimate lambda in that case
    if (FALSE) { # This is hypothetical code for the case where we wish to maximize use of augZXy:
      # In that case .preprocess_augZXy() should be modified to allow augZXy to estimate lambda in the presence of an init phi;
      # phi_by_augZXy should depend on this init; and phi_by_augZXy should not affect proc1$augZXy_cond
      phi_by_augZXy <- ( augZXy_cond && is.null(init.optim$phi) && is.null(user.lower$phi) && is.null(user.upper$phi)) 
      # NOTHING TO DO to processed$augZXy_cond
      # but 'phi_by_augZXy' still separately needed in call below.
    } else { # This is for the case where we decide not to use augZXy_cond for lambda whene we cannot use it for phi:
      # (yet we seem to do that in .makeCovEst1()... as controlled by the "inner" attribute)
      # In that case .preprocess_augZXy() must have set augZXy to FALSE in the presence of an init phi; and:
      phi_by_augZXy <- ( augZXy_cond && is.null(user.lower$phi) && is.null(user.upper$phi)) 
      if (augZXy_cond && ! phi_by_augZXy) {
        proc_it$augZXy_cond <- structure(phi_by_augZXy, inner=attr(proc_it$augZXy_cond, "inner"))
        .do_TRACE(processed)
      }
    }
    # Now create a template for optimization, deciding for outer/inner optimisations:
    if (For=="fitmv") {
      init.optim <- .more_init_optim(proc1=proc_it, 
                                     processed=processed, # critical mv info 
                                     corr_types=corr_types, init.optim=init.optim, 
                                     phi_by_augZXy = phi_by_augZXy)  
    } else init.optim <- .more_init_optim(proc1=proc_it, 
                                          processed=proc_it, # even for multi()
                                          corr_types=corr_types, init.optim=init.optim, 
                                          phi_by_augZXy = phi_by_augZXy)  
  }
  .check_conflict_init_fixed(fixed,init.optim, "given as element of both 'fixed' and 'init'. Check call.")
  .check_conflict_init_fixed(init.HLfit,init.optim, "given as element of both 'init.HLfit' and 'init'. Check call.") ## has quite poor effect on fits
  if (For=="fitmv") { # in that case the ranefs differ across submodels and we want submodel-specific proc1 to be used eg in 
    #   IMRF_pars=attr(attr(attr(processed$ZAlist,"exp_spatial_terms")[[rd]],"type"),"pars")
    moreargs <- .calc_moreargs(processed= proc_it, 
                               # We might need something equiv to .calc_range_info computing a mean(nbUnique) from multi()-processed
                               # But (1) .makeLowUp_stuff_mv() is *the* place for handling processed itself, 
                               #           so does it or can it deal with that ? __FIXME__
                               #     (2) structure of mv-processed differs from that of multi()-processed.
                               #     (3) Even the structure of ranefs of mv-processed$unmerged differs from that of multi()-processed.
                               corr_types=corr_types, fixed=fixed, init.optim=init.optim, control_dist=proc_it$control_dist, 
                               init.HLfit=init.HLfit, corr_info=corr_info, verbose=verbose, lower=lower, upper=upper)  
  } else moreargs <- .calc_moreargs(processed=processed, # possibly a list of environments -> .calc_range_info -> scans then to compute a mean(nbUnique) 
                             corr_types=corr_types, fixed=fixed, init.optim=init.optim, control_dist=proc_it$control_dist, 
                             init.HLfit=init.HLfit, corr_info=corr_info, verbose=verbose, lower=lower, upper=upper)
  fixed <- .expand_hyper(fixed, hyper_info=proc_it$hyper_info, moreargs=moreargs)
  inits <- .calc_inits(init.optim=init.optim, init.HLfit=init.HLfit,
                       ranFix=fixed,  corr_info=corr_info,
                       moreargs=moreargs,
                       user.lower=user.lower, user.upper=user.upper,
                       optim.scale=optim.scale, 
                       For=For, hyper_info=proc_it$hyper_info
  )
  .check_conflict_init_fixed(fixed,inits$init.optim, "present in both 'fixed' and 'inits$init.optim'. Contact the maintainer.")
  .check_conflict_init_fixed(init.HLfit,inits$init.optim, "present in both 'init.HLfit' and 'inits$init.optim'. Contact the maintainer.") 
  ################
  ################
  if (For == "fitmv") {
    # The local moreargs may be defective for adjacency models 
    # (as we try to avoid determining spprec and what depends on it when preprocessing eahc submodel)
    # And then, any merged LowUp will be defective.
    # Instead, .makeLowUp_stuff_mv() will compute mv-LUarglist and mv-LowUp from the global "processed"  and the merged other arguments
    return(list(inits=inits, fixed=fixed, corr_types=corr_types))
  } else {
    LUarglist <- list(canon.init=inits$`init`, 
                      init.optim=inits$init.optim,
                      user.lower=user.lower,user.upper=user.upper,
                      corr_types=corr_types,
                      ranFix=fixed, # inits$ranFix, # Any change in $ranFix would be ignored 
                      optim.scale=optim.scale, 
                      moreargs=moreargs) ## list needed as part of attr(,"optimInfo")
    LowUp <- do.call(".makeLowerUpper",LUarglist) 
    return(list(inits=inits, fixed=fixed, corr_types=corr_types, LUarglist=LUarglist,LowUp=LowUp))
   ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  } 
} 

.is.multi <- function(family) {
  # Where call at the end of fitting fn, 'family' is the evaluated argument of the call, possibly not yet interpreted as final family()
  # i.e,, "negbin" -> character, negbin-> function, negbin() -> family (which is also a list!)
  # EXCEPT that multi(...) -> list but not family
  return(
    is.list(family) && # i.e. value of multi(...), or value of <family>() 
      family$family=="multi"
  )
}

.anyNULL <- function(x) {
  if (is.null(x)) {
    return(TRUE)
  } else if (is.list(x)) {
    if (is.null(anyNULL <- attr(x,"anyNULL"))) anyNULL <- any(lengths(x)==0L)
    return(anyNULL)
  } else return(FALSE) # stop("Unhandled 'x' class in .anyNULL().") # if numeric, not any null...
}

.allNULL <- function(x) {
  if (is.null(x)) {
    return(TRUE)
  } else if (is.list(x)) {
    if (is.null(allNULL <- attr(x,"anyNULL"))) allNULL <- all(lengths(x)==0L)
    return(allNULL)
  } else  return(FALSE) # stop("Unhandled 'x' class in .allNULL().") # if numeric, not all null...
}

.check_suspect_rho <- function(corr_types, ranPars_in_refit, LowUp) {
  locmess <- NULL
  if (any(geostat <- corr_types %in% c("Matern","Cauchy"))) {
    char_rnf <- as.character(which(geostat))
    corlow <- unlist(LowUp$lower$corrPars[char_rnf])
    if (length(corlow)) {
      corpars <- unlist(ranPars_in_refit$corrPars[char_rnf])[names(corlow)]
      corup <- unlist(LowUp$upper$corrPars[char_rnf])
      at_lower_bound <- which((corpars-corlow)/(corup-corlow)<1e-5)
      if ( length(at_lower_bound) && 
           length(which(grepl("trRho", names(corpars))))>1L && ## at least two rho's 
           #  (case with one rho at the bound would deserve another message; NB occurs in the vignette...)
           any(grepl("trRho", names(at_lower_bound)))) {
        locmess <- paste0("The estimate(s) of some 'rho' scale parameter(s) = their lower bound in numerical optimization.\n",
                          "This suggests that the expected value of some ranef(s) should be non-zero,\n ",
                          "but that the model term(s) for intercept(s) do not account for that.")
      }
    }
  }
  locmess
}


