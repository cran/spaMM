# processed$hyper_info <- .preprocess_hyper(.) provides a list used by .calc_inits() -> .calc_inits_hyper() and by .expand_hyper()
# .calc_inits_hyper() fills init.optim$hyper <- structure(hyper,map=hyper_info$map,ranges=hyper_info$rangess, trL_skeleton=trL_skeleton), 
#                           whic his a nested list of params with attr. The nested list  serves as a template for .makeLowerUpper.
# expand_hyper() is used (notably) by HLCor.obj to fill $corrPars and $trLambda from [ranPars derived from ranefParsVec] and processed$hyper_info


MRF <- function(...) {
  canonize <- function(corrPars_rd, cP_type_rd, moreargs_rd, checkComplete, ...) {
    if (!is.null(corrPars_rd$trKappa)) { ## 
      corrPars_rd$kappa <- .kappaInv(corrPars_rd$trKappa,KAPPAMAX=moreargs_rd$KAPPAMAX) 
      corrPars_rd$trKappa <- NULL
      cP_type_rd$kappa <- cP_type_rd$trKappa 
      cP_type_rd$trKappa <- NULL
    }
    kappa <- corrPars_rd$kappa
    if (is.null(kappa) && checkComplete) {
      stop("kappa missing from ranPars (or correlation model mis-identified).")
    }
    return(list(corrPars_rd=corrPars_rd, cP_type_rd=cP_type_rd))
  }
  calc_inits <- function(inits, char_rd, moreargs_rd, user.lower, user.upper, optim.scale, ranFix, init.optim, ...) {
    inits <- .calc_inits_MRF(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                             user.lower=user.lower,user.upper=user.upper,optim.scale=optim.scale,KAPPAMAX=moreargs_rd$KAPPAMAX,char_rd=char_rd)
    # # Nugget: remains NULL through all computations if NULL in init.optim
    if (is.null(.get_cP_stuff(ranFix,"Nugget",which=char_rd))) { ## new spaMM3.0 code
      inits$init$corrPars[[char_rd]] <- .modify_list(inits$init$corrPars[[char_rd]],
                                                     list(Nugget=.get_cP_stuff(init.optim,"Nugget",which=char_rd)))
    }
    return(inits)
  }
  calc_corr_from_dist <- function(ranFix,char_rd,distmat,...) { stop("This should not be called for MRF terms") }
  #
  calc_moreargs <- function(KAPPAMAX, ...) {
    moreargs_rd <- list(KAPPAMAX=KAPPAMAX) 
    return(moreargs_rd)
  }
  make_new_corr_list <- function(object, old_char_rd, control_dist_rd, geonames, newuniqueGeo, olduniqueGeo, which_mats, make_scaled_dist, new_rd) {
    ## hmmm is that useful ?
  }
  makeLowerUpper <- function() {
    ## hmmm is that useful ?
  }
  # : changes the parent.env of all the member functions.
  structure(list(corr_family="MRF",
                 names_for_post_process_parlist= c("kappa"),
                 canonize=canonize,
                 calc_inits=calc_inits,
                 calc_corr_from_dist=calc_corr_from_dist,
                 calc_moreargs=calc_moreargs
                 #,make_new_corr_list=make_new_corr_list
  ),
  class="corr_family")
}

.MRFcrossfactor <- function(xsteps, ysteps=xsteps, margin, kappa) {
  xstwm <- xsteps +2L*margin
  ystwm <- ysteps +2L*margin
  indices <- expand.grid(seq(xstwm), seq(ystwm))
  nc_I <- nrow(indices)
  locpos0 <- seq(nc_I)
  locpos1 <- seq(xstwm*(ystwm-1L))
  locpos2 <- locpos0[-xstwm*seq(ystwm)]
  B <- sparseMatrix(i=c(locpos0,locpos1,locpos2), j=c(locpos0,locpos1+xstwm,locpos2+1L), 
                    x=c(rep(4+kappa^2,nc_I),rep(-1,length(locpos1)+length(locpos2))), symmetric =TRUE)
  # B <- as(B,"symmetricMatrix") # dsCMatrix ## potentially slow and only useful if speeding a crossprod...
  return(B) ## such that Q= crossprod(B)
  # or :
  B <- new("dgCMatrix",i= 0:(nc_I-1L), p=c(0L:(nc_I)), Dim=c(nc_I,nc_I),x=rep(4+kappa^2,nc_I))
  locpos <- seq(xstwm*(ystwm-1L))
  B[cbind(locpos,locpos+xstwm)] <- -1
  locpos <- seq(xstwm*(ystwm-1L))
  B[cbind(locpos+xstwm, locpos)] <- -1
  locpos <- seq(nc_I)[-xstwm*seq(ystwm)]
  B[cbind(locpos, locpos+1L)] <- -1
  B[cbind(locpos+1L, locpos)] <- -1
  #B <- as(B,"symmetricMatrix") # dsCMatrix
}

.Wendland <- function(d) {
  res <- 0*d
  islowd <- d<1
  lowd <- d[islowd]
  res[islowd] <- (35*lowd^2 + 18*lowd + 3)*(1-lowd)^6 /3 
  return(res)
}

.to_grid_coord <- function(uniqueScal, ranges, pars) {
  for (dt in seq(ncol(uniqueScal))) {
    rge <- ranges[[dt]]
    steplen <- diff(rge)/(pars[["nd"]]-1L)
    uniqueScal[,dt] <- uniqueScal[,dt]/steplen - rge[1]/steplen ## from 0 to pars$nd-1L
  }
  return(uniqueScal)
}

.calc_AMatrix_MRF <- function(term, data, dist.method, old_AMatrix_rd=NULL, scale) {
  pars <- attr(attr(term,"type"),"pars")
  blob <- .get_dist_nested_or_not(term, data=data, distMatrix=NULL, uniqueGeo=NULL, 
                                  dist.method = dist.method, needed=c(uniqueGeo=TRUE),  geo_envir=NULL)
  uniqueScal <- blob$uniqueGeo
  rownames(uniqueScal) <- seq(nrow(uniqueScal)) ## hmmmf... this must match dataordered_levels_blob <- .calc_dataordered_levels(...) but...
  if ( ! is.null(old_AMatrix_rd)) {
    ranges <- attr(old_AMatrix_rd,"ranges")
    grid_arglist <- attr(old_AMatrix_rd,"grid_arglist") ## reconstruct the original grid
    scale <- attr(old_AMatrix_rd,"scale")
  } else {
    ranges <- lapply(uniqueScal, range)
    grid_arglist <- vector("list",length(ranges))
    for (dt in seq(ncol(uniqueScal))) grid_arglist[[dt]] <- (-pars$m):(pars$nd+pars$m-1L)
    #
  }
  uniqueScal <- .to_grid_coord(uniqueScal,ranges=ranges,pars=pars)
  grid <- expand.grid(grid_arglist)
  #
  scaldistmat <- proxy::dist(x=uniqueScal,y=grid,method=dist.method)
  if ( ! is.null(old_AMatrix_rd)) {  
    colnames(scaldistmat) <- colnames(old_AMatrix_rd) 
    return(as(.Wendland(scaldistmat[]/scale),"sparseMatrix")) ## dividing the matrix not most economical but clearer # [] convert crossdist to matrix
  } else {
    colnames(scaldistmat) <- apply(grid,1L,paste0,collapse=":") ## provide colnames(ZA), expected at least by ranef.HLfit 
    return(structure(as(.Wendland(scaldistmat[]/scale),"sparseMatrix"), ## dividing the matrix not most economical but clearer # [] convert crossdist to matrix
                     grid_arglist=grid_arglist, ranges=ranges,scale=scale)) ## attributes for quick reconstruction in prediction
  }  
}

.assign_AMatrices_MRF <- function(corr_info, Zlist, exp_barlist, data, control_dist,
                                  scale=2.5) { ## scale value from Nychka et al 2015 http://dx.doi.org/10.1080/10618600.2014.914946, p. 584
  corr_types <- corr_info$corr_types
  isMRF <- (corr_types=="MRF")
  if (any(isMRF, na.rm=TRUE)) {
    if (is.null(corr_info$AMatrices)) corr_info$AMatrices <- vector("list",length(Zlist))
    for (rd in which(isMRF)) {
      corr_info$AMatrices[[rd]] <- .calc_AMatrix_MRF(term=exp_barlist[[rd]], data=data, 
                                                     dist.method=control_dist[[rd]][["dist.method"]], scale=scale)
    }
  }
}

.get_new_AMatrices <- function(object, new_mf_ranef) {
  amatrices <- attr(object$ZAlist,"AMatrices")
  exp_spatial_terms <- attr(object$ZAlist,"exp_spatial_terms")
  isMRF <- (attr(exp_spatial_terms,"type") == "MRF")
  for (rd in which(isMRF)) {
    amatrices[[rd]] <- .calc_AMatrix_MRF(term=exp_spatial_terms[[rd]], data=new_mf_ranef, 
                                         dist.method=.get_control_dist(object,rd)$dist.method, old_AMatrix_rd = amatrices[[rd]])
  }
  return(amatrices)
}

#"multiSAR"


.expand_multiSAR <- function(bar, levels, margin, coarse=10L, ...) {
  nodes_per_dim <- (coarse)*2^seq(levels) # +2*arglist$margin
  bar <- paste0(bar, ", nd=", nodes_per_dim, ", m=",margin)
  as.formula(paste("~ ( MRF(", paste(bar, collapse = ") + MRF("), ") )"))[[2]]
}


# called by .preprocess_formula():

.expand_multiSARs <- local({
  hyper_info <- NULL
  function (term,init=TRUE) {
    if (init) hyper_info <<- numeric(0)
    if (!is.name(term) && is.language(term)) {
      if (term[[1]] == as.name("(")) {
        term[[2]] <- .expand_multiSARs(term[[2]],init=FALSE)
      }
      stopifnot(is.call(term))
      term1 <- as.character(term[[1]])
      if (term1 == "multiSAR") {
        args <- paste(paste(names(term),"=",term)[-c(1,2)], collapse=",")
        arglist <- eval(parse(text = paste("list(",args,")")))
        arglist$bar <- paste0(deparse(term[[2]]))
        hyper_info <<- c(hyper_info,arglist$levels)
        return(do.call(".expand_multiSAR", arglist)) # fn call => uses standard R argument matching 
      }
      term[[2]] <- .expand_multiSARs(term[[2]],init=FALSE)
      if (length(term) == 3) 
        term[[3]] <- .expand_multiSARs(term[[3]],init=FALSE)
    }
    if (init && length(hyper_info)) {
      return(structure(term, hyper_info=hyper_info))
    } else return(term)
  }
})

.calc_inits_hyper <- function(inits, hyper_info) {
  if ( ! is.null(hyper_info)) {
    init.optim <- inits$init.optim
    ranges <- hyper_info$ranges
    hyper <- hyper_info$template
    for (rg_it in seq_along(ranges)) {
      firstrand <- (rd_range <- ranges[[rg_it]])[1L]
      # names (by as.character) required for .modify_list() to work
      hyper[[as.character(rg_it)]] <- list(hy_trL=init.optim$trLambda[firstrand], hy_trK=init.optim$corrPars[[firstrand]]$trKappa) 
      init.optim$trLambda[rd_range] <- NA 
      for (rd in rd_range) { init.optim$corrPars[[rd]]$trKappa <- NULL }
    }
    trL_skeleton <- init.optim$trLambda ## with NA's
    init.optim$trLambda <- init.optim$trLambda[!is.na(init.optim$trLambda)]
    inits[["init"]]$lambda <- inits[["init"]]$lambda[!is.na(init.optim$trLambda)] ## as this is used as a template for lower$trLambda...
    init.optim$hyper <- structure(hyper,map=hyper_info$map,ranges=hyper_info$ranges, trL_skeleton=trL_skeleton)
    inits[["init.optim"]] <- init.optim
  }
  return(inits)
}

.preprocess_hyper <- function(processed) { ## needs $predictor, $ZAlist
  hyper_info <- attr(processed$predictor,"hyper_info") ## a vector of number of levels
  if (len <- length(hyper_info)) {
    ranges <- vector("list",len)
    hyper_map <- rep(NA_integer_,length(processed$ZAlist)) ## maps rd to elements of hyper
    MRF_ranges <- cumsum(c(0,hyper_info))
    which_MRF <- which(attr(processed$ZAlist,"exp_ranef_types") =="MRF")
    for (it in seq_len(len)) { 
      rd_range <- (which_MRF[MRF_ranges[it]+1L]):(MRF_ranges[it+1L])
      hyper_map[rd_range] <- it 
      ranges[[it]] <- rd_range ## inverse map
    }
    template <- structure(vector("list",length(ranges)),names=paste(seq(length(ranges))))
    return(list(map=hyper_map,ranges=ranges, template= template) )
  }
}

.expand_hyper <- function(ranPars, trL_skeleton, processed) {
  if (!is.null(hyper_info <- processed$hyper_info)) {
    ranges <- hyper_info$ranges
    for (rg_it in seq_along(ranges)) {
      hyper_el <- ranPars$hyper[[rg_it]]
      if ("hy_trL" %in% names(hyper_el)) {
        rd_range <- ranges[[rg_it]] 
        hy_lam <- .dispInv(hyper_el$hy_trL)
        lambdas <- hy_lam / (2^(2*seq_len(length(rd_range)))) ## 
        trL_skeleton[rd_range] <- .dispFn(lambdas) ## over-writes temporary NaN
      }
      if ("hy_trK" %in% names(hyper_el)) {
        rd_range <- ranges[[rg_it]] 
        for(rd in rd_range) {
          ranPars$corrPars[[rd]]$trKappa <- hyper_el$hy_trK ## same SAR parameter for all levels
        }
      }
    }
    ranPars$trLambda <- trL_skeleton 
    #ranPars$hyper <- NULL 
    #rpType$hyper <- NULL
  } 
  return(ranPars) ## with trLambda, $corrPars modified
}