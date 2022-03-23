# 'Canned' correlation families

.kappa_range_factor <- function(mesh, IMRF_design, EEV_required) {
  range_factor <- c(diff(range(mesh$loc[,1])),diff(range(mesh$loc[,2])))
  range_factor <- sqrt(sum(range_factor*range_factor)) # avoid computation of all distances (or on the hull only...)
  if (is.null(IMRF_design)) { # alpha not prefixed in the fitted model => IMRF_design not precomputed
    .inla.spde2.matern <- get("inla.spde2.matern", asNamespace("INLA"), inherits = FALSE)
    if (inherits(.inla.spde2.matern,"try-error")) stop("INLA must be installed in order to design the 'MaternIMRFa' covariance.")
    local_design <- .inla.spde2.matern(mesh, alpha=1)
    e1 <- extreme_eig(as(local_design$param.inla$M2,"dgCMatrix"), symmetric=TRUE, required=EEV_required)[1L]
  } else {
    e1 <- extreme_eig(as(IMRF_design$param.inla$M1,"dgCMatrix"), symmetric=TRUE, required=EEV_required)[1L]
    if ( ! is.null(e1)) {
      e1 <- max(e1,
                extreme_eig(as(IMRF_design$param.inla$M2,"dgCMatrix"), symmetric=TRUE, required=EEV_required)[1L])
    }
  }
  if ( ! is.null(e1)) range_factor <- min(range_factor, 0.01/sqrt((e1 - 1e5 * 0)/(-1 + 1e5 + e1 - 1e5 * 0)))
  range_factor
}

.make_new_corr_lists_IMRFs <- function(newLv_env, which_mats, object, ranefs, new_rd, old_rd) {
  ## in this case a new A matrix (by .get_new_AMatrices()) must be computed (elsewhere: it goes in newZAlist)
  ## but the corr matrix between the nodes is unchanged as node positions do not depend on response position.
  ## => In this special case, the correlation between new and old locations is the correlation between old locations.
  if (which_mats$no || which_mats$nn[new_rd]) uuColdold <- .tcrossprod(object$strucList[[old_rd]], perm=TRUE) # perm=TRUE means for ranefs permuted as in ZA
  if (which_mats$no) newLv_env$cov_newLv_oldv_list[[new_rd]] <- 
      structure(uuColdold,
                ranefs=ranefs[[new_rd]]
      ) ## always necess for .compute_ZAXlist(XMatrix = cov_newLv_oldv_list, ZAlist = newZAlist)
  ## correct code but should be useless:
  # warning("Suspect code: covariance matrix computation requested for IMRF term.")
  # if (which_mats$nn[new_rd]) {
  #   newLv_env$cov_newLv_newLv_list[[new_rd]] <- uuColdold
  # } else newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- diag(x=uuColdold)
}

MaternIMRFa <- function(mesh, tpar=c(alpha=1.25,kappa=0.1), fixed=NULL) { 
  
  force(mesh)
  ZAcolnames <- paste0(mesh$loc[,1],":",mesh$loc[,2]) # match between columns of A and rows of Q
  # Arownames need to be recomputed for each new data (fit vs predict)
  
  .inla.spde2.matern <- get("inla.spde2.matern", asNamespace("INLA"), inherits = FALSE)
  if (inherits(.inla.spde2.matern,"try-error")) stop("INLA must be installed in order to design the 'MaternIMRFa' covariance.")
  
  IMRF_design <- NULL
  if ((  ! is.null(fixed) ) && 
      (! is.na(fixed["alpha"]))) IMRF_design <- .inla.spde2.matern(mesh, alpha=fixed["alpha"])
  
  Cf <- function(parvec) {
    alpha <- parvec[1L]
    if (is.null(IMRF_design)) IMRF_design <- .inla.spde2.matern(mesh, alpha=alpha) # it as been precomputed in the parent envir if alpha was fixed
    kappa <- parvec[2L]
    Qmat <- .calc_IMRF_Qmat(pars=list(model=IMRF_design, SPDE_alpha=alpha), 
                            grid_arglist=NULL, 
                            kappa=kappa, test=FALSE)
    rownames(Qmat) <- ZAcolnames
    #list(precision=Qmat) # returns a precision matrix
    Qmat
  }
  
  # .calc_AMatrix_IMRF has typical arguments (term=exp_spatial_terms[[rd]], data=newdata, 
  #                   dist.method=.get_control_dist(object,char_rd)$dist.method, 
  #                  old_AMatrix_rd = amatrices[[char_rd]])
  # 'term' appears as needed in Af() here as in .calc_AMatrix_IMRF
  Af <- function(newdata, 
                 term,  
                 fit.=FALSE,
                 ... # assumed by the call by .assign_AMatrices_corrFamily() 
  ) { # the 'data' locations are potentially distinct from what has been used to create the IMRF_design
    # and generally distinct from the mesh locations!
    # and the positions in the fitted data are the cols of Z=rows of A
    # The positions on the mesh are the cols of A = cols fo ZA = rows of Q =mespos
    if (is.null(newdata)) {
      stop("'newdata' needed to determine the data locations (not for names, but for .spaMM_spde.make.A(.,loc) )") 
    } else {
      blob <- .get_dist_nested_or_not(term, data=newdata, distMatrix=NULL, uniqueGeo=NULL, 
                                      dist.method = NULL, needed=c(uniqueGeo=TRUE),  geo_envir=NULL)
      uniqueScal <- blob$uniqueGeo
      uniqueScal_levels_blob <- .as_factor(txt=paste(colnames(uniqueScal), collapse="+"), 
                                           mf=uniqueScal, type=.spaMM.data$options$uGeo_levels_type) # to match the call in .calc_Zmatrix()
      Arownames <- uniqueScal_levels_blob$factor # not levels(.) which are alphanumerically ordered !  
      # this should match the colnames of the ZMatrix matching the 'newdata'. We could check that, 
      # but it's not clear when this could not be the case (except forgetting type=.spaMM.data$options$uGeo_levels_type...), so the B/C of the  check seems low.
    }
    Amatrix <- .spaMM_spde.make.A(mesh=mesh, points=as.matrix(uniqueScal))
    if (fit. && any(rowSums(Amatrix)<0.99)) message("Some 'data' locations appear out of the mesh. Mismatched inputs, or strong mesh cutoff?") # check this only for the fit.
    # : If all of them are out of the mesh, as the structure of ZA is used to determine initial lambda value, I had a NaN there, 
    #   so no outer init lambda where expected => bug)
    # Amatrix rows are, with the default arguments, ordered as uniqueScal rows
    rownames(Amatrix) <- Arownames # names from coordinates in the data (notably, fitme(., data)...)
    colnames(Amatrix) <- ZAcolnames # names from coordinates of the mesh points
    return(Amatrix)
  }
  
  calc_moreargs <- function(corrfamily, ...) {
    range_factor <- .kappa_range_factor(mesh, IMRF_design, EEV_required=FALSE) # where IMRF_design may still be NULL
    list(
      lower=c(alpha=1,kappa=0.01/range_factor), 
      init=c(alpha=1.15,kappa=1/range_factor), 
      upper=c(alpha=2,kappa=200/range_factor) 
    )
  }
  
  make_new_corr_lists <- function(newLv_env, which_mats, object, ranefs, new_rd, old_rd, ...) {
    .make_new_corr_lists_IMRFs(newLv_env=newLv_env, which_mats=which_mats, object=object, ranefs=ranefs, new_rd=new_rd, old_rd=old_rd)
  }
  
  list(Cf=Cf, tpar=tpar, type="precision", Af=Af, fixed=fixed, calc_moreargs=calc_moreargs, 
       make_new_corr_lists=make_new_corr_lists,
       levels_type=.spaMM.data$options$uGeo_levels_type,
       need_Cnn=FALSE,
       tag="MaternIMRFa") # the mesh is in the environment of the functions
}

gridIMRF <- function(..., fixed=NULL) {  # internal family constructor, not descriptor => could use many argument, does not even need tpar
  
  gridpars <- list(...)
  # as in .process_IMRF_bar()
  if (is.null(gridpars$no)) gridpars$no <- TRUE
  if (is.null(gridpars$ce)) gridpars$ce <- TRUE
  # and per the IMRF() doc, it expects $margin (allowed shorthand: $m) and $nd
  
  tpar <- c(kappa=0.1)
  
  grid_arglist <- precnames <- NULL
  #
  initialize <- function(Zmatrix, ...) { 
    
    Zlevels <- colnames(Zmatrix)
    uniqueScal <- strsplit(Zlevels,":") 
    uniqueScal <- do.call(rbind,uniqueScal)
    class(uniqueScal) <- "numeric"
    
    nc <- ncol(uniqueScal)
    ranges <- diffs <- grid_arglist <- vector("list",nc) ## typically nc=2
    for (dt in seq(nc)) { 
      ranges[[dt]] <- range(uniqueScal[,dt]) 
      diffs[[dt]] <- diff(ranges[[dt]])
    }
    widedim <- which.max(unlist(diffs))
    steplen <- diffs[[widedim]]/(gridpars[["nd"]]-1L)
    origin <- numeric(nc)
    for (dt in seq(nc)) {
      if (gridpars$ce) {
        origin[dt] <- ranges[[dt]][1L] - (diffs[[dt]] %% steplen)/2 ## or using trunc... ## 0 for the wide dimension
      } else origin[dt] <- ranges[[dt]][1L] # ce=FALSE reproduces the LatticeKrig results
      nd_dt <- (diffs[[dt]] %/% steplen) +1L ## central grid node for dim dt
      grid_arglist[[dt]] <- (-gridpars$m):(nd_dt+gridpars$m-1L) ## sequence of nodes in grid coordinates for dim dt
    }
    
    grid_arglist <<- grid_arglist
    
    grid <- expand.grid(grid_arglist) ## INTEGER GRID
    precnames <<-  apply(grid,1L,paste0,collapse=":") ## Should match colnames(ZA) 
    
  }
  
  Cf <- function(parvec) { # a case wher parvec is actually a list
    Qmat <- .calc_IMRF_Qmat(pars=gridpars, grid_arglist=grid_arglist, kappa=unlist(parvec), test=FALSE)
    rownames(Qmat) <- precnames
    Qmat
  }
  
  Af <- function(newdata, 
                 term,  
                 fit.=FALSE,
                 dist.method, 
                 old_AMatrix_rd=NULL, 
                 scale,
                 ... # assumed by the call by .assign_AMatrices_corrFamily() 
  ) {
    # .calc_AMatrix_IMRF has typical arguments (term=exp_spatial_terms[[rd]], data=newdata, 
    #                   dist.method=.get_control_dist(object,char_rd)$dist.method, 
    #                  old_AMatrix_rd = amatrices[[char_rd]])
    # 'term' appears as needed in Af() here as in .calc_AMatrix_IMRF
    .calc_AMatrix_IMRF(term, # its attr(.,"pars") carries grid parameters
                       data=newdata, # only for .get_dist_nested_or_not()  
                       dist.method=dist.method, old_AMatrix_rd=old_AMatrix_rd, 
                       scale=scale,
                       pars=gridpars,
                       fit.=fit.
    )
  }
  
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
  
  calc_inits <- function(inits, char_rd, moreargs_rd, user.lower, user.upper, optim.scale, init.optim, ...) {
    inits <- .calc_inits_IMRF(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                              user.lower=user.lower,user.upper=user.upper,optim.scale=optim.scale,
                              moreargs_rd=moreargs_rd,char_rd=char_rd)
    ### There is no documented Nugget for IMRF models...
    # # # Nugget: remains NULL through all computations if NULL in init.optim
    # if (is.null(.get_cP_stuff(inits$ranFix,"Nugget",which=char_rd))) { 
    #   inits$init$corrPars[[char_rd]] <- .modify_list(inits$init$corrPars[[char_rd]],
    #                                                  list(Nugget=.get_cP_stuff(init.optim,"Nugget",which=char_rd)))
    # }
    
    # inits$init$corrPars[[char_rd]] <- unlist(inits$init$corrPars[[char_rd]]) # $Cf() expects a vector (documented); .calc_inits_IMRF() produced a list. 
    
    return(inits)
  }
  
  calc_moreargs <- function(KAPPAMAX, ...) {
    # * Always on transformed scale;
    # * lists rather than vectors, as for the "IMRF-keyword" ranef
    moreargs_rd <- list(KAPPAMAX=KAPPAMAX,
                        lower=list(trKappa= .kappaFn(kappa=1e-4,KAPPAMAX=KAPPAMAX)),
                        upper=list(trKappa= .kappaFn(kappa=KAPPAMAX-1e-6,KAPPAMAX=KAPPAMAX))) # optimization should not try to approach infinity)) 
    return(moreargs_rd)
  }

  make_new_corr_lists <- function(newLv_env, which_mats, object, ranefs, new_rd, old_rd, ...) {
    .make_new_corr_lists_IMRFs(newLv_env=newLv_env, which_mats=which_mats, object=object, ranefs=ranefs, new_rd=new_rd, old_rd=old_rd)
  }
  
  list(Cf=Cf, tpar=tpar, type="precision", Af=Af, fixed=fixed, initialize=initialize, canonize=canonize, calc_moreargs=calc_moreargs, 
       calc_inits=calc_inits, 
       normIMRF= gridpars$no, # for normalization of ZA
       make_new_corr_lists=make_new_corr_lists,
       tag="gridIMRF") # the mesh is in the environment of the functions
}




