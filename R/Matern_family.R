Matern <- function(...) {
  canonize <- function(corrPars_rd, cP_type_rd, moreargs_rd, checkComplete, ...) {
    if (!is.null(corrPars_rd$trNu)) { ## either we have nu,rho or trNu,trRho 
      corrPars_rd$nu <- .nuInv(corrPars_rd$trNu,NUMAX=moreargs_rd$NUMAX) ## before trRho is removed...
      corrPars_rd$trNu <- NULL
      cP_type_rd$nu <- cP_type_rd$trNu 
      cP_type_rd$trNu <- NULL
    }
    nu <- corrPars_rd$nu
    if (is.null(nu) && checkComplete) {
      stop("nu missing from ranPars (or correlation model mis-identified).")
    }
    if ( ! is.null(corrPars_rd$trRho)) { ## assuming a single trRho with possibly several elements
      corrPars_rd$rho <- .rhoInv(corrPars_rd$trRho,RHOMAX=moreargs_rd$RHOMAX)  
      corrPars_rd$trRho <- NULL
      cP_type_rd$rho <- cP_type_rd$trRho
      cP_type_rd$trRho <- NULL
    } ## else there may simply be rho rather than trRho (including for adjacency model through optim procedure !)
    rho <- corrPars_rd$rho
    if (is.null(rho) && checkComplete) stop("rho missing from ranPars.")
    return(list(corrPars_rd=corrPars_rd, cP_type_rd=cP_type_rd))
  }
  calc_inits <- function(inits, char_rd, moreargs_rd, user.lower, user.upper, optim.scale, init.optim, ...) {
    inits <- .calc_inits_geostat_rho(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                     user.lower=user.lower,user.upper=user.upper,
                                     maxrange=moreargs_rd$maxrange,optim.scale=optim.scale,RHOMAX=moreargs_rd$RHOMAX,char_rd=char_rd)
    inits <- .calc_inits_nu(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                            control_dist_rd=moreargs_rd$control.dist, optim.scale=optim.scale, NUMAX=moreargs_rd$NUMAX,char_rd=char_rd)
    # Nugget: remains NULL through all computations if NULL in init.optim
    if (is.null(.get_cP_stuff(inits$ranFix,"Nugget",which=char_rd))) { ## new spaMM3.0 code
      inits$init$corrPars[[char_rd]] <- .modify_list(inits$init$corrPars[[char_rd]],
                                                     list(Nugget=.get_cP_stuff(init.optim,"Nugget",which=char_rd)))
    }
    return(inits)
  }
  calc_corr_from_dist <- function(ranFix,char_rd,distmat,...) {
    nu <- .get_cP_stuff(ranFix,"nu",which=char_rd)
    #if (is.null(nu)) stop("NULL 'nu' here...") 
    # if (is.null(nu)) nu <- .nuInv(.get_cP_stuff(ranFix,"trNu",which=char_rd),  # not sure test is ever TRUE
    #                               NUMAX =.getPar(attr(ranFix,"moreargs"),"NUMAX"))
    corr_mat <- MaternCorr(nu=nu,
                           Nugget=.get_cP_stuff(ranFix,"Nugget",which=char_rd),
                           d=distmat)  ## so that rho=1 in MaternCorr
    return(corr_mat)
  }
  #
  calc_moreargs <- function(fixed, char_rd, init.optim, processed, rd, control_dist, NUMAX, LDMAX, ...) {
    rho_size <- max(length(fixed$corrPars[[char_rd]][["rho"]]),length(init.optim$corrPars[[char_rd]][["rho"]]))
    range_info_blob <- .calc_range_info(rho_size, processed, rd, control_dist[[char_rd]]) 
    RHOMAX <- 1.1*30*range_info_blob$nbUnique/range_info_blob$maxrange## matches init$rho in calc_inits() 
    control_dist[[char_rd]]$rho.mapping <- range_info_blob$rho_mapping 
    moreargs_rd <- list(maxrange=range_info_blob$maxrange, RHOMAX=RHOMAX,
                                NUMAX=NUMAX, LDMAX=LDMAX,## not variable...
                                nbUnique=range_info_blob$nbUnique, control.dist=control_dist[[char_rd]]## needed for LUarglist, not calc_inits
    ) 
    return(moreargs_rd) ## with control_dist with possibly modified rho_mapping
  }
  # This fn is not used (it is an unfinished attempt), so it does not matter the Cauchy code is "not even wrong".
  calc_cov_info_mat <- function(control.dist, char_rd, spatial_term, corr_type, rho, processed, rd, ranPars) {
    control_dist_rd <- control.dist[[char_rd]]
    txt <- paste(c(spatial_term[[2]][[3]])) ## the RHS of the ( . | . ) # c() to handle very long RHS
    if (length(grep("%in%",txt))) {
      stop(paste0("(!) ",corr_type,"( . | <coord> %in% <grp>) is not yet handled."))
    } 
    msd.arglist <- list(rho = rho)
    msd.arglist$`dist.method` <- control_dist_rd$`dist.method` ## may be NULL
    if (length(rho)>1L) {
      geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("uniqueGeo"), 
                                 dist_method_rd=control_dist_rd$dist.method)
      msd.arglist$uniqueGeo <- geo_envir$uniqueGeo
      msd.arglist$`rho.mapping` <- control_dist_rd$`rho.mapping` ## may be NULL
    } else {
      geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("distMatrix"), 
                                 dist_method_rd=control_dist_rd$dist.method)
      msd.arglist$distMatrix <- geo_envir$distMatrix   
    }
    cov_info_mat <- do.call("make_scaled_dist",msd.arglist)
    ## at this point if a single location, dist_mat should be dist(0) and make_scaled_dist was modified to that effect
    if ( nrow(cov_info_mat)>1 ) { ## >1 locations
      Nugget <- .get_cP_stuff(ranPars,"Nugget",which=char_rd)
      if (corr_type == "Matern") {
        nu <- .get_cP_stuff(ranPars,"nu",which=char_rd)
        if (is.null(nu)) nu <- .nuInv(.get_cP_stuff(ranPars,"trNu",which=char_rd), NUMAX =.getPar(attr(ranPars,"moreargs"),"NUMAX")) # not sure test is ever TRUE
        cov_info_mat <- MaternCorr(nu=nu, Nugget=Nugget, d=cov_info_mat)        
      } else { ## ,"Cauchy"
        longdep <- .get_cP_stuff(ranPars,"longdep",which=char_rd)
        #if (is.null(longdep)) longdep <- .longdepInv(.get_cP_stuff(ranPars,"trLongdep"), LDMAX =.getPar(attr(ranPars,"moreargs"),"LDMAX")) # not sure test is ever TRUE
        shape <- .get_cP_stuff(ranPars,"shape",which=char_rd)
        cov_info_mat <- CauchyCorr(shape=shape, longdep=longdep, Nugget=Nugget, d=cov_info_mat)        
      }
      # no rho because the MaternCorr input will be an already scaled distance 'cov_info_mat'
    } 
    return(cov_info_mat)
  }
  
  make_new_corr_list <- function(object, old_char_rd, control_dist_rd, geonames, newuniqueGeo, olduniqueGeo, which_mats, make_scaled_dist, new_rd) {
    ### rho only used to compute scaled distances
    rho <- .get_cP_stuff(object$ranFix,"rho", which=old_char_rd)
    if ( ! is.null(rho_mapping <- control_dist_rd$rho.mapping) 
         && length(rho)>1L ) rho <- .calc_fullrho(rho=rho,coordinates=geonames,rho_mapping=rho_mapping)
    ## rows from newuniqueGeo, cols from olduniqueGeo:
    msd.arglist <- list(uniqueGeo=newuniqueGeo,uniqueGeo2=olduniqueGeo,
                        rho=rho,return_matrix=TRUE)
    if ( ! is.null(dist.method <- control_dist_rd$dist.method)) msd.arglist$dist.method <- dist.method ## make_scaled_dist does not handle NULL
    resu <- list()
    if (which_mats$no) resu$uuCnewold <- do.call(make_scaled_dist,msd.arglist) ## ultimately allows products with Matrix ## '*cross*dist' has few methods, not even as.matrix
    if (which_mats$nn[new_rd])  {
      msd.arglist$uniqueGeo2 <- NULL
      if (nrow(msd.arglist$uniqueGeo)==1L) {
        resu$uuCnewnew <- matrix(0) ## trivial distance matrix for single point
      } else resu$uuCnewnew <- do.call(make_scaled_dist,msd.arglist) 
    }
    return(resu)
  }
  makeLowerUpper <- function() {
    ## hmmm is that useful ?
  }
  #parent.env(environment(calc_corr_from_dist)) <- environment(spaMM::MaternCorr) ## parent = <environment: namespace:stats>
  # : changes the parent.env of all the member functions.
  structure(list(corr_family="Matern",
                 names_for_post_process_parlist= c("rho","nu","Nugget"),
                 canonize=canonize,
                 calc_inits=calc_inits,
                 calc_corr_from_dist=calc_corr_from_dist,
                 #calc_cov_info_mat=calc_cov_info_mat,
                 #make_new_corr_list=make_new_corr_list,
                 calc_moreargs=calc_moreargs ),
            class="corr_family")
}

Cauchy <- function(...) {
  canonize <- function(corrPars_rd, cP_type_rd, checkComplete, moreargs_rd, ...) {
    if (!is.null(corrPars_rd$trLongdep)) { ## either we have longdep,rho or trLongdep,trRho 
      corrPars_rd$longdep <- .longdepInv(corrPars_rd$trLongdep,LDMAX=moreargs_rd$LDMAX)
      corrPars_rd$trLongdep <- NULL
      cP_type_rd$longdep <- cP_type_rd$trLongdep 
      cP_type_rd$trLongdep <- NULL
    }
    longdep <- corrPars_rd$longdep
    if (is.null(longdep) && checkComplete) {
      stop("longdep missing from ranPars (or correlation model mis-identified).")
    }
    if ( ! is.null(corrPars_rd$trRho)) { ## assuming a single trRho with possibly several elements
      corrPars_rd$rho <- .rhoInv(corrPars_rd$trRho,RHOMAX=moreargs_rd$RHOMAX)  
      corrPars_rd$trRho <- NULL
      cP_type_rd$rho <- cP_type_rd$trRho
      cP_type_rd$trRho <- NULL
    } ## else there may simply be rho rather than trRho (including for adjacency model through optim procedure !)
    rho <- corrPars_rd$rho
    if (is.null(rho) && checkComplete) {
      stop("rho missing from ranPars.")
    }
    return(list(corrPars_rd=corrPars_rd, cP_type_rd=cP_type_rd))
  }
  #
  calc_inits <- function(inits, user.lower, user.upper, moreargs_rd, optim.scale, char_rd, init.optim, ...) {
    inits <- .calc_inits_geostat_rho(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                     user.lower=user.lower,user.upper=user.upper,
                                     maxrange=moreargs_rd$maxrange,optim.scale=optim.scale,RHOMAX=moreargs_rd$RHOMAX,char_rd=char_rd)
    inits <- .calc_inits_cauchy(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                control_dist_rd=moreargs_rd$control.dist, optim.scale=optim.scale, LDMAX=moreargs_rd$LDMAX,char_rd=char_rd)
    # Nugget: remains NULL through all computations if NULL in init.optim
    if (is.null(.get_cP_stuff(inits$ranFix,"Nugget",which=char_rd))) { ## new spaMM3.0 code
      inits$init$corrPars[[char_rd]] <- .modify_list(inits$init$corrPars[[char_rd]],
                                                     list(Nugget=.get_cP_stuff(init.optim,"Nugget",which=char_rd)))
    }
    return(inits)
  }
  #
  calc_corr_from_dist <- function(ranFix, char_rd, distmat,...) {
    #longdep <- .get_cP_stuff(ranFix,"longdep",which=char_rd)
    # #if (is.null(longdep)) longdep <- .longdepInv(.get_cP_stuff(ranFix,"trLongdep"), LDMAX =.getPar(attr(ranFix,"moreargs"),"LDMAX")) # not sure test is ever TRUE
    corr_mat <- CauchyCorr(shape=.get_cP_stuff(ranFix,"shape",which=char_rd),
                           longdep=.get_cP_stuff(ranFix,"longdep",which=char_rd),
                           Nugget=.get_cP_stuff(ranFix,"Nugget",which=char_rd),
                           d=distmat)  ## so that rho=1 in CauchyCorr
    return(corr_mat)
  }
  #
  calc_moreargs <- function(fixed, char_rd, init.optim, processed, rd, control_dist, NUMAX, LDMAX, ...) {
    rho_size <- max(length(fixed$corrPars[[char_rd]][["rho"]]),length(init.optim$corrPars[[char_rd]][["rho"]]))
    range_info_blob <- .calc_range_info(rho_size, processed, rd, control_dist[[char_rd]]) 
    RHOMAX <- 1.1*30*range_info_blob$nbUnique/range_info_blob$maxrange## matches init$rho in calc_inits() 
    control_dist[[char_rd]]$rho.mapping <- range_info_blob$rho_mapping 
    moreargs_rd <- list(maxrange=range_info_blob$maxrange, RHOMAX=RHOMAX,
                        NUMAX=NUMAX, LDMAX=LDMAX,## not variable...
                        nbUnique=range_info_blob$nbUnique, control.dist=control_dist[[char_rd]]## needed for LUarglist, not calc_inits
    ) 
    return(moreargs_rd) ## with control_dist with possibly modified rho_mapping
  }
  #
  structure(list(corr_family="Cauchy",
                 names_for_post_process_parlist=c("rho","shape","longdep","Nugget"),
                 canonize=canonize,calc_inits=calc_inits,
                 calc_corr_from_dist=calc_corr_from_dist,
                 calc_moreargs=calc_moreargs
                ),
  class="corr_family")
}

corrMatrix <- function(...) {
  structure(list(corr_family="corrMatrix",
                 names_for_post_process_parlist=c(),
                 canonize=function(...) {return(NULL)}, ## NULL OK since elements are copied in [[char_rd]], not [[rd]]
                 calc_inits=function(inits, ...) {return(inits)},
                 calc_moreargs= function(...) {return(NULL)} ## NULL OK since it returns in [[char_rd]], not [[rd]]
  ),
  class="corr_family")
}

corrFamily <- function(corrfamily=NULL, ...) {
  if ( ! is.null(corrfamily)) { # true initialization
    structure(list(corr_family="corrFamily",
                   names_for_post_process_parlist=corrfamily$parnames,
                   canonize=function(corrPars_rd, checkComplete, ...) {list(corrPars_rd=corrPars_rd)}, ## no transfo defined yet
                   calc_inits=function(inits, char_rd, 
                                       optim.scale, # currently ignored (not passed)
                                       init.optim,  
                                       ranFix,  
                                       user.lower,
                                       user.upper,
                                       ...) { # possibly add moreargs_rd
                     inits <- .calc_inits_corrFamily(corrfamily=corrfamily, init=inits$init, char_rd=char_rd, optim.scale="", 
                                                     init.optim=inits$init.optim, init.HLfit=inits$init.HLfit, ranFix=inits$ranFix, 
                                                     user.lower=user.lower,user.upper=user.upper)
                     return(inits)
                   },
                   calc_moreargs= function(...) {return(NULL)} ## NULL OK since it returns in [[char_rd]], not [[rd]]
    ),
    class="corr_family")
  } else { # a kind of declaration
    structure(list(corr_family="corrFamily",
                   names_for_post_process_parlist=function(...) {stop("code missing here (1)")},
                   canonize=function(...) {stop("code missing here (2)")}, ## NULL OK since elements are copied in [[char_rd]], not [[rd]]
                   calc_inits=function(inits, ...) {stop("code missing here (3)")},
                   calc_moreargs= function(...) {return(NULL)} ## NULL OK since it returns in [[char_rd]], not [[rd]]
    ),
    class="corr_family")
  }
}

AR1 <- function(...) {
  canonize <- function(corrPars_rd, checkComplete, ...) {
    if (is.null(corrPars_rd$ARphi) && checkComplete) {
      stop("ARphi missing from ranPars.")
    }
    return(list(corrPars_rd=corrPars_rd))
  }
  #
  calc_inits <- function(inits, user.lower, user.upper, char_rd, ...) {
    inits <- .calc_inits_ARphi(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                               user.lower=user.lower,user.upper=user.upper,char_rd=char_rd)
    return(inits)
  }
  #
  calc_corr_from_dist <- function(ranFix, char_rd, distmat, ...) {
    if (methods::.hasSlot(distmat,"x")) { # "dsCDIST"; or dgCMatrix after subsetting in prediction code  
      corr_mat <- distmat
      corr_mat@x <- .get_cP_stuff(ranFix,"ARphi",which=char_rd)^(corr_mat@x)  
      corr_mat@x[is.na(corr_mat@x)] <- 1
      attr(corr_mat,"dsCDIST") <- NULL
    } else {
      corr_mat <- .get_cP_stuff(ranFix,"ARphi",which=char_rd)^distmat   
      corr_mat[distmat==Inf] <- 0 # instead of the NaN resulting from negative rho...
    }
    corr_mat
  }
  #
  structure(list(corr_family="AR1",
                 names_for_post_process_parlist=c("ARphi"),
                 canonize=canonize, calc_inits=calc_inits,
                 calc_corr_from_dist=calc_corr_from_dist,
                 calc_moreargs= function(...) {return(NULL)} ## NULL OK since it returns in [[char_rd]], not [[rd]]
                ),
            class="corr_family")
}

adjacency <- function(...) {
  canonize <- function(corrPars_rd, cP_type_rd, checkComplete, moreargs_rd, ...) {
    if ( ! is.null(corrPars_rd$trRho)) { ## assuming a single trRho with possibly several elements
      corrPars_rd$rho <- .rhoInv(corrPars_rd$trRho,RHOMAX=moreargs_rd$RHOMAX)  
      corrPars_rd$trRho <- NULL
      cP_type_rd$rho <- cP_type_rd$trRho
      cP_type_rd$trRho <- NULL
    } ## else there may simply be rho rather than trRho (including for adjacency model through optim procedure !)
    return(list(corrPars_rd=corrPars_rd, cP_type_rd=cP_type_rd))
  }
  #
  calc_inits <- function(inits, user.lower, user.upper, moreargs_rd, char_rd, For, ...) {
    inits <- .calc_inits_Auto_rho(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                                  user.lower=user.lower,user.upper=user.upper,
                                  rhorange=moreargs_rd$rhorange,For=For,char_rd=char_rd)
    return(inits)
  }
  #
  calc_moreargs <- function(decomp, verbose, lower, upper, char_rd, ...) {
    if (is.null(decomp$eigrange)) { # this may occur when processing submodel info in mv case
      # then range(NULL) is -Inf Inf and rhorange would be 0 0, and by the code below.
      # -> would go into lower, upper -> would constrain any further attempt to define rhorange. Hence, Instead...
      rhorange <- .Machine$double.xmax*c(-0.1,0.1) # diff(rhorange)  should not overflow
      # and .makeLowUp_stuff_mv() should provide the true bounds
    } else {
      rhorange <- sort(1/range(decomp$eigrange)) ## keeping in mind that the bounds can be <>0
    }
    if(verbose["SEM"])  cat(paste("Feasible rho range: ",paste(signif(rhorange,6),collapse=" -- "),"\n"))
    rhorange[1L] <- max(rhorange[1L],lower$corrPars[[char_rd]][["rho"]])
    rhorange[2L] <- min(rhorange[2L],upper$corrPars[[char_rd]][["rho"]])
    return(list(rhorange=rhorange))
  }
  #
  structure(list(corr_family="adjacency",
                 names_for_post_process_parlist=c("rho"),
                 canonize=canonize, calc_inits=calc_inits,
                 calc_moreargs=calc_moreargs
                ),
            class="corr_family")
}

SAR_WWt <- adjacency ## for .calc_inits and .calc_moreargs: implied by older code where only argument values may differ between adjacency and SAR_WWt 

print.corr_family <- function (x, ...) {
  cat("\nCorrelation family:", x$corr_family, "\n")
  #cat("Link function:", x$link, "\n\n")
  invisible(x)
}

