# processed$hyper_info <- .preprocess_hyper(.) provides a list used by .calc_inits() -> .calc_inits_hyper() and by .expand_hyper()
# .calc_inits_hyper() fills init.optim$hyper <- structure(hyper,map=hyper_info$map,ranges=hyper_info$ranges), 
#                           whic his a nested list of params with attr. The nested list  serves as a template for .makeLowerUpper.
# expand_hyper() is used (notably) by HLCor.obj to fill $corrPars and $trLambda from [ranPars derived from ranefParsVec] and processed$hyper_info


IMRF <- function(...) {
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
    inits <- .calc_inits_IMRF(init=inits$init,init.optim=inits$init.optim,init.HLfit=inits$init.HLfit,ranFix=inits$ranFix,
                             user.lower=user.lower,user.upper=user.upper,optim.scale=optim.scale,
                             moreargs_rd=moreargs_rd,char_rd=char_rd)
    # # Nugget: remains NULL through all computations if NULL in init.optim
    if (is.null(.get_cP_stuff(ranFix,"Nugget",which=char_rd))) { ## new spaMM3.0 code
      inits$init$corrPars[[char_rd]] <- .modify_list(inits$init$corrPars[[char_rd]],
                                                     list(Nugget=.get_cP_stuff(init.optim,"Nugget",which=char_rd)))
    }
    return(inits)
  }
  calc_corr_from_dist <- function(ranFix,char_rd,distmat,...) { stop("This should not be called for IMRF terms") }
  #
  calc_moreargs <- function(KAPPAMAX, IMRF_pars, ...) {
    moreargs_rd <- list(KAPPAMAX=KAPPAMAX) 
    if ( ! is.null(spde_info <- IMRF_pars$model)) {
      Cmat <- spde_info$param.inla$M0
      if (IMRF_pars$SPDE_alpha==2) { # nu=1, alpha=nu+d/2=2
        Gmat <- spde_info$param.inla$M1
        Gzut <- Cholesky(Gmat)
        # Here we may see that G is numerically singular with a small negative eigenvalue:
        eigvals_chol_G <- suppressWarnings(diag(as(Gzut,"sparseMatrix"))) 
        posbad <- which.min(eigvals_chol_G)
        badgval <- eigvals_chol_G[posbad]
        cval <- Cmat@x[posbad]
        # minimum value of kappa such that minimal eigen of Qmat is '1e-12' (?)
        moreargs_rd$minKappa <- sqrt((1e-5*sqrt(cval)-badgval)/cval) # solving (k^2 cval+badgval)^2 / cval =1e-10
      } else if (IMRF_pars$SPDE_alpha==1) { # nu=0, alpha=nu+d/2=1
        Gmat <- spde_info$param.inla$M2
        Gzut <- Cholesky(Gmat)
        # Here we may see that G is numerically singular with a small negative eigenvalue:
        eigvals_chol_G <- suppressWarnings(diag(as(Gzut,"sparseMatrix"))) 
        posbad <- which.min(eigvals_chol_G)
        badgval <- eigvals_chol_G[posbad]
        cval <- Cmat@x[posbad]
        moreargs_rd$minKappa <- sqrt(1e-10-badgval) # solving k^2 cval+badgval=1e-10
      } else {
        # more generally there are two matrices G_1 and G_2 in $M1 and $M2
        x_p <- 0.005
        while (TRUE) {
          chol_p <- Cholesky(.calc_IMRF_Qmat(IMRF_pars,kappa=x_p))
          min_p <- min(suppressWarnings(diag(as(chol_p,"sparseMatrix"))))
          if (min_p>1e-05) break
          x_p <- x_p*10
        } 
        x_m <- x_p/(dlog <- 10)
        while (TRUE) {
          chol_m <- Cholesky(.calc_IMRF_Qmat(IMRF_pars,kappa=x_m))
          min_m <- min(suppressWarnings(diag(as(chol_m,"sparseMatrix"))))
          if (min_m > 0) break
          x_m <- x_m * 9
        }
        # assuming a log-log relationship:
        slope <- (log(min_p)-log(min_m))/log(dlog)
        # targeting a min eigenvalue=1e-10 hence 1e-05 for chol factor
        minKappa <- exp((log(1e-5)-log(min_p))/slope)*x_p 
        while ((min(suppressWarnings(diag(as(Cholesky(.calc_IMRF_Qmat(IMRF_pars,kappa=minKappa)),"sparseMatrix")))))<5e-6) {
          minKappa  <- minKappa*exp(log(2)/slope)
        }
        moreargs_rd$minKappa <- minKappa
      }
    }
    return(moreargs_rd)
  }
  make_new_corr_list <- function(object, old_char_rd, control_dist_rd, geonames, newuniqueGeo, olduniqueGeo, which_mats, make_scaled_dist, new_rd) {
    ## hmmm is that useful ?
  }
  makeLowerUpper <- function() {
    ## hmmm is that useful ?
  }
  # : changes the parent.env of all the member functions.
  structure(list(corr_family="IMRF",
                 names_for_post_process_parlist= c("kappa"),
                 canonize=canonize,
                 calc_inits=calc_inits,
                 calc_corr_from_dist=calc_corr_from_dist,
                 calc_moreargs=calc_moreargs
                 #,make_new_corr_list=make_new_corr_list
  ),
  class="corr_family")
}

.IMRFcrossfactor <- function(xstwm, ystwm, kappa) {
  indices <- expand.grid(seq(xstwm), seq(ystwm))
  nc_I <- nrow(indices)
  locpos0 <- seq(nc_I)
  locpos1 <- seq(xstwm*(ystwm-1L))
  locpos2 <- locpos0[-xstwm*seq(ystwm)]
  B <- sparseMatrix(i=c(locpos0,locpos1,locpos2), j=c(locpos0,locpos1+xstwm,locpos2+1L), 
                    x=c(rep(4+kappa^2,nc_I),rep(-1,length(locpos1)+length(locpos2))), symmetric =TRUE)
  # B <- as(B,"symmetricMatrix") # dsCMatrix ## potentially slow and only useful if speeding a crossprod...
  return(B) ## such that Q= crossprod(B)
  ## or :
  # B <- new("dgCMatrix",i= 0:(nc_I-1L), p=c(0L:(nc_I)), Dim=c(nc_I,nc_I),x=rep(4+kappa^2,nc_I))
  # locpos <- seq(xstwm*(ystwm-1L))
  # B[cbind(locpos,locpos+xstwm)] <- -1
  # locpos <- seq(xstwm*(ystwm-1L))
  # B[cbind(locpos+xstwm, locpos)] <- -1
  # locpos <- seq(nc_I)[-xstwm*seq(ystwm)]
  # B[cbind(locpos, locpos+1L)] <- -1
  # B[cbind(locpos+1L, locpos)] <- -1
  #B <- as(B,"symmetricMatrix") # dsCMatrix
}

.Wendland <- function(d) {
  res <- 0*d
  islowd <- d<1
  lowd <- d[islowd]
  res[islowd] <- (35*lowd^2 + 18*lowd + 3)*(1-lowd)^6 /3 
  return(res)
}

.to_grid_coord <- function(uniqueScal, steplen, origin) {
  for (dt in seq(ncol(uniqueScal))) {
    uniqueScal[,dt] <- (uniqueScal[,dt]-origin[dt])/steplen  
  }
  return(uniqueScal)
}

# optimized 2D version of locatePointinvT
.locate_in_tv <- function(points, ## the points to be located
                          tv, meshloc## mesh$graph$tv, and the mesh vertices coordinates mesh$loc
                          ) {
  pmul <- cbind(-1,diag(nrow =2))
  nr <- nrow(points)
  p2m.t <- rep(NA_integer_, nr)
  p2m.b <- matrix(NA,nrow=nr,ncol=3)
  solve_tvM_table <- vector("list", nrow(tv))
  solve_origins_table <- vector("list", nrow(tv))
  origins <- meshloc[tv[,1],1:2]
  for (tv_it in seq(nrow(tv))) {
    simplex <- meshloc[tv[tv_it,],1:2]
    vM <- pmul %*% simplex # is simplex[-1,]-simplex[1,] for each simplex
    detVM <- vM[4]*vM[1]-vM[2]*vM[3]
    if (detVM==0) { ## which should not occur if the mesh is OK
      solve_tvM_table[[tv_it]] <- ginv(t(vM))
    } else solve_tvM_table[[tv_it]] <- matrix(c(vM[c(4,3,2,1)]*c(1,-1,-1,1))/detVM,ncol=2) #solve(t(vM))
    solve_origins_table[[tv_it]] <- solve_tvM_table[[tv_it]] %*% origins[tv_it,]
  }
  solve_tvMs <- do.call(rbind,solve_tvM_table)
  solve_origins <- unlist(solve_origins_table)
  allWeights <- solve_tvMs %*% t(points) - solve_origins
  seqpos <- seq(nrow(tv))*2
  allmW <- pmin(1-allWeights[seqpos-1L,]-allWeights[seqpos,],allWeights[seqpos-1L,],allWeights[seqpos,]) ## minimum coord in each triangle
  for (pt_it in seq_len(nr)) {
    #point <- points[pt_it,] ## drops...
    #Weights <- solve_tvMs %*% point - solve_origins
    #Weights <- allWeights[,pt_it]
    #mW <- pmin(1-Weights[seqpos-1L]-Weights[seqpos],Weights[seqpos-1L],Weights[seqpos]) ## minimum coord in each triangle
    mW <- allmW[,pt_it]
    which_t <- which(mW>0) # If the point is in a triangle there is one such value
    if ( length(which_t) != 1L ) {
      # either the length is 0 (point seems out of tesselation) = we take the best candidate triangle and apply the next test
      #    or langth>1: numerical precision issue => a point appears to be in two triangles. We again take the best candidate.
      # => which.max() in both cases
      which_t <- which.max(mW)
    } 
    if (mW[which_t] > - 1e-12) { # numerical threshold for out-of-mesh
      p2m.t[pt_it] <- which_t
      pos <- seqpos[which_t]
      Weights <- allWeights[,pt_it]
      p2m.b[pt_it,] <- c(1-Weights[pos-1L]-Weights[pos],Weights[pos-1L],Weights[pos])
    } # else (point is out of tesselation) p2m.t remains NA_integer and p2m.b remains NA
  }
  return(list(p2m.t=p2m.t, # vector (rather than INLA's 1-col matrix)
              p2m.b=p2m.b))
}

.spaMM_spde.make.A <- function(mesh, points) { # replacement for inla.spde.make.A()
  locations <- .locate_in_tv(points=points, #[1:2,]+0.01,
                             tv=mesh$graph$tv,
                             meshloc=mesh$loc[,1:2])
  ti <- locations$p2m.t # vector !
  ii <- which(ok <- (! is.na(ti)))
  # rows of A will remain 0 for points out of the mesh:
  A <- sparseMatrix(dims = c(nrow(points), mesh$n), 
                    i = rep(ii, 3), 
                    j = as.vector(mesh$graph$tv[ti[ii], ]), 
                    x = as.vector(locations$p2m.b[ii, ]))
  A <- drop0(A) 
  #list(t = ti, bary = locations$p2m.b, A = A, ok = ok, A = A) # this is what inla.mesh.project.inla.mesh returns
  return(A)
}

.calc_AMatrix_IMRF <- function(term, # its attr(.,"pars") carries grid parameters
                               data, dist.method, old_AMatrix_rd=NULL, 
                               scale) { ## scale is in steps units so does not depend on the step length in users' coordinates
  pars <- attr(attr(term,"type"),"pars")
  blob <- .get_dist_nested_or_not(term, data=data, distMatrix=NULL, uniqueGeo=NULL, 
                                  dist.method = dist.method, needed=c(uniqueGeo=TRUE),  geo_envir=NULL)
  uniqueScal <- blob$uniqueGeo
  uniqueScal_levels_blob <- .as_factor(txt=paste(colnames(uniqueScal), collapse="+"), 
                                       mf=uniqueScal, type=.spaMM.data$options$uGeo_levels_type) # to match the call in .calc_Zmatrix()
  if ( ! is.null(spde_info <- pars$model)) { # F I X M E does not handle nesting
    # Amatrix <- INLA:::inla.spde.make.A(spde_info$mesh, as.matrix(uniqueScal)) 
    if (inherits(spde_info,"inla.spde2")) {
      Amatrix <- NULL
      if (identical(.spaMM.data$options$INLA_A,TRUE)) { # for debugging, use trace(INLA::inla.spde.make.A) to see the fmesher.args args contained in "mesh", etc.
        Amatrix <- .do_call_wrap("inla.spde.make.A", 
                                 arglist=list(mesh=spde_info$mesh, # or simply mesh=spde_info, as doc'ed
                                              loc=as.matrix(uniqueScal)), 
                                 pack="INLA")
      } # Amatrix still NULL if INLA not available.
      if (is.null(Amatrix)) {
        Amatrix <- .spaMM_spde.make.A(mesh=spde_info$mesh, # FIXME ? save the matrix precomputations from .locate_in_tv() ?
                                           points=as.matrix(uniqueScal))
      } else Amatrix <- drop0(Amatrix) ## there is some noise
      # Amatrix rows are, with the default arguments, ordered as uniqueScal rows
      rownames(Amatrix) <- uniqueScal_levels_blob$factor # not levels(.) which are alphanumerically ordered !  
      return(Amatrix)
    } else stop("Unhandled model class for IMRF")
  }
  # ELSE case or regular grid+Wendland
  if ( ! is.null(old_AMatrix_rd)) {
    uniqueScal <- .to_grid_coord(uniqueScal,steplen= attr(old_AMatrix_rd,"steplen"),
                                 origin=attr(old_AMatrix_rd,"origin")) #IXME would be nicer if ti kept row names.
    grid_arglist <- attr(old_AMatrix_rd,"grid_arglist") ## reconstruct the original grid
    scale <- attr(old_AMatrix_rd,"scale")
  } else {
    nc <- ncol(uniqueScal)
    ranges <- diffs <- grid_arglist <- vector("list",nc) ## typically nc=2
    for (dt in seq(nc)) { 
      ranges[[dt]] <- range(uniqueScal[,dt]) 
      diffs[[dt]] <- diff(ranges[[dt]])
    }
    widedim <- which.max(unlist(diffs))
    steplen <- diffs[[widedim]]/(pars[["nd"]]-1L)
    origin <- numeric(nc)
    for (dt in seq(nc)) {
      if (pars$ce) {
        origin[dt] <- ranges[[dt]][1L] - (diffs[[dt]] %% steplen)/2 ## or using trunc... ## 0 for the wide dimension
      } else origin[dt] <- ranges[[dt]][1L] # ce=FALSE reproduces the LatticeKrig results
      nd_dt <- (diffs[[dt]] %/% steplen) +1L ## central grid node for dim dt
      grid_arglist[[dt]] <- (-pars$m):(nd_dt+pars$m-1L) ## sequence of nodes in grid coordinates for dim dt
    }
    #
    uniqueScal <- .to_grid_coord(uniqueScal, origin=origin, steplen=steplen)
  }
  rownames(uniqueScal) <- uniqueScal_levels_blob$factor # these names must be retained in the output (otherwise visible error in .calc_ZAlist())  
  grid <- expand.grid(grid_arglist) ## INTEGER GRID
  #
  scaldistmat <- proxy::dist(x=uniqueScal,y=grid,method=dist.method) ## coordinates of data in grid units
  if ( ! is.null(old_AMatrix_rd)) {  
    # colnames(scaldistmat) <- colnames(old_AMatrix_rd) # not always true as old_AMatrix_rd cols may have been permuted after original call to .calc_AMatrix_IMRF() 
    return(as(.Wendland(scaldistmat[]/scale),"sparseMatrix")) ## dividing the matrix not most economical but clearer # [] convert crossdist to matrix
  } else {
    colnames(scaldistmat) <- apply(grid,1L,paste0,collapse=":") ## provide colnames(ZA), expected at least by ranef.HLfit 
    return(structure(as(.Wendland(scaldistmat[]/scale),"sparseMatrix"), ## dividing the matrix not most economical but clearer # [] convert crossdist to matrix
                     grid_arglist=grid_arglist, origin=origin, steplen=steplen, scale=scale)) ## attributes for quick reconstruction in prediction
  }  
}

.assign_AMatrices_IMRF <- function(corr_info, Zlist, exp_barlist, data, control_dist,
                                  scale=2.5, ## scale value from Nychka et al 2015 http://dx.doi.org/10.1080/10618600.2014.914946, p. 584
                                  centered=TRUE) { 
  corr_types <- corr_info$corr_types
  isIMRF <- (corr_types=="IMRF")
  if (any(isIMRF, na.rm=TRUE)) {
    if (is.null(corr_info$AMatrices)) corr_info$AMatrices <- list()
    for (rd in which(isIMRF)) {
      corr_info$AMatrices[[as.character(rd)]] <- .calc_AMatrix_IMRF(term=exp_barlist[[rd]], # its attr(.,"pars") carries grid parameters
                                                                    data=data, 
                                                     dist.method=control_dist[[rd]][["dist.method"]], scale=scale)
    }
  }
}

.get_new_AMatrices <- function(object, new_mf_ranef) {
  amatrices <- attr(object$ZAlist,"AMatrices")
  exp_spatial_terms <- attr(object$ZAlist,"exp_spatial_terms")
  isIMRF <- (attr(exp_spatial_terms,"type") == "IMRF")
  for (rd in which(isIMRF)) {
    char_rd <- as.character(rd)
    perm <- attr(amatrices[[char_rd]], "perm") # provided by .assign_geoinfo_and_LMatrices_but_ranCoefs()
    amatrices[[char_rd]] <- .calc_AMatrix_IMRF(term=exp_spatial_terms[[rd]], data=new_mf_ranef, 
                                         dist.method=.get_control_dist(object,char_rd)$dist.method, 
                                         old_AMatrix_rd = amatrices[[char_rd]])
    if ( ! is.null(perm)) amatrices[[char_rd]] <- .subcol_wAttr(amatrices[[char_rd]], j=perm, drop=FALSE)
  }
  return(amatrices)
}

#"multIMRF"


.expand_multIMRF <- function(bar, levels, margin, coarse=10L, norm=TRUE, centered=TRUE, ...) {
  #nodes_per_dim <- (coarse)*2^seq(levels) # +2*arglist$margin
  nodes_per_dim <- (coarse-1L)*2^(seq(levels)-1L)+1L 
  bar <- paste0(bar, ", nd=", nodes_per_dim, ", m=",margin, ", no=", norm, ", l=",seq(levels), ", ce=",centered)
  as.formula(paste("~ ( IMRF(", paste(bar, collapse = ") + IMRF("), ") )"))[[2]]
}


# called by .preprocess_formula():

.expand_multIMRFs <- local({ 
  levels <- NULL
  function (term,init=TRUE) {
    if (init) levels <<- numeric(0)
    if (!is.name(term) && is.language(term)) {
      if (term[[1]] == as.name("(")) {
        term[[2]] <- .expand_multIMRFs(term[[2]],init=FALSE)
      }
      stopifnot(is.call(term))
      term1 <- as.character(term[[1]])
      if (term1 == "multIMRF") {
        args <- paste(paste(names(term),"=",term)[-c(1,2)], collapse=",")
        arglist <- eval(parse(text = paste("list(",args,")")))
        arglist$bar <- paste0(deparse(term[[2]]))
        levels <<- c(levels,arglist$levels)
        return(do.call(".expand_multIMRF", arglist)) # fn call => uses standard R argument matching setting defaults for coarse, norm...
      } else if (term1 == "IMRF") levels <<- c(levels,NA) # NA: placeholder for other IMRF.
      term[[2]] <- .expand_multIMRFs(term[[2]],init=FALSE)
      if (length(term) == 3) 
        term[[3]] <- .expand_multIMRFs(term[[3]],init=FALSE)
    }
    if (init && length(levels)) {
      return(structure(term, hyper_info=list(levels=levels)))
    } else return(term)
  }
})

.calc_inits_hyper <- function(inits, hyper_info, fixed, moreargs ) {
  # the "expanded" corrPars and trLambda have been prefilled by .calc_inits_IMRF() for the reason explained there. 
  # We remove them here if hyper parameters are to be used. 
  if ( ! is.null(hyper_info)) {
    init.optim <- inits$init.optim ## currently mixes the canon.init and the user init, with conflicting info. 
    ##                                The following code rsolves such conflicts
    canon.init <- inits$init
    ranges <- hyper_info$ranges
    hyper <- hyper_info$template
    for (char_hyper_it in names(hyper)) {
      char_rd_range <- as.character(ranges[[char_hyper_it]])
      firstrand <- (char_rd_range)[1L]
      fix_cP <- fixed$hyper[[char_hyper_it]]
      if ( is.null(fix_cP$hy_trK) && is.null(fix_cP$hy_kap) ) { # if NOT fixed
        hy_kap <- init.optim$hyper[[char_hyper_it]]$hy_kap # user init
        if (is.null(hy_kap)) {
          hy_trK <- init.optim$corrPars[[firstrand]]$trKappa
          hy_kap <- .kappaInv(hy_trK, moreargs[[firstrand]]$KAPPAMAX)
        } else hy_trK <- .kappaFn(hy_kap, moreargs[[firstrand]]$KAPPAMAX)
        canon.init$hyper[[char_hyper_it]]$kappa <- hy_kap
        hyper[[char_hyper_it]]$hy_trK <- hy_trK
        for (char_rd in char_rd_range) { 
          init.optim$corrPars[[char_rd]]$trKappa <- NULL 
        }
      } 
      for (char_rd in char_rd_range) { ## in all cases because either tracking the init.optim removal, or the fact it is fixed
        canon.init$corrPars[[char_rd]]$kappa <- NULL 
      }
      if (is.null(fix_cP$hy_trL) && is.null(fix_cP$hy_lam)) {
        hy_lam <- init.optim$hyper[[char_hyper_it]]$hy_lam # user init
        if (is.null(hy_lam)) hy_lam <- init.optim$lambda[firstrand]
        canon.init$hyper[[char_hyper_it]]$hy_lam <- hy_lam
        hyper[[char_hyper_it]]$hy_trL <- .dispFn(hy_lam)
      } 
      ## In all cases because either tracking the init.optim$hyper, or the fact it is fixed:
      init.optim$lambda[char_rd_range] <- NA # tag for removal below
      canon.init$lambda[char_rd_range] <- NA # tag for removal below
    } # but for fiXed values, expanded values are still in init.optim hence remove_from_parlist() below
    if (length(unlist(hyper))) init.optim$hyper <- hyper # if some hyper_params are not fixed
    # 'fixed' must have expanded values for use in remove_from_parlist and for all further purposes
    # 'init.optim' must still have expanded values for use in remove_from_parlist
    init.optim <- remove_from_parlist(init.optim, inits$ranFix) # removes fixed 'hyper', for example
    # if hyper still there, we add attributes (which would have been lost by remove_from_parlist())
    if ( ! is.null(init.optim$hyper)) init.optim$hyper <- structure(init.optim$hyper,
                                                                    map=hyper_info$map,ranges=hyper_info$ranges )
    # Now it is important to remove NA's from expanded values:
    init.optim$lambda <- init.optim$lambda[!is.na(init.optim$lambda)]
    inits[["init.optim"]] <- init.optim
    # And parallel removal in canon.init
    canon.init$lambda <- canon.init$lambda[!is.na(canon.init$lambda)] 
    inits[["init"]] <- canon.init ## as this is used as a template for lower$trLambda...
  }
  return(inits)
}

.preprocess_hyper <- function(processed) { ## needs $predictor, $ZAlist
  level_info <- attr(processed$predictor,"hyper_info")$levels ## a vector of number of levels, provided by .expand_multiMRFs()
  hy_levels <- na.omit(level_info) 
  NA_pos <- attr(hy_levels,"na.action") # position of non-multIMRF IMRFs
  if (len <- length(hy_levels)) { # if multIMRF's...
    template <- ranges <- structure(vector("list",len),names=paste(seq(len))) # could we use ranges as template ?
    nrand <- length(processed$ZAlist)
    hyper_map <- rep(NA_integer_, nrand) ## maps rd to elements of hyper
    level_info[NA_pos] <- 1L
    IMRF_ranges <- c(0,cumsum(level_info))
    which_IMRF <- which(attr(processed$ZAlist,"exp_ranef_types") =="IMRF")
    hy_pos <- seq_along(level_info)
    if ( ! is.null(NA_pos)) hy_pos <- hy_pos[- NA_pos]
    for (it in seq_along(hy_pos)) { 
      hy_it <- hy_pos[it]
      rd_range <- which_IMRF[(IMRF_ranges[hy_it]+1L):(IMRF_ranges[hy_it+1L])]
      hyper_map[rd_range] <- it 
      ranges[[as.character(it)]] <- rd_range ## inverse map
    }
    names(hyper_map) <- paste(seq_len(length(hyper_map)))
    blocs <- list()
    rd <- blocit <- 1L
    # Matrix::bdiag ignores names
    rowN <- colN <- list()
    while (rd < nrand+1L) {
      if (is.na(hyper_it <- hyper_map[rd])) {
        char_rd <- as.character(rd)
        blocs[[blocit]] <- 1
        rowN[[blocit]] <- colN[[blocit]] <- char_rd 
        rd <- rd+1L
      } else {
        bllen <- length(ranges[[hyper_it]])
        blocs[[blocit]] <- rep(1,bllen)
        seq_rd <- rd-1L + seq(bllen)
        rowN[[blocit]] <- paste(seq_rd)
        colN[[blocit]] <- paste0(seq_rd[1L],":",seq_rd[bllen]) 
        rd <- rd+bllen
      }
      blocit <- blocit+1L
    }
    summingMat <- as.matrix(Matrix::bdiag(blocs))
    rownames(summingMat) <- unlist(rowN)
    colnames(summingMat) <- unlist(colN)
    return(list(map=hyper_map,ranges=ranges, template= template, summingMat=summingMat) )
  }
}

.expand_hyper <- function(ranPars, hyper_info, moreargs) { # called to build inits, in HLCor.obj, and for refit
  # DO NOT assume any type attribute here => no attributes management here
  if (!is.null(hyper_info)) {
    ranges <- hyper_info$ranges ## e.g list("1"=1:3)
    trLam <- hyper_info$map # template
    trLam[] <- NA
    trLam[seq_along(ranPars$trLambda)] <- ranPars$trLambda ## initialization with value for other ranefs
    ## ranPars$trLambda typically NULL from user input (in $lambda), but no longer so after ranPars_in_refit <- ... 
    ## FIXME ideally ranPars$trLambda has char_rd names... but I'm not sure it's always the case.
    for (rg_it in seq_along(ranges)) { # then: 1
      hyper_el <- ranPars$hyper[[rg_it]]
      char_rd_range <- as.character(ranges[[rg_it]]) 
      if (is.null(hy_lam <- hyper_el$hy_lam)) {
        hy_trL <- hyper_el$hy_trL
        if ( ! is.null(hy_trL)) hy_lam <- .dispInv(hyper_el$hy_trL) # else hy_lam remains NULL
      } 
      if (! is.null(hy_lam)) {
        lam_fac <- c(1, 2^(-2*(seq_len(length(char_rd_range)-1L)))) # 1, 1/4, 1/16
        lambdas <- hy_lam * lam_fac/sum(lam_fac) 
        trLam[char_rd_range] <- .dispFn(lambdas) ## over-writes temporary NaN ## names needed for remove_from_parlist() later
      }
      if (is.null(hy_trK <- hyper_el$hy_trK)) {
        hy_kap <- hyper_el$hy_kap
        if (!is.null(hy_kap)) {
          KAPPAMAX <- moreargs[[char_rd_range[1L]]]$KAPPAMAX
          if ( hy_kap > KAPPAMAX - 5e-7) {
            warning("Fixed value of IMRF's kappa appears above allowed upper value. I reduce it to the allowed maximum. ")
            hy_kap <- KAPPAMAX - 1e-6
          }
          hy_trK <- .kappaFn(hy_kap,KAPPAMAX=KAPPAMAX)
        }
      }
      if (! is.null(hy_trK)) {
        # this syntax does not work:
        #ranPars$corrPars[char_rd_range]$trKappa <- hyper_el$hy_trK ## same SAR parameter for all levels
        for (char_rd in char_rd_range) {
          ranPars$corrPars[[char_rd]]$trKappa <- hy_trK ## same SAR parameter for all levels
        }
      } 
    }
    ranPars$trLambda <- trLam # [ ! is.na(trLam)] ... removing NA's would be a problem for .modify_list()...
  } 
  return(ranPars) ## with trLambda, $corrPars modified with TRANSFORMED values
}


.strip_IMRF_args <- function (term) { ## compare to lme4:::nobars_ ; 'term is formula or any term, recursively
  if (!("|" %in% all.names(term))) 
    return(term)
  if (term[[1]] == as.name("IMRF") || term[[1]] == as.name("multIMRF")) {
    return(term[c(1,2)]) # (lhs|rhs) without the pars 
  }
  if (length(term) == 2) {
    nb <- .strip_IMRF_args(term[[2]])
    if (is.null(nb)) 
      return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- .strip_IMRF_args(term[[2]])
  nb3 <- .strip_IMRF_args(term[[3]])
  if (is.null(nb2)) 
    return(nb3)
  if (is.null(nb3)) 
    return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

str.inla.spde2 <- function(object, verbose=("package:INLA" %in% search()), ...) {
  if (verbose) {
    class(object) <- setdiff(class(object), "inla.spde2") # to call default str() methods
    str(object, ...)
  } else {
    cat("<inla.spde2 object>; load INLA, or use str(.,verbose=TRUE) if you want to str() this object.\n")
    invisible()
  }
}
