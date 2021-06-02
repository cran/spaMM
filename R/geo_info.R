.get_geo_info <- function(processed, which_ranef, which="", dist_method_rd=processed$control_dist[[which_ranef]]$dist.method) {
  geo_envir <- processed$geo_info[[which_ranef]] ## should be an environment, possibly empty. Check:
  if (is.null(geo_envir)) stop("The required environment was not created during preprocessing.")
  needed <- c("distMatrix"=("distMatrix" %in% which && is.null(geo_envir$distMatrix)),
              "uniqueGeo"=("uniqueGeo" %in% which && is.null(geo_envir$uniqueGeo)),
              "nbUnique"=("nbUnique" %in% which && is.null(geo_envir$nbUnique)),
              "notSameGrp"=("notSameGrp" %in% which && is.null(geo_envir$notSameGrp)) )
  if (any(needed)) {
    blob <- .get_dist_nested_or_not(spatial_term=attr(processed$ZAlist,"exp_spatial_terms")[[which_ranef]], 
                                data=processed$data, distMatrix=geo_envir$distMatrix, 
                                uniqueGeo=geo_envir$uniqueGeo, ## precomputed only for AR1
                                dist.method = dist_method_rd, needed=needed,
                                geo_envir=geo_envir)
    ## Avoid overwriting preexisting ones with possible NULL's:
    if (is.null(geo_envir$distMatrix)) geo_envir$distMatrix <- blob$distMatrix 
    if (is.null(geo_envir$uniqueGeo)) geo_envir$uniqueGeo <- blob$uniqueGeo  
    if (is.null(geo_envir$nbUnique)) geo_envir$nbUnique <- blob$nbUnique 
    if (is.null(geo_envir$notSameGrp)) geo_envir$notSameGrp <- blob$notSameGrp 
    if (is.null(geo_envir$coordinates)) geo_envir$coordinates <- blob$coordinates 
  }
  return(geo_envir) ## some elements may still be NULL
}

.inla_spde2_theta2phi0 <- function (spde, theta) {
  if (spde$n.theta > 0) {
    return(exp(spde$param.inla$B0[, 1, drop = TRUE] + 
                 spde$param.inla$B0[,-1, drop = FALSE] %*% theta))
  } else return(exp(spde$param.inla$B0[, 1, drop = TRUE]))
}

.inla_spde2_theta2phi1 <- function (spde, theta) {
  if (spde$n.theta > 0) {
    return(exp(spde$param.inla$B1[, 1, drop = TRUE] + 
                 spde$param.inla$B1[, -1, drop = FALSE] %*% theta))
  } else return(exp(spde$param.inla$B1[, 1, drop = TRUE]))
}

.inla_spde2_theta2phi2 <- function (spde, theta) {
  if (spde$n.theta > 0) {
    phi = ((spde$param.inla$B2[, 1, drop = TRUE] + 
              spde$param.inla$B2[, -1, drop = FALSE] %*% theta))
  } else {
    phi = ((spde$param.inla$B2[, 1, drop = TRUE]))
  }
  if (spde$param.inla$transform == "identity") {
    return(phi)
  } else if (spde$param.inla$transform == "logit") {
    return(cos(pi/(1 + exp(-phi))))
  } else if (spde$param.inla$transform == "log") {
    return(2 * exp(phi) - 1)
  } else {
    warning(paste("Unknown link function '", spde$param.inla$transform, 
                  "' phi2.  Using identity link instead.", sep = ""))
    return(phi)
  }
}

.inla_spde_precision_inla_spde2 <- function (spde, theta = NULL, 
                                            phi0 = .inla_spde2_theta2phi0(spde, theta), 
                                            phi1 = .inla_spde2_theta2phi1(spde, theta), 
                                            phi2 = .inla_spde2_theta2phi2(spde, theta), ...) 
{
  if (spde$f$model != "spde2") {
    stop("spaMM only supports some internal inla models 'spde2'")
  }
  # matches formula for Q p.5 of https://www.jstatsoft.org/article/view/v063i19
  D0 = Diagonal(spde$n.spde, phi0) # T
  D1 = Diagonal(spde$n.spde, phi1) # K^2
  D12 = Diagonal(spde$n.spde, phi1 * phi2) # K^2 ....
  # for $MO = C, $M1 = G_1, $M2 = G_2
  Q = (D0 %*% (D1 %*% spde$param.inla$M0 %*% D1 + D12 %*% spde$param.inla$M1 + 
                 t(spde$param.inla$M1) %*% D12 + spde$param.inla$M2) %*% D0)
  return(Q)
}

.calc_IMRF_Qmat <- function(pars, grid_arglist, kappa, test=FALSE) {
  if (is.null(spde_info <- pars$model)) { 
    crossfac_Q <- .IMRFcrossfactor(xstwm=length(grid_arglist[[1]]), ystwm=length(grid_arglist[[2]]),kappa=kappa)
    sparse_Qmat <- crossprod(crossfac_Q) # same relation as for t_chol_Q ## automatically dsCMatrix
  } else {
    if (inherits(spde_info,"inla.spde2")) { # F I X M E only for d=2 
      sparse_Qmat <- NULL
      # spacial cases for sparse_Qmat:
      # test for # inla.spde.matern with default parameters B.tau and B.kappa
      theta_system <- pars$model$param.inla$BLC[c("tau.1","kappa.1"),] 
      is_default_spde2_matern <- (diff(range((theta_system - matrix(c(0,0,1,0,0,1),ncol=3))))==0)
      if (test || is_default_spde2_matern) { 
        if (pars$SPDE_alpha==2) { # nu=1, alpha=nu+d/2=2
          Cmat <- spde_info$param.inla$M0
          Gmat <- spde_info$param.inla$M1
          Kmat <- kappa^2 * Cmat +Gmat
          if (isDiagonal(Cmat)) {
            # solvesqC <- Cmat
            # solvesqC@x <- 1/sqrt(solvesqC@x)
            # crossfac_Q <- Kmat %*% solvesqC
            tcrossfac_Q <- .Matrix_times_Dvec(Kmat, 1/sqrt(Cmat@x))
            sparse_Qmat <- tcrossprod(tcrossfac_Q) # same relation as for t_chol_Q ## automatically dsCMatrix
          } else {stop("Cmat does not appear to be diagonal")}
        } else if (pars$SPDE_alpha==1) { # nu=0, alpha=nu+d/2=1
          Cmat <- spde_info$param.inla$M0
          Gmat <- spde_info$param.inla$M2
          sparse_Qmat <- as(kappa^2 * Cmat +Gmat,"symmetricMatrix")
        }      
      }
      if (is.null(sparse_Qmat)) { # that should be the general method
        # 'theta_system' represents the system for log(tau) and log(kappa) given p.5 of https://www.jstatsoft.org/article/view/v063i19.
        # The following code inverts this system so that .inla_spde_precision_inla_spde2() always returns the same sparse_Qmat
        # irrespective of the value of theta_system; and thus it always gives the same sparse_Qmat as a call of 
        # inla.spde2.matrix() with default B.tau and B.kappa, which we know to match well the Matern() results.
        theta <- solve(theta_system[,2:3],c(0,log(kappa))-theta_system[,1]) # theta reduces to c(0,log(kappa)) if (is_default_spde2_matern)
        sparse_Qmat <- .inla_spde_precision_inla_spde2(pars$model, theta=theta)
        sparse_Qmat <- as(sparse_Qmat,"symmetricMatrix")
      }
    } else stop("Unhandled model class for IMRF")
  }
  return(sparse_Qmat)
}

.ZA_update <- function(rd, Q_CHMfactor, processed, Amat) {
  tPmat <- t(as(Q_CHMfactor,"pMatrix"))
  if (any(tPmat@perm!=seq_along(tPmat@perm))) {
    levelnames <- colnames(processed$ZAlist[[rd]])
    RRsP <- sort.list(tPmat@perm) 
    colnames(tPmat) <- levelnames[RRsP]
    if (is.null(Amat)) {
      rownames(tPmat) <- colnames(processed$ZAlist[[rd]])
      # : when there is an A matrix, .calc_normalized_ZAlist() checks its names 
      attr(processed$ZAlist,"AMatrices")[[as.character(rd)]] <- structure(tPmat, permuted_Q=TRUE)
    } else {
      Amat <- .subcol_wAttr(Amat,j=RRsP, drop=FALSE) # Amat %*% tPmat
      attr(Amat,"perm") <- RRsP # used by .get_new_AMatrices()
      attr(processed$ZAlist,"AMatrices")[[as.character(rd)]] <- structure(Amat, permuted_Q=TRUE)
    }
    processed$ZAlist[[rd]] <- processed$ZAlist[[rd]] %*% tPmat
    .assign_ZAfix(processed)
  } # else ignore identity tPmat
}

.assign_Lunique_from_Q_CHM <- function(processed, rd, Q_CHMfactor, corr_type, spatial_term, type) {
  ## fitme sparse_precision has an incomplete symSVD=> corr matrix not computed, 
  ##    and try(mat_sqrt(symSVD=symSVD)) fails. Instead use code always valid:
  use_ZA_L <- .spaMM.data$options$use_ZA_L
  if (
    (
      identical(use_ZA_L,TRUE) ||
      (is.null(use_ZA_L) && processed$augZXy_cond) # default is to use_ZA_L for augZXy_cond
    )
    && ! processed$AUGI0_ZX$vec_normIMRF[rd]) { # solve may be OK when IMRF is normalized (otherwise more code needed in .normalize_IMRF) 
    Lunique <- Q_CHMfactor # representation of Lunique, not Lunique itself
    # note that this has an attribute "type" !  
  } else {
    Lunique <- solve(Q_CHMfactor,system="Lt") #  L_Q^{-\top}=LMatrix_correlation 
    ## Lunique is a dtCMatrix; it is for single correlated effect, and still used in HLfit_body; 
    #  Keeping it as dtCMatrix might be faster for some operations (but not .tcrossprod) (F I X M E test)
    ## whether it is used or not in MME_method (sXaug_...), a lot of other code still expects it
    attr(Lunique,"Q_CHMfactor") <- Q_CHMfactor 
    attr(Lunique, "type") <- type
  } 
  attr(Lunique,"corr.model") <- corr_type
  attr(Lunique,"ranefs") <- paste(c(spatial_term)) ## essentiel pour la construction de ZAL! ## paste(c()) handles very long RHS
  processed$AUGI0_ZX$envir$LMatrices[[rd]] <- Lunique
  attr(processed$AUGI0_ZX$envir$LMatrices,"is_given_by")[rd] <- "AUGI0_ZX$envir" 
}

# replaces in HLCor_body() the direct call to .init_AUGI0_ZX_envir_spprec_info() from other _body() functions
.assign_geoinfo_and_LMatrices_but_ranCoefs <- function(processed, corr_types, spatial_terms, 
                                                       ranPars, control.dist, argsfordesignL) {
  # * assigns geo_envir <- .get_geo_info(...)
  # * modifies processed$AUGI0_ZX$envir by .init_AUGI0_ZX_envir_spprec_info(...) 
  # * computes processed$AUGI0_ZX$envir$LMatrices except for ranCoefs (the latter being filled in HLfit_body by .update_AUGI0_ZX_envir_ranCoefs_info())
  if (processed$is_spprec) .init_AUGI0_ZX_envir_spprec_info(processed) 
  AUGI0_ZX_envir <- processed$AUGI0_ZX$envir
  for (rd in seq_along(corr_types)) {
    corr_type <- corr_types[rd]
    symSVD <- NULL ## reinitialize as rd is tested within the loop
    msd.arglist <- NULL
    if (!is.na(corr_type)) {
      char_rd <- as.character(rd)
      rho <- .get_cP_stuff(ranPars,"rho",which=char_rd)
      spatial_term <- spatial_terms[[rd]] 
      adj_rho_is_inner_estimated <- (corr_type=="adjacency"
                                     && is.null(rho) ) ## can occur in direct call of HLCor 
      #
      ###
      # Provide cov_info_mat, sparse_Qmat (except for inner-estimated adjacency rho)
      if (corr_type == "adjacency") {
        adjMatrix <- processed$corr_info$adjMatrices[[rd]] # dsC whether spprec or not
        symSVD <- attr(adjMatrix,"symSVD")
        if ( ! adj_rho_is_inner_estimated) { # typical fitme() call
          sparse_Qmat <- - rho * adjMatrix
          sparse_Qmat <- .dsCsum(sparse_Qmat, attr(adjMatrix,"dsCdiag")) #diag(sparse_Qmat) <- diag(sparse_Qmat)+1 # a maybe-not-yet-perfect solution to the poor perf of 'diag<-' which is not doc'ed for dsCMatrix
          if ( ! processed$is_spprec && is.null(symSVD)) {
            if (ncol(sparse_Qmat>200L)) { # a bit of an ad hoc guess based on ohio...
              cov_info_mat <- as.matrix(solve(sparse_Qmat))
            } else cov_info_mat <- solve(as.matrix(sparse_Qmat))
          } 
          AUGI0_ZX_envir$updateable[rd]=(rho!=0)
        } else { # inner estimation of adjacency rho
          if (is.null(symSVD)) {
            ## Direct call of HLCor (SEM or not), 
            ## I also wrote that this could occur in fitme/corrHLfit if(list(processed)) "bc symSVD added to proc1 (i.e. not to all proc's)"
            ##   but there is not good tests of this case...
            symSVD <- .provide_AR_factorization_info(
              adjMatrix, 
              sparse_precision=processed$is_spprec, # this must be FALSE for inner estimation of adjacency rho (and there is bug-catching code for this)
              corr.model=corr_type)
            attr(processed$corr_info$adjMatrices[[rd]],"symSVD") <- symSVD
            #
            # typically an *initial* value:
            rho <- attr(ranPars,"init.HLfit")$corrPars[[char_rd]]$rho # a bit unclear, but we will compute Lmatrix from this 'initial value' below
            if (is.null(rho)) stop("'rho' missing. Contact the maintainer") # protection to avoid programming errors resulting in segfaults.
          } 
          # symSVD may be modified below
        }
      } else if (corr_type =="IMRF") {
        kappa <- .get_cP_stuff(ranPars,"kappa",which=char_rd)
        pars <- attr(attr(spatial_term,"type"),"pars")
        sparse_Qmat <- .calc_IMRF_Qmat(pars, 
                                       grid_arglist=attr(attr(processed$ZAlist,"AMatrices")[[char_rd]],"grid_arglist"), # promise for spde case
                                       kappa)
        if ( ! processed$is_spprec) cov_info_mat <- chol2inv(chol(sparse_Qmat)) # solve(sparse_Qmat) 
      } else if (corr_type =="SAR_WWt") { 
        adjMatrix <- processed$corr_info$adjMatrices[[rd]]
        sparse_Qmat <- - rho * adjMatrix
        sparse_Qmat <- .dsCsum(sparse_Qmat, attr(adjMatrix,"dsCdiag")) #diag(sparse_Qmat) <- diag(sparse_Qmat)+1 # a maybe-not-yet-perfect solution to the poor perf of 'diag<-' which is not doc'ed for dsCMatrix
        UDU. <- attr(adjMatrix,"UDU.")
        if (is.null(UDU.)) {
          # cf comments on symSVD above
          cov_info_mat <- as.matrix(solve(sparse_Qmat)) 
        } else {
          cov_info_mat <- .ZWZt(UDU.$u, 1/(1-rho*UDU.$d)) # UDU.$u %*% sweep(UDU.$u.,MARGIN=1,1/(1-rho*UDU.$d),`*`) 
        }
        cov_info_mat <- .tcrossprodCpp(cov_info_mat,NULL)
      }  else if (corr_type=="AR1" && ! processed$is_spprec) {
        geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("distMatrix"), 
                                   dist_method_rd=control.dist[[char_rd]]$dist.method)
        cov_info_mat <- .get_cP_stuff(ranPars,"ARphi",which=char_rd)^(geo_envir$distMatrix)  
        cov_info_mat[geo_envir$distMatrix==Inf] <- 0 ## should not be necess, but is.
      } else  if (corr_type %in% c("Matern","Cauchy")) {
        control_dist_rd <- control.dist[[char_rd]]
        msd.arglist <- list(rho = rho)
        msd.arglist$`dist.method` <- dist_method_rd <- control_dist_rd$`dist.method` ## may be NULL
        # .get_geo_info() is *the* code that accounts for grouping but never uses rho (hence never returns a scaled distMatrix)
        # So for non-trivial rho we cannot directly use the distMatrix, and make_scaled_dist() won't have access to grouping info in distMatrix.
        # So for non-trivial rho plus grouping we need both the uniqueGeo and the distMatrix
        if (length(rho)>1L) {
          txt <- paste(c(spatial_term[[2]][[3]])) ## the RHS of the ( . | . ) # c() to handle very long RHS
          if (length(grep("%in%",txt))) {
            needed <- c("uniqueGeo", "notSameGrp")
          } else needed <- "uniqueGeo"
          geo_envir <- .get_geo_info(processed, which_ranef=rd, which=needed, 
                                     dist_method_rd=dist_method_rd)
          coord_within <- .extract_check_coords_within(spatial_term=spatial_term)
          msd.arglist$uniqueGeo <- geo_envir$uniqueGeo[,coord_within,drop=FALSE]
          if ( ! is.null(dist_method_rd) && dist_method_rd %in% c("Earth","EarthChord") ) {
            if (grepl("lat", coord_within[1])) warning("Hmm... the first coordinate should be longitude, but seems to be latitude.")
            if (max(abs(msd.arglist$uniqueGeo[,1])-180>1e-14)) warning("Hmm... max(abs(longitude)) should be <= 180, but is not.")
            if (max(abs(msd.arglist$uniqueGeo[,1])-90>1e-14)) warning("Hmm... max(abs(latitude)) should be <= 90, but is not.")
          }
          msd.arglist$`rho.mapping` <- control_dist_rd$`rho.mapping` ## may be NULL
          cov_info_mat <- do.call("make_scaled_dist",msd.arglist)
          #
          if (length(grep("%in%",txt))) { # hmmmf
            cov_info_mat <- as.matrix(cov_info_mat)
            cov_info_mat[geo_envir$notSameGrp] <- Inf
            cov_info_mat <- proxy::as.simil(cov_info_mat)
          }
        } else {
          geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("distMatrix"), 
                                     dist_method_rd=dist_method_rd)
          msd.arglist$distMatrix <- geo_envir$distMatrix   
          cov_info_mat <- do.call("make_scaled_dist",msd.arglist)
        }
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
      } else if (corr_type== "corrMatrix") { ## could cover all specifications of a constant 'cov' Matrix without correlation parameter
        if ( ! processed$is_spprec) cov_info_mat <- processed$corr_info$cov_info_mats[[rd]] ## correlation or precision...
      } 
      ###
      #
      ## Provide Lunique if not already available (and optionally additional stuff)
      if (processed$is_spprec) {
        .init_AUGI0_ZX_envir_spprec_info(processed)   
        if (corr_type=="AR1") { 
          types <- AUGI0_ZX_envir$finertypes
          ARphi <- .get_cP_stuff(ranPars,"ARphi",which=char_rd)
          ## the dim of Q *HAS* to match that of ZA
          AR1_block_n_u_h_s <- attr(processed$ZAlist[[rd]],"AR1_block_n_u_h_s") # block size[S for single AR1 term if nested]
          ilist <- jlist <- xlist <- vector("list",length(AR1_block_n_u_h_s))
          cum_ncol <- 0L
          for (bloc in seq_along(AR1_block_n_u_h_s)) {
            triplets <- .calc_AR1_t_chol_Q_block(AR1_block_n_u_h_s[[bloc]], ARphi=ARphi) ## L block(s) of given size(s)
            ilist[[bloc]] <- cum_ncol + triplets$i
            jlist[[bloc]] <- cum_ncol + triplets$j 
            xlist[[bloc]] <- triplets$x
            cum_ncol <- cum_ncol + AR1_block_n_u_h_s[[bloc]]
          }
          t_chol_Q <- sparseMatrix(dims = c(cum_ncol,cum_ncol), i=unlist(ilist), j=unlist(jlist),
                                   x = unlist(xlist), triangular=TRUE)
          AUGI0_ZX_envir$precisionFactorList[[rd]] <- list(chol_Q=t(t_chol_Q), # Linv
                                                           Qmat=crossprod(t_chol_Q))
          AUGI0_ZX_envir$updateable[rd]=(ARphi!=0) # to be specifically used for G_CHMfactor
          Lunique <- solve(t_chol_Q)
          attr(Lunique, "type") <- "from_AR1_specific_code"
          attr(Lunique,"corr.model") <- corr_type
          attr(Lunique,"ranefs") <- paste(c(spatial_term)) ## essentiel pour la construction de ZAL! ## paste(c()) handles very long RHS
          AUGI0_ZX_envir$LMatrices[[rd]] <- Lunique
          attr(AUGI0_ZX_envir$LMatrices,"is_given_by")[rd] <- "AUGI0_ZX$envir" 
        } else if (corr_type=="corrMatrix") {
          #.init_AUGI0_ZX_envir_spprec_info() has done the work                #   sparse_Qmat <- drop0(cov_info_mat[["matrix"]])
        } else { ## General code that should be correct for AR1 and corrMatrix too
          ## (1) Provide sparse_Qmat for all parameteric correlation models 
          if (corr_type %in% c("Matern","Cauchy")) {
            ## at this point cov_info_mat is a dist object !
            cov_info_mat <- proxy::as.matrix(cov_info_mat, diag=1)
            ## If the order of columns of ZAlist[[it]] had been permuted by .calc_Zmatrix(), we ould need something like:
            # ZAnames <- colnames(processed$ZAlist[[rd]])
            # cov_info_mat <- cov_info_mat[ZAnames,ZAnames]
            # but this would require controlling the names of Z to be those of cov_info_mat 
            # (spMMFactorList_locfn() does not care for that but .preprocess() does something similar for the Z of precision matrices)
            sparse_Qmat <- as_precision(cov_info_mat)$matrix
          } else if (corr_type=="adjacency" 
                             && ! adj_rho_is_inner_estimated) {
            ## Implies call from fitme_body with outer rho estim.
            # sparse_Qmat already computed before in this fn
            ## Cholesky gives proper LL' (think LDL')  while chol() gives L'L...
          } else if (corr_type== "IMRF") {
            # Remember that we need dtCMatrix'es 'chol_Q' so that bdiag() gives a dtCMatrix
            # Hence use next general code to produce precisionFactorList[[rd]]
          } else stop("Some error occurred (inner estimation of adjacency rho with requested sparse precision ?)")     
          ## solve(sparse_Qmat) gives the correlation matrix
          ## (2) Builds from sparse_Qmat, providing Provides precisionFactorList[[rd]] as expected by .reformat_Qmat_info()
          if (is.null(template <- AUGI0_ZX_envir$precisionFactorList[[rd]]$template)) { ## occurs if $update_CHM is FALSE, OR first comput. of CHM, OR not yet all $updateable[rd]
            get_template <- .spaMM.data$options$update_CHM && AUGI0_ZX_envir$updateable[rd]
            if (get_template) {
              ###  ONE-TIME CODE
              #  I once had a problem with perm_Q=TRUE in test-predVar-Matern-corrMatrix -> predict(f2,.) so make sure to check this when changing the code
              Amat <- attr(processed$ZAlist,"AMatrices")[[as.character(rd)]]
              if (is.null(perm_Q <- .spaMM.data$options$perm_Q)) {
                perm_Q <- ( ! is.null(Amat) || corr_type=="adjacency") # 1st condition => to get IMRFs notably
                # remaining cases are Matern,Cauchy where any permuted precision matrix is presumably full without any useful pattern of zeros 
                ## __F I X M E__ also depend on expected sparsity ? on ZA being identity?
              } 
              Q_CHMfactor <- Cholesky(sparse_Qmat,LDL=FALSE,perm=perm_Q) 
              AUGI0_ZX_envir$precisionFactorList[[rd]]$template <- Q_CHMfactor
              if (perm_Q) .ZA_update(rd, Q_CHMfactor, processed, Amat) # else ZA is not permuted, there is no permuted_Q attribute, so the G matrix will be constructed from unpermuted ZA and unpermuted sparse_Qmat
              ### END OF ONE-TIME CODE
            } else Q_CHMfactor <- Cholesky(sparse_Qmat,LDL=FALSE,perm=FALSE) 
            # perm=TRUE without saving the result in a template should be a problem: 
            #  as(Q_CHMfactor, "sparseMatrix") should be permuted (with possibly different permutations over calls), but then ZA should be permuted
            # And .ZA_update() is valid for mermuting ZA only once but not repeteably (original unpermutate ZA not kept)
          } else { # template for updating already exists
            Q_CHMfactor <- Matrix::.updateCHMfactor(template, parent=sparse_Qmat, mult=0) 
          }
          #
          permuted_Q <- attr(attr(processed$ZAlist,"AMatrices")[[as.character(rd)]],"permuted_Q") # clumsy but need to get info from one-time code
          # $Qmat <- sparse_Qmat will be used together with ZA independently from the CHM to construct the Gmat
          # If we use permuted Chol AND a template, ZA has been permuted (once). Then we must permute each new sparse_Qmat, by 
          if (identical(permuted_Q,TRUE)) sparse_Qmat <- tcrossprod(as(Q_CHMfactor,"sparseMatrix")) 
          AUGI0_ZX_envir$precisionFactorList[[rd]]$Qmat <- sparse_Qmat # should by *dsC*  
          AUGI0_ZX_envir$precisionFactorList[[rd]]$chol_Q <- as(Q_CHMfactor, "sparseMatrix") # Linv
          #
          .assign_Lunique_from_Q_CHM(processed, rd, Q_CHMfactor, corr_type, spatial_term, type="from_Q_CHMfactor")  ## assigns to AUGI0_ZX$envir
        }
      } else { ## sparse or dense CORRELATION algorithms
        # __F I X M E__ it would be nice to extract a function similar to .assign_Lunique_from_Q_CHM() but there are many arguments used here:
        # this is called through HLCor, hence there must be a dense correlation structure
        ## densePrecision cases with outer rho estimation, and some residual inner rho estimation cases
        if ( ! is.null(symSVD)) {
          symSVD$d <- 1/(1-rho*symSVD$adjd) ## from adjMatrix to correlation matrix
          # outer optim, not spprec -> LMatrix recomputed from this for each rho  
          Lunique <- mat_sqrt(symSVD=symSVD) ## using $d not $adjd : $d=NaN should catch erroneous calls to mat_sqrt.
          if (inherits(Lunique,"try-error")) { ## mat_sqrt has try()'s and can return "try-error"'s
            print("correlation parameter was rho=:",rho,quote=FALSE) ## makes sense if mat_sqrt already issued some warning
            stop()
          }
          if ( adj_rho_is_inner_estimated ) { # contrived way of construction Lunique with the correct attributes.
            Lunique[] <- attr(Lunique,"symsvd")$u   ## "[] <- " keeps attributes... except for Matrix...
          }
        } else { ## cov_info_mat must exist
          if (processed$HL[1L]=="SEM") argsfordesignL$try.chol <- FALSE
          if (inherits(cov_info_mat,"dist")) {cov_info_mat <- proxy::as.matrix(cov_info_mat, diag=1)} ## else full matrix may be a COV matrix with non-unit diag
          Lunique <- do.call(mat_sqrt,c(list(m=cov_info_mat),argsfordesignL)) ## mat_sqrt has try()'s and can return "try-error"'s
          # Lunique is a tcrossfac: tcrossprod(mat_sqrt(X))=X
          if (inherits(Lunique,"try-error")) { 
            print("correlation parameters were:",quote=FALSE) ## makes sense if mat_sqrt already issued some warning
            print(unlist(argsfordesignL))    
            stop()
          }
        }
        attr(Lunique,"msd.arglist") <- msd.arglist ## NULL except for Matern, Cauchy
        attr(Lunique,"corr.model") <- corr_type
        attr(Lunique,"ranefs") <- paste(c(spatial_term)) ## essentiel pour la construction de ZAL! ## paste(c()) handles very long RHS
        AUGI0_ZX_envir$LMatrices[[rd]] <- Lunique
        attr(AUGI0_ZX_envir$LMatrices,"is_given_by")[rd] <- "AUGI0_ZX$envir" 
      }
    }
    if ( ! is.null(symSVD$adjd)) { 
      ## we should not replace attr(...$adjMatrices[[rd]], "symSVD") which has a given $d distinct from the new symSVD$d
      ## Moreover, currently we use only he new symSVD's $adjd
      attr(processed$corr_info$adjMatrices[[rd]],"adjd") <- symSVD$adjd ## 
    } 
  }
}

.normalize_IMRF <- function(ZAlist, # already ZA in non-IMRF input (or even IMRF input not to be normalized), 
                                    #   in contrast to .calc_normalized_ZAlist()
                            vec_normIMRF, 
                            Zlist=attr(ZAlist,"Zlist"),
                            strucList) {
  AMatrices <- attr(ZAlist,"AMatrices")
  for (rd in  seq_len(length(ZAlist))) { 
    if (vec_normIMRF[rd]) { 
      char_rd <- as.character(rd)
      colnams <- colnames(Zlist[[char_rd]]) 
      if ( ! setequal(rownames(AMatrices[[char_rd]]), colnams)) {
        stop(paste0("Any 'A' matrix must have row names that match the levels of the random effects\n",
                    "(i.e. the colnames of the 'Z' design matrix)"))
      } # ELSE:       
      AL <- AMatrices[[char_rd]] %*% strucList[[rd]]
      invnorm <- 1/sqrt(rowSums(AL^2)) # diag(tcrossprod...)
      normAL <- .Dvec_times_Matrix(invnorm, AMatrices[[char_rd]])
      ZAlist[[char_rd]] <- Zlist[[char_rd]] %id*% normAL[colnams,]
    }
  }
  return(ZAlist) ## with unchanged attributes
}

.calc_normalized_ZAlist <- function(Zlist, # creates ZA from Z and A, even for non-IMRF
                            AMatrices,
                            vec_normIMRF, 
                            strucList) {
  if (length(Zlist) && length(AMatrices)) {
    for (char_rd in  names(Zlist)) { # critically uses names here; Zlist and AMatrices are incomplete lists
      if ( ! is.null(Amatrix <- AMatrices[[char_rd]])) {
        colnams <- colnames(Zlist[[char_rd]])
        if (length(setdiff(colnams,rownames(Amatrix)))) {
          stop(paste0("Any 'A' matrix must have row names that match the levels of the random effects\n", 
                      "(i.e. the colnames of the 'Z' design matrix)"))
        } ## ELSE
        rd <- as.integer(char_rd) # I cannot yet assume strucList[[char_rd]] (nor vec_normIMRF[char_rd])
        if (vec_normIMRF[rd]) { 
          AL <- Amatrix %*% strucList[[rd]]
          invnorm <- 1/sqrt(rowSums(AL^2)) # diag(tcrossprod...)
          normAL <- .Dvec_times_Matrix(invnorm, Amatrix)
          Zlist[[char_rd]] <- Zlist[[char_rd]] %id*% normAL[colnams,,drop=FALSE]
        } else {
          is_incid <- attr(Zlist[[char_rd]],"is_incid")
          if (inherits(Amatrix,"pMatrix")) {
            # subsetting by rownames does not generally work on permutation matrices
            Amatrix <- as(Amatrix,"ngTMatrix")
          } else is_incid <- NULL
          Zlist[[char_rd]] <- Zlist[[char_rd]] %*% Amatrix[colnams,,drop=FALSE]
          attr(Zlist[[char_rd]],"is_incid") <- is_incid
        }
      }
    }
    attr(Zlist, "AMatrices") <- AMatrices
  } ## else attr(Zlist, "AMatrices") remains nULL
  return(Zlist) ## with other attributes unchanged
}

.init_assign_geoinfo <- function(processed, ZAlist, For=processed$For, corr_info=processed$corr_info, corr_types=corr_info$corr_types, 
                             exp_barlist, nrand=length(ZAlist), distMatrix, 
                             sparse_precision=processed$is_spprec, uniqueGeo) {
  if ( length(corr_types[ ! is.na(corr_types)])) {
    
    if (For=="HLfit") {
      ranef_string <- attr(ZAlist,"exp_ranef_strings")[ ! is.na(corr_types)][1L]
      stop(paste("Term",ranef_string,"not allowed in HLfit(). Try another fitting function such as fitme()."))
    } else {
      ## Cannot be unlisted bc cf ?unlist " Non-vector elements ... are not coerced, and so a list containing ... remains a list":
      attr(ZAlist,"exp_spatial_terms") <- exp_spatial_terms <- .findSpatial(barlist=exp_barlist, nomatch=NA, expand=TRUE) ## match ZAlist for predict()
      if (For != "fitmv") {
        geo_info <- structure(vector("list",nrand), names=seq_len(nrand)) ## each element will be an environment with typical elements $distMatrix, $uniqueGeo, $nbUnique
        cov_info_mats <- vector("list",nrand)
        for (it in seq_along(corr_types)) {
          corr_type <- corr_types[it]
          if ( ! is.na(corr_type)) {
            if (corr_type== "corrMatrix") { ## could cover all specifications of a constant 'cov' Matrix without correlation parameter
              ## .check_subset_corrMatrix has the effect that the corr or prec mat used in later computation is a permutation of the inputed one
              #  according to the order of columns of ZAlist[[it]] 
              ## In the spprec case corr_info$corrMatrices[[it]] already contains a precision matrix
              cov_info_mats[[it]] <- .check_subset_corrMatrix(corrMatrix=corr_info$corrMatrices[[it]],ZA=ZAlist[[it]]) ## correlation or precision...
              if ( inherits(cov_info_mats[[it]],"precision")) {
                # fast test by verif3 in test pedigree (no new cols) and longer test by test-Matern-spprec.R#9 (new cols)
                precmat <- cov_info_mats[[it]][["matrix"]]
                precnames <- colnames(precmat)
                is_incid <- attr(ZAlist[[it]],"is_incid") 
                ZAlist[[it]] <- .addrightcols_Z(Z=ZAlist[[it]], precnames)
                # now we are sure that they have the same names, only the orders are uncertain, so we can test order by any( != )
                ZAnames <- colnames(ZAlist[[it]])
                if (any(precnames!=ZAnames)) {
                  if (FALSE) { 
                    ZAlist[[it]] <- ZAlist[[it]][ ,precnames]
                  } else cov_info_mats[[it]] <- cov_info_mats[[it]][ZAnames,ZAnames]
                }
                attr(ZAlist[[it]],"is_incid") <- is_incid 
              }
            } else { ## all cases where geo_info (even empty) is needed 
              geo_info[[it]] <- new.env(parent=emptyenv())
              if ( ! is.null(distMatrix)) { ## user-provided distMatrix )
                if (is.list(distMatrix)) { ## spaMM3.0 extended syntax
                  geo_info[[it]]$distMatrix <- distMatrix[[it]] 
                } else geo_info[[it]]$distMatrix <- distMatrix 
              } # geo_info[[it]]$distMatrix may be modified below for Matern(LHS|.)
              if (corr_type == "AR1") {
                ugeo <-  attr(ZAlist[[it]],"uniqueGeo")
                if ( ! sparse_precision &&
                     ! is.null(dataordered_unique_levels <- attr(ZAlist[[it]],"dataordered_unique_levels"))) {
                  # The above test implies "if not already ! sparse when ZAlist was first evaluated"
                  # Subsetting ZAlist drops useful attributes => "uniqueGeo" must be secured in a more consistent place
                  rownames(ugeo) <- apply(ugeo,1L,paste0,collapse=":")
                  geo_info[[it]]$uniqueGeo <- ugeo[(dataordered_unique_levels), ,drop=FALSE]
                  #
                  is_incid <- attr(ZAlist[[it]],"is_incid") 
                  ZAlist[[it]] <- ZAlist[[it]][ ,(dataordered_unique_levels)]
                  # so that nice attributes are lost: "is_incid" "LHS_levels" "namesTerm" "dataordered_unique_levels" "AR1_block_n_u_h_s" "uniqueGeo"  
                  attr(ZAlist[[it]],"is_incid") <- is_incid 
                } else {
                  for (nam in names(ugeo)) if (is.factor(fac <- ugeo[[nam]])) ugeo[[nam]] <- as.character(levels(fac))[fac]
                  geo_info[[it]]$uniqueGeo <- ugeo
                }
                # } else if (corr_type %in% c("Matern","Cauchy") ) {
                # in that case we either
                # (1) have columns in ZAlist[[it]] not reordered: no further issues. This is the current design
                #  or
                # (2) have columns in ZAlist[[it]] not reordered in .calc_Zmatrix Then
                #     the cov_info_mats[[it]] info is not yet available at preprocessing time...
                #     See further comment in .assign_geoinfo_and_LMatrices_but_ranCoefs()
              } else if (corr_type %in% c("Matern","Cauchy")) { 
                if (attr(ZAlist[[it]],"namesTerm") != "(Intercept)") { # Matern(LHS|.)
                  geo_info[[it]]$activelevels <- activelevels <- which(colSums(ZAlist[[it]]!=0L)>0L)
                  geo_info[[it]]$distMatrix <- geo_info[[it]]$distMatrix[activelevels,activelevels, drop=FALSE]
                  is_incid <- attr(ZAlist[[it]],"is_incid") 
                  ZAlist[[it]] <- ZAlist[[it]][ ,activelevels,drop=FALSE]
                  attr(ZAlist[[it]],"is_incid") <- is_incid 
                } 
              }
              # not "else" because Matern/Cauchy must be combined with this!
              if ( ! is.null(uniqueGeo)) { ## user-provided uniqueGeo (no example anywhere! :-) )
                if (is.list(uniqueGeo)) { ## spaMM3.0 extended syntax
                  geo_info[[it]]$uniqueGeo <- uniqueGeo[[it]] 
                } else geo_info[[it]]$uniqueGeo <- uniqueGeo #
              } 
            }
          }
        }
        corr_info$cov_info_mats <- cov_info_mats
        processed$geo_info <- geo_info
      }
    }
  } else { ## no true_corr_types
    if ( For %in% c("HLCor","corrHLfit")) {
      stop("Correlation model not specified in 'formula': was valid in version 1.0 but not later.")
      ## Matern or corrMatrix were allowed without a tag then
    }
  }
  return(ZAlist)
}
