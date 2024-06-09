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
    if (is.null(geo_envir$uniqueGeo)) geo_envir$uniqueGeo <- blob$uniqueGeo  # 'typically' a data.frame if not provided to geo_envir 
                                                                             # in some other format earlier (see AR1 - time series code) 
    if (is.null(geo_envir$nbUnique)) geo_envir$nbUnique <- blob$nbUnique 
    if (is.null(geo_envir$notSameGrp)) geo_envir$notSameGrp <- blob$notSameGrp 
    if (is.null(geo_envir$coordinates)) geo_envir$coordinates <- blob$coordinates 
  }
  return(geo_envir) ## some elements may still be NULL
}

.inla_spde2_theta2phi0 <- function (spde, theta) {
  if (spde$n.theta > 0L) {
    return(exp(spde$param.inla$B0[, 1, drop = TRUE] + 
                 spde$param.inla$B0[,-1, drop = FALSE] %*% theta))
  } else return(exp(spde$param.inla$B0[, 1, drop = TRUE]))
}

.inla_spde2_theta2phi1 <- function (spde, theta) {
  if (spde$n.theta > 0L) {
    return(exp(spde$param.inla$B1[, 1, drop = TRUE] + 
                 spde$param.inla$B1[, -1, drop = FALSE] %*% theta))
  } else return(exp(spde$param.inla$B1[, 1, drop = TRUE]))
}

.inla_spde2_theta2phi2 <- function (spde, theta) {
  if (spde$n.theta > 0L) {
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

# called for IMRF() but not for MaternIMRFa
.inla_spde_precision_inla_spde2 <- function (spde, theta = NULL, 
                                            phi0 = .inla_spde2_theta2phi0(spde, theta), 
                                            phi1 = .inla_spde2_theta2phi1(spde, theta), 
                                            phi2 = .inla_spde2_theta2phi2(spde, theta), ...) 
{
  ## Matches formula for Q p.5 of https://www.jstatsoft.org/article/view/v063i19 (Lindgren & Rue 2015)
  if (spde$f$model != "spde2") {
    stop("spaMM only supports some internal inla models 'spde2'")
  }
  # if (FALSE) {
  #   D0 <- Diagonal(spde$n.spde, phi0) # T
  #   D1 <- Diagonal(spde$n.spde, phi1) # K^2
  #   D12 <- Diagonal(spde$n.spde, phi1 * phi2) # K^2 ....
  #   # for $MO = C, $M1 = G_1, $M2 = G_2
  #   Q <- (D0 %*% (D1 %*% spde$param.inla$M0 %*% D1 + D12 %*% spde$param.inla$M1 + 
  #                   t(spde$param.inla$M1) %*% D12 + spde$param.inla$M2) %*% D0)
  #   # 
  #   # => Starting with Matrix 1.5-2, the above Q is dgT, no longer dgC... argh
  #   return(as(Q,"CsparseMatrix"))
  #   # This probably comes fom  $MO, $MA, $M2 as produced by INLA code 
  #   # But the idea now is to preprocess them in .process_IMRF_bar() by .inla.param_dgT2dgC().
  # } else {
    ## This assumes that the matrices are dgC, otherwise see above. 
    D1M0D1 <- .Dvec_times_Matrix(phi1,spde$param.inla$M0)
    D1M0D1 <- .Matrix_times_Dvec(D1M0D1,phi1)
    D12M1 <- .Dvec_times_Matrix(phi1 * phi2,spde$param.inla$M1)
    Q <- D1M0D1 + D12M1 + t(D12M1) + spde$param.inla$M2
    Q <- .Dvec_times_Matrix(phi0, Q)
    Q <- .Matrix_times_Dvec(Q, phi0)
    ## Currently the phi vectors appear to be constant vectors 
    ##  so this could be further optimized, but...
    Q
  # }
}

# called for MaternIMRFa
.calc_IMRF_Qmat <- function(pars, grid_arglist, kappa, test=FALSE) {
  if (is.null(spde_info <- pars$model)) { 
    crossfac_Q <- .IMRFcrossfactor(xstwm=length(grid_arglist[[1]]), ystwm=length(grid_arglist[[2]]),kappa=kappa)
    sparse_Qmat <- crossprod(crossfac_Q) # same relation as for t_chol_Q ## automatically dsCMatrix
  } else {
    if (inherits(spde_info,"inla.spde2")) { # F I X M E only for dim=2 
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
    if (.hasSlot(tPmat,"margin") && tPmat@margin==2L) { # change introduced in Matrix 1.6.0 (argh)
      RRsP <- tPmat@perm
    } else RRsP <- sort.list(tPmat@perm) # older Matrix versions without indMatrix class
    colnames(tPmat) <- levelnames[RRsP]
    if (is.null(Amat)) {
      rownames(tPmat) <- colnames(processed$ZAlist[[rd]])
      # : when there is an A matrix, .calc_normalized_ZAlist() checks its names 
      processed$corr_info$AMatrices[[as.character(rd)]] <- structure(tPmat, permuted_Q=TRUE)
    } else {
      Amat <- .subcol_wAttr(Amat,j=RRsP, drop=FALSE) # Amat %*% tPmat
      attr(Amat,"perm") <- RRsP # used by .get_new_AMatrices()
      processed$corr_info$AMatrices[[as.character(rd)]] <- structure(Amat, permuted_Q=TRUE)
    }
    ZA <- processed$ZAlist[[rd]] %*% tPmat
    attr(ZA,"is_incid") <- attr(processed$ZAlist[[rd]],"is_incid")
    attr(ZA,"RHS_info") <- attr(processed$ZAlist[[rd]],"RHS_info") # cannot be modified by tPmat in composite case at least
    processed$ZAlist[[rd]] <- ZA
    .assign_ZAfix(processed)
  } # else ignore identity tPmat
}

.Lunique_info_from_Q_CHM <- function(processed, rd, Q_CHMfactor, 
                                     corr_type=processed$corr_info$corr_types[[rd]], 
                                     spatial_term=attr(processed$ZAlist,"exp_spatial_terms")[[rd]], 
                                     type) {
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
    # Q_CHMfactor of corr matrix: no direct relation to the Xaug nor to processed$AUGIO_ZX$I...
    Lunique <- solve(Q_CHMfactor,system="Lt", b=.sparseDiagonal(n=ncol(Q_CHMfactor), shape="g")) #  L_Q^{-\top}=LMatrix_correlation 
    ## Lunique is a dtCMatrix; it is for single correlated effect, and still used in HLfit_body; 
    #  Keeping it as dtCMatrix might be faster for some operations (but not .tcrossprod) (F I X M E test)
    ## whether it is used or not in MME_method (sXaug_...), a lot of other code still expects it
    attr(Lunique,"Q_CHMfactor") <- Q_CHMfactor 
    attr(Lunique, "type") <- type
  } 
  attr(Lunique,"corr.model") <- processed$corr_info$corr_types[[rd]]
  ####  attr(Lunique,"ranefs") <- paste(c(spatial_term)) ## essentiel pour la construction de ZAL! ## paste(c()) handles very long RHS
  Lunique # expanded, in composite case
}


.calc_Lunique_for_correl_algos <- function(processed, symSVD, rho, adj_rho_is_inner_estimated, argsfordesignL, 
                                           cov_info_mat, condnum=1) {
  if ( ! is.null(symSVD)) {
    if (! is.numeric(rho)) { # protection from crashing bug...
      stop("Invalid rho value passed to .calc_Lunique_for_correl_algos()")
    }
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
  } else if (inherits(cov_info_mat,"precision")) { ## cov_info_mat must exist
    ##  warning("Presumably inefficient computation in .calc_Lunique_for_correl_algos(). PLease contact the maintainer.")
    # This never occurred in the long tests [until *], presumably bc precision matrices are converted to corr mat
    # before reaching this point. 
    # But again, this does not happen in the tests...
    # [ * ... until a test has been added, forcing ARp to be fitted by CORR algos. => cf test-composite.R]
    Lunique <- solve(Matrix::chol(cov_info_mat$matrix)) # tcrossfac; no need for an environment storing chol_Q since the sparse matrix already does this.
  } else {
    if (processed$HL[1L]=="SEM") argsfordesignL$try.chol <- FALSE
    if (inherits(cov_info_mat,"dist")) {cov_info_mat <- proxy::as.matrix(cov_info_mat, diag=1)} ## else full matrix may be a COV matrix with non-unit diag
    #Lunique <- do.call(.sym_sqrt,c(list(A=cov_info_mat),argsfordesignL)) ## .sym_sqrt has try()'s and can return "try-error"'s
    Lunique <- do.call(mat_sqrt,c(list(m=cov_info_mat),argsfordesignL)) ## mat_sqrt has try()'s and can return "try-error"'s
    # Lunique is a tcrossfac: tcrossprod(mat_sqrt(X))=X
    if (inherits(Lunique,"try-error")) { 
      print("correlation parameters were:",quote=FALSE) ## makes sense if .sym_sqrt already issued some warning
      print(unlist(argsfordesignL))    
      stop()
    }
  }
  Lunique
}

..get_Q_CHMfactor__assign_template__perm_ZA <- function(processed, rd, sparse_Qmat, AUGI0_ZX_envir,
                                                        corr_type=processed$corr_info$corr_types[[rd]],
                                                        levels_type=processed$corr_info$levels_types[[rd]]) {
  #  I once had a problem with perm_Q=TRUE in test-predVar-Matern-corrMatrix -> predict(f2,.) so make sure to check this when changing the code
  Amat <- processed$corr_info$AMatrices[[as.character(rd)]]
  if (is.null(perm_Q <- .spaMM.data$options$perm_Q)) {
    ## ! Two distinct concepts: Whether the CHM is updateable and whether it is permuted. 
    ## But updating is useful ~ only when permutations are used.
    ## We reach here only if it is updateable,
    ## but perm_Q is here TRUE only in possibly more restrictive conditions *for user's perm_Q default=NULL*
    ## perm_Q is FALSE in remaining cases: Matern,Cauchy where any permuted precision matrix is presumably full without any useful pattern of zeros 
    ## _F I X M E__ could also depend on expected sparsity ? on ZA being identity?
    #
    ## for AR1 the 'AR1_block_u_h_ranges' info was lost in .ZA_update() [which permutes ZA cols according to Q_perm]
    ## But it seems that AR1_block_u_h_ranges can be used [it's used to create the unpermuted Q 
    ## from which CHM updates are computed]. The fit is OK but the predict() fails
    ## if perm_Q is forced to TRUE. Cf tests composite-extra
    # p1s <- predict(fit1sp, newdata=fit1$data)   # AR1
    # pps <- predict(fitpsp, newdata=fit1$data)   # ARp
    ##
    ## I previously had distinct conditions if not composite, identified by:
    # not_composite <- ! processed$ranCoefs_blob$is_composite[[rd]]
    ## In composite case, the sparse_Qmat is already the Kronecker product here,
    ## (whether adjacency, IMRF or MaternIMRFa).
    ##
    ## => perm_Q=TRUE not OK for AR1 and ARp at least as explained above (___F I X M E____)
    ## ANd the following condition isquite speculative for yet-untested combinations of Amat and corr_type
    perm_Q <- (
      corr_type %in% c("adjacency") ||
        ( ! is.null(Amat) && # Amat condition => to get IMRFs notably
            levels_type != "time_series" 
        )
    )
  } 
  Q_CHMfactor <- Cholesky(sparse_Qmat,LDL=FALSE,perm=perm_Q) 
  AUGI0_ZX_envir$precisionFactorList[[rd]]$template <- Q_CHMfactor
  if (perm_Q) .ZA_update(rd, Q_CHMfactor, processed, Amat) # else ZA is not permuted, there is no permuted_Q attribute, so the G matrix will be constructed from unpermuted ZA and unpermuted sparse_Qmat
  Q_CHMfactor
}

.get_Q_CHMfactor__assign_template__perm_ZA <- function(AUGI0_ZX_envir, rd, processed, sparse_Qmat) {
  if (is.null(template <- AUGI0_ZX_envir$precisionFactorList[[rd]]$template)) { ## occurs if $update_CHM is FALSE, OR first comput. of CHM, OR not yet all $updateable[rd]
    get_template <- .spaMM.data$options$update_CHM && AUGI0_ZX_envir$updateable[rd]
    if (get_template) {
      Q_CHMfactor <-  ..get_Q_CHMfactor__assign_template__perm_ZA(processed, rd, sparse_Qmat, AUGI0_ZX_envir)
    } else Q_CHMfactor <- Cholesky(sparse_Qmat,LDL=FALSE,perm=FALSE) 
    # perm=TRUE without saving the result in a template should be a problem: 
    #  as(Q_CHMfactor, "sparseMatrix") might have different permutations over calls, while 
    # .ZA_update() is valid for permuting ZA only once but not repeatedly (original unpermuted ZA not kept)
  } else { # template for updating already exists
    Q_CHMfactor <- Matrix::.updateCHMfactor(template, parent=sparse_Qmat, mult=0) 
  }
}

.get_t_chol_Q_AR1 <- function(AUGI0_ZX_envir, ARphi, char_rd, processed, rd) {
  types <- AUGI0_ZX_envir$finertypes
  ## the dim of Q *HAS* to match that of ZA
  AR1_block_u_h_ranges <- attr(processed$ZAlist[[rd]],"RHS_info")$AR1_block_u_h_ranges # block size[S for single AR1 term if nested]
  ilist <- jlist <- xlist <- vector("list",length(AR1_block_u_h_ranges))
  cum_ncol <- 0L
  for (bloc in seq_along(AR1_block_u_h_ranges)) {
    n_u_h_block <- diff(AR1_block_u_h_ranges[[bloc]])+1L
    triplets <- .calc_AR1_t_chol_Q_block(n_u_h_block, ARphi=ARphi) ## L block(s) of given size(s)
    ilist[[bloc]] <- cum_ncol + triplets$i
    jlist[[bloc]] <- cum_ncol + triplets$j 
    xlist[[bloc]] <- triplets$x
    cum_ncol <- cum_ncol + n_u_h_block
  }
  t_chol_Q <- sparseMatrix(dims = c(cum_ncol,cum_ncol), i=unlist(ilist), j=unlist(jlist),
                           x = unlist(xlist), triangular=TRUE)
  t_chol_Q
}

.get_cov_info_geostat <- function(control.dist, char_rd, rho, spatial_term, processed, rd, ranPars) {
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
                               dist_method_rd=dist_method_rd)  # wrapper for .get_dist_nested_or_not()
    # For e.g. Matern(1|longitude+latitude %in% grp) the returned $uniqueGeo may be 'e_uniqueGeo within 3 cols
    msd.arglist$uniqueGeo <- geo_envir$uniqueGeo[,attr(geo_envir$uniqueGeo,"coord_within"),drop=FALSE] 
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
    msd.arglist$distMatrix <- geo_envir$distMatrix   # may be a "dsCDIST" (dsC-with-NA-on-diag; 0's being interpreted differently) 
    # "dsCDIST" must be handled by make_scaled_dist() (nothing specific yet needed);
    # "dsCDIST" must be handled by $calc_corr_from_dist() -> MaternCorr() -> MaternCorr.dsCMatrix() which returns a standard dsC whithout NA 
    cov_info_mat <- do.call("make_scaled_dist",msd.arglist)
  }
  if ( nrow(cov_info_mat)>1 ) { ## >1 locations
    cov_info_mat <- processed$corr_info$corr_families[[rd]]$calc_corr_from_dist(ranFix=ranPars, char_rd=char_rd, distmat=cov_info_mat)
  } # else a single location => dist_mat should be dist(0) and make_scaled_dist was modified to that effect
  list(cov_info_mat=cov_info_mat, msd.arglist=msd.arglist)
}


.precision2cov <- function(sparse_Qmat) {
  if (ncol(sparse_Qmat)>200L) { # a bit of an ad hoc guess based on ohio...
    cov_info_mat <- as.matrix(solve(sparse_Qmat)) # not clearly distinct from chol2inv(chol(sparse_Qmat)) (ohio)
  } else cov_info_mat <- solve(as.matrix(sparse_Qmat)) # faster than chol2inv(chol(sparse_Qmat)) ! (ohio)
  cov_info_mat
}

.AR1_assign_geoinfo_and_LMatrices_but_ranCoefs <- function(AUGI0_ZX_envir, rd, processed, sparse_Qmat, ranPars, char_rd, corr_type, spatial_term) {
  Q_CHMfactor <- try(.get_Q_CHMfactor__assign_template__perm_ZA(AUGI0_ZX_envir, rd, processed, sparse_Qmat))
  if (inherits(Q_CHMfactor, "try-error")) { 
    cP <- unlist(ranPars$corrPars[[char_rd]])
    mess <- paste(names(cP),"=",cP, collapse=", ")
    stop(paste0("A numerical problem occurred for random effect ",char_rd," with parameter(s) " ,mess,". Try to restrict the parameter range."))
  }
  #
  # If we use permuted Chol AND a template, ZA has been permuted (once). Then we must permute each new sparse_Qmat, by 
  chol_Q <- as(Q_CHMfactor,"CsparseMatrix")
  permuted_Q <- attr(processed$corr_info$AMatrices[[as.character(rd)]],"permuted_Q") # clumsy but need to get info from one-time code
  if (identical(permuted_Q,TRUE)) sparse_Qmat <- tcrossprod(chol_Q) 
  #
  AUGI0_ZX_envir$precisionFactorList[[rd]]$Qmat <- sparse_Qmat # should by *dsC*  
  # $Qmat <- sparse_Qmat will be used together with ZA independently from the CHM to construct the Gmat
  AUGI0_ZX_envir$precisionFactorList[[rd]]$chol_Q <- chol_Q # Linv
  AUGI0_ZX_envir$LMatrices[[rd]] <- .Lunique_info_from_Q_CHM(
    processed=processed, rd=rd, Q_CHMfactor=Q_CHMfactor, 
    corr_type=corr_type, spatial_term=spatial_term, # presumably the default args but they are available 
    type="from_Q_CHMfactor")
  attr(AUGI0_ZX_envir$LMatrices,"is_given_by")[rd] <- "from_Q_CHMfactor" 
} # no used return.


# Only called in HLCor_body() where it replaces the direct call to .init_AUGI0_ZX_envir_spprec_info() from other _body() functions
.assign_geoinfo_and_LMatrices_but_ranCoefs <- function(processed, corr_types, spatial_terms, 
                                                       ranPars, control.dist, argsfordesignL,
                                                       corr_families=processed$corr_info$corr_families) {
  # * writes into geo_envir <- .get_geo_info(...)
  # * modifies processed$AUGI0_ZX$envir by .init_AUGI0_ZX_envir_spprec_info(...) 
  # * computes processed$AUGI0_ZX$envir$LMatrices except for ranCoefs (the latter being filled in <HLfit_body> 
  #   by .wrap_precisionFactorize_ranCoefs())
  if (processed$is_spprec) .init_AUGI0_ZX_envir_spprec_info(processed) 
  AUGI0_ZX_envir <- processed$AUGI0_ZX$envir
  corr_info <- processed$corr_info
  for (rd in seq_along(corr_types)) {
    corr_type <- corr_types[[rd]]
    symSVD <- NULL ## reinitialize as rd is tested within the loop
    msd.arglist <- NULL
    if ( ! is.na(corr_type)) { # including composite ones, otherwise need:   && ! AUGI0_ZX_envir$finertypes[rd]=="ranCoefs") {
      char_rd <- as.character(rd)
      spatial_term <- spatial_terms[[rd]] 
      ###
      # Provide cov_info_mat and/or sparse_Qmat (with many exceptions) and sometimes symSVD or msd.arglist
      if (corr_type == "adjacency") {
        adjMatrix <- corr_info$adjMatrices[[rd]] # dsC whether spprec or not
        symSVD <- attr(adjMatrix,"symSVD")
        rho <- .get_cP_stuff(ranPars,"rho",which=char_rd)
        adj_rho_is_inner_estimated <- (is.null(rho) ) ## can occur in direct call of HLCor 
        if ( ! adj_rho_is_inner_estimated) { # typical fitme() call
          sparse_Qmat <- - rho * adjMatrix
          sparse_Qmat <- .dsCsum(A=sparse_Qmat, B=attr(adjMatrix,"dsCdiag"), keep_names = TRUE) #diag(sparse_Qmat) <- diag(sparse_Qmat)+1 # a maybe-not-yet-perfect solution to the poor perf of 'diag<-' which is not doc'ed for dsCMatrix
          cov_info_mat <- list(matrix=sparse_Qmat)
          class(cov_info_mat) <- c("precision", class(cov_info_mat))
          if ( ! processed$is_spprec && is.null(symSVD)) cov_info_mat <- .precision2cov(sparse_Qmat=sparse_Qmat)
          AUGI0_ZX_envir$updateable[rd]=(rho!=0)
        } else { # inner estimation of adjacency rho => only symSVD
          if (is.null(symSVD)) {
            ## Direct call of HLCor (SEM or not), 
            ## I also wrote that this could occur in fitme/corrHLfit if(list(processed)) "bc symSVD added to proc1 (i.e. not to all proc's)"
            ##   but there is not good tests of this case...
            symSVD <- .provide_AR_factorization_info(
              adjMatrix, 
              sparse_precision=processed$is_spprec, # this must be FALSE for inner estimation of adjacency rho (and there is bug-catching code for this)
              corr.model=corr_type)
            attr(corr_info$adjMatrices[[rd]],"symSVD") <- symSVD
            #
            # typically an *initial* value:
            rho <- attr(ranPars,"init.HLfit")$corrPars[[char_rd]]$rho # a bit unclear, but we will compute Lmatrix from this 'initial value' below
            if (is.null(rho)) stop("'rho' missing. Contact the maintainer") # protection to avoid programming errors resulting in segfaults.
          } 
          # symSVD may be modified below
        }
      } else if (corr_type =="IMRF") {
        # sparse_Qmat
        sparse_Qmat <- .calc_IMRF_Qmat(pars=attr(attr(spatial_term,"type"),"pars"), 
                                       grid_arglist=attr(processed$corr_info$AMatrices[[char_rd]],"grid_arglist"), # promise for spde case
                                       kappa=.get_cP_stuff(ranPars,"kappa",which=char_rd))
        # cov_info_mat
        if ( ! processed$is_spprec) cov_info_mat <- .precision2cov(sparse_Qmat) # solve(sparse_Qmat) # ~ never used ( ~ always spprec);
      } else if (corr_type =="SAR_WWt") { 
        rho <- .get_cP_stuff(ranPars,"rho",which=char_rd)
        adjMatrix <- corr_info$adjMatrices[[rd]]
        # sparse_Qmat
        sparse_Qmat <- - rho * adjMatrix
        sparse_Qmat <- .dsCsum(sparse_Qmat, attr(adjMatrix,"dsCdiag")) #diag(sparse_Qmat) <- diag(sparse_Qmat)+1 # a maybe-not-yet-perfect solution to the poor perf of 'diag<-' which is not doc'ed for dsCMatrix
        # cov_info_mat
        if ( ! processed$is_spprec) {
          UDU. <- attr(adjMatrix,"UDU.")
          if (is.null(UDU.)) {
            # cf comments on symSVD above
            cov_info_mat <- list(matrix=sparse_Qmat)
            class(cov_info_mat) <- c("precision", class(cov_info_mat))
          } else {
            cov_info_mat <- .ZWZt(UDU.$u, 1/(1-rho*UDU.$d)) # UDU.$u %*% sweep(UDU.$u.,MARGIN=1,1/(1-rho*UDU.$d),`*`) 
          }
          cov_info_mat <- .tcrossprodCpp(cov_info_mat,NULL)
        }
      }  else if ( identical(corr_families[[rd]]$levels_type,"time_series") ) { # AR1, ARp, ARMA 
        if (corr_type=="AR1") {
          # sparse_Qmat never needed
          # cov_info_mat
          if (processed$is_spprec) {
            # there is ad hoc code below where t_chol_Q is computed first stored in AUGI0_ZX_envir$precisionFactorList[[rd]])
            #    and other matrices may be deduced from it.
            # the formals of a wrapper for such code would be (ranPars, char_rd, AUGI0_ZX_envir, processed, rd, corr_type, spatial_term)
          } else {
            geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("distMatrix"), 
                                       dist_method_rd=control.dist[[char_rd]]$dist.method)
            # only the non-diagonal elements as a dist object: (tailored to AR1 case)
            cov_info_mat <- corr_families[[rd]]$calc_corr_from_dist(ranFix=ranPars, char_rd=char_rd, 
                                                                    distmat=geo_envir$distMatrix)
          } 
        } else { ## other time_series
          # the cov_info_mat is reduced (! spprec) or not (spprec) according to distances present in the data
          # => .corrfamily2Matrix() when not reduced, else $calc_corr_from_dist(., distmat=geo_envir$distMatrix)
          if ( processed$is_spprec) {
            cov_info_mat <- .corrfamily2Matrix(corr_info$corr_families[[rd]], ranPars$corrPars[[char_rd]], AUGI0_ZX_envir, rd)
          } else {
            geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("distMatrix"),
                                       dist_method_rd=control.dist[[char_rd]]$dist.method)
            # NOT a dist object:
            cov_info_mat <- corr_families[[rd]]$calc_corr_from_dist(ranFix=ranPars, char_rd=char_rd,
                                                                    distmat=geo_envir$distMatrix)
          }
        }
      } else  if (corr_type %in% c("Matern","Cauchy")) {
        rho <- .get_cP_stuff(ranPars,"rho",which=char_rd)
        # sparse_Qmat not provided, will be deduced from cov_info_mat later => cov_info_mat needed for spprec too
        # cov_info_mat (for spprec too)
        cov_info_geostat <- .get_cov_info_geostat(control.dist, char_rd, rho, spatial_term, processed, rd, ranPars)
        cov_info_mat <- cov_info_geostat$cov_info_mat
        msd.arglist <- cov_info_geostat$msd.arglist
      } else if (corr_type== "corrMatrix") { 
        if (processed$is_spprec) {
          #   .init_assign_geoinfo() has set a "blob" attribute to cov_info_mats[[rd]], with promises for Lunique, etc
          # (if composite ranCoefs) {
          #   the promises are there to be used.  
          # } else .init_AUGI0_ZX_envir_spprec_info() runs precisionFactorList[[rd]] <- .calc_corrMatrix_precisionFactor__assign_Lunique(processed, rd)
        } else cov_info_mat <- corr_info$cov_info_mats[[rd]] ## correlation or precision or dist...
      } else if (corr_type== "corrFamily") { # fallback case
        # see comments in [time_series, not AR1] case.
        # Only .corrfamily2Matrix() has been considered for a long time. This was 
        # presumably inappropriate in some then-untested cases, but we cannot revise 
        # the code in a way that would break uses of corrFamilies without $calc_corr_from_dist. 
        # e.g. (diallel_fit <- fitme(z ~1 +diallel(1|id1+id2), data=dyaddf)).
        # Ideally (___F I X M E____) one could generate  reasonable default $calc_corr_from_dist
        # a by redefining $Cf so that it can handle newlevels, and then try:
        # calc_corr_from_dist <- function(ranFix, char_rd, distmat, ...) { # The AR1 code use distance matrices to handle the nested AR1 case...
        #   parvec <- .fill_parvec(parvec=ranFix$corrPars[[char_rd]], fixed=fixed, npar=p)
        #   dimnams <- dimnames(distmat)
        #   levelrange <- range(as.integer(.unlist(dimnams)))
        #   Qmat <- <redefined Cf>(parvec=parvec, newlevels=seq(levelrange[1L],levelrange[2L]))
        #   corr_mat <- chol2inv(chol(Qmat)) # .precision2cov(Qmat) 
        #   if (inherits(distmat,"dist")) {
        #     corr_mat[dimnams,dimnams] 
        #   } else {
        #     corr_mat[dimnams[[1]],dimnams[[2]]]
        #   }
        # }
        if ( processed$is_spprec ||
             inherits(corr_families[[rd]]$calc_corr_from_dist,"stub")) {
          cov_info_mat <- .corrfamily2Matrix(corr_info$corr_families[[rd]], ranPars$corrPars[[char_rd]], AUGI0_ZX_envir, rd)
        } else {
          geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("distMatrix"),
                                     dist_method_rd=control.dist[[char_rd]]$dist.method)
          # Should NOT be a dist object (currently this would break a possible call to .calc_denseness()):
          cov_info_mat <- corr_families[[rd]]$calc_corr_from_dist(ranFix=ranPars, char_rd=char_rd,
                                                                  distmat=geo_envir$distMatrix)
        }
        #
        # if ( ( ! processed$is_spprec) && inherits(cov_info_mat,"precision")) {
        ## This would occur only when an explicit algebra=<"spcorr"|"decorr"> has been requested, as otherwise 
        ## .determine_spprec() selects spprec when any(.unlist(lapply(corr_info$corr_families,`[[`, "type"))=="precision").
        ## Hence the following warning does not seem useful. 
        #   mess <- paste("Inefficient code: corrFamily's $f() returns a precision matrix, but selected algebraic method is",
        #           if (processed$QRmethod=="sparse") {"'spcorr'"} else {"'decorr'"})
        #   warning(mess)
        #   cov_info_mat <- .precision2cov(cov_info_mat$matrix) # and this is useless...
        # }
        #
      }
      ###
      #
      ## Provide Lunique if not already available (and optionally additional stuff)
      if (processed$is_spprec) { 
        #### SPPREC
        if (corr_type=="AR1") {
          # AR1 = special case where Qmat is deduced from a crossfac, rather than a (t)crossfac deduced from Qmat.
          ARphi <- .get_cP_stuff(ranPars,"ARphi",which=char_rd)
          t_chol_Q <- .get_t_chol_Q_AR1(AUGI0_ZX_envir, ARphi, char_rd, processed, rd)
          sparse_Qmat <- crossprod(t_chol_Q)
          AUGI0_ZX_envir$updateable[rd]=(ARphi!=0) # to be specifically used for G_CHMfactor
          if (AUGI0_ZX_envir$finertypes[rd]=="ranCoefs") { # composite case! 
            # .process_ranCoefs() will build the info for the full composite model, including the LMatrix
            AUGI0_ZX_envir$precisionFactorList[[rd]]$RHS_Qmat <- sparse_Qmat # should by *dsC*  
            processed$ranCoefs_blob$new_compos_info[rd] <- TRUE
          } else #if (TRUE) {
            .AR1_assign_geoinfo_and_LMatrices_but_ranCoefs(AUGI0_ZX_envir, rd, processed, sparse_Qmat, ranPars, char_rd, corr_type, spatial_term)
          # } else {
          ## Non-CHM update version
          #   AUGI0_ZX_envir$precisionFactorList[[rd]] <- list(chol_Q=t(t_chol_Q), # Linv
          #                                                    Qmat=sparse_Qmat)
          #   Lunique <- solve(t_chol_Q) # that is super slow compared to using .Lunique_info_from_Q_CHM with represents the LMatrix as a CHMfactor !!
          #   attr(Lunique, "type") <- "from_AR1_specific_code"
          #   attr(Lunique,"corr.model") <- corr_type
          #   ####   attr(Lunique,"ranefs") <- paste(c(spatial_term)) ## essentiel pour la construction de ZAL! ## paste(c()) handles very long RHS
          #   AUGI0_ZX_envir$LMatrices[[rd]] <- Lunique
          #   attr(AUGI0_ZX_envir$LMatrices,"is_given_by")[rd] <- "AUGI0_ZX$envir" 
          # }
        } else if (corr_type=="corrMatrix") { # both composite and non-composite one !
          #.init_AUGI0_ZX_envir_spprec_info() -> .calc_corrMatrix_precisionFactor__assign_Lunique() has done the work
        } else { ## General spprec code that should be correct for AR1 and corrMatrix too
          #### (1) Provide sparse_Qmat for all parametric correlation models (solve(sparse_Qmat) gives the correlation matrix)
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
          } else if (corr_type=="corrFamily") {
            if (is.null(cov_info_mat)) { # e.g., ranGCA...
              nc <- ncol(processed$ZAlist[[rd]])
              sparse_Qmat <- .symDiagonal(n=nc)
              # chol_Q=new("dtCMatrix",i= 0:(nc-1L), p=0:(nc), Dim=c(nc,nc),x=rep(1,nc)) )
            } else if (inherits(cov_info_mat,"precision")) {
              sparse_Qmat <- cov_info_mat$matrix
            } else sparse_Qmat <- as_precision(cov_info_mat)$matrix 
          } else if (corr_type== "IMRF") {
            # Remember that we need dtCMatrix'es 'chol_Q' so that bdiag() gives a dtCMatrix
            # Hence use next general code to produce precisionFactorList[[rd]]
          } else stop("Some error occurred (inner estimation of adjacency rho with requested sparse precision ?)") 
          #
          #### (2) (usually) build FROM sparse_Qmat, providing precisionFactorList[[rd]] as expected by .reformat_Qmat_info(),
          #        and the LMatrix.     
          if (AUGI0_ZX_envir$finertypes[rd]=="ranCoefs") { # composite case! 
            # .process_ranCoefs() will build the Q_CHMfactor info for the full composite model, including the LMatrix
            AUGI0_ZX_envir$precisionFactorList[[rd]]$RHS_Qmat <- sparse_Qmat # should by *dsC*  
            processed$ranCoefs_blob$new_compos_info[rd] <- TRUE
          } else {
            Q_CHMfactor <- try(.get_Q_CHMfactor__assign_template__perm_ZA(AUGI0_ZX_envir, rd, processed, sparse_Qmat))
            if (inherits(Q_CHMfactor, "try-error")) { 
              cP <- unlist(ranPars$corrPars[[char_rd]])
              mess <- paste(names(cP),"=",cP, collapse=", ")
              mess <- paste0("A numerical problem occurred for random effect ",char_rd,
                     " with parameter(s) " ,mess,".\n Try to restrict the parameter range")
              if (corr_type=="IMRF" && 
                  "kappa" %in% names(cP) && 
                  environment(.try_RSpectra)$RSpectra_warned
              ) mess <- paste0(mess, " (installing the 'RSpectra' package may also help)")
              stop(paste0(mess, "."))
            }
            #
            # If we use permuted Chol AND a template, ZA has been permuted (once). Then we must permute each new sparse_Qmat, by 
            permuted_Q <- attr(processed$corr_info$AMatrices[[as.character(rd)]],"permuted_Q") # clumsy but need to get info from one-time code
            if (identical(permuted_Q,TRUE)) sparse_Qmat <- tcrossprod(as(Q_CHMfactor,"CsparseMatrix")) 
            #
            # $Qmat <- sparse_Qmat will be used together with ZA independently from the CHM to construct the Gmat
            AUGI0_ZX_envir$precisionFactorList[[rd]]$Qmat <- sparse_Qmat # should by *dsC*; possibly permuted version 
            # of the factorized matrix. Info about the permutation can be found in the template used for updating.  
            AUGI0_ZX_envir$precisionFactorList[[rd]]$chol_Q <- as(Q_CHMfactor, "CsparseMatrix") # Linv
            
            processed$AUGI0_ZX$envir$LMatrices[[rd]] <- .Lunique_info_from_Q_CHM(
              processed=processed, rd=rd, Q_CHMfactor=Q_CHMfactor, 
              corr_type=corr_type, spatial_term=spatial_term, # presumably the default args but they are available 
              type="from_Q_CHMfactor")
            # processed$AUGI0_ZX$envir$LMatrices[[rd]] <- .Lunique_info_from_Q_CHM is for all correlated ranefs ./.
            # not specifically for composite ones. For the latter is as the dimension of the kronecker product.
            attr(processed$AUGI0_ZX$envir$LMatrices,"is_given_by")[rd] <- "from_Q_CHMfactor" 
          }
        }
      } else { 
        #### sparse or dense CORRELATION algorithms. For fixed "corrMatrix" => this IS executed but only once.
        if (processed$ranCoefs_blob$is_composite[rd]) { ## ~ if(AUGI0_ZX_envir$finertypes[rd]=="ranCoefs") 
          # => provide info to be used by .process_ranCoefs().
          # In that specific case, (CORREL algo, composite ranef), the RHS Lunique is stored 
          # in attr(corr_info$cov_info_mats[[rd]],"blob")$Lunique, 
          # not in AUGI0_ZX_envir$LMatrices[[rd]] where 'full' Lunique are stored.
          #
          NOT_cM_cT <- (corr_type != "corrMatrix")
          if (NOT_cM_cT) corr_info$cov_info_mats[[rd]] <- "'cov_info_mat' not (yet) stored" # non-NULL to allow attributes:
          if (NOT_cM_cT || 
              is.null(attr(corr_info$cov_info_mats[[rd]],"blob")$Lunique) # TRUE only the 1st time this block is reached 
             ) {
              # it would be nice to use a delayedAssign as for spprec. But how to deal with the extra arguments?
              Lunique <- .calc_Lunique_for_correl_algos(processed, symSVD, rho, adj_rho_is_inner_estimated, argsfordesignL, 
                                                        cov_info_mat)
              attr(corr_info$cov_info_mats[[rd]],"blob")$Lunique <- Lunique 
              attr(AUGI0_ZX_envir$kron_Y,"is_given_by")[rd] <- "AUGI0_ZX$envir"
              processed$ranCoefs_blob$new_compos_info[rd] <- TRUE
          }
        } else { # NOT composite. --- Block including the final if (ok) {...}
          if ( corr_type=="corrFamily" && # test TRUE for both un- and registered corrFamily's
               ( is.null(cov_info_mat) || # bc in that case we don't want the alternative to be tried
                 ( ! inherits(cov_info_mat,"precision") &&
                   .calc_denseness(cov_info_mat,relative=TRUE) < 0.15 )
               )
          ) {
            if ( ok <- ! is.null(cov_info_mat)) { # quite contrived. ranGCA and ARMA are test cases with opposite bugs...
              ## Allowing Cholesky updates for SPARSE cov_info_mat
              corr_CHMfactor <- .provide_CHMfactor(cov_info_mat, CHMtemplate=corr_info$corr_families[[rd]]$corr_CHM_template)
              r <- Matrix::expand(corr_CHMfactor) # return type of expand not formally doc'ed, but example from doc => $P is a 'pMatrix' (=> .crossprod() probably not suitable) and $L a 'dtCMatrix'
              Lunique <- crossprod(r$P,r$L) # Lunique no longer triangular ! 
              if (processed$AUGI0_ZX$envir$updateable[rd] && # : initialized to FALSE by .init_AUGI0_ZX() for correlation algos, and
                  # checked/updated for each new sparse Cf(tpar), so this is relatively safe (but see comment on check).
                  # For dense cov_info_mat this block is not run given the .calc_denseness() check. 
                  is.null(corr_info$corr_families[[rd]]$corr_CHM_template)) {
                corr_info$corr_families[[rd]]$corr_CHM_template <- corr_CHMfactor
                if (processed$QRmethod=="sparse" && .calc_denseness(r$L,relative=TRUE) > 0.15) { 
                  cat(crayon::green("Sparse-matrix methods have been selected, but it seems that setting algebra='decorr' would be better.\n"))
                  ## If r$L has high denseness, the slow step by "spcorr" will be the Matrix::qr(<sXaug>) step, not the Cholesky updates for corr_CHMfactor 
                  ## The follwing fix does not really work bc previous $QRmethod has been used by 
                  ## .preprocess -> .init_AUGI0_ZX(.,as_mat=.eval_as_mat_arg(processed)) with some irreversible effects (see e.g. comments on
                  ## type of Xscal <- .make_Xscal(...) in .solve_IRLS_as_ZX()).
                  ## so we would have to reevaluate this preprocessing step (which seems risky (but when do we first need it ?)).
                  #
                  # cat(crayon::green("Using algebra='decorr' instead\n"))
                  # processed$QRmethod <- "dense"
                  # processed$as_matrix <- NULL # will be recomputed by  .eval_as_mat_arg()
                } else if (processed$QRmethod=="dense" && .calc_denseness(r$L,relative=TRUE) < 0.15) {
                  cat(crayon::green("Dense-matrix methods have been selected, but it seems that setting algebra='spcorr' would be better.\n"))
                }
              }
              # Other cases typically not dsC, (except hacked type dscDIST). But any 'atypical" case worth considering ? 
            }
          } else if (ok <- (corr_type != "corrMatrix")) { # include corrFamilies
            Lunique <- .calc_Lunique_for_correl_algos(processed, symSVD, rho, adj_rho_is_inner_estimated, argsfordesignL, 
                                                      cov_info_mat)
          } else if (ok <- is.null(AUGI0_ZX_envir$LMatrices[[rd]])) { # *once* for corrMatrix !
            Lunique <- .calc_Lunique_for_correl_algos(processed, symSVD, rho, adj_rho_is_inner_estimated, argsfordesignL, 
                                                      cov_info_mat) # wraps mat_sqrt() as cov_info_mat is still dgC, not dsC
          }
          if (ok) {
            attr(Lunique,"msd.arglist") <- msd.arglist ## NULL except for Matern, Cauchy
            attr(Lunique,"corr.model") <- corr_type
            #### attr(Lunique,"ranefs") <- paste(c(spatial_term)) ## essentiel pour la construction de ZAL! ## paste(c()) handles very long RHS
            AUGI0_ZX_envir$LMatrices[[rd]] <- Lunique
            attr(AUGI0_ZX_envir$LMatrices,"is_given_by")[rd] <- "AUGI0_ZX$envir" 
          }
        }
      }
    } #  condition on corr_type...
    if ( ! is.null(symSVD$adjd)) { 
      ## we should not replace attr(...$adjMatrices[[rd]], "symSVD") which has a given $d distinct from the new symSVD$d
      ## Moreover, currently we use only he new symSVD's $adjd
      attr(corr_info$adjMatrices[[rd]],"adjd") <- symSVD$adjd ## 
    } 
  } # rd loop
}

# Recomputes ZA with A modified as function of L such that AL is tcrossfac of correlation matrix 
.normalize_IMRF_ZA <- function(Z, A, L, colnams=NULL) {
  if (is.null(Z)) return(NULL) #this occurs in mv fits .calc_ZAlist_newdata_mv() -> .calc_normalized_ZAlist(Zlist ...) -> here
                               # where Zlist may have some NULL elements.
  # : it's no longer clear when this L is the tcross factor (with the Q_CHMfactor as attribute...) or is the Q_CHMfactor
  # Maybe it is the Q_CHMfactor post-fit.
  if (inherits(L,"dCHMsimpl")) L <- solve(L, system="Lt", b=.sparseDiagonal(n=ncol(L), shape="g"))  
  AL <- A %*% L
  invnorm <- 1/sqrt(rowSums(AL^2)) # diag(tcrossprod...)
  normAL <- .Dvec_times_Matrix(invnorm, A)
  if (is.null(colnams)) {
    ZA <- Z %id*% normAL
  } else ZA <- Z %id*% normAL[colnams,,drop=FALSE]
  attr(ZA,"is_incid") <- FALSE
  ZA
}

#call by HLCor_body
.normalize_IMRF <- function(processed, # with $ZAList already ZA in non-IMRF input (or even IMRF input not to be normalized), 
                                    #   in contrast to .calc_normalized_ZAlist()
                            vec_normIMRF, 
                            Zlist=attr(ZAlist,"Zlist"),
                            strucList) {
  ZAlist <- processed$ZAlist
  AMatrices <-processed$corr_info$AMatrices
  for (rd in  seq_len(length(ZAlist))) { 
    if (vec_normIMRF[rd]) { 
      char_rd <- as.character(rd)
      colnams <- colnames(Zlist[[char_rd]]) 
      Amatrix <- AMatrices[[char_rd]]
      if ( ! setequal(rownames(Amatrix), colnams)) {
        # col Z must be = rows of A
        stop(paste0("Any 'A' matrix must have row names that match all the levels of the random effects\n",
                    "(i.e. the colnames of the 'Z' design matrix)"))
      } # ELSE:       
      ZAlist[[char_rd]] <- .normalize_IMRF_ZA(Z=Zlist[[char_rd]], A=Amatrix, L=strucList[[rd]], colnams=colnams)
    }
  }
  return(ZAlist) ## with unchanged attributes
}

..calc_normalized_ZA <- function(Z_, normIMRF, Amatrix, L) {
  if (normIMRF) {
    Z_ <- .normalize_IMRF_ZA(Z=Z_, A=Amatrix, L=L, colnams=NULL)
  } else {
    mostAttrs <- attributes(Z_)
    is_incid <- mostAttrs[["is_incid"]]
    if (inherits(Amatrix,"pMatrix")) {
      # ... colnams-using code removed.
    } else if ( ! is.null(is_incid)) {
      if (is_incid) is_incid <- attr(Amatrix,"is_incid") # .spaMM_spde.make.A() provides this attr. Otherwise, may be NULL, in which case ./.
      # ./. a later correct message may occur ("'is_incid' attribute missing, which suggests inefficient code in .calc_new_X_ZAC().)
    }           
    Z_ <- Z_ %*% Amatrix 
    attr(Z_,"is_incid") <- is_incid
    mostAttrs <- mostAttrs[setdiff(names(mostAttrs), c("class","is_incid", slotNames(Z_)))]
    for (st in names(mostAttrs)) attr(Z_,st) <- mostAttrs[[st]] # "is_incid", etc. 
  }
  Z_
}

# this version uses 'colnams' which is a potential bag of bugs
# and de facto does not work for composite mv at least (failed test adjacency mv)
..calc_normalized_ZA_with_subsetting <- function(Z_, normIMRF, Amatrix, L, colnams) {
  if (normIMRF) {
    Z_ <- .normalize_IMRF_ZA(Z=Z_, A=Amatrix, L=L, colnams=colnams)
  } else {
    mostAttrs <- attributes(Z_)
    is_incid <- mostAttrs[["is_incid"]]
    if (inherits(Amatrix,"pMatrix")) {
      # subsetting by rownames does not generally work on per, mutation matrices
      Amatrix <- as(as(Amatrix, "nMatrix"), "TsparseMatrix") # => ngTMatrix 
    } else if ( ! is.null(is_incid)) {
      if (is_incid) is_incid <- attr(Amatrix,"is_incid") # .spaMM_spde.make.A() provides this attr. Otherwise, may be NULL, in which case ./.
      # ./. a later correct message may occur ("'is_incid' attribute missing, which suggests inefficient code in .calc_new_X_ZAC().)
    }          
    Z_ <- Z_ %*% Amatrix[colnams,,drop=FALSE]
    attr(Z_,"is_incid") <- is_incid
    mostAttrs <- mostAttrs[setdiff(names(mostAttrs), c("class","is_incid", slotNames(Z_)))]
    for (st in names(mostAttrs)) attr(Z_,st) <- mostAttrs[[st]] # "is_incid", etc. 
  }
  Z_
}


# Called post fit: 
.calc_normalized_ZAlist <- function(Zlist, # creates ZA from Z and A, even for non-IMRF
                                    AMatrices,
                                    vec_normIMRF, 
                                    strucList) {
  if (length(Zlist) && length(AMatrices)) {
    for (char_rd in  names(Zlist)) { # 
      if ( ! is.null(Amatrix <- AMatrices[[char_rd]])) {
        Z_ <- Zlist[[char_rd]]
        # colnams <- colnames(Z_) # may be coordinates in the form "-5.3469:36.1291", or gridcode for adjacency ranefs
        rd <- as.integer(char_rd) # I cannot yet assume strucList[[char_rd]] (nor vec_normIMRF[char_rd])
        # rownams <- rownames(Amatrix)
        # `==` recycle its arguments... which may actually be useful. But no longer used  here.
        # if ( ! all(colnams == rownams)) { # Ultimately it would be nice to remove this check 
        #   # But for user-provided A at least, the rows may not be ordered as the Z cols, so names woudl be needed. 
        #   warning("(!) Possible problem with Amatrix permutation.", immediate.=TRUE)
        #   Zlist[[char_rd]] <- ..calc_normalized_ZA_with_subsetting(
        #     Z_=Z_, normIMRF=vec_normIMRF[rd], Amatrix=Amatrix, L=strucList[[rd]], colnams=colnams)
        # } else {
          ## Now the check on names is in preprocessing -> .ZxA_with_attrs();
          ## There, the ZA's were built for each submodel. Post-fit we need to 'makelong' the A matrix:
          # nblocks <- length(colnams) %/% length(rownams)
          nblocks <- ncol(Z_) %/% nrow(Amatrix)
          if (nblocks>1L) Amatrix <- .bdiag_Amatrix(Amatrix,2L)
          Zlist[[char_rd]] <- ..calc_normalized_ZA(
            Z_=Z_, normIMRF=vec_normIMRF[rd], Amatrix=Amatrix, L=strucList[[rd]])
        # }
      }
    }
  } 
  return(Zlist) ## with other attributes unchanged
}


.init_promises_spprec_compos_corrMatrix <- function(cov_info_mat) {
  blob <- attr(cov_info_mat,"blob")
  delayedAssign("Q_CHMfactor", {
    sparse_Qmat <- drop0(cov_info_mat[["matrix"]])
    Cholesky(sparse_Qmat,LDL=FALSE,perm=FALSE) # automatically saves the dCHMsimpl in the sparse_Qmat
  }, assign.env = blob ) 
  delayedAssign("kron_Y_chol_Q", {drop0(as(blob$Q_CHMfactor,"CsparseMatrix"),1e-17)}, assign.env = blob ) 
  # the drop0 levels affects "4th decimals" and computation time 
  delayedAssign("Lunique", {
    Lunique <- drop0(solve(blob$Q_CHMfactor,system="Lt", b=.sparseDiagonal(n=ncol(blob$Q_CHMfactor), shape="g")),1e-17) # ___TAG___ note that .Lunique_info_from_Q_CHM() does not back permute
    # so if we want code consistent with that fn...  but potential divergent uses of this Lunique?
    colnames(Lunique) <- rownames(Lunique) <- colnames(cov_info_mat[["matrix"]])
    attr(Lunique,"Q_CHMfactor") <- blob$Q_CHMfactor 
    attr(Lunique, "type") <- "from_Q_CHMfactor" # used in post-fit code for matrix inversion
    Lunique
  }, assign.env = blob ) 
}

.get_rownames_corrMatrix <- function(corrMatrix) {
  if (inherits(corrMatrix,"dist")) {
    corrnames <- labels(corrMatrix) ## unclear
  } else if (inherits(corrMatrix,"precision")) {
    corrnames <- rownames(corrMatrix[["matrix"]])
  } else if (inherits(corrMatrix,c("matrix","Matrix"))) {
    corrnames <- rownames(corrMatrix)
  } else stop("Unhandled class of corrMatrix object.")
  corrnames
}

.get_nestednames_ZA <- function(ZA) {
  if (is.null(ZAnames <- attr(ZA,"RHS_nesting_info")$nested_LEVELS)) ZAnames <- colnames(ZA)
  ZAnames
}



# This fn is called by .init_assign_geoinfo() directly. 
.check_rownames_corrMatrix <- function(corrMatrix, ZAnames, For, corrnames=.get_rownames_corrMatrix(corrMatrix)) { # 
  if (is.null(ZAnames)) {
    stop("NULL colnames in (a block of) the design matrix for random effects. Some mishandling of 'AMatrices'?")
  }
  if (is.null(corrnames)) {
    if (For=="corrFamily$Cf") {
      mess <- paste("corrFamily$Cf() did not provide rownames so the match with",
                    "\n   levels of the random effect (",
                    paste0(ZAnames[1L:min(5L,length(ZAnames))], collapse=" "),if(length(ZAnames)>5L){"...)"} else{")"},
                    "\n   cannot be checked. Provide names.")
      stop(mess)
    } else {
      mess <- paste("(!) corrMatrix without labels or row names: the grouping levels, in order",
                    paste0(ZAnames[1L:min(5L,length(ZAnames))], collapse=" "),if(length(ZAnames)>5L){"...,"} else{","},
                    "\n are matched in this order to rows and columns of corrMatrix, without further check.",
                    "\n This may cause later visible errors (notably, wrongly dimensioned matrices)",
                    "\n or even silent errors. See help(\"corrMatrix\") for a safer syntax.")
      warning(mess, immediate. = TRUE)
    }
  } else if (length(extraZAnames <- setdiff(ZAnames,corrnames))) { # There are ZAnames without matching corrnames. 
    if  (For=="corrFamily$Cf") {
      ## i.e. if not all ZAnames in corrnames
      mess <- paste("The row names provided by corrFamily$Cf() do not include all levels",
                    "\n   of the random effect (",
                    paste0(ZAnames[1L:min(5L,length(ZAnames))], collapse=" "),if(length(ZAnames)>5L){"...)."} else{")."},
                    "\n: check the definition of <corrFamily>$Cf() and possibly $Af().")
      stop(mess)
    } else if (inherits(corrMatrix,"precision")) {
      ## Uses a strict approach bc it's already complicated enough, and later code will again compare the names.
      stop("Some levels of the grouping variable are missing from the row names of the precision matrix\n: check corrMatrix dimensions and/or names.")
    } else {
      ## For true covariance matrix spaMM tries to accomodate: it tries to match matrices by row order,
      # but need identical dimensions for such a match.
      uZAnames <- unique(ZAnames)
      if ( length(corrnames)!=length(uZAnames)){ 
        stop("The corrMatrix does not match the levels of the grouping variable: different levels (names) and different dimensions.")
      } else { ## same dimensions, but names do not match
        message(paste0("spaMM is not able to match levels of the random effect to the names of corrMatrix,\n",
                       " and matches levels to rows of the matrix by their respective orders.\n",
                       " See help(\"corrMatrix\") for a safer syntax."))
      }
    }
  }
}

.uniqueGeo_from_ulevels <- function(unique_levels, uGeo_colnames) { # unique_levels expected to have a "colnames" attribute
  if (is.character(unique_levels)) { # seq_levelrange when it is character vec from apply(.,paste0,collapse=":"), or dataordered_unique_levels  := unique(as.character(RHS_info$factor)
                                     # RHS_info$factor was itself provided as a factor from int or factor or char... its split may have 1 or two elements
    uniqueGeo <- do.call(rbind,strsplit(unique_levels,":"))
    uniqueGeo <- as.data.frame(uniqueGeo) # 1-col or 2-col data frame
    uniqueGeo[,1] <- as.integer(uniqueGeo[,1]) # 1st col is int, optional 2nd col is char from the strsplit
    colnames(uniqueGeo) <- uGeo_colnames
    rownames(uniqueGeo) <- unique_levels
  } else { # only for seq_levelrange being integer vector
    uniqueGeo <- unique_levels
    dim(uniqueGeo) <- c(length(uniqueGeo),1L) # 1-col matrix from 'seq_levelrange' when it is integer vec.
    colnames(uniqueGeo) <- uGeo_colnames
  }
  uniqueGeo # 1-col matrix from integer vector, or 1-col data frame (int), or 2-col data frame (int, char)
}

.redef_Cf_4_nested_ranef <- function(oldCf, RHS_nesting_info) {
  if (is.null( oldCf )) { 
    stop("corrFamily API requires a 'Cf' element.")
    # i.e., no code handling RHS nesting in private case without $Cf
  } else {
    # redefine Cf
    newCf <- function(parvec) {
      cov_info_mat <- oldCf(parvec) # calling the original $Cf
      nested_Zcols <- RHS_nesting_info$nested_Zcols
      nbks <- length(nested_Zcols)
      blocks <- fullnames <- vector("list", nbks)
      blocknames <- names(nested_Zcols)
      for (blk in seq_len(nbks)) {
        indices <- nested_Zcols[[blk]]
        in_levels <- RHS_nesting_info$nested_LEVELS[indices]
        bloc <- cov_info_mat[in_levels,in_levels]
        blocks[[blk]] <- bloc
        fullnames[[blk]] <- paste0(in_levels,":",blocknames[blk])
      }
      cov_info_mat <- Matrix::bdiag(blocks) # ignores names
      rownames(cov_info_mat) <- colnames(cov_info_mat) <- .unlist(fullnames)
      cov_info_mat <- cov_info_mat[RHS_nesting_info$full_LEVELS,RHS_nesting_info$full_LEVELS] # is ZA order
      # str(cov_info_mat)
      # print(diag(cov_info_mat))
      cov_info_mat
    }
    newCf
  }
}


.init_assign_geoinfo <- function(processed, ZAlist, For=processed$For, corr_info=processed$corr_info, corr_types=corr_info$corr_types, 
                             exp_barlist, nrand=length(ZAlist), distMatrix, 
                             sparse_precision=processed$is_spprec) {
  if ( length(corr_types[ ! is.na(corr_types)])) {
    
    if (For=="HLfit") {
      ranef_string <- attr(ZAlist,"exp_ranef_strings")[ ! is.na(corr_types)][1L]
      stop(paste("Term",ranef_string,"not allowed in HLfit(). Try another fitting function such as fitme()."))
    } else {
      ## Cannot be unlisted bc cf ?unlist " Non-vector elements ... are not coerced, and so a list containing ... remains a list":
      attr(ZAlist,"exp_spatial_terms") <- exp_spatial_terms <- .findSpatial(barlist=exp_barlist, nomatch=NA, expand=TRUE) ## match ZAlist for predict()
      if (For != "fitmv") {
        #### In the preprocessing of fitmv models, 
        #     .preprocess is first called on each sub model with For="fitmv"
        #      So this block is not run
        #      Then .init_assign_geoinfo is called on processed=merged  with For="fitme", so this block is run on the merged info
        #      Finally there is a step unmerged[[mv_it]]$geo_info <- merged$geo_info[rd_in_mv]
        #      so submodel info is deduced from merged info, rather than the reverse.
        #      On the other hand the 'merged' AMatrices is a merging of A matrices for submodels
        #      Hence the A matrices have been computed while this block has not been run, as is the case also for fitme calls().
        #   => If I want to make geo_info available for AMatrix computation in fitme(), I need to move its declaration
        #      out of this function and before AMatrix computations.
        #   BUT there is no way to create a merged$geo_info 'before' AMatrix omputations for each submodel, unless the whole
        #      AMatrix code for .preprocess_MV is revised. Hmmmm.
        #   => geo_info <- processed$geo_info        not feasible here
        # 
        geo_info <- structure(vector("list",nrand), names=seq_len(nrand)) ## each element will be an environment with typical elements $distMatrix, $uniqueGeo, $nbUnique
        cov_info_mats <- vector("list",nrand)
        for (it in seq_along(corr_types)) {
          corr_type <- corr_types[[it]]
          if ( ! is.na(corr_type)) {
            if (corr_type== "corrMatrix") { 
              ## For CORREL algos .subset_corrMatrix has the effect that the correl matrix used in later computations 
              # is a subset & permutation of the user-level one according to the order of columns of ZAlist[[it]]
              #
              ## In the SPPREC case corr_info$corrMatrices[[it]] already contains a precision matrix
              # it cannot be subsetted. Nothing is done at this point(in particular, it would be pointless to modify the ZA argument locally)
              # Instead, the Z matrix will be extended by .addrightcols_Z(), then the cov_info_mat will be reordered 
              # according to the columns of this augmented Z matrix.
              # This implies that the order cov_info_mat will then differ from that of the $corrMatrices[[]]
              # and that the latter should not be used (=> _F I X M E__ remove this corrMatrices[[]] object?)
              corrMatrix <- corr_info$corrMatrices[[it]]
              corrnames <- .get_rownames_corrMatrix(corrMatrix)
              NESTED_ZAnames <- .get_nestednames_ZA(ZAlist[[it]]) # rather than simply colnames(ZAlist[[it]])
              .check_rownames_corrMatrix(corrMatrix=corrMatrix, ZAnames=NESTED_ZAnames, For="", corrnames=corrnames)
              cov_info_mat <- .subset_corrMatrix(corrMatrix=corrMatrix, ZAnames=NESTED_ZAnames, corrnames=corrnames) ## correlation or precision...
              if ( inherits(cov_info_mat,"precision")) {
                # fast test by verif3 in test pedigree (no new cols) and longer test by test-Matern-spprec.R#9 (new cols)
                # Also Gryphon... caught bug  mv  RHS in .addrightcols_Z() 
                precmat <- cov_info_mat[["matrix"]]
                precnames <- colnames(precmat) # NESTED levels
                mostAttrs <- attributes(ZAlist[[it]])
                mostAttrs <- mostAttrs[setdiff(names(mostAttrs), c("class",slotNames(ZAlist[[it]])))]
                ZAlist[[it]] <- .addrightcols_Z(Z=ZAlist[[it]], precnames) # Handles LHS and nested RHS; 
                # now we are sure that they have the same names, only the orders are uncertain, so we can test order by any( != )
                uZAnames <- unique(colnames(ZAlist[[it]]))
                if (any(precnames!=uZAnames)) cov_info_mat <- cov_info_mat[uZAnames,uZAnames] # only a permutation, using `[.precision`
                for (st in names(mostAttrs)) attr(ZAlist[[it]],st) <- mostAttrs[[st]] # "is_incid", etc. 
                # AFTER the subsetting !
                attr(cov_info_mat,"blob") <- new.env(parent=emptyenv()) 
                .init_promises_spprec_compos_corrMatrix(cov_info_mat) # "blob" environment gets promises for Q_CHMfactor, kron_Y_chol_Q, Lunique
                # : this envir was conceived for composite ranefs, but could be used beyond?
              } else attr(cov_info_mat,"blob") <- new.env(parent=emptyenv()) # "blob" without promises 
              cov_info_mats[[it]] <- cov_info_mat
            } else { ## all cases where geo_info (even empty) is needed 
              geo_info[[it]] <- new.env(parent=emptyenv())
              if ( ! is.null(distMatrix)) { ## user-provided distMatrix )
                if (is.list(distMatrix)) { ## spaMM3.0 extended syntax
                  geo_info[[it]]$distMatrix <- distMatrix[[it]] 
                } else geo_info[[it]]$distMatrix <- distMatrix 
              } # geo_info[[it]]$distMatrix may be modified below for Matern(LHS|.) => # not "else" because Matern/Cauchy must be combined with this!
              RHS_info <- attr(ZAlist[[it]], "RHS_info")
              if (corr_type == "AR1" || 
                  .get_levels_type(corr_info=corr_info, it=it, default="")=="time_series") { # includes ARp...
                if ( ! sparse_precision) {
                  # The above test implies "if not already ! sparse when ZAlist was first evaluated"
                  dataordered_unique_levels <- RHS_info$dataordered_unique_levels
                  geo_info[[it]]$uniqueGeo <- .uniqueGeo_from_ulevels(unique_levels=dataordered_unique_levels, # $uniqueGeo is here 1-col data frame (int), or 2-col data frame (int, char)
                                                                      uGeo_colnames=RHS_info$splt)
                  ### Column deletion: 
                  mostAttrs <- attributes(ZAlist[[it]])
                  mostAttrs <- mostAttrs[setdiff(names(mostAttrs), c("class",slotNames(ZAlist[[it]])))]
                  if (length(mostAttrs$namesTerm)>1L) { # AR1(multiple-column LHS|.)
                    ## In Matern(multiple-column LHS|.) we opted for not trying to remove columns.
                    ## But for AR1 this would conflict with the code from .assign_geoinfo_and_LMatrices_but_ranCoefs() 
                    # geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("distMatrix"), 
                    #                           dist_method_rd=control.dist[[char_rd]]$dist.method)
                    ## which provides a 'reduced' geo_envir$distMatrix to
                    # cov_info_mat <- corr_families[[rd]]$calc_corr_from_dist(ranFix=ranPars, char_rd=char_rd, 
                    #                                                        distmat=geo_envir$distMatrix)
                    ## so we do reduce columns here. The which() is intended to handle repeated names.
                    ZAlist[[it]] <- ZAlist[[it]][ ,which(colnames(ZAlist[[it]]) %in% 
                                                           dataordered_unique_levels), drop=FALSE]
                  } else ZAlist[[it]] <- ZAlist[[it]][ ,(dataordered_unique_levels), drop=FALSE]
                  ## subsetting dropped attributes 
                  ## "is_incid"      "namesTerm"     "LHS_levels"    "RHS_info"      "Z_levels_type"
                  for (st in names(mostAttrs)) attr(ZAlist[[it]],st) <- mostAttrs[[st]]
                  ###
                } else {
                  seq_levelrange <- .seq_levelrange(RHS_info)
                  geo_info[[it]]$uniqueGeo <- .uniqueGeo_from_ulevels(unique_levels=seq_levelrange, # $uniqueGeo is here 1-col matrix from integer vector, or 1-col data frame (int), or 2-col data frame (int, char)
                                                                      uGeo_colnames=RHS_info$splt)
                }
                # } else if (corr_type %in% c("Matern","Cauchy") ) {
                # in that case we either
                # (1) have columns in ZAlist[[it]] not reordered: no further issues. This is the current design
                #  or
                # (2) have columns in ZAlist[[it]] not reordered in .calc_Zmatrix Then
                #     the cov_info_mats[[it]] info is not yet available at preprocessing time...
                #     See further comment in .assign_geoinfo_and_LMatrices_but_ranCoefs()
              } else if (corr_type %in% c("Matern","Cauchy")) { 
                # Here I tried to set up the 'activelevels' info so that levels can be removed 
                # from corr matrices for nested random effects. But this never worked. 
                # * For non-composite ranefs, the attempt failed too (activelevels being indexed as columns of ZA but
                # conceived to remove rows of uniqueGeo that only refer to levels of the spatial variable)
                # * For composite ranefs, the problem would probably be worse, as ZAlist cols should presumably match 
                # those of the kronecker product; e.g. 54 distinct spatial positions, and LHS -> factor with 18 levels. 
                # ZAlist may have 972 cols. Again, the distinct attempts to set up 'activelevels' in that case did not work.
                # => I remove this code from version 4.1.34.
                # The corrFamily code (next) may be interesting for any future attempt.
                # even simple random-coeff terms, say (female|grp) may have inactive levels.
              } else if (corr_type== "corrFamily") { 
                # for time-series, this is here treated as the 'corrMatrix' case, rather than the AR1 case: 
                # in the SPPREC case, extention of Z; in the he CORR cases, .subset_corrFamily() instead of .subset_corrMatrix()
                corrfamily <- corr_info$corr_families[[it]]
                
                #### PBBLY not the best place for the next block of code, which could be run earlier ? (_F I X M E__?)
                if (! is.null(RHS_nesting_info <- attr(ZAlist[[it]],"RHS_nesting_info"))) { 
                  if (ncol(corrfamily$template)!=ncol(ZAlist[[it]])) {
                    corrfamily$"Cf" <- .redef_Cf_4_nested_ranef(oldCf=corrfamily$"Cf", RHS_nesting_info)
                    corrfamily$template <- corrfamily$"Cf"(corrfamily$tpar)
                    if (inherits(corrfamily$template,"CsparseMatrix")) corrfamily$sp_chk <- length(corrfamily$template@i)
                  } # else the user defined a $Cf that already has the extended dimension: Cf last example in corrFamily-design:
                }
                ####
                
                template <- corrfamily$template
                if ( ! is.null(template)) { # allow null template
                  corrnames <- rownames(template)
                  ZAnames <- colnames(ZAlist[[it]])
                  # template provide info for ZAlist; but may not be used beyond that, notably if $"Cf" exists...
                  .check_rownames_corrMatrix(corrMatrix=template, ZAnames=ZAnames, For="corrFamily$Cf", corrnames=corrnames)
                  # the template for the code is the corrMatrix code since it also handles a pr?defined matrix
                  if (processed$is_spprec) {
                    # spprec => possible extension/permutation of *ZA*
                    mostAttrs <- attributes(ZAlist[[it]])
                    mostAttrs <- mostAttrs[setdiff(names(mostAttrs), c("class",slotNames(ZAlist[[it]])))]
                    ZAlist[[it]] <- .addrightcols_Z(Z=ZAlist[[it]], corrnames, 
                                                    # warning about additional levels in the prec mat is not always pertinent;
                                                    # this inhibitsthe warning for time_series:
                                                    verbose=( ! identical(corrfamily$levels_type,"time_series")))
                    # now we are sure that they have the same names, only the orders are uncertain, so we can test order by any( != )
                    uZAnames <- unique(colnames(ZAlist[[it]])) # ! updated names
                    if (any(corrnames!=uZAnames)) .corrfamily_permute(corrfamily, perm=uZAnames)
                    for (st in names(mostAttrs)) attr(ZAlist[[it]],st) <- mostAttrs[[st]] # "is_incid", etc. 
                  } else { 
                    # sp|de corr => possible reduction/permutation of *corrFamily template*
                    .subset_corrFamily(corrfamily, ZAnames, corrnames) # calls .corrfamily_permute(), which ma redefine $"Cf" 
                  }
                }
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
