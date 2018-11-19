.get_geo_info <- function(processed, which_ranef, which="", dist.method) {
  geo_envir <- processed$geo_info[[which_ranef]] ## should be an environment, possibly empty. Check:
  if (is.null(geo_envir)) stop("The required environment was not created during preprocessing.")
  needed <- c("distMatrix"=("distMatrix" %in% which && is.null(geo_envir$distMatrix)),
              "uniqueGeo"=("uniqueGeo" %in% which && is.null(geo_envir$uniqueGeo)),
              "nbUnique"=("nbUnique" %in% which && is.null(geo_envir$nbUnique)) )
  if (any(needed)) {
    blob <- .get_dist_nested_or_not(spatial_term=attr(processed$ZAlist,"exp_spatial_terms")[[which_ranef]], 
                                data=processed$data, distMatrix=geo_envir$distMatrix, 
                                uniqueGeo=geo_envir$uniqueGeo, ## precomputed only for AR1
                                dist.method = dist.method,needed=needed)
    ## Avoid overwriting preexisting ones with possible NULL's:
    if (is.null(geo_envir$distMatrix)) geo_envir$distMatrix <- blob$distMatrix 
    if (is.null(geo_envir$uniqueGeo)) geo_envir$uniqueGeo <- blob$uniqueGeo  
    if (is.null(geo_envir$nbUnique)) geo_envir$nbUnique <- blob$nbUnique 
    if (is.null(geo_envir$coordinates)) geo_envir$coordinates <- blob$coordinates 
  }
  return(geo_envir) ## some elements may still be NULL
}

.assign_geoinfo_and_LMatrices_but_ranCoefs <- function(processed, corr_types, spatial_terms, 
                                                       ranPars, control.dist, argsfordesignL) {
  # * assigns geo_envir <- .get_geo_info(...)
  # * modifies processed$AUGI0_ZX$envir by .init_precision_info(...) 
  # * computes processed$AUGI0_ZX$envir$LMatrices except for ranCoefs (the latter being filled in HLfit_body)
  envir <- processed$AUGI0_ZX$envir
  if (processed$sparsePrecisionBOOL) .init_precision_info(processed,NULL) ## modifies processed$AUGI0_ZX$envir  
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
      if (corr_type == "adjacency") {
        adjMatrix <- processed$corr_info$adjMatrices[[rd]]
        symSVD <- attr(adjMatrix,"symSVD")
        if (is.null(symSVD) && adj_rho_is_inner_estimated) { 
          ## Can occur in direct call of HLCor, or more genrally if(list(processed)),
          #  bc symSVD added to proc1 in fitme/corrHLfit as side effect of determining rhorange
          # We tray to avoid the computation when its cost/benefit ratio is low.
          # Nevertheless if computed here, it is added to the proc envir hence computed only once per environment
          rho <- attr(ranPars,"init.HLfit")$corrPars[[char_rd]]$rho
          if (isSymmetric(adjMatrix)) {
            symSVD <- eigen(adjMatrix, symmetric=TRUE)
            svdnames <- names(symSVD)
            svdnames[svdnames=="values"] <- "d"
            svdnames[svdnames=="vectors"] <- "u"
            names(symSVD) <- svdnames
          }             
        }
        if (is.null(symSVD)) {
          cov_info_mat <- as.matrix(solve(diag(nrow(adjMatrix))-rho*(adjMatrix))) 
        } else {
          symSVD$adjd <- symSVD$d
          # fixme remove $d from symSVD to avoid later confusions ?
        }
        # symSVD may be modified below
      }  else if (corr_type =="SAR_WWt") { 
        adjMatrix <- processed$corr_info$adjMatrices[[rd]]
        UDU. <- attr(adjMatrix,"UDU.")
        if (is.null(UDU.)) {
          # cf comments on symSVD above
          cov_info_mat <- as.matrix(solve(diag(nrow(adjMatrix))-rho*(adjMatrix))) 
        } else {
          cov_info_mat <- .ZWZt(UDU.$u, 1/(1-rho*UDU.$d)) # UDU.$u %*% sweep(UDU.$u.,MARGIN=1,1/(1-rho*UDU.$d),`*`) 
        }
        cov_info_mat <- .tcrossprodCpp(cov_info_mat,NULL)
      }  else if (corr_type=="AR1" && ! processed$sparsePrecisionBOOL) {
        geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("distMatrix"), 
                                   dist.method=control.dist[[char_rd]]$dist.method)
        cov_info_mat <- .get_cP_stuff(ranPars,"ARphi",which=char_rd)^(geo_envir$distMatrix)  
        cov_info_mat[geo_envir$distMatrix==Inf] <- 0 ## should not be necess, but is.
      } else  if (corr_type %in% c("Matern","Cauchy")) {
        control_dist_rd <- control.dist[[char_rd]]
        txt <- paste(c(spatial_term[[2]][[3]])) ## the RHS of the ( . | . ) # c() to handle very long RHS
        if (length(grep("%in%",txt))) {
          stop(paste0("(!) ",corr_type,"( . | <coord> %in% <grp>) is not yet handled."))
        } 
        msd.arglist <- list(rho = rho)
        msd.arglist$`dist.method` <- control_dist_rd$`dist.method` ## may be NULL
        if (length(rho)>1L) {
          geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("uniqueGeo"), 
                                     dist.method=control_dist_rd$dist.method)
          msd.arglist$uniqueGeo <- geo_envir$uniqueGeo
          msd.arglist$`rho.mapping` <- control_dist_rd$`rho.mapping` ## may be NULL
        } else {
          geo_envir <- .get_geo_info(processed, which_ranef=rd, which=c("distMatrix"), 
                                     dist.method=control_dist_rd$dist.method)
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
      }
      if (corr_type== "corrMatrix") { ## could cover all specifications of a constant 'cov' Matrix without correlation parameter
        cov_info_mat <- processed$corr_info$cov_info_mats[[rd]] ## correlation or precision...
      } 
      ## Provide Lunique if not already available (and optionally additional stuff)
      if (processed$sparsePrecisionBOOL) {
        .init_precision_info(processed,NULL) ## modifies processed$AUGI0_ZX$envir  
        ## (1) Provide sparse_Qmat
        if (corr_type=="corrMatrix" && inherits(cov_info_mat,"precision")) {
          sparse_Qmat <- drop0(cov_info_mat[["matrix"]])
        } else if (corr_type=="adjacency" 
                   && ! adj_rho_is_inner_estimated) {
          ## Implies call from fitme_body with outer rho estim.
          sparse_Qmat <- - rho * drop0(adjMatrix)
          diag(sparse_Qmat) <- diag(sparse_Qmat)+1
          ## Cholesky gives proper LL' (think LDL')  while chol() gives L'L...
        } else if (corr_type=="AR1") { 
          types <- processed$AUGI0_ZX$envir$finertypes
          ARphi <- .get_cP_stuff(ranPars,"ARphi",which=char_rd)
          seq_n_u_h <- diff(processed$cum_n_u_h)
          names(seq_n_u_h) <- types
          ## the dim of Q *HAS* to match that of ZA
          AR1_block_n_u_h_s <- attr(processed$ZAlist[[rd]],"AR1_block_n_u_h_s") # block size[S for single AR1 term if nested]
          t_chol_Q <- sapply(AR1_block_n_u_h_s, .calc_AR1_t_chol_Q_block,ARphi=ARphi) ## L block(s) of given size(s)
          ## bdiag required if nested AR1 term:
          Lunique <- Matrix::bdiag(lapply(t_chol_Q, solve )) # sparse triangular solve 
          attr(Lunique, "type") <- "from_AR1_specific_code"
          sparse_Qmat <- Matrix::bdiag(lapply(t_chol_Q, crossprod ))
          processed$AUGI0_ZX$envir$precisionFactorList[[rd]] <- list(chol_Q=Matrix::bdiag(lapply(t_chol_Q, t )), # Linv
                                                                     Qmat=sparse_Qmat)
        } else if (corr_type %in% c("Matern","Cauchy")) {
          ## at this point cov_info_mat is a dist object !
          cov_info_mat <- as.matrix(cov_info_mat)
          diag(cov_info_mat) <- 1
          ## If the order of columns of ZAlist[[it]] had been permuted by .calc_Zmatrix(), we ould need something like:
          # ZAnames <- colnames(processed$ZAlist[[rd]])
          # cov_info_mat <- cov_info_mat[ZAnames,ZAnames]
          # but this would require controlling the names of Z to be those of cov_info_mat (spMMFactorList_locfn() does not care for that)
          sparse_Qmat <- (as(solve(cov_info_mat),"sparseMatrix")) # even though it's not mathematically sparse
        } else {stop("Some error occurred (inner estimation of adjacency rho with requested sparse precision ?)")}     
        ## (2) Builds from sparse_Qmat if not already available
        if (corr_type != "AR1") { ## General code that is correct for AR1 too
          ## I have tried to use Q_CHMfactor uniformly but Cholesky may poorly handle nearly-singular correlation matrices
          ## solve(sparse_Qmat) gives the correlation matrix
          Q_CHMfactor <- Cholesky(sparse_Qmat,LDL=FALSE,perm=FALSE) ## called for each corrPars
          ## fitme sparse_precision has an incomplete symSVD=> corr mattrix not computed, 
          ##    and try(<mat_sqrt_fn>(symSVD=symSVD)) fails. Instead use code always valid:
          Lunique <- solve(Q_CHMfactor,system="Lt") #  L_Q^{-\top}=LMatrix_correlation 
          ## Lunique is a dtCMatrix; it is for single correlated effect, and still used in HLfit_body; 
          ## whether it is used or not in MME_method (sXaug_...), a lot of other code still expects it
          attr(Lunique,"Q_CHMfactor") <- Q_CHMfactor ## FIXME remove ? limited use that saves only optional sparse triangular solves
          attr(Lunique, "type") <- "from_Q_CHMfactor"
          processed$AUGI0_ZX$envir$precisionFactorList[[rd]] <- list(chol_Q=as(Q_CHMfactor, "sparseMatrix"), # Linv
                                                                     Qmat=sparse_Qmat)
        }
      } else {
        # this is called through HLCor, hence there must be a dense correlation structure
        ## densePrecision cases with outer rho estimation, and some residual inner rho estimation cases
        if ( ! is.null(symSVD)) {
          symSVD$d <- 1/(1-rho*symSVD$adjd) ## from adjMatrix to correlation matrix
          # outer optim -> LMatrix recomputed from this for each rho  
          Lunique <- do.call(.spaMM.data$options$mat_sqrt_fn,list(symSVD=symSVD)) ## using $d not $adjd
          if (inherits(Lunique,"try-error")) { ## mat_sqrt_fn has try()'s and can return "try-error"'s
            print("correlation parameter was rho=:",rho,quote=FALSE) ## makes sense if mat_sqrt_fn already issued some warning
            stop()
          }
          if ( adj_rho_is_inner_estimated ) { # contrived way of construction Lunique with the correct attributes.
            Lunique[] <- attr(Lunique,"symsvd")$u   ## "[] <- " keeps attributes... except for Matrix...
          }
        } else { ## cov_info_mat must exist
          if (processed$HL[1L]=="SEM") argsfordesignL$try.chol <- FALSE
          if (inherits(cov_info_mat,"dist")) {
            cov_info_mat <- as.matrix(cov_info_mat)
            diag(cov_info_mat) <- 1L ## IF diag missing in input corrMatrix THEN assume a correlation matrix
          } ## else full matrix may be a COV matrix with non-unit diag
          Lunique <- do.call(.spaMM.data$options$mat_sqrt_fn,c(list(m=cov_info_mat),argsfordesignL)) ## mat_sqrt_fn has try()'s and can return "try-error"'s
          if (inherits(Lunique,"try-error")) { 
            print("correlation parameters were:",quote=FALSE) ## makes sense if mat_sqrt_fn already issued some warning
            print(unlist(argsfordesignL))    
            stop()
          }
        }
        attr(Lunique,"msd.arglist") <- msd.arglist ## NULL except for Matern, Cauchy
      }
      attr(Lunique,"corr.model") <- corr_type
      attr(Lunique,"ranefs") <- paste(c(spatial_term)) ## essentiel pour la construction de ZAL! ## paste(c()) handles very long RHS
      envir$LMatrices[[rd]] <- Lunique
      attr(envir$LMatrices,"is_given_by")[rd] <- "AUGI0_ZX$envir" 
    }
    if ( ! is.null(symSVD$adjd)) { 
      ## we should not replace attr(...$adjMatrices[[rd]], "symSVD") which has a given $d distinct from the new symSVD$d
      ## Moreover, currently we use only he new symSVD's $adjd
      attr(processed$corr_info$adjMatrices[[rd]],"adjd") <- symSVD$adjd ## 
    } 
  }
}
