.get_geo_info <- function(processed, which_ranef, which="", dist.method) {
  geo_envir <- processed$geo_info[[which_ranef]] ## should be an environment, possibly empty. Check:
  if (is.null(geo_envir)) stop("The required environment was not created during preprocessing.")
  needed <- c(distMatrix=("distMatrix" %in% which && is.null(geo_envir$distMatrix)),
              uniqueGeo=("uniqueGeo" %in% which && is.null(geo_envir$uniqueGeo)),
              nbUnique=("nbUnique" %in% which && is.null(geo_envir$nbUnique)) )
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
                                                       ranPars, control.dist, trueCorrpars, argsfordesignL) {
  # * assigns geo_envir <- .get_geo_info(...)
  # * modifies processed$AUGI0_ZX$envir by .init_precision_info(...) 
  # * computes processed$AUGI0_ZX$envir$LMatrices except for ranCoefs (the latter being filled in HLfit_body)
  envir <- processed$AUGI0_ZX$envir
  rho <- ranPars$rho ## more general code needed for spaMM 3.0 !! 
  if (processed$sparsePrecisionBOOL) .init_precision_info(processed,NULL) ## modifies processed$AUGI0_ZX$envir  
  for (it in seq_along(corr_types)) {
    corr_type <- corr_types[it]
    if (!is.na(corr_type)) {
      spatial_term <- spatial_terms[[it]] 
      symSVD <- NULL
      adj_rho_is_inner_estimated <- (corr_type=="adjacency"
                                     #                               && ! is.null(attr(ranPars,"type")) ## through corrHLfit (or fitme !?) or some direct HLCor call
                                     && identical(attr(ranPars,"type")$rho,"var") ## can occur in direct call of HLCor 
      )    
      if (corr_type == "adjacency") {
        adjMatrix <- processed$corr_info$adjMatrices[[it]]
        symSVD <- attr(adjMatrix,"symSVD")
        if (is.null(symSVD) && adj_rho_is_inner_estimated) { 
          ## Can occur in direct call of HLCor, or more genrally if(list(processed)),
          #  bc symSVD added to proc1 in fitme/corrHLfit as side effect of determining rhorange
          # We tray to avoid the computation when its cost/benefit ratio is low.
          # Nevertheless if computed here, it is added to the proc envir hence computed only once per environment
          if (isSymmetric(adjMatrix)) {
            symSVD <- sym_eigen(adjMatrix)
            attr(processed$corr_info$adjMatrices[[it]],"symSVD") <- symSVD  
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
        adjMatrix <- processed$corr_info$adjMatrices[[it]]
        UDU. <- attr(adjMatrix,"UDU.")
        if (is.null(UDU.)) {
          # cf comments on symSVD above
          cov_info_mat <- as.matrix(solve(diag(nrow(adjMatrix))-rho*(adjMatrix))) 
        } else {
          cov_info_mat <- UDU.$u %*% sweep(UDU.$u.,MARGIN=1,1/(1-rho*UDU.$d),`*`) 
        }
        cov_info_mat <- .tcrossprodCpp(cov_info_mat,NULL)
      }  else if (corr_type=="AR1" && ! processed$sparsePrecisionBOOL) {
        geo_envir <- .get_geo_info(processed, which_ranef=it, which=c("distMatrix"), 
                                   dist.method=control.dist$dist.method)
        cov_info_mat <- trueCorrpars$ARphi^(geo_envir$distMatrix)  
        cov_info_mat[geo_envir$distMatrix==Inf] <- 0 ## should not be necess, but is.
      } else  if (corr_type =="Matern") {
        txt <- paste(c(spatial_term[[2]][[3]])) ## the RHS of the ( . | . ) # c() to handle very long RHS
        if (length(grep("%in%",txt))>0) {
          stop("(!) Matern( . | <coord> %in% <grp>) is not yet handled.")
        } 
        msd.arglist <- list(rho = rho)
        msd.arglist$`dist.method` <- control.dist$`dist.method` ## may be NULL
        if (length(rho)>1L) {
          geo_envir <- .get_geo_info(processed, which_ranef=it, which=c("uniqueGeo"), 
                                     dist.method=control.dist$dist.method)
          msd.arglist <- c(msd.arglist,list(uniqueGeo=geo_envir$uniqueGeo))
          msd.arglist$`rho.mapping` <- control.dist$`rho.mapping` ## may be NULL
        } else {
          geo_envir <- .get_geo_info(processed, which_ranef=it, which=c("distMatrix"), 
                                     dist.method=control.dist$dist.method)
          msd.arglist$distMatrix <- geo_envir$distMatrix   
        }
        cov_info_mat <- do.call("make_scaled_dist",msd.arglist)
        ## at this point if a single location, dist_mat should be dist(0) and make_scaled_dist was modified to that effect
        if ( nrow(cov_info_mat)>1 ) { ## >1 locations
          norho <- trueCorrpars 
          norho$rho <- NULL ## because the MaternCorr input will be an already scaled distance 'cov_info_mat'
          cov_info_mat <- do.call(MaternCorr,args=c(norho,list(d=cov_info_mat)))        
        } 
      }
      if (corr_type== "corrMatrix") { ## could cover all specifications of a constant 'cov' Matrix without correlation parameter
        cov_info_mat <- processed$corr_info$cov_info_mats[[it]] ## correlation or precision...
      } 
      if (processed$verbose["trace"] && length(trueCorrpars)>0) print(unlist(trueCorrpars))
      ## call designL.from.Corr if Lunique not available
      ## Lunique can be available here either bc the user provided it explictly (will not occur in public usage)
      ##   or if we add an explicit calculation above
      #sparse_Qmat <- NULL
      ## Provide Lunique (and optionally additional stuff)
      #  The following defines two exclsive types of ranef types. In spaMM 3.0 it would be nice to fill some of the gaps 
      if (processed$sparsePrecisionBOOL) {
        .init_precision_info(processed,NULL) ## modifies processed$AUGI0_ZX$envir  
        ## (1) Provide sparse_Qmat
        if (corr_type=="corrMatrix" && inherits(cov_info_mat,"precision")) {
          sparse_Qmat <- Matrix::drop0(cov_info_mat[["matrix"]])
        } else if (corr_type=="adjacency" 
                   && ! adj_rho_is_inner_estimated) {
          ## Implies call from fitme_body with outer rho estim.
          sparse_Qmat <- - rho * Matrix::drop0(adjMatrix)
          diag(sparse_Qmat) <- diag(sparse_Qmat)+1
          ################################L_Q <- .designL.from.Qmat(Qmat) ## solve(tcrossprod(LMatrix)) = Qmat or
          ## Cholesky gives proper LL' (think LDL')  while chol() gives L'L...
        } else if (corr_type=="AR1") { 
          types <- processed$AUGI0_ZX$envir$finertypes
          ARphi <- ranPars$ARphi
          seq_n_u_h <- diff(processed$cum_n_u_h)
          names(seq_n_u_h) <- types
          ## the dim of Q *HAS* to match that of ZA
          AR1_block_n_u_h_s <- attr(processed$ZAlist[[which(types=="AR1")]],"AR1_block_n_u_h_s")
          t_chol_Q <- sapply(AR1_block_n_u_h_s, .calc_AR1_t_chol_Q_block,ARphi=ARphi)
          t_chol_Q <- Matrix::bdiag(t_chol_Q)
          ## bc we need the CHMfactor per se, we do not directly provide chol_Q. This is not elegant.
          sparse_Qmat <- crossprod(t_chol_Q) 
          ## print(solve(sparse_Qmat)) shows the correlation matrix
          ## as(Matrix::Cholesky(sparse_Qmat,perm=FALSE,LDL=FALSE),"sparseMatrix") gives back Linv...
          ## older code used: 
          #Lunique <- as.matrix(solve(tLinv/sqrt(1-ARphi^2))) ## corrmat is tcrossprod(Lunique): we keep tcrossprod but L' is tri.sup ! 
          #attr(Lunique,"type") <- "cholR_RRt" ## not equivalent to base::chol() which whould produce cholR_tRR 
        }      
        ## (2) Builds from sparse_Qmat
        #if ( ! is.null(sparse_Qmat)) {
        #processed$AUGI0_ZX$envir$Qmat <- sparse_Qmat ## fixme do we need this copy ?
        Q_CHMfactor <- Cholesky(sparse_Qmat,LDL=FALSE,perm=FALSE) ## called for each corrPars
        processed$AUGI0_ZX$envir$precisionFactorList[[it]] <- list(Qmat=sparse_Qmat,
                                                                   chol_Q=as(Q_CHMfactor, "sparseMatrix"))
        #}
        ## fitme sparse_precision has an incomplete symSVD=> corr mattrix not computed, 
        ##    and try(designL.from.Corr(symSVD=symSVD)) fails. Instead use code always valid:
        Lunique <- solve(Q_CHMfactor,system="Lt") #  L_Q^{-\top}=LMatrix_correlation 
        attr(Lunique, "type") <- "from_Q_CHMfactor"
        attr(Lunique,"Q_CHMfactor") <- Q_CHMfactor ## FIXME store directly the Q_CHMfactor rather than Lunique ?
        # Lunique is a dtCMatrix; it is for single correlated effect, and still used in HLfit_body; 
        ## whether it is used or not in MME_method (sXaug_...), a lot of other code still expects it
        
        ## failed concept:
        # if (types[it]=="AR1") {
        #   colnames(Lunique) <- rownames(Lunique) <- seq(levelrange[1L],levelrange[2L])
        #   ## then the non-square LMatrix: remove rows of L that won't match cols of ZA
        #   Lunique <- Lunique[colnames(processed$ZAlist[[it]]),,drop=FALSE]
        #  
        #  # BUT such an Lunique is more like an AMatrix and we should have created ZA before computing cum_n_u_h
        #  All the code (fit_as_sparse precision included) uses cum_n_u_h
        # 
        # }
      } else {
        # this is HLCor, hence there must be a dense correlation structure
        ## densePrecision cases with outer rho estimation, and some residual inner rho estimation cases
        if ( ! is.null(symSVD)) {
          symSVD$d <- 1/(1-rho*symSVD$adjd) ## from adjMatrix to correlation matrix
          # outer optim -> LMatrix recomputed from this for each rho  
          Lunique <- try(designL.from.Corr(symSVD=symSVD)) ## usin $d not $adjd
          # Lunique is num, with attr(*, "type")= chr "symsvd" ....
          if ( adj_rho_is_inner_estimated ) { # contrived way of construction Lunique with the correct attributes.
            Lunique[] <- attr(Lunique,"symsvd")$u   ## "[] <- " keeps attributes... except for Matrix...
          }
        } else { ## cov_info_mat must exist
          if (processed$HL[1L]=="SEM") argsfordesignL$try.chol <- FALSE
          if (inherits(cov_info_mat,"dist")) {
            cov_info_mat <- as.matrix(cov_info_mat)
            diag(cov_info_mat) <- 1L ## IF diag missing in input corrMatrix THEN assume a correlation matrix
          } ## else full matrix may be a COV matrix with non-unit diag
          Lunique <- try(do.call(designL.from.Corr,c(list(m=cov_info_mat),argsfordesignL)))
        }
        if (inherits(Lunique,"try-error")) { 
          print("correlation parameters were:",quote=FALSE) ## makes sense if designL.from.Corr already issued some warning
          print(unlist(trueCorrpars))    
          stop()
        }
        if (corr_type=="Matern") attr(Lunique,"msd.arglist") <- msd.arglist
      }
      attr(Lunique,"corr.model") <- corr_type
      attr(Lunique,"ranefs") <- paste(c(spatial_term)) ## essentiel pour la construction de ZAL! ## paste(c()) handles very long RHS
      envir$LMatrices[[it]] <- Lunique
      attr(envir$LMatrices,"is_given_by")[it] <- "AUGI0_ZX$envir" 
    }
  }
  envir$adj_symSVD <- symSVD ## may be NULL ## not yet spaMM 3.0 compatible; (like adjMatrix)
}
