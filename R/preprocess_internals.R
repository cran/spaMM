# there is Matrix::nnzero() that may be useful in more general contexts but which is slower.
.nonzeros <- function(spm) {
  if (inherits(spm,"ddiMatrix") && spm@diag=="U") {
    return(ncol(spm))
  } else if (inherits(spm,"sparseMatrix")) {
    nz <- length(spm@x)
    if ( methods::.hasSlot(spm, "diag")) nz <- nz+ncol(spm)
    return(nz) 
  } else return(sum(spm !=0)) ## AMatrices reaches here
}

.calc_QR_diagnosis <- function(corr_info, exp_ranef_types, ZAlist) {
  is_cF <- corr_info$corr_types=="corrFamily"
  #
  possiblydense <- exp_ranef_types %in% c("adjacency", "IMRF", "Matern","Cauchy", "corrMatrix", "AR1",
                                          "diallel", "MaternIMRFa", "ARMA", "ARp" # ad-hoc code for some corrFamilies...
  )
  for (rd in which(is_cF)) possiblydense[rd] <- corr_info$corr_families[[rd]]$possiblyDenseCorr 
  #
  actuallysparse <- (
    possiblydense & (
      grepl("%in%", attr(ZAlist, "exp_ranef_strings"), fixed=TRUE) | # nested models otherwise with possibly dense *corr*
        sapply(ZAlist, .is01col) # This should detect <corr_type>( <0/1 LHS> | RHS)
    )
  ) 
  stillpossiblydense <- possiblydense & ( ! actuallysparse)
  #
  # sparseprecs <- exp_ranef_types %in% c("adjacency", "IMRF", "AR1", "MaternIMRFa", "ARp")
  # for (rd in which(is_cF)) sparseprecs[rd] <- corr_info$corr_families[[rd]]$sparsePrec  
  
  list(stillpossiblydense=stillpossiblydense,
       actuallysparse=actuallysparse)
}

.provide_QR_diagnosis <- function(corr_info, ZAlist, exp_ranef_types=attr(ZAlist, "exp_ranef_types")) {
  if (is.null(corr_info$QR_diagnosis)) {
    corr_info$QR_diagnosis <- .calc_QR_diagnosis(corr_info=corr_info,
                                                 ZAlist=ZAlist, 
                                                 exp_ranef_types=exp_ranef_types)
  }
  return(corr_info$QR_diagnosis)
}


## Fast assessment of relative SPPREC vs CORREL cost, by comparing the denseness ZAX(corrlist) products "noAR" and "ZL", 
##   where low noAR denseness relative to ZL favors SPPREC. The corrlist elements are NOT actual L model matrices.  
## For *corrMatrix* case, the corrlist element is first set to a drop0(<precmat>) for noAR denseness computation,
##   then it is set to a dummy full triangular matrix for ZL denseness computation.
## For *corrFamily* case, this is similar to corrMatrix case, except as noted in the relevant block.
## For *Matern, Cauchy*, full dense dummy factor in noAR, triangular dummy factor in ZL (hence disfavoring SPPREC)
## For *ranCoefs* (composite or not -- note that $is_composite is not yet available to this fn nor elsewhere for simple syntaxes),
##   we further perform kronecker products (distinct for noAR and ZL) with LHS a full triangular matrix 
##     (except for NULL-for-noAR RHS, which are left NULL-for-noAR). 
## For (.|.)...) the elements for noAR and ZL are left NULL. 
## For other ranef types (AR1, IMRF) the elements for noAR are left NULL,  Dummy dense triangular matrices are build for ZL (if the preclist was NULL)
## The unexpected time-consuming step may then be compute_ZAL on potentially huge dummy dense triangular matrices. (_F I X M E__)

# **************** This function is no longer called by default ******************** 
.calc_G_diagnosis_fast <- function(corr_info, ZAlist) { # OLD method
  preclist <- corr_info$adjMatrices
  corrlist <- corr_info$corrMatrices
  corrFamlist <- corr_info$corr_families
  # Don't try to drop0(corrMatrices) here before evaluating various densenesses: these densenesses would not hold for the
  # corr_info$cov_info_mats used to compute QMat. 
  # This is why the drop0() must happen earlier, in .assign_cov_matrices__from_covStruct()
  nrand <- length(ZAlist)
  exp_ranef_types <- attr(ZAlist, "exp_ranef_types")
  which_nested_RHS <- grep("%in%", names(attr(ZAlist,"exp_ranef_terms")))
  vec_nc_ori <- integer(nrand) # two loops are run here so we need this but not in the alternative function.
  for (rd in seq_len(nrand)) {
    vec_nc_ori[rd] <- ncol(ZAlist[[rd]]) # save this info before any .addRightcols()
    # ncol(ZA) LHS-induced times RHS-induced dimension* say LHS-dim * spatial * RHS-nesting
    LHS_nlev <- length(attr(ZAlist[[rd]],"namesTerm"))
    type_rd <- exp_ranef_types[rd]
    # preclist, corrlist should only have dim of the RHS of a potential kronecker product i.e. spatial * RHS-nesting
    if (type_rd=="corrMatrix" ) {
      # If we specified "adjMatrix", we suspect it's better (but spprec is not automatic); this block 
      # is not executed and corrlist[[rd]] is (locally) NULL
      # Alternatively, we may have specified a precision matrix by "corrMatrix"... 
      ## cov_info_mats elements may be correlation matrices, or they may be lists...
      if (is.list(corrlist[[rd]])) { 
        preclist[[rd]] <- corrlist[[rd]]$matrix
        corrlist[rd] <- list(NULL) # no non-trivial alternative available? 
        # It will be replaced to compute the denseness of ZL. Anyway, this NULL means that
        # the noAR denseness will be low, favoring spprec, which makes sense we provided a precision matrix.
      } else {
        template <- corrlist[[rd]]
        if (inherits(template,"dist")) {
          template <- proxy::as.matrix(template, diag=1)
        }  
        if (got_chol <- (is.matrix(template) || is(template,"Matrix"))) { # testing dim() is not appropriate bc spaMM has a dim.precision method
          cholcorr <- try(chol(template)) # base::chol or Matrix::chol 
          if (inherits(cholcorr,"try-error")) stop(cli::format_error("A correlation matrix is (nearly) singular. Check the correlation model and/or see {.topic [sparse_precision](spaMM::sparse_precision)}.")) 
          preclist[[rd]] <- drop0(chol2inv(cholcorr), tol = .Machine$double.eps)
          ZAlist[[rd]] <- .addrightcols_Z(Z=ZAlist[[rd]], colnames(template), verbose=FALSE)
          # drop0(cholcorr...) might be useful here if also performed in the relevant .calc_Lunique_for_correl_algos() -> mat_sqrt() -> .wrap_Ltri_t_chol() call...
          corrlist[[rd]] <- drop0(cholcorr,tol=0) # more logical than (previously)  preclist[[rd]], apart from speed 
          # (could
          # (1) corrlist[[rd]] <- drop0(cholcorr,tol=0)
          # (2) preclist[[rd]] <- drop0(chol2inv(cholcorr), tol = .Machine$double.eps) 
          # be more efficient ?)
          # cf alternative fn: corrlist[[rd]] <- .fill_denseL(nc_ori, nc=ncol(ZAlist[[rd]]), LHS_nlev, existing=drop0(cholcorr,tol=0)) 
          ## I previously assumed that "the chol factor of a relatedness matrix appears" 
          # "as sparse as its precmat" but this is only roughly true if one compares only drop0() values.
          # The Gryphon example is one where denseness is low for chol2inv(cholcorr) and for drop0(cholcorr)
          # but not for cholcorr itself. So to assess density of noAR, 
          # one may use here either preclist[[rd]] or drop0(cholcorr,tol=0):
          # In that Gryphon case the user-level input corrMatrix is a dense matrix with 0's and
          # corr_info$corrMatrices[[rd]] is a dsC matrix already without drop0(.,tol=0)-able elements.
          #
          # To assess ZL density, ## : corrlist[[rd]] will be later set to a dummy fully-dense 
          # toeplitz lower triangle accounting for any .addrightcols_Z()-ed empty cols. 
        } 
      }
    } else if (type_rd == c("corrFamily") ) {
      if ( ! is.null(corrFamlist[[rd]][["f"]])) { # otherwise structure is ignored
        template <- corrFamlist[[rd]]$template
        if (inherits(template,"precision")) {
          preclist[[rd]] <- template$matrix
          corrlist[rd] <- list(NULL)   ## consistently with s=corrMatrix code, explained above
        } else {
          cholcorr <- try(chol(template)) # base::chol or Matrix::chol 
          if (inherits(cholcorr,"try-error")) stop("The matrix template does not look like a symmetric positive definite matrix. Check the 'tpar' argument..") 
          preclist[[rd]] <- drop0(chol2inv(cholcorr), tol = .Machine$double.eps)
          ZAlist[[rd]] <- .addrightcols_Z(Z=ZAlist[[rd]], colnames(template), verbose=FALSE)
          corrlist[[rd]] <- preclist[[rd]]  ## consistently with s=corrMatrix code, explained above
        }
      }
    } else if (type_rd %in% c("Matern","Cauchy")  ) {
      nc <- vec_nc_ori[rd]/LHS_nlev # assuming the later kronecker  product 
      if ( rd %in% which_nested_RHS ) { # more memory efficient than alternative code
        RHS_nesting_info <- attr(ZAlist[[rd]],"RHS_nesting_info")
        blcks <- lapply(RHS_nesting_info$nested_Zcols, function(v) {
          nc <- length(v) 
          matrix(TRUE,ncol=nc,nrow=nc)
        })
        preclist[[rd]] <- Matrix::bdiag(blcks) # slow in nested-Matern example
        blcks <- lapply(blcks, lower.tri, diag=TRUE)
        corrlist[[rd]] <- Matrix::bdiag(blcks)
        # Add names for for check in .compute_ZAXlist()
        # The precise order of names here may not be essential, but this could be a template for other algorithms where the oreder would matter.
        perm_ZAnames <- RHS_nesting_info$full_LEVELS[.unlist(RHS_nesting_info$nested_Zcols)] # must be (permuted) ZAnames
        rownames(corrlist[[rd]]) <- colnames(corrlist[[rd]]) <- perm_ZAnames
        corrlist[[rd]] <- corrlist[[rd]][RHS_nesting_info$full_LEVELS,RHS_nesting_info$full_LEVELS] # now in ZA order
      } else {
        preclist[[rd]] <- matrix(TRUE,ncol=nc,nrow=nc) # as(allTRUE,"lgCMatrix") #new("lgCMatrix",i=rep(c(0L,seq_len(nc-1L)),nc),p=c(0L,seq_len(nc))*nc,x=rep(TRUE,nc^2),Dim=c(nc,nc)) 
        corrlist[[rd]] <- lower.tri(preclist[[rd]],diag = TRUE) # logi
      } 
    }
    if (LHS_nlev>1L ) { # before .compute_ZAL()
      if ( ! is.null(corrlist[[rd]])) {
        rownames_kron_Y <- rownames(corrlist[[rd]])
        corrlist[[rd]] <- kronecker(lower.tri(matrix(1,ncol=LHS_nlev,nrow=LHS_nlev),diag=TRUE), corrlist[[rd]]) # make.dimnames does not seem to work with Matrix'es
        if (! is.null(attr(ZAlist[[rd]],"RHS_nesting_info"))) {
          rownames(corrlist[[rd]]) <- rep(rownames_kron_Y,LHS_nlev)
        }
      }
      # It is was NULL it remains NULL, which looks OK for noAR
    } 
  }
  # suppressMessages here and below as .provide_G_diagnosis() is not the right context for messages.
  noAR <- suppressMessages( .compute_ZAL(XMatrix=corrlist,ZAlist,as_matrix = FALSE) )# without the implied corrlist cost of terms with sparse precision structures 
  
  
  ## next we fill more matrices in corrlist for ZL denseness.
  for (rd in seq_len(nrand)) {
    ## much better in spprec: adjacency-long and some in test AR1 (long):
    ## fitar1 <- corrHLfit(obs ~ 1+AR1(1|age),family=poisson(),data=fake,verbose=c(TRACE=TRUE)) 
    ## fit_ar1nested <- ... also
    ## correlation algos can still be selected and are better (ohio small and many scotlip tests)
    ## SPPREC and CORREL roughly as fast for fitNF in test-devel-predVar-AR1 (test-determine_spprec.R version is different)
    updated <- FALSE
    
    nc <- ncol(ZAlist[[rd]])
    LHS_nlev <- length(attr(ZAlist[[rd]],"namesTerm"))
    if (type_rd %in% c("AR1", "IMRF", "adjacency") ) { # If the user provided e.g. a huge adjmatrix it's risky to compute a huge inverse
      if (LHS_nlev>1L) nc <- nc %/% LHS_nlev 
      corrlist[[rd]] <- lower.tri(matrix(TRUE,ncol=nc,nrow=nc),diag = TRUE) # template matrix created in faster way than diag()
      updated <- TRUE
    } else if (type_rd %in% c("corrMatrix","corrFamily")) {
      corrlist[[rd]] <- .fill_denseL(nc_ori=vec_nc_ori[rd], nc=nc, LHS_nlev=LHS_nlev, existing=corrlist[[rd]])
      updated <- TRUE
    } else  if ( is.null(corrlist[[rd]]) &&  ! is.null(preclist[[rd]])) { 
      # presumably for the never used .provide_G_diagnosis(., fast=FALSE), case;
      # not clear whether this can occur for fast=TRUE. But then computing the mat_sqrt, next solve(t(.)) may be inefficient
      warning("Possibly inefficient code in .provide_G_diagnosis(). Please contact the maintainer.") # _F I X M E__
      tcrossfac_adj <- mat_sqrt(preclist[[rd]]) # tcrossfac hence t(solve()) is the tcrossfac of the corr mat (<=> Lunique) which is the following backsolve
      if ( attr(tcrossfac_adj,"type")=="cholL_LLt") {
        if (inherits(tcrossfac_adj,"dtCMatrix")) {
          corrlist[[rd]] <- solve(t(tcrossfac_adj)) 
        } else corrlist[[rd]] <- .backsolve(tcrossfac_adj,upper.tri = FALSE, transpose=TRUE)
      } else corrlist[[rd]] <- t(solve(tcrossfac_adj)) # quick patch when another facto used.  
      updated <- TRUE
    }
    if (updated && LHS_nlev>1L ) {
      corrlist[[rd]] <- kronecker(lower.tri(matrix(1,ncol=LHS_nlev,nrow=LHS_nlev),diag=TRUE), corrlist[[rd]])
    } 
  }
  
  ZL <- suppressMessages( .compute_ZAL(XMatrix=corrlist,ZAlist,as_matrix = FALSE) )  ## > qq s for large ZA 
  
  if (is.logical(ZL)) {# logical corr by identity Z gives logi ZL, not handled by .crossprod
    cross_ZL <- crossprod(ZL)
  } else cross_ZL <- .crossprod(ZL) ## forces a call to forceSymmetric => result is Matrix either dsy or sparse.
  corr_info$G_diagnosis <- list(rel_denseness_noAR=.calc_denseness(noAR, relative=TRUE),
                                rel_denseness_cross_ZL=.calc_denseness(cross_ZL, relative=TRUE), # crossprod ideally dsC except if ZL is really dense, 
                                dimZL=dim(ZL),
                                crossZL_is_dsy=inherits(cross_ZL,"dsyMatrix"),
                                fast=TRUE)
  # if there are only "(.|.)", we compare noAR to crossprod(ZL) =crossprod(noAR) and the latter may be less dense!
  # corr_info is an environment so no need for return value.
}

.fill_denseL <- function(nc_ori, nc, LHS_nlev, existing) {
  if (LHS_nlev>1L) {
    nc <- nc %/% LHS_nlev 
    nc_ori <- nc_ori %/% LHS_nlev 
  } # so that it will return a block of the final matrix
  if (nc > nc_ori) { # nc after .addRightcols()
    toeplitz(
      c(rep(FALSE,nc-nc_ori), rep(TRUE,nc_ori)),
      rep(FALSE,nc)
    )
  } else if (inherits(existing,"dtCMatrix")) {
    existing # from corrlist[[rd]] <- cholcorr. Leave it as is.
  } else lower.tri(matrix(TRUE,ncol=nc,nrow=nc),diag = TRUE) # template matrix created in faster way than diag()
  
}

.calc_G_diagnosis_new <- function(corr_info, ZAlist) { 
  preclist <- corr_info$adjMatrices
  corrlist <- corr_info$corrMatrices
  corrFamlist <- corr_info$corr_families
  # Don't try to drop0(corrMatrices) here before evaluating various densenesses: these densenesses would not hold for the
  # corr_info$cov_info_mats used to compute QMat. 
  # This is why the drop0() must happen earlier, in .assign_cov_matrices__from_covStruct()
  nrand <- length(ZAlist)
  exp_ranef_types <- attr(ZAlist, "exp_ranef_types")
  which_nested_RHS <- grep("%in%", names(attr(ZAlist,"exp_ranef_terms")))
  for (rd in seq_len(nrand)) {
    nc_ori <- ncol(ZAlist[[rd]]) # save this info before any .addRightcols()
    LHS_nlev <- length(attr(ZAlist[[rd]],"namesTerm"))
    type_rd <- exp_ranef_types[rd]
    # ncol(ZA) LHS-induced times RHS-induced dimension* say LHS-dim * spatial * RHS-nesting
    # preclist, corrlist should only have dim of the RHS of a potential kronecker product i.e. spatial * RHS-nesting
    LHS_done <- FALSE
    if (type_rd=="corrMatrix" ) {
      # If we specified "adjMatrix", we suspect it's better (but spprec is not automatic); this block 
      # is not executed and corrlist[[rd]] is (locally) NULL
      # Alternatively, we may have specified a precision matrix by "corrMatrix"... 
      ## cov_info_mats elements may be correlation matrices, or they may be lists...
      if (is.list(corrlist[[rd]])) { 
        preclist[[rd]] <- corrlist[[rd]]$matrix
        corrlist[[rd]] <- .fill_denseL(nc_ori, nc=nc_ori, LHS_nlev, existing=NULL) 
      } else {
        template <- corrlist[[rd]]
        if (inherits(template,"dist")) {
          template <- proxy::as.matrix(template, diag=1)
        }  
        if (got_chol <- (is.matrix(template) || is(template,"Matrix"))) { # testing dim() is not appropriate bc spaMM has a dim.precision method
          cholcorr <- try(chol(template)) # base::chol or Matrix::chol 
          if (inherits(cholcorr,"try-error")) stop(cli::format_error("A correlation matrix is (nearly) singular. Check the correlation model and/or see {.topic [sparse_precision](spaMM::sparse_precision)}.")) 
          preclist[[rd]] <- drop0(chol2inv(cholcorr), tol = .Machine$double.eps)
          ZAlist[[rd]] <- .addrightcols_Z(Z=ZAlist[[rd]], colnames(template), verbose=FALSE)
          corrlist[[rd]] <- .fill_denseL(nc_ori, nc=ncol(ZAlist[[rd]]), LHS_nlev, 
                                         existing=drop0(cholcorr,tol=0)) 
          # The Gryphon example is one where denseness is low for chol2inv(cholcorr) and for drop0(cholcorr)
          # but not for cholcorr itself. 
          # In that Gryphon case the user-level input corrMatrix is a dense matrix with 0's and
          # corr_info$corrMatrices[[rd]] is a dsC matrix already without drop0(.,tol=0)-able elements.
          #
          # To assess ZL density, ## : corrlist[[rd]] will be later set to a dummy fully-dense 
          # toeplitz lower triangle accounting for any .addrightcols_Z()-ed empty cols. 
        } 
      } 
    } else  if (type_rd == c("adjacency")) { # preclist[[rd]] presumably has 0s on diagonal
      nc <- ncol(preclist[[rd]])
      corrlist[[rd]] <- lower.tri(matrix(TRUE,ncol=nc,nrow=nc),diag = TRUE) # template matrix created in faster way than diag()
      rowmax <- max(rowSums(preclist[[rd]]))
      preclist[[rd]] <- preclist[[rd]] + .symDiagonal(n=nc) # ,x=rowmax+1) # make it diagonally dominant
    } else if (type_rd %in% c("AR1", "IMRF") ) { 
      nc <- ncol(ZAlist[[rd]])
      if (LHS_nlev>1L) nc <- nc %/% LHS_nlev 
      corrlist[[rd]] <- lower.tri(matrix(TRUE,ncol=nc,nrow=nc),diag = TRUE) # template matrix created in faster way than diag()
      if (is.null(preclist[[rd]])) preclist[[rd]] <- 
        Matrix::toeplitz(as(c(1, 1/2, rep(0, nc-2L)), "sparseVector"))
    } else if (type_rd == c("corrFamily") ) {
      if ( ! is.null(corrFamlist[[rd]][["f"]])) { # otherwise structure is ignored
        template <- corrFamlist[[rd]]$template
        if (inherits(template,"precision")) {
          preclist[[rd]] <- template$matrix
          corrlist[[rd]] <- .fill_denseL(nc_ori, nc=nc_ori, LHS_nlev, existing=NULL)   ## consistently with s=corrMatrix code, explained above
        } else {
          cholcorr <- try(chol(template)) # base::chol or Matrix::chol 
          if (inherits(cholcorr,"try-error")) stop("The matrix template does not look like a symmetric positive definite matrix. Check the 'tpar' argument..") 
          preclist[[rd]] <- drop0(chol2inv(cholcorr), tol = .Machine$double.eps)
          ZAlist[[rd]] <- .addrightcols_Z(Z=ZAlist[[rd]], colnames(template), verbose=FALSE)
          corrlist[[rd]] <- .fill_denseL(nc_ori, nc=ncol(ZAlist[[rd]]), LHS_nlev, existing=drop0(cholcorr,tol=0)) ## consistently with s=corrMatrix code, explained above
        }
      }
    } else if (type_rd %in% c("Matern","Cauchy")  ) {
      if ( rd %in% which_nested_RHS ) { 
        RHS_nesting_info <- attr(ZAlist[[rd]],"RHS_nesting_info")
        blcks <- lapply(RHS_nesting_info$nested_Zcols, function(v) {
          nc <- length(v) 
          matrix(TRUE,ncol=nc,nrow=nc)
        })
        preclist[[rd]] <- Matrix::bdiag(blcks) # slow in nested-Matern example
        blcks <- lapply(blcks, lower.tri, diag=TRUE)
        corrlist[[rd]] <- Matrix::bdiag(blcks)
        # Add names for for check in .compute_ZAXlist()
        # The precise order of names here may not be essential, but this could be a template for other algorithms where the oreder would matter.
        perm_ZAnames <- RHS_nesting_info$full_LEVELS[.unlist(RHS_nesting_info$nested_Zcols)] # must be (permuted) ZAnames
        rownames(corrlist[[rd]]) <- colnames(corrlist[[rd]]) <- perm_ZAnames
        corrlist[[rd]] <- corrlist[[rd]][RHS_nesting_info$full_LEVELS,RHS_nesting_info$full_LEVELS] # now in ZA order
      } else {
        nc <- ncol(ZAlist[[rd]])
        preclist[[rd]] <- matrix(TRUE,ncol=nc,nrow=nc) # as(allTRUE,"lgCMatrix") #new("lgCMatrix",i=rep(c(0L,seq_len(nc-1L)),nc),p=c(0L,seq_len(nc))*nc,x=rep(TRUE,nc^2),Dim=c(nc,nc)) 
        corrlist[[rd]] <- lower.tri(preclist[[rd]],diag = TRUE) # logi
      } 
    } else  if ( is.null(corrlist[[rd]]) &&  ! is.null(preclist[[rd]])) { 
      warning("Possibly inefficient code in .calc_G_diagnosis_new(). Please contact the maintainer.") # _F I X M E__
      tcrossfac_adj <- mat_sqrt(preclist[[rd]]) # tcrossfac hence t(solve()) is the tcrossfac of the corr mat (<=> Lunique) which is the following backsolve
      if ( attr(tcrossfac_adj,"type")=="cholL_LLt") {
        if (inherits(tcrossfac_adj,"dtCMatrix")) {
          corrlist[[rd]] <- solve(t(tcrossfac_adj)) 
        } else corrlist[[rd]] <- .backsolve(tcrossfac_adj,upper.tri = FALSE, transpose=TRUE)
      } else corrlist[[rd]] <- t(solve(tcrossfac_adj)) # quick patch when another facto used.  
    } 
    if (LHS_nlev>1L ) { 
      if ( ! is.null(corrlist[[rd]])) { # it may already have an LHS_nlev-adjusted size, so:
        if (ncol(ZAlist[[rd]])> nrow(corrlist[[rd]])) {
          rownames_kron_Y <- rownames(corrlist[[rd]])
          corrlist[[rd]] <- kronecker(lower.tri(matrix(1,ncol=LHS_nlev,nrow=LHS_nlev),diag=TRUE), 
                                      corrlist[[rd]]) # make.dimnames does not seem to work with Matrix'es
          if (! is.null(attr(ZAlist[[rd]],"RHS_nesting_info"))) {
            rownames(corrlist[[rd]]) <- rep(rownames_kron_Y,LHS_nlev)
          } # bc dimnames may be needed deep in .compute_ZAL() -> .compute_ZAXlist()?
        }
      }
    } 
    if ( is.null(preclist[[rd]])) { # "(.|.)" 
      ## tnb <- fitme(resp~1+(1|ID), data=lll,family=Tnegbin(2)) is a test case where spprec is clearly slower
      preclist[[rd]] <- .symDiagonal(TRUE,n=ncol(ZAlist[[rd]])/LHS_nlev) # .symDiagonal(TRUE,n=ncol(ZAlist[[rd]]))  ## (null corrlist[[rd]] must mean the same thing) 
    }
    fac <- ncol(ZAlist[[rd]])/ncol(preclist[[rd]])
    if (fac>1L) { # ie if LHS_nlev not yet accounted for precision matrix
      preclist[[rd]] <- kronecker(lower.tri(matrix(1,ncol=fac,nrow=fac),diag=TRUE), 
                                  preclist[[rd]]) 
    }
    preclist[[rd]] <- as(forceSymmetric(preclist[[rd]]),"CsparseMatrix")
    
  }
  # suppressMessages here and below as .provide_G_diagnosis() is not the right context for messages.
  ZL <- suppressMessages( .compute_ZAL(XMatrix=corrlist,ZAlist,as_matrix = FALSE) )  ## > qq s for large ZA 
  dimZL <- dim(ZL)
  
  locQ <- do.call(Matrix::bdiag, list(preclist)) # dsC
  Z_ <- suppressMessages( .compute_ZAL(XMatrix=NULL,ZAlist,as_matrix = FALSE, bind.=TRUE, force_bindable=FALSE) )
  locG <- .crossprod(Z_)+locQ # ideally dsC except if Z_ is really dense
  c_sparse <- 10000/ncol(locG)^3 
  # I cannot really use my earlier 5/3 pow as ref since corss_ZL was used
  ZLdim_fac <- (min(dimZL)^2 *max(dimZL)/(ncol(locG)^3))
  algfacs <- .spaMM.data$options$algfacs
  costs <- c("spprec"= algfacs[["spprec"]]*
               (c_sparse+ .calc_denseness(locG, relative=TRUE)),
             "spcorr"= algfacs[["spcorr"]]*
                          .calc_denseness(ZL, relative=TRUE)*ZLdim_fac^(2/3) ,
             "decorr"=                                       ZLdim_fac^(3/4)
             )
  G_diagnosis <- list(costs=costs, fast=FALSE)
  # print(unlist(G_diagnosis))
  corr_info$G_diagnosis <- G_diagnosis
  
    # corr_info is an environment so no need for return value.
}

.provide_G_diagnosis <- local({
  time_warned <- FALSE
  diagnosis_time <- 0
  function(corr_info, ZAlist, fast=.spaMM.data$options$fast_G_diagnosis) { # default is fast=FALSE...
    # .assign_geoinfo_and_LMatrices_but_ranCoefs() may serve as a template, but we don't want actual matrices except to assess computation costs
    if (is.null(corr_info$G_diagnosis)) {
      if ( ! time_warned) time1 <- Sys.time()
      if (fast) {
        corr_info$G_diagnosis <- .calc_G_diagnosis_fast(corr_info=corr_info, ZAlist=ZAlist)
      } else {
        corr_info$G_diagnosis <- .calc_G_diagnosis_new(corr_info=corr_info, ZAlist=ZAlist)
      }
      if ( ! time_warned) {
        time1 <- .timerraw(time1)
        if (time1>1) diagnosis_time <<- time1
      }
    }
    return(corr_info$G_diagnosis)
  }
})



# Does not set 'algebra' (see .set_augX_methods() instead); 
# sets spprec boolean depending on two optional user inputs, 'algebra' and 'sparse_precision'
.wrap_determine_spprec <- function(control.HLfit, ZAlist, processed, X.pv) {
  # .spaMM.data$options[[".algebra"]] <- # ugly but convenient to check user input in .ad_hoc_dsy_warning()
    algebra <- control.HLfit$algebra
  if (is.null(algebra)) {
    sparse_precision <- control.HLfit$sparse_precision
    if (is.null(sparse_precision)) {
      sparse_precision <- 
        .determine_spprec(ZAlist=ZAlist, processed=processed, # uses $For, $corr_info, $init_HLfit
                          X.pv=X.pv) ## possibly writes $corr_info$G_diagnosis! .../...
    }
  } else sparse_precision <- (algebra=="spprec")
  sparse_precision
}

if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") {
  
  .is01col <- function(ZA) {
    is01col <- attr(ZA,"is_incid")
    if (is.null(is01col)) { # But presumably is_incid has been checked previously to .is01col() beign called?
      warning("(_LOCAL_TESTS_-check only:) 'is_incid' attribute missing.\n Might be correct in some cases but is it here?",
              immediate. = TRUE)
      return(FALSE)
    }
    if (is01col) {
      is01col <- attr(is01col,"is01col") 
      if (is.null(is01col)) {
        warning("(_LOCAL_TESTS_-check only:) Possible inefficiency: 'is01col' attribute missing.", immediate. = TRUE)
        return(FALSE)
      }
    }
    is01col
  }
  
} else {
  
  .is01col <- function(ZA) {
    is01col <- attr(ZA,"is_incid")
    if (is01col) {
      is01col <- attr(is01col,"is01col") 
      if (is.null(is01col)) return(FALSE) # __F I X M E___ See alternative, _LOCAL_TESTS_ definition for devel checks.
    }
    is01col
  }
  
}

# OLD, conceived for use with 'fast' G diagnosis
.old_QR_from_fast_G_diagnosis <- function(corr_info, ZAlist, is_spprec, G_diagnosis, 
                                          sp_alg_thresholds, nrand) {
  
  QR_diagnosis <- .provide_QR_diagnosis(corr_info=corr_info, ZAlist=ZAlist)
  stillpossiblydense <- QR_diagnosis$stillpossiblydense
  actuallysparse <- QR_diagnosis$actuallysparse
  if (( ! is_spprec) && all(stillpossiblydense)) { ## simple subcase of the next case => exclude SPCORR, but not SPPREC
    ## LMatrices are not available, and it may be better to use the density of the correlation matrices anyway:
    ## for maximally sparse Z, ZL's denseness is that of the retained rows of L. This suggests that ZL could be made sparser 
    ## by reordering the levels of the correlation matrix so that the most represented levsl come first in a triangular L factor. 
    ## But this would not affect the denseness of .crossprod(ZW) in .get_absdiagR_blocks(), 
    ## and this leads to use "dense" whenever the correlation matrix is dense.
    ## Gryphon example is useful here, L is sparse but the "dense" QRmethod still appears as good as the "sparse" one.
    # Next test is good for vullioud2018 otherwise decorr is selected
    if (G_diagnosis$rel_denseness_cross_ZL < sp_alg_thresholds[["spcorr"]]) {
      return("sparse")
    } else return("dense") # "that WON'T prevent SPPREC from being possibly selected!" (but cf test '! spprec')
    # But for the smallest networks in the CAR_timings, for which spprec is not used, decorr is thus used.
    # But even if .provide_G_diagnosis() were run, decorr would still be selected. The reason is that the true denseness via ZL is not assessed;
    # instead L is assumed maximally dense. For those smallest networks, this is not optimal bc the actual ZL would be quite sparse.
    # but these simulated adjacency matrices are highly unusual and a bit silly (most of their eigenvalues are zero). sparse spatial correlation matrices should be unusual.
    # So there would be little interest in trying to assess the sparsity of the correl matrix, even if a fast & simple algo were available for that.
  } else if (any(stillpossiblydense) || any(actuallysparse)) {
    if (G_diagnosis$crossZL_is_dsy) { ## sufficient, but loose, condition for using dense
      return("dense")
    } else {
      if (any(actuallysparse)) {
        if (G_diagnosis$rel_denseness_noAR < sp_alg_thresholds[["spcorr"]]) { 
          return("sparse")
        } else {
          return("dense") ## ZAlist actually not so sparse
        }
      } else { # some stillpossiblydense, no actuallysparse  # the difference of decision rule in this case appeared more clearly after tidying the code, but is old
        totdim <- 0L
        for (it in seq_along(ZAlist)) totdim <- totdim + dim(ZAlist[[it]])
        if (totdim[1L]>4L*totdim[2L]) {
          return("sparse")
        } else return("dense")
      }
    }
  } else if (nrand==1L && .is_identity(ZAlist[[1]])) { ## test pertinent slmt pour non-spatial models !
    return("sparse") ## special case for poisson or binomial with saturated ranef
  } else { ## several block effects...
    if ((G_diagnosis$rel_denseness_noAR < sp_alg_thresholds[["spcorr"]])) { 
      return("sparse")
    } else {
      return("dense") ## ZAlist actually not so sparse
    }
  } 
}

.calc_default_QR_method <- 
  function(corr_info, ZAlist, is_spprec,
           G_diagnosis=.provide_G_diagnosis(corr_info=corr_info, ZAlist=ZAlist),
           sp_alg_thresholds=.spaMM.data$options$sp_alg_thresholds
  ) {
  nrand <- length(ZAlist)
  if (nrand>0L) {
    # adjacency speed to be tested on 2nd example from test-spaMM.R
    #
    if ( ! .spaMM.data$options$fast_G_diagnosis) {
      if (which.min(G_diagnosis$costs[c("spcorr","decorr")])==1L) {
        return("sparse")
      } else return("dense")
    } else {
      .old_QR_from_fast_G_diagnosis(
        corr_info=corr_info, ZAlist= ZAlist, is_spprec=is_spprec, 
        G_diagnosis=G_diagnosis, sp_alg_thresholds=sp_alg_thresholds, nrand=nrand)
    }
  } else return("dense") 
}


# Does *not* set 'algebra' (see .set_augX_methods() instead): 
# sets dense/sparse for the case where "corr"elation method will be chosen.
# Return value when is_spprec is TRUE (as previously set by .wrap_determine_spprec()) may be ignored.
# When is_spprec is FALSE the function addresses the issue that 
# even when the Z's are sparse postmultplication by LMatrix leads some of the ZAL's to dgeMatrix (dense).
.choose_QRmethod <- function(ZAlist, corr_info, is_spprec, processed, control.HLfit) {
  if (is_spprec) return("sparse") # 08/2021: currently QRmethod operates only though .eval_as_mat_arg() 
  # which ignores QRmethod when is_spprec is TRUE (so returning NULL or NaN should have the same effect).
  if ( ! is.null(algebra <- control.HLfit$algebra)) {
    if (algebra=="decorr") return("dense")
    if (algebra=="spcorr") return("sparse")
  }
  # ELSE
  if ( is.null(QRmethod <- .spaMM.data$options$QRmethod) ) {
    QRmethod <- .calc_default_QR_method(corr_info=corr_info, ZAlist=ZAlist, is_spprec=is_spprec)
  }
  QRmethod
}

.set_augX_methods <- function(processed) {
  # sets algebra and possibly *two* methods...
  if (.eval_as_mat_arg(processed)) {
    algebra <- "decorr"
    processed$sXaug_method <- .spaMM.data$options$matrix_method
  } else {
    algebra <- "spcorr"
    if (processed$how$obsInfo) {
      processed$sXaug_method <- .spaMM.data$options$Hobs_Matrix_method
    } else {
      processed$sXaug_method <- .spaMM.data$options$Matrix_method
    }
  } # NOT alternative to first block (in spprec we still need  processed$sXaug_method for .MakecovEst1)
  if (processed$is_spprec) {
    algebra <- "spprec" 
    processed$spprec_method <- .spaMM.data$options$spprec_method 
  }
  processed$how$algebra <- algebra
  algebra
}

.check_time_G_diagnosis <- function(.provide_G_diagnosis, processed, algebra) {
  if (( ! environment(.provide_G_diagnosis)$time_warned) && 
      (time1 <- environment(.provide_G_diagnosis)$diagnosis_time)) {
    message(cli::format_message(paste0("(One-time message:) Choosing matrix methods took ",time1," s.\n",
                   "  If you perform many similarly costly fits, setting the method\n",
                   '  by control.HLfit=list(algebra=<"spprec"|"spcorr"|"decorr">) may be useful,\n',
                   '  see {.topic [algebra](spaMM::algebra)}. "',algebra,'" has been selected here.'))) 
    environment(.provide_G_diagnosis)$time_warned <- TRUE
  }
}

.any_gaussian_inverse <- function(processed) {
  if (is.null(family <- processed$family)) {
    families <- processed$families
    for (mv_it in seq_along(families)) {
      fam <- families[[mv_it]]
      if (fam$link=="inverse" && fam$family=="gaussian") return(TRUE) 
    }
    return(FALSE)
  } else return(family$link=="inverse" && family$family=="gaussian")
}

.preprocess_LevM <- function(user_LM, processed, nrand) {
  if (attr(processed[["models"]],"LMMbool")) user_LM <- FALSE  # (_F I X M E_) removing this and forcing LevM creates an error in the tests, meaning that LevM does not handle LMMs (OK)
  if (is.null(user_LM)) user_LM <- .spaMM.data$options$LevenbergM
  if (is.list(user_LM)) stop("is.list(LevenbergM)")
  # we may want: 
  # no Levenberg: user's LevenbergM=FALSE
  # Levenberg from start: user's LevenbergM=TRUE
  # optional LevenbergM, with start as decided by following code: user's LevenbergM=NULL
  # full control overriding code below: e.g. user's LevenbergM=c(user_LM=TRUE, LM_start=TRUE) (not API)
  if (length(setdiff(c("user_LM","LM_start"), names(user_LM)))) {
    if ( ! is.logical(user_LM)) user_LM <- NA # handles default case where user_LM is NULL 
    if (is.na(user_LM)) { 
      if (processed$bin_all_or_none ) { 
        if (nrand>0 && processed$HL[1L]==0L) {
          ## PQL/L + LevenbergM combine safely and relatively fast.... for small data
          # bigranefs -> PQL/L+LevM much smaller than ML!
          ## adjlg has 1000 levels and is faster without LevM 
          LM_start <- (tail(processed$cum_n_u_h,n=1)<500L)[[1]] # drop the automatic name from cum_n_u_h...
        } else LM_start <- FALSE  ## BINARYboot test to assess effect on timings
      } else if (.any_gaussian_inverse(processed)) {
        LM_start <- TRUE # maximum safety in numerically difficult case. test-devel-LLM example shows this is necessary.
      } else LM_start <- FALSE
    } else LM_start <- user_LM[[1L]] # important to drop name else the names of the vector are wrong
    return(c(user_LM=user_LM[[1L]], LM_start=LM_start) ) 
  } else return(user_LM) # allows a full vector to be user-provided, for full control
}

.subset_corrMatrix <- function(corrMatrix, ZAnames, corrnames=.get_rownames_corrMatrix(corrMatrix)) {
  # if (is.null(ZAnames)) { ## set by .calc_Zlist() or .calc_ZAlist(), with two cases for corrMatrix # with repeated names for ranCoefs sensu lato. 
  #   stop("NULL colnames in (a block of) the design matrix for random effects. Some mishandling of 'AMatrices'?")
  # }
  # 
  # if (is.null(colnames(corrMatrix))) { # copy row names to col names
  #   if (inherits(corrMatrix, c("matrix", "Matrix"))) {
  #     colnames(corrMatrix) <- corrnames 
  #   }
  #   else if (inherits(corrMatrix, "precision")) {
  #     colnames(corrMatrix[["matrix"]]) <- corrnames
  #   }
  # }
  # 
  # In the nested RHS case the ZAnames have additional :1, :2 in comparison to the corrnames. So the next test is FALSE and no subsetting/reordering is performed.
  if ( ! length(extraZAnames <- setdiff(ZAnames,corrnames))) { ## i.e. if all ZAnames in corrnames
    if ( inherits(corrMatrix,"precision")) { 
      # do not subset a precmat => nothing here but the corresponding Z matrix may be modified (.addrightcols) by .init_assign_geoinfo() 
      # and then the precision matrix may be reordered according to the cols of this possibly augmented Z. 
    } else {
      uZAnames <- unique(ZAnames)
      if ( length(setdiff(corrnames,uZAnames)) || any(corrnames!=uZAnames) ) { # reordering and subsetting. corr_info$corrMatrices still contains the full matrix
        perm <- pmatch(uZAnames, corrnames)
        if (inherits(corrMatrix,"dist")) {
          corrMatrix <- (proxy::as.matrix(corrMatrix,diag=1)[perm,perm]) ## IF diag missing in input corrMatrix THEN assume a correlation matrix
          ## it's not useful to convert back to dist (either uglily by as.dist(); or package 'seriation' has (permute.dist-> C code)
        } else corrMatrix <- corrMatrix[perm,perm]  
      } ## else orders already match
    }
  } 
  return(corrMatrix)
}

.addrightcols_Z <- function(Z, precnames, verbose=TRUE) {
  # We have tested in .init_assign_geoinfo() -> .check_rownames_corrMatrix() whether all ZAnames were in precnames 
  # so (1) There is no need to check names again before calling this fn in .provide_G_diagnosis()
  #    (2) the only possible difference between sets of names is additional names in precnames
  nestedZAnames <- .get_nestednames_ZA(Z)
  supplevels <- setdiff(precnames,nestedZAnames)
  if (length(supplevels)) {
    next_nestednames <- c(nestedZAnames, supplevels)
    ncol_prec <-  length(precnames)
    ncol_Z <- ncol(Z)
    LHS_nlev <- length(attr(Z,"namesTerm"))
    nlev_RHS <- ncol_Z/LHS_nlev
    if (suppcols <- ncol_prec-nlev_RHS) {
      if (LHS_nlev>1L ) { # suppcols are per block; possibly poorly tested code
        if (verbose) message(paste("Note: Precision matrix has", suppcols, 
                                   "more levels than there are in the data.")) # and <0 values are a bug...
        # add cols of zeros one the right; cols to be reordered next.
        Zp <- Z@p
        colindices <- matrix(seq(ncol_Z),ncol=LHS_nlev)
        colindices <- rbind(colindices,
                            matrix(rep(colindices[nlev_RHS,],suppcols),ncol=LHS_nlev,byrow = TRUE))
        dim(colindices) <- NULL
        Z@p <-  c(0L,Zp[colindices+1L])
        Z@Dim[2L] <- next_ncol <- ncol_prec*LHS_nlev
        nextcolnames <- matrix(colnames(Z),ncol=LHS_nlev)
        nextcolnames <- rbind(nextcolnames,
                              matrix( # paste0(supplevels,":dummy"), # ":dummy" for easy tracking in case of problem .. # NO, bc:
                                supplevels, # set of names must be the same as precnames to allow permutation of the cov_info_mat
                                nrow=suppcols, ncol=LHS_nlev))
        dim(nextcolnames) <- NULL 
        Z@Dimnames[[2L]] <- nextcolnames
        suppcols <- nlev_RHS+seq(suppcols)
      } else { # (names=levels) in precmat but not in the data
        if (verbose) message(paste("Note: Precision matrix has", suppcols, 
                                   "more levels than there are in the data.")) # and <0 values are a bug...
        # add cols of zeros on the right; cols to be reoredered next.
        if ( inherits(Z,"ddiMatrix")) Z <- .as_ddi_dgC(Z)
        Zp <- Z@p
        Z@p <- c(Zp, Zp[length(Zp)] + rep(0L,suppcols))
        Z@Dim[2L] <- ncol_prec
        Z@Dimnames[[2L]] <- next_nestednames
        suppcols <- ncol_Z+seq(suppcols)
      }
      if ( ! is.null(RHS_nesting_info <- attr(Z,"RHS_nesting_info"))) {
        RHS_nesting_info$full_LEVELS <- nextcolnames
        RHS_nesting_info$nested_LEVELS <- next_nestednames
        RHS_nesting_info$nested_Zcols[["dummy"]] <- nlev_RHS+suppcols
        attr(Z,"RHS_nesting_info") <- RHS_nesting_info
      }
    }
  }
  Z
}

.preprocess_pw <- function(subs_p_weights, nobs, model_frame) {
  if (is.null(subs_p_weights)) {
    prior.weights <- structure(rep(1L,nobs),unique=TRUE, is_unit=TRUE) ## <- 1L prevented by glm -> model.frame(... prior.weights)
  } else if (inherits(subs_p_weights,"call") && subs_p_weights[[1L]] == "quote")  {## 'prior.weights' is a quoted expression
    prior.weights <- structure(subs_p_weights, unique=FALSE, is_unit=FALSE)
  } else {
    prior.weights <- as.vector(stats::model.weights(model_frame)) ## as.vector as in say lm() protects against array1d
    if ( ! is.numeric(prior.weights)) 
      stop("'weights' must be a numeric vector")
    if (any(prior.weights < 0)) 
      stop("negative weights not allowed")
    is_unique <- length(unique(prior.weights))==1L
    prior.weights <- structure(prior.weights, unique=is_unique, is_unit=(is_unique && prior.weights[1]==1))
    #attr(prior.weights,"only1") <- all(upw==1L)
  }   
  return(prior.weights)
}

################################################################################
# for a random slope term, ie v= v_1+x v_2 , the x went into the general ZAL matrix 
# (construction of Zlist by .calc_Zlist(), and
# we are still estimating the lambda's using a X_lamres with 0/1's only
# unless there is a non-trivial model for the lambdas
################################################################################
.calc_X_lamres <- function(processed, models=processed$models, ZAlist=processed$ZAlist, nrand=length(ZAlist)) {
  if (all(models[["lambda"]]=="lamScal")) { ## all mixed models handled in 06/2015 (random slope, adjacency...) hence currently always TRUE
    Xi_cols <- attr(ZAlist,"Xi_cols")
    if (any(Xi_cols>1 & ! processed$lcrandfamfam=="gaussian")) {
      stop("(!) random slope models implying correlated non-gaussian random effects are not fitted.")
    }
    cum_Xi_cols <- cumsum(c(0, Xi_cols)) ## if two ranef,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
    n_u_h <- rep(0, sum(Xi_cols))
    #for (i in 1:nrand) n_u_h[(cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]] <- ncol(ZAlist[[i]]) ##  nlevels(Subject[[i]])
    # if 18 responses in a random slope model ncol(ZAlist[[i]]) is 36 while nlevels(Subject[[i]]) was 18
    for (i in seq_len(nrand)) n_u_h[cum_Xi_cols[i]+seq(Xi_cols[i])] <- ncol(ZAlist[[i]])/Xi_cols[i]
    # h_u_h not n_u_h ...
    cum_h_u_h <- cumsum(c(0, n_u_h)) ## if two "Intercept" ranefs,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
    ## if (1+X|...) +(1|...),  with n_u_h=(3,4), this is 0,3,6,10. cum_h_u_h[sum(Xi_cols)+1] is then 10, the total # of realizations
    X_lamres <- matrix(0,cum_h_u_h[sum(Xi_cols)+1L],sum(Xi_cols))
    colnames(X_lamres) <- unlist(attr(ZAlist,"namesTerms"))
    for (i in seq_len(nrand)) {
      for (j in (cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]) {
        X_lamres[(cum_h_u_h[j]+1L):cum_h_u_h[j+1L],j] <- 1L ## this maps the deviance residuals to the lambda's to be estimated from them. None of the random-slope columns is a constant full column because each dev res is used for estimating only one lambda. Nevertheless, the first col will be called "(Intercept)", and this makes a valid output.
      }
    }
    return(X_lamres)
  } else {  ## linear predictor for variance of random effects (lambda) (lamGLM or lamHGLM) 
    if (any(models[["lambda"]]=="lamHGLM")) { ##need distinct calls... to fit each lambda model  
      if (length(formulaLambda)==2) formulaLambda <- as.formula(paste('"lambda"',paste(formulaLambda,collapse=" ")))
      if (!is.null(.parseBars(formulaLambda))) {  ## lamHGLM
        models[["lambda"]] <- list("lamHGLM")
      } else models[["lambda"]] <- list("lamGLM")  
      colnames(X_lambda) <- colnames(fr_lambda$X) ## but code not effective, fr_lambda not computed
    } else { ## can use a single design matrix for all random effects, which is convenient.
      stop("LIKELY missing code to handle linear predictor for lambda.")
      # la suite c'est dexu residus de code a assembler: il faut une liste de terms_info du type
      fr_lambda <- .get_terms_info(formula=formulaLambda,data="data", famfam="famfam") ## but the "data" should probably be distinct data here, with nrow=number of reals of ranefs 
      # (pobablement calculee en amont pour determiner lamScal aussi...) ensuite extraire les design matrices
      #X_lamres ? Xi_cols ?
    }
  } 
}

.preprocess_resid <- function(preprocess_arglist) {
  residProcessed <- do.call(.preprocess,preprocess_arglist) ## cf verbose explicitly set to NULL 
  # we add ".phi" to attr(residProcessed$predictor - for summary() only ? But then same operation on version with hyper-ranefs
  fullform <-  .preprocess_formula(as.formula(paste(".phi",.DEPARSE(residProcessed$predictor))))
  mostattributes(fullform) <- attributes(residProcessed$predictor)
  if ( ! is.null(hy_form <- attr(fullform,"hyper_info")$formula)) attr(fullform,"hyper_info")$formula <- as.formula(paste(".phi",.DEPARSE(hy_form)))
  residProcessed$predictor <- fullform
  return(residProcessed)
}

.check_identifiability_LMM <- function(processed, nobs) {
  ## identifiability checks cf modular.R -> checkNlevels() in lmer:
  vec_n_u_h <- diff(processed$cum_n_u_h)
  if (any(vec_n_u_h<2L) && is.null(processed$phi.Fix)) {
    problems <- which(vec_n_u_h<2L) 
    for (rd in problems) {
      mess <- paste0("Only ",vec_n_u_h[rd]," level for random effect ",
                     attr(processed$ZAlist,"exp_ranef_strings")[rd],
                     ";\n   its variance may be confounded with phi.")
      warning(mess, immediate.=TRUE)
    }
  }
  if (any(vec_n_u_h==nobs) && processed$models[["phi"]] %in% c("phiScal","phiGLM")) { 
    if (attr(processed$residModel$formula,"has_intercept")) { ## there is an intercept in the resid.model formula
      # ideally for random-coefficients models we should compare the design columns... 
      ## FR->FR cf isNested check as in https://github.com/lme4/lme4/blob/master/R/utilities.R, 
      problems <- which(vec_n_u_h==nobs) 
      for (rd in problems) {
        term_ranef <- attr(processed$ZAlist,"exp_ranef_strings")[rd]
        if (substr(term_ranef, 1, 1)=="(" ## excludes spatial (and more generally 'keyword') ranefs 
            && ! is.numeric(processed$lambda.Fix[rd])
        ) {
          mess <- paste0("Number of levels = number of observations for random effect ", term_ranef,
                         ";\n   this model cannot be fitted unless phi is fixed, or the variance",
                         "\n   of this effect is fixed, or a non-trivial correlation matrix is given.") 
          stop(mess)
        }          
      }
    }
  }
}

.check_phi_Fix <- function(phi.Fix, family) {
  if ( ! (constr_fit <- ! is.null(phi.Fix))) {
    if (constr_fam <- ! family$family %in% c("gaussian","Gamma","tweedie")) {
      phi.Fix <- 1 
    } # else if (var(y)==0) phi.Fix <- .spaMM.data$options$min_disp
  } else if (any(phi.Fix==0)) stop("phi cannot be fixed to 0.")
  if ( ! is.null(phi.Fix)) phi.Fix <- structure(phi.Fix,
                                                constr_fit= constr_fit,
                                                constr_phi= (constr_fit || constr_fam) ) 
  # so that the fitobject$phi gets these attributes.
  return(phi.Fix)
}

.get_residProcessed <- function(processed, control.HLfit, resid.model, HLmethod, 
                                resid.formula, data, control.glm, print_phiHGLM_info) {
  preprocess_arglist <- list(control.HLfit=control.HLfit, ## constrained
                             ranFix=resid.model$fixed, 
                             HLmethod=HLmethod, ## constrained
                             predictor=resid.formula, ## obvious
                             resid.model=resid.model$resid.model, # potentially allows nested resid.model's... 
                             REMLformula=NULL, # constrained
                             data=data, # obvious (?) 
                             family=resid.model$family, # obvious, gaussian
                             BinomialDen=NULL, # obviously no binomial response
                             rand.families=resid.model$rand.family, # (NULL not handled by preprocess); 
                             #   outer preprocess calls *receive* a default value from formals(HLfit)
                             etaFix=resid.model$etaFix, ## not constrained, but should rather use 'resid.model$fixed'
                             prior.weights=NULL, ## currently defined  dynamically using lev_phi...
                             control.glm=control.glm, ## constrained
                             verbose=c(print_phiHGLM_info=print_phiHGLM_info), ## TRACE would be overriden by the final do_TRACE call of the parent .preprocess()
                             For="fitme", ## constrained: preprocess must allow spatial and non-spatial models
                             init.HLfit=as.list(resid.model$init.HLfit), ## converts NULL to list() as exp'd by .preprocess()
                             CONTROL=list(ppc_reactvt_warn=TRUE, dyndyn=FALSE)
                             )
  ## preprocess formal arguments that were ignored up to v.2.4.30 14/05/2018:
  other_preprocess_args <- setdiff(names(formals(.preprocess)),names(preprocess_arglist))
  preprocess_arglist[other_preprocess_args] <- resid.model[other_preprocess_args]
  .preprocess_resid(preprocess_arglist)
}

.preprocess_phi_model <- function(processed, models, resid.model, control.HLfit, HLmethod, data, 
                                  control.glm, family) {
  residFrames <- NULL
  resid.model$fixed <- .preprocess_fixed(resid.model$fixed)
  resid.formula <- resid.model$formula
  is_mixed <- ! is.null(.parseBars(resid.formula)) # will be identified as "phiHGLM"
  has_etaFix <- ! is.null(resid.model$etaFix) # NOT to be confused with variable etaFix used for outer optim.
  mainfamfam <- family$family
  if ( is.null(processed$phi.Fix) ) {
    if (is_mixed && is.null(resid.model$rand.family)) resid.model$rand.family <- gaussian() # avoids rand.families being NULL in .preprocess_resid() -> .preprocess()
    
    if (is_mixed || # phiHGLM
        has_etaFix # resid (phiGLM) with etaFix: we use the residProcessed interface as a quick implementation of this case
                   # Such a fixed etaFix can, and should better, be fitted by an offset (alternative below).
        ) { 
      processed$residProcessed <- .get_residProcessed(processed=processed,
        control.HLfit=control.HLfit, resid.model=resid.model, HLmethod=HLmethod, 
        resid.formula=resid.formula, data=data, control.glm=control.glm, 
        print_phiHGLM_info=processed$verbose["phifit"] && is_mixed)
      # Message is a nuisance for bootstraps:
      # if (processed$verbose["phifit"] && # should be FALSE in bootstraps bc the message is a nuisance in that context
      #     identical(names(resid.model$fixed$phi),"default")) message("'phi' of residual dispersion model set to 1 by default")
    } 
    
    if (is_mixed) {
      models[["phi"]] <- "phiHGLM" 
      p_phi <- NA
    } else {
      residFrames <- .get_terms_info(formula=resid.formula, data=data, famfam=resid.model$family$family)
      attr(resid.formula,"off") <- off <- model.offset(residFrames$mf) ## only for summary.HLfit() (and below)
      attr(resid.formula,"has_intercept") <- (attr(residFrames$fixef_off_terms,"intercept")!=0L) ## check_identifiability_LMM looks for it
      resid.model$formula <- resid.formula  ## put it back after attributes have been added (no equivalent if phiHGLM)
      if (has_etaFix) {
        models[["phi"]] <- "phiGLM" 
        p_phi <- .old_NCOL(processed$residProcessed$AUGI0_ZX$X.pv) # without the etaFix ones; info for .more_init_optim() to eval allPhiScalorFix which must be non-NA
      } else { # NOT mixed-effect model NOR etaFix. Still allows a fixed offset that should give result equivalent to a fixed etaFix (but with different X)
        ## if formula= ~1 and data is an environment, there is no info about nobs, => fr_disp$X has zero rows, which is a problem later 
        X <- residFrames$X
        p_phi <- .old_NCOL(X)
        namesX_disp <- colnames(X)
        if (p_phi==1L && namesX_disp[1]=="(Intercept)"
            && is.null(off) ## added 06/2016 (bc phiScal does not handle offset in a phi formula) 
        ) {
          models[["phi"]] <- "phiScal"
        } else if (p_phi==0L) { # resid.formula has only an offset term.
          # set phi.Fix so that it is used by fitting functions, instead of running phiGLM code : 
          #   leverages, dev.res , .calc_dispGammaGLM() -> model.frame(), model.matrix()... to find that there is nothing to fit!
          processed$phi.Fix <- resid.model$family$linkinv(model.offset(residFrames$mf)) # fitting fns see this as phi.Fix
          models[["phi"]] <- "phiGLM" # meaningful: see how new offset values are predicted in .calcResidVar()
        } else { 
          models[["phi"]] <- "phiGLM"
          disp_env <- family$resid.model # envir distinct from the processed$residModel *list*
          disp_env$resid.formula <- resid.formula 
          if ( ! is.null(off)) disp_env$off <- off
          disp_env$colnames_X <- namesX_disp
          disp_env$X <- X
        }
      } 
    }
  } else { # Typically when phi.Fix was not NULL. In particular for rdisPars it is presumably 1. models[["phi"]] remains ""
    if ( # ( ! mainfamfam %in% c("gaussian","Gamma")) && 
         .DEPARSE(resid.formula) != "~1") {
      if ( ! mainfamfam %in% c("beta_resp", "betabin", "negbin1","negbin2", "gaussian", "Gamma", "tweedie")) 
        warning(paste0("resid.model may be ignored in ",mainfamfam,"-response models"), immediate. = TRUE)
      resid.formula <- resid.model$formula
      if ( ! is.null(.parseBars(resid.formula))) stop("Random effects are not allowed in model for family parameter.")
      # NOT mixed-effect model NOR etaFix. Still allows a fixed offset that should give result equivalent to a fixed etaFix
      residFrames <- .get_terms_info(formula=resid.formula, data=data, famfam="")
      off <- model.offset(residFrames$mf)
      # attr(resid.formula,"has_intercept") <- (attr(residFrames$fixef_off_terms,"intercept")!=0L) ## for identifiability checks
      ## if formula= ~1 and data is an environment, there is no info about nobs, => fr_disp$X has zero rows, which is a problem later
      X <- residFrames$X
      p_phi <- .old_NCOL(X)
      if (p_phi==1L && colnames(X)[1]=="(Intercept)"
          && is.null(off) ## added 06/2016 (bc phiScal does not handle offset in a phi formula) 
      ) {
        models[["rdispar"]] <- "rdiScal"
      } else {
        disp_env <- family$resid.model # virgin envir distinct from the processed$residModel list
        disp_env$resid.formula <- resid.formula
        if ( ! is.null(off)) disp_env$off <- off
        if (p_phi==0L) { # resid.formula has only an offset term.
          # processed$phi.Fix <- resid.model$family$linkinv(model.offset(residFrames$mf)) # fitting fns see this as phi.Fix
          models[["rdispar"]] <- "rdiOff" # meaningful: see how new offset values are predicted in .calcResidVar()
          # resid.model <- list(formula=resid.formula, off=off)  ## put it back after attributes have been added (no equivalent if phiHGLM has been detected?)
        } else { 
          models[["rdispar"]] <- "rdiForm"
          disp_env$colnames_X <- colnames(X)
          disp_env$X <- X
          # disp_env$init <- resid.model$init
          # disp_env$etaFix <- resid.model$etaFix # finally NOT used by .init_rdisPars().
        }
      }
    } else p_phi <- 0L
  }
  processed$p_fixef_phi <- p_phi # no X_disp is saved in processed
  processed$residModel <- resid.model 
  models
}

.warn_augZXy_scaling_once <- local({
  mess_scaling <- FALSE
  function() {
    if (! mess_scaling) {
      message("'y-augmented' algorithm: in TRACE displays, variable lambda values are shown relative to phi values.")
      mess_scaling <<- TRUE
    }
  } 
})

# This fn first defines the 'tracer' and 'exit' functions, then call trace() using them.
.do_TRACE <-   function(processed) { ## no need for an 'unTRACE' call at the end of each fit since the next fit will cope with everything.
  ## trouble when called from another package while not attached (bboptim example)
  # THe syntax spaMM::HLfit_body, where=spaMM::fitme does not stop() in that case, 
  #     but HLfit_body is not effectively traced when using spaMM directly attached (standard library(spaMM))
  # The syntax spaMM::HLfit_body without where=also does not trace when using spaMM directly attached
  
  level <- processed$verbose["TRACE"]
  if ("package:spaMM" %in% search()) {
    if ( level ) { # may be 0.5...
      if (level >= 1L ) {
        # 'tracer' function:
        if (processed$augZXy_cond) { # .HLfit_body_augZXy() does not have an etaFix argument: See comments in its source file.
          .warn_augZXy_scaling_once()
          tracing_op <- quote(try(.TRACE_fn(fixed, processed=processed))) # the closure of the traced function must have a 'fixed' variable
        } else tracing_op <- quote(try(.TRACE_fn(fixed, etaFix, processed))) # the closure of the traced function must have both 'fixed' and 'etaFix' variable
        if (is.null(processed$REMLformula)) { ## default REML case
          if (processed$HL[[1L]]=='SEM')  {
            objLik <- "logLapp"
          } else objLik <- "p_bv"
        } else {
          if (attr(processed$REMLformula,"isML")) {  
            objLik <- "p_v"
          } else { ## if nontrivial REML formula was used...
            objLik <- "p_bv"
          }
        }
        # 'exit' function
        exit_op <- substitute({
          aphl <- res$APHLs[[objLik]]
          if (is.null(aphl)) {
            print("(objective not found)",quote=FALSE)
          } else print(paste0(objLik,"= ",.prettify_num(aphl,nsmall=4)),quote=FALSE)
        }, list(objLik=objLik))
      } else { ## e.g. level=0.5 : will print only the "progress bars"  
        tracing_op <- quote({})
        exit_op <- quote({})
      }
      
      suppressMessages(trace(processed$HLfit_body_fn, where=asNamespace("spaMM"), print=FALSE, 
                             tracer=tracing_op, # shows the parameters
                             exit=exit_op)) # shows the objective fn
      
      if ( ! is.null(processed$augZXy_body_fn)) {
        suppressMessages(trace(processed$augZXy_body_fn, where=asNamespace("spaMM"), print=FALSE, 
                               tracer=tracing_op, # shows the parameters
                               exit=exit_op)) # shows the objective fn
      }
      for (method_st in c("matrix_method","Matrix_method","spprec_method","Hobs_Matrix_method")) {
        fn <- paste("get_from_MME",strsplit(spaMM.getOption(method_st),"def_")[[1L]][2],sep=".") 
        if (level<4L) {
          suppressMessages(untrace(fn, where=asNamespace("spaMM")))
        } else if (level==4L) {
          suppressMessages(trace(fn,print=FALSE,tracer=quote(cat(which))))
        } else { # level>=5
          suppressMessages(trace(fn,print=TRUE, exit=quote({cat(which); str(resu)})))
        }
      }
    } else { # TRACE=0
      .silent_M_E(untrace(processed$HLfit_body_fn, where=asNamespace("spaMM")))   
      if ( ! is.null(processed$augZXy_body_fn)) {
        .silent_M_E(untrace(processed$augZXy_body_fn, where=asNamespace("spaMM")))   
      }
      for (method_st in c("matrix_method","Matrix_method","spprec_method","Hobs_Matrix_method")) {
        fn <- paste("get_from_MME",strsplit(spaMM.getOption(method_st),"def_")[[1L]][2],sep=".") 
        suppressMessages(untrace(fn, where=asNamespace("spaMM")))
      } 
    }
  } else if (level) {warning("The 'spaMM' package must be *attached* for verbose(TRACE=...) tracing to fully operate",
                             immediate.=TRUE)}
} 

.preprocess_init.HLfit <- function(init.HLfit, corr_info) {
  if ( ! is.null(rho <- init.HLfit$rho)) {
    init.HLfit$corrPars <- list()
    if (length(adj_rd <- which(corr_info$corr_types=="adjacency"))) {
      init.HLfit$corrPars[[as.character(adj_rd)]][["rho"]] <- rho
      init.HLfit$rho <- NULL
    } else stop("Invalid ambiguous 'init.HLfit' argument: single 'rho' but not single adjacency random-effect term.")
  }
  init.HLfit
}


.calc_Binomial_Den <- function(Y, famfam, nobs) {
  if (famfam %in% c("binomial","betabin") && NCOL(Y)>1) {
    BinomialDen <- rowSums(Y)
    if (any(chk <- BinomialDen == 0L)) {
      # if (For_fitmv) {
        return(chk)
      # } else stop("please remove missing data (i.e. for which binomial sample size is 0).")
    }
    ## It's not really possible to remove data at this stage as this may not match the dimension of the distance matrices
    ## moreover one cannot simply remove rows of a matrix "root"...
  } else {
    BinomialDen <- rep(1,nobs)
  }
  BinomialDen
}

.check_y <- function(family, y, BinomialDen) {
  famfam <- family$family
  if (famfam %in% c("binomial","betabin")) {
    if (length(y)==1L || (var(y)==0 && var(BinomialDen)==0) ) { warning("var(response) = 0, which may cause errors.") }  
    bin_all_or_none <- all(pmin(y,BinomialDen-y)==0L)
  } else { 
    bin_all_or_none <- FALSE
    if ( ! is.null(y)) { ## y may be NULL in evaluation of residProcessed
      if ( length(y)==1L || var(y)==0) { # (~1, family=poisson, data=<single response> ) can be fitted
        if (famfam %in% c("gaussian", "Gamma","tweedie")) 
          warning("var(response) = 0, which may cause errors.", immediate. = TRUE) 
      } else if (var(y)<1e-3 && famfam=="gaussian") {
        warning("The variance of the response is low, which may lead to imperfect estimation of variance parameters.\n Perhaps rescale the response?")
      }
    }  
  } # (e1071::svm should fail when var response=0)
  bin_all_or_none
}

.calc_iter_mean_dispFix <- function(control.HLfit, family, y) {
  iter_mean_dispFix <- control.HLfit$iter.mean.dispFix ## documented
  if (is.null(iter_mean_dispFix)) iter_mean_dispFix <- control.HLfit$max.iter.mean ## public control
  if (is.null(iter_mean_dispFix)) {
    if (family$family=="Gamma" && family$link=="log") {
      if (is.null(y)) { 
        iter_mean_dispFix <- NaN # will be replaced when residProcessed$y is known.
      } else iter_mean_dispFix <- max(200L,ceiling(100* mean(abs(log2(y)))))
    } else iter_mean_dispFix <- 200L 
  }
  iter_mean_dispFix  
}

.preprocess_spaMM_tol <- function(bin_all_or_none, control.HLfit) {
  spaMM_tol <- spaMM.getOption("spaMM_tol") 
  if ( ! is.list(spaMM_tol)) stop("spaMM_tol must be a list")
  if (bin_all_or_none) {
    spaMM_tol$rescue_thr <- .spaMM.data$options$LevM_HL11_method$rescue_thr_AoN
  } else spaMM_tol$rescue_thr <- .spaMM.data$options$LevM_HL11_method$rescue_thr_null
  if (spaMM_tol$rescue_thr["rescue"]) {
    spaMM_tol$v_pot_tol <- spaMM_tol$v_pot_tol_rescue
  } else spaMM_tol$v_pot_tol <- spaMM_tol$v_pot_tol_noresc
  # then user's explicit conv.threshold controls
  if ( ! is.null(conv.threshold <- control.HLfit$conv.threshold)) spaMM_tol[["Xtol_rel"]] <- conv.threshold
  # then user's explicit spaMM_tol controls
  if ( ! is.null(user_spaMM_tol <- control.HLfit$spaMM_tol)) spaMM_tol[names(user_spaMM_tol)] <- user_spaMM_tol
  spaMM_tol$fpot_tol <- .spaMM.data$options$fpot_tol # kept separately in the options for easier access
  spaMM_tol$fpot_cond <- spaMM_tol$fpot_tol>0
  spaMM_tol
}

.calc_iter_mean_dispVar <- function(control.HLfit, family, y) {
  iter_mean_dispVar <- control.HLfit$iter.mean.dispVar ## documented
  if (is.null(iter_mean_dispVar)) iter_mean_dispVar <- control.HLfit$max.iter.mean ## public control 
  if (is.null(iter_mean_dispVar)) {
    if (family$family=="Gamma" && family$link=="log") {
      if (is.null(y)) { 
        iter_mean_dispVar <- NaN
      } else iter_mean_dispVar <- max(50L,ceiling(100* mean(abs(log2(y)))))
    } else iter_mean_dispVar <- 50L 
  } 
  iter_mean_dispVar  
}

.preprocess_HL_REMLformula <- function(HLmethod, processed, BinomialDen, nobs, control.HLfit, y, REMLformula) {
  HL <- eval(str2lang(paste0("c",substr(HLmethod,3,100)))) ## extracts the (...) part into a vector
  if (length(HL)==2) HL <- c(HL,1)
  processed$HL <- HL ## ! this may be modified locally !!
  if (HL[1L]=="SEM") processed$SEMargs <- .preprocess_SEMargs(BinomialDen, nobs, control.HLfit, y)
  if (HL[1L]==0L) {processed$p_v_obj <-"hlik"} else processed$p_v_obj <-"p_v" ## objective for beta(_v) estim only: != outer obj 
  if (substr(HLmethod,0,2)=="ML") { # && HL[1]!="SEM") { ## FR->FR c'est bizarre d'exclure le SEM là... p_bv est il vraiment utilisé ?
    if ( ! is.null(REMLformula)) { # input REMLformula should be NULL ***BUT*** processed$REMLformula is not
      stop(paste0("Confusing combination of arguments: 'HLmethod=ML(...)' with non-null 'REMLformula'.\n",
                  "  Make sure what you mean and simplify the arguments."))
    }
    # We do not need the LHS !
    # if (length(predictor)==3) {
    #   lhs <- paste(predictor)[[2]] ## extract response, either cbind or not
    # } else lhs <- "" ## occurs for predictor of phi if dhglm ML fit.
    # REMLformula <- as.formula(paste(lhs,"~ 0")) 
    REMLformula <- ~ 0
    attr(REMLformula,"isML") <- TRUE
  } ## else do nothing: keeps input REMLformula, which may be NULL or a non-trivial formula
  # REMLformula <- .preprocess_formula(REMLformula)
  # emptyenv() insufficient for .preprocess -> ... -> assign_X.Re_objective -> .get_terms_info(formula = REMLformula, data = data, famfam = "") 
  if ( ! is.null(REMLformula)) environment(REMLformula) <- new.env(parent=baseenv())
  processed$REMLformula <- REMLformula  
}

.as_mat_wAttr <- function(X) {
  if (inherits(X,"sparseMatrix"))  { # if X.pv was found to be sparse but decorr was ultimately selected.
    # sparse X.pv will contaminaes cbind()'s and in not compatible with "decorr" method.
    # So we have converted X.pv to sparse format and we convert it back...
    # Not easy to avoi that given that:
    # * rankinfo was determined by fast sparse QR facto; 
    # * X.pv sparsity may be used in at least some subcases for choosing the "algebra". 
    #   One would have to achieve the same effect while keeping X.pv storage dense to avoid a possible back conversion
    # Overall this back-conversion is quite rare: its first occurrence was in pathological case
    # (formula with Matern() + location factor as fixed effect...)
    Xattr <- attributes(X)
    X <- as.matrix(X)
    names_lostattrs <- setdiff(names(Xattr), c(names(attributes(X)),"Dim","Dimnames","i","p","x","factors","class"))
    attributes(X)[names_lostattrs] <- Xattr[names_lostattrs] # as in .subcol_wAttr(). 
  }
  X
}

.as_spMat_wAttr <- function(X) {
  if (is.matrix(X))  { # if X.pv was found to be sparse but decorr was ultimately selected.
    # sparse X.pv will contaminaes cbind()'s and in not compatible with "decorr" method.
    # So we have converted X.pv to sparse format and we convert it back...
    # Not easy to avoi that given that:
    # * rankinfo was determined by fast sparse QR facto; 
    # * X.pv sparsity may be used in at least some subcases for choosing the "algebra". 
    #   One would have to achieve the same effect while keeping X.pv storage dense to avoid a possible back conversion
    # Overall this back-conversion is quite rare: its first occurrence was in pathological case
    # (formula with Matern() + location factor as fixed effect...)
    Xattr <- attributes(X)
    X <- as(X, "sparseMatrix")
    names_lostattrs <- setdiff(names(Xattr), c(names(attributes(X)),"dim","dimnames","class"))
    attributes(X)[names_lostattrs] <- Xattr[names_lostattrs] # as in .subcol_wAttr(). 
  }
  X
}

.init_AUGI0_ZX <- function(X.pv, vec_normIMRF, ZAlist, nrand, n_u_h, sparse_precision, as_mat) {
  if (nrand) {
    pforpv <- ncol(X.pv)
    if (as_mat) {
      X.pv <- .as_mat_wAttr(X.pv)
      AUGI0_ZX <- list2env( list(I=diag(nrow=n_u_h),ZeroBlock= matrix(0,nrow=n_u_h,ncol=pforpv), X.pv=X.pv) )
    } else {
      AUGI0_ZX <- list2env( list(I=.sparseDiagonal(n=n_u_h, shape="g"), ## to avoid repeated calls to as() through rbind2...; previously used .trDiagonal()
                                 ZeroBlock= Matrix(0,nrow=n_u_h,ncol=pforpv), X.pv=X.pv) )
      # delayedAssign("Ilarge", .trDiagonal(n=ncol(I)+ncol(X.pv), unitri = FALSE), eval.env = AUGI0_ZX, assign.env = AUGI0_ZX) # hmf: .trDiagonal  ~ 8e-4 s. (but delayedA is 500 times faster)
    } ## $ZAfix added later   and   X.pv scaled below  !!
    AUGI0_ZX$vec_normIMRF <- vec_normIMRF
    AUGI0_ZX$envir <- list2env(list(finertypes=attr(ZAlist,"exp_ranef_types"), ## to be modified later
                                    LMatrices=structure(vector("list",nrand),
                                                        is_given_by=rep("",nrand)), ## to be modified later
                                    kron_Y=structure(vector("list",nrand),
                                                     is_given_by=rep("",nrand)), ## to be modified later
                                    updateable=attr(ZAlist, "exp_ranef_types")=="(.|.)" ## to be modified later
                                    ),    
                               parent=environment(.preprocess))
    # updateable introduced here 2022/01/30. For CORR methods by .init_AUGI0_ZX(), it is updated for some correlated ranefs in .assign_geoinfo_and_LMatrices_but_ranCoefs()
    # That latter code would deserve to be checked (whether some FALSE could become TRUE) since it has been used only by aug_ZXy's .get_absdiagR_blocks() procedure .
    # For SPPREC, it is rebuilt by .init_AUGI0_ZX_envir_spprec_info()
    # For CORR methods, this paves the way for extended usage, but won't solve what prevents using permuted Cholesky for which="R_scaled_v_h_blob" in QRP_CHM. 
    if (sparse_precision) {
      AUGI0_ZX$trDiag <- ..trDiagonal(n=n_u_h) # -> ... .sparseDiagonal(n=n_u_h, shape="t")
    } else {
      if (inherits(AUGI0_ZX$ZeroBlock,"sparseMatrix")) {
        AUGI0_ZX$trDiag <- ..trDiagonal(n=n_u_h+pforpv)
        AUGI0_ZX$Zero_X <- rbind2(AUGI0_ZX$ZeroBlock, as(AUGI0_ZX$X.pv,"CsparseMatrix"))
      } else AUGI0_ZX$Zero_X <- rbind2(AUGI0_ZX$ZeroBlock, AUGI0_ZX$X.pv)
    }
    AUGI0_ZX <- .add_ZAfix_info(AUGI0_ZX, ZAlist, sparse_precision, as_mat=as_mat)
  } else AUGI0_ZX <- list(X.pv=X.pv)
  AUGI0_ZX
}

.preprocess_X_XRe_off <- function(main_terms_info, predictor, processed, X.pv, etaFix, data, objective, nobs) {
  off <- model.offset(main_terms_info$mf) ## look for offset from (ori)Formula 
  if ( ! is.null(off) ) { ## offset (ori)Formula term. Check attribute (Predictor ensures offset is not both in formula and attr)
    no_offset <- .stripOffset(predictor) ## note that .stripOffset()ing is not performed in .preprocess for the resid.formula 
    off <- pmax(log(.Machine$double.xmin),off) ## handling log(0) ## but if input off were NULL, output off would be is numeric(0) where it should remain NULL
    processed$predictor <- structure(predictor, no_offset=no_offset)
  } else processed$predictor <- structure(predictor, no_offset=predictor)
  ## Extract columns for fixed oefficients (involves offset; ! can conflict with .rankTrim() results)
  rownames(X.pv) <- NULL
  XReinput <- X.pv ## may be overwritten from etaFix$beta code before the .assign_X.Re_objective() call..
  ## reimplementation of etaFix$beta (2015/03)
  if ( length(betaFix <- etaFix$beta)>0 ) {
    namesbetafix <- names(betaFix)
    if (is.null(namesbetafix)) {
      message("The elements of etaFix$beta should be named and the names should match the column names of the design matrix.")
    }
    if (length(setdiff(namesbetafix,colnames(X.pv)))==0L) { ## if no incorrect name
      offFromEtaFix <- drop(X.pv[ ,namesbetafix,drop=FALSE] %*% betaFix) # must be vector not matrix
      namesbetavar <- setdiff(colnames(X.pv),namesbetafix)
      X.pv <- .subcol_wAttr(X.pv, j=namesbetavar, drop=FALSE)
      if (is.null(off)) {
        off <- offFromEtaFix
      } else off <- off + offFromEtaFix
      ## TRUE by default:
      if ( is.null( keepInREML <- attr(betaFix,"keepInREML") ) ||  ( ! keepInREML) ) {
        XReinput <- X.pv # If X.pv is modified, XReinput is too
      } else { # (keepInREML TRUE)=> XReinput and X.pv are now different, XReinput being the X.pv before .subcol_wAttr()ing
        if (is.null(processed$REMLformula)) processed$REMLformula <- predictor # fix 2022/12/08; 
        # without that .assign_X.Re_objective does not set up the correct design matrix for REML correction.
        # Test is by numInfo(<REML fit>, which=..."beta") not failing the gradient check and replicating the likelihood (use verbose to check that)
      }
    } else {
      stop("The names of elements of etaFix$beta should all match column names of the design matrix.")
    }
  } 
  .assign_X.Re_objective(processed, XReinput=XReinput, processed$REMLformula, data, X.pv, objective) ##] assigns processed$X.Re and processed$objective
  
  if (is.null(off)) { ## model.frame.default(formula = locform, offset = off,...) expects a vector....
    processed$off <- rep(0,nobs) ## long form expected by spaMM_glm.fit() [as by glm.fit()] and then possibly assumed by further code 
  } else {
    processed$off <- off
  }
  X.pv
}

.preprocess_arglists <- function(processed) {
  processed$u_h_v_h_from_v_h <- .def_u_h_v_h_from_v_h(processed)
  processed$updateW_ranefS <- .def_updateW_ranefS(processed)
  arglists <- list() 
  nrand <- length(processed$rand.families)
  rand_list <- vector("list", nrand) # to avoid repeating this in .initialize_v_h(). Same idea below.
  repNAnrand <- rep(NA_real_, nrand) 
  names(rand_list) <- names(repNAnrand) <- seq_len(nrand)
  arglists$rand_list <- rand_list
  arglists$repNAnrand <- repNAnrand
  arglists
}

.eval_check_vec_n_u_h <- function(ZAlist, nobs, processed) {
  vec_n_u_h <- unlist(lapply(ZAlist,ncol)) ## nb cols each design matrix = nb realizations each ranef
  if (min(nobs/vec_n_u_h)<2) {
    if (processed$family$family == "binomial" && processed$bin_all_or_none && processed$HL[1]==1L) {
      if ( ! identical(spaMM.getOption("PQL_warned"),TRUE)) {
        message("Fits using Laplace approximation may diverge for (nearly) all-or-none binomial data:\n check PQL or PQL/L methods in that case.")
        .spaMM.data$options$PQL_warned <- TRUE
      }
    }
  }
  vec_n_u_h
}

.vecdisneeded <- function(pforpv, family, processed) {
  if (pforpv) {
    coef12needed <- ! attr(processed$models,"unit_Hobs_weights") 
    coef3needed <- any(attr(processed$rand.families,"unique.psi_M")!=0) # psi_M happening to be 0 only for gaussian(identity) ranef
    vecdisneeded <- c( coef12needed, coef12needed, coef3needed )
  } else vecdisneeded <- rep(FALSE,3L)
  vecdisneeded
}

.need_obsAlgo <- function(HLmethod, family, canonicalLink=family$flags$canonicalLink, 
                          nrand, opt=.spaMM.data$options$obsInfo) {
  # not GLM -> always obsInfo
  # GLM fam with canonical link -> this is indifferent. here obsInfo set to false, 
  if ( ! family$flags$exp ) { # no expected-info capacity; not GLM
    obsInfo <- TRUE
  } else if (canonicalLink) {
    obsInfo <- 0L # let's try to have two kinds of false; used in how.HLfit()
  } else if (length(HLmethod)==2L) {
    obsInfo <- switch( HLmethod[[2L]],
                       "obs" = TRUE,
                       "exp" = FALSE,
                       stop("Unhandled second specifier in 'method' argument"))
  } else {
    obsInfo <- (
      HLmethod[[1L]]!="SEM" && 
        opt ## had previously : && (nrand || opt>1L) # so that default was still "exp" in fixed-effect models
    ) # use .spaMM.data's default method is not SEM
  }
  obsInfo
}
