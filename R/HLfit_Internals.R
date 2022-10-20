
## declared "arglist" to print a clean summary instead of very long list
## the print(summary()) is required if the "arglist" is a member (eg as.list(<call>));
## summary alone would print nothing
print.arglist <- function(x,...) {
  if (inherits(x,"environment")) { ## not clear that this case should occur, but it occurs notably in HLCor.args of refit
    str(x) ## bc summary(<env>) is a bug
  } else print(summary(x,...))
}
##

.Dvec_times_m_Matrix <- function(Dvec, X) {
  ## Should be consistent with R's diag(Dvec) %*% X that keeps only the X colnames:
  if (inherits(X,"Matrix")) {
    return(.Dvec_times_Matrix(Dvec=Dvec, X=X))
  } else return(.Dvec_times_matrix(Dvec=Dvec, X=X))
}

.Dvec_times_matrix <- function(Dvec, X) { ## for *m*atrix input
  if (nrow(X)!=length(Dvec)) {
    stop("nrow(X)!=length(Dvec) ") ## fatal error for eigen code...
  } else if (ncol(X)==0L) {
    return(X)
  } else {
    # "Error in .Rcpp_sweepZ1W(X, Dvec) : Wrong R type for mapped matrix" may signal an integer matrix...
    res <- .Rcpp_sweepZ1W(X,Dvec) ## Rcpp (fast !) version of sweep ( MARGIN=1L ) which is also X * Dvec
    ## Consistent with R's diag(Dvec) %*% X that keeps only the X colnames:
    colnames(res) <- colnames(X)
    return(res)
  }
}

.m_Matrix_times_Dvec <- function(X, Dvec) {
  if (inherits(X,"ZAXlist")) {
    return(.ZAXlist_times_Dvec(X=X, Dvec=Dvec))
  } else if (inherits(X,"Kronfacto")) {
    return(.m_Matrix_times_Dvec(X=X@BLOB$long, Dvec=Dvec))
  } else if (inherits(X,"Matrix")) {
    return(.Matrix_times_Dvec(X=X, Dvec=Dvec))
  } else return(sweep(X, MARGIN=2L, Dvec, `*`))
}

.Matrix_times_Dvec <- function(X,Dvec, 
                               keep_dsC=TRUE, # subcase of check_dsC...
                               check_dsC=TRUE # to be able to bypass the check in special cases
                               ) { # try to keep the type except when input is dsC and output is presumably not symmetric.
  if (inherits(X,"ddiMatrix")) {
    if (X@diag=="U") { ## diag + unitary => identity matrix
      X <- Diagonal(x=Dvec)
    } else X@x <- X@x * Dvec ## raw diag matrix
  } else if (inherits(X,c("dgCMatrix", "dtCMatrix"))) {
    X@x <- X@x * rep(Dvec,diff(X@p))    
    ## a triangular matrix (dtC) with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
    if ( methods::.hasSlot(X, "diag") && X@diag=="U") Matrix::diag(X) <- Dvec
  } else if (inherits(X,"dsCMatrix")) { 
    if (check_dsC) {
      if (generally_dgC <- (length(unique(Dvec))>1L)) {
        # The correct result is presumably NOT symmetric so direct manipulation of slots of X-dsCMatrix is not appropriate. 
        # We go through dgC (correct result in all cases) and may try converting back to dsC afterwards.
        if (.spaMM.data$options$Matrix_old) {
          X <- as(X,"dgCMatrix")
        } else X <- as(X,"generalMatrix")
      }
      X@x <- X@x * rep(Dvec,diff(X@p))    
      if (generally_dgC) {
        if (keep_dsC) {
          warning("Suspect call: return value presumably not symmetric, but dsCMatrix requested.")
          X <- as(X,"symmetricMatrix")
        } # else warning(".Matrix_times_Dvec(<dsCMatrix>,...) returns a dgCMatrix.") # we explicitly used non-default keep_dsC, so warning not useful
      } # else we already have a correct dsC
    } else { # we assume that keeping dsC is OK: risky, non-default, but for precisionBlocks it's really hat we want (so no warning).
      # warning("call to .Matrix_times_Dvec(<dsCMatrix>, ..., check_dsC=FALSE) is possibly inefficient code.")
      X@x <- X@x * rep(Dvec,diff(X@p))    
    }
  } else if (inherits(X,c("dgTMatrix"))) {
    X@x <- X@x * Dvec[X@j+1L]    
  } else {
    warning("inefficient code in .make_Xscal or .Matrix_times_Dvec") ## eg dgeMatrix: dense matrix in the S4 Matrix representation
    # but also dtC for which the warning is not warranted
    X <- X %*% Diagonal(x=Dvec)
  } 
  if (methods::.hasSlot(X, "factors")) X@factors <- list()
  return(X)
}

if (FALSE) {
  .m_Matrix_urblock_times_Dvec <- function(X, Dvec) { ## hmmm but it's in the other corner...
    max_row <- length(Dvec)
    prev_col <- ncol(X)-length(Dvec)
    if (inherits(X,"ddiMatrix")) {
      warning("Suspect call to .Matrix_ZAL_ulblock_times_Dvec()") # the matrix is square hence X.pv not NULL, and ZAL is a zero block !
    } else if (inherits(X,c("dgCMatrix","dsCMatrix", "dtCMatrix"))) { 
      which_i_affected_rows <- X@i<max_row
      col_indices <- rep(1L:(ncol(X)),diff(X@p))
      col_indices <- col_indices-prev_col
      which_affected_col_inds <- col_indices>0L
      which_affected_x <- (which_i_affected_rows & which_affected_col_inds)
      X@x[which_affected_x ] <- X@x[which_affected_x]* Dvec[col_indices[which_affected_x]]    
      ## A triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
      #   but we dont expec't a triangular matrix here.
      if ( methods::.hasSlot(X, "diag") && X@diag=="U") warning("Suspect call to .Matrix_ZAL_ulblock_times_Dvec()")
    } else { # dense matrix either _m_atrix or eg dgeMatrix in the S4 Matrix representation
      whichrows <- 1L:max_row
      whichcols <- (prev_col+1L):ncol(X)
      X[whichrows,whichcols] <- .m_Matrix_times_Dvec(X[whichrows,whichcols],Dvec)
    } 
    return(X)
  }
  
  .m_Matrix_llblock_times_Dvec <- function(X, Dvec) {
    n_u_h <- length(Dvec)
    if (inherits(X,"ddiMatrix")) {
      warning("Suspect call to .Matrix_ZAL_ulblock_times_Dvec()") # the matrix is square hence X.pv not NULL, and ZAL is a zero block !
    } else if (inherits(X,c("dgCMatrix","dsCMatrix", "dtCMatrix"))) { 
      which_i_affected_rows <- X@i>(n_u_h-1L) # couldn't we use a precomputed which_i_affected_rows as elesewhere?
      col_indices <- rep(1L:(ncol(X)),diff(X@p)) ## from 1!
      which_affected_col_inds <- col_indices<=n_u_h
      which_affected_x <- (which_i_affected_rows & which_affected_col_inds)
      X@x[which_affected_x ] <- X@x[which_affected_x]* Dvec[col_indices[which_affected_x]]    
      ## A triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
      #   but we dont expec't a triangular matrix here.
      if ( methods::.hasSlot(X, "diag") && X@diag=="U") warning("Suspect call to .Matrix_ZAL_ulblock_times_Dvec()")
    } else { # dense matrix either _m_atrix or eg dgeMatrix in the S4 Matrix representation
      whichrows <- (n_u_h+1L):nrow(X)
      whichcols <- seq_len(n_u_h)
      X[whichrows,whichcols] <- .m_Matrix_times_Dvec(X[whichrows,whichcols],Dvec)
    } 
    return(X)
  }
} ## end if (FALSE)

.Dvec_times_Matrix <- function(Dvec,X) {
  if (inherits(X,"ddiMatrix")) {
    if (X@diag=="U") { ## diag + unitary => identity matrix
      X <- Diagonal(x=Dvec)
    } else X@x <- X@x * Dvec ## raw diag matrix
  } else if (inherits(X,c("dgCMatrix","dsCMatrix", "dtCMatrix"))) { 
    X@x <- X@x * Dvec[X@i+1L] 
    ## a triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
    if ( methods::.hasSlot(X, "diag") && X@diag=="U") Matrix::diag(X) <- Dvec
  } else if (inherits(X,c("dgeMatrix"))) { ## dense Matrix
    X <- X * Dvec 
  } else {
    # dgeMatrix occurs in Matern fitted by sparse corr methods (occurs if other ranef sparsifies stuff !?)
    warning("Possibly inefficient code in .Dvec_times_Matrix") ## eg dgeMatrix: dense matrix in the S4 Matrix representation
    ## warning if a dgeMatrix has been created = possibly inefficient code
    ## Other matrix formats are possibly not excluded
    ## if it is really a dgeMatrix then Dvec * X appears correct. (FIXME double check and implement ?)
    X <- Diagonal(x=Dvec) %*% X
  } 
  if (methods::.hasSlot(X, "factors")) X@factors <- list()
  return(X)
}

.Dvec_times_Matrix_lower_block <- function(Dvec,X,min_row) { # min_row from 0
  if (inherits(X,"ddiMatrix")) {
    if (X@diag=="U") { ## diag + unitary => identity matrix
      X <- Diagonal(x=c(rep(1,min_row),Dvec))
    } else X@x <- X@x * c(rep(1,min_row),Dvec) ## raw diag matrix
  } else if (inherits(X,c("dgCMatrix", "dtCMatrix"))) { 
    Xi <- X@i
    which_i_affected_rows <- Xi>(min_row-1L) 
    X@x[which_i_affected_rows] <- X@x[which_i_affected_rows]*Dvec[Xi[which_i_affected_rows]-min_row+1L]
    ## a triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
    if ( methods::.hasSlot(X, "diag") && X@diag=="U") Matrix::diag(X) <- c(rep(1,min_row),Dvec)
  } else {
    warning("inefficient code in .Dvec_times_Matrix_lower_block")
    X <- Diagonal(x=c(rep(1,min_row),Dvec)) %*% X
  } 
  return(X)
}

.Dvec_times_m_Matrix_upper_block <- function(Dvec,X) {
  more_rows <- nrow(X)-length(Dvec)
  if (inherits(X,"Matrix")) {
    if (inherits(X,"ddiMatrix")) {
      if (X@diag=="U") { ## diag + unitary => identity matrix
        X <- Diagonal(x=c(Dvec,rep(1,more_rows)))
      } else X@x <- X@x * c(Dvec,rep(1,more_rows)) ## raw diag matrix
    } else if (inherits(X,c("dgCMatrix", "dtCMatrix"))) { 
      which_i_affected_rows <- X@i<length(Dvec) 
      X@x[which_i_affected_rows] <- X@x[which_i_affected_rows]*Dvec[X@i[which_i_affected_rows]+1L]
      ## a triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
      if ( methods::.hasSlot(X, "diag") && X@diag=="U") Matrix::diag(X) <- c(Dvec,rep(1,more_rows))
    } else {
      warning("inefficient code in .Dvec_times_m_Matrix_upper_block")
      X <- Diagonal(x=c(Dvec,rep(1,more_rows))) %*% X
    } 
  } else X <- diag(x=c(Dvec,rep(1,more_rows))) %*% X
  return(X)
}





## the following fns try to keep the input class in output, but are called with dense matrices (except first tested case).
# les Matrix::(t)crossprod  paraissent (parfois au moins) remarquablement inefficaces !!
# idem pour Diagonal()
.ZWZtwrapper <- function(ZAL,w,as_sym=TRUE) { 
  if (inherits(ZAL,"Matrix")) {
    if (inherits(ZAL,"ddiMatrix")) { ## if diagonal... returns dsC
      if ( ZAL@diag=="U") {
        return(.symDiagonal(x=w)) ## diagonal and unitary => identity. But not ddi: ""returns an object of class dsCMatrix or lsCMatrix"
      } else return(.symDiagonal(x = w * diag(ZAL)^2)) ## this case may never occur
    } else {
      sqrtW <- sqrt(w)
      ZW <- .Matrix_times_Dvec(ZAL,sqrtW)
      return(.tcrossprod( ZW, as_sym=as_sym)) ## dsCMatrix if as_sym=TRUE (default) if input matrix was dgC. if it was dge, result is dpo
    }
  } else if (inherits(ZAL,"ZAXlist")) {
    sqrtW <- sqrt(w)
    ZW <- .m_Matrix_times_Dvec(ZAL,sqrtW)
    return(.tcrossprod( ZW,as_sym=as_sym)) ## dsCMatrix if as_sym=TRUE (default)
  } else return(.ZWZt(ZAL,w))
} ## so seems always dsC

#.ZWZtbase <- function(ZAL,w) tcrossprod(sweep(ZAL,MARGIN=2L,w,`*`),ZAL) # slower than Rcpp version, with no obvious benefits
#  !  not the same as     ZAL %*% sweep(ZAL,MARGIN=1L,w,`*`) 

.ZtWZwrapper <- function(ZAL,w, allow_as_mat=TRUE, sqrtW=sqrt(w)) { ## used in several contexts
  if (inherits(ZAL,"ZAXlist")) {
    tcrossfac <- .crossprod(ZAL,Diagonal(x=sqrtW)) # t() %*% sqrtw as dsCMatrix # by ad-hoc 
    return(.tcrossprod(tcrossfac))  #
  } else if (inherits(ZAL,"Matrix")) {
    # but Diagonal() is slow and it's better to produce matrices of symmetric type by construction so that forceSymmetric is not needed
    # Matrix::.symDiagonal better for addition with another symmetric matrix, see its doc.
    if (inherits(ZAL,"ddiMatrix")) { ## if diagonal
      if ( ZAL@diag=="U") {
        return(.symDiagonal(x=w)) ## diagonal and unitary => identity
      } else return(.symDiagonal(x = w * diag(ZAL)^2)) ## this case may never occur
    } else  {
      DZAL <- ZAL
      DZAL@x <- DZAL@x * sqrtW[DZAL@i+1L] ## W Z
      ## a triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
      if (  methods::.hasSlot(DZAL, "diag") &&  DZAL@diag=="U") Matrix::diag(DZAL) <- sqrtW
      return(.crossprod(DZAL,allow_as_mat=allow_as_mat))  
    }
  } else return(.ZtWZ(ZAL,w)) ## but this uses sqrt() hence assumes w>0
} # Matrix result is always of symmetric type (dsC or ddi)

.sapply <- function(liste, FUN, ...) {
  resu <- lapply(X = liste, FUN = FUN, ...)
  .unlist(resu)
}

.sapply_null_or_NA <- function(liste) {
  resu <- logical(length(liste))
  for (it in seq_along(liste)) resu[it] <- length(liste[[it]])<2L # length is 0 fo NULL, 1 for NA; 
  # is.null(liste[[it]]) || is.na(liste[[it]]) is inappropriate for any vector 1,0,1 or Na,0,NA
  resu
}

.sapply_anyNA <- function(liste) {
  resu <- logical(length(liste))
  for (it in seq_along(liste)) resu[it] <- anyNA(liste[[it]])
  resu
}

.get_rC_longLv_phantom_map <- function(envir, Xi_ncol) {
  if (is.null(envir$rC_longLv_phantom_map)) {
    #
    map_X <- matrix(seq(Xi_ncol^2),nrow=Xi_ncol,ncol=Xi_ncol)
    phantom_X <- upper.tri(map_X, diag=TRUE)
    map_X <- as(map_X*phantom_X,"dtCMatrix")
    #
    phantom_Y <- as(envir$Lunique,"ltCMatrix") # sparse logical matrix...
    LHS <- .makelong(map_X, kron_Y=phantom_Y) 
    #
    phantom_X <- as(phantom_X,"ltCMatrix")
    envir$rC_longLv_phantom_map <- list(
      phantom_X=phantom_X,
      LHS_x= as.integer(LHS@x), # OK bc with a tolerance of .Machine$double.eps
      RHS=.makelong(phantom_X, kron_Y=envir$Lunique) 
    )
  }
  envir$rC_longLv_phantom_map
}


.wrap_makelong_outer_composite <- function(processed, rt, latentL_blob) { # called by process_ranCoefs()
  # template <- processed$ranCoefs_blob$longLv_templates[[rt]]
  cov_info_mat <- processed$corr_info$cov_info_mats[[rt]]
  if (use_corrMatrix_fast_code <- TRUE) { # ____TAG___ OK bc unpermuted Chol for corrMatrix
    # only makelongs and ranCoefs's compact-matrix operations 
    if (FALSE) { # safe code using slow, repetitive Matrix::kronecker() cals 
      rC_blob_longLvMatrix <- .makelong(solve(t(latentL_blob$compactchol_Q)),#longsize=ncol(processed$ZAlist[[rt]]),
                                        # template=template,
                                        kron_Y=attr(cov_info_mat,"blob")$Lunique) # $Lunique is upper tri
    } else { # code with precomputation of two "phantoms" of the final kronecker product kronecker(X,Y),
      # the 'RHS' with indices of the elements of X that will factor with elements of Y 
      # the 'LHS' with copies of the actual Y (say copies of Lunique)
      # then kronecker product = RHS@x * X[LHS@x] is the only computation for each new X
      phantasmapic <- .get_rC_longLv_phantom_map(envir=attr(cov_info_mat,"blob"), Xi_ncol=ncol(latentL_blob$compactchol_Q))
      kron_X <- .rawsolve(t(latentL_blob$compactchol_Q)) 
      rC_blob_longLvMatrix <- phantasmapic$RHS
      kron_X <- as.matrix(kron_X) # as.matrix() before subsetting.
      rC_blob_longLvMatrix@x <- rC_blob_longLvMatrix@x * kron_X[phantasmapic$LHS_x] 
      # if (any(kron_X[phantasmapic$phantom_X]==0)) rC_blob_longLvMatrix <- drop0(rC_blob_longLvMatrix) 
      ## : the $phantom_X serves here to test whether the kron_X factor is even sparser than phantom_X allows 
      ## but costs may outweight benefits
    }
  } else { # allowing permuted cholesky:
    cov_long <- .makelong(tcrossprod(latentL_blob$compactchol_Q),#longsize=ncol(processed$ZAlist[[rt]]), 
                          # template=template,
                          kron_Y=cov_info_mat$matrix # should be dsC 
    )
    # calling/updating Cholesky for any new compact chol Q:
    long_Q_CHMfactor <- 
      .get_Q_CHMfactor__assign_template__perm_ZA(processed$AUGI0_ZX$envir, rt, processed, cov_long) 
    #
    rC_blob_longLvMatrix <- .Lunique_info_from_Q_CHM(
      processed=processed, rd=rt, Q_CHMfactor=long_Q_CHMfactor, 
      # with default corr_type 
      # with default spatial_term=attr(processed$ZAlist,"exp_spatial_terms")[[rt]], 
      type="from_Q_CHMfactor") # result is either the (long)Q_CHMfactor or has the (long) Q_CHMfactor as attribute
  }
  rC_blob_longLvMatrix
}


.process_ranCoefs <- function(processed, ranCoefs, trRanCoefs, use_tri_CORREL) {
  ranCoefs_blob <- processed$ranCoefs_blob
  ZAlist <- processed$ZAlist 
  exp_ranef_types <- attr(ZAlist,"exp_ranef_types")
  # For composite ranef:
  # exp_ranef_types[] is the is the <corrMatrix|...> type (vs (.|.) for a pure ranCoefs)
  # processed$corr_info$corr_types[] is the <corrMatrix|...> type ( vs NA for a pure ranCoef)
  # finertypes[] is 'ranCoefs'
  # corr.model[] is "random-coef"
  nrand <- length(ZAlist)
  cum_n_u_h <- processed$cum_n_u_h
  if (nrand && is.null(ranCoefs_blob)) { ## in call from .preprocess(), allows simplified user input as numeric vector
    Xi_cols <- attr(ZAlist, "Xi_cols")
    isRandomSlope <- Xi_cols>1L ## FIXME seems oK for later code but semantically sloppy, cf (X-1|id) terms
    hasRandomSlope <- any(isRandomSlope)
    longLv_templates <- vector("list", length(Xi_cols))
    #constraints <- rep(list(NULL), length(Xi_cols))
    if (hasRandomSlope) {
      for(rd in seq_len(length(Xi_cols))) {
        if ((Xi_ncol <- Xi_cols[rd])>1L) {
          templa <- matrix(seq(Xi_ncol^2),nrow=Xi_ncol,ncol=Xi_ncol)
          longLv_templates[[rd]] <- .makelong(templa, ncol(ZAlist[[rd]]), ## templa is (compact) dense array of indices, here serving to build the sparse long template 
                                              kron_Y=NULL)         
        } # in corrMatrix(mv(1,2)), kron_Y mauy be NULL in the two .preprocess() calls, then non-NULL in the .merge... call  
      }
    }
    ranCoefs_blob <- list(isRandomSlope=isRandomSlope, 
                          is_composite=(isRandomSlope & exp_ranef_types != "(.|.)"), 
                          is_set=(isRandomSlope & FALSE), # vector
                          new_compos_info=(isRandomSlope & FALSE), # vector 
                          longLv_templates=longLv_templates, one_time_check_not_done=TRUE)
    
    if (hasRandomSlope && ! is.null(ranCoefs)) { ## ranCoefs in preprocess: fixed values
      if (is.numeric(ranCoefs)) {
        if (sum(isRandomSlope)==1L) {
          ranCoefs <- vector("list", nrand)
          ranCoefs[[isRandomSlope]] <- ranCoefs
        } else stop("Cannot match numeric ranCoefs to a single random-coefficient term.")
      } else if (length(ranCoefs) < nrand ){
        rand_list <- processed$reserve$rand_list
        rand_list[names(ranCoefs)] <- ranCoefs
        ranCoefs <- rand_list
      }
      not_constr_or_set <- .sapply_null_or_NA(ranCoefs) # tests all vectors => result is vector
      newly_set <- ! (not_constr_or_set | .sapply_anyNA(ranCoefs))  
      # newly_set contains some TRUE (unless the user provided an uninformative ranCoefs) =>  no return yet ! 
      ranCoefs_blob$is_set <- newly_set
    } else {
      ranCoefs_blob$is_diag_family <- rep(FALSE, length(Xi_cols))
      return(ranCoefs_blob) # still preprocess case
    }
  } else if (is.null(ranCoefs) || # after preprocess: occurs if nrand=0
             ! any(isRandomSlope <- ranCoefs_blob$isRandomSlope)) {
    return(ranCoefs_blob)
  } else { ## both ranCoefs_blob and ranCoefs non-NULL, && any(isRandomSlope)
    ## assume 'incomplete' list
    nrand <- length(isRandomSlope)
    if (length(ranCoefs) < nrand ){
      rand_list <- processed$reserve$rand_list
      rand_list[names(ranCoefs)] <- ranCoefs
      ranCoefs <- rand_list
    }
    newly_set <- ! .sapply_null_or_NA(ranCoefs) 
    ##### Handling of partially fixed ranCoefs:
    # in post-processing, we reach here from various fitting functions (fitme, HLfit...)
    # If from fitme, ranCoefs should no longer have NA
    # If from HLfit, ranCoefs still have its NA's
    # But implementing partial constraints in inner optim is complicated:
    # One would have to hack the compactcovmat from .makeCovEst1() -> .objfn() -> .calc_cov_from_trRancoef(),
    #   including the compactcovmat's 'chol_crossfac' attribute. Possible complications without a clear interest => not a priority
    # One would have also to provide access to the constraint in .makeCovEst1. This is apparently not currently the case
    #   although this might be done by adding the constraint in processed$ranCoefs_blob
    # Hence:
    if (ranCoefs_blob$one_time_check_not_done && any( newly_set &  .sapply_anyNA(ranCoefs))) {
      stop("Partial constraints on ranCoefs are not handled by this function. Use fitme() instead.")
    }
    ranCoefs_blob$one_time_check_not_done <- FALSE
    #####
    newly_set <- newly_set & ( ! ranCoefs_blob$is_set)
    if ( any(newly_set | ranCoefs_blob$new_compos_info)) {
      Xi_cols <- attr(ZAlist, "Xi_cols")
      ##  & continue final block of code untilfinal return(), skipping immediate return 
    } else return(ranCoefs_blob)
  }
  ## newly_set or ranCoefs_blob$new_compos_info contains some TRUE
  spprecBool <- processed$is_spprec
  nrand <- length(ZAlist)
  if (is.null(LMatrices <- ranCoefs_blob$LMatrices)) LMatrices <- vector("list", nrand)
  if (is.null(lambda_est <- ranCoefs_blob$lambda_est)) lambda_est <- numeric(cum_n_u_h[nrand+1L])
  min_lam <- .spaMM.data$options$tol_ranCoefs_inner["tol"]
  ## fitme/fitmv
  # For an composite ranef:
  # if ranCoefs is fixed: the call from .[merge_]processed has newly_set=TRUE (the ranCoefs ARE newly available) 
  #                       and new_compos_info=FALSE (never TRUE at this point). The latentL_blob must then be created 
  #                       hence the condition below is (newly_set | ...) ( using '&' leads to a NaN loglik)
  #                       The first call from HLfit_body builds the full info (newly_set=FALSE but 
  #                         new_compos_info=TRUE).
  # if ranCoefs is variable: the call from .[merge_]processed has both tests FALSE
  #                       the call from  HLfit_body has both tests TRUE
  # => In neither case new_compos_info is true before newly_set
  ## Now if we try refit = TRUE... 
  for (rt in which(isRandomSlope & (newly_set | 
                                    ( ranCoefs_blob$is_set & ranCoefs_blob$new_compos_info)))) {
    Xi_ncol <- Xi_cols[rt]
    
    compactcovmat <- .C_calc_cov_from_ranCoef(ranCoef=ranCoefs[[rt]], Xi_ncol=Xi_ncol) ## hfff made a test in test-poly relative as Eigen remains light on numerical precision
    latentL_blob <- .calc_latentL(compactcovmat, use_tri_CORREL=use_tri_CORREL, spprecBool=spprecBool)
    #  with(latentL_blob,.ZWZt(design_u,d)) = compactcovmat hence design_u is more of a tcrossprod factor
    ## we have a repres in terms of ZAL and of a diag matrix of variances; the latter affects only the hlik computation
    #
    if (ranCoefs_blob$new_compos_info[rt]) {
      if (spprecBool) { 
        if (TRUE) {
          rC_blob_longLvMatrix <- .def_Kronfacto(lhs=.rawsolve(t(latentL_blob$compactchol_Q)), 
                                                 rhs=attr(processed$corr_info$cov_info_mats[[rt]],"blob")$Lunique)
        } else rC_blob_longLvMatrix <- .wrap_makelong_outer_composite(processed, rt=rt, latentL_blob=latentL_blob) 
      } else {
        rC_blob_longLvMatrix <- .makelong(latentL_blob$design_u,longsize=ncol(ZAlist[[rt]]),# the variances are taken out in $d
                                          template=ranCoefs_blob$longLv_templates[[rt]],
                                          kron_Y=attr(processed$corr_info$cov_info_mats[[rt]],"blob")$Lunique) # CORREL algos
      }
      ranCoefs_blob$new_compos_info[rt] <- FALSE # might be reset to TRUE by some external code
    } else if (newly_set[rt]) {
      rC_blob_longLvMatrix <- .makelong(latentL_blob$design_u,longsize=ncol(ZAlist[[rt]]),# the variances are taken out in $d
                                        template=ranCoefs_blob$longLv_templates[[rt]],
                                        kron_Y=NULL,
                                        drop0template = FALSE ) # ___F I X M E___ can we save time in ZAL <- ... by using a Kronfacto even when kron_Y is NULL? not clear 
    } 
    attr(rC_blob_longLvMatrix,"trRancoef") <- .ranCoefsFn(ranCoefs[[rt]], rC_transf=.spaMM.data$options$rC_transf_inner) # used for init_trRancoef in .MakeCovEst()
    attr(rC_blob_longLvMatrix,"Ltype") <- "cov"
    
    attr(rC_blob_longLvMatrix,"latentL_blob") <- latentL_blob ## kept for updating in next iteration, strucList, predVar, and output
    ranef_info <- attr(ZAlist,"exp_ranef_strings")[rt]
    attr(ranef_info, "type") <- exp_ranef_types[rt]## type is "(.|.)" if LMatrix is for pure random slope ## ajout 2015/06
    #### attr(rC_blob_longLvMatrix,"ranefs") <-  ranef_info
    attr(rC_blob_longLvMatrix, "corr.model") <- "random-coef"
    LMatrices[[rt]] <- rC_blob_longLvMatrix
    n_levels <- ncol(ZAlist[[rt]])/Xi_ncol 
    u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
    lambda_est[u.range] <- 1
    lambda_est[lambda_est<min_lam] <- min_lam ## arbitrarily small eigenvalue is possible for corr=+/-1 even for 'large' parvec
  }
  # with ranCoefs_blob$new_compos_info we can reach with is_set TRUE and newly_set FALSE
  # so previous code ranCoefs_blob$is_set <- newly_set is clearly wrong
  # but it was also suspect for mixed types of ranCoefs?
  ranCoefs_blob$is_set <- (ranCoefs_blob$is_set | newly_set) 
  ranCoefs_blob$LMatrices <- LMatrices
  ranCoefs_blob$lambda_est <- lambda_est
  return(ranCoefs_blob)
}

.make_long_lambda <- function(d, n_levels, Xi_ncol=length(d)) { ## n_levels as given by rhs of ranef term
  lambda <- rep(d,rep(n_levels,Xi_ncol)) 
  return(lambda)
}


.calcRanefPars <- function(HLfit_corrPars=list(),
                          lev_lambda,
                          ranefEstargs, # mostly list of arg for .makeCovEst1
                          ranCoefs_blob,
                          lam_fix_or_outer_or_NA, ## distinctly used for CARdispGammaGLM
                          rand.families,
                          psi_M,
                          verbose,
                          control,
                          iter, ## ajustement gamma(identity...)
                          maxLambda
) {
  ## Build pseudo response for lambda GLM/HGLM
  glm_lambda <- NULL
  next_LMatrices <- prev_LMatrices <- ranefEstargs$prev_LMatrices
  if ( ! is.null(prev_LMatrices) && ! is.list(prev_LMatrices)) prev_LMatrices <- list(dummyid=prev_LMatrices)
  ranefs <- attr(ranefEstargs$ZAlist,"exp_ranef_strings")
  nrand <- length(ranefEstargs$ZAlist) ## Cf notes du 3/6/2015
  done <- rep(FALSE,nrand)
  anyVanishLam <- NULL
  u_h <- ranefEstargs$u_h
  cum_n_u_h <- ranefEstargs$cum_n_u_h
  resp_lambda <- matrix(0,cum_n_u_h[nrand+1L],1L)
  next_lambda_est <- numeric(length(u_h)) 
  #########################
  isRandomSlope <- ranCoefs_blob$isRandomSlope
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  if (any(isRandomSlope)) {
    is_set <- ranCoefs_blob$is_set
    done[is_set] <- TRUE
    var_ranCoefs <- ( isRandomSlope & ! is_set ) 
    if (any(var_ranCoefs)) {
      ## handling correlation in random slope models # slmt pr gaussian ranefs, verif dans preprocess
      ranefEstargs$var_ranCoefs <- var_ranCoefs 
      LMatricesBlob <- do.call(".makeCovEst1", ranefEstargs)  ## ranCoefs estimation
      HLfit_corrPars$info_for_conv_rC <- LMatricesBlob$info_for_conv_rC
      next_LMatrices[var_ranCoefs] <- LMatricesBlob$updated_LMatrices[var_ranCoefs] # update next_LMatrices only if (any(var_ranCoefs))
      ## : updated_LMatrices is list of matrices where only random-slope elements are non NULL
      done[ var_ranCoefs ] <- TRUE
    }
    for (it in which(isRandomSlope)) { 
      u.range <- (cum_n_u_h[it]+1L):cum_n_u_h[it+1L]
      if (var_ranCoefs[it]) {
        next_lambda_est[u.range] <- LMatricesBlob$next_lambda_est[u.range] 
      } else next_lambda_est[u.range] <- ranCoefs_blob$lambda_est[u.range]
      ## resp_lambda only for testing convergence: 
      resp_lambda[u.range] <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=1) ## must give d1 in table p 989 de LeeN01
    }
  } else { ## only 'declarations' for all further code
    HLfit_corrPars$info_for_conv_rC <- NULL
  }
  ### next the other LMatrix models
  
  for (rd in (seq_len(nrand)[ ! done]) ) {
    char_rd <- as.character(rd)
    if ( ! is.null(HLfit_corrPars[[char_rd]][["rho"]])) {
      adjd <- attr(ranefEstargs$processed$corr_info$adjMatrices[[rd]],"adjd") ## may be NULL
      if (is.null(adjd)) stop("is.null(adjd)")
      locdf <- data.frame(adjd=adjd) ## $adjd, not $d which is (1/(1-rho * $adjd)): adj, not corr
      u.range <- (cum_n_u_h[rd]+1L):cum_n_u_h[rd+1L]
      locdf$resp <- resp_lambda[u.range] <- u_h[u.range]^2
      ## here CAR allows REML contrary to the SEM CAR, hence leverages
      glm_lambda <- .calc_CARdispGammaGLM(data=locdf, lambda_rd=lam_fix_or_outer_or_NA[rd], lev=lev_lambda[u.range],control=control)
      attr(glm_lambda,"whichrand") <- rd
      next_lambda_est[u.range] <- fitted(glm_lambda) ## prediction of heteroscedastic variances
      coeffs <- coefficients(glm_lambda)
      if (is.na(lam_fix_or_outer_or_NA[rd])) { 
        rho <- - coeffs[["adjd"]]/ coeffs[["(Intercept)"]]
      } else {
        rho <- - coeffs[1]*lam_fix_or_outer_or_NA
      }
      HLfit_corrPars[[char_rd]]$rho <- rho 
      done[rd] <- TRUE
    }
  }
  ### next the (no L matrices) or (Matern model and other fixed L matrix cases)
  for (it in which( ! done )) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    if (is.na(unique_lambda <- lam_fix_or_outer_or_NA[it])) {
      prior_lam_fac <- rand.families[[it]]$prior_lam_fac # prior_lam_fac is the 'design' for non-ranCoef (wei-1|.)
      if (is.null(prior_lam_fac)) wt <- 1 else wt <- 1/prior_lam_fac
      safe_dev.res_it <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=wt) ## must give d1 in table p 989 de LeeN01
      denom <- 1-lev_lambda[u.range]
      unique_lambda <- sum(safe_dev.res_it)/sum(denom) ## NOT in linkscale 
      if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { ## Gamma(identity)
        # unique_lambda <- pmin(unique_lambda,1-1/(2^(iter+1)))  ## impose lambda<1 dans ce cas 
        while(unique_lambda>1) {
          localmax <- max(safe_dev.res_it)*0.9999/unique_lambda
          safe_dev.res_it <- pmin(safe_dev.res_it,localmax) # to upper-bound the lambda near max lambda=1
          unique_lambda <- sum(safe_dev.res_it)/sum(denom) ## NOT in linkscale 
        }
      } else {
        if (unique_lambda<(minLambda <- .spaMM.data$options$minLambda)) { # minLambda indep of family link
          # Different type of correction to avoid infinite loop when min dev.res=0...
          safe_dev.res_it <- safe_dev.res_it+denom*(minLambda-unique_lambda) # to lower-bound the lambda near 1e-8
          unique_lambda <- sum(safe_dev.res_it)/sum(denom) ## NOT in linkscale 
          anyVanishLam <- TRUE
        } else while (unique_lambda>maxLambda[it]) { # processed$maxLambda depending on response-family link
          localmax <- max(safe_dev.res_it)*0.9999*maxLambda[it]/unique_lambda
          safe_dev.res_it <- pmin(safe_dev.res_it,localmax) # to upper-bound the lambda near max lambda=1
          unique_lambda <- sum(safe_dev.res_it)/sum(denom) ## NOT in linkscale 
        }
      }
      resp_lambda[u.range] <- safe_dev.res_it # In all cases, save the final corrected dev.res for use in Gamma GLM so that the latter gives the same lambda 
    } 
    if ( ! is.null(prior_lam_fac <- rand.families[[it]]$prior_lam_fac)) {
      next_lambda_est[u.range] <- unique_lambda*prior_lam_fac
    } else next_lambda_est[u.range] <- rep(unique_lambda,length(u.range))
  }
  if (verbose["trace"]) { 
    ## this is cryptic but a better output would require a lot of reformatting as at the end of HLfit_body. 
    print(paste("unique(next_lambda_est)=",paste(signif(unique(next_lambda_est),4),collapse=" ")),quote=FALSE)
    print_corr_est <- unlist(HLfit_corrPars) 
    if ( ! is.null(print_corr_est)) print(paste("corr_est=",paste(signif(print_corr_est,4),collapse=" ")),quote=FALSE)
  }
  return(list(next_LMatrices=next_LMatrices,
              resp_lambda=resp_lambda, ## for final glm...
              HLfit_corrPars=HLfit_corrPars,
              next_lambda_est=next_lambda_est, ## heterosc
              glm_lambda=glm_lambda, ## potentially to be replaced by a list of glms later
              isRandomSlope=isRandomSlope,
              anyVanishLam=anyVanishLam
              ))
}

.process_resglm_list <- function(resglm_lambdaS, ## les 2 autres args for handling errors
                                nrand) {
  lambda_seS <- as.list(rep(NA,nrand)) # return value will be a list of length nrand
  coefficients_lambdaS <- as.list(rep(NA,nrand)) ## idem
  linkS <- list() # list of same length as resglm_lambdaS
  linkinvS <- list() # list of same length as resglm_lambdaS 
  warnmesses <- list() 
  for (glmit in seq_len(length(resglm_lambdaS))) {
    glm_lambda <- resglm_lambdaS[[glmit]]
    # next line for SEM
    coeff_lambdas <- coefficients(glm_lambda)
    coeffs_substitute <- glm_lambda$coeffs_substitute ## convoluted but avoids permutations by using the names:
    ## code for SEM 09/2015:
    ##The coefficients have different names whether a precomputed design matrix was used in ~X-1 or the data=<data.frame> syntax
    if (inherits(glm_lambda$data,"environment")) {
      locform <- glm_lambda$formula
      if (deparse(locform[[length(locform)]]) == "X - 1") {
        names(coeff_lambdas) <- colnames(glm_lambda$model$X)
      } else stop(paste("code missing for names(coeff_lambdas) with formula:",deparse(glm_lambda$formula)))
    } ## else names are already OK (and model$X does not exist)
    ## FR->FR mais je n'ai encore aucun controle sur l'identite des noms de params
    if ( ! is.null(coeffs_substitute)) coeff_lambdas <- coeffs_substitute[names(coeff_lambdas)]
    coefficients_lambdaS[[attr(glm_lambda,"whichrand")]] <- coeff_lambdas
    #
    p_lambda <- length(coeff_lambdas) ## only for the local resglm
    lambda_seS[[attr(glm_lambda,"whichrand")]] <- summary(glm_lambda,dispersion=1)$coefficients[(p_lambda+1L):(2L*p_lambda)]
    linkS[[glmit]] <- glm_lambda$family$link
    linkinvS[[glmit]] <- glm_lambda$family$linkinv
    rf <- attr(resglm_lambdaS,"rand.families")
    for (it in seq_len(length(rf))) { ## iteration over ranefs only for the local resglm
      if (tolower(rf[[it]]$family)=="gamma" && rf[[it]]$link=="identity" && coeff_lambdas[it]>0) {
        message("lambda failed to converge to a value < 1 for gamma(identity) random effects.")
        message("This suggests that the gamma(identity) model with lambda < 1 is misspecified, ")
        message("and Laplace approximations are unreliable for lambda > 1. ")
      }
    }
    warnmesses[[glmit]] <- glm_lambda$warnmess
  }
  return( list(lambda_seS =lambda_seS, #list
             coefficients_lambdaS = coefficients_lambdaS, #list
             linkS = linkS, # list of same length as resglm_lambdaS
             linkinvS = linkinvS, # list of same length as resglm_lambdaS 
             warnmesses = warnmesses  ))
}

..calcPHI <- function(processed,
                      dev.res,
                      # mv case:
                      phi.Fix,
                      #
                      lev_phi, phimodel, verbose, ## using verbose "TRACE", "trace" and "phifit"
                      # glm case:
                      method="glm", control.HLfit,
                      # hglm case:
                      residProcessed=processed$residProcessed, 
                      residModel=processed$residModel, 
                      iter, prev_PHIblob
) {
  if (phimodel == "phiHGLM") { ## random effect(s) in predictor for phi
    phifitarglist <- .update_phifitarglist(processed, 
                                           residProcessed=residProcessed,
                                           residModel=residModel,
                                           dev.res= dev.res, 
                                           lev_phi=lev_phi,
                                           iter=iter, verbose=verbose, phifit=prev_PHIblob$phifit) 
    # It's the parent 'processed' which bears the TRACE info (cf comment in residProcessed <- .preprocess(.) arguments)
    if (verbose["TRACE"]) {cat(paste("\nBegin tracing residual dispersion fit for iter=",iter,":\n"))}
    phifit <- do.call("fitme_body",phifitarglist)
    phifit$how$fnname <- "fitme_body" # for get_inits_from_fit()
    #if ( ! is.null(processed$port_env$port_fit_values)) undebug(glm.fit)
    if (verbose["TRACE"]) {cat(paste("... end tracing residual dispersion fit for iter=",iter,"."))}
    prevmsglength <- .overcat_phifit_progress(phifit, verbose, processed, iter, prev_PHIblob$prevmsglength)
    next_phi_est <- .sanitize_phi_est(phifit$fv, control.HLfit)
    return(list(next_phi_est=next_phi_est,  #low phi values are handled in calc_APHLs...
                phifit=phifit,
                prevmsglength=prevmsglength))
  } else {
    residFamily <- residModel$family
    if ( identical(deparse(residModel$formula),"~1")) { ## one case where we can easily avoid an explicit call to a glm (but one will be used to compute SEs later) 
      next_phi_est <- sum(dev.res)/sum(1-lev_phi) ## NOT in linkscale
      beta_phi <- c("(Intercept)"=residFamily$linkfun(next_phi_est)) ## linkscale value
      glm_phi <- NULL
      # to obtain the phi estimate given by summary.glm(), one must use
      #  dev.res <- wt * ((y-mu)/family$mu.eta(eta))^2 ## wt * EQL residuals
      # but the logLik differs from that given by logLik.glm(). See Details in ?HLfit
    } else { ## which means that calcPHI cannot be used for final estim phi
      glm_phi <- .calc_dispGammaGLM(formula=residModel$formula, dev.res=dev.res,
                                    data=processed$data,lev=lev_phi, family=residFamily,
                                    etastart=control.HLfit$`control.phi`$etastart, ## private, availabble, but NULL so far
                                    control=processed[["control.glm"]])
      beta_phi <- coefficients(glm_phi) ## numeric(0) if phi fixed by some offset
      next_phi_est <- fitted(glm_phi)
      if (residFamily$link!="log" && any(next_phi_est<=0)) { stop("Gamma GLM for dispersion yields negative phi estimates.") }
    }
    if (verbose["trace"] && length(beta_phi)) {cat("str(phi_est): ",str(next_phi_est))}
    next_phi_est <- .sanitize_phi_est(next_phi_est, control.HLfit)
    return(list(next_phi_est=next_phi_est,  #low phi values are handled in calc_APHLs...
                glm_phi=glm_phi,
                beta_phi=beta_phi ## used at least to initiate final GLM in "~1" case
    )) ## 
  }
}

.calcPHIs <- function(processed,
                      y, mu, wt, phi.Fix,
                      #
                      lev_phi, phimodels, verbose, ## using verbose "TRACE", "trace" and "phifit"
                      # glm case:
                      method="glm", control.HLfit,
                      # hglm case:
                      residProcessed=processed$residProcessed, 
                      residModel=processed$residModel, 
                      iter, prev_PHIblob
) {
  if ( ! is.null(vec_nobs <- processed$vec_nobs)) { # fitmv case
    cum_nobs <- c(0L,cumsum(vec_nobs))
    multiPHI <- next_phi_est <- structure(vector("list",length(vec_nobs)), names=seq_len(length(vec_nobs))) # names used in call to .Modify_list()
    for (mv_it in seq_along(vec_nobs)) if (is.null(phi.Fix[[mv_it]])) { 
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      # to obtain the phi estimate given by summary.glm(), one must use
      #  dev.res <- wt * ((y-mu)/family$mu.eta(eta))^2 ## wt * EQL residuals
      # but the logLik differs from that given by logLik.glm(). See Details in ?HLfit
      multiPHI[[mv_it]] <- ..calcPHI(processed,
                                    dev.res= processed$families[[mv_it]]$dev.resids(y[resp_range],mu[resp_range],wt[[mv_it]]), 
                                    residProcessed=processed$residProcesseds[[mv_it]],
                                    residModel=processed$residModels[[mv_it]],
                                    lev_phi=lev_phi[resp_range],
                                    phimodel=phimodels[[mv_it]], 
                                    verbose=verbose,
                                    control.HLfit=control.HLfit,
                                    #
                                    iter=iter, prev_PHIblob=prev_PHIblob$multiPHI[[mv_it]] # _FIXME_ does not contain the correct prevmsglength if several phiHGLMs...
      )
      next_phi_est[[mv_it]] <- multiPHI[[mv_it]]$next_phi_est # value of *phi* (not phi_i:= phi/prior.weights as pw are used in GLMweights, not here)
    } else next_phi_est[[mv_it]] <- phi.Fix[[mv_it]]
    return(list(multiPHI=multiPHI, next_phi_est=next_phi_est))
  } else ..calcPHI(processed,
                   dev.res= processed$family$dev.resids(y,mu,wt=wt), # eta-family !
                   #: times pw to be an estimate of same phi across level of response
                   # but not same phi as when there is no pw !
                   # double pw => double phi_est so that phi_est_i :=phi_est/pw_i is unchanged
                   lev_phi=lev_phi, phimodel=phimodels, verbose=verbose, 
                   # glm:
                   control.HLfit=control.HLfit,
                   # hglm:
                   iter=iter, prev_PHIblob=prev_PHIblob) # => list(next_phi_est=next_phi_est, glm_phi=glm_phi, beta_phi=beta_phi)
}

.addPhiGLMwarning <- function(PHIblob, models, warningList) {
  if ( ! is.null(multiPHI <- PHIblob$multiPHI)) { # fitmv case
    for (mv_it in seq_along(multiPHI)) {
      if (models[["phi"]][mv_it]=="phiHGLM") {
      } else {
        if (! is.null(locw <- multiPHI[[mv_it]]$glm_phi$warnmess)) warningList$innerPhiGLM <- locw
      }
    }
  } else if (models[["phi"]]=="phiHGLM") {
  } else {
    if (! is.null(locw <- PHIblob$glm_phi$warnmess)) warningList$innerPhiGLM <- locw
  }
  warningList
}

# always retrun a list in mv obsInfo case.
.calc_w_resid <- function(GLMweights,phi_est, obsInfo) { ## One should not correct this phi_est argument by prior.weights (checked)
  if (inherits(GLMweights,"mvlist")) { # an mvlist of lists
    mv_res <- vector("list",length(GLMweights))
    islist <- logical(length(GLMweights))
    for (mv_it in seq_along(mv_res)) {
      mv_res[[mv_it]] <- .calc_w_resid(GLMweights[[mv_it]],phi_est[[mv_it]]) # a single vector, 
                      #                or a list with the elements of the GLMweights[[mv_it]] list, plus $w_resid
      islist[mv_it] <- is.list(mv_res[[mv_it]])
    }
    if (any(islist)) { # one of the families of the mv model has zero_truncated info for Hexp
      w_resid <- WU_WT <- vector("list",length(GLMweights))
      for (mv_it in seq_along(w_resid)) {
        if (islist[mv_it]) {
          w_resid[[mv_it]] <- mv_res[[mv_it]]$w_resid
          WU_WT[[mv_it]] <- mv_res[[mv_it]][["WU_WT"]]  
        } else {
          w_resid[[mv_it]] <- mv_res[[mv_it]]
          WU_WT[[mv_it]] <- rep(1, length( mv_res[[mv_it]]))  
        }
      }
      if ( ! obsInfo) {
        w_resid <- .unlist(w_resid)
        attr(w_resid,"unique") <- attr(w_resid,"is_unit") <- FALSE 
      }
      return(list(w_resid=w_resid, WU_WT=.unlist(WU_WT), mvlist=mv_res)) # mvlist being a list with elements of two possible types as described above
      # each element of 'mvlist' has all the info in the expected info for a single response
    } else {
      if (obsInfo) {
        return(list(w_resid="non", mvlist=mv_res))
      } else {
        w_resid <- .unlist(mv_res) # vector
        attr(w_resid,"unique") <- attr(w_resid,"is_unit") <- FALSE 
        # .calc_H_global_scale() assumes the "unique" attribute is here, for mv or not.
        # is_unit expected to be true only for augZXy hence not for mv 
        return(w_resid) # mvlist being a list with elements of tow possible types as described above
      }
    }
    # in mv model w_resid are not 'unique' although blocks may be so, for gaussian submodels.
  }
  
  if (is.character(GLMweights)) return(GLMweights) # can still be a list so don"t test is.numeric()
  
  phi_est[phi_est<1e-12] <- 1e-11 ## 2014/09/04 local correction, cf comment in calc_APHLS...
  if (has_ZT_Hexp_info <- inherits(GLMweights,"list")) { ## has zero_truncated info for Hexp
    res <- GLMweights ## includes WU_WT, d3logMthdth3 .... so that d3logMthdth3 ... becomes an element of w.resid...
    GLMweights <- GLMweights$truncGLMweights
  }
  #if (any(is.nan(GLMweights/phi_est))) stop("any(is.nan(GLMweights/phi_est))")
  w.resid <- as.vector(GLMweights/phi_est)
  attr(w.resid,"unique") <- (attr(GLMweights,"unique") && length(phi_est)==1L)
  attr(w.resid,"is_unit") <- FALSE
  if (has_ZT_Hexp_info) {
    res[["w_resid"]] <- w.resid 
    return(res) # a list with the elements of the 'GLMweights' list, plus $w_resid
  } else {
    return(w.resid) # vector
  }
}

.calc_H_global_scale <- function(w.resid) { # where w.resid depends on phi_est and on prior weights
  # small adjustment added 2019/11/11 resolve some numerical pb in extreme binary fits (mu \approx y, w.resid->0)
  if (is.list(w.resid)) w.resid <- w.resid$w_resid
  # is_unique <- attr(w.resid,"unique")
  # if (is.null(is_unique)) {
  #   warning('Presumably inefficient code: missing attr(w.resid,"unique")')
  #   u_w_resid <- unique(w.resid)
  #   is_unique <- (length(u_w_resid)==1L)
  # } else if (is_unique) u_w_resid <- w.resid[1L]
  H_scale_regul <- .spaMM.data$options$H_scale_regul
  if (attr(w.resid,"unique")) { # __FIXME__ more efficient code could be conceived for mv models
    u_w_resid <- w.resid[1L]
    H_global_scale <- 1/u_w_resid # so that weight_X is rep(1,...) 
    # unbiased in w.resid=1 (if lop1p used,  difference between lop1p(x) and log(1+x) could create a 'bias' => don't use log1p):
    if (H_global_scale>1/H_scale_regul) H_global_scale <- exp(log(H_scale_regul+1) - log(H_scale_regul+u_w_resid)) # test '871.8141' affected if this is not conditional
  } else {
    if  ( ! is.null(attr(w.resid,"signs"))) w.resid <- abs(w.resid) # quite local
    H_global_scale <- exp(- mean(log(w.resid)))
    if (H_global_scale>1/H_scale_regul) H_global_scale <- exp(log(H_scale_regul+1) - mean(log(H_scale_regul+w.resid))) # test '871.8141' affected if this is not conditional
  }
  return(H_global_scale) 
}

.D_dlWdmu_Deta <- function(family,mu) { # mu is muFREQS;      detailed computations in Hessian_weights.nb
  # .D_dlWdmu_Deta() has no fallback so it is the limiting factor in implementing obsInfo.
  switch(family$family,
         # "Gamma" log not needed but others missing ?
         "negbin2" = family$d_dlWdmu_detafun(mu), # using the shape parameter
         "poisson" = family$d_dlWdmu_detafun(mu), 
         "binomial" = {
           switch(family$link,
                  "logit" = (1 - 2*mu + 2*mu^2)/((-1 + mu)*mu),  # Hypothetical since not needed
                  "cloglog" = { 
                    log1mu <- log(1 - mu)
                    (-2 *mu^2 - 2 *mu^2 *log1mu + (1 - 2*mu)* log1mu^2)/((-1 + mu)* mu^2 * log1mu)
                  }, 
                  "probit" = {
                    eta <- qnorm(mu)
                    denom <- ((mu*(1-mu))^2)
                    dnorm_eta <- dnorm(eta)
                    dnorm_eta * (1-2*mu+2*mu^2+ 2*(2*mu^3-mu^4-mu^2-denom*eta^2)/dnorm_eta^2 ) / denom
                  },
                  "cauchit" = {
                    -(-1 + 2*mu - 16*mu^3 *pi^2 + 8*mu^4 *pi^2 + 
                         mu^2 *(-2 + 8*pi^2) + (1 - 2*mu + 2*mu^2)*cos(2*mu*pi))/(
                           2*(-1 + mu)^2 *mu^2 *pi)
                  },
                  stop("d(dlWdmu)deta cannot yet be computed for this link.") 
           )
         }, 
         "gaussian" = {
           switch(family$link,
                  "log" = -2/mu
           )
         },
         stop("d(dlWdmu)deta cannot yet be computed for this family.")
  )
}

.DdetadmuDeta <- function(link, mu) { # in Hessian_weights.nb; mu is muFREQS
  switch( # d (deta/dmu) / deta # Cf Hessian_weights.nb
    link,
    "log" =  -1/mu,
    "inverse" = -2/mu, # in the Gamma(inverse) sense
    "identity" = 0,
    "sqrt" = -1/(2*mu),
    "cloglog" = {(1 + 1/log(1 - mu))/(1 - mu)}, # Hypothetical, not needed in binomial(logit) fits
    "logit" = {(1 - 2 *mu)/((-1 + mu)* mu)},
    "probit" = {
      eta <- qnorm(mu)
      eta/dnorm(eta)
    },
    "cauchit" = -2*pi/tan(pi*mu),
    # loglambda not necess bc appears only as canonical link
    stop("This cannot yet be computed for this family link.") 
  )
}

# .D2detadmuDeta2 is defined further below



# .dlW_Hexp__dmu used in obsInfo by NON-LLgeneric algo  to compute Hobs_w.resid.
# But note that the therein-called function dlW_Hexp__detafun() is also used for expInfo 
#
# Similar to  ..calc_dlW_deta(obsInfo=FALSE), but not needing to handle mv, and as fn of mu, not eta.
# Indeed calls ..calc_dlW_deta() when there is no ad-hoc code.
#
# There are calculations of special cases in Hessian_weights.nb
.dlW_Hexp__dmu <- function(family, # always needed
                     ### (possibly) needed for predefined cases: (mu always, and one of dmudeta and dmu.deta)
                     mu, # mu or muFREQS
                     # non binomial case:
                     dmudeta, # = muetablob$dmudeta, # could use muetablob below
                     # binomial case
                     dmu.deta=family$mu.eta(family$linkfun(mu)), # for binomial, differ from dumdeta=dmu.deta*BinomialDen
                     ### (possibly) needed for generic fallback code.
                     muetablob, # always needed for generic fallback code. 
                     BinomialDen # May be needed for generic fallback code, which should use it only in binomial cases.
                     ) {
  dlW_dmu <- switch( # dlW_deta * [deta.dmu written as 1/dmu.deta|_eta=eta(mu)
    family$family,
    "Gamma" = {
      switch(family$link,
             "log" = 0, 
             ..calc_dlW_deta(muetablob, family, calcCoef1=FALSE, w.resid=NULL, 
                             processed=NULL, Hratio_factors=NULL,
                             obsInfo=FALSE)$dlW_deta/dmudeta      
      )
    } 
    ,
    "negbin2" = family$dlW_Hexp__detafun(mu)/dmudeta, # family member fn for access to the shape parameter... 
    "poisson" = family$dlW_Hexp__detafun(mu)/dmudeta, # Poisson family member. Maybe I could move the code here 
    #      (design decision to avoid replacing poisson with Poisson unless for truncated which oes not use this code)
    "binomial" = {
      switch(family$link,
      "logit" = (1-2*mu)/dmu.deta,  # or (1-2*mu)/(mu*(1-mu)) # Hypothetical, not needed in binomial(logit) fits
      "cloglog" = (1/mu + 2/log(1 - mu))/(-1 + mu), # or (2 + log(1 - mu)/mu)/dmu.deta
      "cauchit" = (1 - 2*mu + 4*(-1 + mu)*mu*pi/tan(mu*pi))/((-1 + mu)*mu),
      # "probit" in notebook, is ugly
      ..calc_dlW_deta(muetablob, family, calcCoef1=FALSE, w.resid=NULL, 
                      BinomialDen=BinomialDen, processed=NULL, Hratio_factors=NULL,
                      obsInfo=FALSE)$dlW_deta/dmu.deta      # note difference from .../dmudeta
      )
    },
    "gaussian" = {
      switch(family$link,
             "log" = 2/mu
             )
    },
    ..calc_dlW_deta(muetablob, family, calcCoef1=FALSE, w.resid=NULL, 
                    #BinomialDen=BinomialDen,
                    processed=NULL, Hratio_factors=NULL,
                    obsInfo=FALSE)$dlW_deta/dmudeta 
  )
  dlW_dmu
}


..calc_H_w.resid <- function(w.resid, Hratio_factors, partials=Hratio_factors$partials, yMmu=Hratio_factors$yMmu, BinomialDen) {
  if (is.list(w.resid)) w.resid <- w.resid$w_resid # has a "unique" attribute
  w.resid* (1 - partials * yMmu/BinomialDen) 
}

.calc_H_w.resid <- function(w.resid, 
                            muetablob, processed, # as providers of next variables... 
                                                  # but in fact muetablob distinctly needed for obsInfo
                                                  # would be a problem for a postfit call of .calc_H_w.resid().
                            mu=muetablob$mu, y=processed$y,
                            obsInfo=processed$how$obsInfo,
                            family=processed$family,
                            families=processed$families,
                            BinomialDen=processed$BinomialDen
                            ) {
  if (obsInfo &&  ! is.null(cum_nobs <- attr(families,"cum_nobs"))) {
    Hobs_w.resid <- numeric(length(muetablob$mu))
    for (mv_it in seq_along(families)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      Hobs_w.resid_it <- .calc_H_w.resid(w.resid="ignore", 
                                         muetablob=muetablob$mv[[mv_it]],
                                         y=y[resp_range],
                                         obsInfo=obsInfo,
                                         family=families[[mv_it]],
                                         BinomialDen=BinomialDen[resp_range],
                                         families=NULL
      )
      Hobs_w.resid[resp_range] <- Hobs_w.resid_it
    }
    if (any(Hobs_w.resid<0)) { # absence of attribute meaning that all signs are >0;
      attr(Hobs_w.resid,"signs") <- as.vector(sign(Hobs_w.resid)) # as.vector drops the attribute own attributes, as sign() keeps them (...)
    }
    attr(Hobs_w.resid,"unique") <- FALSE # presumably
    attr(Hobs_w.resid,"is_unit") <- FALSE # presumably; used by spprec
    return(Hobs_w.resid)
  }
  # 
  if (obsInfo) {
    Hobs_w.resid <-  muetablob$Md2logcLdeta2 # present for LLF... but no longer only for them: truncated GLM families too.
    if (is.null(Hobs_w.resid)) { # obsInfo but not LLF
      # we seek MINUS the derivative, derivative in with the first term is -1 from D(y MINUS mu)...
      # ergo     -  ( -1 + (y-mu) ...)
      if (family$flags$canonicalLink) { # This is reached in the mv case mixing canonical and non-canonical links
        if (is.list(w.resid)) w.resid <- w.resid$w_resid # has a "unique" attribute
        Hobs_w.resid <- w.resid
        rep0L <- rep(0L, length(w.resid))
        attr(Hobs_w.resid,"Hratio_factors") <- list("yMmu"=rep0L, # not true but should not matter
                                                    "partials"= rep0L) # true...
      } else {
        yMmu <- drop(y-mu)
        if (is.list(w.resid)) yMmu <- yMmu -w.resid$dlogMthdth # truncated family: Tpoisson or Tnegbin 
                                #but neither pure LLF nor in the "hybrid" family negbin2 (which only fits by Hobs)
                                # so this code not useful for comparison of computations in negbin2 case.
                                # we still want to allow Hexp truncated fits so we cannot always use LLF-like algos.
        if (family$family=="binomial") {
          muFREQS <- mu/BinomialDen
          partials <- .DdetadmuDeta(family$link, mu=muFREQS) + 
            .dlW_Hexp__dmu(family=family, mu=muFREQS, 
                     ## for generic fallback code:
                     muetablob=muetablob, BinomialDen=BinomialDen)
          ## Hobs_w.resid <- w.resid* (1 - drop(y-mu)*(DdetadmuDeta + dlW_dmu)) :
          ## When there is a BinomialDen, this is
          # GLMweights * (BinomialDem- (y-muCOUNT)*partials)
          # =  w.resid* (1 - (y-muCOUNT)*partials/BinomialDen)  bc w.resid=BinomialDen * GLMweights. Cf gradientsHobs_cloglog.nb
          Hobs_w.resid <- ..calc_H_w.resid(w.resid, partials=partials, yMmu=yMmu,BinomialDen=BinomialDen)
        } else {
          partials <- .DdetadmuDeta(family$link, mu=mu) + 
            .dlW_Hexp__dmu(family=family, mu=mu, dmudeta = muetablob$dmudeta,
                     ## for generic fallback code:
                     # BinomialDen=BinomialDen, # should not be needed
                     muetablob=muetablob)
          if (is.list(w.resid)) partials <- partials * w.resid$WU_WT # again Tpoisson or Tnegbin (only: see further comments on negbin2 above)
          Hobs_w.resid <- ..calc_H_w.resid(w.resid, partials=partials, yMmu=yMmu,BinomialDen=1)
        }
        # if (any(Hobs_w.resid<0)) browser()
        attr(Hobs_w.resid,"Hratio_factors") <- list("yMmu"=yMmu, "partials"= partials)
      }
    } else { # generic LLF procedure
      if (any(Hobs_w.resid<0)) { # absence of attribute meaning that all signs are >0;
        attr(Hobs_w.resid,"signs") <- as.vector(sign(Hobs_w.resid)) # as.vector drops the attribute own attributes, as sign() keeps them (...)
      }
    }
    attr(Hobs_w.resid,"unique") <- FALSE # presumably
    attr(Hobs_w.resid,"is_unit") <- FALSE # presumably; for spprec
    return(Hobs_w.resid)
  } else {
    if (is.list(w.resid)) w.resid <- w.resid$w_resid # has a "unique" attribute
    w.resid
  }
}

.calc_weight_X <- function(Hobs_w.resid, H_global_scale, obsInfo) {
  scaled_H_w.resid  <- as.vector(H_global_scale*Hobs_w.resid) # signed ! ; as.vector todrop attributes (but attributes of "H_w.resid" attribute below must remain) 
  if ( ! is.null(attr(Hobs_w.resid,"signs"))) {
    weight_X <- sqrt(abs(scaled_H_w.resid)) # keeps the "signs" and other attributes
  } else weight_X <- sqrt(scaled_H_w.resid)
  attr(weight_X,"H_w.resid") <- Hobs_w.resid # # signed & UNscaled ! non always required
  weight_X
}

## spaMM_Gamma() fixes Gamma()$dev.resids(1e10+2,1e10,1) is < 0
# dev.resids() must be >0 for computation deviance_residual in fitting Gamma GLMM, and also for $aic() computation.
spaMM_Gamma <- local({
  link_warned <- FALSE
  function (link = "inverse") {
  mc <- match.call()
  if (is.null(mc$link) && ! link_warned) {
    message("Gamma family's default link is 'inverse', not 'log'")
    link_warned <- TRUE
  }
  linktemp <- substitute(link) ## does not evaluate
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp) ## converts to char the unevaluated expression
  okLinks <- c("inverse", "log", "identity")
  if (linktemp %in% okLinks) 
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link) ## evals expression converted to char (with  $name in particular); but returns link=linktemp, not link=stats$name
    # problem is that the families fns return link=linktemp, which seems weird: better is   
    linktemp <- stats$name ## line not in Gamma() [and different in binomial()], which absence prevents programming with link argument...  
  } else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name)) 
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for gamma family; available links are %s", 
                    linktemp, paste(sQuote(okLinks), collapse = ", ")), 
           domain = NA)
    }
  }
  variance <- function(mu) mu^2 ## problem for extreme mu values
  validmu <- function(mu) all(mu > 0)
  dev.resids <- function(y, mu, wt) { # mu<0 -> Inf; mu=0 -> NaN; some huge mu's that could yield dev_res<0 => fixed
    if (any(mu < 0)) return(Inf) ## 2015/04/27; shortcut; without it next line warns about NaNs produced and NaN's are returned
    dev_res <- -2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
    dev_res[dev_res<.Machine$double.eps] <- .Machine$double.eps ##  ## FR: added this (*would* ignore NaN's)
    ## Once in Gamma GLM, y=25.75563, y-mu=-5.996906e-08, yielded negative dev.resid
    return(dev_res)
  }
  aic <- function(y, n, mu, wt, dev) {
    n <- sum(wt)
    disp <- dev/n
    -2 * sum(dgamma(y, 1/disp, scale = mu * disp, log = TRUE) * 
               wt) + 2
  }
  initialize <- expression({
    if (any(y <= 0)) stop("non-positive values not allowed for the gamma family") 
    n <- rep.int(1, nobs)
    mustart <- y
  })
  simfun <- function(object, nsim) {
    wts <- object$prior.weights
    if (any(wts != 1)) 
      message("using weights as shape parameters")
    ftd <- fitted(object)
    shape <- MASS::gamma.shape(object)$alpha * wts ## explicit MASS:: is in source of stats::Gamma() !
    resu <- rgamma(nsim * length(ftd), shape = shape, rate = shape/ftd)
    if (nsim>1L) resu <- matrix(resu,ncol=nsim)
    resu
  }
  # linkinv <- function (eta) pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax) 
  # : permet des plantages severes dans glm.fit ensuite (R CMD check en detecte) 
  ## all closures defined here have parent.env the environment(spaMM_Gamma) ie <environment: namespace:spaMM>
  ## changes the parent.env of all these functions (aic, dev.resids, simfun, validmu, variance): 
  # as.list(environment(aic)) ## this has an unexplained effet on saveSize!
  parent.env(environment(aic)) <- environment(stats::Gamma) ## parent = <environment: namespace:stats>
  ## That _does_ reduce the size of the fitted objects using spaMM_Gamma (eg in a phi.object)
  ## That does not eliminate an hidden environment shared among member functions 
  #    _after_ compiling the package:  
  ##  spaMM:::.saveSize(attr(attr(spaMM_Gamma()$aic,"srcref"),"srcfile")) grows
  ## compared to the non-compiled version
  ## It _is_ an environment : ls(attr(attr(spaMM_Gamma()$aic,"srcref"),"srcfile")) lists it.  
  ## But its size is not explained by its contents...
  structure(list(family =  structure("Gamma",patch="spaMM_Gamma"), 
                 link = linktemp, linkfun = stats$linkfun, 
                 linkinv = stats$linkinv, 
                 variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun), 
            class = "family")
}
  })


.get_clik_fn <- function(family) {
  # return value of each fn must be a vector if y is a vector
  switch(family$family,
         gaussian = function(theta,y,nu) {nu*(theta*y-(theta^2)/2)- ((y^2)*nu+log(2*pi/nu))/2}, 
         poisson = function(theta,y,nu) { ## matches different families : stats::poisson, Poisson, Tpoisson...
           ## with theta = log(mu)
           # res <- nu*(theta*y-attr(theta,"mu"))   -  lfactorial(y)
           # res[theta== -Inf & y==0] <- 1
           # res
           ## uses C code:
           -family$aic(y=y, mu=exp(theta), wt=1)/2 ## handles truncation from given untruncated mu
         },
         binomial = function(theta,freqs,sizes,nu) {
           #nu*sizes*(freqs*theta-log(1+exp(theta))) +lchoose(sizes,round(sizes*freqs)) :
           ## uses stats:: C code:
           muFREQS <- plogis(theta)
           nu * dbinom(round(sizes * freqs), sizes, prob=muFREQS, log = TRUE)
          },
         # gamma = function(theta,y,nu) {nu*(y*theta+log(-theta))+nu*(log(nu*y))-lgamma(nu)-log(y)} ## mean mu=-1/th, **** var = mu^2 / vu ****
         # same by using ad hoc C library...
         Gamma = function(theta,y,nu) {
           dgamma(y, shape=nu , scale = attr(theta,"mu")/nu, log = TRUE) 
         },
         COMPoisson = function(theta, y, nu, muetaenv=NULL) {
           COMP_nu <- environment(family$aic)$nu # using the (variable nu) from the family captured in this fn's definition envir 
           .CMP_clik_fn(theta, y, nu, COMP_nu, muetaenv=muetaenv)
         },
         negbin = function(theta,y,nu) { ## theta is the canonical param, -log(1+shape/mu)
           if (is.null(mu <- attr(theta,"mu"))) {
             NB_shape <- environment(family$aic)$shape
             mu <- NB_shape/(exp(-theta)-1) # note that NB_shape/(exp(log(1 + NB_shape/x)) - 1) numerically fails to be x for large x
           }
           -family$aic(y=y, mu=mu, wt=1)/2 ## handles truncation from given untruncated mu
         },
         negbin2 = function(theta,y,nu) { ## theta is the canonical param, -log(1+shape/mu)
           if (is.null(mu <- attr(theta,"mu"))) {
             NB_shape <- environment(family$aic)$shape
             mu <- NB_shape/(exp(-theta)-1) # note that NB_shape/(exp(log(1 + NB_shape/x)) - 1) numerically fails to be x for large x
           }
           -family$aic(y=y, mu=mu, wt=1)/2 ## handles truncation from given untruncated mu
         },
         negbin1 = function(theta,y,nu) { ## theta is a pseudo-canonical parameter... defined as mu (pseudo canonical link def as identity)
           -family$aic(y=y, mu=theta, wt=1)/2 ## Should handle truncation from given untruncated mu
         },
         beta_resp= function(theta,y,pw) # note how prior weights are passed by .calc_clik -> clik_fn(theta, y, pw=eval(prior.weights))
           -family$aic(y=y, mu=theta, wt=pw)/2,
         stop("code missing for this family")
  )
}


`.theta.mu.canonical` <- function(mu,family) { # essential for all families; and LLF-GLM ones must has ad-hoc case in the switch
  ## the (fixed) canonical link between theta and mu, not the family link between eta and mu 
  if (inherits(family,"family")) {
    famfam <- family$family
  } else famfam <- family
  switch(famfam,
         gaussian = mu ,
         poisson = { # structure(log(mu),mu=mu)
           th <- log(mu)
           # attr(th,"mu") <- mu # mu # presumably not used
           th
         },
         binomial = .logit(mu), #                
         # or slow access to C_logit_link:         
         #   make.link("logit")$linkfun(mu), # wraps C_logit_link, but make.link calls structure()... 
         ## if any numerical issue, use 
         #                 { 
         #                    theta <- .logit(mu)
         #                    theta[theta>27.6310] <- 27.6310 ## mu>1-1e-12
         #                    theta[theta < -27.6310] <- -27.6310 
         #                    theta
         #                 },
         Gamma = { # structure(-1/mu,mu=mu), ## "-": McC&N p. 290
           th <- -1/mu
           attr(th,"mu") <- mu # for clik_fn, as attr(theta,"mu")
           th
         },
         COMPoisson = { 
           if (is.null(lambda <- attr(mu,"lambda"))) {
             lambda <-  family$mu2lambda(mu) 
           }  
           th <- log(lambda) 
           attr(th,"mu") <- mu # structure(log(lambda),mu=mu) with mu used by .CMP_loglambda_linkinv
           th
         },
         negbin = { #structure(-log(1+(environment(family$aic)$shape)/mu),mu=mu)
           th <- -log(1+(environment(family$aic)$shape)/mu)
           attr(th,"mu") <- mu ## keep mu b/c useful for clik_fn, as attr(theta,"mu")
           th
         },
         negbin2 = { #structure(-log(1+(environment(family$aic)$shape)/mu),mu=mu)
           th <- -log(1+(environment(family$aic)$shape)/mu)
           attr(th,"mu") <- mu ## keep mu b/c useful for clik_fn, as attr(theta,"mu")
           th
         },
         # For 'LLM' families outside the exponential class, let the 'canonical' link be the identity link => theta=mu  
         if ( ! family$flags$exp ) {
           mu
         } else stop("code missing for this family") 
  )
} ## returns values for given mu

# For families which do not have specific Md*logcLdeta* functions for obsInfo code (cf call in ..calc_dlW_deta() excluding case where these fns are available)
# In particular for the Hexp case of standard GLMs
.thetaMuDerivs <- function(mu, family, 
                           muFREQS, # for Binomial 
                           muetablob # for COMPoisson
                           ) { ## used for non-canonical links
  familyfam <- family$family
  ## these definitions depend only on the canonical link
  Dtheta.Dmu <- switch(familyfam,
                       gaussian = rep(1,length(mu)) ,
                       poisson = 1/mu ,
                       binomial = 1/(muFREQS*(1-muFREQS)),
                       Gamma = 1/mu^2,
                       negbin = family$Dtheta.Dmu(mu), # 1/(mu*(1+mu/NB_shape)),
                       COMPoisson = muetablob$thetaMuDerivs_2$Dtheta.Dmu,  
                       # nothing for non-GLM families
                       stop("code missing") # Dtheta.Dmu deriv(body(<family(<canonical link>)>$linkfun),"mu") but see binomial() or Gamma() for various issues
  ) ## values for given mu
  # if (familyfam=="binomial") Dtheta.Dmu <- Dtheta.Dmu/BinomialDen removed 2022/05/01. No bug anytime, but change of design:
  # up to that time the calling code used the present results in combination with muetablob$dmudeta which has a BinomialDen factor. 
  # The calling code has now been modified. 
  D2theta.Dmu2 <- switch(familyfam,
                         gaussian = rep(0,length(mu)) ,
                         poisson = -1/mu^2 ,
                         binomial = -(1-2*muFREQS)/(muFREQS*(1-muFREQS))^2,
                         Gamma = -2/mu^3,
                         negbin = family$D2theta.Dmu2(mu), # -(1+2*mu/NB_shape)/(mu*(1+mu/NB_shape))^2,
                         COMPoisson = muetablob$thetaMuDerivs_2$D2theta.Dmu2,  
                         stop("code missing") 
  ) ## values for given mu
  # if (familyfam=="binomial") D2theta.Dmu2 <- D2theta.Dmu2/(BinomialDen^2) removed 2022/05/01. Cf comment above.
  return(list(Dtheta.Dmu=Dtheta.Dmu,D2theta.Dmu2=D2theta.Dmu2))
}

.sanitize_eta_log_link <- function(eta, max, min=-max,y, nu=1, warn_neg_y=TRUE) { # cf reasons below to always provide y 
  ### Before implementation of y argument:
  ## 'max' used to sanitize eta should be higher than the max used to sanitize the response in .calc_dispGammaGLM()
  ## Otherwise we would never have sanitized eta = sanitized y in cases where this would be the correct solution.
  ## More generally the LevM algo may then be stuck.
  ### BUT: additional concept in v.3.0.49; extend bound according to y (for both sanitized eta & sanitized y)
  #  test code with huge .calc_dispGammaGLM() response:
  # {
  #   spaMM.options(example_maxtime=60)
  #   Infusion::Infusion.options(spaMM.options("example_maxtime"))
  #   options(warn=2)
  #   options(error=recover)
  #   #Infusion.options(example_maxtime=20)  ## allows basic example
  #   #Infusion.options(example_maxtime=120)  ## allows refine()
  #   example("example_raw",package="Infusion",ask=FALSE) 
  # }
  # max=40 means mu < exp(40)=2.353853e+17
  if (any(is.infinite(eta))) {
    isinf <- which(is.infinite(eta))
    eta[isinf] <- sign(eta[isinf]) * .Machine$double.xmax  # following code further sanitizes eta if within +/- double.xmax
  }
  tmp1 <-  6.5663667753507884 # log(1+log(.Machine$double.xmax))
  tmpmax <- (max - tmp1)*nu # when nu<1, this shrinks both the maximum and the range over which a correction is applied
  tmpmin <- min+tmp1 # nu appears to have little effect on mu when eta is small
  if ( ! is.null(y)) {
    if (warn_neg_y) {
      logy <- log(y) 
    } else logy <- suppressWarnings(log(y)) # for gaussian(log) y can be negative.
    tmpmax <- pmax(logy, tmpmax) ## if log(y)>tmpmax corrected eta is log(y)+(< log(1+double.xmax)=tmp1)*nu i.e. corrected eta is < log(y)+tmp1*nu
    #   otherwise  corrected eta is < (max-tmp1)*nu+(< log(1+double.xmax)=tmp1)*nu    <   max*nu  
    #  and "-tmpmax" ensures continuity in eta=tmpmax) as log(1+log(1+...))=0 there
    high_eta <- ( (!is.nan(tmpmax)) & eta>tmpmax) # added the is.nan(tmpmax) condition to allow some negative y when fitting gaussian(log)
    eta[high_eta] <- tmpmax[high_eta] + log(1+log(1+eta[high_eta]-tmpmax[high_eta]))*nu
    ## smooth correction proved useful in spaMM_glm.fit (at least), presumably to avoid a plateau of objective fn
    #
    tmpmin <- pmin(tmpmin,logy)
    low_eta <- ( (!is.nan(tmpmax)) & eta< tmpmin)
    eta[low_eta] <- tmpmin[low_eta] - log(1+log(1-eta[low_eta]+tmpmin[low_eta]))
  }  else {
    eta[eta>tmpmax] <- tmpmax + log(1+log(1+eta[eta>tmpmax]-tmpmax))*nu
    eta[eta< tmpmin] <- tmpmin - log(1+log(1-eta[eta< tmpmin]+tmpmin))## symmetric correction for negative eta
  }
  return(eta)
} 

# using y is important in cases where a comparison btwn eta (or mu) and y is important for convergence of an algo
# In that case we sanitize only to the extant that individual values of y allow.
# If this sanitization is not enough, then it is y which must be sanitized.
.sanitize_eta <- function(eta, y=NULL,family,max=.spaMM.data$options$sanitize_eta["otherlog"],
                          bin_mu_tol=.spaMM.data$options$bin_mu_tol) {
  if (family$link =="log") {
    if (family$family=="gaussian") {
      eta <- .sanitize_eta_log_link(eta, max=.spaMM.data$options$sanitize_eta["gauslog"],  y=y, warn_neg_y = FALSE) 
    } else if (family$family=="COMPoisson") {
      eta <- .sanitize_eta_log_link(eta, max=.spaMM.data$options$sanitize_eta["COMPlog"],  y=y) 
    } else eta <- .sanitize_eta_log_link(eta, max=max,  y=y)
  } else if (family$family == "COMPoisson" && family$link =="loglambda") {
    # this should be consistent with poisson(log) when nu=1, except that we may wish to avoid the computational burden of large eta values
    COMP_nu <- environment(family$aic)$nu 
    eta <- .sanitize_eta_log_link(eta, max=max, y=y, nu=COMP_nu) ## will use log(mu) ~ eta/nu for large eta and small nu
  } else if (family$link=="inverse") {
    if (family$family=="Gamma") {
      etamin <- sqrt(.Machine$double.eps)
      neg_eta <- (eta<etamin)
      eta[neg_eta] <- etamin ## both eta and mu must be >0
      attr(eta,"any_neg_eta") <- any(neg_eta) ## actually 'not strictly positive'
    } else if (family$family=="gaussian") {
      absetamax <- 1/.Machine$double.eps # GLMweights ~ -1/eta^3...
      huge_eta <- (abs(eta)>absetamax)
      eta[huge_eta] <- sign(eta[huge_eta]) * absetamax # a smooth correction as for log link would be better...
      absetamin <- sqrt(.Machine$double.eps)
      tiny_eta <- (abs(eta)<absetamin)
      eta[tiny_eta] <- sign(eta[tiny_eta]) * absetamin ## both eta and mu must be !=0
      attr(eta,"any_tiny_eta") <- any(tiny_eta) ## actually 'not strictly positive'
    }
  } else if (family$link=="sqrt") {
    etamax <- sqrt(1/.Machine$double.eps)
    big_eta <- (eta>etamax)
    eta[big_eta] <- etamax 
    etamin <- 1/etamax
    neg_eta <- (eta<etamin)
    eta[neg_eta] <- etamin 
  } else if (family$family == "binomial") { ## y ignored here: important as there would be ambiguity whether y is counts or frequencies in calling function.
    #for binomial cases stats::binomial(...)$linkinv corrects eta for all links so that mu is always .Machine$double.eps from 0 or 1
    # This may not be enough... or inconsistent with corrections used elsewhere, so we overcome the stats:: correction
    # by correcting eta differently here and by computing mu <- .binomial_raw_linkinv(eta,link=family$link) (in .muetafn and possibly elsewhere)
    eta <- .binomial_corr_eta(eta,link=family$link, tol=bin_mu_tol)
  }
  return(eta)
}


.DlogM_Tnegbin_mpfr <- function(mu, shape, theta=NULL) {
  if (is.null(theta)) theta <- - log1p(shape/mu)
  umeth <- -expm1(theta)
  eth <- 1-umeth
  gmu <- .do_call_wrap("mpfr",arglist=list(x=mu, precBits=128),pack = "Rmpfr")
  gsh <- .do_call_wrap("mpfr",arglist=list(x=shape, precBits=128),pack = "Rmpfr")
  gth <- - log1p(gsh/gmu)
  umeth <- -expm1(gth)
  eth <- 1-umeth
  g2 <- -(umeth^gsh-gsh*umeth+gsh-1)*((eth*umeth^(gsh-2)*gsh)/(umeth^gsh-1)^2)
  fac <- (1+eth)*(-1+umeth^gsh)^2 + 3*gsh*eth*(-1+umeth^gsh) + (gsh*eth)^2*(1+umeth^gsh)
  g3 <- -(fac)*((eth*umeth^(gsh-3)*gsh)/(umeth^gsh-1)^3)
  return(list(d2logMthdth2=as.numeric(g2), d3logMthdth3=as.numeric(g3)))
}

.muetafn_truncated_GLM <- local({
  Rmpfr_warned <- FALSE
  mumax <- log(.Machine$double.xmax*1e-8)/2 # so that expmu ^2 does not overflow 
  function(family, mu, GLMweights, Vmu) {
    if (family$family=="poisson") { 
      ## D[Log[1 - E^-E^theta], {theta, 2}] /. {theta -> Log[mu]} // Simplify
      mu <- pmin(mu,mumax)
      expmu <- exp(mu)
      p0 <- 1/expmu
      dlogMthdth <- -mu/(1-expmu) ## useful to correct rhs
      d2logMthdth2 <- -(1 + expmu * (mu-1))* mu/((expmu-1)^2)  # (d2logMthdth2 -> 0 for Infinite mu) 
      exp2mu <- expmu^2
      d3logMthdth3 <- mu * (1 + 3 * mu *(expmu - exp2mu) + (mu^2) *(expmu + exp2mu) + exp2mu - 2 * expmu )/(expmu-1)^3
    } else if (family$family=="negbin2") { ## computations in TruncatedNegBin.nb
      ## D[Log[1 - E^-E^theta], {theta, 2}] /. {theta -> Log[mu...]} // Simplify
      shape <- .get_family_par(family)
      p0 <- .negbin2_p0(mu,shape) ## (1-p)^r
      dlogMthdth <- (mu * p0)/(1-p0) ## family-specific relation btwn correction term M(th) for truncation and p0
      # d2logMthdth2 <- -((mu * p0 *(shape *(-1 + p0) + mu *(-1 + shape + p0)))/(shape *(-1 + p0)^2))
      # d3logMthdth3 <- -((mu * p0 * (shape^2 *(-1 + p0)^2 + 3 * mu * shape *(-1 + p0) * (-1 + shape + p0) + 
      #                                 mu^2 *(3 *shape *(-1 + p0) + 2 *(-1 + p0)^2 + shape^2 *(1 + p0))))/(
      #                                   shape^2 *(-1 + p0)^3))
      theta <- - log1p(shape/mu)
      umeth <- -expm1(theta)
      eth <- 1-umeth
      # num is -((mu * p0 *(shape *(-1 + p0) + mu *(-1 + shape + p0)))/(shape) and denom is (-1 + p0)^2
      d2logMthdth2 <- -(umeth^shape-shape*umeth+shape-1)*((eth*umeth^(shape-2)*shape)/(umeth^shape-1)^2)
      lowmu <- ((mu^2 * (1+shape)/(shape))<1e-11) # derived from second order term in Normal[Series[(shape/(mu + shape))^shape, {mu, 0, 2}]]
      if (any(lowmu)) {
        lmu <- mu[lowmu]
        # see section on approximations in the notebook
        d2logMthdth2[lowmu] <- -(lmu^3 *(lmu^2 + (2 + (-2 + lmu)* lmu)* shape)* (lmu + (-1 + lmu) *shape + shape^2))/(4 *shape^3 *umeth[lowmu]^shape-1)^2
      }
      #
      fac <- (1+eth)*(-1+umeth^shape)^2 + 3*shape*eth*(-1+umeth^shape) + (shape*eth)^2*(1+umeth^shape)
      d3logMthdth3 <- -(fac)*((eth*umeth^(shape-3)*shape)/(umeth^shape-1)^3)
      #
    } 
    truncGLMweights <- GLMweights*(1+d2logMthdth2/Vmu) 
    if (any(wrong <- (truncGLMweights<=0))) {
      if (family$family=="negbin2") { 
        has_Rmpfr <- suppressWarnings(do.call("require",list(package="Rmpfr", quietly = TRUE))) # given it's not in DESCRIPTION
        if (has_Rmpfr) {
          dnlogMthdthn <- .DlogM_Tnegbin_mpfr(mu[wrong],shape)
          d2logMthdth2[wrong] <- dnlogMthdthn$d2logMthdth2
          d3logMthdth3[wrong] <- dnlogMthdthn$d3logMthdth3
          truncGLMweights[wrong] <- GLMweights*(1+d2logMthdth2[wrong]/Vmu) 
        } else {
          if ( ! Rmpfr_warned) {
            message("If the 'Rmpfr' package were installed, better numerical precision would be possible in some Tnegbin computation.")
            Rmpfr_warned <<- TRUE
          }
          truncGLMweights <- pmax(truncGLMweights, .Machine$double.eps)
          # not clear what to do about d3logMthdth3 but it may be quite close to d2logMthdth2
        }
      }
    }
    WU_WT <- GLMweights/truncGLMweights 
    return(list(mu=mu, p0=p0,
                GLMweights=list(truncGLMweights=truncGLMweights,WU_WT=WU_WT,dlogMthdth=dlogMthdth, # this where the w.resid list is created
                                d2logMthdth2=d2logMthdth2, d3logMthdth3=d3logMthdth3)))
  }
})

.calc_GLMweights <- function(LMbool, processed, pw, dmudeta, Vmu) {
  if (LMbool || attr(processed$models,"unit_GLMweights") ) { # really needing unit weights and not simply constant ones here.
    GLMweights <- pw
    GLMweights[] <- eval(pw) # a way of keeping attributes of pw in GLMweights, including "is_unit"
  } else {
    # in Gamma() (inverse) case dmudeta=Vmu => dmudeta^2 /Vmu= Vmu
    # this can diverge: add a warning/protection ? (__F I X M E__ need examples)
    GLMweights <- eval(pw) * dmudeta^2 /Vmu ## must be O(n) in binomial cases
    attr(GLMweights, "unique") <- FALSE ## might actually be true sometimes
    attr(GLMweights, "is_unit") <- FALSE
  }
  GLMweights
}

# elements to be added to the .muetafn output for LLgeneric method:
.add_Md_logcLdeta_terms <- function(muetablob, family, y, mu, pw, dmudeta, eta, phi) { # obsInfo code:standard code for LLgeneric, but also used for truncated non-generic code  

  dlogLdmu <- family$DlogLDmu(mu=mu,y=y,wt=pw, phi=phi)
  muetablob$dlogcLdeta <- dlogLdmu*dmudeta  #                                                               1st
  d2logLdmu2 <- family$D2logLDmu2(mu=mu,y=y,wt=pw, phi=phi) 
  d2mudeta2 <- family$D2muDeta2(eta)
  muetablob$Md2logcLdeta2 <- - (d2logLdmu2*dmudeta^2+dlogLdmu*d2mudeta2) #                                  2nd 
  d3logLdmu3 <- family$D3logLDmu3(mu=mu,y=y,wt=pw, phi=phi) 
  d3mudeta3 <- family$D3muDeta3(eta)
  # high precision necessary for d3logLdmu3 as the dmudeta^3 weight may be largest 
  muetablob$Md3logcLdeta3 <- - drop(d3logLdmu3*dmudeta^3 + 3* d2logLdmu2*d2mudeta2*dmudeta + dlogLdmu*d3mudeta3)
  muetablob
}

.binom_add_Md_logcLdeta_terms <- function(muetablob, family, y, mu, 
                                          dmudeta, # muFREQS !
                                          eta, BinomialDen) { 
  
  y <- drop(y)
  muFREQS <- mu/BinomialDen 
  dlogLdmu <- family$DlogLDmu(muCOUNT=mu, muFREQS=muFREQS , y=y,BinomialDen=BinomialDen)
  muetablob$dlogcLdeta <- dlogLdmu*dmudeta #                                                                1st
  d2logLdmu2 <- family$D2logLDmu2(muFREQS=muFREQS, y=y,BinomialDen=BinomialDen) 
  d2mudeta2 <- family$D2muDeta2(eta)
  muetablob$Md2logcLdeta2 <- - (d2logLdmu2*dmudeta^2+dlogLdmu*d2mudeta2) #                                  2nd                
  d3logLdmu3 <- family$D3logLDmu3(muFREQS=muFREQS, y=y,BinomialDen=BinomialDen) 
  d3mudeta3 <- family$D3muDeta3(eta)
  muetablob$Md3logcLdeta3 <- - drop(d3logLdmu3*dmudeta^3 + 3* d2logLdmu2*d2mudeta2*dmudeta + dlogLdmu*d3mudeta3)
  muetablob
}


.COMP_add_Md_logcLdeta_terms <- function(muetaenv, family, y, mu=muetaenv$mu, 
                                         dmudeta=muetaenv$dmudeta, eta=muetaenv$sane_eta) {

  derivs <- muetaenv$thetaMuDerivs_3 # 3rd deriv needed in family$D3logLDmu3
  dlogLdmu <- family$DlogLDmu(mu=mu,y=y, thetaMuDerivs=derivs)                                        
  muetaenv$dlogcLdeta <- dlogLdmu*dmudeta #  rename as dclikdeta ?                                          1st 
  d2logLdmu2 <- family$D2logLDmu2(mu=mu,y=y, thetaMuDerivs=derivs) 
  d2mudeta2 <- family$D2muDeta2(eta)
  muetaenv$Md2logcLdeta2 <- - (d2logLdmu2*dmudeta^2+dlogLdmu*d2mudeta2) #                                   2nd
  d3logLdmu3 <- family$D3logLDmu3(mu=mu,y=y, thetaMuDerivs=derivs) 
  d3mudeta3 <- family$D3muDeta3(eta)
  # high precision necessary for d3logLdmu3 as the dmudeta^3 weight may be largest 
  muetaenv$Md3logcLdeta3 <- - drop(d3logLdmu3*dmudeta^3 + 3* d2logLdmu2*d2mudeta2*dmudeta + dlogLdmu*d3mudeta3)
  muetaenv
}

.muetafn_LLgeneric <- function(family, LMbool, processed, pw, dmudeta, eta, mu, BinomialDen, y, phi_est) {
  # The <stats:: family>2LLF functions have set LLgeneric to TRUE, so they must have added the member functions required by .add_Md_logcLdeta_terms() 
  
  # if ( family$family %in% c("gaussian","Gamma")  ) { 
  #   GLMweights <- .calc_GLMweights(LMbool, processed, pw, dmudeta, Vmu)
  #   # this was previously needed for w.resid gradient weights used in .hatvals2std_lev -> dh0deta, thus for GLMs with phi param
  #   # but now I use dlogcLikdeta there in all distinct obsInfo cases.
  # } else 
    GLMweights <- "Not supposed to use GLMweights"  # (or GLMs without phi param) for devel : negbin2() may be declared as not using GLM methods. =
    #   # _F I X M E_ try to use unit_GLMweights concept for dlogcLikdeta ? would be limited though (not even Gamma(log))
    
  if ( family$family=="binomial") {
    mu <- mu * BinomialDen
    muetablob <- list(mu=mu, 
                      dmudeta=dmudeta * BinomialDen, # muCOUNTS here
                      GLMweights=GLMweights, sane_eta=eta)
    muetablob <- .binom_add_Md_logcLdeta_terms(muetablob=muetablob, family, y, mu, 
                                               dmudeta, # muFREQS here 
                                               eta, BinomialDen=BinomialDen)
  } else {
    muetablob <- list(mu=mu, dmudeta=dmudeta, GLMweights=GLMweights, sane_eta=eta)
    # p0 computation as in .muetafn_truncated_GLM:
    if (identical(family$zero_truncated,TRUE)) { 
      if (family$family=="poisson") { 
        ## D[Log[1 - E^-E^theta], {theta, 2}] /. {theta -> Log[mu]} // Simplify
        muetablob$p0 <- 1/exp(pmin(mu,log(.Machine$double.xmax*1e-8)/2))
      } else if (family$family=="negbin2") { ## computations in TruncatedNegBin.nb
        ## D[Log[1 - E^-E^theta], {theta, 2}] /. {theta -> Log[mu...]} // Simplify
        shape <- .get_family_par(family)
        muetablob$p0 <- .negbin2_p0(mu,shape) ## (1-p)^r
      } 
    } 
    muetablob <- .add_Md_logcLdeta_terms(muetablob=muetablob, family, y, mu, pw, dmudeta, eta, phi = phi_est)
  }
  # With the derivatives added by the ...add_Md_logcLdeta_terms(),  ..calc_dlW_deta() uses its LLgeneric code and ignores the LLadhoc code 
  return(muetablob) 
  
}


.muetafn_expInfo <- function(family, eta, Vmu, BinomialDen, mu, dmudeta, LMbool, processed, pw) {
  if (family$family=="binomial") {
    if ( family$link=="probit") { ## _F I X M E_ other links need correction too (same Vmu function of mu), but eta-mu link not even symmetric for cloglog
      islarge <- abs(eta)>8
      etalarge <- eta[islarge]
      Vmu[islarge] <- pnorm(etalarge,lower.tail = FALSE)*pnorm(etalarge)
    }
    Vmu <- Vmu * BinomialDen 
    mu <- mu * BinomialDen
    dmudeta <- dmudeta * BinomialDen
  } 
  GLMweights <- .calc_GLMweights(LMbool, processed, pw, dmudeta, Vmu)
  if (identical(family$zero_truncated,TRUE)) { 
    muetablob <- .muetafn_truncated_GLM(family, mu, GLMweights, Vmu) # list(mu, p0, GLMweights=list(truncGLMweights, WU_WT, dlogMthdth, d2logMthdth2, d3logMthdth3))
    muetablob$dmudeta <- dmudeta
    muetablob$sane_eta <- eta
  } else muetablob <- list(mu=mu,dmudeta=dmudeta,GLMweights=GLMweights, sane_eta=eta)
  return(muetablob) 
  
}

.muetafn <-   function(eta, BinomialDen, processed, family=processed$family, 
                       #LMbool=.is_LM(family),  # reevaluates in a backward compatible way is needed # but not called post-fit ?
                       LMbool=family$flags$LMbool,
                       pw=processed$prior.weights, 
                       vec_nobs=processed$vec_nobs,
                       y=processed$y,
                       phi_est # added in v3.13.2 for Gamma and gaussian obsInfo. (-> .add_Md_logcLdeta_terms() -> D*logLdmu* functions). Not quite elegant. 
                       ) { ## note outer var BinomialDen 
  ### if ( ! is.null(names(eta))) stop(" ! is.null(names(eta))")
  # names(eta) <- NULL ## no longer useful because rownames(X.pv) <- NULL and rownames(ZAL) <- NULL
  if (! is.null(vec_nobs)) {
    cum_nobs <- attr(processed$families,"cum_nobs")
    mu <- dmudeta <- sane_eta <- numeric(tail(cum_nobs,1L))
    mv <- vector("list",length(vec_nobs))
    GLMweights <- structure(vector("list",length(vec_nobs)), class="mvlist")
    any_neg_eta <- any_tiny_eta <- NULL
    for (mv_it in seq_along(vec_nobs)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      fam_it <- processed$families[[mv_it]]
      mv[[mv_it]] <- muetablob <- .muetafn(eta[resp_range],BinomialDen[resp_range],processed, family=fam_it, 
                                           LMbool=fam_it$flags$LMbool, pw=pw[[mv_it]], y=y[resp_range],
                                           vec_nobs=NULL, phi_est=phi_est[[mv_it]])
      mu[resp_range] <- muetablob$mu
      dmudeta[resp_range] <- muetablob$dmudeta
      # basic GLMweights may be a list. And  it must have a 'unique' attr. Cf .calc_w.resid()
      GLMweights[[mv_it]] <- muetablob$GLMweights 
      sane_eta[resp_range] <- sane_eta_it <- muetablob$sane_eta
      any_neg_eta <- c(any_neg_eta, attr(sane_eta_it,"any_neg_eta"))
      any_tiny_eta <- c(any_tiny_eta, attr(sane_eta_it,"any_tiny_eta"))
    }
    attr(sane_eta,"any_neg_eta") <- any(any_neg_eta)
    attr(sane_eta,"any_tiny_eta") <- any(any_tiny_eta)
    return(list(mu=mu,dmudeta=dmudeta,GLMweights=GLMweights, sane_eta=sane_eta, mv=mv))
    # : a redundant object with a list 'mv' of sub-muetablob's added to a synthetic muetablob with itself a list of sub-GLMweights'
  } 
  ############### single family code:
  eta <- .sanitize_eta(eta, family=family) #, bin_mu_tol=processed$envir$bin_mu_tol)
  if (family$family=="COMPoisson") {
    muetaenv <- .CMP_muetaenv(family, pw, eta)
    if (processed$how$obsInfo) muetaenv <- .COMP_add_Md_logcLdeta_terms(muetaenv=muetaenv, family=family, y=y)
    return(muetaenv) # possibly including obsInfo extra definitions 
  } 
  
  if (family$family == "binomial") { ## 
    #for binomial cases stats::binomial(...)$linkinv corrects eta for all links so that mu is always .Machine$double.eps from 0 or 1
    # This may not be enough... or inconsistent with corrections used elsewhere, so we overcome the stats:: correction
    # by computing eta <- .binomial_corr_eta(eta,link=family$link, tol=.spaMM.data$options$bin_mu_tol) in .sanitize_eta() and by:
    mu <- .binomial_raw_linkinv(eta,link=family$link)
  } else if (family$family == "poisson") { ## same troubles as for binomial
    mu <- .poisson_raw_linkinv(eta,link=family$link)
  } else mu <- family$linkinv(eta) ## linkinv(eta) is FREQS for binomial, COUNTS for poisson...
  # if (any(is.infinite(mu))) stop()
  dmudeta <- family$mu.eta(eta) 
  Vmu <- family$variance(mu) 
  
  if (processed$how$obsInfo) {
    if (family$flags$LLgeneric) { # routine obsInfo case. 
      # The <stats:: family>2LLF functions have set LLgeneric to TRUE, so they must have added the member functions required by .add_Md_logcLdeta_terms() 
      # With the derivatives added by .muetafn_LLgeneric -> ...add_Md_logcLdeta_terms(),  
      #       ..calc_dlW_deta() uses its LLgeneric code and ignores the LLadhoc code 
      muetablob <- .muetafn_LLgeneric(family, LMbool, processed, pw, dmudeta, eta, mu, BinomialDen, y, phi_est)
    } else { #  LL  ad-hoc case. To maintain functionality of code for experimental obsInfo via H_ratio factors
      muetablob <- .muetafn_expInfo(family, eta, Vmu, BinomialDen, mu, dmudeta, LMbool, processed, pw) 
    }
  } else { # 
    muetablob <- .muetafn_expInfo(family, eta, Vmu, BinomialDen, mu, dmudeta, LMbool, processed, pw)
  }
  return(muetablob)
} ## end def .muetafn

.updateWranef <- function(rand.family,lambda,u_h,v_h) {
  dudv <- rand.family$mu.eta(v_h) ## general cf InvGamma with log link rand.family$mu.eta(v_h) = exp(v) =u is du/d(log u)   
  ## compute w.ranef := - d^2 log dens(v)/dv^2 := 1/Sigma^2_v (= 1/lambda for LMM). See Appendix 3 of LeeN01 + my notes
  ## computed either directly or as (dudv/V_M)*(dudv/lambda)
  ## compute dlogWran_dv_h := d log w.ranef/dv
  if (rand.family$family=="gaussian") {
    if (rand.family$link=="identity") {
      V_M <- rand.family$variance(u_h) ##rep(1,length(u_h)) ## GLMMs in general
      dlogWran_dv_h <- rep(0L,length(u_h))
    }
  } else if (rand.family$family=="Gamma") { 
    if (rand.family$link=="log") {
      V_M <- u_h ## V(u), canonical conjugate Gamma as in canonical Poisson Gamma HGLM
      dlogWran_dv_h <- rep(1L,length(u_h))
    } else if (rand.family$link=="identity") { ## gamma(identity)
      w.ranef <- as.numeric((1-lambda)/(lambda * u_h^2)) ## vanishes for lambda=1 and negative above... (in which case the Laplace approx is bad anyway)
      dlogWran_dv_h <- -2/as.numeric(u_h)
      return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))  ###### return here !
    } 
  } else if (rand.family$family=="inverse.Gamma") { ## for Gamma HGLM 
    ## the canonical form gives the density of theta(u)
    if (rand.family$link=="log") {
      w.ranef <- as.numeric(1/(u_h * lambda)) ## W1/lambda, W1 computation shown in appendix 3 of LeeN01; also in Noh and Lee's code.
      dlogWran_dv_h <- rep(-1L,length(u_h)) ## v=log u, dlogW/dv= dlog(1/u)/dv=-1
      return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))  ###### return here !
    } else if (rand.family$link=="-1/mu") {
      ## D[\[Nu] (\[Theta][u] - (-Log[-\[Theta][u]])), {\[Theta][u], 2}]
      V_M <- rand.family$variance(u_h) ## u_h^2 ## V(u), canonical conjugate HGLM 
      dlogWran_dv_h <- 2 * u_h ## no independent check 
    }
  } else if (rand.family$family=="Beta") {
    if (rand.family$link=="logit") {
      V_M <- rand.family$variance(u_h) ##  u_h*(1-u_h) ## canonical conjugate HGLM
      dlogWran_dv_h <- 1 - 2 * u_h ## D[Log[u (1 - u)] /. u -> 1/(1 + E^-v), v] /. v -> Log[u/(1 - u)] ; no independent check
    }
  }
  ## dudv/V_M may be 1 as both diverge: 
  w.ranef <- as.numeric((dudv/V_M)*(dudv/lambda)) ## semble valide quand v=g(u) = th(u): not previous return()
  w.ranef[w.ranef>1e10] <- 1e10 ## patch useful to avoid singular d2hdv2 in PLoG model
  return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))
}


.updateW_ranefS <- function(cum_n_u_h, rand.families, lambda, u_h, v_h, 
                            w.ranef_list=vector("list", length(rand.families)), 
                            dlogWran_dv_h_list=vector("list", length(rand.families)), 
                            dvdu_list=vector("list", length(rand.families))) {
  nrand <- length(rand.families)
  for (it in seq_len(nrand)) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    blob <- .updateWranef(rand.family=rand.families[[it]],lambda[u.range],u_h[u.range],v_h[u.range])
    w.ranef_list[[it]] <- blob$w.ranef
    dlogWran_dv_h_list[[it]] <- blob$dlogWran_dv_h
    dvdu_list[[it]] <- blob$dvdu
  }
  resu <- list(w.ranef=.unlist(w.ranef_list),dlogWran_dv_h=.unlist(dlogWran_dv_h_list),dvdu=.unlist(dvdu_list))
  ## the test is invalid for ranCoefs:
  # if (nrand==1L && rand.families[[1L]]$family=="gaussian") resu$unique_w.ranef <- w.ranef[[1L]] # used in sparsePrecision code
  #if (length(unique_w.ranef <- unique(w.ranef))==1L) resu$unique_w.ranef <- unique_w.ranef # used in sparsePrecision code
  resu
}

.calc_d2mudeta2 <- function(link,mu=NULL,eta=NULL,muFREQS=NULL) { ## d2 MuCOUNTS d etaFREQS^2 
  # BinomialDen removed for same reason as in .thetaMuDerivs()
  switch(link,
         identity = 0,
         log = mu, 
         inverse = 2 * mu^3 , ## canonical for Gamma()
         ## next three make sense for Binomial response data
         logit = d2muFREQS <- muFREQS*(1-muFREQS)*(1-2*muFREQS),
         probit = -eta * dnorm(eta), # assuming eta corrected as in binomial(probit)$mu.eta
         cloglog = {
           eta <- pmin(eta,700) ## as in binomial(cloglog)$mu.eta
           ## in particular d2mudeta2/dmudeta is then numerically OK
           exp(eta-exp(eta))*(1-exp(eta)) ## D[1 - E^-E^\[Eta], {\[Eta], 2}]
          }, 
         cauchit = -2 *eta/(pi * (1+eta^2)^2),
         sqrt = 2, 
         stop(paste("unhandled link'",link,"'in .calc_d2mudeta2()"))
  )
} 

#derivatives of GLM weights wrt eta 
.CMP_calc_dlW_deta_locfn <- function(i,lambdas,mu,COMP_nu) { # (only used for canonical link case)
  lambdai <- lambdas[[i]]
  if (lambdai==0) {
    return(c(dmudeta=0, d2mudeta2=0))
  }
  ## ELSE
  moments <- .COMP_Pr_moments(lambda=lambdai,nu=COMP_nu,moments=c("2","3"))
  mui <- mu[[i]]
  dmu.dlogLambda <- moments[["2"]] - mui^2 # =family$mu.eta() without repeating some computations    # ...that's the family $variance()...
  d2mu.dlogLambda2 <- moments[["3"]]-mui*moments[["2"]]-2*mui*dmu.dlogLambda
  return(c(dmudeta=dmu.dlogLambda, d2mudeta2=d2mu.dlogLambda2))
}

# returns bounded eta
.binomial_corr_eta <- function(eta, link=family$link, tol=.Machine$double.eps) {
  switch(link,
         logit= { # .Call(C_logit_linkinv, eta) de facto corrects eta has follows before computing mu:
           thresh <- -qlogis(tol)
           .bracket(eta, -thresh, thresh) #pmin(pmax(eta, -thresh), thresh)
         },
         probit= {
           ### .muetafn has called binomial(probit)$linkinv to compute mu(FREQS), which first corrected eta as follows:
           thresh <- -qnorm(tol)
           .bracket(eta, -thresh, thresh) # pmin(pmax(eta, -thresh), thresh) # so that abs(dnorm(eta)) is bounded by .Machine$double.eps   
         },
         cloglog= {
           lo <- log(-log1p(-tol))
           up <- log(-log(tol))
           .bracket(eta, lo, up) #pmin(pmax(eta, lo), up) 
          },
         cauchit= {
           thresh <- -qcauchy(tol)
           .bracket(eta, -thresh, thresh) # pmin(pmax(eta, -thresh), thresh)
          },
         stop("Unhandled link")
  )
}
## _F I X M E_ disjunction between previous and next function may cause pbs
.binomial_raw_linkinv <- function(eta, link=family$link) { # without control of eta extreme values 
  switch(link,
         logit= 1/(1+exp(-eta)),
         probit= pnorm(eta),
         cloglog= -expm1(-exp(eta)),
         cauchit= pcauchy(eta),
         stop("Unhandled link")
  )
}

.poisson_raw_linkinv <- function(eta, link=family$link) { # without control of eta extreme values 
  switch(link,
         log=exp(eta),
         identity= eta,
         sqrt= eta^2,
         stop("Unhandled link")
  )
}

.D2detadmuDeta2 <- function(link, mu) { # in Hessian_weights.nb
  switch( #  d (d (deta/dmu) / deta2
    link,
    "log" =  1/mu,
    "inverse" = -2, # in the Gamma(inverse) sense 
    "identity" = 0,
    "sqrt" = 1/mu^(1.5),
    "cloglog" = {(1 + log(1 - mu) + log(1 - mu)^2)/((-1 + mu)*log(1 - mu))}, 
    "logit" = {(1 - 2*mu + 2*mu^2)/(mu - mu^2)},
    "probit" = {
      eta <- qnorm(mu)
      (1+eta^2)/dnorm(eta)
    },
    "cauchit" = 2*pi,
    # loglambda not necess bc appears only as canonical link
    stop("This cannot yet be computed for this family link.") 
  )
}

.dlog_HratioDeta <- function(Hratio_factors, mu, family, 
                             # binomial:
                             BinomialDen, muFREQS=mu/BinomialDen) {
  yMmu <- Hratio_factors$yMmu
  partials <- Hratio_factors$partials
  if (family$family=="binomial") {
    # this blcok should be correct in the more general case, but then with useless divisins by vectors of 1 
    dmu.deta <- family$mu.eta(eta=family$linkfun(muFREQS)) # as function of this fn argument, mu, and distinct from muetablob$dmudeta
    DfactorsDeta  <- .D2detadmuDeta2(family$link, mu=muFREQS) + .D_dlWdmu_Deta(family, mu=muFREQS) # 
    # Role of BinomialDen as explained in .calc_H_w.resid()
    yMmuFREQS <- yMmu/BinomialDen
    dHratio <- (dmu.deta * partials- DfactorsDeta * yMmuFREQS)
    Hratio <- 1 - partials * yMmuFREQS 
  } else {
    dmudeta <- family$mu.eta(eta=family$linkfun(mu)) # as function of mu (would muetablob be accessible?)
    DfactorsDeta  <- .D2detadmuDeta2(family$link, mu=mu) + .D_dlWdmu_Deta(family, mu=mu)  
    dHratio <- (dmudeta * partials- DfactorsDeta * yMmu)
    Hratio <- 1 - partials * yMmu
  }
  dHratio / Hratio
} # we will need to add this to the dlW_deta term in coef2

## Called for LLM- and GLM-family objects. Returns derivatives of W_Hobs when the objective is Laplace-obs, 
## whether for LLM or for GLM with obsInfo implemented through Hratio factors, for untruncated only.
#### meaning of code:
# if (LLgeneric) {
#   dlW_deta <- is directly dlW_obs/deta
# } else if ( ! is.null(dlW_Hexp__detafun <- family$dlW_Hexp__detafun)) {
#   dlW_deta <- dlW_Hexp__detafun gives dlW_exp/deta 
#   if (obsInfo) {
#     dlW[obs]_deta <- dlW[exp]_deta + correction
#   } else no obsInfo => dlW[exp]_deta is used
# }
..calc_dlW_deta <- function(muetablob, family, calcCoef1, w.resid, BinomialDen, processed, Hratio_factors,
                            #
                            dmudeta=muetablob$dmudeta, mu=muetablob$mu, eta=muetablob$sane_eta,
                            muFREQS=mu/BinomialDen, y=processed$y,
                            obsInfo=processed$how$obsInfo,
                            Hobs_w.resid=muetablob$Hobs_w.resid # added for LLM case
                            ) {
  coef1 <- NULL
  
  # Note that LLFs are always treated as non-canonical link. negbin2 show why this may not always be optimal. 
  # But optimizatiosn are of limited interest in the truncated case.
  # if negbin() is faster than negbin2(), we could arrage for negbin2() to return negbin()? -- Would make debugging more difficult? 
  
  if ( ! is.null(muetablob$Md3logcLdeta3)) {  # obsInfo by LLgeneric method
    res <- list(dlW_deta=muetablob$Md3logcLdeta3/muetablob$Md2logcLdeta2) # coef2 
    if (calcCoef1) res$coef1 <- res$dlW_deta/muetablob$Md2logcLdeta2 # Hobs_w.resid <-  muetablob$Md2logcLdeta2
    return(res) # and this means no $WU_WT is included in the result, contrary to Hexp truncated case
  }  
  
  ## We first handle the canonical link cases, where comput. of coef1 depends only on the link  
  ## here w=dmudeta; d1=dwdmu dmudeta /w^2 = dlogwdeta/w = (d2mu/deta2)/(dmu/deta) /w =
  ##      (d2mu/deta2)/(dmu/deta)^2 = (d(dmudeta)/dmu)/dmudeta where d(dmudeta)/dmu is the numerator as detailed:
  if (family$flags$canonicalLink) {
    #dlW_deta <- d2mudeta2 / dmudeta or :
    if (family$family=="gaussian") {
      if (calcCoef1) coef1 <- rep(0L,length(mu))
      dlW_deta <- rep(0L,length(mu))
    } else if (family$family=="poisson") {
      ## numerator is D[D[E^\[Eta], \[Eta]] /. {E^\[Eta] -> \[Mu]}, \[Mu]] =1 
      if (identical(family$zero_truncated,TRUE)) { ## 
        d_tildeW_deta <- mu + w.resid$d3logMthdth3 ## MolasL10 p 3307
        dlW_deta <- d_tildeW_deta/(w.resid$w_resid) ## it really is dlog_tildeW_deta
        if (calcCoef1) coef1 <- dlW_deta/dmudeta ## coef1 = dl_tildeW_deta/ W_untrunc
      } else { ## coef1 = dlW_deta/ W
        dlW_deta <- rep(1L,length(mu))
        if (calcCoef1) coef1 <- 1/dmudeta
      }
    } else if (family$family=="binomial") {
      dlW_deta <- (1-2*muFREQS)  
      if (calcCoef1) coef1 <- dlW_deta/dmudeta # F I X M E there's a pattern
      ## numerator is D[D[1/(1 + E^-\[Eta]), \[Eta]] /. {E^-\[Eta]->(1-\[Mu])/\[Mu]} ,\[Mu]]=1-2 mu 
    } else if (family$family=="Gamma") { ## link= "inverse" !
      ## numerator is D[D[-1/\[Eta], \[Eta]] /. {\[Eta] -> -1/\[Mu]}, \[Mu]] =2 mu 
      dlW_deta <- 2*mu
      if (calcCoef1) coef1 <- dlW_deta /dmudeta
    } else if (family$family=="COMPoisson") { 
      if (TRUE) { # using the muetablob environment
        dmudeta <- muetablob$dmudeta # = # ...that's the family $variance()...
        d2mudeta2 <- muetablob$EX3-muetablob$EX*muetablob$EX2-2*muetablob$EX*dmudeta
        dlW_deta <- d2mudeta2 / dmudeta
        dlW_deta[is.nan(dlW_deta)] <- 0 ## quick patch for cases that should have low lik 
        if (calcCoef1) {
          coef1 <- dlW_deta / dmudeta
          coef1[is.nan(coef1)] <- 0 ## idem
        }
      } else {
        COMP_nu <- environment(family$aic)$nu
        lambdas <- exp(eta) ## pmin(exp(eta),.Machine$double.xmax) ##  FR->FR lambdas missing as mu attribute here ?
        blob <- sapply(seq(length(lambdas)), .CMP_calc_dlW_deta_locfn,lambdas=lambdas,mu=mu,COMP_nu=COMP_nu)
        dlW_deta <- blob["d2mudeta2",] / blob["dmudeta",]
        dlW_deta[is.nan(dlW_deta)] <- 0 ## quick patch for cases that should have low lik
        if (calcCoef1) {
          coef1 <- dlW_deta / blob["dmudeta",]
          coef1[is.nan(coef1)] <- 0 ## idem
        }
      }
    } 
  } else if (family$family=="binomial" && family$link=="probit") { ## ad hoc non canonical case 
    pmax_dnorm_eta <- dnorm(eta)
    VmuFREQS <- pnorm(eta,lower.tail = FALSE)*pnorm(eta)  # more accurate than muFREQS*(1-muFREQS)
    dlW_deta <- -2*eta - pmax_dnorm_eta*(1-2*muFREQS)/VmuFREQS
    if (obsInfo) {
      dlW_deta <- dlW_deta + .dlog_HratioDeta(Hratio_factors, mu, family, 
                                              muFREQS=muFREQS, BinomialDen=BinomialDen) # and now this is the d log w_Hobs 
      if (calcCoef1) {
        Hobs_w.resid <- ..calc_H_w.resid(w.resid, Hratio_factors=Hratio_factors,BinomialDen=BinomialDen)
        coef1 <- dlW_deta/Hobs_w.resid # (with BinomialDen in denom)
      }
    } else if (calcCoef1) {
      coef1 <- dlW_deta *(VmuFREQS)/ (BinomialDen * pmax_dnorm_eta^2) 
      # or coef1 <- (-1+2*(mu-VmuFREQS*eta/pmax_dnorm_eta))/(BinomialDen*pmax_dnorm_eta)
    }
  } else if (family$family=="Gamma" && family$link=="log") { ## ad hoc non canonical case 
    if (obsInfo) {
      # special case not calling Hratio_factors:
      dlW_deta <- rep(-1,length(mu)) 
      if (calcCoef1) coef1 <- - muetablob$mu/drop(w.resid*y) # should be -1/wi= -1/(w.resid*y/mu)
    } else {
      dlW_deta <- rep(0L,length(mu)) ## because they both involve dW.resid/dmu= 0
      if (calcCoef1) coef1 <- dlW_deta
    }
  } else if (family$family %in% c("negbin2","poisson")) { # 
    dlW_deta <- family$dlW_Hexp__detafun(mu) 
    if (obsInfo) { # by Hratio-factors method
      # UNtruncated => has dlW_Hexp__detafun available. Use it instead of the generic fallback code.
      dlW_deta <- dlW_deta + .dlog_HratioDeta(Hratio_factors, mu, family) # and now this is the d log w_Hobs 
      if (calcCoef1) {
        Hobs_w.resid <- ..calc_H_w.resid(w.resid, Hratio_factors=Hratio_factors,BinomialDen=1)
        coef1 <- dlW_deta/Hobs_w.resid # (with BinomialDen in denom)
      }
    } else { # Hexp ! 
      if (calcCoef1) coef1 <- family$coef1fun(mu)
    }
    # Function not present for truncated models
  } else { # Has been general code for GLM H_exp or Hobs. LLF-family never reach this point, see first lines of the fn.
    # Should become only Hexp for GLMs 
    #
    ## we need to update more functions of mu...
    # Now (2022/05/01) in terms of derivatives of functions of eta and muFREQS so there is no BinomialDen factor in them.
    # Instead yMmu/BinomialDen is used in obsInfo case and in denom of coef1 generally.
    ### (1) compute dlW[H_exp] 
    dmu.deta <- dmudeta/BinomialDen
    d2mu.deta2 <- .calc_d2mudeta2(link=family$link,mu=mu,eta=eta, muFREQS=muFREQS)
    tmblob <- .thetaMuDerivs(mu, family=family, muFREQS=muFREQS, muetablob=muetablob)
    Dtheta.Dmu <- tmblob$Dtheta.Dmu 
    D2theta.Dmu2 <- tmblob$D2theta.Dmu2 
    D2theta.Deta2_Dtheta.Deta <- (D2theta.Dmu2 * dmu.deta^2 + Dtheta.Dmu * d2mu.deta2)/(Dtheta.Dmu * dmu.deta)
    dlW_deta <- d2mu.deta2 / dmu.deta + D2theta.Deta2_Dtheta.Deta
    if (obsInfo) {
      if (identical(family$zero_truncated,TRUE)) {
        stop("This GLM family object cannot be used to fit a zero-truncated model observed-Hessian Laplace approximation.")
      } 
      ### (2) compute dlW[H_obs] 
      dlW_deta <- dlW_deta + .dlog_HratioDeta(Hratio_factors, mu=mu, family=family, 
                                              # next two only for binomial, and provided by the ..calc_dlW_deta() promises. 
                                              muFREQS=muFREQS, BinomialDen=BinomialDen) # and now this is the d log w_Hobs 
      if (calcCoef1) {
        Hobs_w.resid <- ..calc_H_w.resid(w.resid, Hratio_factors=Hratio_factors,BinomialDen=BinomialDen)
        coef1 <- dlW_deta/Hobs_w.resid # (with BinomialDen in denom)
      }
    } else if (calcCoef1) {
      if (identical(family$zero_truncated,TRUE)) {  
        # at this point we have the untruncated dlW_deta. We evaluate the truncated dlW_deta dlog_tildeW_deta
        dlW_deta <- dlW_deta + w.resid$WU_WT *
          (D2theta.Dmu2/Dtheta.Dmu * w.resid$d2logMthdth2 +  Dtheta.Dmu * w.resid$d3logMthdth3) * Dtheta.Dmu* dmu.deta
      } 
      ## in truncated case we have coef1 = (truncated dlW_deta)/ W_untrunc so the following is still correct
      coef1 <- dlW_deta / (BinomialDen* Dtheta.Dmu * dmu.deta^2) ## note that coef2 is indep of the BinomialDen, but coef1 depends on it 
    }
  }
  
  res <- list(dlW_deta=dlW_deta,## dlW_deta equiv coef2
              coef1=coef1)
  if (identical(family$zero_truncated,TRUE)) res$WU_WT <- w.resid$WU_WT # again, LLF-family never reach this point.
                                                                        # For LLF families, this w.resid is not a list.
  return(res) 
}


.calc_dlW_deta <- function(muetablob, 
                           processed, # note use not only as provided of next arguments, in observed-Hessian case
                           calcCoef1=FALSE, 
                           w.resid, # potentially the list with $w_resid element, etc.  
                           Hratio_factors=NULL,
                           family=processed$family, 
                           families=processed$families, 
                           BinomialDen=processed$BinomialDen 
) {
  if ( ! is.null(cum_nobs <- attr(families,"cum_nobs"))) {
    coef1 <- dlW_deta <- numeric(length(muetablob$mu))
    for (mv_it in seq_along(families)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      # in the mv case w.resid is a list only when there are zero-truncated families
      if (is.list(w.resid)) {
        w.resid_it <- w.resid[["mvlist"]][[mv_it]]
      } else w.resid_it <- w.resid[resp_range] 
      dlW_detablob <- ..calc_dlW_deta(muetablob=muetablob$mv[[mv_it]],
                                     family=families[[mv_it]],
                                     BinomialDen=BinomialDen[resp_range],
                                     y=processed$y[resp_range],
                                     processed=processed,
                                     calcCoef1=calcCoef1,
                                     w.resid=w.resid_it,
                                     Hratio_factors=list(Hratio_factors$yMmu[resp_range], 
                                                         partials=Hratio_factors$partials[resp_range])# for obsInfo 
      )
      dlW_deta[resp_range] <- dlW_detablob$dlW_deta
      if (is.null(dlW_detablob$coef1)) {
        coef1[resp_range] <- 0
      } else coef1[resp_range] <- dlW_detablob$coef1
    }
    return(list(dlW_deta=dlW_deta, coef1=coef1))
  } else ..calc_dlW_deta(muetablob=muetablob, family=family, calcCoef1=calcCoef1, w.resid=w.resid, BinomialDen=BinomialDen, 
                  processed=processed, Hratio_factors=Hratio_factors)
}

.safesolve_qr_vector <- function(qr.a,b,silent=TRUE) { ## solve.qr with fall-back; qr.a should be a qr object, b must be a vector
  if (inherits(qr.a,"sparseQR")) { ## pas de 'essai' en var locale !
    ## there was a 'Matrix' subcode prior to 10/03/2013; another try on 11/2013
    res <- qr.coef(qr.a,b)
  } else {
    res <- try(solve.qr(qr.a,b),silent=silent)
    if (inherits(res,"try-error")) {   ## then some weird code, but...
      ## we try to solve(<original 'a' matrix>) using the QR decomp... this may work when solve.qr fails !
      ## The following is equivalent to solveA <- try(solve(qr.R(qr.a)[,pivI])%*%t(qr.Q(qr.a)),silent=silent) then return solveA %*% b 
      ## but uses only backward mutliplication of vectors and transpose of vectors  
      res <- t(crossprod(b, (qr.Q(qr.a)))) ## not yet res
      #pivI <- sort.list(qr.a$pivot) ## inverse perm such as pivI[$pivot]=$pivot[pivI]= identity
      #res <- try(solve(qr.R(qr.a)[, pivI]) %*% res, silent=silent)
      res <- try(backsolve(qr.R(qr.a),res[qr.a$pivot]))
      if (inherits(res,"try-error")) {
        stop("backsolve() failed in .safesolve_qr_vector().") ## perhaps recover A by qr.X and solve(A) ?
      } 
    }  
  }
  return(res)
}

# Matrix needs some help: for
# Z <- Matrix::.symDiagonal(n=10000)
# system.time(replicate(10000, {if (spaMM:::.is_identity(Z)) {Z}})) is *much* faster than
# system.time(replicate(10000, Z %*% Z))
# So the function is useful, as long as isDiagonal() is not called.
# Commented line may be useful to detect inefficient code
.is_identity <- function( x, matrixcheck=FALSE, tol=1e-8 ) {
  if (inherits(x,"Matrix")) {
    if (inherits(x,"ddiMatrix") ) return(x@diag=="U")
    if (matrixcheck) { # which is never TRUE...)
      if (! isDiagonal( x ) ) return( FALSE ) # isDiagonal is terribly slow
      ## hence now diagonal:
      # if (max(abs(range(diag(x))-1))<tol) stop("ICI")
      return(max(abs(range(diag(x))-1))<tol)
    } else return(FALSE)
  } else {
    if (matrixcheck) {
      return(ncol(x)==nrow(x) && max(abs(x- diag(ncol(x))))<tol)
    } else return(FALSE)  
  }
}

.wrap_Ltri_t_chol <- function(mat, use_eigen=.spaMM.data$options$USEEIGEN) { 
  ## returns lower tri as transpose(base::chol) in all cases; 
  if (inherits(mat,"sparseMatrix")) {
    return( t(Matrix::chol(mat))) ## matches Cholesky !
  } else if (inherits(mat,"Matrix")) {
    mat <- as.matrix(mat) 
  }
  if (use_eigen) {
    chol <- .RcppChol(mat) 
    if ( chol$Status==1L) { 
      return(chol$L) 
    } else stop("chol$Status !=1L") ## best used in combination with try()
  } else return(t(chol(mat))) 
} 

.Cholwrap <- .wrap_Ltri_t_chol

.Utri_chol_by_qr <- function(mat) { # transpose of .wrap_Ltri_t_chol
  esys <- .eigen_sym(mat) 
  eigvals <- esys$values
  if (min(eigvals)<0) { # from removed function .regularize_Wattr():
    target_min_d <- max(eigvals)/1e100 ## so that corrected condition number is at most the denominator:     
    # (yes, but with high 'condnum' and large negative esys$values, target_min_d may not be numerically 
    # distinct from target_min_d-esys$values and then output matrix has a zero eigenvalue... hence add small regul_ranCoefs[1L] here:
    regul_ranCoefs <- .spaMM.data$options$regul_ranCoefs
    d_corr <- max(c(0,(target_min_d-eigvals)*(1+regul_ranCoefs[1L])))
    ## ... but then target_min_d loses its verbatim meaning...
    eigvals <- eigvals + d_corr # all diag is corrected => added a constant diagonal matrix to mat
  } ## .smooth_regul() implements another approach
  crossfac <- .Dvec_times_matrix(sqrt(eigvals), t(esys$vectors)) # crossfac but not upper triangular.
  qrblob <- qr(crossfac)
  crossfac <- qr.R(qrblob) # applying .lmwithQR() systematically (badly) affects numerical precision
  if (! all(unique(diff(qrblob$pivot))==1L)) { # eval an unpermuted triangular R
    crossfac <- .lmwithQR(crossfac[, sort.list(qrblob$pivot)] ,yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled # not permuted in contrast to qr(crossfac) # upper tri crossfac
  } 
  crossfac <- .Dvec_times_matrix(sign(.diagfast(crossfac)), crossfac) ## upper.tri crossfac
  return(crossfac)
}

.wrap_Utri_chol <- function(mat, use_eigen=.spaMM.data$options$USEEIGEN) { # transpose of .wrap_Ltri_t_chol
  ## returns upper tri crossfactor as base::chol in all cases
  if (inherits(mat,"sparseMatrix")) {
    return(chol(mat)) ## Matrix::chol
  } else if (inherits(mat,"Matrix")) {
    ## this was the typical case when crossr22 in .calc_r22 was a dsyMatrix (dense symmetric)
    # but we now avoid costly Matrix operations that yields crossr22 as a dsyMatrix
    mat <- as.matrix(mat) 
  } 
  if (use_eigen) {
    chol <- .Rcpp_chol_R(mat) 
    if (chol$Status==1L) { 
      return(chol$R) 
    } else { ## tested by tests_private/test-for-scaled-spprec.R
      ## OK for small matrices such as crossr22
      crossfac <- .Utri_chol_by_qr(mat) # not permuted in contrast to qr(crossfac) # upper tri crossfac
      return(crossfac)
    }
  } else try(chol(mat)) 
} 




#LogAbsDetWrap <- function(...) .LogAbsDetWrap(...) ## 

.LogAbsDetWrap <- function(mat,logfac=0) { ## M or m
  if (ncol(mat)==0) return(0) ## GLM fitted by ML: d2hdbv2 is 0 X 0 matrix 
  # un piege est que mat/(2*pi) conserve les attributes de mat (telle qu'une décomp QR de mat...)
  # il nefaut  dont pas demander LogAbsDetWrap(mat/(2*pi))
  if (inherits(mat,"Matrix")) {
    lad <- Matrix::determinant(mat)$modulus[1]
  } else if (.spaMM.data$options$USEEIGEN) {
    lad <- .LogAbsDetCpp(mat)
  } else lad <- determinant(mat)$modulus[1]
  # pb general est cert eigenvalues peuvent -> +inf et d'autres -inf auquel cas logabsdet peut être innocuous mais pas estimaable précisément   
  if (is.nan(lad) || is.infinite(lad)){## because of determinant of nearly singular matrix
    zut <- abs(eigen(mat,only.values = TRUE)$values) 
    zut[zut<1e-12] <- 1e-12
    lad <- sum(log(zut)) 
  }
  lad <- lad + nrow(mat)*logfac
  return(lad)
}

# currently not used:
.Matrix_tcrossprod <- function(x,y, as_sym=TRUE) { # x Matrix; y Matrix or NULL
  if (inherits(x,"dgCMatrix")) {
    if (is.null(y)) {
      xyt <- .dgCtcrossprod(x,NULL) 
      if (as_sym) {
        return(forceSymmetric(xyt)) ## not as(.,"sparseMatrix") which checks -> slow
      } else return(xyt)    
    } else if (inherits(y,"dgCMatrix")) {
      return(.dgCtcrossprod(x,y)) ## return(Matrix::tcrossprod(x,y)) ## 
    } else return(Matrix::tcrossprod(x,y))
  } else return(Matrix::tcrossprod(x,y))
}

.tcrossprod <-  function(x, y=NULL, chk_sparse2mat=TRUE, as_sym=TRUE, 
                         perm # for the dCHMsimpl case
                         ) { # XX' tcrossprod, or more general concept for dCHMsimpl
  # this function can be used to compute a cov matrix from its tcrossfac;
  # but also to compute a permuted cov matrix from the dCHMsimpl of an unpermuted precision matrix.
  # Hence the ad hoc code in the latter case, and the absence of default value of 'perm' argument, as a default could cause hard to track bugs. 
  if (is.null(x)) return(NULL) ## allows lapply(,.tcrossprod) on a listof (matrix or NULL)
  if (inherits(x,"ZAXlist")) {
    return(tcrossprod(x,y)) # ZAXlist method
  } else if (inherits(x,"dCHMsimpl")) { # This case is tricky as we do not strictly want a tcrossprod...
    if (perm) {
      if (is.null(y)) {
        resu <- .tcrossprod(solve(x, system="Lt", b=.sparseDiagonal(n=ncol(x), shape="g"))) # cov matrix for permuted ranefs from L=chol_Q  # or Matrix::chol2inv(t(as(x,"sparseMatrix"))) 
      } else resu <- t( solve(x, b=t(y), system="L") ) # possibly never called 
    } else { # previous code:
      if (is.null(y)) {
        resu <- solve(x, system="A") # cov matrix for unpermuted ranefs from (L=chol_Q,P) factorization of its inverse
      } else resu <- t( solve(x, b=as(x,"pMatrix") %*% t(y), system="L") ) # old code, # possibly never called, effectively assuming 'perm'=FALSE
    }
  } else if (inherits(x,"Matrix") || inherits(y,"Matrix")) {
    if (is.null(y)) {
      if (inherits(x,"dgCMatrix")) { 
        if (as_sym) {
          # xxt <- .dgCtcrossprod(x,y)  # might be faster than Matrix::tcrossprod() in some cases, 
          # but not consistently so + we need forceSymmetric() + less memory efficient.
          # return(forceSymmetric(xxt)) ## not as(.,"sparseMatrix") which checks -> slow
          return(Matrix::tcrossprod(x)) # directly dsC
        } else return(.dgCtcrossprod(x,y)) # cf case where we want to avoid "...BLOB$Gmat + dampdG is dsC + dsC..."  ... except that this may return dsC below
      } else if (chk_sparse2mat && inherits(x,"sparseMatrix") && (maxd <- prod(dim(x)))>1L ) { 
        x_denseness <- .calc_denseness(x)/maxd 
        if (x_denseness>0.35) {
          resu <- .tcrossprodCpp(as.matrix(x),y) ##  so much faster (simulate(mrf,...) )-- poor Matrix:: coding ?
          if (as_sym) resu <- as(resu,"symmetricMatrix") # but the as(, ) may produce a quite dense dsC (or dsy)
          #                                                             
          colnames(resu) <- rownames(resu) <- rownames(x)
          return(resu)
        } else return(Matrix::tcrossprod(x)) # .... may well be dsC
      } else return(Matrix::tcrossprod(x))
#    } else if (FALSE && inherits(x,"dgCMatrix") &&  inherits(y,"dgCMatrix") ) { ## FASTer when (FALSE &&)
#      xyt <- .dgCtcrossprod(x,y) 
#      return(xyt)
    } else if (chk_sparse2mat) { 
      if (inherits(x,"sparseMatrix") && (maxd <- prod(dim(x)))>1L) { # trivial >1L bc the alternatives to Matrix::tcrossprod are always useful
        x_denseness <- .calc_denseness(x)/maxd 
        if (x_denseness>0.35) x <- as.matrix(x)
      } #else x_denseness <- 0 ## ie if the test is FALSE, no explicit conversion is attempted
      if (inherits(y,"sparseMatrix") && (maxd <- prod(dim(y)))>1L) { # idem
        y_denseness <- .calc_denseness(y)/maxd
        if (y_denseness>0.35) y <- as.matrix(y)
      } #else y_denseness <- 0 ## idem
      if (is.matrix(x) && is.matrix(y)) { # sparsity overhead is ~ 1/0.35
        return(.tcrossprodCpp(x, y)) 
      } else return(Matrix::tcrossprod(x,y)) ## both <0.35
    } else return(Matrix::tcrossprod(x,y))
  } else {
    if (is.integer(x)) storage.mode(x) <- "double"
    if (is.integer(y)) storage.mode(y) <- "double"
    resu <- .tcrossprodCpp(x,y)
    # Tried calling directly .tcrossprodCpp and ignoring the names for dense hatval_Z but gains are negligible
    if (is.null(y)) {
      colnames(resu) <- rownames(resu) <- rownames(x)
    } else {
      rownames(resu) <- rownames(x)
      colnames(resu) <- rownames(y)
    }
    return(resu)
  }
}

.Matcrossprod <- function(x,y=NULL,return_mat=FALSE) { # wraps an inefficient bit of code, to be avoided 
  if (is.null(y)) {
    return(Matrix::crossprod(x)) 
  } else if (return_mat) {
    return(as.matrix(Matrix::crossprod(x,y))) # that's the ineff code, not called in the long tests.
  } else {
    return(Matrix::crossprod(x,y)) # perfectly OK and routine () ...except...
    # the time of Matrix::crossprod(x,<dge>, y) seems eqivalent to that of crossprod(as.matrix(x),y) so the question is whether as.matrix(x) should be prevcomputed or not
  }
}

## .crossprod() designed to handle all cases. It may call:
#  .crossprodCpp_d(dense) which is remarkably fast when y=NULL at least (what about y!=NULL?) so method of choice for as_mat=TRUE
#  .Rcpp_crossprod() for most cases except vector y. .Rcpp_crossprod calls itself .crossprod_not_dge() after conversion of dge
#  For remaining vector y and unanticipated cases, .dgCcrossprod is called (for dgC matrices), or 
#         the above wrapper for Matrix::crossprod(x,y), or .crossprodCpp_d() again.  



# as_mat=TRUE and as_sym=TRUE may return inconsistent results.
# It is distinctly better to call directly .Rcpp_crossprod() if we know in advance that it will be called by .crossprod 
.crossprod <- function(x, y=NULL, 
                       allow_as_mat=TRUE, # controls chk_sparse2mat AND an automatic conversion.
                       chk_sparse2mat=allow_as_mat, # means that .calc_denseness may be called.on a sparseMatrix
                       # If so and if the Matrix is really sparse and we waste the tiem of checking that;
                       # is not chk_sparse2mat, Matrix::crossprod is called. So as much as possible we should JOINTLY
                       # let x be sparseMatrix only when it is truely sparse; and set chk_sparse2mat=FALSE  
                       as_sym=( ! as_mat && is.null(y)), use_Rcpp_cp=.spaMM.data$options$Rcpp_crossprod, 
                       eval_dens=TRUE, as_mat=FALSE,
                       keep_names=TRUE # for .Rcpp_crossprod call: this fn keeps names by default, contrary to others, so for consistency with others, it may be useful to reverse the default.
                       ) {
  #if (inherits(x,"dgeMatrix")) stop("Possibly inefficient use of 'dgeMatrix' caught by .calc_denseness.") 
  if (is.null(x)) return(NULL) ## allows lapply(,.tcrossprod) on a list of (matrix or NULL)
  if (inherits(x,"ZAXlist")) {
    return(crossprod(x,y)) # ZAXlist method
  } else if (as_mat && is.null(y)) {
    return(.crossprodCpp_d(as.matrix(x),yy=NULL)) # fast relative to .Rcpp_crossprod(sparse)... probably not general !
  } else if (use_Rcpp_cp) {
    if ( inherits(x,c("dgCMatrix","dgeMatrix","matrix"))
         && inherits(y,c("dgCMatrix","dgeMatrix", "matrix","NULL"))) { ## all cases handled by .Rcpp_crossprod()
      txy <- .Rcpp_crossprod(x, y, eval_dens, as_mat=as_mat, keep_names=keep_names) # returns 'matrix' or 'dgCMatrix' (the latter for sparse input not converted to dens after eval_dens check); 
      #                                         next forceSymmetric() can produce dgC (sp) OR dsy (dense)
      if (as_sym) { txy <- forceSymmetric(txy) } # absolutely necessary for use as .updateCHMfactor(., parent=txy) with a (t)crossprod 'parent'; see ?`CHMfactor-class`
      #                                            BUT:  computing the (t)crossfac in the first place is not needed.
      if ( ( ! allow_as_mat) && ( ! inherits(txy,"sparseMatrix"))) txy <- as(txy,"sparseMatrix") # fortunately rarely needed but occurs in HLCor(lh~1 +AR1(1|time) example
      return(txy)
    } else { ## show unhandled classes (ddi...) and continue (but vector y reaches here - OK - and eval_dens is ignored - OK - )
      if (use_Rcpp_cp>1L) {
        print(c("-",class(x),class(y))) 
        utils::flush.console()
      }
    }
  }
  # vector y reaches here; would also operate if ! use_Rcpp_cp and in other unanticipated cases
  if (inherits(x,"Matrix") || inherits(y,"Matrix")) { 
    if (is.null(y)) { # if ! use_Rcpp_cp and in other unanticipated cases
      if (TRUE && inherits(x,"dgCMatrix")) { ## FAST when (TRUE &&)
        # .dgCcrossprod+forceSymmetric manages to be faster than Matrix::crossprod() 
        txx <- .dgCcrossprod(x,y, as_mat=as_mat) 
        if (as_sym) {
          return(forceSymmetric(txx)) ## not as(.,"sparseMatrix") which checks -> slow
        } else return(txx)
      } else if (chk_sparse2mat && (maxd <- prod(dim(x)))>1L ) { 
        x_denseness <- .calc_denseness(x)/maxd 
        if (x_denseness>0.35) {
          ## Matrix::crossprod() is inefficient on math-dense matrices, even for dgeMatrix ! We try to maintain its return type, typically dsCMatrix:
          resu <- .crossprodCpp_d(as.matrix(x),y)
          if (as_sym) {
            return(forceSymmetric(resu)) 
          } else return(resu)
        } else return(.Matcrossprod(x, return_mat=as_mat)) 
      } else return(.Matcrossprod(x, return_mat=as_mat)) 
    } else if (FALSE && inherits(x,"dgCMatrix") &&  inherits(y,"dgCMatrix") ) { ## FASTer when (FALSE &&)
      txy <- .dgCcrossprod(x,y, as_mat=as_mat) 
      return(txy)
    } else if (chk_sparse2mat) {
      # .calc_denseness() is fast on Csparse but slow on dge
      if (inherits(x,"sparseMatrix") && (maxd <- prod(dim(x)))>1L) { 
        x_denseness <- .calc_denseness(x)/maxd 
        if (x_denseness>0.35) x <- as.matrix(x)
      } #else x_denseness <- 0 ## ie if the test is FALSE, no explicit conversion is attempted
      if (inherits(y,"sparseMatrix") && (maxd <- prod(dim(y)))>1L) { # idem
        y_denseness <- .calc_denseness(y)/maxd
        if (y_denseness>0.35) y <- as.matrix(y)
      }
      if (is.matrix(x) && is.matrix(y)) { # sparsity overhead is ~ 1/0.35
        return(.crossprodCpp_d(x, y)) # .crossprodCpp_d() slightly better than crossprod()
      } else return(.Matcrossprod(x, y, return_mat=as_mat)) 
    } else return(.Matcrossprod(x, y, return_mat=as_mat)) 
  } else { # *m*atrix case
    if (is.integer(x)) storage.mode(x) <- "double"
    if (is.integer(y)) storage.mode(y) <- "double"
    resu <- .crossprodCpp_d(x,y) ## faster than base::crossprod, but the matrix cannot be integer
    if (is.null(y)) {
      colnames(resu) <- rownames(resu) <- colnames(x)
    } else {
      rownames(resu) <- colnames(x)
      colnames(resu) <- colnames(y)
    }
    return(resu)
  }
}

.is_spprec_fit <- function(object) {
  if (is.null(MME_method <- object$how$MME_method)) MME_method <- object$MME_method
  length(intersect(c("AUGI0_ZX_sparsePrecision","AUGI0_ZX_spprec"), MME_method))>0L
}

# Post-fit fn
.get_tcrossfac_beta_w_cov <- function(fit) { 
  if (is.null(fit$envir$tcrossfac_beta_w_cov)) { ## tcrossfac...w stored under $envir while tcrossfac...v is under $beta_v_cov
    if (.is_spprec_fit(fit)) {
      tcrossfac_beta_v_cov <- .get_tcrossfac_beta_v_cov(fit$X.pv, fit$envir) ## typically permuted from triangular
      pforpv <- ncol(fit$X.pv)
      lhs <- fit$envir$chol_Q
      if (pforpv) {
        lhs_Xpv <- new("dtCMatrix",i= 0:(pforpv-1L), p=0:(pforpv), Dim=c(pforpv,pforpv),x=rep(1,pforpv))
        lhs <- Matrix::bdiag(lhs_Xpv, lhs )
      } 
      fit$envir$tcrossfac_beta_w_cov <- drop0(lhs %*% tcrossfac_beta_v_cov, tol=1e-16) 
    } else {
      if (is.matrix(beta_cov_info <- fit$envir$beta_cov_info) || ## matrix is old format, should be a list now
          is.null(tcrossfac_beta_w_cov <- beta_cov_info$tcrossfac_beta_v_cov)) {
        tcrossfac_beta_w_cov <- .get_beta_cov_info(fit)$tcrossfac_beta_v_cov 
      }    
      invL <- .get_invL_HLfit(fit) ## correlation matrix of ranefs is solve((t(invL)%*%(invL)))
      # invL is currently a single matrix for allranefs. de facto a block matrix when several ranefs
      if ( ! is.null(invL) && ! .is_identity(invL)) {
        pforpv <- ncol(fit$X.pv)
        v.range <- pforpv+seq(ncol(invL))
        if (inherits(tcrossfac_beta_w_cov,"sparseMatrix")) { # cf "dgC_good" code; we keep the sparseness as determined there
          fit$envir$tcrossfac_beta_w_cov <- rbind(tcrossfac_beta_w_cov[-v.range,],
                                                  .crossprod(invL, tcrossfac_beta_w_cov[v.range,], allow_as_mat = FALSE, eval_dens=FALSE))
        } else {
          tcrossfac_beta_w_cov[v.range,] <- .crossprod(invL, tcrossfac_beta_w_cov[v.range,], as_mat=TRUE) # invL may be sparse
          fit$envir$tcrossfac_beta_w_cov <- tcrossfac_beta_w_cov
        }
      } else fit$envir$tcrossfac_beta_w_cov <- tcrossfac_beta_w_cov
    }
  }
  return(fit$envir$tcrossfac_beta_w_cov) ## repeated calls to predVar suffices to make it useful to save it.
}

.eval_as_mat_arg <- function(processed) {
  if (is.null(processed$as_matrix)) {
    processed$as_matrix <- (
      processed$HL[1L]=="SEM" || ## SEM code does not yet handle sparse as it uses a dense Sig matrix
        ( ! processed$is_spprec && processed$QRmethod!="sparse" )
    )
  }
  processed$as_matrix
}

.eval_as_mat_arg.HLfit <- function(object) { 
  (
    object$HL[1L]=="SEM" || ## SEM code does not yet handle sparse as it uses a dense Sig matrix
      ! identical(object$QRmethod,"sparse")
  )
}



.get_invColdoldList <- function(res,regul.threshold=1e-7, control) {
  ## the validity of this fn is tested by checking that Evar (in predVar) is null in the case where explicit newdata=ori data
  if (identical(is.nan(control$fix_predVar),TRUE)) { # documented!
    res$envir$invColdoldList <- NULL
    control$fix_predVar <- NULL
  }
  if (is.null(invColdoldList <- res$envir$invColdoldList)) { 
    ## returns a list of inv(Corr) from the LMatrix
    strucList <- res$strucList
    fix_predVar_NULL <- FALSE
    if ( ! is.null(strucList)) {
      cum_n_u_h <- attr(res$lambda,"cum_n_u_h")
      vec_n_u_h <- diff(attr(res$lambda,"cum_n_u_h")) 
      invColdoldList <- structure(lapply(vec_n_u_h,Diagonal), names=seq_along(vec_n_u_h))
      if (res$spaMM.version < "2.2.116") {
        ranefs <- attr(res$ZAlist,"ranefs") 
      } else ranefs <- attr(res$ZAlist,"exp_ranef_strings") 
      for (Lit in seq_len(length(strucList))) {
        lmatrix <- strucList[[Lit]]
        if ( ! is.null(lmatrix)) {
          ## end of mat_sqrt implies either type is cholL_LLt or we have decomp $design_u and $d 
          type <-  attr(lmatrix,"type")
          invCoo <- NULL
          if ( ! is.null(latentL_blob <- attr(lmatrix,"latentL_blob"))) { ## from .process_ranCoefs
            if ( ! is.null(kron_Y <- res$ranef_info$sub_corr_info$kron_Y_LMatrices[[Lit]])) {
              # for CORREL kron_Y_LMatrices[[]] is the result of mat_sqrt(); 
              # for SPPREC it is the result of solve(blob$Q_CHMfactor,system="Lt") and the next solve() is not maximally efficient
              ## but in both case this is a *t*crossfac so 
              kron_Y <- crossprod(.inv_Lmatrix(kron_Y, regul.threshold=regul.threshold) ) # errors here should be easily detectable by 'newdata' tests
            }
            if ( is.null(invC <- latentL_blob$gmp_compactprecmat)) {
              if (is.null(invC <- latentL_blob$compactprecmat)) {
                invC <- solve(latentL_blob$compactcovmat) 
              } else invC <- as.matrix(invC)
            } # If kron_Y is NULL invC should be bigq or numeric matrix but not Matrix because:
            invCoo <- .makelong(invC,longsize=ncol(lmatrix),kron_Y=kron_Y) # (no template arg && NULL kron_Y) => .makelong_bigq() or [.C_makelong() => invC must be numeric matrix]
          } else if (inherits(lmatrix,"dCHMsimpl")) { # # before any test on type, because dCHMsimpl has a @type slot
            invCoo <- tcrossprod(as(lmatrix,"sparseMatrix")) # assuming 'lmatrix' is an unpermuted CHM factor of the precision factor (as for attr(lmatrix,"Q_CHMfactor"))
          } else if (type == "from_AR1_specific_code")  {
            invCoo <- crossprod(solve(lmatrix)) # cost of a sparse triangular solve.
          } else if (type == "from_Q_CHMfactor")  {
            invCoo <- tcrossprod(as(attr(lmatrix,"Q_CHMfactor"),"sparseMatrix")) ## correct but requires the attribute => numerical issues in computing Q_CHMfactor
          } else if (type %in% c("cholL_LLt","t(Matrix::chol)"))  { # Rmatrix upper tri in both case => suitable for chol2inv
            Rmatrix <- t(lmatrix)
            if (attr(lmatrix,"need_gmp")) {
              # control$fix_predVar controls what do do when there is a near singularity identified by "need_gmp".
              # In some cases it would be OK to use gmp, but not for large matrices => the default is thus to use ginv().
              # There is further control whether to output a message or not 
              # In the context of external packages, the user can only modify those packages' options. So to allow user control, 
              # the control$fix_predVar values used by these packages should be taken from these packages' options. 
              # That would also avoid checking the call stack by .check_frames(); more generally, for this reason,
              # calling predict() with explicit control$fix_predVar is recommended in programming.
              # An effective example is given by the objectivefn() definition in Infusion::MSL()
              if (is.null(fix_predVar <- control$fix_predVar)) {
                predVar_exceptions <- .spaMM.data$options$fix_predVar
                # Either 
                # .spaMM.data$options$fix_predVar lists a calling function in its NA|TRUE|FALSE components => 
                #   this block will not be called for all predict() calls, only once for any problematic fit 
                #    so the impact .check_frames() on performance should be negligible. 
                # or 
                # .spaMM.data$options$fix_predVar does NOT list a calling function 
                #   (e.g., blackbox::sampleByResp was forgotten at some stage) => repeated warnings
                # as fix (invCoo <- ginv(crossprod(Rmatrix))) only local since save in envir occurs only if ( ! fix_predVar_NULL)
                if (.check_frames(which=predVar_exceptions$"NA")) {
                  fix_predVar <- NA # single message + ginv + saved in envir  = (FALSE+warning)
                } else if (.check_frames(which=predVar_exceptions$"TRUE")) {
                  fix_predVar <- TRUE # NO warning + gmp + saved in envir; gmp would be a problem for large matrices, and so TRUE cannot be a default 
                } else if (.check_frames(which=predVar_exceptions$"FALSE")) {
                  fix_predVar <- FALSE # NO warning + ginv + saved in envir  = (NA - warning)
                } # else fix_predVar remains NULL => # repeated warnings + ginv + NOT saved in envir  
              }
              if (is.null(fix_predVar)) { # when user can handle this
                warning("A nearly-singular correlation matrix was fitted; see help('fix_predVar') for handling this", immediate. = TRUE)
                fix_predVar_NULL <- TRUE
              } else if (is.na(fix_predVar)) { # When user can do nothing and gmp cannot be assumed: warning+fallback 'if (is.null(invCoo))' computation below 
                message("A nearly-singular correlation matrix was fitted; ginv() will be used.")
              } else if (fix_predVar) { 
                invCoo <- gmp::tcrossprod(gmp::solve.bigq(gmp::as.bigq(Rmatrix))) 
                #invCoo <- Rmpfr::mpfr(invCoo,128)
              } # (else other cases including fix_predVar=FALSE:) fallback computation of  invCoo below.
            } # (else ! need_gmp: fallback computation of  invCoo below.
          } else { ## Rcpp's symSVD, or R's eigen() => LDL (also possible bad use of R's svd, not normally used)
            condnum <- kappa(lmatrix,norm="1")
            if (condnum<1/regul.threshold) {
              decomp <- attr(lmatrix,attr(lmatrix,"type")) ## of corr matrix !
              invCoo <-  try(.ZWZt(decomp$u,1/decomp$d),silent=TRUE) ## the idea is to use ginv() if this fails, rather than 
              #                                                         to attempt regularization, given that condnum has been tested as not large.
              if (inherits(invCoo,"try-error")) invCoo <- NULL
            }
            if (is.null(invCoo)) Rmatrix <- .lmwithQR(t(lmatrix),yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled # no pivoting compared to qr.R(qr(t(lmatrix))) 
          }
          # The fallback computation:
          if (is.null(invCoo)){ 
            invCoo <- try(chol2inv(Rmatrix),silent=TRUE)
            if (inherits(invCoo,"try-error") || max(abs(range(invCoo)))> 1e12) {
              invCoo <- ginv(.crossprod(Rmatrix, as_mat=TRUE))
            }
          }
          invColdoldList[[Lit]] <- invCoo
        } 
      }
      if ( ! fix_predVar_NULL) res$envir$invColdoldList <- invColdoldList
    } else return(NULL)
  } 
  return(invColdoldList)
}

.calcD2hDv2 <- function(ZAL,w.resid,w.ranef) { ## wrapper for a "ZtWZ minus diag" computation
  ## Si Gamma(identity) avec lambda >1 et w.ranef approche de de -1e6, et si on dit phi <- 1e-06, w.resid = 1e6 
  #    d2hdv2 peut etre une diag matrix with zome 0 elements => logabsdet=log(0)
  if (inherits(ZAL,"ddiMatrix") && ZAL@diag=="U") {
    d2hdv2 <- Diagonal(x= - w.resid - w.ranef)
  } else if (inherits(ZAL,"Matrix")) {
    if (attr(w.resid,"unique")) {
      crossprodZAL <- .crossprod(ZAL)
      d2hdv2 <- - w.resid[1L] * crossprodZAL
    } else d2hdv2 <- - .ZtWZwrapper(ZAL,w.resid) 
    if (inherits(d2hdv2,"dsCMatrix")) {
      d2hdv2 <- .dsCsum(d2hdv2,.symDiagonal(x=-w.ranef)) 
    } else { # dsy
      nc <- ncol(d2hdv2)
      diagPos <- seq.int(1L,nc^2,nc+1L)
      d2hdv2[diagPos] <- d2hdv2[diagPos] - w.ranef 
    }
  } else if (attr(w.resid,"unique")) {
    crossprodZAL <- .crossprod(ZAL)
    d2hdv2 <- - w.resid[1L] * crossprodZAL
    nc <- ncol(d2hdv2)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    d2hdv2[diagPos] <- d2hdv2[diagPos] - w.ranef 
  } else d2hdv2 <- .Rcpp_d2hdv2(ZAL,w.resid,w.ranef) 
  return(d2hdv2)
}

.gmp_solve <- function(X, as_numeric=TRUE) { ## "garbage in, garbage out:" the main problem may be the imprecision in the input X
  gmp_X <- gmp::as.bigq(X)
  invX <- gmp::solve.bigq(gmp_X) ##  with this invX gmp::`%*%`(gmp_X, invX) is exactly Identity
  if (as_numeric) invX <- gmp::asNumeric(invX) ## unfortunately with this invX X %*% invX can widely differ from Identity
  rownames(invX) <- colnames(invX) <- colnames(X)
  return(invX)
}

# .force_solve seems OK for a small matrix such as logdispInfo. Otherwise, ginv should be tried.
.force_solve <- function(X, 
                         try_gmp=identical(spaMM.getOption("try_gmp"),TRUE) ## default is FALSE
                         ) {
  condnum <- kappa(X)
  if (condnum==Inf) {
    ## logdispInfo can be exactly singular ! cf a twolambda example
    invX <- ginv(X) ## works exactly in onelambda/twolambda test!
    warning(paste("The matrix looks exactly singular."))
  } else {
    ## Regularization (minimal for gmp case: we need to control the sign of eigenvalues at least)
    esys <- eigen(X, only.values = TRUE)
    evalues <- esys$values
    #if (abs(evalues[1L])*abs(evalues[length(evalues)])<0L) {} ## serious inaccuracy as the matrix should pos or neg-definite
    negpos <- sign( abs(evalues[1L])-abs(evalues[length(evalues)]) )
    if (try_gmp) {threshold <- 1e100} else threshold <- 1e14 ## try_gmp case not tested after change
    if (negpos>0) { # large posive eigenvalue: aim for positive-def matrix
      min_d <- evalues[1L]/threshold ## so that corrected condition number is at most the denominator
      diagcorr <- max(c(0,min_d-evalues)) # SINGLE SCALAR
    } else  { # large negative eigenvalue: aim for negative-def matrix
      min_d_is_neg <- evalues[length(evalues)]/threshold ## same logic on negative
      diagcorr <- min(c(0,min_d_is_neg-evalues)) # SINGLE SCALAR
    } 
    diag(X) <- diag(X) + diagcorr ## # all diag is corrected => added a constant diagonal matrix to compactcovmat
    #
    if (try_gmp) {
      #     if ("package:gmp" %in% search()) { ## would be OK if default use-gmp were TRUE. 
      if (requireNamespace("gmp",quietly=TRUE)) { ## a *single* non-default act by the user is required: set use_gmp=TRUE
        invX <- .gmp_solve(X)
      } else {
        warning(#"  If the 'gmp' package were loaded, the inaccuracy could be substantially reduced."))
          "  If the 'gmp' package were installed, the inaccuracy could be substantially reduced.")
        invX <- ginv(X)
      }
    } else invX <- solve(X) ## assuming that the regularizaion was sufficient
  }
  return(invX)
}

.old_calc_d2hdv2_info <- function(object, ZAL) {
  envir <- object$envir
  if (is.null(factor_inv_Md2hdv2 <- envir$factor_inv_Md2hdv2)) { ## old code not using promises
    if ( ! is.null(envir$G_CHMfactor)) {
      if (length(object$y) > ncol(envir$chol_Q)) { # We precompute inverse(d2hdv2) so that .calc_dvdlogphiMat_new() has ONE costly (r*r) %*% (r*n).
        if ( is.null(envir$invL_G.P)) { # allowing for old objects that did not have this element
          LL <- Matrix::solve(envir$G_CHMfactor, as(envir$G_CHMfactor,"pMatrix") %*% envir$chol_Q,system="L") ## dgCMatrix 
        } else  LL <- envir$invL_G.P %*% envir$chol_Q ## dgCMatrix 
        envir$factor_inv_Md2hdv2 <- Matrix::drop0(LL,tol=.Machine$double.eps)
        d2hdv2_info <- - .crossprod(envir$factor_inv_Md2hdv2) 
      } else d2hdv2_info <- envir # then .calc_dvdlogphiMat_new() has TWO  (r*r) %*% (r*n) (plus an efficient solve())
    } else {
      d2hdv2 <- .calcD2hDv2(ZAL,.get_H_w.resid(object),object$w.ranef) 
      if (inherits(d2hdv2,"sparseMatrix")) {
        d2hdv2_info <- suppressWarnings(try(Cholesky( - d2hdv2,LDL=FALSE,perm=TRUE ), silent=TRUE )) #  '-' ! 
        if (inherits(d2hdv2_info, "CHMfactor")) {
          #d2hdv2_info <- structure(d2hdv2_info, BLOB=list2env(list(), parent=environment(.solve_CHM)))
          rank <- ncol(d2hdv2)
        } else {
          # oldMDCopt <- options(Matrix.warnDeprecatedCoerce = 0) # qr(<dsC>) problem in Matrix v1.4.2
          d2hdv2_info <- qr(d2hdv2, tol=spaMM.getOption("qrTolerance")) # sparseQR
          # options(oldMDCopt)
          rank <- sum(abs(diag(qrR(d2hdv2_info,backPermute=FALSE)))>1e-7)
        }
      } else { # d2hdv2 was a matrix or a Matrix in dense format (dgeMatrix...)
        d2hdv2_info <- qr(d2hdv2, tol=spaMM.getOption("qrTolerance")) # qr 
        rank <- d2hdv2_info$rank
      }
      if (rank<ncol(d2hdv2)) {
        d2hdv2_info <- .force_solve(d2hdv2) # inverse(d2hdv2)
      } # else we keep the QR facto.
    }
  } else d2hdv2_info <- - .crossprodCpp_d(as.matrix(factor_inv_Md2hdv2), yy=NULL) # inverse(d2hdv2) # much faster than Matrix::crossprod(dsCMatrix...)
  #                       and still faster than any of the as_mat variants, eg .Rcpp_crossprod(factor_inv_Md2hdv2, BB = NULL, as_mat=TRUE) 
  #                       bc .crossprodCpp_d(dense) is fast relative to .Rcpp_crossprod(sparse), irrespective of as_mat
  return(d2hdv2_info)
}

.calc_d2hdv2_info <- function(object, ZAL) {
  if (TRUE) return(.old_calc_d2hdv2_info(object, ZAL)) # rediscovered ong afterwards... so ./.
  ## The code below first checks for a "factor_inv_Md2hdv2" promise. But 
  ##   there is no such promise in the current object bc HLfit_body() has run  
  ##   .init_promises_spprec(sXaug, non_X_ones=FALSE, nullify_X_ones =TRUE) 
  ##   and non_X_ones=FALSE implies that there is no "factor_inv_Md2hdv2" promise.
  #
  # ELSE
  envir <- object$envir
  if ( ! is.null(envir$G_CHMfactor)) { ## spprec code
    if (! .is_evaluated("factor_inv_Md2hdv2", envir)) {
      if (length(object$y) > ncol(envir$chol_Q)) { # We precompute inverse(d2hdv2) so that .calc_dvdlogphiMat_new() has ONE costly (r*r) %*% (r*n).
        d2hdv2_info <- - .crossprod(envir$factor_inv_Md2hdv2) 
      } else d2hdv2_info <- envir # then .calc_dvdlogphiMat_new() has TWO  (r*r) %*% (r*n) (plus an efficient solve())
    } else d2hdv2_info <- - .crossprodCpp_d(as.matrix(envir$factor_inv_Md2hdv2), yy=NULL)
  } else {
    if (is.null(factor_inv_Md2hdv2 <- envir$factor_inv_Md2hdv2)) { 
      d2hdv2 <- .calcD2hDv2(ZAL,object$w.resid,object$w.ranef) 
      if (inherits(d2hdv2,"sparseMatrix")) {
        d2hdv2_info <- suppressWarnings(try(Cholesky( - d2hdv2,LDL=FALSE,perm=TRUE ), silent=TRUE )) #  '-' ! 
        if (inherits(d2hdv2_info, "CHMfactor")) {
          #d2hdv2_info <- structure(d2hdv2_info, BLOB=list2env(list(), parent=environment(.solve_CHM)))
          rank <- ncol(d2hdv2)
        } else {
          d2hdv2_info <- qr(d2hdv2, tol=spaMM.getOption("qrTolerance")) # sparseQR
          rank <- sum(abs(diag(qrR(d2hdv2_info,backPermute=FALSE)))>1e-7)
        }
      } else { # d2hdv2 was a matrix or a Matrix in dense format (dgeMatrix...)
        d2hdv2_info <- qr(d2hdv2, tol=spaMM.getOption("qrTolerance")) # qr 
        rank <- d2hdv2_info$rank
      }
      if (rank<ncol(d2hdv2)) {
        d2hdv2_info <- .force_solve(d2hdv2) # inverse(d2hdv2)
      } # else we keep the QR facto.
    } else d2hdv2_info <- - .crossprodCpp_d(as.matrix(envir$factor_inv_Md2hdv2), yy=NULL) # inverse(d2hdv2) # much faster than Matrix::crossprod(dsCMatrix...)
    #                       and still faster than any of the as_mat variants, eg .Rcpp_crossprod(factor_inv_Md2hdv2, BB = NULL, as_mat=TRUE) 
    #                  
  }
  return(d2hdv2_info)
}

.get_glm_phi <- function(fitobject, mv_it=NULL) { 
  # gets always a single glm fit, not the list of the mv case
  if ( ! is.null(mv_it)) { # fitmv case. Then glm_phi may be (1) an element of phi.object[[mv_it]] (2) envir$glmS_phi[[mv_it]]
    glm_phi <- fitobject$phi.object[[mv_it]][["glm_phi"]] 
    if (is.null(glm_phi)) glm_phi <- fitobject$envir$glmS_phi[[mv_it]]
    if (is.null(glm_phi)) {
      vec_nobs <- fitobject$vec_nobs
      if (is.null(fitobject$envir$glmS_phi)) fitobject$envir$glmS_phi <- vector("list", length(vec_nobs)) 
      cum_nobs <- c(0L,cumsum(vec_nobs))
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      validrownames <- attr(fitobject$data,"validrownames")[[mv_it]]
      if (length(fitobject$phi.object[[mv_it]]$phi_outer)>1L) {
        nobs <- length(resp_range)
        glm_phi_args <- list(formula=.get_phiform(fitobject, mv_it), 
                             lev=rep(0,nobs), dev.res=rep(1,nobs), control=list(), # dummy values
                             data=fitobject$data[validrownames,,drop=FALSE], 
                             family= .get_phifam(fitobject, mv_it))
      } else glm_phi_args <- c(fitobject$phi.object[[mv_it]]$glm_phi_args, 
                        list(formula=.get_phiform(fitobject, mv_it), 
                             lev=fitobject$lev_phi[resp_range], data=fitobject$data[validrownames,,drop=FALSE], 
                             family= .get_phifam(fitobject, mv_it))
      )
      fitobject$envir$glmS_phi[[mv_it]] <- glm_phi <- do.call(".calc_dispGammaGLM", glm_phi_args)
    }
  } else { # univariate case.  Then glm_phi may be (1) an element of phi.object (2) envir$glm_phi
    glm_phi <- fitobject$phi.object[["glm_phi"]] 
    if (is.null(glm_phi)) glm_phi <- fitobject$envir$glm_phi
    if (is.null(glm_phi)) { 
      if (length(fitobject$phi.object$phi_outer)>1L) {
        nobs <- nrow(fitobject$X.pv)
        glm_phi_args <- list(formula=.get_phiform(fitobject), # formula with offset only: see explanations in .calcResidVar()
                             lev=rep(0,nobs), dev.res=rep(1,nobs), control=list(), # dummy values
                             data=fitobject$data, 
                             family= .get_phifam(fitobject))
      } else glm_phi_args <- c(fitobject$phi.object$glm_phi_args, 
                        list(formula=.get_phiform(fitobject),
                             lev=fitobject$lev_phi, data=fitobject$data, 
                             family= .get_phifam(fitobject))
      )
      fitobject$envir$glm_phi <- glm_phi <- do.call(".calc_dispGammaGLM", glm_phi_args)
    } 
  }
  return(glm_phi)
}

.calc_wAugX <- function(XZ_0I,sqrt.ww) {
  return(.Dvec_times_m_Matrix(sqrt.ww, XZ_0I)) # sweep(TT,MARGIN=1,sqrt.ww,`*`) # rWW %*%TT
}

.calc_XZ_0I <- function(AUGI0_ZX,ZAL) { ## single use to create old augX order from AUGI0_ZX 
  if (inherits(ZAL,"Matrix")) {
    XZ_0I <- suppressMessages(cbind2(
      rbind2(AUGI0_ZX$X.pv,AUGI0_ZX$ZeroBlock), # XZ_0I order
      rbind2(ZAL, AUGI0_ZX$I)
    )) 
  } else {
    XZ_0I <- cbind(
      rbind(AUGI0_ZX$X.pv,AUGI0_ZX$ZeroBlock), 
      rbind(ZAL, AUGI0_ZX$I)
    ) 
  }
  return(XZ_0I) ## aug design matrix 
}

# mostly post-fit
.calc_Md2hdvb2_info_spprec_by_QR <- function(envir, which="tcrossfac_v_beta_cov") {## sparse qr is fast
  qr_R <- qrR(envir$qrXa, backPermute=FALSE) # the function was called on the condition that envir$qrXa had been evaluated 
  X_scale <- envir$X_scale # pre v.2.7.34
  if (is.null(X_scale)) X_scale <- attr(envir$sXaug$AUGI0_ZX$X.pv,"scaled:scale") # from version 3.2.19
  if (is.null(X_scale)) X_scale <- attr(envir$sXaug,"scaled:scale") # from v2.7.34 to 3.2.18 
  if ( ! is.null(X_scale)) { ## X_scale impacts post-fit computations. 
    Xnames <- names(X_scale) 
    qr_R[,Xnames] <- .m_Matrix_times_Dvec(qr_R[,Xnames,drop=FALSE], X_scale)
  }
  if (which=="tcrossfac_v_beta_cov") { ## 
    R_invMd2hdvb2 <- solve(qr_R) ## such that tcrossprod(R_invMd2hdvb2) (indep of back-permutation)= solve(Md2hdvb2) = solve(tcrossprod(Xscal))
    return(R_invMd2hdvb2) ## not necess. triangular
    # For other uses than a tcross factor, we may need:
    # qI <- sort.list(envir$qrXa@q)
    # return(R_invMd2hdvb2[,qI]) ## non-long triangular
  } else if (which=="R_Md2hdbv2") {
    stop("need to account for permutation")
    return(qr_R) ## such that tcrossprod(R_invMd2hdvb2) = solve(Md2hdbv2) = solve(crossprod( this ))
  } else { ##
    v_beta_cov <- Matrix::chol2inv(qr_R) ## true indep of back-permutation
    return(v_beta_cov) ## v_beta_cov
  }
}

# Extractor mostly post fit, although using compatibility with spprec structure, with envir=sXaug$BLOB,
# to use it in an spprec function called "mostly post fit" but possibly in-fit.
.get_H_w.resid <- function(object, envir=object$envir) {
  ## this is called in contexts where the fit object is not accessible, so 
  H_w.resid <- envir$H_w.resid # new HLfit objects , with structured "w.resid" attribute
  if (is.null(H_w.resid)) {
    H_w.resid <- attr(envir$sXaug,"w.resid") # old MM
    if (is.null(H_w.resid)) H_w.resid <- object$w.resid  # old GLM and SEM
  }
  if ( is.null(attr(H_w.resid,"w.resid")) && # old object
       ! missing(object))  {
    attr(H_w.resid,"w.resid") <- object$w.resid # back compat, so that there is always something 
    # for old objects this is the vector so may not be appropriate in cases where old spaMM failed (e.g., hatvalues from old mv fit)
  }
  H_w.resid
}

# mostly post-fit but not only
.calc_Md2hdvb2_info_spprec_by_r22 <- function(X.pv, # currently the unscaled one (the object's $X.pv)
                                              envir, which="tcrossfac_v_beta_cov") {
  # called e.g. for "tcrossfac_v_beta_cov" for which r12,r22 needed
  # same if it is called for "R_Md2hdbv2" (as through .get_tcrossfac_beta_w_cov())  
  # if we used it for v_beta_cov, could we simplify the code by only computing crossprod_r22 using $invG_ZtWX ? What about r_12 ?
  
  if (FALSE) { # if ( ! object$spaMM.version < "3.7.10")  except that "object" not easily accessible
    #if (is.null(tcrossfac_v_beta_cov <- envir$tcrossfac_v_beta_cov)) # No strong reason to save it. The only re-use is for cAIC_pd
    # scaled_X <- envir$sXaug$AUGI0_ZX$X.pv # !! use the scaled X.pv to compute quantities consistent with other promises !!
    # # !!! but then R_invMd2hdvb2 etc will refer to a scaled X !!! 
    {
      if (which=="R_Md2hdbv2") { ## not called in the source, and the single occurrence of solve(factor_inv_Md2hdv2)
        if (ncol(X.pv)) {
          r22 <- as(envir$r22,"sparseMatrix") # not sure why
          R_Md2hdbv2 <- rbind(cbind(t(solve(envir$factor_inv_Md2hdv2)), envir$r12), # cbind(chol_Md2hdv2,r12),
                              cbind(Matrix(0,nrow=ncol(X.pv),ncol=ncol(envir$chol_Q)), r22)) ## R_a in the working doc
        } else R_Md2hdbv2 <- t(solve(envir$factor_inv_Md2hdv2))
        return(R_Md2hdbv2) ## such that tcrossprod(R_invMd2hdvb2) = solve(Md2hdbv2) = solve(crossprod(R_Md2hdbv2))
      } else {
        if (ncol(X.pv)) {
          inv11 <- t(envir$factor_inv_Md2hdv2)# solve(chol_Md2hdv2)
          inv22 <- solve(envir$r22)
          inv12 <- -1* inv11 %*% (envir$r12 %*% inv22) ## does not seem to understand unary '-' (late PS: in old version of Matrix?)
          R_invMd2hdvb2 <- rbind(cbind(inv11,inv12),
                                 cbind(Matrix(0,nrow=ncol(X.pv),ncol=ncol(envir$chol_Q)), inv22))
        } else R_invMd2hdvb2 <- t(envir$factor_inv_Md2hdv2)# only the inv11 block, solve(chol_Md2hdv2)
        if (.spaMM.data$options$perm_G) {
          tcrossfac_v_beta_cov <- as(R_invMd2hdvb2,"sparseMatrix") # not necess. triangular
        } else tcrossfac_v_beta_cov <- as(R_invMd2hdvb2,"triangularMatrix") ## not sure that being triangular would be of any use.
      }
    }
  } else {
    if (is.null(factor_inv_Md2hdv2 <- envir$factor_inv_Md2hdv2)) {
      if ( is.null(envir$invL_G.P)) { # allowing for old objects that did not have this element
        LL <- Matrix::solve(envir$G_CHMfactor, as(envir$G_CHMfactor,"pMatrix") %*% envir$chol_Q,system="L") ## dgCMatrix 
      } else  LL <- envir$invL_G.P %*% envir$chol_Q ## dgCMatrix 
      envir$factor_inv_Md2hdv2 <- factor_inv_Md2hdv2 <- Matrix::drop0(LL,tol=.Machine$double.eps)
      # such that chol_Md2hdv2 = t(solve(factor_inv_Md2hdv2))
    }
    if (is.null(r22 <- envir$r22)){ # => post-fit, presumably
      if (ncol(X.pv)) {
        # 'object' is not accessible here so its $spaMM.version isn't. This comes at the end of a long stack of calls...
        H_w.resid <- .get_H_w.resid(envir=envir)
        if (is.null(envir$ZtW)) envir$ZtW <- t(.Dvec_times_m_Matrix(H_w.resid, envir$sXaug$AUGI0_ZX$ZAfix))
        r12 <- as(Matrix::solve(envir$G_CHMfactor, as(envir$G_CHMfactor,"pMatrix") %*% envir$ZtW %*% X.pv,system="L"),"sparseMatrix") ## solve(as(envir$G_CHMfactor,"sparseMatrix"), envir$ZtW %*% AUGI0_ZX$X.pv)
        r22 <- .calc_r22_postfit(X.pv,H_w.resid,r12, XtWX=envir$XtWX) ## both lines as explained in working doc
        r22 <- as(r22,"sparseMatrix") 
      } else r22 <- diag(nrow=0L) # don't leave it NULL; + we use its ncol
    } else r12 <- envir$r12 
    if (which=="R_Md2hdbv2") { ## not called in the source, and the single occurrence of solve(factor_inv_Md2hdv2)
      if (ncol(r22)) {
        R_Md2hdbv2 <- rbind(cbind(t(solve(factor_inv_Md2hdv2)), r12), # cbind(chol_Md2hdv2,r12),
                            cbind(Matrix(0,nrow=ncol(X.pv),ncol=ncol(envir$chol_Q)), r22)) ## R_a in the working doc
      } else R_Md2hdbv2 <- t(solve(factor_inv_Md2hdv2))
      return(R_Md2hdbv2) ## such that tcrossprod(R_invMd2hdvb2) = solve(Md2hdbv2) = solve(crossprod(R_Md2hdbv2))
    } else {
      if (ncol(r22)) {
        inv11 <- t(factor_inv_Md2hdv2)# solve(chol_Md2hdv2)
        inv22 <- solve(r22, b= .sparseDiagonal(n=ncol(r22), shape="t"))
        inv12 <- -1* inv11 %*% (r12 %*% inv22) ## does not seem to understand unary '-'
        R_invMd2hdvb2 <- rbind(cbind(inv11,inv12),
                               cbind(Matrix(0,nrow=ncol(X.pv),ncol=ncol(envir$chol_Q)), inv22))
      } else R_invMd2hdvb2 <- t(factor_inv_Md2hdv2)# only the inv11 block, solve(chol_Md2hdv2)
      if (.spaMM.data$options$perm_G) {
        tcrossfac_v_beta_cov <- as(R_invMd2hdvb2,"sparseMatrix") # not necess. triangular
      } else tcrossfac_v_beta_cov <- as(R_invMd2hdvb2,"triangularMatrix") ## not sure that being triangular would be of any use.
    } 
  }
  if (which=="tcrossfac_v_beta_cov") {
    return(tcrossfac_v_beta_cov)  
  } else { 
    v_beta_cov <- .tcrossprod(tcrossfac_v_beta_cov) ## it happens that the beta,beta block is solve(crossprod(r22)), a property used in .calc_inv_beta_cov() but not here.
    return(v_beta_cov)
  }
}



.calc_Md2hdvb2_info_spprec <- function(X.pv,envir, which="tcrossfac_v_beta_cov") { 
  # called by post-fit functions, but else by .get_tcrossfac_beta_v_cov which *may* be called in-fit 
  # Hence from spprec there *may* be an envir$qrXa promise. 
  if ( .is_evaluated("qrXa",envir)) { # If we have a QR factorization of augmented X, we can use it efficiently
    # This is tested for example by get_predVar(fitsparse,newdata=newXandZ) on 'landsat' data
    .calc_Md2hdvb2_info_spprec_by_QR(envir, which)
  } else {
    .calc_Md2hdvb2_info_spprec_by_r22(X.pv, envir, which)
  }
}

# called by .get_tcrossfac_beta_w_cov() post-fit; or by .calc_beta_cov_info_spprec(), in-fit or post-fit:
# ie this *may* be called during a fit, by  .calc_beta_cov_info_spprec(... envir = sXaug$BLOB ...) -> .calc_beta_cov_info_spprec -> .get_tcrossfac_beta_v_cov
#    in the diagnostic code for poor breakcond.
.get_tcrossfac_beta_v_cov <- function(X.pv, envir) {
  if (is.null(envir$beta_cov_info$tcrossfac_beta_v_cov)) {
    ## Using sparse matrices as much as possible for computation:
    tcrossfac_v_beta_cov <- .calc_Md2hdvb2_info_spprec(X.pv, envir, which="tcrossfac_v_beta_cov")
    # but the v_beta_cov matrix is (mathematically) dense whatever its format
    pforpv <- ncol(X.pv)
    n_u_h <- ncol(envir$chol_Q)
    seqp <- seq_len(pforpv)
    perm <- c(n_u_h+seqp,seq_len(n_u_h))
    envir$beta_cov_info$tcrossfac_beta_v_cov <- tcrossfac_v_beta_cov[perm,,drop=FALSE] # useful to keep it for predVar computations?
  }
  return(envir$beta_cov_info$tcrossfac_beta_v_cov) ## Used for beta_v_cov itself, then for (beta_w_cov for predVar). Seem useful to  keep it
}

# may be called in-fit or post-fit
.calc_beta_cov_info_spprec <- function(X.pv, envir) { 
  tcrossfac_beta_v_cov <- .get_tcrossfac_beta_v_cov(X.pv, envir)
  if (is.null(beta_cov <- envir$beta_cov_info$beta_cov)) { # if stripHLfit()ed
#    beta_v_cov <- .tcrossprod(tcrossfac_beta_v_cov) # if the L differs among method the meaning of v differs too, and beta_v_cov too  
    # ./. but we can check that is is = solve(crossprod(wAugX)) by reconstructing wAugX with the matching L as in .get_LSmatrix()
    pforpv <- ncol(X.pv)
    seqp <- seq_len(pforpv)
#    beta_cov <- beta_v_cov[seqp,seqp,drop=FALSE]
    beta_cov <- .tcrossprod(tcrossfac_beta_v_cov[seqp,,drop=FALSE]) # if the L differs among method the meaning of v differs too, and beta_v_cov too  
    colnames(beta_cov) <- rownames(beta_cov) <- colnames(X.pv)
    envir$beta_cov_info$beta_cov <- as.matrix(beta_cov) # make_beta_table() expects a *m*atrix
  }
  return(envir$beta_cov_info)
}

## returns a list !!
## input XMatrix is either a single LMatrix which is assumed to be the spatial one, or a list of matrices 
.compute_ZAXlist <- function(XMatrix, ZAlist, force_bindable=FALSE) {
  ## ZAL is nobs * (# levels ranef) and ZA too
  ## XMatrix is (# levels ranef) * (# levels ranef) [! or more generally a list of matrices!]
  ## the levels of the ranef must match each other in multiplied matrices
  ## the only way to check this is to have the levels as rownames and colnames and to check these
  if (is.null(ZAlist)) return(list())
  ## ELSE
  ZAX <- ZAlist
  bindable <- TRUE
  if ( ! is.null(XMatrix) && length(ZAlist) ) {
    if ( ! inherits(XMatrix,"list")) XMatrix <- list(dummyid=XMatrix)
    LMlen <- length(XMatrix)
    for (Lit in seq_len(LMlen)) {
      xmatrix <- XMatrix[[Lit]]
      if (inherits(xmatrix,"dCHMsimpl")) { ## occurs more or less depending on use_ZA_L
        ## if use_ZA_L is TRUE, in spprec and with the exception of AR1 ; 
        ## or if use_ZA_L is NULL, with further condition of using augZXy (e.g. test-spde-perm_Q.R: IMRF LMM):
        ## In .assign_geoinfo_and_LMatrices_but_ranCoefs():
        ## if (
        ##   (
        ##      identical(use_ZA_L,TRUE) ||
        ##      (is.null(use_ZA_L) && processed$augZXy_cond) # default is to use_ZA_L for augZXy_cond
        ##   )
        ##   && ! processed$AUGI0_ZX$vec_normIMRF[rd]) { # solve may be OK when IMRF is normalized (otherwise more code needed in .normalize_IMRF) 
        ##   Lunique <- Q_CHMfactor # representation of Lunique, not Lunique itself
        ##   # note that this has an attribute "type" !  
        ## }
        ## indeed runs Lunique <- Q_CHMfactor. Then .compute_ZAL(), which calls this fn, by default 
        ##    returns an object of S4 class 'ZAXlist'. 
        if (force_bindable) {
          xmatrix <- as(xmatrix,"pMatrix") %*% solve(xmatrix, system="Lt", b=.sparseDiagonal(n=ncol(xmatrix), shape="g"))
        } else bindable <- FALSE
      } else if (inherits(xmatrix,"bigq") ) {
        xmatrix <- .mMatrix_bigq(xmatrix) # this does not seem to be used for the Evar computation so we may drop precision
      } 
      if ( is.null(xmatrix)) {
        # do nothing, ZAX[[Lit]] remains equal to ZA[[Lit]]
      } else {
        ZA <- ZAlist[[Lit]]
        if (inherits(xmatrix,"Kronfacto")) {
          if (force_bindable) {
            xmatrix <- xmatrix@BLOB$long 
          } else if (.is_identity(ZA)) {
            ZAX[[Lit]] <- xmatrix 
          } else ZAX[[Lit]] <- structure(list(ZA=ZA, Kronfacto=xmatrix), class=c("ZA_Kron","list")) 
        } else if (.is_identity(ZA)) {
          if (inherits(xmatrix,"dCHMsimpl")) {
            ZAX[[Lit]] <- structure(list(ZA=ZA, Q_CHMfactor=xmatrix), class=c("ZA_QCHM","list")) 
          } else ZAX[[Lit]] <- xmatrix ## matrix may replace Matrix here...          
        } else { ##
          # Including the case of composite nested ranefs i.e  <keyword>(RHS|nested LHS) 
          # composite => there is an x matrix, we reach here
          # nested LHS => not much difference from non nested one ?
          # RHS => in provide G diagnosis, xmatrix is a kronecker product locnr>locnc ; some care is needed there beforethe .compute_ZAXlist() call.  
          locnc <- ncol(ZA)
          locnr <- nrow(xmatrix)
          if ( locnc %% locnr !=0L) {
            # then names are needed to match the row and columns
            if (length(setdiff(colnames(ZA),rownames(xmatrix)))) { # some names missing from xmatrix 
              mess <- "Cannot match the correlation matrix:\n"
              mess <- paste0(mess, "The number of levels of the grouping variable in random term ", 
                             attr(ZAlist,"exp_ranef_strings")[Lit])
              mess <- paste0(mess,"\n  is not a multiple of the dimension of the correlation matrix;") ## by distMatrix checking in corrHLfit or no.info check somewhere...
              mess <- paste0(mess,"\n  and levels of the grouping variable cannot all be matched by name to dimnames of the correlation matrix.")
              stop(mess) # for spprec,  this may mean something wrong occurred in .init_assign_geoinfo() when cols where added to Z (eg non-unique names in input)
            } else {
              xmatrix <- xmatrix[colnames(ZA),colnames(ZA),drop=FALSE]
              locnr <- locnc
            }
          } 
          #
          # removed (=>3.6.1) here code beginning with nblocks <- locnc %/% locnr, maybe for an obsolete version of nested AR1. 
          #
          ## With a proxy::dist or 'crossdist' matrix, it is likely that ZA was = I and we don't reach this code;
          ## However, exceptions can occur: cf Infusion with CIpoint = MLE => replicate in points where MSEs are to be estimated
          ## Then the xmatrix must have been converted from proxy style to a matrix.
          if (inherits(xmatrix,"dCHMsimpl")) {
            zax <- structure(list(ZA=ZA, Q_CHMfactor=xmatrix), class=c("ZA_QCHM","list")) 
            # F I X M E not sure about colnames handling here. Can this occur for the xmatrix of a random-coef?
          } else if (inherits(ZA,"dgCMatrix") &&  inherits(xmatrix,"dgCMatrix") ) {
            # In random slope models, xmatrix can be a Matrix
            #   matmult <- getMethod(`%*%`, c("CsparseMatrix", "CsparseMatrix"))
            #   zax <- matmult(ZA, xmatrix)
            zax <- .dgCprod(ZA, xmatrix) 
            #zax@Dimnames[2] <- xmatrix@Dimnames[2] # now in cpp code
            ## : colnames needed for predict(., newdata) -> match_old_new_levels(), and maybe elsewhere. 
          } else {
            zax <- ZA %*% xmatrix # measurably slow... 
            #if (identical(attr(xmatrix,"corr.model"),"random-coef")) colnames(zax) <- colnames(ZA) ## presumably obsolete. xmatrix must provide colnames
          }
          ZAX[[Lit]] <- zax
        }
      }
    }
  }
  if ( ! bindable) class(ZAX) <- c("list","notBindable") ## keep this order otherwise the ZAXlist constructor does not .../...
  #                                                   recognize a list since "notBindable" is not declared as extending "list"
  return(ZAX)
}

## replace missing (as(<ddi>, "dGCMatrix)) method
# also handles missing 'Z' but explicit 'ncol_Z' and 'x'
.as_ddi_dgC <- function(Z=NULL, ncol_Z=ncol(Z), x=Z@x) { 
  if ( ! length(x)) x <- rep(1,ncol_Z)
  Z <- new("dgCMatrix",i= 0:(ncol_Z-1L), p=c(0L:(ncol_Z)), Dim=c(ncol_Z,ncol_Z),x=x)
  return(Z)
}


# even though the Z's were sparse postmultplication by LMatrix leads some of the ZAL's to dgeMatrix (dense) 
#                                 [and replacement by LMatrix may give a *m*atrix !]
.ad_hoc_cbind <- function(ZALlist, as_matrix ) {
  nrand <- length(ZALlist)
  if ( as_matrix ) {
    for (rd in seq_len(nrand)) {
      if (inherits(ZALlist[[rd]],"Kronfacto")) ZALlist[[rd]] <- ZALlist[[rd]]@BLOB$long
      ZALlist[[rd]] <- as.matrix(ZALlist[[rd]]) 
    }
    if (nrand>1L) {ZAL <- do.call(cbind,ZALlist)} else ZAL <- ZALlist[[1L]]
  } else {
    for (rd in seq_len(nrand)) {
      if (inherits(ZALlist[[rd]],"Kronfacto")) ZALlist[[rd]] <- ZALlist[[rd]]@BLOB$long
      if ( # is.matrix(ZALlist[[rd]]) || ## seems to work but at a cost for speed.
        inherits(ZALlist[[rd]],"dgeMatrix")) ZALlist[[rd]] <- as(ZALlist[[rd]],"CsparseMatrix")
    }
    ## but leave diagonal matrix types unchanged 
    if (nrand>1L) { ## ZAL <- suppressMessages(do.call(cbind,ZALlist))
      ZAL <- ZALlist[[1L]]
      for (rd in 2L:nrand) {
        if (inherits(ZAL,"dgCMatrix") &&  inherits(ZALlist[[rd]],"dgCMatrix") ) {
          ZAL <- .cbind_dgC_dgC(ZAL, ZALlist[[rd]]) 
        } else ZAL <- cbind(ZAL,ZALlist[[rd]])
      }
      
    } else ZAL <- ZALlist[[1L]]
  } 
  return(ZAL)
}

.compute_ZAL <- function(XMatrix, ZAlist, as_matrix, bind.=TRUE, force_bindable=bind.) { # ideally force_bindable should be ( ! processed$is_spprec)
  ZALlist <- .compute_ZAXlist(XMatrix, ZAlist, force_bindable=force_bindable) # force_bindable=TRUE to avoid Kronfacto in result 
  if ( bind. && ! inherits(ZALlist,"notBindable")) {
    ZAL <- .ad_hoc_cbind(ZALlist, as_matrix )
    return(ZAL)
  } else return(new("ZAXlist", LIST=ZALlist)) # __F I X M E__ could add warning if bind. was TRUE
}


if (FALSE) { # that's not used.
  ## cette fonction marche que si on a fixed effect + un terme aleatoire....
  .eval_corrEst_args <- function(family,rand.families,predictor,data,X.Re,
                                 REMLformula,ranFix,
                                 term=NULL,
                                 Optimizer) {
    ## ici on veut une procedure iterative sur les params de covariance
    #  HLCor.args$processed <- processed ## FR->FR dangerous in early development
    corrEst.args <- list(family=family,rand.family=rand.families) ## but rand.families must only involve a single spatial effect 
    loc.lhs <- paste(predictor)[[2]]
    ## build formula, by default with only spatial effects
    if (is.null(term)) term <- .findSpatial(predictor)
    corrEst.form <-  as.formula(paste(loc.lhs," ~ ",paste(term)))
    corrEst.args$data <- data ## FR->FR way to use preprocess ???                    
    # if standard ML: there is an REMLformula ~ 0; ____processed$X.Re is 0-col matrix, not NULL____
    # if standard REML: REMLformula is NULL: processed$X.Re is NULL
    # non standard REML: other REMLformula: processed$X.Re may take essentially any value
    if (is.null(X.Re) ) { ## processed$X.Re should be NULL => actual X.Re=X.pv => standard REML 
      corrEst.args$REMLformula <- predictor ## standard REML 
    } else corrEst.args$REMLformula <- REMLformula ## _ML_ _or_ non-standard REML
    if (NCOL(X.Re)) { ## some REML correction (ie not ML)
      corrEst.args$objective <- "p_bv" ## standard or non-standard REML
    } else corrEst.args$objective <- "p_v" ## ML
    corrEst.args$ranFix <- ranFix ## maybe not very useful
    corrEst.args$control.corrHLfit$optimizer<- Optimizer ## (may be NULL) 
    corrEst.args$control.corrHLfit$optim$control$maxit <- 1 
    corrEst.args$control.corrHLfit$optimize$tol <- 1e10 
    return(list(corrEst.args=corrEst.args,corrEst.form=corrEst.form))
  }
  
}

.corr_notEQL_lambda <- function(nrand,cum_n_u_h,lambda_est,lcrandfamfam) {  
  ## d h/ d !log! lambda correction (nul for gaussian ranef)
  ## ! correction for not using the deviance residuals as approx for the distribution of the random effects. It's not specifically ReML !
  ## this is a trick for still using deviances residuals in the Gamma GLM
  notEQL <- vector("list", nrand)
  for (it in seq_len(nrand)) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    loclambda <- lambda_est[u.range]
    notEQL[[it]] <- switch(lcrandfamfam[it], 
                   gaussian=rep(0,length(u.range)),
                   gamma=1+2*(log(loclambda)+digamma(1/loclambda))/loclambda,## cf notes on p. 89 of the book
                   "inverse.gamma"=1+2*(log(loclambda)-loclambda+digamma(1+(1/loclambda)) )/loclambda, ## appears to be the same as for the gamma case [digamma(1+x)=digamma(x)+1/x]... 
                   beta=1-2*(digamma(1/loclambda)/loclambda)+2*(digamma(1/(2*loclambda))/loclambda)+log(4)/loclambda
    ) ## consistent with HGLMMM
  }
  return(.unlist(notEQL))
}

.initialize_v_h <- function(processed, etaFix, init.HLfit) {
  v_h <- etaFix$v_h
  if (is.null(v_h) ) v_h <- processed$port_env$port_fit_values$v_h
  if (is.null(v_h) ) v_h <- init.HLfit$v_h
  if (is.null(v_h) ) {
    psi_M <- processed$psi_M
    cum_n_u_h <- processed$cum_n_u_h
    rand.families <- processed$rand.families
    v_h_list <- processed$reserve$rand_list
    for (it in seq_len(length(rand.families))) {
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      v_h_list[[it]] <- rand.families[[it]]$linkfun(psi_M[u.range]) ## v as link(mean(u)) 
    }
    v_h <- .unlist(v_h_list)
  }
  v_h
}

.cosy_names_lambdas <- function(namesnames) { # attempt to standardize previously messy things
  # these are names that users can see from lambda.object and were used in examples or programming
  for (it in seq_len(length(namesnames))) if (nchar(namesnames[it])>10) namesnames[it] <- paste0(substr(namesnames[it],0,9),".")
  make.unique(namesnames,sep=".")
}

.make_lambda_object <- function(nrand, lambda_models, cum_n_u_h, lambda_est,                        
                               process_resglm_blob, rand.families, ZAlist, next_LMatrices, lambdaType) {
  ## redefine names(namesTerms) and namesTerms elements
  print_namesTerms <- attr(ZAlist,"namesTerms") ## a list, which names correspond to the grouping variable, and elements are the names of the coefficients fitted
  names(print_namesTerms) <- .cosy_names_lambdas(names(print_namesTerms)) ## makes group identifiers unique (names of coeffs are unchanged); using base::make.unique
  # these names are used e.g. in VarCorr. There a ranCoefs has an Intercept and a <predictor>  coefficient, 
  #  and we create a conflict if "(Intercept)" is replaced by <predictor> (ie, by a namesname)
  # OTOH
  # now print_namesTerms is a list with unique names but possibly with repeated elements "(Intercept)" which is bad as only the elements will be used
  # => we keep here an informative lambda_list with redundant names and extract elsewhere a cosy $lambda by an ad hoc function .print_lambda().
  lambda_list <- vector("list",nrand) 
  Xi_cols <- attr(ZAlist,"Xi_cols")
  for (rd in seq_len(nrand)) {
    plam_rd <- process_resglm_blob$lambda_pred_list[[rd]] # from .calc_dispGammaGLM our CAR dispGammaGLM
    if (anyNA(plam_rd)) { ## those for which no glm was available, such as fixed lambdas...
      Xi_ncol <- Xi_cols[rd]
      plam_rd <- numeric(Xi_ncol)
      u.range <- (cum_n_u_h[rd]+1L):(cum_n_u_h[rd+1L])
      blocksize <- length(u.range)/Xi_ncol
      for (colit in seq_len(Xi_cols[rd])) {
        u.subrange <- u.range[((colit-1L)*blocksize+1L):(colit*blocksize)]
        plam_rd[colit] <- unique(lambda_est[u.subrange])
      } # With "chol" transfo pour ranCoefs, all the lambda's are 1, so subrange'ing is essential # all 1 ? not sure for inner-estimated ones.
    }
    # lambdaType does not distinguish ranCoefs
    if (lambdaType[rd]!="inner" && ## bc for inner, plam_rd is $lambda_pred_list[[rd]] and is already the single lambda param
        ! is.null(prior_lam_fac <- rand.families[[rd]]$prior_lam_fac)) {
      plam_rd <- plam_rd[1L]/prior_lam_fac[1L]
    } ## F I X M E It would be better to use a lambda_par attribute for lambda_est
    lambda_list[[rd]] <- structure(plam_rd, names=print_namesTerms[[rd]])
  }
  names(lambda_list) <- names(print_namesTerms)
  attr(lambda_list,"cum_n_u_h") <- cum_n_u_h
  lambda.object <- list(lambda_est = lambda_est,  ## full vector for simulate() calc_logdisp_cov()
                        lambda_list=lambda_list, ## nrand-elements list, multiple uses beyond unlisting it into hlfit's $lambda. 
                        type = lambdaType)  
  if (any(lambdaType=="inner")) { ## modifies default namesTerms
    coefficients_lambdaS <- process_resglm_blob$coefficients_lambdaS
    for (it in seq_len(length(coefficients_lambdaS))) { ## detect exceptions
      if (lambdaType[it]=="inner") {
        coefficients <- names(coefficients_lambdaS[[it]]) 
        if ("adjd" %in% coefficients) print_namesTerms[[it]] <- coefficients
      }
    }
    lambda.object <- c(lambda.object,
                       list(coefficients_lambdaS=coefficients_lambdaS,  
                            rand_to_glm_map=process_resglm_blob$rand_to_glm_map,
                            lambda_se=unlist(process_resglm_blob$lambda_seS),
                            linkS = process_resglm_blob$linkS,
                            linkinvS = process_resglm_blob$linkinvS ) 
    )
    attr(lambda.object,"warning") <- unlist(process_resglm_blob$warnmesses) ## may be NULL
  } 
  lambda.object$print_namesTerms <-  print_namesTerms 
  return(lambda.object)
}

.timerraw <- function(time1) {
  return(round(as.numeric(difftime(Sys.time(), time1, units = "secs")), 2))
}

.post_process_v_h_LMatrices <- function(next_LMatrices, v_h, u_h, processed, ZAlist=processed$ZAlist, ranCoefs_blob, 
                                    cum_n_u_h=processed$cum_n_u_h) { 
  strucList <- vector("list", length(ZAlist))
  for (rd in seq_len(length(ZAlist))) {
    lmatrix <- next_LMatrices[[rd]]
    ## keep in mind that str(S4...) does not show extra attributes
    ## lmatrix is an S4 object...
    ## ranefs itself has attribute type="(.|.)"
    if ( ! is.null(lmatrix)) {
      corr.model <- attr(lmatrix, "corr.model") 
      if (is.null(corr.model)) warning('attr(next_LMatrix, "corr.model") is NULL')
      if (corr.model=="random-coef") {
        latentL_blob <- attr(lmatrix,"latentL_blob")
        design_u <- latentL_blob$design_u
        # Old versions fitted a model with L (->strucList) based on design$u and lambda_est containing the $d
        # if the d's were not 1,  some code here previously aimed to put them back into the compact L. 
        # BUT then, the long L (=>strucList) and the lambda could be inconsistent with the compact L. 
        # For example, the lambda's are distinctly needed in (predVar with new data) -> .get_logdispObject() -> (dvdloglamMat_needed)
        # This bug led to negligible imprecisions in tests (see successive versiosn of test-devel-predVar-rC_spprec.R)
        # until composite corrMatrix models were tested. 
        # Instead of fixing a lot of stuff, it appeared better to use latent_d=1 even for spprec.
        # Finally (v3.9.19) $d was removed from the latentL_blob. Only post-fit code may need to check its existence and then use it, for back compatibility. 
        if ((kappa(design_u)>1e06)) { # occurs in HLfit3 example; also ./. 
          # testthat check for fit_small in test-devel-predVar-ranCoefs)
          latentL_blob$gmp_compactcovmat <- gmp::as.bigq(latentL_blob$compactcovmat)
          # the gmp bug [str() failing] reported in Jan 2020 is still there in version 0.6.2.
          # MRE:  str(gmp::as.bigq(matrix(seq(4),ncol=2)))
          latentL_blob$gmp_compactprecmat <- gmp::solve.bigq(latentL_blob$gmp_compactcovmat)
          latentL_blob$compactchol_Q_w <- t(gmp::asNumeric(gmp::solve.bigq(gmp::as.bigq(design_u)))) # returns bigq matrix 
        } else if (is.null(latentL_blob$compactchol_Q)) {
          latentL_blob$compactchol_Q_w <- t(solve(design_u))
        } else {
          latentL_blob$compactchol_Q_w <- latentL_blob$compactchol_Q
        }
        latentL_blob$compactchol_Q <- NULL
        latentL_blob$design_u <- NULL
        attr(lmatrix,"latentL_blob") <- latentL_blob 
      } else if (! is.null(msd.arglist <- attr(lmatrix,"msd.arglist"))) { ## Matern, Cauchy
        ## remove potentially large distMatrix to save space, but stores its "call" attribute for getDistMat()
        if ( ! is.null(dM <- msd.arglist$distMatrix) ) { 
          if ( ! is.null(distcall <- attr(dM,"call"))) {
            msd.arglist$distcall <- distcall ## eg language proxy::dist(x = uniqueGeo, method = dist.method)
            msd.arglist$distMatrix <- NULL ## removes the big matrix
          }
          attr(lmatrix,"msd.arglist") <- msd.arglist
        }
      }
      # need_gmp attribute does not control use of above gmp object; rather it control use of gmp in other cases
      # 'given' that I cannot modify strucList in post_fit code, the need_gmp attribute must be assigned now.  
      attr(lmatrix,"need_gmp") <- (
        ( ! inherits(lmatrix,"dCHMsimpl") ) && ## before any test on type, because dCHMsimpl has a @type slot; very lately caught "issue "length > 1 in coercion to logical"
          attr(lmatrix,"type")=="cholL_LLt" && corr.model %in% c("Matern","Cauchy") && kappa(lmatrix)>1e05
      ) # in other cases, we might still 'need' it, but .get_invColdoldList) might not be ready to handle it.
      ## Now fill strucList[[it]]
      #
      if (corr.model=="random-coef") { 
        #
        if (inherits(lmatrix,"dCHMsimpl")) { 
          #   # ranef() assumes strucList[[it]] is a mMatrix so we cannot simply copy the dCHMsimpl.
          lmatrix_ <- as(lmatrix,"pMatrix") %*% solve(lmatrix, system="Lt", b=.sparseDiagonal(n=ncol(lmatrix), shape="g")) 
        } else {
          lmatrix_ <- next_LMatrices[[rd]]
        }
        #
        if (isS4(lmatrix_)) {
          extras_attr <- setdiff(names(attributes(lmatrix)),c("class",slotNames(lmatrix)))
        } else extras_attr <- setdiff(names(attributes(lmatrix)),names(attributes(strucList[[rd]])))
        #
        if (inherits(lmatrix_,"Kronfacto")) lmatrix_ <- lmatrix_@BLOB$long # because some the prediction code does not (yet) handle Kronfacto: test_composite l. 181 or 182
        #
        # Add descriptor of the correlation model necessary to predict with new data:
        #
        if (processed$ranCoefs_blob$is_composite[rd]) {
          # nothing to do yet
          ## # ___TAG___ code potentially needed here  to extend composite ranefs beyond corrMatrix
          
          ### for corrMatrix:
          # * processed$corr_info$cov_info_mats[[.]] contains (~user-level) info for corrMatrix, plus a handy "blob" attribute
          # that should be used whether the info is a list (SPPREC) or not (CORREL).
          # It holds the "short" $Lunique and possibly other elements
          
          ## For CORREL algos:
          # * next_LMatrices and ranCoefs_blob$LMatrices (but not processed$ranCoefs_blob...) should contain an rC_blob_longLvMatrix;
          # * processed$ranCoefs_blob$LMatrices contains a useless matrix (until something is changed in the code);
          # * processed$AUGI0_ZX$envir$LMatrices has a NULL element
          #      ( BUT for standard corrMatrix, the SPPREC code -> .Lunique_info_from_Q_CHM() + assign Lunique in $AUGI0_ZXenvir$LMatrices).;
        }
        attributes(lmatrix_)[extras_attr] <- attributes(lmatrix)[extras_attr] 
        #if (all(latent_d ==1)) attr(lmatrix_,"latentL_blob")[["d"]] <- NULL # so that we know they are all 1 without testing them
        attr(lmatrix_,"latentL_blob")[["d"]] <- NULL # so that we know they are all 1 without testing them
        strucList[[rd]] <- lmatrix_ 
      } else strucList[[rd]] <- lmatrix  
    } 
  }
  return(list(strucList = structure(strucList, isRandomSlope=ranCoefs_blob$isRandomSlope),## isRandomSlope is boolean vector
              v_h = structure(v_h, u_h= structure(u_h, cum_n_u_h=cum_n_u_h)))) ## cum_n_u_h  duplicates info in lambda object. But it's handy.)) 
}


.nothing_to_fit <- function(phi.Fix, off, models, etaFix, rand.families, cum_n_u_h, 
                            lambda.Fix, vec_n_u_h, n_u_h, fixed_adjacency_info, ZAL, BinomialDen, processed) {
  ## nothing to fit. We just want a likelihood
  ### a bit the same as max.iter<1 ... ?
  phi_est <- phi.Fix
  eta <- off
  if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
    ## we need u_h in calc_APHLS...() and v_h here for eta...
    v_h <- etaFix$v_h
    u_h <- etaFix$u_h
    if (is.null(u_h)) u_h <- processed$u_h_v_h_from_v_h(v_h)
    lambda_est <- .resize_lambda(lambda.Fix,vec_n_u_h,n_u_h, adjacency_info=fixed_adjacency_info)
    eta <- eta + drop(ZAL  %id*id%  etaFix$v_h) ## updated at each iteration
  } ## FREQS
  ## conversion to mean of response variable (COUNTS for binomial)
  muetablob <- .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed, phi_est=phi_est) 
  w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est, obsInfo=processed$how$obsInfo) ## 'weinu', must be O(n) in all cases
  if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
    wranefblob <- processed$updateW_ranefS(u_h=u_h,v_h=v_h,lambda=lambda_est)
    # at this point w.resid is always the result of .calc_w_resid()
    # and when it is a list with info about mv model it has a complete vector $w_resid.
    H_w.resid <- .calc_H_w.resid(w.resid, muetablob=muetablob, processed=processed) # for LLF w.resid is not generally defined.
    H_global_scale <- .calc_H_global_scale(H_w.resid)
    weight_X <- .calc_weight_X(H_w.resid, H_global_scale, obsInfo=processed$how$obsInfo) ## sqrt(s^2 W.resid)  # -> .... sqrt(w.resid * H_global_scale)
    ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
    Xscal <- .make_Xscal(ZAL, ZAL_scaling, processed=processed)
    sXaug <- do.call(processed$corr_method,
                     list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
  } else sXaug <- NULL 
  res <- list(APHLs=.calc_APHLs_from_ZX(processed=processed, which="p_v", sXaug=sXaug, phi_est=phi_est, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, muetablob=muetablob))
  return(res)
}

.calc_APHLs_from_params <- function(muetablob, phi_est, lambda_est, u_h,v_h, ZAL, processed, which, ad_hoc_corrPars) {
  w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est, obsInfo=processed$how$obsInfo) ## 'weinu', must be O(n) in all cases
  wranefblob <- processed$updateW_ranefS(u_h=u_h,v_h=v_h,lambda=lambda_est)
  # at this point w.resid is always the result of .calc_w_resid()
  # and when it is a list with info about mv model it has a complete vector $w_resid.
  H_w.resid <- .calc_H_w.resid(w.resid, muetablob=muetablob, processed=processed) # for LLF w.resid is not generally defined.
  H_global_scale <- .calc_H_global_scale(H_w.resid)
  weight_X <- .calc_weight_X(H_w.resid, H_global_scale, obsInfo=processed$how$obsInfo) ## sqrt(s^2 W.resid)  # -> .... sqrt(w.resid * H_global_scale)
  ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
  if (processed$is_spprec) {
    sXaug <- do.call(processed$spprec_method, # ie, def_AUGI0_ZX_spprec
                     list(AUGI0_ZX=processed$AUGI0_ZX, corrPars=ad_hoc_corrPars, 
                          cum_n_u_h=processed$cum_n_u_h, w.ranef=wranefblob$w.ranef, H_w.resid=H_w.resid))
  } else {
    Xscal <- .make_Xscal(ZAL, ZAL_scaling, processed=processed, as_matrix=.eval_as_mat_arg(processed))
    sXaug <- do.call(processed$corr_method, 
                          list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale)) 
  }
  APHLs <- .calc_APHLs_from_ZX(processed=processed, which=which, sXaug=sXaug, phi_est=phi_est, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, muetablob=muetablob)
  return(APHLs)
}

.adhoc_rbind_dgC_dvec <- function(X, dvec) { 
  Ilen <- X@Dim[2L]
  newlen <- Ilen+length(X@x) # ne @x length
  Iseq <- seq_len(Ilen)
  Ip <- c(0L,Iseq)
  newp <- Ip+X@p
  Ipos <- newp[-1L]
  #
  newx <- numeric(newlen)
  newx[Ipos] <- dvec
  newx[-Ipos] <- X@x
  newi <- integer(newlen)
  newi[Ipos] <- X@Dim[1L]-1L+Iseq
  newi[-Ipos] <- X@i
  #
  X@i <- newi
  X@x <- newx
  X@p <- newp
  X@Dim[1L] <- Ilen+X@Dim[1L]
  if ( ! is.null(X@Dimnames[[1L]])) X@Dimnames[[1L]] <- c(X@Dimnames[[1L]],rep("",Ilen)) 
  return(X)
}

.adhoc_rbind_dtC_dvec <- function(X, dvec) { ## this works but is not used 
  Ilen <- X@Dim[2L]
  newlen <- Ilen+length(X@x) # ne @x length
  Iseq <- seq_len(Ilen)
  Ip <- c(0L,Iseq)
  newp <- Ip+X@p
  Ipos <- newp[-1L]
  #
  newx <- numeric(newlen)
  newx[Ipos] <- dvec
  newx[-Ipos] <- X@x
  newi <- integer(newlen)
  newi[Ipos] <- Ilen-1L+Iseq
  newi[-Ipos] <- X@i
  #
  Dimnames <- X@Dimnames
  if ( ! is.null(Dimnames[[1L]])) Dimnames[[1L]] <- c(Dimnames[[1L]],rep("",Ilen))
  X <- new("dgCMatrix",i=newi, p=newp, x=newx, Dim=X@Dim+c(Ilen,0L), Dimnames=Dimnames)
  return(X)
}

.XDtemplate <- function(X, upperTri) {
  if (inherits(X,"dtCMatrix")) {
    #XDtemplate <- .adhoc_rbind_dtC_dvec(X,dvec=rep(1,ncol(X))) # equivalent but with less ad-hoc unsafe code:
    if (.spaMM.data$options$Matrix_old) { # ugly... but such versions do not handle as(, "dMatrix"))
      X <- as(X, "dgCMatrix") 
    } else X <- as(X,"generalMatrix") 
    res <- .adhoc_rbind_dgC_dvec(X,dvec=rep(1,ncol(X)))
  } else if (inherits(X,"dgCMatrix")) {
    res <- .adhoc_rbind_dgC_dvec(X,dvec=rep(1,ncol(X)))
  } else if (inherits(X,"Matrix")) {
    res <- .adhoc_rbind_dgC_dvec(X,dvec=rep(1,ncol(X)))
  # } else if (upperTri) {
  #   # test whether one can get advantage of sparse QR in .damping_to_solve() (=> No)
  #   res <- .adhoc_rbind_dgC_dvec(as(X,"dgCMatrix"),dvec=rep(1,ncol(X)))
  } else res <- rbind(X, diag(nrow=ncol(X)))
  attr(res,"upperTri") <- upperTri
  res
}

.damping_to_solve <- function(X, XDtemplate=NULL, dampDpD, rhs=NULL,method="QR", .drop=TRUE) { ## cf my notes on Mor\'e 1977 
  if (.drop) rhs <- drop(rhs) ## 1-col m/Matrix to vector ## affects indexing below but the result of the chol2inv line is always Matrix
  if (method=="QR") { ## seems always true
    if (is.null(XDtemplate)) {
      warning("Possibly inefficient call to .damping_to_solve() without precomputed XDtemplate")
      XDtemplate <- .XDtemplate(X, upperTri=FALSE)
    }
    nr <- nrow(XDtemplate)-length(dampDpD)
    if (inherits(XDtemplate,"Matrix")) { # both cases occur in routine use # presumably spprec, ou sparse-correl
      XD <- .Dvec_times_Matrix_lower_block(Dvec=sqrt(dampDpD),X=XDtemplate,min_row=nr)
      RP <- qr(XD) ## i.e. Matrix::qr
      RRsP <- sort.list(RP@q) 
      # solve(crossprod(XD)) = chol2inv(Matrix::qrR(RP,backPermute = FALSE))[RRsp,RRsP]
      # solve(crossprod(XD), rhs) = (Matrix::chol2inv(Matrix::qrR(RP, backPermute = FALSE)) %*% rhs[sort.list(RRsP)])[RRsP]
      if (is.null(rhs)) {
        return(list(inv=Matrix::chol2inv(qrR(RP,backPermute = FALSE)),Rperm=RP@q+1L,RRsP=RRsP))
      } else {
        # test cloglog (with large enough maxtime) and test dhglm reaches here
        # test DHGLM ncol(XD)=21 or 23 (v or v_b) is already faster by .backsolve()
        # test cloglog ncol(XD)=144 or 192 is indistinct
        # Other tests would be 'bigranefs' with forced LevM; twinR
        # apparently Matrix::chol2inv -> tcrossprod(solve(as(x, "triangularMatrix")))  # cf getMethod("chol2inv",signature=c(x="sparseMatrix"))
        ## very slow tcrossprod for large matrices,
        ## not specifically bc of the tcrossprod implementation (compared to crossprod) but because the result is *much* denser (bigranefs example) 
        ## It's then much better to avoid chol2inv(), at least when there is a RHS!
        ## PLUS 05/2020 : more efficient .backsolve avoids backsolve -> as.matrix(). 
        ## BUT 01/2022 : the naive solve(., solve(t(.), ..)) on dtC avoids any conversion => much faster
        if (TRUE) {
          qrr <- qrR(RP,backPermute = FALSE) # must be dtC
          resu <- solve(qrr, solve(t(qrr), rhs[RP@q + 1])) # Matrix::solve for dtC # faster than mess below, test-nloptr's simuland example.
          if (.drop) { # rhs has been converted to vector so is not matrix
            resu <- resu[RRsP,1]
          } else resu <- resu[RRsP,]
        } else { # all attempts compared
          # if (FALSE) { 
          #   if (is.matrix(rhs)) {
          #     resu <- (Matrix::chol2inv(qrR(RP,backPermute = FALSE)) %*% rhs[RP@q+1L,])[RRsP,]
          #   } else { ## if rhs is one col 
          #     resu <- (Matrix::chol2inv(qrR(RP,backPermute = FALSE)) %*% rhs[RP@q+1L])[RRsP]
          #   }
          # } else {
          #   if (TRUE) {
          #     qrr <- qrR(RP,backPermute = FALSE) # must be dtC
          #     if (TRUE) {
          #       resu <- solve(qrr, solve(t(qrr), rhs[RP@q + 1])) # Matrix::solve for dtC # faster than mess below, test-nloptr's simuland example.
          #       if (.drop) { # rhs has been converted to vector so is not matrix
          #         resu <- resu[RRsP,1]
          #       } else resu <- resu[RRsP,]
          #     } else {
          #       qrr <- as(qrr,"dgCMatrix") ## avoid conversions in each call to .backsolve()
          #       #if (TRUE) { # Alternatives compared while debugging a large negbin GLMM (from twins ms), => same results but faster by .backsolve
          #       # BUT ... the singla as(.,"dgCMatrix") is much more costly...
          #       resu <- .backsolve(qrr, .backsolve(qrr, rhs[RP@q + 1], transpose = TRUE))
          #       if (is.matrix(rhs)) {
          #         resu <- resu[RRsP,]
          #       } else resu <- resu[RRsP,1]
          #       # } else {
          #       #   resu <- backsolve(qrr, backsolve(qrr, rhs[RP@q + 1], transpose = TRUE))
          #       #   if (is.matrix(rhs)) {
          #       #     resu <- resu[RRsP,]
          #       #   } else resu <- resu[RRsP]
          #       # }
          #     }
          #   } else {
          #     solveR <- solve(qrR(RP,backPermute = FALSE))
          #     if (is.matrix(rhs)) {
          #       resu <- (solveR %*% .crossprod(solveR,rhs[RP@q+1L,]))[RRsP,]
          #     } else { ## if rhs is one col 
          #       resu <- (solveR %*% .crossprod(solveR,rhs[RP@q+1L]))[RRsP]
          #     }
          #   }
          # }
          # 
        }
      }
      ## assuming rhs is 'slim', this minimizes permutations. For computation of dv_h, it is square rather than slim...
    } else { # decorr
      XD <- .Dvec_times_matrix(Dvec=c(rep(1,nr),sqrt(dampDpD)), X=XDtemplate)
      if (TRUE) { 
        # on test-levM.R arabidopsis example, .lmwithQR faster than .lmwithQRP faster than a firsr .update_R_in_place_essai() (code in givens.R) 
        # Finally a revised .update_R_in_place() may be marginally faster than .lmwithQR()
        # check: test-spaMM.R#8: try(fitme(I(1 + cases) ~ I(prop.ag/10)...
        if (attr(XDtemplate,"upperTri")) { # The update does not work when ".lmwithQRP() has been used.
          R_scaled <- .update_R_in_place(XD) # bypass QR factorization code to get new R (for unknown new Q) with XD modified in place.
        } else R_scaled <- .lmwithQR(XD,yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled ## Eigen QR OK since we don't request Q
         # dVscaled <- chol2inv(RP$R_scaled)[RRsP,RRsP] %*% rhs
        if (is.null(rhs)) {
          return(list(inv=chol2inv(R_scaled),Rperm=seq(ncol(R_scaled)),RRsP=seq(ncol(R_scaled))))
        } else {
          # here again I tested the speed of .Rcpp_chol2solve(R_scaled, rhs) and its slower; possible reason is that a new 
          # matrix is allocated (solve in place not possible) while this is avoided by:
          resu <- .Rcpp_backsolve(R_scaled, .Rcpp_backsolve(R_scaled,rhs,transpose=TRUE))
          # which seems to bring a slight gain over
          # resu <- backsolve(R_scaled, backsolve(R_scaled,rhs,transpose=TRUE))
        }
      } else {
        RP <- .lmwithQRP(XD,yy=NULL,returntQ=FALSE,returnR=TRUE) ## Eigen QR OK since we don't request Q
        RRsP <- sort.list(RP$perm)
        # dVscaled <- chol2inv(RP$R_scaled)[RRsP,RRsP] %*% rhs
        if (is.null(rhs)) {
          return(list(inv=chol2inv(RP$R_scaled),Rperm=RP$perm+1L,RRsP=RRsP))
        } else {
          if (is.matrix(rhs)) {
            ## case never run, hence not tested.
            # resu <- (chol2inv(RP$R_scaled) %*% rhs[RP$perm+1L,])[RRsP,]  ## assuming rhs is 'slim', this minimizes permutations   
            # resu <- .Rcpp_chol2solve(RP$R_scaled,rhs[RP$perm+1L,])[RRsP,] 
            resu <- backsolve(RP$R_scaled, backsolve(RP$R_scaled,rhs[RP$perm+1L,],transpose=TRUE))[RRsP,]
          } else { ## if rhs is one col 
            ## test-spaMM replicated 20 times ->  .Rcpp_chol2solve has no obvious benefit
            # resu <- .chol2solve(RP$R_scaled,rhs[RP$perm+1L])[RRsP] 
            resu <- backsolve(RP$R_scaled, backsolve(RP$R_scaled,rhs[RP$perm+1L],transpose=TRUE))[RRsP] 
          }
        }
      }
    }
    if (.drop) resu <- drop(resu)
    return(resu)
  } else { # other 'method'; not used.
    diag(X) <- diag(X)+dampDpD
    if (inherits(X,"Matrix")) {
      if (is.vector(rhs)) {
        return(drop(Matrix::solve(X,rhs)))
      } else return(Matrix::solve(X,rhs))
    } else {      
      return(solve(X,rhs))
    }
  }
}


.new_phifit_init_corrPars <- function(phifit, innershift) { ## innershfit has remarkable impact (even cyclical behaviour for 1e-4)
  corrPars <- get_ranPars(phifit,which="corrPars")
  type <- attr(corrPars,"type")
  notouter <- which((u_list <- unlist(type))!="outer")
  corrPars <- structure(  .remove_from_cP(corrPars, u_names = names(notouter)),
                             type=.remove_from_cP(type, u_list=u_list,u_names = names(notouter)) )
  # correction corrPars to handle numerical problems; presumably nothing req for IMRFs bc better conceived control of upper kappa
  if (length(corrPars)) {
    ## FR->FR init must be in user scale, optr output is 'unreliable' => complex code
    LUarglist <- attr(phifit,"optimInfo")$LUarglist ## might be NULL 
    isMatern <- ("Matern" == LUarglist$corr_types)
    isCauchy <- ("Cauchy" == LUarglist$corr_types)
    if ( any(isMatern | isCauchy)) {
      LowUp <- do.call(".makeLowerUpper",LUarglist)
      for (rd in which(isMatern)) {
        char_rd <- as.character(rd)
        ## need to redefine NUMAX and RHOMAX
        if ( !is.null(nu <- corrPars[[char_rd]][["nu"]]) ) {
          if (innershift>0) {
            NUMAX <- LUarglist$moreargs[[char_rd]][["NUMAX"]] 
            initnu <- .nuFn(nu,NUMAX=NUMAX) 
            mini <- LowUp$lower$corrPars[[char_rd]]$trNu
            maxi <- LowUp$upper$corrPars[[char_rd]]$trNu
            margin <- innershift *(maxi-mini)/2
            initnu <- max(min(initnu,maxi-innershift),mini+innershift)
            corrPars[[char_rd]][["nu"]] <- .nuInv(initnu,NUMAX=NUMAX)
          } else corrPars[[char_rd]][["nu"]] <- nu
        }
        if ( !is.null(rho <- corrPars[[char_rd]][["rho"]]) ) {
          if (innershift>0) {
            RHOMAX <- LUarglist$moreargs[[char_rd]][["RHOMAX"]]  
            initrho <- .rhoFn(rho,RHOMAX=RHOMAX) 
            mini <- LowUp$lower$corrPars[[char_rd]]$trRho
            maxi <- LowUp$upper$corrPars[[char_rd]]$trRho
            margin <- innershift *(maxi-mini)/2
            initrho <- max(min(initrho,maxi-innershift),mini+innershift)
            corrPars[[char_rd]][["rho"]] <- .rhoInv(initrho,RHOMAX=RHOMAX)
          } else corrPars[[char_rd]][["rho"]] <- rho
        }
      }
      for (rd in which(isCauchy)) {
        char_rd <- as.character(rd)
        if ( !is.null(longdep <- corrPars[[char_rd]][["longdep"]]) ) {
          if (innershift>0) {
            LDMAX <- LUarglist$moreargs[[char_rd]][["LDMAX"]] 
            initlongdep <- .longdepFn(longdep,LDMAX=LDMAX) 
            mini <- LowUp$lower$corrPars[[char_rd]]$trLongdep
            maxi <- LowUp$upper$corrPars[[char_rd]]$trLongdep
            margin <- innershift *(maxi-mini)/2
            initlongdep <- max(min(initlongdep,maxi-innershift),mini+innershift)
            corrPars[[char_rd]][["longdep"]] <- .longdepInv(initlongdep,LDMAX=LDMAX)
          } else corrPars[[char_rd]][["longdep"]] <- longdep
        }
        if ( !is.null(rho <- corrPars[[char_rd]][["rho"]]) ) {
          if (innershift>0) {
            RHOMAX <- LUarglist$moreargs[[char_rd]][["RHOMAX"]]  
            initrho <- .rhoFn(rho,RHOMAX=RHOMAX) 
            mini <- LowUp$lower$corrPars[[char_rd]]$trRho
            maxi <- LowUp$upper$corrPars[[char_rd]]$trRho
            margin <- innershift *(maxi-mini)/2
            initrho <- max(min(initrho,maxi-innershift),mini+innershift)
            corrPars[[char_rd]][["rho"]] <- .rhoInv(initrho,RHOMAX=RHOMAX)
          } else corrPars[[char_rd]][["rho"]] <- rho
        }
      }
    }
  }
  return(corrPars)
}

.HLfit_finalize_init_lambda <- function(models, init.lambda, processed, ZAL, cum_n_u_h, vec_n_u_h, n_u_h, ranCoefs_blob) {
  ## we typically need  both init.lambda and lambda_est... 
  if ( anyNA(init.lambda) ) { ## some inner may remain undetermined
    NA_in_compact <- is.na(init.lambda)
    prev_lambda_est <- processed$port_env$port_fit_values$lambda_est
    if (is.null(prev_lambda_est)) {
      stillNAs <- which(NA_in_compact)
      ## if reset is TRUE init.lambda is recomputed. Otherwise it is computed only once and then 'got'
      if (is.null(processed$envir$init_lambda_guess)) { 
        ## FIXME: condition not clearly pertinent as init.lambda should generally depend on LMatrix. But see comment in .eval_init_lambda_guess
        ## which occurs either if reset is TRUE or if $get_init_lambda has not yet been called
        processed$envir$init_lambda_guess <- .eval_init_lambda_guess(processed, stillNAs=stillNAs, ZAL=ZAL, cum_n_u_h=cum_n_u_h, For="iter")
        ## HLfit call -> default init_lambda_guess is >= 1e-4, bc of explicit code for For!= "optim" in .eval_init_lambda_guess()
        ## fitme call -> even non-default init values are >=1e-4 bc of .calc_inits_dispPars()
      }
      init.lambda[stillNAs] <- processed$envir$init_lambda_guess[stillNAs]
      lambda_est <- .resize_lambda(init.lambda,vec_n_u_h,n_u_h, adjacency_info=NULL)
    } else { ## do not replace current fixed values !
      lambda_est <- .resize_lambda(init.lambda,vec_n_u_h,n_u_h, adjacency_info=NULL)
      NA_in_expanded <- .resize_lambda(NA_in_compact,vec_n_u_h,n_u_h, adjacency_info=NULL)
      lambda_est[NA_in_expanded] <- prev_lambda_est[NA_in_expanded]
    } 
  } else lambda_est <- .resize_lambda(init.lambda,vec_n_u_h,n_u_h, adjacency_info=NULL)
  for (rd in seq_len(length(vec_n_u_h))) {
    if ( ! is.null(prior_lam_fac <- processed$rand.families[[rd]]$prior_lam_fac)) {
      u.range <- (cum_n_u_h[rd]+1L):(cum_n_u_h[rd+1L])
      lambda_est[u.range] <- lambda_est[u.range]*prior_lam_fac
    }
  }
  for (rd in which(ranCoefs_blob$is_set)) { 
    u.range <- (cum_n_u_h[rd]+1L):cum_n_u_h[rd+1L]
    lambda_est[u.range] <- ranCoefs_blob$lambda_est[u.range]
  }
  attr(lambda_est,"init.lambda") <- init.lambda ## keep init.lambda provisionally bc SEMwrap expects it (FIXME)
  return(lambda_est) 
}

.calc_initial_init_lambda <- function(lam_fix_or_outer_or_NA, nrand, processed, ranCoefs_blob, init.HLfit, fixed) {
  init.lambda <- lam_fix_or_outer_or_NA # already the right size 'nrand' with NA's or non-fixed ones"
  if (is.null(lambdaType <- processed$envir$lambdaType)) { # _F I X M E_ that's midway from pushing this to .preprocess. But...
    # the hyper stuff cannot be fully included in this block and I will wait a better opportunity to dissect what could be done about it. 
    lambdaType <- rep("",nrand) 
    lambdaType[processed$ranCoefs_blob$is_set] <- "fix_ranCoefs" 
    lambdaType[lambdaType=="" & ! is.na(processed$lambda.Fix)] <- "fixed" 
    processed$envir$lambdaType <- lambdaType
  }
  
  if ( ! is.null((hyper_info <- processed$hyper_info)$map)) {
    hyper_type <- attr(fixed,"type")$hyper
    hy_lam_types <- character(length(hyper_type))
    names(hy_lam_types) <- hy_names <- names(hyper_type)
    for (char_hyper_it in hy_names) { 
      hy_lam_type <- hyper_type[[char_hyper_it]]$hy_trL # if outer optim
      if (is.null(hy_lam_type)) hy_lam_type <- "fix" # hyper_type[[char_hyper_it]]$hy_lam # if no outer optim is necessary
      #if (is.null(hy_lam_type)) stop("is.null(hy_lam_types[char_hyper_it])")
      hy_lam_types[char_hyper_it] <- hy_lam_type
    }
    hy_lam_types[] <- paste0(hy_lam_types,"_hyper")
    pos <- .unlist(hyper_info$ranges[hy_names])  
    lambdaType[pos] <- hy_lam_types[hyper_info$map[pos]]
  }
  lambdaType[lambdaType=="" & ranCoefs_blob$is_set] <- "outer_ranCoefs" 
  ## "fixed" (not "fix" -- confusing) is tested eg to compute dfs p_lambda, in summary, in calc_logdisp_cov.R, in .get_logdispObject.
  lambdaType[(lambdaType=="" & ! (is.na(lam_fix_or_outer_or_NA)))] <- "outer" # outer charac. by difference between processed$lambda.Fix and lambda.Fix
  lambdaType[whichInner <- (lambdaType=="")] <- "inner"  
  if (! is.null(init.HLfit$lambda)) init.lambda[whichInner] <- init.HLfit$lambda[whichInner]
  if (! is.null(init.lambda)) attr(init.lambda,"type") <- lambdaType
  return(init.lambda)
}

.post_process_resid_corrPars <- function(phifit_init, resid_corrPars) {
  if ( ! is.null(resid_corrPars)) {
    ## The corrPars (considering only outer-optimized ones so far: minor fixme)
    rPtype <- attr(resid_corrPars,"type") ## hinges on the fact that the input ranFix of HLfit_body must have info about "outer" corrPars!
    # but (FIXME) using attr(get_ranPars(phifit,which="corrPars"),"type"), similarly to lambda, could be cleaner. 
    ## through corrHLfit or fitme call: ranPars inherits values from <'corrfitme'> (...,init.HLfit(...))
    u_rPtype <- unlist(rPtype)
    is_outer <- (u_rPtype=="outer")
    ## init.HLfit must recover elements from ranPars! (bug detected by AIC( <SEM fit> ) in test-CAR.R where it must get rho...
    if ( ! is.null(names(which(is_outer)))) { ## can be NULL for corrMatrix case => not $ranFix    
      not_outer_names <- names(which( ! is_outer))
      resid_corrPars <- structure(.remove_from_cP(resid_corrPars,u_names=not_outer_names), ## loses attributes
                                  type=.remove_from_cP(rPtype,u_rPtype, u_names=not_outer_names) ) 
      phifit_init$corrPars <- .modify_list(phifit_init$corrPars, 
                                           resid_corrPars) ## loses attributes
    }
  } 
  return(phifit_init)
}

.update_port_fit_values_residModel <- function(residProcessed, phimodel, residModel, PHIblob) {
  if (length(phimodel)>1L) {
    for (mv_it in seq_along(phimodel)) {
      if (phimodel[mv_it]=="phiHGLM") .update_port_fit_values_residModel(residProcessed[[mv_it]], # RHS  residProcessed being actually the 'residProcesseds' list
                                                                         phimodel=phimodel[mv_it],
                                                                         residModel=residModel[[mv_it]], # RHS residModel being actually the 'residModels' list
                                                                         PHIblob= PHIblob[[mv_it]] # RHS  PHIblob being actually the 'PHIblob$multiPHI' list
      )
    }
    return()   # leaves the fn but no return value, port_env environments of each residProcessed modified
  }
  phifit <- PHIblob$phifit
  port_env <- residProcessed$port_env
  # phifit (fitme -> HLfit) has already written in processed$residProcessed$port_env$port_fit_values (and this may be erased by the phifitarglist code): 
  # but the condition was on its own APHLs and it may actually have erased the info 
  ## To be used only in iter=0; for iter>0 see parallel code in HLfit_body()
  phifit_init_HLfit <- list(v_h=phifit$v_h)
  if (is.null(residModel$fixed$fixef)) phifit_init_HLfit$fixef <- fixef(phifit) ## : Use an init if the parameters are not fixed. 
  port_env$port_fit_values[["init_HLfit"]] <- phifit_init_HLfit
  port_env$port_fit_values$corrPars <- .new_phifit_init_corrPars(phifit,innershift=0) 
  ## $corrPars : canonical and with full info about types including "outer"
  port_env$port_fit_values$lambda.object <- phifit$lambda.object ## disp fit lambda info.
  # no return value, port_env environment modified
}


.update_port_fit_values <- function(old_obj, new_obj, port_fit_values, models, processed, control.HLfit,
                                    lambda_est, #we only need lambda_est to initiate the next inner optim of the mean fit
                                    PHIblob # we need more info to initiate both inner and outer optim of the next phifit
                                    ) {
  if ( new_obj> old_obj) {
    abs_d_obj <- abs(old_obj-new_obj)
    ## Use FALSE to inhibit all port_env usage:
    if ( ! identical(control.HLfit$write_port_env,FALSE)) assign("objective",new_obj,envir=processed$port_env) 
    if (abs_d_obj<length(port_fit_values$fixef) ||  # previously '5' based on some old Loaloa IMRF example; changed in v3.11.10, 04/2022
        processed$LevenbergM["LM_start"]) { # LM_start condition appears to improve the fit for test negbin difficult.
      # Update port_fit_values and assign it: 
      if (models[["eta"]]=="etaHGLM") {
        port_fit_values$lambda_est <- lambda_est ## mean fit lambda; to be used by .HLfit_finalise_init_lambda
        # 
      }
      processed$port_env$port_fit_values <- port_fit_values # NOT only for phiHGLM
      #
      if ( ! is.null(processed$residModels)) {
        .update_port_fit_values_residModel(residProcessed=processed$residProcesseds,
                                           phimodel=models[["phi"]],
                                           residModel=processed$residModels,
                                           PHIblob=PHIblob$multiPHI)
      } else if (models[["phi"]]=="phiHGLM") .update_port_fit_values_residModel(residProcessed=processed$residProcessed,
                                                                                phimodel=models[["phi"]],
                                                                                residModel=processed$residModel,
                                                                                PHIblob=PHIblob
      )
    }
  } else { ## keep objective value; two cases for init values:
    if (abs(old_obj-new_obj)<1) {
      ## small decrease: Keep old port_env_values
    } else {
      processed$port_env$port_fit_values <- NULL ## remove starting value, not useful for large variations 
      if ( ! is.null(processed$residModels)) {
        for (mv_it in seq_along(models[["phi"]])) if (models[["phi"]][[mv_it]]=="phiHGLM") processed$residProcesseds[[mv_it]]$port_env$port_fit_values <- NULL
      } else processed$residProcessed$port_env$port_fit_values <- NULL
    }  
  }
}

# .unscale_beta_cov_info def removed in [ v2.6.52

.scale <- function(X, beta=NULL) {
  if (is.null(beta)) {
    Xattr <- attributes(X)
    nc <- ncol(X)
    if (nrow(X)==1L) {
      scale <- rep(1, nc)
    } else {
      scale <- numeric(nc)
      for (jt in seq_len(nc)) {
        x <- X[,jt]
        v <- var(x)
        if (v==0) { ## there should be only one constcol, "(Intercept)", or something that happens to play the same role
          scale[jt] <- x[1L]
        } else if ((am <- abs(mean(x)))>10) { ## threshold ad hoc but the sensitive tests are eg HLfit3 and fitme(passengers ~ month * year + AR1(1|time)...)
          scale[jt] <- sqrt(v)*(1+sqrt(am)) # not using sqrt(v) in this case affects adjacency-long... X1 has large mean and var
        } else scale[jt] <- sqrt(v)
      }
    }
    names(scale) <- colnames(X)
    X <- .m_Matrix_times_Dvec(X,1/scale)  # .m_Matrix_times_Dvec() kindly preserve attributes
    attr(X,"scaled:scale") <- scale ## same attr as in base::scale
    names_lostattrs <- setdiff(names(Xattr), names(attributes(X)))
    attributes(X)[names_lostattrs] <- Xattr[names_lostattrs] ## not mostattributes which messes S4 objects ?!
    return(X) ## center=FALSE keeps sparsity
  } else {
    return(beta * attr(X,"scaled:scale"))
  }
}

.unscale <- function(X, beta=NULL) {
  if (is.null(beta)) {
    Xattr <- attributes(X)
    X <- .m_Matrix_times_Dvec(X, attr(X,"scaled:scale")) 
    names_lostattrs <- setdiff(names(Xattr), names(attributes(X)))
    attributes(X)[names_lostattrs] <- Xattr[names_lostattrs] ## not mostattributes hich messes S4 objects ?!
    attr(X,"scaled:scale") <- NULL ## otherwise there is a non-trivial scale on an unscaled matrix
    attr(X,"info") <- "scale removed after unscaling" ## otherwise ther is a non-trivial scale on an unscaled matrix
    return(X)
  } else {
    return(beta / attr(X,"scaled:scale"))
  }
  # if (is.null(beta)) {
  #   X <- scale(X, center=FALSE,scale=1/attr(X,"scaled:scale"))
  #   return(scale(X, center=-attr(X,"scaled:center"),scale=FALSE))
  # } else {
  #   beta <- beta / attr(X,"scaled:scale")
  #   if ( ! is.null(intercept <- beta["(Intercept)"])) beta["(Intercept)"] <- intercept - sum(beta * attr(X,"scaled:center"))
  #   return(beta)
  # }
}

## as(.,"dtCMatrix) is quite slow => this attempt. But even the new() therein is slow.
.as_dtCMatrix_cholL <- function(cholL,template=NULL) { ## lower triangular matrix
  if (inherits(cholL,"dtCMatrix")) {
    return(cholL)
  } else {
    nc <- ncol(cholL) 
    indices <- matrix(0:(nc - 1L),ncol=nc,nrow=nc)
    ltri <- lower.tri(cholL,diag=TRUE)
    if (is.null(template)) {
      return(
        new("dtCMatrix", x= cholL[ltri], # '<sparse>[ <logic> ] : .M.sub.i.logical() maybe inefficient' if cholL is Matrix
            i=indices[ltri], p=c(0L,cumsum(rev(seq(nc)))), Dim= as.integer(c(nc,nc)),uplo="L")
      )
    } else {
      template@x <- cholL[ltri][template@x]
      return(template)
    }
  }
}

#### Old comment, presumably false: For sparse 'r', .Rcpp_backsolve has a definite advantage over backsolve() and possibly also solve().
# actually for dtC, the process was to convert to dgC (as seen here), for which .Rcpp_backsolve() is fast; but conversion is too slow.
# Instead one should use solve() for dtCMatrix
# Overall this makes .backsolve() useless. It was sparsely used. We keep it as a wrapper to check arguments for .Rcpp_backsolve
.backsolve <- function(r, x=NULL, upper.tri = TRUE, transpose = FALSE) {
  if (inherits(r,"Matrix")) { ## then dgCMatrix or dtCMatrix expected. 
    if ( ! inherits(r,"dgCMatrix")) { # .Rcpp_backsolve requires dgCMatrix # it would be nice to be able to pass a dtCMatrix to Eigen
      warning("*possibly inefficient call to .backsolve().") # in particular, for dtCMatrix, Matrix::solve() is more efficient (even with transposition) 
      if (.spaMM.data$options$Matrix_old) { # this block appears to evade the long tests
        r <- as(r,"dgCMatrix")
      } else r <- as(as(r,"generalMatrix") ,"CsparseMatrix") 
    }
  }
  if (inherits(x,"Matrix")) { 
    stop(".backsolve does not handle case inherits(x,'Matrix').")
    # .Rcpp_backsolve_M_M() does not work.
  } else return(.Rcpp_backsolve(r=r, x=x, upper_tri = upper.tri, transpose=transpose))
}

# Called by .add_phi_returns():
.get_phi_object <- function(phi.Fix, PHIblob, dev_res, prior.weights, phi.preFix, nobs=length(dev_res), control) {
  if (is.null(phi.Fix)) { # (ie not "fix" nor "outer") In mv case, this is called for each submodel, phi.Fix is not an mv-list
    beta_phi <- PHIblob$beta_phi 
    # I deleted some long-obsolete beta_phi renaming code and comments from 3.13.22 -> .23 here
    phi.object <- list(fixef=beta_phi, glm_phi=PHIblob$glm_phi, fittedPars=beta_phi)
    if (is.null(phi.object[["glm_phi"]])) {
      # delays computation of glm_phi
      glm_phi_args <- list(dev.res=dev_res*prior.weights,
                           control=control,
                           etastart=rep(PHIblob$beta_phi,nobs)) ## no glm <=> formula was ~1
      phi.object <- c(phi.object, list(glm_phi_args=glm_phi_args ) )
    } 
  } else { ## "fix" nor "outer" => important distinction for (summary, df or LRTs:
    if (is.null(phi.preFix)) { ## absent from original call
      phi.object <- list(phi_outer=structure(phi.Fix,type="var"), 
                         fittedPars=phi.Fix) 
    } else phi.object <- list(phi_outer=structure(phi.Fix,type="fix"))
  }
  return(phi.object)
}

## not used.
# .get_multi_phi_object <- function(phi.Fix, multiPHI, dev_res, prior.weights, processed, vec_nobs) {
#   multi_phi_object <- vector("list",length(vec_nobs))
#   cum_nobs <- c(0L,cumsum(vec_nobs))
#   for (mv_it in seq_along(vec_nobs)) {
#     resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
#     multi_phi_object[[mv_it]] <- .get_phi_object(phi.Fix[[mv_it]], PHIblob=multiPHI[[mv_it]], dev_res[resp_range], 
#                                               prior.weights[[mv_it]], 
#                                               phi.preFix=processed$phi.Fix[[mv_it]], 
#                                               vec_nobs[mv_it],
#                                               control=processed[["control.glm"]])
#   }
#   return(multi_phi_object)
# }

.update_phifitarglist <- function(processed, 
                                  residProcessed,
                                  residModel,
                                  dev.res, lev_phi, iter, verbose, phifit) {
  # Since we need to evaluate the phi-responses, using leverages, outer optimization would not avoid leverage computations.
  residProcessed$prior.weights <- structure((1-lev_phi)/2,unique=FALSE) # expected structure in 'processed'.
  # uses of prior weights matches that in input of calcPHI -> dispGammaGLM 
  residProcessed$data$.phi <- dev.res/(1-lev_phi) 
  residProcessed$y <- residProcessed$data$.phi
  residProcessed$iter_mean_dispFix <- max(200L,ceiling(100* mean(abs(log2(residProcessed$y)))))
  residProcessed$iter_mean_dispVar <- max(50L,ceiling(100* mean(abs(log2(residProcessed$y)))))
  residProcessed$main_terms_info$Y <- as.matrix(residProcessed$data$.phi)
  #  residProcessed$port_fit_values allows a previous 'close' [as defined by condition on logLik change] outer fit 
  #   to be used to initialize the 1st inner fitme of the new HLfit.
  #  residProcessed$fixed, until modified below, contains preprocessed resid.model's fixed lambda, and phi.
  #   and info in residModel$fixed must match (duplicate!) that in residProcessed:
  #   <processed=residProcessed>$lambda.Fix contains globally fixed lambda for the resid-fit.
  #  Finally, residProcessed$envir$ranPars may be used for outer-optimization at the mean-fit level 
  #   ('outer_optim_resid' option) of the resid-fit parameters,
  #   in which case residProcessed$envir$ranPars contains the objective's arguments that are resid-fit parameters,
  #   which have been written in (by [HLCor|HLfit].obj()).
  #   and should now be copied in phifitarglist$fixed:
  phifitarglist <- residModel
  # phifitarglist["init.HLfit"], formula, resid.model, family and rand.family have been <-NULL-ified and hoould not be put back. 
  ## $residModel: including formula, control.dist, family, fixed, rand.family... but most of the info used is in residProcessed!
  # it is tempting to select the arguments in formals(fitme_body) but this includes '...'
  phifitarglist["init.HLfit"] <- NULL ## for clarity (but such arguments must be kept in processed$residModel )
  phifitarglist["formula"] <- NULL ## for clarity
  phifitarglist["resid.model"] <- NULL ## for clarity; processed must have all relevant info
  phifitarglist["family"] <- NULL ## for clarity
  phifitarglist["rand.family"] <- NULL ## for clarity
  phifitarglist$fixed <- structure(.modify_list(phifitarglist$fixed, residProcessed$envir$ranPars ),
                                   type=.modify_list(attr(phifitarglist$fixed,"type"), attr(residProcessed$envir$ranPars,"type") ))
  phifit_init <- as.list(phifitarglist[["init"]]) ## NULL-> list(); same on next line
  lambda_source <- NULL 
  ##### (1) get info to add to phifit_init
  if (iter==0L) { 
    # Late bug fixed here (several instance below): reference to processed$residModel: this is NULL for multivariate-resp fits (instead having processed$residModel*s*)
    if (is.null(residModel$fixed[["phi"]])) { ## ie if user set fixed$phi=NA, processed in .reformat_resid_model()
      phifit_init$phi <- 1 ## A case can be made for fixing it to 1; see working doc, end of 'Gamma GLM formulation of REML estimators'
    } # and  keeps other init values as given by user.
    if ( ! is.null(processed$port_env$objective) && verbose["phifit"]) cat("\n")
    #
    # It is not sure that the following line is useful if residProcessed$port_fit_values are available, 
    #  but otherwise it makes sense:
    residProcessed$envir$inits_by_glm <- NULL ## recompute them if needed; do not use glm results on a priori distinct response
    #
    # The condition for using residProcessed$.$port_fit_values should be 
    # that the parent's processed$.$port_fit_values are available, suggesting convergence of the joint DHGLM:
    if ( ! is.null(processed$port_env$port_fit_values)) { 
      phifit_init_HLfit <- as.list(residProcessed$port_env$port_fit_value[["init_HLfit"]])
      # phifit_init$lambda default value may be NULL which means that phifit's fitme is likely to use outer optimisation
      # but it might also contain NA's. This is controlled by residModel's non-defaults
      lambda_source <- residProcessed$port_env$port_fit_values
    } else {
      residProcessed$port_env$port_fit_values <- NULL ## overrides the decision that HLfit made 'locally' within the phifit
      phifit_init_HLfit <- list()
      resid_lambda_object <- NULL
    }
    resid_corrPars <- residProcessed$port_env$port_fit_values$corrPars
  } else { ## iter>0
    #phifit_init_HLfit created at iteration 0 and stored here:
    phifit_init_HLfit <- residProcessed[["init_HLfit"]]
    if (is.null(residModel$fixed$fixef)) phifit_init_HLfit$fixef <- fixef(phifit) ## : Use an init if the parameters are not fixed. 
    phifit_init_HLfit$v_h <- phifit$v_h
    # FIXME one could *update* a corrPars[[.]]$rho element in init_HLfit... not done currently... instead
    # use ranef parameters from previous 'iter' to initiate new optimization (corrPars likely to be very slightly shifted) 
    phifit_init$corrPars <- .new_phifit_init_corrPars(phifit, innershift=1e-7) 
    lambda_source <- phifit
    if (is.null(residModel$fixed$phi)) phifit_init$phi <- phifit$phi 
    resid_corrPars <- residProcessed$port_env$port_fit_values$corrPars
  } ## end if iter==0L else ...
  # _F I X M E_ I should think about merging this with get_inits_from_fit() but the source is phifit or port_env so there must be differences.
  ###### (2) ADD elements to phifit_init and phifit_init_HLfit
  phifit_init <- .post_process_resid_corrPars(phifit_init, resid_corrPars) # synthesis with phifit_init$corrPars
  if ( ! is.null(lambda_source)) { # lambda_source is fit from previous iter, or is $port_fit_values
    # The default 'isRandomSlope' value of .get_lambdas_notrC_from_hlfit() assumes that hlfit is indeed an HLfit object
    # if lambda_source is the $port_fit_values, the 'isRandomSlope'  default value is (incorrectly) NULL.
    # we override the default with a value which is correct in all cases:
    isRandomSlope <- residProcessed$ranCoefs_blob$isRandomSlope
    phifit_init[["lambda"]] <- .get_lambdas_notrC_from_hlfit(lambda_source, type="outer",
                                                             isRandomSlope=isRandomSlope) 
    phifit_init_HLfit[["lambda"]] <- .get_lambdas_notrC_from_hlfit(lambda_source, type="inner",
                                                                   isRandomSlope=isRandomSlope)
  }
  # ###### (3) copy everything in phifitarglist
  phifitarglist[["init"]] <- phifit_init
  residProcessed[["init_HLfit"]] <- phifit_init_HLfit ## global change since residProcessed is an environment in $processed : .../...
  phifitarglist$processed <- residProcessed   ## that means that the $processed of the next HLfit call contains the phifit_init_HLfit of the last iter
  return(phifitarglist)
}

.overcat_phifit_progress <- function(phifit, verbose, processed, iter, prevmsglength) {
  la <- phifit$lambda  # not using a formal extractor, dirtily & in fear that there is any inner-estimated lambda in the phifit
  if (length(la)) names(la) <- paste0(seq_len(length(la)),".lam")
  outerEsts <- get_inits_from_fit(phifit, to_fn="fitme_body")[["init"]]
  cpla <- c(unlist(outerEsts[["corrPars"]]),la)
  cpla <- c(unlist(outerEsts[["ranCoefs"]]),cpla) # both lines so that the "corrPars" and "ranCoefs" do not get (repeated) in the display. 
  if (verbose["phifit"]) {
    prevmsglength <- overcat(
      paste0(processed$port_env$prefix,
             "phi fit's iter=",iter+1L,
             ", .phi[1]=",signif(phifit$y[1],5),", ",
             paste0(names(cpla),"=", signif(cpla,6), collapse=", "),
             ";           "),
      prevmsglength)
  }
  return(prevmsglength)
}

.sanitize_phi_est <- function(next_phi_est, control.HLfit) {
  if (any(next_phi_est<1e-12)) {
    if (is.null(min_phi <- control.HLfit$min_phi)) {
      .hack_options_error(message=paste0("Low (<1e-12) fitted residual variance (phi):\n",
                                         "this may be a genuine result for data without appropriate replicates\n",
                                         "                             and a model that allows overfitting, but\n",
                                         "(1) this may also point to problems in the data (duplicated response values?);\n",
                                         "(2) this may have led to the present error from a later computation.\n",
                                         "You may overcome this by setting control.HLfit$min_phi\n",
                                         "    to 1e-10 or some other low, but not too low, value.\n",
                                         "Still, the computed likelihood maximum wrt all parameters may be inaccurate.\n"))
      # crash later computations: ::Cholesky(wd2hdv2w) problem
    } else next_phi_est[next_phi_est<min_phi] <- min_phi  
  } else .hack_options_error(message=NULL)
  return(next_phi_est)
}

.mu_U2fv <- function(mu, BinomialDen, muetablob, processed, family=processed$family,
                         families=processed$families, data=processed$data) {
  if (is.null(families)) {
    if (family$family=="binomial") { 
      fv <- mu/BinomialDen ## cf glm(binomial): fitted values are frequencies 
    } else if ( ! is.null(muetablob$p0)) {
      fv <- structure(mu/(1-muetablob$p0), ## mu Truncated from mu Untruncated.
                      mu_U=mu,
                      p0=muetablob$p0) ## more efficient info for simulation.
    } else {fv <- mu} ## fitted values may be counts (cf poisson), or reals
    names(fv) <- rownames(data) ## in some old context it got names given by .ULI(), with repeats... 
  } else {
    cum_nobs <- attr(families,"cum_nobs")
    fvs <- vector("list",length(families)) # clik_fn returns a single value when it uses family()$aic  (e.g., poisson) (so summand=FALSE fails) and a vector otherwise
    for (mv_it in seq_along(fvs)) {
      resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
      fvs[[mv_it]] <- .mu_U2fv(mu=mu[resp_range], BinomialDen=BinomialDen[resp_range], 
                                   muetablob=muetablob$mv[[mv_it]], family=families[[mv_it]],
                                   families=NULL, data=processed$unmerged[[mv_it]]$data) 
    }
    fv <- structure(unlist(fvs, recursive=FALSE, use.names=TRUE), mv=fvs)
  }
  return(fv)
}

.dfs_ranCoefs <- function(types) {
  # $ranCoefs is present only for (partially) fixed ranCoefs...
  # its type is present only for *partially* fixed  ranCoefs!
  # => wonderful ad-hockery.
  #p_trRanCoefs <- length(which(unlist(attr(res$CorrEst_and_RanFix,"type")$trRanCoefs, use.names = FALSE) == "outer")) 
  #p_partiallyfixed_ranCoefs <- length(which(unlist(attr(res$CorrEst_and_RanFix,"type")$ranCoefs, use.names = FALSE) == "fix")) 
  #
  # Even better with isDiagFamily, where ranCoef has more elements than trRanCoef... hence:
  p_rC <- 0L
  ty_trRanCoefs <- types$trRanCoefs
  ty_ranCoefs <- types$ranCoefs
  for (char_rd in names(ty_trRanCoefs)) {
    if (is.null(ty_ranCoefs[[char_rd]])) { # not partially fixed
      p_rC <- p_rC + length(which(ty_trRanCoefs[[char_rd]]=="outer"))
    } else { # partially fixed: trRanCoefs cannot be used bc all types are "outer" and their number differ whether isDiagFamily or not
      p_rC <- p_rC + length(which(ty_ranCoefs[[char_rd]]=="outer"))
    }
  }
  p_rC
}

.get_MME_method <- function(auglinmodblob, HL) {
  if (HL[1]=="SEM") {
    MME_method <- "stochastic EM"
  } else {
    if (is.null(auglinmodblob)) { # GLM with fixed fixef
      MME_method <- "no model matrix"
    } else {
      get_from <- attr(auglinmodblob$sXaug,"get_from")
      if (is.null(get_from)) { # Not sure this still happens: is next comment obsolete? 
        MME_method <- setdiff(class(auglinmodblob$sXaug),c("list")) # "dgCMatrix" for _Matrix_QRP_CHM bc the class of this S4 object cannot be modified... hence the alternative code  
      } else {
        MME_method <- unique(c(substring(get_from,14), setdiff(class(auglinmodblob$sXaug),c("list"))))
      }
    }
  }
  MME_method
}

.post_process_warningList <- function(warningList, processed, maxit.mean,  pforpv, innerj, 
                                      nonSPD,
                                      conv_logL, iter, conv.lambda, conv.phi, conv.corr,
                                      #
                                      HL=processed$HL, models=processed$models, LMMbool=attr(models,"LMMbool"), max.iter=processed$max.iter
                                      ) {
  ## translation of warnings in user-more friendly form 
  if (identical(nonSPD, TRUE)) {
    warningList$regularizedHess <- "The negative-Hessian was not positive definite at the final estimates. Likelihood may not be maximized."
  }
  if ( ! is.null(warningList$anyVanishLam) && warningList$anyVanishLam) {
    warningList$anyVanishLam <- "Some lambda estimates may actually be zero (lower bound controlled by option 'minLambda')"
  }
  if (! is.null(warningList$leveLam1) && warningList$leveLam1) {
    warningList$leveLam1 <- paste("lambda leverages numerically 1 were replaced by 1-",
                                  .spaMM.data$options$regul_lev_lambda,"(as controlled by option 'regul_lev_lambda')")
  }
  if (! is.null(locw <- warningList$innerPhiGLM)) {
    warningList$innerPhiGLM <- paste0("'",locw,"' in some sub-final iteration(s) of phi estimation;")
  }
  if ( HL[1L]!="SEM" && maxit.mean>1 
       && ( models[[1L]]=="etaHGLM" || pforpv>0L) ## cases where iterations are needed : including pforpv>0L for GLM
       && innerj==maxit.mean ) {
    warningList$innerNotConv <- paste0("linear predictor estimation did not converge;",
                                       if ( models[[1L]]=="etaHGLM" && ! LMMbool && ! processed$LevenbergM["LM_start"] ) {
                                         " try control.HLfit=list(LevenbergM=TRUE), or"
                                       },
                                       " increase 'control.HLfit$max.iter.mean' above ",maxit.mean)
  }
  if ( (! is.na(conv_logL)) && iter==max.iter) {
    maxitmess <- paste0("Estimates did not converge;",
                        if ( ! LMMbool && ! processed$LevenbergM["LM_start"] ) {
                          " try control.HLfit=list(LevenbergM=TRUE), or"
                        },
                        " increase 'max.iter' above ",max.iter,
                        "\n (see help('HLfit') for details about 'max.iter')")
    if (models[["eta"]]=="etaHGLM") {
      if (conv_logL  && ! conv.lambda) {
        mainNotConv <- paste0("p_v apparently converged but lambda estimates apparently did not.",
                              "\n This may indicate that some lambda estimate(s) should be zero.",
                              "\n Otherwise try increasing 'max.iter' above ",max.iter,
                              "\n (see help(HLfit) for details about 'max.iter')")          
      } else mainNotConv <- maxitmess        
    } else mainNotConv <- maxitmess        
    attr(mainNotConv,"diagnostics") <- c( conv.phi=conv.phi , conv.lambda=conv.lambda, 
                                          conv.corr=conv.corr )
    warningList$mainNotConv <- mainNotConv
  }
  if ( ((! is.null(warningList$innerNotConv))||(! is.null(warningList$mainNotConv))) && 
       ! is.null(locw <- processed$envir$inits_by_glm$conv_info)) {
    warningList$inits_by_glm <- locw ## added 2019/04/03
  }
  warningList
}

.verbose_warnings <- function(verbose, warningList) {
  if (verbose["all_objfn_calls"]) {
    seriousWarnings <- warningList[intersect(c("innerNotConv","mainNotConv"),names(warningList))]
    if (length(seriousWarnings) ) { 
      warningsss <- paste0("In HLfit :\n",unlist(seriousWarnings))
      abyss <- sapply(warningsss, warning, call.=FALSE) 
      warningList[setdiff(names(warningList),c("innerNotConv","mainNotCov"))] <- NULL
    }
  } 
  if (verbose["trace"]) {
    if (length(warningList) ) {
      warningsss <- paste0(unlist(warningList),"\n")
      abyss <- sapply(warningsss,cat) 
    }
  }
  # no used return value
}

.add_phi_returns <- function(res, processed, PHIblob, phi.Fix, dev_res, phi_est,
                             #
                             models=processed$models, vec_nobs=processed$vec_nobs) {
  if ( ! is.null(vec_nobs) ) {
    res$residModels <- vector("list", length(vec_nobs))
    cum_nobs <- attr(res$families,"cum_nobs")
    for (mv_it in seq_along(vec_nobs)) {
      res$residModels[[mv_it]] <- list(formula=processed$residModels[[mv_it]]$formula,
                                       family=attr(processed$residModels[[mv_it]]$family,"quoted"))
      if (is.null(res[["resid_fits"]])) res[["resid_fits"]] <- vector("list", length(vec_nobs))
      if (models[["phi"]][mv_it]=="phiHGLM") {
        res[["resid_fits"]][[mv_it]] <- PHIblob$multiPHI[[mv_it]]$phifit
      } else {
        resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
        res[["phi.object"]][[mv_it]] <- .get_phi_object(phi.Fix[[mv_it]], PHIblob=PHIblob$multiPHI[[mv_it]], dev_res[resp_range], 
                                                        prior.weights=res$prior.weights[[mv_it]], 
                                                        phi.preFix=processed$phi.Fix[[mv_it]], nobs=vec_nobs[mv_it], 
                                                        control=processed[["control.glm"]])
      }
    }
    res[["phi"]] <- phi_est 
  } else {
    res$residModel <- list(formula=processed$residModel$formula,
                           family=attr(processed$residModel$family,"quoted"))
    if (models[["phi"]]=="phiHGLM") {
      res[["resid_fit"]] <- PHIblob$phifit
      res[["phi"]] <- phi_est 
    } else {
      res[["phi.object"]] <- .get_phi_object(phi.Fix, PHIblob, dev_res, prior.weights=res$prior.weights, 
                                             phi.preFix=processed$phi.Fix, control=processed[["control.glm"]])
      # phi_est comes from PHIblob$next_phi_est, not from final glm,  hence is in minimal form
      if (models[["phi"]]=="phiScal") {res[["phi"]] <- phi_est[1L]} else res[["phi"]] <- phi_est 
    }
  }
  res
}

.broom_corrMatrices <- function(sub_corr_info, processed) { # the *aim* is to keep here only full corrMatrices with more levels than levels that in ZA. 
  for (rd in seq_along(sub_corr_info$corrMatrices)) {
    if ( ! is.null(fullcorrmat <- sub_corr_info$corrMatrices[[rd]]) && # ___F I X M E___ test surely imperfect
         ncol(processed$ZAlist[[rd]]) == ncol(fullcorrmat)){ sub_corr_info$corrMatrices[rd] <- list(NULL)
    } else if (inherits(fullcorrmat,"dist")) {
      sub_corr_info$corrMatrices[[rd]] <- proxy::as.matrix(fullcorrmat, diag=1) 
    } 
  }
  sub_corr_info
}

.add_ranef_returns <- function(res, processed, wranefblob, 
                               lambda_est, process_resglm_blob, LMatrices, init.lambda, v_h, u_h, ranCoefs_blob,
                               #
                               models=processed$models, nrand=length(processed$ZAlist), cum_n_u_h=processed$cum_n_u_h) {
  res$ZAlist <- processed$ZAlist ## needed for prediction variance
  #
  sub_corr_info <- mget(c("corr_families","corr_types", "AMatrices", "corrMatrices"),  processed$corr_info)
  sub_corr_info <- .broom_corrMatrices(sub_corr_info, processed)
  kron_Y_LMatrices <- vector("list", nrand)
  for (it in seq_len(nrand)) kron_Y_LMatrices[it] <- list(attr(processed$corr_info$cov_info_mats[[it]],"blob")$Lunique)
  sub_corr_info$kron_Y_LMatrices <- kron_Y_LMatrices
  # 
  res$ranef_info <- list(sub_corr_info=sub_corr_info, hyper_info=processed$hyper_info,
                         vec_normIMRF=processed$AUGI0_ZX$vec_normIMRF)
  res$w.ranef <- wranefblob$w.ranef ## useful for .get_info_crits() and get_LSmatrix()
  #
  res$lambda.object <- .make_lambda_object(nrand, lambda_models=models[["lambda"]], cum_n_u_h, lambda_est, 
                                           process_resglm_blob, rand.families=processed$rand.families,
                                           ZAlist=processed$ZAlist, LMatrices, lambdaType=attr(init.lambda,"type"))
  
  strucBlob <- .post_process_v_h_LMatrices(next_LMatrices=LMatrices, v_h=v_h, u_h=u_h,
                                           processed=processed, ranCoefs_blob=ranCoefs_blob) 
  res$strucList <- strucBlob$strucList
  res$v_h <- strucBlob$v_h
  # $lambda is a vector that NO LONGER contains lambdas for ranCoefs (too confusing, the more so as they are often 1)
  res[["lambda"]] <- structure(.get_lambdas_notrC_from_hlfit(res, type="adhoc"), cum_n_u_h=cum_n_u_h) ## redundant but very convenient except for programming
  res
}

.wrap_IRLS <- function(nrand, intervalInfo, processed, beta_eta,
                       ZAL, y, lambda_est, muetablob, maxit.mean, etaFix, 
                       wranefblob, u_h, v_h, w.resid, H_w.resid, phi_est, pforpv, verbose,
                       ad_hoc_corrPars) {
  if (nrand) { # (models[["eta"]]=="etaHGLM") {
    auglinmodblob <- processed$solve_IRLS_fn(ZAL=ZAL, y=y, 
                                 lambda_est=lambda_est, # homoscedastic for CAR! confusing!
                                 muetablob=muetablob,
                                 maxit.mean=maxit.mean, etaFix=etaFix,
                                 ## supplement for LevenbergM
                                 beta_eta=beta_eta,
                                 ## supplement for ! GLMM
                                 wranefblob=wranefblob, u_h=u_h, v_h=v_h, w.resid=w.resid, H_w.resid=H_w.resid,
                                 phi_est=phi_est,
                                 # for_init_z_args=list(nrand=nrand, psi_M=processed$psi_M), 
                                 for_intervals=intervalInfo,
                                 ##
                                 processed=processed,
                                 corrPars=ad_hoc_corrPars # ignored by .solve_IRLS_as_ZX
    )
  } else if (pforpv>0L  && maxit.mean) {
    ## resultat nomme' auglinmodblob pour avancer vers simplif du code, mais je ne peux suppser qu'auglinmodblob est toujours calculé
    auglinmodblob <- processed$etaxLM_fn( # .calc_etaLLMblob or .calc_etaGLMblob
      processed=processed,  muetablob=muetablob, old_beta_eta=beta_eta, 
      w.resid=w.resid, # ignored when .calc_etaxLMblob = .calc_etaLLMblob
      phi_est=phi_est, maxit.mean=maxit.mean, 
      verbose=verbose, for_intervals=intervalInfo)
  } else auglinmodblob <- NULL 
  auglinmodblob
}

.eval_conv.phi <- function(phi_est, next_phi_est, spaMM_tol) {
  if (is.list(phi_est)) {
    conv.phi <- TRUE
    for (mv_it in seq_along(phi_est)) {
      conv.phi <- all(abs(next_phi_est[[mv_it]]-phi_est[[mv_it]]) < 
                        spaMM_tol$Xtol_rel*phi_est[[mv_it]] + spaMM_tol$Xtol_abs) 
      if ( ! conv.phi) break
    }
  } else conv.phi <- all(abs(next_phi_est-phi_est) < spaMM_tol$Xtol_rel*phi_est + spaMM_tol$Xtol_abs) ## 'weak convergence'... 
  conv.phi
}
