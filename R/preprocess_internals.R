.nonzeros <- function(spm) {
  if (inherits(spm,"ddiMatrix") && spm@diag=="U") {
    return(ncol(spm))
  } else if (inherits(spm,"sparseMatrix")) {
    nz <- length(spm@x)
    if ( methods::.hasSlot(spm, "diag")) nz <- nz+ncol(spm)
    return(nz) 
  } else return(sum(spm !=0)) ## AMatrices reaches here
}

.provide_G_diagnosis <- function(corr_info, ZAlist) {
  if (is.null(corr_info$G_diagnosis)) {
    loccorrlist <- corr_info$cov_info_mats ## to handle corrMatrix
    exp_ranef_types <- attr(ZAlist, "exp_ranef_types")
    for (rd in which(exp_ranef_types==c("corrMatrix") )) {
      ## cov_info_mats elements may be correlation matrices, or they may be lists...
      if (inherits(loccorrlist[[rd]],"dist")) {
        loccorrlist[[rd]] <- proxy::as.matrix(loccorrlist[[rd]], diag=1)
      } else if (is.list(loccorrlist[[rd]])) {
        loccorrlist[rd] <- list(NULL) ## hmmm that is a quick patch... (F I X M E ?) but this suggests we aim at spprec
        prosparse <- TRUE
      }
    }
    for (rd in which(exp_ranef_types %in% c("Matern","Cauchy"))) {
      nc <- ncol(ZAlist[[rd]])
      locmat <- matrix(0,ncol=nc,nrow=nc)
      locmat[lower.tri(locmat,diag=TRUE)] <- 1
      loccorrlist[[rd]] <- locmat
    }
    # Accounting only for the cost of the non-sparse structures (and of "(.|.)"):
    noAR <- .compute_ZAL(XMatrix=loccorrlist,ZAlist,as_matrix = FALSE)
    cross_noAR <- .crossprod(noAR) ## forces a call to forceSymmetric => result is Matrix either dsy or sparse.
    if (any(exp_ranef_types %in% c("adjacency", "IMRF", "AR1"))) { # then we need to evaluate the cost of the sparse structures 
      which_nested <- grep("%in%", names(attr(ZAlist,"exp_ranef_terms")))
      for (rd in which(exp_ranef_types %in% c("adjacency", "IMRF", "AR1"))) {
        nc <- ncol(ZAlist[[rd]])
        if ( ! (rd %in% which_nested)) {
          locmat <- matrix(0,ncol=nc,nrow=nc)
          locmat[lower.tri(locmat,diag=TRUE)] <- 1 ## for nested AR1 his overestimates the computational cost.
        } else {
          locmat <- diag(nrow=nc)
          diag(locmat[-1,]) <- 1 # heuristic compromise
        }
        loccorrlist[[rd]] <- locmat
      }
      # 
      ZC <- .compute_ZAL(XMatrix=loccorrlist,ZAlist,as_matrix = FALSE)
      cross_ZC <- .crossprod(ZC) ## forces a call to forceSymmetric => result is Matrix either dsy or sparse.
      crossZC_is_dsy <- inherits(cross_ZC,"dsyMatrix") 
      len_crossZC <- length(cross_ZC@x)
    } else {
      crossZC_is_dsy <- inherits(cross_noAR,"dsyMatrix") 
      len_crossZC <- length(cross_noAR@x)
    } 
    corr_info$G_diagnosis <- list(len_crossZC=len_crossZC, len_noAR=length(cross_noAR@x), crossZC_is_dsy=crossZC_is_dsy)
  }
  return(corr_info$G_diagnosis)
}

# even though the Z's were sparse postmultplication by LMatrix leads some of the ZAL's to dgeMatrix (dense)
.choose_QRmethod <- function(ZAlist, predictor, corr_info, trySparse=TRUE) {
  if ( is.null(QRmethod <- .spaMM.data$options$QRmethod) ) { ## user setting. The code should NOT write into it. 
    nrand <- length(ZAlist)
    if (trySparse && nrand>0L) {
      # adjacency speed to be tested on 2nd example from test-spaMM.R
      densecorrs <- attr(ZAlist,"exp_ranef_types") %in% c("adjacency", "IMRF", "Matern","Cauchy", "corrMatrix") 
      if (all(densecorrs)) { ## simple subcase of the next case
        ## LMatrices are not available, and it may be better to use the density of the correlation matrices anyway:
        ## for maximally sparse Z, ZL's denseness is that of the retained rows of L. This suggests tht ZL could be made sparser 
        ## by reordering the levels of the correlation matrix so that the most represented levsl come first in a triangular L factor. 
        ## But this would not affect the denseness of .crossprod(ZW) in .get_absdiagR_blocks(), 
        ## and this leads to use "dense" whenever the correlation matrix is dense.
        return("dense")
      } else if (any(densecorrs)) {
        G_diagnosis <- .provide_G_diagnosis(corr_info=corr_info, ZAlist=ZAlist)
        if (G_diagnosis$crossZC_is_dsy) { ## sufficient, but loose, condition for using dense
          return("dense")
        } else {
          totdim <- colSums(do.call(rbind,lapply(ZAlist,dim)))
          if (totdim[1L]>4L*totdim[2L]) {
            return("sparse")
          } else return("dense")
        }
      } else if (nrand==1L && .is_identity(ZAlist[[1]])) { ## test pertinent slmt pour non-spatial models !
        return("sparse") ## special case for poisson or binomial with saturated ranef
      } else { ## several block effects...
        # could use .crossprod() here too to assess sparsity.
        totdim <- colSums(do.call(rbind,lapply(ZAlist,dim)))
        totsize <- prod(totdim)
        nonzeros <- sum(unlist(lapply(ZAlist, .nonzeros)))          
        if (nonzeros/totsize < .spaMM.data$options$sparsity_threshold) { 
          return("sparse")
        } else {
          return("dense") ## ZAlist actually not so sparse
        }
      }
    } else return("dense") 
  }
  return(QRmethod)
}

.preprocess_LevM <- function(user_LM, processed) {
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
        if (processed$HL[1L]==0L) {
          ## PQL/L + LevenbergM combine safely and relatively fast.... for small data
          # bigranefs -> PQL/L+LevM much smaller than ML!
          LM_start <- (tail(processed$cum_n_u_h,n=1)<500L) ## adjlg has 1000 levels and is faster without LevM 
        } else LM_start <- FALSE  ## BINARYboot test to assess effect on timings
      } else LM_start <- FALSE
    } else LM_start <- user_LM[[1L]] # important to drop name else the names of the vector are wrong
    return(c(user_LM=user_LM[[1L]], LM_start=LM_start) ) 
  } else return(user_LM) # allows a full vector to be user-provided, for full control
}

