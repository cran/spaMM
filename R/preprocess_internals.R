.nonzeros <- function(spm) {
  if (inherits(spm,"ddiMatrix") && spm@diag=="U") {
    return(ncol(spm))
  } else if (inherits(spm,"sparseMatrix")) {
    nz <- length(spm@x)
    if ( methods::.hasSlot(spm, "diag")) nz <- nz+ncol(spm)
    return(nz) 
  } else return(sum(spm !=0)) ## AMatrices reaches here
}

.provide_G_diagnosis <- function(corr_info, ZAlist,thr=0.05, fast=TRUE) {
  # .assign_geoinfo_and_LMatrices_but_ranCoefs() may serve as a template, but we don't want actual matrices except to assess computation costs
  if (is.null(corr_info$G_diagnosis)) {
    adjlist <- corr_info$adjMatrices
    corrlist <- corr_info$corrMatrices
    exp_ranef_types <- attr(ZAlist, "exp_ranef_types")
    which_nested <- grep("%in%", names(attr(ZAlist,"exp_ranef_terms")))
    for (rd in seq_along(exp_ranef_types)) {
      if (exp_ranef_types[rd]==c("corrMatrix") ) {
        ## cov_info_mats elements may be correlation matrices, or they may be lists...
        if (inherits(corrlist[[rd]],"dist")) {
          corrlist[[rd]] <- proxy::as.matrix(corrlist[[rd]], diag=1)
        } else if (is.list(corrlist[[rd]])) {
          corrlist[rd] <- list(NULL) ## hmmm that is a quick patch... (F I X M E ?) but this suggests we aim at spprec
        } 
        if (is.matrix(corrlist[[rd]])) {
          ## using the corrMatrix rather than a chol factor overestimates the computational weight, 
          # and the chol factor of a relatedness matrix appears as sparse as its precmat. So its worth evaluating it.
          ## Example of that is Gryphon, better by correlation algo
          corrlist[[rd]] <- mat_sqrt(corrlist[[rd]]) # as used for Lunique which is tcrossfac (tcrossprod(mat_sqrt(X))=X)
          ZAlist[[rd]] <- .addrightcols_Z(Z=ZAlist[[rd]], colnames(corrlist[[rd]]))
        }
      } else if (exp_ranef_types[rd]%in% c("Matern","Cauchy")  ) {
        nc <- ncol(ZAlist[[rd]])
        adjlist[[rd]] <- matrix(TRUE,ncol=nc,nrow=nc) # as(allTRUE,"lgCMatrix") #new("lgCMatrix",i=rep(c(0L,seq_len(nc-1L)),nc),p=c(0L,seq_len(nc))*nc,x=rep(TRUE,nc^2),Dim=c(nc,nc)) 
        corrlist[[rd]] <- lower.tri(adjlist[[rd]],diag = TRUE) # logi
      } 
    }
    # suppressMessages here and below as .provide_G_diagnosis() is not the right context for messages.
    noAR <- suppressMessages( .compute_ZAL(XMatrix=corrlist,ZAlist,as_matrix = FALSE) )# without the cost of the sparse structures
    #
    for (rd in seq_along(exp_ranef_types)) {
      if (exp_ranef_types[rd] %in% c("AR1", "IMRF") ) { # for "adjacency" we assume adjlist[[rd]] is filled
        if (is.null(adjlist[[rd]])) {
          nc <- ncol(ZAlist[[rd]])
          locmat <- diag(nrow=nc)
          if ( ! (rd %in% which_nested)) {
            diag(locmat[-1,]) <- 1 ## for nested AR1 this overestimates the computational cost.
          } else {
            diag(locmat[-1,]) <- 1 
          }
          adjlist[[rd]] <- (locmat+t(locmat))/2
        }
      } else if (exp_ranef_types[rd] == c("adjacency")) { # adjlist[[rd]] presumably has 0s on diagonal
        rowmax <- max(rowSums(adjlist[[rd]]))
        adjlist[[rd]] <- adjlist[[rd]] + Diagonal(n=ncol(adjlist[[rd]]),x=rowmax+1) # make it diagonally dominant
      }
      if ( ! is.null(corrlist[[rd]])) {
        if ( is.null(adjlist[[rd]])) adjlist[[rd]] <- chol2inv(chol(corrlist[[rd]]))
      } else if ( ! is.null(adjlist[[rd]])) {
        if ( is.null(corrlist[[rd]])) {
          if (exp_ranef_types[rd] %in% c("AR1", "IMRF","adjacency")) {
            # much better in spprec: adjacency-long and some in test AR1 (long):
            # fitar1 <- corrHLfit(obs ~ 1+AR1(1|age),family=poisson(),data=fake,verbose=c(TRACE=TRUE)) 
            # fit_ar1nested <- ... also
            # If the user provided a huge adjmatrix it's risky to compute a huge dense corrMatrix... 
            # correlation algos can still be selected and are better (ohio small and many scotlip tests)
            # spprec and ! roughly as fast for fitNF in test-devel-predVar-AR1; which has crit=43
            nc <- ncol(adjlist[[rd]])
            corrlist[[rd]] <- lower.tri(matrix(TRUE,ncol=nc,nrow=nc),diag = TRUE) # template matrix created in faster way than diag()
          } else { # not clear when this can occur
            tcrossfac_adj <- mat_sqrt(adjlist[[rd]]) # tcrossfac hence t(solve()) is the tcrossfac of the corr mat (<=> Lunique) which is the following backsolve
            if ( attr(tcrossfac_adj,"type")=="cholL_LLt") {
              corrlist[[rd]] <- .backsolve(tcrossfac_adj,upper.tri = FALSE, transpose=TRUE) 
            } else corrlist[[rd]] <- t(solve(tcrossfac_adj)) # quick patch when another facto used.  
          }
        } 
      } else {# "(.|.)" 
        ## tnb <- fitme(resp~1+(1|ID), data=lll,family=Tnegbin(2)) is a test case where spprec is clearly slower
        adjlist[[rd]] <- .symDiagonal(TRUE,n=ncol(ZAlist[[rd]])) # .symDiagonal(TRUE,n=ncol(ZAlist[[rd]]))  ## (null corrlist[[rd]] must mean the same thing) 
      }
    }
    ### > qq s for large ZA 
    ZL <- suppressMessages( .compute_ZAL(XMatrix=corrlist,ZAlist,as_matrix = FALSE) )  
    if (is.logical(ZL)) {# logical corr by identity Z gives logi ZL, not handled by .crossprod
      cross_ZL <- crossprod(ZL)
    } else cross_ZL <- .crossprod(ZL) ## forces a call to forceSymmetric => result is Matrix either dsy or sparse.
    denseness_via_ZL <- .calc_denseness(cross_ZL) # crossprod ideally dsC except if ZL is really dense
    crossZL_is_dsy <- inherits(cross_ZL,"dsyMatrix") 
    ###
    if (fast) {
      corr_info$G_diagnosis <- list(denseness_noAR=.calc_denseness(noAR), 
                                    crossZL_is_dsy=crossZL_is_dsy,
                                    denseness_via_ZL=denseness_via_ZL)
      # if there are only "(.|.)", we compare noAR to crossprod(ZL) =crossprod(noAR) and the latter may be less dense!
    } else {
      for (rd in seq_along(exp_ranef_types)) adjlist[[rd]] <- as(forceSymmetric(adjlist[[rd]]),"dsCMatrix")
      locQ <- do.call(Matrix::bdiag, adjlist) # dsC
      Z_ <- suppressMessages( .compute_ZAL(XMatrix=NULL,ZAlist,as_matrix = FALSE) )
      locG <- .crossprod(Z_)+locQ # ideally dsC except if Z_ is really dense
      denseness_G <- .calc_denseness(locG)
      dens_G_rel_ZL <- denseness_G/denseness_via_ZL
      density_G <- denseness_G/prod(dim(locG))
      crit <- (dens_G_rel_ZL<1 && density_G*dens_G_rel_ZL<thr)
      #
      corr_info$G_diagnosis <- list(denseness_noAR=.calc_denseness(noAR), 
                                    crossZL_is_dsy=crossZL_is_dsy,
                                    denseness_via_ZL=denseness_via_ZL, 
                                    # supplements experimentaux
                                    dens_G_rel_ZL=dens_G_rel_ZL, density_G=density_G,
                                    crit=crit)
    }
  }
  return(corr_info$G_diagnosis)
}

# even though the Z's were sparse postmultplication by LMatrix leads some of the ZAL's to dgeMatrix (dense)
.choose_QRmethod <- function(ZAlist, predictor, corr_info, trySparse=TRUE,is_spprec) {
  if ( is.null(QRmethod <- .spaMM.data$options$QRmethod) ) { ## user setting. The code should NOT write into it. 
    nrand <- length(ZAlist)
    if (trySparse && nrand>0L) {
      # adjacency speed to be tested on 2nd example from test-spaMM.R
      densecorrs <- attr(ZAlist,"exp_ranef_types") %in% c("adjacency", "IMRF", "Matern","Cauchy", "corrMatrix","AR1") 
      sparseprecs <- attr(ZAlist,"exp_ranef_types") %in% c("adjacency", "IMRF", "AR1")
      if (is_spprec && all(sparseprecs)) {
        totdim <- colSums(do.call(rbind,lapply(ZAlist,dim)))
        if (totdim[2L]>1000L) { # a bit a hoc (ohio/adjacency-long/large IMRF)
          return("sparse")
        } else return("dense")
      } else if (( ! is_spprec) && all(densecorrs) ) { ## simple subcase of the next case
        ## LMatrices are not available, and it may be better to use the density of the correlation matrices anyway:
        ## for maximally sparse Z, ZL's denseness is that of the retained rows of L. This suggests that ZL could be made sparser 
        ## by reordering the levels of the correlation matrix so that the most represented levsl come first in a triangular L factor. 
        ## But this would not affect the denseness of .crossprod(ZW) in .get_absdiagR_blocks(), 
        ## and this leads to use "dense" whenever the correlation matrix is dense.
        ## Gryphon example is useful here, L is sparse but the "dense" QRmethod still appears as good as the "sparse" one.
        return("dense")
      } else if (any(densecorrs)) {
        G_diagnosis <- .provide_G_diagnosis(corr_info=corr_info, ZAlist=ZAlist)
        if (G_diagnosis$crossZL_is_dsy) { ## sufficient, but loose, condition for using dense
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

.preprocess_LevM <- function(user_LM, processed, nrand) {
  if (processed$LMMbool) user_LM <- FALSE  # (_F I X M E_) removing this and forcing LevM creates an error in the tests.
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
          LM_start <- (tail(processed$cum_n_u_h,n=1)<500L) ## adjlg has 1000 levels and is faster without LevM 
        } else LM_start <- FALSE  ## BINARYboot test to assess effect on timings
      } else LM_start <- FALSE
    } else LM_start <- user_LM[[1L]] # important to drop name else the names of the vector are wrong
    return(c(user_LM=user_LM[[1L]], LM_start=LM_start) ) 
  } else return(user_LM) # allows a full vector to be user-provided, for full control
}

.check_subset_corrMatrix <- function(corrMatrix,ZA) {
  ZAnames <- colnames(ZA) ## set by .calc_Zlist() or .calc_ZAlist(), with two cases for corrMatrix 
  if (is.null(ZAnames)) {
    stop("NULL colnames in (a block of) the design matrix for random effects. Some mishandling of 'AMatrices'?")
  }
  if (inherits(corrMatrix,"dist")) {
    corrnames <- labels(corrMatrix) ## unclear
  } else if (inherits(corrMatrix,c("matrix","Matrix"))) {
    corrnames <- rownames(corrMatrix)
  } else if ( inherits(corrMatrix,"precision")) {
    corrnames <- rownames(corrMatrix[["matrix"]])
  } else stop("Unhandled class of corrMatrix object.")
  if (is.null(corrnames)) {
    mess <- paste("(!) corrMatrix without labels or row names: the grouping levels, in order",
                  paste0(ZAnames[1L:min(5L,length(ZAnames))], collapse=" "),if(length(ZAnames)>5L){"...,"} else{","},
                  "\n are matched in this order to rows and columns of corrMatrix, without further check.",
                  "\n This may cause later visible errors (notably, wrongly dimensioned matrices)",
                  "\n or even silent errors. See help(\"corrMatrix\") for a safer syntax.")
    warning(mess, immediate. = TRUE)
  } else if (is.null(colnames(corrMatrix))) {
    if (inherits(corrMatrix, c("matrix", "Matrix"))) {
      colnames(corrMatrix) <- corrnames
    }
    else if (inherits(corrMatrix, "precision")) {
      colnames(corrMatrix[["matrix"]]) <- corrnames
    }
  }
  if ( inherits(corrMatrix,"precision")) { 
    # do not subset a precmat => nothing here but the corresponding Z matrix may be modified by .preprocess() 
  } else if ( length(setdiff(ZAnames,corrnames)) ==0L ) { ## i.e. all ZAnames in corrnames
    ## : should be the case when generator = "as.factor"
    if ( length(setdiff(corrnames,ZAnames)) || any(corrnames!=ZAnames) ) { # reordering and subsetting
      if (inherits(corrMatrix,"dist")) {
        corrMatrix <- (proxy::as.matrix(corrMatrix,diag=1)[ZAnames,ZAnames]) ## IF diag missing in input corrMatrix THEN assume a correlation matrix
        ## it's not useful to convert back to dist (either uglily by as.dist(); or package 'seriation' has (permute.dist-> C code)
      } else corrMatrix <- corrMatrix[ZAnames,ZAnames]  
    } ## else orders already match
  } else {
    if ( ! is.null(corrnames)) {
      if ( length(corrnames)!=length(ZAnames)){ 
        stop("The dimension of corrMatrix does not match the number of levels of the grouping variable.")
      } else { ## no clear reordering
        message(paste0("spaMM is not able to match levels of the random effect to the names of corrMatrix,\n",
                       " and matches levels to rows of the matrix by their respective orders.\n",
                       " See help(\"corrMatrix\") for a safer syntax."))
      }
    }
  }
  return(corrMatrix)
}

.addrightcols_Z <- function(Z, precnames) {
  ZAnames <- colnames(Z)
  ncol_prec <-  length(precnames)
  ncol_Z <- ncol(Z)
  # We have tested in .check_subset_corrMatrix() whether all ZAnames were in precnames 
  # so the only possible difference between sets of names is additional names in precnames
  if (suppcols <- ncol_prec-ncol_Z) { # (names=levels) in precmat but not in the data
    message(paste("Note: Precision matrix has", suppcols, "more levels than there are in the data.")) # and <0 vaues are a bug...
    supplevels <- setdiff(precnames,ZAnames)
    # add cols of zeros one the right; cols to be reoredered next.
    if ( inherits(Z,"ddiMatrix")) Z <- .as_ddi_dgC(Z)
    Zp <- Z@p
    Z@p <- c(Zp, Zp[length(Zp)] + rep(0L,suppcols))
    Z@Dim[2L] <- ncol_prec
    Z@Dimnames[[2L]] <- c(ZAnames,supplevels) 
  } 
  Z
}

.preprocess_pw <- function(subs_p_weights, nobs, validdata) {
  if (is.null(subs_p_weights)) {
    prior.weights <- structure(rep(1L,nobs),unique=TRUE,is_call=FALSE) ## <- 1L prevented by glm -> model.frame(... prior.weights)
  } else if ( ! (inherits(subs_p_weights,"call") && subs_p_weights[[1L]] == "quote") )  {
    prior.weights <- as.vector(stats::model.weights(validdata)) ## as.vector as in say lm() protects against array1d
    if ( ! is.numeric(prior.weights)) 
      stop("'weights' must be a numeric vector")
    if (any(prior.weights < 0)) 
      stop("negative weights not allowed")
    prior.weights <- structure(prior.weights, unique= length(unique(prior.weights))==1L, is_call=FALSE)
    #attr(prior.weights,"only1") <- all(upw==1L)
  } else {## 'prior.weights' is a quoted expression
    attr(prior.weights,"unique") <- FALSE 
    attr(prior.weights,"is_call") <- TRUE
  }   
  return(prior.weights)
}


