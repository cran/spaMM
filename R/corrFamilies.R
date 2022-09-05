# The parvec are the partial autocorrelation coefficients  corr(X_s,X_s+t | X_s+1,..., X_s+t-1), in the hypercube (-1,1)^p,
# They are converted to the phi parameters of the autoregression (Barndorff-Nielsen & Schou 1973 eq.10 and various papers cited in Jones 1987...)
# Then the inverse covariance matrix is built from the phi coeffs (Verbyla 1985; seems to have Verbyla seems to have mixed i and j in the G_i,i+j terms)
# Finally this matrix is normalized as an inverse correlation matrix by correcting by the stationary variance, given by B-N & S eq 7 (with unit innovation variance)
# Possibly nonzero elements of the matrix are built in the 'xvec'  vector. Direct operations on the matrix in dense format would look like:
# Qmat <- diag(nc)
# for (jt in seq_len(jmax2)) {
#   sub_diagPos <- diagPos[jt+seq(nc-2*jt)]  # E pos
#   Qmat[sub_diagPos] <- Qmat[sub_diagPos] + phivec[jt]^2
# }
# for (jt in seq_len(jmax2)) {
#   subdiag_Pos <- (diagPos+jt)[-(nc+1L-seq(jt))] # F_j pos
#   Qmat[subdiag_Pos] <- Qmat[subdiag_Pos] - phivec[jt]
# }
# for (jt in seq_len(min(ord-1L,jmax1))) {
#   for (it in seq_len(ord-jt)) {
#     subdiag_Pos <- (diagPos+jt)[-(nc+1L-seq(jt))] # F_j pos
#     sub_subdiag_Pos <- subdiag_Pos[it+seq(length(subdiag_Pos)-2*it)] # G pos
#     Qmat[sub_subdiag_Pos] <- Qmat[sub_subdiag_Pos] + phivec[it]*phivec[it+jt]
#   }
# } 
# The chol crossfac is almost Toeplitz(c(1,-phivec)) except for the lower right block. Is there a way to write it directly? 
# There is a strand of litt ignoring this, but then the cov mat is not toeplitz


.condcorr2canonARpar <- function(condcorr) {
  phivec <- condcorr
  for (kt in seq_along(condcorr)[-1L]) {
    its <- 1L:(kt-1L)
    phivec[its] <- phivec[its] - condcorr[kt]*phivec[kt-its]
  }
  phivec
}

ARp <- function(p=1L, fixed=NULL, corr=TRUE, tpar=1/(1+seq(p))) { 
  force(corr)
  
  oldZlevels <- oldZrange <- postfit_phivec <- NULL
  initialize <- function(Zmatrix, ...) {
    oldZrange <<- range(as.integer(colnames(Zmatrix)))
    oldZlevels <<- paste(oldZrange[1L]:oldZrange[2L])
  }
  
  calc_Qmat_ARp <- function(parvec, newlevels) {    # Conversion to the phi parameters of the autoregression 
    phivec <-.condcorr2canonARpar(condcorr=parvec)
    # Build precision matrix from the phi coeffs
    nc <- length(newlevels)
    ord <- length(parvec)
    jmax1 <- (nc-1L)%/%2
    jmax2 <- min(ord,jmax1)
    diagvals <- rep(1,nc)
    subdiags <- vector("list", jmax2)
    for (jt in seq_len(jmax2)) {
      pos_in_diag <- jt+seq(nc-2*jt)
      diagvals[pos_in_diag] <- diagvals[pos_in_diag] + phivec[[jt]]^2
    }
    for (jt in seq_len(jmax2)) subdiags[[jt]] <- rep( - phivec[jt], nc-jt)
    for (jt in seq_len(min(ord-1L,jmax1))) {
      for (it in seq_len(ord-jt)) {
        pos_in_subdiag <- it+seq(length(subdiags[[jt]])-2*it)
        subdiags[[jt]][pos_in_subdiag] <- subdiags[[jt]][pos_in_subdiag] + phivec[it]*phivec[it+jt]
      }
    } 
    Qmat <- bandSparse(nc,k=c(0,seq(ord)), diagonals=c(list(diagvals), subdiags), symmetric = TRUE)
    if (corr) { # to return inverse correlation matrix
      sta_var <- 1/prod(1-parvec^2) 
      Qmat <- Qmat*sta_var 
    } # else return inverse covariance matrix for unit innovation variance.
    rownames(Qmat) <- newlevels 
    Qmat
  }
  
  
  Cf <- function(parvec) calc_Qmat_ARp(parvec=parvec, newlevels=oldZlevels )

  calc_moreargs <- function(corrfamily, ...) {
    np <- length(corrfamily$parnames)
    init <- rep(0.01,np)
    lower <- rep(-0.999,np)
    upper <- rep(0.999,np)
    names(init) <- names(lower) <- names(upper) <- corrfamily$parnames

    list(init=init, lower=lower, upper=upper)
  }
  
  ..calc_corr_from_dist <- function(ranFix, char_rd, distmat, ...) { # The AR1 code use distance matrices to handle the nested AR1 case...
    # function not used, but may be a useful template
    parvec <- ranFix$corrPars[[char_rd]]
    levelrange <- range(as.integer(.unlist(dimnames(distmat))))
    Qmat <- calc_Qmat_ARp(parvec=parvec, newlevels=seq(levelrange[1L],levelrange[2L]))
    corr <- .precision2cov(Qmat) 
    corr[rownames(distmat),colnames(distmat)] 
  }
  
  make_new_corr_lists <- function(newLv_env, which_mats, ranFix, newZAlist, new_rd, old_rd, ...) { 
    parvec <- ranFix$corrPars[[as.character(new_rd)]] 
    newlevels <- colnames(newZAlist[[new_rd]])
    newrange <- range(as.integer(newlevels))
    levelrange <- range(c(oldZrange,newrange))   
    # Qmat <- calc_Qmat_ARp(parvec=parvec, newlevels=newlevels)
    # corr <- .precision2cov(Qmat)  # correct but so ugly. 
    if (is.null(postfit_phivec)) postfit_phivec <<- .condcorr2canonARpar(condcorr=parvec)
    ACF <- stats::ARMAacf(lag.max=diff(levelrange),ar=postfit_phivec,ma=NULL) 
    corr <- stats::toeplitz(ACF)
    colnames(corr) <- rownames(corr) <- seq(levelrange[1L],levelrange[2L])

    newLv_env$cov_newLv_oldv_list[[new_rd]] <- corr[newlevels,oldZlevels, drop=FALSE]
    if (which_mats$nn[new_rd]) {
      newLv_env$cov_newLv_newLv_list[[new_rd]] <- corr[newlevels,newlevels, drop=FALSE]
    } else { 
      newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- rep(1,length(newlevels)) # allowing subsetting
    }
  }
  
  list(Cf=Cf, tpar=tpar, type="precision", initialize=initialize, fixed=fixed, calc_moreargs=calc_moreargs, 
       levels_type="time_series",
        #calc_corr_from_dist=calc_corr_from_dist,
       make_new_corr_lists=make_new_corr_lists,
       sparsePrec=TRUE, possiblyDenseCorr=TRUE,
       tag="ARp")
}


ARMA <- function(p=1L, q=1L, fixed=NULL, tpar=c(1/(1+seq_len(p)),1/(1+seq_len(q)))) { 
  force(p)
  force(q)
  
  fullparnames <- c(if(p>0L) {paste0("p",seq_len(p))}, if(q>0L) {paste0("q",seq_len(q))})
  names(tpar) <- fullparnames # only way to enforce specific names is to do assign them in the constructor, in tpar or in its body
  
  oldZrange <- oldZlevels <- postfit_phivec <- NULL
  initialize <- function(Zmatrix, ...) {
    oldZrange <<- range(as.integer(colnames(Zmatrix)))
    oldZlevels <<- seq(oldZrange[1L],oldZrange[2L])
  }
  
  calc_Cmat_from_dist <- function(parvec, levelrange, phivec=NULL) {
    arpos <- seq_len(p)
    if (is.null(phivec)) phivec <- .condcorr2canonARpar(condcorr=parvec[arpos]) 
    ACF <- stats::ARMAacf(lag.max=diff(levelrange),ar=phivec,ma=parvec[-arpos]) 
    Cmat <- stats::toeplitz(ACF)
    rownames(Cmat) <- colnames(Cmat) <- seq(levelrange[1L],levelrange[2L])
    Cmat
  }
  
  Cf <- function(parvec) calc_Cmat_from_dist(parvec=parvec, levelrange=oldZrange )
  
  calc_moreargs <- function(corrfamily, ...) {
    varpos <- which (corrfamily$parnames %in% fullparnames) 
    init <- rep(0.01,p+q)[varpos]
    lower <- c(rep(-0.999,p),rep(-Inf,q))[varpos]
    upper <- c(rep(0.999,p),rep(Inf,q))[varpos]
    names(init) <- names(lower) <- names(upper) <- corrfamily$parnames
    
    list(init=init, lower=lower, upper=upper)
  }

  ..calc_corr_from_dist <- function(ranFix, char_rd, distmat, ...) {
    # function not used, but may be a useful template
    parvec <- ranFix$corrPars[[char_rd]]
    levelrange <- range(as.integer(.unlist(dimnames(distmat))))
    corr <- calc_Cmat_from_dist(parvec=parvec, levelrange)
    corr[rownames(distmat),colnames(distmat)] 
  }
  
  make_new_corr_lists <- function(newLv_env, which_mats, ranFix, newZAlist, new_rd, old_rd, ...) { 
    parvec <- ranFix$corrPars[[as.character(new_rd)]] 
    newlevels <- colnames(newZAlist[[new_rd]])
    levelrange <- range(c(oldZrange,as.integer(newlevels)))  
    if (is.null(postfit_phivec)) postfit_phivec <<- .condcorr2canonARpar(condcorr=parvec[seq(p)])
    corr <- calc_Cmat_from_dist(parvec=parvec, levelrange, phivec=postfit_phivec)

    newLv_env$cov_newLv_oldv_list[[new_rd]] <- corr[newlevels,oldZlevels, drop=FALSE]
    if (which_mats$nn[new_rd]) {
      newLv_env$cov_newLv_newLv_list[[new_rd]] <- corr[newlevels,newlevels, drop=FALSE]
    } else { 
      newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- rep(1,length(newlevels)) # allowing subsetting
    }
  }
  
  list(Cf=Cf, tpar=tpar, initialize=initialize, fixed=fixed, calc_moreargs=calc_moreargs, 
       levels_type="time_series",
       # calc_corr_from_dist=calc_corr_from_dist,
       make_new_corr_lists=make_new_corr_lists,
       sparsePrec=FALSE, possiblyDenseCorr=TRUE,
       tag="ARMA")
}

.ranGCA_Z2A <- function(Zlevels) {
  Z2A <- strsplit(Zlevels,":") 
  Z2A <- do.call(rbind,Z2A)
  which_sort <- which(Z2A[,1]>Z2A[,2])
  Z2A[which_sort,c(1L,2L)] <- Z2A[which_sort,c(2L,1L)]
  Afactor <- factor(unique(unlist(Z2A)))
  Alevels <- levels(Afactor)
  Z2A <- data.frame(ID1=factor(Z2A[,1], levels=Alevels),
                    ID2=factor(Z2A[,2], levels=Alevels))
  list(Afactor=Afactor, Z2A=Z2A)
}


ranGCA <- function() {
  
  oldZlevels <- NULL

  initialize <- function(Zmatrix, ...) {
    oldZlevels <<- colnames(Zmatrix)
  }
  
  Af <- function(newdata, 
                 term,
                 fit.=FALSE,
                 ...) {
    if (fit.) {
      Zlevels <- oldZlevels # provided by initialize()
    } else {
      rhs <- term[[2L]][[3L]]
      txt <- .DEPARSE(rhs) ## should be the rhs of (|) cleanly converted to a string by terms(formula,data) in .get_terms_info()
      RHS_info <- .as_factor(txt=txt,mf=newdata, type=parent.env(environment())$levels_type) # 'levels_type' as provided there by .preprocess_corrFamily
      Zlevels <- levels(RHS_info$factor)
    }
    Z2Ablob <- .ranGCA_Z2A(Zlevels)
    Z2A <- Z2Ablob$Z2A
    Afactor <- Z2Ablob$Afactor
    Alevels <- levels(Z2Ablob$Afactor)
    if (length(Alevels)==1L) { # that must be a single 'homozygote'
      Amatrix <- .trivial_incidMat[,rep(1L,2L), drop=FALSE] 
    } else {
      Amatrix <- sparse.model.matrix(~ID1-1,Z2A)+sparse.model.matrix(~ID2-1,Z2A) # *sum* of the two individual random effects
    }
    colnames(Amatrix) <- Alevels
    rownames(Amatrix) <- Zlevels # ... matches the levels of Z with reciprocal ordered pairs i:j and j:i...
    # Now using the "global" or the local Afactor depending on 'fit.':
    attr(Amatrix,"is_incid") <- FALSE
    
    Amatrix
  }

  Cf <- function(tpar) NULL 
  
  make_new_corr_lists <- function(newLv_env, which_mats, newZAlist, new_rd, ...) { 
    # see code for (.|.) ranefs: (twhere the spaceial case for fix_info does not seem to require any adaptation here)
    newLv_env$cov_newLv_oldv_list[new_rd] <- list(NULL) # .calc_sub_diagmat_cov_newLv_oldv() will be called as for (.|.)
    if (which_mats$nn[new_rd]) {
      newLv_env$cov_newLv_newLv_list[new_rd] <- list(NULL) 
    } else { 
      newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- rep(1,ncol(newZAlist[[new_rd]])) # allowing subsetting
    }
  }
  
  list(Af=Af, tpar=numeric(0L), Cf=Cf, initialize=initialize, 
       make_new_corr_lists=make_new_corr_lists,
       levels_type="data_order",
       sparsePrec=TRUE, possiblyDenseCorr=FALSE,
       tag="ranGCA")
}

if (FALSE) {
  fit_rho <- fitme(dista ~ dad + mom + diff.mean.age.photo.st + msex + corrFamily(1|bb1+bb2), data=for1stLMM, 
                   covStruct=list(corrFamily=ranGCA()), 
                   control.HLfit=list(algebra="decorr"), # spares 20s otherwise needed for spaMM to select this method
                   verbose=c(TRACE=1), prior.weights=nbphoto, method="ML")
}

########  diallel  #########################################################################

.calc_splitord <- function(Zlevels) {
  splitord <- strsplit(Zlevels,":") 
  splitord <- do.call(rbind,splitord)
  which_sort <- which(splitord[,1]>splitord[,2]) # hmff. comparison on strings...
  splitord[which_sort,c(1L,2L)] <- splitord[which_sort,c(2L,1L)]
  splitord
}

.make_Cor_template_diallel <- function(usplitall) {
  # large but integer matrix, fast operations:
  CorNA <- (sapply(usplitall[,1L],`==`,y=usplitall[,1L])+
              (id12 <- sapply(usplitall[,1L],`==`,y=usplitall[,2L]))+ t(id12) +
              sapply(usplitall[,2L],`==`,y=usplitall[,2L]))
  CorNA <- Matrix::forceSymmetric(as(CorNA,"sparseMatrix"))
  CorNA
}


diallel <- function(tpar=0.25, fixed=NULL, public=NULL) {
  
  force(public)
  CorNA <- rhopos <- oldAfactor <- newAfactor <- oldZlevels <- newZlevels <- usplitold <- usplitnew <-  NULL

  Af <- function(newdata, 
                 term,
                 fit.=FALSE,
                 ...) {
    if (fit.) {
      Zlevels <- oldZlevels # provided by initialize()
      Afactor <- oldAfactor # provided by initialize()
    } else {
      rhs <- term[[2L]][[3L]]
      txt <- .DEPARSE(rhs) ## should be the rhs of (|) cleanly converted to a string by terms(formula,data) in .get_terms_info()
      RHS_info <- .as_factor(txt=txt,mf=newdata, type=parent.env(environment())$levels_type) # 'levels_type' as provided there by .preprocess_corrFamily
      newZlevels <<- Zlevels <- levels(RHS_info$factor)
      
      splitord <- .calc_splitord(Zlevels)
      Afactor <- paste0(splitord[,1],":",splitord[,2])
      #  At this point Afactor must allow duplicate rows for reciprocal pairs so that  rownames(Amatrix) <- Zlevels will work
      newAfactor <<- Afactor <- factor(Afactor,levels=unique(Afactor)) # i.e. keep the order determined by cols of Z.   
      
      usplitnew <<- unique(splitord) # saved for prediction
    }
    #
    #
    Alevels <- levels(Afactor)
    if (length(Alevels)==1L) {
      Amatrix <- .trivial_incidMat[rep(1L,length(Afactor)),, drop=FALSE] # maybe not most efficient, but avoids warning for possible inefficiency elsewhere
      # a general approach as in .calc_ZMatrix is to construct  a 'modmat' for the LHS, 
      #   an 'im' matrix for the RHS (the grouping levels; with special case when only one),
      #   and finally call .calc_raw_ZA(incidMat=im, modmat).
      # For Amatrices the 'modmat' is always trivial, but a special case for only one new A level is still needed.
    } else {
      Amatrix <- sparse.model.matrix(~Afactor-1,data.frame(Afactor=Afactor))
    }
    colnames(Amatrix) <- Alevels
    rownames(Amatrix) <- Zlevels # ... matches the levels of Z with reciprocal ordered pairs i:j and j:i...
    # Now using the "global" or the local Afactor depending on 'fit.':
    attr(Amatrix,"is_incid") <- TRUE
    Amatrix
  }
  
  initialize <- function(Zmatrix, ...) { # called before $Af()
    
    oldZlevels <<- oldZlevels <- colnames(Zmatrix)
    
    splitord <- .calc_splitord(oldZlevels)
    
    Afactor <- paste0(splitord[,1],":",splitord[,2])
    #  At this point Afactor must allow duplicate rows for reciprocal pairs so that  rownames(Amatrix) <- Zlevels will work
    oldAfactor <<- Afactor <- factor(Afactor,levels=unique(Afactor)) # i.e. keep the order determined by cols of Z.   
    
    usplitold <<- usplitold <- unique(splitord) # saved for prediction
    
    if (is.null(public$CorNA)) {
      # define template for correlation matrix
      CorNA <- .make_Cor_template_diallel(usplitold)
      x <- CorNA@x
      idpos <- which(x==2)
      
      rhopos <<- rhopos <- which(x==1)
      
      x[idpos] <- 1
      x[rhopos] <- NA_real_
      CorNA@x <- x
      colnames(CorNA) <- rownames(CorNA) <- levels(Afactor)
      CorNA <<- CorNA
      
      if (is.environment(public)) {
        public$CorNA <- CorNA
      }
    } else {
      rhopos <<- rhopos <- which(is.na(CorNA@x))
    }
  }
  
  Cf <- function(rho) {
    CorNA@x[rhopos] <- rho
    CorNA
  }
  
  make_new_corr_lists <- function(newLv_env, which_mats, ranFix, ranefs, newZAlist, new_rd, old_rd, ...) {
    
    usplitall <- unique(rbind(usplitnew,usplitold))
    allCorNA <- .make_Cor_template_diallel(usplitall)
    x <- allCorNA@x
    idpos <- which(x==2)
    rhopos <<- rhopos <- which(x==1)
    
    x[idpos] <- 1
    x[rhopos] <- ranFix$corrPars[[as.character(old_rd)]]
    allCorNA@x <- x
    
    allClevels <- paste0(usplitall[,1],":",usplitall[,2])
    rownames(allCorNA) <- colnames(allCorNA) <- allClevels
    
    newClevels <- colnames(newZAlist[[new_rd]]) # paste0(usplitnew[,1],":",usplitnew[,2])
    
    if (which_mats$no) newLv_env$cov_newLv_oldv_list[[new_rd]] <- structure(allCorNA[newClevels, 
                                                                                     rownames(CorNA), # reordering relative to allClevels
                                                                                     drop=FALSE])
    if (which_mats$nn[new_rd]) {
      newLv_env$cov_newLv_newLv_list[[new_rd]] <- allCorNA[newClevels, newClevels ,drop=FALSE]
    } else {
      newLv_env$diag_cov_newLv_newLv_list[[new_rd]] <- rep(1,length(newClevels))
    }
    # no return value, the newLv_env has been modified
  }
  
  

  calc_moreargs <- function(corrfamily, ...) {
    lower <- setNames(0, nm=corrfamily$parnames) # matrix easily not diag-dominant for negative values...
    upper <- setNames(0.49999, nm=corrfamily$parnames)
    list(lower=lower,upper=upper)
  }
  
  list(tpar=tpar, Cf=Cf, Af=Af, initialize=initialize, fixed=fixed, calc_moreargs=calc_moreargs, 
       make_new_corr_lists=make_new_corr_lists,
       levels_type="data_order",
       sparsePrec=FALSE, possiblyDenseCorr=TRUE,
       tag="diallel", public=public )
}

if (FALSE) {
  # match with the additive case.  Decimals of fixed matter.
  # the variance of the single ranef must be the double of the variances of the individual random effect.
  fit_rho <- fitme(dista ~ dad + mom + diff.mean.age.photo.st + msex + corrFamily(1|bb1+bb2), data=for1stLMM, 
                   covStruct=list(corrFamily=diallel(tpar=0.42,fixed=c(p1=0.499999))), 
                   control.HLfit=list(algebra="decorr"), # spares 20s otherwise needed for spaMM to select this method
                   init=list(phi=6.54159, lambda=3.25445),
                   verbose=c(TRACE=1), prior.weights=nbphoto, method="ML")
  
  
  fit_rho <- fitme(dista ~ dad + mom + diff.mean.age.photo.st + msex + corrFamily(1|bb1.bb2), data=for1stLMM, 
                   covStruct=list(corrFamily=diallel(tpar=0.42)), 
                   control.HLfit=list(algebra="decorr"), # spares 20s otherwise needed for spaMM to select this method
                   verbose=c(TRACE=1), prior.weights=nbphoto, method="ML", 
                   lower=list(corrPars=list("1"=c(p1=0.4))), 
                   upper=list(corrPars=list("1"=c(p1=0.45))), 
                   init=list(corrPars=list("1"=c(p1=0.42))))
  
}




