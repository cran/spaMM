
## declared "arglist" to print a clean summary instead of very long list
## the print(summary()) is required if the "arglist" is a member (eg as.list(<call>));
## summary alone would print nothing
print.arglist <- function(x,...) {
  if (inherits(x,"environment")) { ## not claer that this case shoudl occur, but it occurs notably in HLCor.args of refit
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
    col_indices <- rep(1L:(ncol(X)),diff(X@p))
    X@x <- X@x * Dvec[col_indices]    
    ## a triangular matrix (dtC) with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
    if ( methods::.hasSlot(X, "diag") && X@diag=="U") Matrix::diag(X) <- Dvec
  } else if (inherits(X,"dsCMatrix")) { 
    if (check_dsC) {
      if (generally_dgC <- (length(unique(Dvec))>1L)) {
        # The correct result is presumably NOT symmetric so direct manipulation of slots of X-dsCMatrix is not appropriate. 
        # We go through dgC (correct result in all cases) and may try converting back to dsC afterwards.
        X <- as(X,"dgCMatrix")
      }
      col_indices <- rep(1L:(ncol(X)),diff(X@p))
      X@x <- X@x * Dvec[col_indices]
      if (generally_dgC) {
        if (keep_dsC) {
          warning("Suspect call: return value presumably not symmetric, but dsCMatrix requested.")
          X <- as(X,"symmetricMatrix")
        } # else warning(".Matrix_times_Dvec(<dsCMatrix>,...) returns a dgCMatrix.") # we explictly used non-default keep_dsC, so warning not useful
      } # else we already have a correct dsC
    } else { # we assume that keeping dsC is OK: risky, non-default, but for precisionBlocks it's really hat we want (so no warning).
      # warning("call to .Matrix_times_Dvec(<dsCMatrix>, ..., check_dsC=FALSE) is possibly inefficient code.")
      col_indices <- rep(1L:(ncol(X)),diff(X@p))
      X@x <- X@x * Dvec[col_indices]
    }
  } else {
    warning("inefficient code in .make_Xscal or .Matrix_times_Dvec") ## eg dgeMatrix: dense matrix in the S4 Matrix representation
    # but also dtC for which the warning is not warranted
    X <- X %*% Diagonal(x=Dvec)
  } 
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
} ## end if (FALSE)

.m_Matrix_llblock_times_Dvec <- function(X, Dvec) {
  n_u_h <- length(Dvec)
  if (inherits(X,"ddiMatrix")) {
    warning("Suspect call to .Matrix_ZAL_ulblock_times_Dvec()") # the matrix is square hence X.pv not NULL, and ZAL is a zero block !
  } else if (inherits(X,c("dgCMatrix","dsCMatrix", "dtCMatrix"))) { 
    which_i_affected_rows <- X@i>(n_u_h-1L)
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
  return(X)
}

.Dvec_times_Matrix_lower_block <- function(Dvec,X,min_row) { # min_row from 0
  if (inherits(X,"ddiMatrix")) {
    if (X@diag=="U") { ## diag + unitary => identity matrix
      X <- Diagonal(x=c(rep(1,min_row),Dvec))
    } else X@x <- X@x * c(rep(1,min_row),Dvec) ## raw diag matrix
  } else if (inherits(X,c("dgCMatrix", "dtCMatrix"))) { 
    which_i_affected_rows <- X@i>(min_row-1L)
    X@x[which_i_affected_rows] <- X@x[which_i_affected_rows]*Dvec[X@i[which_i_affected_rows]-min_row+1L]
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
      return(.tcrossprod( ZW, as_sym=as_sym)) ## dsCMatrix if as_sym=TRUE (default)
    }
  } else if (inherits(ZAL,"ZAXlist")) {
    sqrtW <- sqrt(w)
    ZW <- .m_Matrix_times_Dvec(ZAL,sqrtW)
    return(.tcrossprod( ZW,as_sym=as_sym)) ## dsCMatrix if as_sym=TRUE (default)
  } else return(.ZWZt(ZAL,w))
} ## so seems always dsC

#.ZWZtbase <- function(ZAL,w) tcrossprod(sweep(ZAL,MARGIN=2L,w,`*`),ZAL) # slower that Rcpp version, with no obvious benefits
#  !  not the same as     ZAL %*% sweep(ZAL,MARGIN=1L,w,`*`) 

.ZtWZwrapper <- function(ZAL,w) { ## used in several contexts
  # if (ncol(ZAL)==0L) {
  #   stop(".ZtWZwrapper called with ncol(ZAL)=0") ## temporary devel code since all calls are in principle protected 
  # } else 
  if (inherits(ZAL,"ZAXlist")) {
    tcrossfac <- .crossprod(ZAL,Diagonal(x=sqrt(w))) # t() %*% sqrtw as dsCMatrix # by ad-hoc 
    return(.tcrossprod(tcrossfac))  #
  } else if (inherits(ZAL,"Matrix")) {
    # but Diagonal() is slow and it's better to produce matrices of symmetric type by construction so that forceSymmetric is not needed
    # Matrix::.symDiagonal better for addition with another symmetric matrix, see its doc.
    if (inherits(ZAL,"ddiMatrix")) { ## if diagonal
      if ( ZAL@diag=="U") {
        return(.symDiagonal(x=w)) ## diagonal and unitary => identity
      } else return(.symDiagonal(x = w * diag(ZAL)^2)) ## this case may never occur
    } else  {
      sqrtW <- sqrt(w) # this uses sqrt() hence assumes w>0
      DZAL <- ZAL
      DZAL@x <- DZAL@x * sqrtW[DZAL@i+1L] ## W Z
      ## a triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
      if (  methods::.hasSlot(DZAL, "diag") &&  DZAL@diag=="U") Matrix::diag(DZAL) <- sqrtW
      return(.crossprod(DZAL))  
    }
  } else return(.ZtWZ(ZAL,w)) ## but this uses sqrt() hence assumes w>0
} # Matrix result is always of symmetric type (dsC or ddi)

.process_ranCoefs <- function(processed, ranCoefs, trRanCoefs, invLs=NULL, use_tri_Nspprec, need_longLv=TRUE) {
  if (! is.null(invLs)) ranCoefs <- trRanCoefs
  ranCoefs_blob <- processed$ranCoefs_blob
  ZAlist <- processed$ZAlist 
  cum_n_u_h <- processed$cum_n_u_h
  if (is.null(ranCoefs_blob)) { ## in call from .preprocess(), allows simplified user input as numeric vector
    Xi_cols <- attr(ZAlist, "Xi_cols")
    isRandomSlope <- Xi_cols>1L ## FIXME seems oK for later code but semantically sloppy, cf (X-1|id) terms
    hasRandomSlope <- any(isRandomSlope)
    longLv_templates <- vector("list", length(Xi_cols))
    if (hasRandomSlope) {
      for(rd in seq_len(length(Xi_cols))) {
        if ((Xi_ncol <- Xi_cols[rd])>1L) {
          templa <- matrix(seq(Xi_ncol^2),nrow=Xi_ncol,ncol=Xi_ncol)
          longLv_templates[[rd]] <- .makelong(templa, ncol(ZAlist[[rd]])) ## templa is dense array of indices
        }
      }
    }
    ranCoefs_blob <- list(isRandomSlope=isRandomSlope, is_set=(isRandomSlope & FALSE),
                          longLv_templates=longLv_templates)
    if (hasRandomSlope && ! is.null(ranCoefs)) { ## ranCoefs in preprocess: fixed values
      nrand <- length(isRandomSlope)
      if (is.numeric(ranCoefs)) {
        if (sum(isRandomSlope)==1L) {
          ranCoefs <- vector("list", nrand)
          ranCoefs[[isRandomSlope]] <- ranCoefs
        } else stop("Cannot match numeric ranCoefs to a single random-coefficient term.")
      } else if (length(ranCoefs) < nrand ){
        posnames <- setdiff(seq(nrand), names(ranCoefs))
        ranCoefs <- c(ranCoefs,structure(vector("list",length(posnames)),names=posnames))
        ranCoefs <- ranCoefs[as.character(seq(nrand))]
      }
      newly_set <- ! (sapply(ranCoefs,is.null) | sapply(ranCoefs,identical,y=NA)) ## is.na would test each element of each ranCoefs element  
      # newly_set contains some TRUE (unless the user provided an uninformative ranCoefs) => final return() 
      # no return yet ! => final return() 
    } else return(ranCoefs_blob)
  } else if (is.null(ranCoefs) || ! any(isRandomSlope <- ranCoefs_blob$isRandomSlope)) {
    return(ranCoefs_blob)
  } else { ## both ranCoefs_blob and ranCoefs non-NULL, && any(isRandomSlope)
    ## assume 'incomplete' list
    nrand <- length(isRandomSlope)
    if (length(ranCoefs) < nrand ){
      posnames <- setdiff(seq(nrand), names(ranCoefs))
      ranCoefs <- c(ranCoefs,structure(vector("list",length(posnames)),names=posnames))
      ranCoefs <- ranCoefs[as.character(seq(nrand))]
    }
    newly_set <- ! (sapply(ranCoefs,is.null) | sapply(ranCoefs,identical,y=NA)) ## is.na would test each element of each ranCoefs element  
    newly_set <- newly_set & ( ! ranCoefs_blob$is_set)
    if ( any(newly_set)) {
      Xi_cols <- attr(ZAlist, "Xi_cols")
      ##  & continue final block of code untilfinal return(), skipping immediate return 
    } else return(ranCoefs_blob)
  }
  ## newly_set contains some TRUE
  spprecBool <- processed$sparsePrecisionBOOL
  nrand <- length(ZAlist)
  if (is.null(LMatrices <- ranCoefs_blob$LMatrices)) LMatrices <- vector("list", nrand)
  if (is.null(lambda_est <- ranCoefs_blob$lambda_est)) lambda_est <- numeric(cum_n_u_h[nrand+1L])
  min_lam <- .spaMM.data$options$tol_ranCoefs_inner["tol"]
  for (rt in which(isRandomSlope & newly_set)) {
    Xi_ncol <- Xi_cols[rt]
    if (need_longLv || is.null(invLs)) { ## only the econd case appears to be true sometimes, and only in develcode...
      compactcovmat <- .calc_cov_from_ranCoef(ranCoef=ranCoefs[[rt]], Xi_ncol=Xi_ncol)
      latentL_blob <- .calc_latentL(compactcovmat, use_tri_Nspprec=use_tri_Nspprec, spprecBool=spprecBool, try_chol= ! spprecBool)
      #  with(latentL_blob,.ZWZt(design_u,d)) = compactcovmat hence design_u is more of a tcrossprod factor
      ## we have a repres in terms of ZAL and of a diag matrix of variances; only the latter affects hlik computation
      longLv <- .makelong(latentL_blob$design_u,longsize=ncol(ZAlist[[rt]]),
                          template=ranCoefs_blob$longLv_templates[[rt]]) ## the variances are taken out in $d
      attr(longLv,"trRancoef") <- .ranCoefsFn(ranCoefs[[rt]], rC_transf=.spaMM.data$options$rC_transf_inner) # used for init_trRancoef in .MakeCovEst()
      attr(longLv,"Ltype") <- "cov"
    } else {
      char_rt <- as.character(rt)
      latentL_blob <- list(crossfac_Q=invLs[[char_rt]], Xi_ncol=Xi_ncol,
                           d=rep(1,Xi_ncol) ## used a few lines below
      )
      longLv <- .makelong(latentL_blob$crossfac_Q, longsize=ncol(ZAlist[[rt]]), 
                          template=ranCoefs_blob$longLv_templates[[rt]])  
      attr(longLv,"trRancoef") <- trRanCoefs[[char_rt]]
      attr(longLv,"Ltype") <- "prec"
    }
    attr(longLv,"latentL_blob") <- latentL_blob ## kept for updating in next iteration, strucList, predVar, and output
    attr(longLv,"ranefs") <-  structure(attr(ZAlist,"exp_ranef_strings")[rt], 
                                              type= attr(ZAlist,"exp_ranef_types")[rt] )
    attr(longLv, "corr.model") <- "random-coef"
    LMatrices[[rt]] <- longLv
    n_levels <- ncol(ZAlist[[rt]])/Xi_ncol 
    u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
    lambda_est[u.range] <- .make_long_lambda(latentL_blob$d, n_levels, Xi_ncol=Xi_ncol)
    lambda_est[lambda_est<min_lam] <- min_lam ## arbitrarily small eigenvalue is possible for corr=+/-1 even for 'large' parvec
  }
  ranCoefs_blob$is_set <- newly_set
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
                          ranefEstargs,
                          ranCoefs_blob,
                          lambda.Fix, ## distinctly used for CARdispGammaGLM
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
      glm_lambda <- .calc_CARdispGammaGLM(data=locdf, lambda.Fix=lambda.Fix[rd], lev=lev_lambda[u.range],control=control)
      attr(glm_lambda,"whichrand") <- rd
      next_lambda_est[u.range] <- fitted(glm_lambda) ## prediction of heteroscedastic variances
      coeffs <- coefficients(glm_lambda)
      if (is.na(lambda.Fix[rd])) { 
        rho <- - coeffs[["adjd"]]/ coeffs[["(Intercept)"]]
      } else {
        rho <- - coeffs[1]*lambda.Fix
      }
      HLfit_corrPars[[char_rd]]$rho <- rho 
      done[rd] <- TRUE
    }
  }
  ### next the (no L matrices) or (Matern model and other fixed L matrix cases)
  for (it in which( ! done )) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    if (is.na(unique.lambda <- lambda.Fix[it])) {
      prior_lam_fac <- rand.families[[it]]$prior_lam_fac
      if (is.null(prior_lam_fac)) wt <- 1 else wt <- 1/prior_lam_fac
      resp_lambda[u.range] <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=wt) ## must give d1 in table p 989 de LeeN01
      unique.lambda <- sum(resp_lambda[u.range])/sum(1-lev_lambda[u.range]) ## NOT in linkscale 
      unique.lambda <- max(unique.lambda,1e-8) # FR->FR still corrected
      unique.lambda <- min(unique.lambda,maxLambda)  
      if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { ## Gamma(identity)
        unique.lambda <- pmin(unique.lambda,1-1/(2^(iter+1)))  ## impose lambda<1 dans ce cas 
      }
    } 
    if ( ! is.null(prior_lam_fac <- rand.families[[it]]$prior_lam_fac)) {
      next_lambda_est[u.range] <- unique.lambda*prior_lam_fac
    } else next_lambda_est[u.range] <- rep(unique.lambda,length(u.range))
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
              isRandomSlope=isRandomSlope
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

.calcPHI <- function(oriFormula, ## with offset
                    dev.res,family,data,
                    lev_phi,
                    phimodel,verbose,method="glm",control.phi=list(),
                    control) {
  if (phimodel != "phiHGLM") {  
    if ( identical(deparse(oriFormula),"~1")) { ## one case where we can easily avoid an explicit call to a glm (but one will be used to compute SEs later) 
      next_phi_est <- sum(dev.res)/sum(1-lev_phi) ## NOT in linkscale
      beta_phi <- c("(Intercept)"=family$linkfun(next_phi_est)) ## linkscale value
      glm_phi <- NULL
    } else { ## which means that calcPHI cannot be used for final estim phi
      glm_phi <- .calc_dispGammaGLM(formula=oriFormula, dev.res=dev.res,
                                   data=data,lev=lev_phi, family=family,
                                   etastart=control.phi$etastart, ## private, availabble, but NULL so far
                                   control=control)
      beta_phi <- coefficients(glm_phi) ## numeric(0) if phi fixed by some offset
      next_phi_est <- fitted(glm_phi)
      if (family$link!="log" && any(next_phi_est<=0)) { stop("Gamma GLM for dispersion yields negative phi estimates.") }
    }
    if (verbose["trace"] && length(beta_phi)) {cat("str(phi_est): ",str(next_phi_est))}
  } else { ## random effect(s) in predictor for phi
    stop("This function should not be called when a residual model with random effects is fitted.")
  } 
  return(list(next_phi_est=next_phi_est,  #low phi values are handled in calc_APHLs...
              glm_phi=glm_phi,
              beta_phi=beta_phi ## used at least to initiate final GLM in "~1" case
  )) ## 
}


.calc_w_resid <- function(GLMweights,phi_est) { ## One should not correct this phi_est argument by prior.weights (checked)
  phi_est[phi_est<1e-12] <- 1e-11 ## 2014/09/04 local correction, cf comment in calc_APHLS...
  if (ilg <- is.list(GLMweights)) { ## zero_truncated family
    res <- GLMweights ## includes WU_WT, d3logLthdth3....
    GLMweights <- GLMweights$truncGLMweights
  }
  #if (any(is.nan(GLMweights/phi_est))) stop("any(is.nan(GLMweights/phi_est))")
  w.resid <- structure(as.vector(GLMweights/phi_est),unique= (attr(GLMweights,"unique") && length(phi_est)==1L))
  if (ilg) {
    res[["w_resid"]] <- w.resid
    return(res)
  } else {
    return(w.resid)
  }
}

.calc_H_global_scale <- function(w.resid) { # where w.resid depends on phi_est and on prior weights
  # small adjustment added 2019/11/11 resolve some numerical pb in extreme binary fits (mu \approx y, w.resid->0)
  if (is.list(w.resid)) w.resid <- w.resid$w_resid
  u_w_resid <- unique(w.resid)
  if (length(u_w_resid)==1L) {
    H_global_scale <- 1/u_w_resid # so that weight_X is rep(1,...) 
  } else {
    H_global_scale <- exp(- mean(log(w.resid)))
  }
  # unbiased in w.resid=1 (if lop1p used,  difference between lop1p(x) and log(1+x) could create a 'bias' => don't use log1p):
  H_scale_regul <- .spaMM.data$options$H_scale_regul
  if (H_global_scale>1/H_scale_regul) H_global_scale <- exp(log(H_scale_regul+1) - mean(log(H_scale_regul+w.resid))) # test '871.8141' affected if this is not conditional
  return(H_global_scale) 
}

.calc_weight_X <- function(w.resid, H_global_scale) {
  if (is.list(w.resid)) {
    return(sqrt(H_global_scale*w.resid$w_resid))
  } else return(sqrt(H_global_scale*w.resid))
}

## spaMM_Gamma() fixes Gamma()$dev.resids(1e10+2,1e10,1) is < 0
# dev.resids() must be >0 for computation deviance_residual in fitting Gamma GLMM, and also for $aic() computation.
spaMM_Gamma <- function (link = "inverse") {
  mc <- match.call()
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
         COMPoisson = function(theta,y,nu) { ## theta = log(lambda) 
           COMP_nu <- environment(family$aic)$nu
           logLs <- numeric(length(y))
           for (i in seq(length(y))) {
             comp_z <- .COMP_Z(lambda=exp(theta[i]),nu=COMP_nu)
             logLs[i] <- nu[i]*(theta[i]*y[i]-comp_z[[1]]-log(comp_z[[2]])) - COMP_nu * lfactorial(y[i])
           }
           logLs[theta== -Inf & y==0] <- 1
           logLs
         },
         negbin = function(theta,y,nu) { ## theta is the canonical param, -log(1+shape/mu)
           NB_shape <- environment(family$aic)$shape
           if (is.null(mu <- attr(theta,"mu"))) mu <- NB_shape/(exp(-theta)-1) # note that NB_shape/(exp(log(1 + NB_shape/x)) - 1) numerically fails to be x for large x
           -family$aic(y=y, mu=mu, wt=1)/2 ## handles truncation from given untruncated mu
         },
         stop("code missing for this family")
  )
}


`.theta.mu.canonical` <- function(mu,family) { 
  ## the (fixed) canonical link between theta and mu, not the family link between eta and mu 
  if (inherits(family,"family")) {
    famfam <- family$family
  } else famfam <- family
  switch(famfam,
         gaussian = mu ,
         poisson = structure(log(mu),mu=mu) ,
         binomial = make.link("logit")$linkfun(mu),  
         ## if this does no work, use 
         #                 { 
         #                    theta <- logit(mu)
         #                    theta[theta>27.6310] <- 27.6310 ## mu>1-1e-12
         #                    theta[theta < -27.6310] <- -27.6310 
         #                    theta
         #                 },
         Gamma = structure(-1/mu,mu=mu), ## "-": McC&N p. 290
         COMPoisson = {
           if (is.null(lambda <- attr(mu,"lambda"))) {
             lambda <-  family$linkfun(mu,log=FALSE) ## valid bc only the canonical link is implemented 
           }  
           structure(log(lambda),mu=mu)
         },
         negbin = structure(-log(1+(environment(family$aic)$shape)/mu),mu=mu), ## keep mu b/c useful for clik_fn
         stop("code missing for this family") 
  )
} ## returns values for given mu

.thetaMuDerivs <- function(mu,BinomialDen,family) { ## used for non-canonical links
  familyfam <- family$family
  if (familyfam=="binomial") muFREQS <- mu/BinomialDen
  if (familyfam == "negbin") NB_shape <- environment(family$aic)$shape
  ## these definitions depend only on the canonical link
  Dtheta.Dmu <- switch(familyfam,
                       gaussian = rep(1,length(mu)) ,
                       poisson = 1/mu ,
                       binomial = 1/(muFREQS*(1-muFREQS)),
                       Gamma = 1/mu^2,
                       negbin = 1/(mu*(1+mu/NB_shape)),
                       stop("code missing")
                       # COMPoisson has no implemented non-canonical link in which case this point should not be reached
  ) ## values for given mu
  if (familyfam=="binomial") Dtheta.Dmu <- Dtheta.Dmu/BinomialDen
  D2theta.Dmu2 <- switch(familyfam,
                         gaussian = rep(0,length(mu)) ,
                         poisson = -1/mu^2 ,
                         binomial = -(1-2*muFREQS)/(muFREQS*(1-muFREQS))^2,
                         Gamma = -2/mu^3,
                         negbin = -(1+2*mu/NB_shape)/(mu*(1+mu/NB_shape))^2,
                         stop("code missing")
                         # COMPoisson again has no implemented non-canonical link
  ) ## values for given mu
  if (familyfam=="binomial") D2theta.Dmu2 <- D2theta.Dmu2/(BinomialDen^2)
  return(list(Dtheta.Dmu=Dtheta.Dmu,D2theta.Dmu2=D2theta.Dmu2))
}

.sanitize_eta_log_link <- function(eta, max, min=-max,y, nu=1) { # cf reasons below to always provide y 
  ### Before implementation of y argument:
  ## 'max' used to sanitize eta should be higher than the max used to sanitize the response in .calc_dispGammaGLM()
  ## Otherwise we would never have sanitized eta = sanitized y in cases where this would be the correct solution.
  ## More generally the LevM algo may then be stuck.
  ### BUT: additional concept in v.3.0.49; extend bound according to y (for both sanitized eta & sanitized y)
  #  test code with huge .calc_dispGammaGLM() response:
  # {
  #   spaMM.options(example_maxtime=70)
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
  if ( ! is.null(y)) {
    tmpmax <- pmax(log(y), tmpmax) ## if log(y)>tmpmax corrected eta is log(y)+(< log(1+double.xmax)=tmp1)*nu i.e. corrected eta is < log(y)+tmp1*nu
    #   otherwise  corrected eta is < (max-tmp1)*nu+(< log(1+double.xmax)=tmp1)*nu    <   max*nu  
    #  and "-tmpmax" ensures continuity in eta=tmpmax) as log(1+log(1+...))=0 there
    high_eta <- (eta>tmpmax)
    eta[high_eta] <- tmpmax[high_eta] + log(1+log(1+eta[high_eta]-tmpmax[high_eta]))*nu
  }  else eta[eta>tmpmax] <- tmpmax + log(1+log(1+eta[eta>tmpmax]-tmpmax))*nu
  ## smooth correction proved useful in spaMM_glm.fit (at least), presumably to avoid a plateau of objective fn
  tmpmin <- min+tmp1 # nu appears to have little effect on mu when eta is small
  if ( ! is.null(y)) {
    tmpmin <- pmin(tmpmin,log(y))
    low_eta <- (eta< tmpmin)
    eta[low_eta] <- tmpmin[low_eta] - log(1+log(1-eta[low_eta]+tmpmin[low_eta]))
  } else eta[eta< tmpmin] <- tmpmin - log(1+log(1-eta[eta< tmpmin]+tmpmin))## symmetric correction for negative eta
  return(eta)
} 

# using y is important in cases where a comparison btwn eta (or mu) and y is important for convergence of an algo
# In that case we sanitize only to the extant that individual values of y allow.
# If this sanitization is not enough, then it is y than must be sanitized.
.sanitize_eta <- function(eta, y=NULL,family,max=.spaMM.data$options$sanitize_eta["otherlog"]) {
  if (family$link =="log") {
    if (family$family=="gaussian") {
      eta <- .sanitize_eta_log_link(eta, max=.spaMM.data$options$sanitize_eta["gauslog"],  y=y) 
    } else eta <- .sanitize_eta_log_link(eta, max=max,  y=y)
  } else if (family$family == "COMPoisson" && family$link =="loglambda") {
    # this should be consistent with poisson(log) when nu=1
    COMP_nu <- environment(family$aic)$nu 
    eta <- .sanitize_eta_log_link(eta, max=max, y=y, nu=COMP_nu) ## will use log(mu) ~ eta/nu for large eta and small nu
  } else if (family$link=="inverse" && family$family=="Gamma") {
    etamin <- sqrt(.Machine$double.eps)
    neg_eta <- (eta<etamin)
    eta[neg_eta] <- etamin ## both eta and mu must be >0
    attr(eta,"any_neg_eta") <- any(neg_eta) ## actually 'not strictly positive'
  } else if (family$family == "binomial") { ## 
    #for binomial cases stats::binomial(...)$linkinv corrects eta for all links so that mu is always .Machine$double.eps from 0 or 1
    # This may not be enough... or inconsistent with corrections used elsewhere, so we overcome the stats:: correction
    # by correcting eta differently here and by computing mu <- .binomial_raw_linkinv(eta,link=family$link) (in .muetafn and possibly elsewhere)
    eta <- .binomial_corr_eta(eta,link=family$link, tol=.spaMM.data$options$bin_mu_tol)
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

.muetafn_truncated <- local({
  Rmpfr_warned <- FALSE
  function(family, mu, GLMweights, Vmu) {
    if (family$family=="poisson") { 
      ## D[Log[1 - E^-E^theta], {theta, 2}] /. {theta -> Log[mu]} // Simplify
      expmu <- exp(mu)
      p0 <- 1/expmu
      dlogMthdth <- -mu/(1-expmu) ## useful to correct rhs
      d2logMthdth2 <- -(1 + expmu * (mu-1))* mu/((expmu-1)^2)
      exp2mu <- expmu^2
      d3logMthdth3 <- mu * (1 + 3 * mu *(expmu - exp2mu) + (mu^2) *(expmu + exp2mu) + exp2mu - 2 * expmu )/(expmu-1)^3
    } else if (family$family=="negbin") { ## computations in TruncatedNegBin.nb
      ## D[Log[1 - E^-E^theta], {theta, 2}] /. {theta -> Log[mu...]} // Simplify
      shape <- .get_family_par(family, "shape")
      p0 <- (shape/(mu+shape))^shape ## (1-p)^r
      dlogMthdth <- (mu * p0)/(1-p0) ## family-specific relation btwn correction term M(th) for truncation and p0
      # d2logMthdth2 <- -((mu * p0 *(shape *(-1 + p0) + mu *(-1 + shape + p0)))/(shape *(-1 + p0)^2))
      # d3logMthdth3 <- -((mu * p0 * (shape^2 *(-1 + p0)^2 + 3 * mu * shape *(-1 + p0) * (-1 + shape + p0) + 
      #                                 mu^2 *(3 *shape *(-1 + p0) + 2 *(-1 + p0)^2 + shape^2 *(1 + p0))))/(
      #                                   shape^2 *(-1 + p0)^3))
      theta <- - log1p(shape/mu)
      umeth <- -expm1(theta)
      eth <- 1-umeth
      # num is -((mu * p0 *(shape *(-1 + p0) + mu *(-1 + shape + p0)))/(shape) an denome is (-1 + p0)^2
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
    if (any(wrong <- truncGLMweights<=0)) {
      if (family$family=="negbin") { 
        has_Rmpfr <- suppressWarnings(do.call("require",list(package="Rmpfr", quietly = TRUE))) # given it's nto in DESCRIPTION
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
                GLMweights=list(truncGLMweights=truncGLMweights,WU_WT=WU_WT,dlogMthdth=dlogMthdth, 
                                d2logMthdth2=d2logMthdth2, d3logMthdth3=d3logMthdth3)))
  }
})


.muetafn <-   function(eta,BinomialDen,processed) { ## note outer var BinomialDen 
  names(eta) <- NULL
  family <- processed$family
  eta <- .sanitize_eta(eta, family=family)
  if (family$family == "binomial") { ## 
    #for binomial cases stats::binomial(...)$linkinv corrects eta for all links so that mu is always .Machine$double.eps from 0 or 1
    # This may not be enough... or inconsistent with corrections used elsewhere, so we overcome the stats:: correction
    # by computing eta <- .binomial_corr_eta(eta,link=family$link, tol=.spaMM.data$options$bin_mu_tol) in .sanitize_eta() and by:
    mu <- .binomial_raw_linkinv(eta,link=family$link)
  } else if (family$family == "poisson") { ## same troubles as for binomial
    mu <- .poisson_raw_linkinv(eta,link=family$link)
  } else mu <- family$linkinv(eta) ## linkinv(eta) is FREQS for binomial, COUNTS for poisson...
  dmudeta <- family$mu.eta(eta) ## aberrant at hoc code for cloglog 'elsewhere'...
  Vmu <- family$variance(mu) 
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
  if (processed$LMMbool) {
    GLMweights <- eval(processed$prior.weights) ## with attr(.,"unique")
    attr(GLMweights,"unique") <- attr(processed$prior.weights,"unique") ## might actually be true sometimes
  } else  if (family$family=="Gamma" && family$link=="log") {
    GLMweights <- eval(processed$prior.weights)  ## with attr(.,"unique")
    attr(GLMweights,"unique") <- attr(processed$prior.weights,"unique") ## might actually be true sometimes
  } else {
    # in Gamma() (inverse) case dmudeta=Vmu => dmudeta^2 /Vmu= Vmu
    GLMweights <- eval(processed$prior.weights) * dmudeta^2 /Vmu ## must be O(n) in binomial cases
    attr(GLMweights,"unique") <- FALSE ## might actually be true sometimes
  }
  if (identical(family$zero_truncated,TRUE)) {
    resu <- .muetafn_truncated(family, mu, GLMweights, Vmu)
    resu$dmudeta <- dmudeta
    resu$sane_eta <- eta
    return(resu)
  } else return(list(mu=mu,dmudeta=dmudeta,GLMweights=GLMweights, sane_eta=eta))
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
  return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))
}

.updateW_ranefS <- function(cum_n_u_h,rand.families,lambda,u_h,v_h) {
  nrand <- length(rand.families)
  w.ranef <- dlogWran_dv_h <- dvdu <- vector("list",nrand)
  for (it in seq_len(nrand)) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    blob <- .updateWranef(rand.family=rand.families[[it]],lambda[u.range],u_h[u.range],v_h[u.range])
    w.ranef[[it]] <- blob$w.ranef
    dlogWran_dv_h[[it]] <- blob$dlogWran_dv_h
    dvdu[[it]] <- blob$dvdu
  }
  w.ranef <- unlist(w.ranef)
  w.ranef[w.ranef>1e10] <- 1e10 ## patch useful to avoid singular d2hdv2 in PLoG model
  dlogWran_dv_h <- unlist(dlogWran_dv_h)
  dvdu <- unlist(dvdu)
  resu <- list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=dvdu)
  ## the test is invalid for ranCoefs:
  # if (nrand==1L && rand.families[[1L]]$family=="gaussian") resu$unique_w.ranef <- w.ranef[[1L]] # used in sparsePrecision code
  #if (length(unique_w.ranef <- unique(w.ranef))==1L) resu$unique_w.ranef <- unique_w.ranef # used in sparsePrecision code
  return(resu)
}

.calc_d2mudeta2 <- function(link,mu=NULL,eta=NULL,BinomialDen=NULL) { ## d2 MuCOUNTS d etaFREQS^2
  switch(link,
         identity = 0,
         log = mu, 
         inverse = 2 * mu^3 , ## canonical for Gamma()
         ## next three make sense for Binomial response data
         logit = {muFREQS <- mu/BinomialDen;
                  d2muFREQS <- muFREQS*(1-muFREQS)*(1-2*muFREQS);
                  d2muFREQS * BinomialDen
         },
         probit = -eta * dnorm(eta) * BinomialDen, # assuming eta corrected as in binomial(probit)$mu.eta
         cloglog = {
           eta <- pmin(eta,700) ## as in binomial(cloglog)$mu.eta
           ## in particular d2mudeta2/dmudeta is then numerically OK
           exp(eta-exp(eta))*(1-exp(eta)) * BinomialDen ## D[1 - E^-E^\[Eta], {\[Eta], 2}]
          }, 
         cauchit = -2 *eta/(pi * (1+eta^2)^2),
         stop(paste("unhandled link'",link,"'in .calc_d2mudeta2()"))
  )
} 

#derivatives of GLM weights wrt eta 
.CMP_calc_dlW_deta_locfn <- function(i,lambdas,mu,COMP_nu) {
  lambdai <- lambdas[[i]]
  if (lambdai==0) {
    return(c(dmudeta=0, d2mudeta2=0))
  }
  ## ELSE
  moments <- .COMP_Pr_moments(lambda=lambdai,nu=COMP_nu,moments=c("2","3"))
  mui <- mu[[i]]
  dmu.dlogLambda <- moments[["2"]] - mui^2 # =family$mu.eta() without repeating some computations
  d2mu.dlogLambda2 <- moments[["3"]]-mui*moments[["2"]]-2*mui*dmu.dlogLambda
  return(c(dmudeta=dmu.dlogLambda, d2mudeta2=d2mu.dlogLambda2))
}

# returns bounded eta
.binomial_corr_eta <- function(eta, link=family$link, tol=.Machine$double.eps) {
  switch(link,
         logit= { # .Call(C_logit_linkinv, eta) de facto corrects eta has follows before computing mu:
           thresh <- -qlogis(tol)
           pmin(pmax(eta, -thresh), thresh)
         },
         probit= {
           ### .muetafn has called binomial(probit)$linkinv to compute mu(FREQS), which first corrected eta as follows:
           thresh <- -qnorm(tol)
           pmin(pmax(eta, -thresh), thresh) # so that abs(dnorm(eta)) is bounded by .Machine$double.eps   
         },
         cloglog= {
           lo <- log(-log1p(-tol))
           up <- log(-log(tol))
           pmin(pmax(eta, lo), up) 
          },
         cauchit= {
           thresh <- -qcauchy(tol)
           pmin(pmax(eta, -thresh), thresh)
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

.calc_dlW_deta <- function(dmudeta,family,mu,eta,BinomialDen,canonicalLink,calcCoef1=FALSE,w.resid) {
  coef1 <- NULL
  ## We first handle the canonical link cases, where comput. of coef1 depends only on the link  
  ## here w=dmudeta; d1=dwdmu dmudeta /w^2 = dlogwdeta/w = (d2mu/deta2)/(dmu/deta) /w =
  ##      (d2mu/deta2)/(dmu/deta)^2 = (d(dmudeta)/dmu)/dmudeta where d(dmudeta)/dmu is the numerator as detailed:
  if (canonicalLink) {
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
      dlW_deta <- (1-2*mu/BinomialDen)  
      if (calcCoef1) coef1 <- dlW_deta/dmudeta # F I X M E there's a pattern
      ## numerator is D[D[1/(1 + E^-\[Eta]), \[Eta]] /. {E^-\[Eta]->(1-\[Mu])/\[Mu]} ,\[Mu]]=1-2 mu 
    } else if (family$family=="Gamma") { ## link= "inverse" !
      ## numerator is D[D[-1/\[Eta], \[Eta]] /. {\[Eta] -> -1/\[Mu]}, \[Mu]] =2 mu 
      dlW_deta <- 2*mu
      if (calcCoef1) coef1 <- dlW_deta /dmudeta
    } else if (family$family=="COMPoisson") { 
      COMP_nu <- environment(family$aic)$nu
      lambdas <- exp(eta) ## pmin(exp(eta),.Machine$double.xmax) ##  FR->FR lambdas missing as mu attribute here ?
      blob <- sapply(seq(length(lambdas)), .CMP_calc_dlW_deta_locfn,lambdas=lambdas,mu=mu,COMP_nu=COMP_nu)
      dlW_deta <- blob["d2mudeta2",] / blob["dmudeta",]
      if (family$family=="COMPoisson") dlW_deta[is.nan(dlW_deta)] <- 0 ## quick patch for cases that should have low lik
      if (calcCoef1) {
        coef1 <- dlW_deta / blob["dmudeta",]
        if (family$family=="COMPoisson") coef1[is.nan(coef1)] <- 0 ## idem
      }
    } 
  } else if (family$family=="binomial" && family$link=="probit") { ## ad hoc non canonical case 
    pmax_dnorm_eta <- dnorm(eta)
    muFREQS <- mu/BinomialDen
    VmuFREQS <- pnorm(eta,lower.tail = FALSE)*pnorm(eta)  # more accurate than muFREQS*(1-muFREQS)
    dlW_deta <- -2*eta - pmax_dnorm_eta*(1-2*muFREQS)/VmuFREQS
    if (calcCoef1) {
      coef1 <- dlW_deta *(VmuFREQS)/ (BinomialDen * pmax_dnorm_eta^2) 
      # or coef1 <- (-1+2*(mu-VmuFREQS*eta/pmax_dnorm_eta))/(BinomialDen*pmax_dnorm_eta)
    }
  } else if (family$family=="Gamma" && family$link=="log") { ## ad hoc non canonical case 
    dlW_deta <- rep(0L,length(mu)) ## because they both involve dW.resid/dmu= 0
    if (calcCoef1) coef1 <- dlW_deta
  } else { ## "negbin" only called with non-canonical links always ends here 
    ## we need to update more functions of mu...
    tmblob <- .thetaMuDerivs(mu,BinomialDen,family)
    Dtheta.Dmu <- tmblob$Dtheta.Dmu # calcul co fn de muFREQS puis / BinomialDen
    D2theta.Dmu2 <- tmblob$D2theta.Dmu2 # calcul co fn de muFREQS puis / BinomialDen ^2
    d2mudeta2 <- .calc_d2mudeta2(link=family$link,mu=mu,eta=eta,BinomialDen=BinomialDen)
    ## ... to compute this:
    D2theta.Deta2_Dtheta.Deta <- (D2theta.Dmu2 * dmudeta^2 + Dtheta.Dmu * d2mudeta2)/(Dtheta.Dmu * dmudeta)
    dlW_deta <- d2mudeta2 / dmudeta + D2theta.Deta2_Dtheta.Deta
    if (identical(family$zero_truncated,TRUE)) {  
      # at this point we have the untruncated dlW_deta. We evaluate the truncated dlW_deta dlog_tildeW_deta
      dlW_deta <- dlW_deta + w.resid$WU_WT *
        (D2theta.Dmu2/Dtheta.Dmu * w.resid$d2logMthdth2 +  Dtheta.Dmu * w.resid$d3logMthdth3) * Dtheta.Dmu* dmudeta
    } 
    ## in truncated case we have coef1 = (truncated dlW_deta)/ W_untrunc so the following is still correct
    if (calcCoef1) coef1 <- dlW_deta / (Dtheta.Dmu * dmudeta^2) ## note that coef2 is indep of the BinomialDen, but coef1 depends on it 
  }
  res <- list(dlW_deta=dlW_deta,## dlW_deta equiv coef2
              coef1=coef1)
  if (identical(family$zero_truncated,TRUE)) {
    res$WU_WT <- w.resid$WU_WT
  }
  return(res) 
}

.safesolve_qr_vector <- function(qr.a,b,silent=TRUE,stop.on.error=TRUE) { ## solve.qr with fall-back; qr.a should be a qr object, b must be a vector
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
        if (stop.on.error) {
          stop("backsolve() failed in .safesolve_qr_vector().") ## perhaps recover A by qr.X and solve(A) ?
        } else return(res) ## passes control to calling function
      } 
    }  
  }
  return(res)
}

.is_identity <- function( x, matrixcheck=FALSE, tol=1e-8 ) {
  #if (is.null(is_id <- attr(x,"is_identity"))) {
    if (inherits(x,"Matrix")) {
      if (inherits(x,"ddiMatrix") ) return(x@diag=="U")
      if (! isDiagonal( x ) ) return( FALSE ) # F I X M E isDiagonal is slow (tests on spatial predictions ?)
      ## hence now diagonal:
      return(max(abs(range(diag(x))-1))<tol)
    } else {
      if (matrixcheck) {
        return(ncol(x)==nrow(x) && max(abs(x- diag(ncol(x))))<tol)
      } else return(FALSE)  
    }
  #} else return(is_id)
}

.CHOL <- function(x, CHOLtol=1e-8, ...) {
  res <- try(chol(x, ...), silent=TRUE)
  if (inherits(res,"try-error")) {
    cholpert <- CHOLtol*diag(ncol(x))
    res <- 2*chol(x+cholpert, ...) - chol(x+2*cholpert, ...)
  }
  return(res)
}

.Cholwrap <- .wrap_Ltri_t_chol <- function(mat, use_eigen=.spaMM.data$options$USEEIGEN) { 
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
    ## this was the typical case for crossr22: a dsyMatrix (dense symmetric)
    # so in fact it's better to avoid costly Matrix operations in .calc_r22 that yields crossr22 as a dsyMatrix
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

.tcrossprod <-  function(x, y=NULL, allow_as_mat=TRUE, as_sym=TRUE) {
  if (is.null(x)) return(NULL) ## allows lapply(,.tcrossprod) on a listof (matrix or NULL)
  if (inherits(x,"ZAXlist")) {
    return(tcrossprod(x,y)) # ZAXlist method
  } else if (inherits(x,"dCHMsimpl")) {
    if (is.null(y)) {
      resu <- solve(x, system="A") # cov matrix from chol_Q
    } else resu <- t( solve(x, b=as(x,"pMatrix") %*% t(y), system="L") ) 
  } else if (inherits(x,"Matrix") || inherits(y,"Matrix")) {
    if (is.null(y)) {
      if (inherits(x,"dgCMatrix")) { 
        xxt <- .dgCtcrossprod(x,y) 
        if (as_sym) {
          return(forceSymmetric(xxt)) ## not as(.,"sparseMatrix") which checks -> slow
        } else return(xxt)
      } else if (allow_as_mat && (maxd <- prod(dim(x)))>1L ) { 
        x_denseness <- .calc_denseness(x)/maxd 
        if (x_denseness>0.35) {
          resu <- .tcrossprodCpp(as.matrix(x),y) ##  so much faster (simulate(mrf,...) )-- poor Matrix:: coding ?
          if (as_sym) resu <- as(resu,"symmetricMatrix") # but the as(, ) may produce a quite dense dsC (or dsy)
          #                                                             
          colnames(resu) <- rownames(resu) <- rownames(x)
          return(resu)
        } else return(Matrix::tcrossprod(x))
      } else return(Matrix::tcrossprod(x))
    } else if (FALSE && inherits(x,"dgCMatrix") &&  inherits(y,"dgCMatrix") ) { ## FASTer when (FALSE &&)
      xyt <- .dgCtcrossprod(x,y) 
      return(xyt)
    } else if (allow_as_mat) {
      if (inherits(x,"Matrix") && (maxd <- prod(dim(x)))>1L) { # trivial >1L bc the alternatives to Matrix::tcrossprod are always useful
        x_denseness <- .calc_denseness(x)/maxd 
        if (x_denseness>0.35) x <- as.matrix(x)
      } #else x_denseness <- 0 ## ie if the test is FALSE, no explicit conversion is attempted
      if (inherits(y,"Matrix") && (maxd <- prod(dim(y)))>1L) { # idem
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
    if (is.null(y)) {
      colnames(resu) <- rownames(resu) <- rownames(x)
    } else {
      rownames(resu) <- rownames(x)
      colnames(resu) <- rownames(y)
    }
    return(resu)
  }
}

# currently not used:
.Matrix_crossprod <- function(x,y, as_sym=TRUE) { # x Matrix; y Matrix or NULL
  if (inherits(x,"dgCMatrix")) {
    if (is.null(y)) {
      txy <- .dgCcrossprod(x,NULL) 
      #if (class(Matrix::crossprod(x)) != class(forceSymmetric(.dgCcrossprod(x,NULL)))) browser()
      if (as_sym) {
        return(forceSymmetric(txy)) ## not as(.,"sparseMatrix") which checks -> slow
      } else return(txy)    
    } else if (inherits(y,"dgCMatrix")) {
      #txy <- Matrix::crossprod(x,y); if (!inherits(txy,"dgCMatrix")) browser(); 
      return(.dgCcrossprod(x,y)) ## return(Matrix::crossprod(x,y)) ## 
    } else return(Matrix::crossprod(x,y))
  } else return(Matrix::crossprod(x,y))
}

.crossprod <- function(x, y=NULL, allow_as_mat=TRUE, as_sym=is.null(y), use_Rcpp_cp=.spaMM.data$options$Rcpp_crossprod, 
                       eval_dens=TRUE) {
  if (is.null(x)) return(NULL) ## allows lapply(,.tcrossprod) on a list of (matrix or NULL)
  if (inherits(x,"ZAXlist")) {
    return(crossprod(x,y)) # ZAXlist method
  } else if (use_Rcpp_cp) {
    if ( inherits(x,c("dgCMatrix","dgeMatrix","matrix"))
         && inherits(y,c("dgCMatrix","dgeMatrix", "matrix","NULL"))) {
      if (use_Rcpp_cp>2L) {
        print(c("+",class(x),class(y))) ## ddi...
        utils::flush.console()
      }
      txy <- .Rcpp_crossprod(x, y, eval_dens) # returns dge; next forceSymmetric() can produce dgC (sp) OR dsy (dense)
      if (as_sym) { txy <- forceSymmetric(txy) } # absolutely necessary for use as update(., parent=txy); see ?`CHMfactor-class`
      if ( ( ! allow_as_mat) && ( ! inherits(txy,"sparseMatrix"))) txy <- as(txy,"sparseMatrix") # fortunately rarely needed but occurs in HLCor(lh~1 +AR1(1|time) example
      return(txy)
    } else { ## show unhandled classes and continue
      if (use_Rcpp_cp>1L) {
        print(c("-",class(x),class(y))) ## ddi...
        utils::flush.console()
      }
    }
  }
  if (inherits(x,"Matrix") || inherits(y,"Matrix")) { 
    # it's OK to use Matrix::crossprod() rather than .crossprod() when x has a known quite-sparse structure such as chol_Q
    if (is.null(y)) {
      if (TRUE && inherits(x,"dgCMatrix")) { ## FAST when (TRUE &&)
        txx <- .dgCcrossprod(x,y) 
        if (as_sym) {
          return(forceSymmetric(txx)) ## not as(.,"sparseMatrix") which checks -> slow
        } else return(txx)
      } else if (allow_as_mat && (maxd <- prod(dim(x)))>1L ) { 
        x_denseness <- .calc_denseness(x)/maxd 
        if (x_denseness>0.35) {
          ## Matrix::crossprod() is inefficient on math-dense matrices, even for dgeMatrix ! We try to maintain its return type, typically dsCMatrix:
          resu <- .crossprodCpp(as.matrix(x),y)
          if (as_sym) {
            return(forceSymmetric(resu)) 
          } else return(resu)
        } else return(Matrix::crossprod(x)) ## return(.Matrix_crossprod(x,y=NULL,as_sym=as_sym)) ## 
      } else return(Matrix::crossprod(x)) ## return(.Matrix_crossprod(x,y=NULL,as_sym=as_sym)) ## 
    } else if (FALSE && inherits(x,"dgCMatrix") &&  inherits(y,"dgCMatrix") ) { ## FASTer when (FALSE &&)
      txy <- .dgCcrossprod(x,y) 
      return(txy)
    } else if (allow_as_mat) {
      if (inherits(x,"Matrix") && (maxd <- prod(dim(x)))>1L) { # trivial >1L bc the alternatives to Matrix::tcrossprod are always useful
        x_denseness <- .calc_denseness(x)/maxd 
        if (x_denseness>0.35) x <- as.matrix(x)
      } #else x_denseness <- 0 ## ie if the test is FALSE, no explicit conversion is attempted
      if (inherits(y,"Matrix") && (maxd <- prod(dim(y)))>1L) { # idem
        y_denseness <- .calc_denseness(y)/maxd
        if (y_denseness>0.35) y <- as.matrix(y)
      }
      if (is.matrix(x) && is.matrix(y)) { # sparsity overhead is ~ 1/0.35
        return(.crossprodCpp(x, y)) # .crossprodCpp() slightly better than crossprod()
      } else return(Matrix::crossprod(x,y)) ## return(.Matrix_crossprod(x,y)) ## either is <0.35 so Matrix was retained. but dgeMatrix-y ends here
    } else return(Matrix::crossprod(x,y)) ## return(.Matrix_crossprod(x,y)) ## ! allow_as_mat
  } else {
    if (is.integer(x)) storage.mode(x) <- "double"
    if (is.integer(y)) storage.mode(y) <- "double"
    resu <- .crossprodCpp(x,y) ## faster than base::crossprod, but the matrix cannot be integer
    if (is.null(y)) {
      colnames(resu) <- rownames(resu) <- colnames(x)
    } else {
      rownames(resu) <- colnames(x)
      colnames(resu) <- colnames(y)
    }
    return(resu)
  }
}


.get_tcrossfac_beta_w_cov <- function(fit) { 
  if (is.null(fit$envir$tcrossfac_beta_w_cov)) { ## tcrossfac...w stored under $envir while tcrossfac...v is under $beta_v_cov
    if ("AUGI0_ZX_sparsePrecision" %in% fit$MME_method) {
      tcrossfac_beta_v_cov <- .get_tcrossfac_beta_v_cov(fit$X.pv, fit$envir, fit$w.resid) ## typically permuted from triangular
      pforpv <- ncol(fit$X.pv)
      lhs <- Matrix::bdiag(new("dtCMatrix",i= 0:(pforpv-1L), p=0:(pforpv), Dim=c(pforpv,pforpv),x=rep(1,pforpv)),
                           fit$envir$chol_Q)
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
          tcrossfac_beta_w_cov[v.range,] <- as.matrix(.crossprod(invL, tcrossfac_beta_w_cov[v.range,])) # invL may be sparse
          fit$envir$tcrossfac_beta_w_cov <- tcrossfac_beta_w_cov
        }
      } else fit$envir$tcrossfac_beta_w_cov <- tcrossfac_beta_w_cov
    }
  }
  return(fit$envir$tcrossfac_beta_w_cov) ## repeated calls to predVar suffices to make it useful to save it.
}

.eval_as_mat_arg <- function(object) { 
  (
    object$HL[1L]=="SEM" || ## SEM code does not yet handle sparse as it uses a dense Sig matrix
    ! identical(object$QRmethod,"sparse") ## => conversion to matrix if object$QRmethod is NULL
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
      invColdoldList <- lapply(vec_n_u_h,Diagonal)
      if (res$spaMM.version < "2.2.116") {
        ranefs <- attr(res$ZAlist,"ranefs") 
      } else ranefs <- attr(res$ZAlist,"exp_ranef_strings") 
      for (Lit in seq_len(length(strucList))) {
        lmatrix <- strucList[[Lit]]
        if (!is.null(lmatrix)) {
          ## end of mat_sqrt implies either type is cholL_LLt or we have decomp $design_u and $d 
          type <-  attr(lmatrix,"type")
          invCoo <- NULL
          if ( ! is.null(latentL_blob <- attr(lmatrix,"latentL_blob"))) { ## from .process_ranCoefs
            if ( ! is.null(gmp_design_u <- latentL_blob$gmp_design_u)) {
              #cat(crayon::red("YES"))
              invC <- gmp::solve.bigq(gmp::tcrossprod(gmp_design_u))
            } else {
              invC <- solve(.tcrossprod(latentL_blob$design_u ))
            }
            invCoo <- .makelong(invC,longsize=ncol(lmatrix))
          } else if (inherits(lmatrix,"dCHMsimpl")) { # before any test on type...
            invCoo <- tcrossprod(as(lmatrix,"sparseMatrix")) # assuming 'lmatrix' is an unpermuted CHM factor of the precision factor (as for attr(lmatrix,"Q_CHMfactor"))
          } else if (type == "from_AR1_specific_code")  {
            invCoo <- crossprod(solve(lmatrix)) # cost of a sparse triangular solve.
          } else if (type == "from_Q_CHMfactor")  {
            invCoo <- tcrossprod(as(attr(lmatrix,"Q_CHMfactor"),"sparseMatrix")) ## correct but requires the attribute => numerical issues in computing Q_CHMfactor
          } else if (type == "cholL_LLt")  {
            Rmatrix <- t(lmatrix)
            if (attr(lmatrix,"need_gmp")) {
              if (is.null(fix_predVar <- control$fix_predVar)) { # best avoided in programming, to avoid the .check_frames() stuff
                # Hence calling predict() with explicit control$fix_predVar is recommended in programming.
                # Yet without explicit control$fix_predVar, this block will not be called for all predict() calls, only once for any problematic fit 
                # so the impact .check_frames() on performance should be negligible. 
                predVar_exceptions <- .spaMM.data$options$fix_predVar
                if (.check_frames(which=predVar_exceptions$"NA")) {
                  fix_predVar <- NA
                } else if (.check_frames(which=predVar_exceptions$"TRUE")) {
                  fix_predVar <- TRUE
                } else if (.check_frames(which=predVar_exceptions$"FALSE")) {
                  fix_predVar <- FALSE
                } 
              }
              if (is.null(fix_predVar)) { # when user can handle this
                warning("A nearly-singular correlation matrix was fitted; see help('fix_predVar') for handling this", immediate. = TRUE)
                fix_predVar_NULL <- TRUE
              } else if (is.na(fix_predVar)) { # When user can do nothing and gmp cannot be assumed: warning+fallback computation below 
                warning("A nearly-singular correlation matrix was fitted", immediate. = TRUE)
              } else if (fix_predVar) { 
                invCoo <- gmp::tcrossprod(gmp::solve.bigq(gmp::as.bigq(Rmatrix))) 
                #invCoo <- Rmpfr::mpfr(invCoo,128)
              } # (else other cases including fix_predVar=FALSE:) fallback computation of  invCoo below.
            }
          } else { ## Rcpp's symSVD, or R's eigen() => LDL (also possible bad use of R's svd, not normally used)
            condnum <- kappa(lmatrix,norm="1")
            if (condnum<1/regul.threshold) {
              decomp <- attr(lmatrix,attr(lmatrix,"type")) ## of corr matrix !
              invCoo <-  try(.ZWZt(decomp$u,1/decomp$d),silent=TRUE) ## the idea is to use ginv() if thisfails, rather than 
              #                                                         to attempt regularization, given that condnum has been tested as not large.
              if (inherits(invCoo,"try-error")) invCoo <- NULL
            }
            if (is.null(invCoo)) Rmatrix <- .lmwithQR(t(lmatrix),yy=NULL,returntQ=FALSE,returnR=TRUE)$R_scaled # no pivoting compared to qr.R(qr(t(lmatrix))) 
          }
          if (is.null(invCoo)){ 
            invCoo <- try(chol2inv(Rmatrix),silent=TRUE)
            if (inherits(invCoo,"try-error") || max(abs(range(invCoo)))> 1e12) {
              invCoo <- ginv(crossprod(Rmatrix))
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
  } else if (attr(w.resid,"unique")) {
    crossprodZAL <- .crossprod(ZAL)
    d2hdv2 <- - w.resid[1L] * crossprodZAL
    nc <- ncol(d2hdv2)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    d2hdv2[diagPos] <- d2hdv2[diagPos] - w.ranef 
  } else if (inherits(ZAL,"Matrix")) {
    d2hdv2 <- - .ZtWZwrapper(ZAL,w.resid)     
    nc <- ncol(d2hdv2)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    d2hdv2[diagPos] <- d2hdv2[diagPos] - w.ranef 
  } else d2hdv2 <- .Rcpp_d2hdv2(ZAL,w.resid,w.ranef) 
  #cat("\n new d2hdv2")
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

.calc_d2hdv2_info <- function(object, ZAL) {
  if (is.null(factor_inv_Md2hdv2 <- object$envir$factor_inv_Md2hdv2)) {
    if ( ! is.null(object$envir$G_CHMfactor)) {
      if (length(object$y) > ncol(object$envir$chol_Q)) { # We precompute inverse(d2hdv2) so that .calc_dvdlogphiMat_new() has ONE costly (r*r) %*% (r*n).
        if (is.null(factor_inv_Md2hdv2 <- object$envir$factor_inv_Md2hdv2)) {
          LL <- with(object$envir, Matrix::solve(object$envir$G_CHMfactor, as(object$envir$G_CHMfactor,"pMatrix") %*% chol_Q,system="L")) ## dgCMatrix 
          object$envir$factor_inv_Md2hdv2 <- Matrix::drop0(LL,tol=.Machine$double.eps)
        }
        d2hdv2_info <- - .crossprod(object$envir$factor_inv_Md2hdv2) 
      } else d2hdv2_info <- object$envir # then .calc_dvdlogphiMat_new() has TWO  (r*r) %*% (r*n) (plus an efficient solve())
    } else {
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
    }
  } else d2hdv2_info <- - .crossprodCpp(as.matrix(factor_inv_Md2hdv2), yy=NULL) # inverse(d2hdv2) # much faster than Matrix::crossprod(dsCMatrix...)
  return(d2hdv2_info)
}

.get_glm_phi <- function(fitobject) {
  if (is.null(fitobject$envir$glm_phi)) { 
    glm_phi_args <- c(fitobject$phi.object$glm_phi_args, 
                      list(formula=fitobject$resid.predictor,
                           lev=fitobject$lev_phi, data=fitobject$data,  
                           family= fitobject$resid.family)
    )
    fitobject$envir$glm_phi <-  do.call(".calc_dispGammaGLM", glm_phi_args)
  } 
  return(fitobject$envir$glm_phi)
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

.calc_Md2hdvb2_info_spprec_by_QR <- function(envir, which="tcrossfac_v_beta_cov") {## sparse qr is fast
  qr_R <- qrR(envir$qrXa, backPermute=FALSE) 
  X_scale <- envir$X_scale
  if ( ! is.null(X_scale)) {
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

.calc_Md2hdvb2_info_spprec_by_r22 <- function(X.pv, envir, w.resid, which="tcrossfac_v_beta_cov") {
  # currently called for "tcrossfac_v_beta_cov" for which r12,r22 needed
  # same if it is called for "R_Md2hdbv2" (as through .get_tcrossfac_beta_w_cov())  
  # if we used it for v_beta_cov, could we simplify the code by only computing crossprod_r22 using $invG_ZtWX ? What about r_12 ?
  
  #if (is.null(tcrossfac_v_beta_cov <- envir$tcrossfac_v_beta_cov)) # No strong reason to save it. The only re-use is for cAIC_pd
  {
    if (is.null(factor_inv_Md2hdv2 <- envir$factor_inv_Md2hdv2)) {
      LL <- with(envir, Matrix::solve(G_CHMfactor, as(G_CHMfactor,"pMatrix") %*% chol_Q,system="L")) ## dgCMatrix 
      envir$factor_inv_Md2hdv2 <- factor_inv_Md2hdv2 <- Matrix::drop0(LL,tol=.Machine$double.eps)
      # such that chol_Md2hdv2 = t(solve(factor_inv_Md2hdv2))
    }
    if (is.null(r22 <- envir$r22)){
      if (ncol(X.pv)==0L) { 
        r12 <- matrix(nrow=nrow(factor_inv_Md2hdv2),ncol=0L)
        r22 <- diag(nrow=0L)
      } else {
        r12 <- as(Matrix::solve(envir$G_CHMfactor, as(envir$G_CHMfactor,"pMatrix") %*% envir$ZtW %*% X.pv,system="L"),"sparseMatrix") ## solve(as(envir$G_CHMfactor,"sparseMatrix"), envir$ZtW %*% AUGI0_ZX$X.pv)
        r22 <- .calc_r22(X.pv,w.resid,r12, sXaug=NULL) ## both lines as explained in working doc
      }
    } else r12 <- envir$r12 
    r22 <- as(r22,"sparseMatrix") 
    n_u_h <- ncol(envir$chol_Q)
    pforpv <- ncol(X.pv)
    if (which=="R_Md2hdbv2") { ## not called in the source, and the single occurrence of solve(factor_inv_Md2hdv2)
      R_Md2hdbv2 <- rbind(cbind(t(solve(factor_inv_Md2hdv2)),r12), # cbind(chol_Md2hdv2,r12),
                          cbind(Matrix(0,nrow=pforpv,ncol=n_u_h),r22)) ## R_a in the working doc
      return(R_Md2hdbv2) ## such that tcrossprod(R_invMd2hdvb2) = solve(Md2hdbv2) = solve(crossprod(R_Md2hdbv2))
    } else {
      # weirdly solve(chol_Md2hdv2,system="L") (with/out b) is slower that solving the sparseMatrix.
      inv11 <- t(factor_inv_Md2hdv2)# solve(chol_Md2hdv2)
      inv22 <- solve(r22)
      inv12 <- -1* inv11 %*% (r12 %*% inv22) ## does not seem to understand unary '-'
      R_invMd2hdvb2 <- rbind(cbind(inv11,inv12),
                             cbind(Matrix(0,nrow=pforpv,ncol=n_u_h),inv22))
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


.calc_Md2hdvb2_info_spprec <- function(X.pv,envir, w.resid, which="tcrossfac_v_beta_cov") { ## this is called post-fit
  if ( ! is.null(envir$qrXa) ) { # If we have a qr facto of augmented X, wecan use it efficiently
    # This is tested for example by get_predVar(fitsparse,newdata=newXandZ) on 'landsat' data
    .calc_Md2hdvb2_info_spprec_by_QR(envir, which)
  } else {
    .calc_Md2hdvb2_info_spprec_by_r22(X.pv, envir, w.resid, which)
  }
}

.get_tcrossfac_beta_v_cov <- function(X.pv, envir, w.resid) {
  if (is.null(envir$beta_cov_info$tcrossfac_beta_v_cov)) {
    ## Using sparse matrices as much as possible for computation:
    tcrossfac_v_beta_cov <- .calc_Md2hdvb2_info_spprec(X.pv, envir, w.resid, which="tcrossfac_v_beta_cov")
    # but the v_beta_cov matrix is (mathematically) dense whatever its format
    pforpv <- ncol(X.pv)
    n_u_h <- ncol(envir$chol_Q)
    seqp <- seq_len(pforpv)
    perm <- c(n_u_h+seqp,seq_len(n_u_h))
    envir$beta_cov_info$tcrossfac_beta_v_cov <- tcrossfac_v_beta_cov[perm,,drop=FALSE] # useful to keep it for predVar computations?
  }
  return(envir$beta_cov_info$tcrossfac_beta_v_cov) ## Used for beta_v_cov itself, then for (beta_w_cov for predVar). Seem useful to  keep it
}

.calc_beta_cov_info_spprec <- function(X.pv, envir, w.resid) { 
  tcrossfac_beta_v_cov <- .get_tcrossfac_beta_v_cov(X.pv, envir, w.resid)
  if (is.null(beta_cov <- envir$beta_cov_info$beta_cov)) { # if stripHLfit()ed
    beta_v_cov <- .tcrossprod(tcrossfac_beta_v_cov) # if the L differs among method the meaning of v differs too, and beta_v_cov too  
    # ./. but we can check that is is = solve(crossprod(wAugX)) by reconstructing wAugX with the matching L as in .get_LSmatrix()
    pforpv <- ncol(X.pv)
    seqp <- seq_len(pforpv)
    beta_cov <- beta_v_cov[seqp,seqp,drop=FALSE]
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
      if (inherits(xmatrix,"dCHMsimpl")) {
        if (force_bindable) {
          xmatrix <- as(xmatrix,"pMatrix") %*% solve(xmatrix, system="Lt")
        } else bindable <- FALSE
      } else if (inherits(xmatrix,"bigq") ) {
        xmatrix <- .mMatrix_bigq(xmatrix) # this does not seem to be used for the Evar computation so we may drop precision
      } 
      if ( ! is.null(xmatrix)) {
        ZA <- ZAlist[[Lit]]
        if (.is_identity(ZA)) {
          if (inherits(xmatrix,"dCHMsimpl")) {
            ZAX[[Lit]] <- structure(list(ZA=ZA, Q_CHMfactor=xmatrix), class=c("ZA_QCHM","list")) 
          } else ZAX[[Lit]] <- xmatrix ## matrix may replace Matrix here...          
        } else {
          locnc <- ncol(ZA)
          locnr <- nrow(xmatrix)
          if ( locnc %% locnr !=0) {
            mess <- paste0("The number of levels of the grouping variable in random term ", 
                           attr(ZAlist,"exp_ranef_strings")[Lit])
            mess <- paste0(mess,"\n  is not a multiple of the dimension of the correlation matrix.") ## by distMatrix checking in corrHLfit or no.info check somewhere...
            stop(mess)
          }         
          nblocks <- locnc %/% locnr 
          if (nblocks>1) {
            # no longer clear when it occurs. Nested AR1 now uses full-dim chol_Q (see "from_AR1_specific_code")
            # and currently AR1 does not use ZAXlist, even when TRY_ZAX=TRUE (F I X M E?)
            # Any example of nested corrMatrix ?
            if (inherits(xmatrix,"dCHMsimpl")) {
              warning("Possibly inefficient code here (nested random effect+TRY_ZAX).") 
              xmatrix <- as(xmatrix,"pMatrix") %*% solve(xmatrix, system="Lt") 
            }
            locZA <- ZA
            for (bt in 1:nblocks) 
              locZA[,locnr*(bt-1)+(1:locnr)] <- locZA[,locnr*(bt-1)+(1:locnr)] %*% xmatrix
            ZAX[[Lit]] <- locZA
          } else {
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
  }
  if ( ! bindable) class(ZAX) <- c("list","notBindable") ## keep this order otherwise the ZAXlist constructor does not .../...
  #                                                   recognize a list since "notBindable" is not declared as extending "list"
  return(ZAX)
}

# even though the Z's were sparse postmultplication by LMatrix leads some of the ZAL's to dgeMatrix (dense) 
#                                 [and replacement by LMatrix may give a *m*atrix !]
.ad_hoc_cbind <- function(ZALlist, as_matrix ) {
  nrand <- length(ZALlist)
  if ( as_matrix ) {
    for (rd in seq_len(nrand)) ZALlist[[rd]] <- as.matrix(ZALlist[[rd]]) 
    if (nrand>1L) {ZAL <- do.call(cbind,ZALlist)} else ZAL <- ZALlist[[1L]]
  } else {
    for (it in seq_len(nrand)) {
      if ( # is.matrix(ZALlist[[it]]) || ## seems to work but at a cost for speed.
        inherits(ZALlist[[it]],"dgeMatrix")) ZALlist[[it]] <- as(ZALlist[[it]],"dgCMatrix")
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

.compute_ZAL <- function(XMatrix, ZAlist, as_matrix, force_bind=FALSE) {
  ZALlist <- .compute_ZAXlist(XMatrix,ZAlist, force_bindable=force_bind)
  if (inherits(ZALlist,"notBindable")) {
    return(new("ZAXlist", LIST=ZALlist))
  } else {
    ZAL <- .ad_hoc_cbind(ZALlist, as_matrix )
    return(ZAL)
  }
}



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
  return(unlist(notEQL))
}

.initialize_v_h <- function(psi_M,etaFix,init.HLfit,cum_n_u_h,rand.families,port_env) {
  v_h <- etaFix$v_h
  if (is.null(v_h) ) v_h <- port_env$port_fit_values$v_h
  if (is.null(v_h) ) v_h <- init.HLfit$v_h
  if (is.null(v_h) ) {
    v_h <- vector("list", length(rand.families))
    for (it in seq_len(length(rand.families))) {
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      v_h[[it]] <- rand.families[[it]]$linkfun(psi_M[u.range]) ## v as link(mean(u)) 
    }
    v_h <- unlist(v_h)
  }
  return(v_h)
}

## u_h and v_h in box constraints fro m unconstrainted v_h
.u_h_v_h_from_v_h <- function(v_h,rand.families,cum_n_u_h,lcrandfamfam,lower.v_h,upper.v_h) {
  if(!is.null(lower.v_h)) {v_h[v_h<lower.v_h] <- lower.v_h}
  if(!is.null(upper.v_h)) {v_h[v_h>upper.v_h] <- upper.v_h}
  nrand <- length(rand.families)
  u_list <- vector("list", nrand)
  for (it in seq_len(nrand)) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    u_list[[it]] <- rand.families[[it]]$linkinv(v_h[u.range])
    if (any(is.infinite(u_list[[it]]))) {
      warning("infinite random values ('u_h') were constrained to finite range.") 
      u_list[[it]] <- pmin(.Machine$double.xmax, pmax(-.Machine$double.xmax,u_list[[it]]) )
    }
  }
  u_h <- unlist(u_list)
  ## if there were box constr, v_h may have been modified, we put it in return value
  if ( ! (is.null(lower.v_h) && is.null(upper.v_h))) attr(u_h,"v_h") <- v_h
  return(u_h)
}

.make_lambda_object <- function(nrand, lambda_models, cum_n_u_h, lambda_est,                        
                               process_resglm_blob, rand.families, ZAlist, next_LMatrices, lambdaType) {
  ## redefine names(namesTerms) and namesTerms elements
  print_namesTerms <- attr(ZAlist,"namesTerms") ## a list, which names correspond to the grouping variable, and elements are the names of the coefficients fitted
  namesnames <- names(print_namesTerms)
  for (it in seq_len(length(namesnames))) if (nchar(namesnames[it])>10) namesnames[it] <- paste0(substr(namesnames[it],0,9),".")
  names(print_namesTerms) <- make.unique(namesnames,sep=".") ## makes group identifiers unique (names of coeffs are unchanged); using base::make.unique
  print_lambda <- vector("list",nrand) 
  Xi_cols <- attr(ZAlist,"Xi_cols")
  for (rd in seq_len(nrand)) {
    plam_rd <- process_resglm_blob$print_lambdas[[rd]]
    if (anyNA(plam_rd)) { ## those for which no glm was available, such as fixed lambdas...
      Xi_ncol <- Xi_cols[rd]
      plam_rd <- numeric(Xi_ncol)
      u.range <- (cum_n_u_h[rd]+1L):(cum_n_u_h[rd+1L])
      blocksize <- length(u.range)/Xi_ncol
      for (colit in seq_len(Xi_cols[rd])) {
        u.subrange <- u.range[((colit-1L)*blocksize+1L):(colit*blocksize)]
        plam_rd[colit] <- unique(lambda_est[u.subrange])
      } # With "chol" transfo pour ranCoefs, all the lambda's are 1, so subrange'ing is essential
    }
    if (lambdaType[rd]!="inner" && ## bc for inner, plam_rd is $print_lambdas[[rd]] and is already the single lambda param
        ! is.null(prior_lam_fac <- rand.families[[rd]]$prior_lam_fac)) {
      plam_rd <- plam_rd[1L]/prior_lam_fac[1L]
    } ## F I X M E It would be better to use a lambda_par attribute for lambda_est
    print_lambda[[rd]] <- structure(plam_rd, names=print_namesTerms[[rd]])
  }
  attr(print_lambda,"cum_n_u_h") <- cum_n_u_h
  lambda.object <- list(lambda_est = lambda_est,  ## full vector for simulate() calc_logdisp_cov()
                        lambda=print_lambda)  ## nrand-elements list, in output -> used by simulate, useful for init another fit, may substitute to the coefficients_lambdaS when the latter have not bee computed from  glm, etc.
  lambda.object$type <- lambdaType
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
  return(round(as.numeric(difftime(Sys.time(), time1, units = "secs")), 1))
}

.post_process_LMatrices <- function(next_LMatrices, ZAlist, ranCoefs_blob) { 
  strucList <- vector("list", length(ZAlist)) 
  for (it in seq_len(length(ZAlist))) {
    lmatrix <- next_LMatrices[[it]]
    ## keep in mind that str(S4...) does not show extra attributes
    ## lmatrix is an S4 object...
    ## ranefs itself has attribute type="(.|.)"
    if ( ! is.null(lmatrix)) {
      corr.model <- attr(lmatrix, "corr.model") 
      if (is.null(corr.model)) warning('attr(next_LMatrix, "corr.model") is NULL')
      if (corr.model=="random-coef") {
        design_u <- attr(lmatrix,"latentL_blob")$design_u
        if ( ( ! is.null(design_u)) && kappa(design_u)>1e06) attr(lmatrix,"latentL_blob")$gmp_design_u <- gmp::as.bigq(design_u)
        # singularity issues not handled here in spprec (crossfac_Q representation)
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
      # 'given' that I cannot modify strucList in post_fit code, the need_gmp attribute must be assigned now.  
      if (attr(lmatrix,"type")=="cholL_LLt" && corr.model %in% c("Matern","Cauchy")) {
        attr(lmatrix,"need_gmp") <- (kappa(lmatrix)>1e05)
      } else attr(lmatrix,"need_gmp") <- FALSE ## it's not that we may not need it, it's that we may not handle it
      if (FALSE && inherits(lmatrix,"dCHMsimpl")) { # inhibited => dCHMsimpl goes into strucList
        extras <- setdiff(names(attributes(lmatrix)),c("class",slotNames(lmatrix)))
        lmatrix_ <- as(lmatrix,"pMatrix") %*% solve(lmatrix, system="Lt") 
        for (st in extras) attr(lmatrix_,"st") <- attr(lmatrix,"st")
        attr(lmatrix_,"Q_CHMfactor") <- lmatrix 
        attr(lmatrix_, "type") <- "from_Q_CHMfactor"
        strucList[[it]] <- lmatrix_ 
      } else strucList[[it]] <- lmatrix  
    } 
  }
  return(structure(strucList, isRandomSlope=ranCoefs_blob$isRandomSlope)) ## isRandomSlope is boolean vector
}


.nothing_to_fit <- function(phi.Fix, off, models, etaFix, rand.families, cum_n_u_h, lcrandfamfam, 
                            lambda.Fix, vec_n_u_h, n_u_h, fixed_adjacency_info, ZAL, BinomialDen, processed) {
  ## nothing to fit. We just want a likelihood
  ### a bit the same as max.iter<1 ... ?
  phi_est <- phi.Fix
  eta <- off
  if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
    ## we need u_h in calc_APHLS...() and v_h here for eta...
    v_h <- etaFix$v_h
    u_h <- etaFix$u_h
    if (is.null(u_h)) {u_h <- .u_h_v_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,
                                                lcrandfamfam=lcrandfamfam,lower.v_h=NULL,upper.v_h=NULL)}
    lambda_est <- .resize_lambda(lambda.Fix,vec_n_u_h,n_u_h, adjacency_info=fixed_adjacency_info)
    eta <- eta + drop(ZAL  %id*id%  etaFix$v_h) ## updated at each iteration
  } ## FREQS
  ## conversion to mean of response variable (COUNTS for binomial)
  muetablob <- .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
  w.resid <- .calc_w_resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
  if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
    wranefblob <- .updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## no fit, likelihood computation
    H_global_scale <- .calc_H_global_scale(w.resid)
    ZAL_scaling <- 1/sqrt(wranefblob$w.ranef*H_global_scale) ## Q^{-1/2}/s
    Xscal <- .make_Xscal(ZAL, ZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX)
    weight_X <- sqrt(H_global_scale*w.resid) ## sqrt(s^2 W.resid)
    if (inherits(Xscal,"Matrix")) {
      mMatrix_method <- .spaMM.data$options$Matrix_method
    } else mMatrix_method <- .spaMM.data$options$matrix_method
    sXaug <- do.call(mMatrix_method,
                     list(Xaug=Xscal, weight_X=weight_X, w.ranef=wranefblob$w.ranef, H_global_scale=H_global_scale))
  } else sXaug <- NULL 
  res <- list(APHLs=.calc_APHLs_from_ZX(processed=processed, which="p_v", sXaug=sXaug, phi_est=phi_est, 
                                        lambda_est=lambda_est, dvdu=wranefblob$dvdu, u_h=u_h, mu=muetablob$mu))
  return(res)
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
  newi[Ipos] <- Ilen-1L+Iseq
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

.XDtemplate <- function(X) {
  if (inherits(X,"dtCMatrix")) {
    #XDtemplate <- .adhoc_rbind_dtC_dvec(X,dvec=rep(1,ncol(X))) # equivalent but with less ad-hoc unsafe code:
    return(.adhoc_rbind_dgC_dvec(as(X,"dgCMatrix"),dvec=rep(1,ncol(X))))
  } else if (inherits(X,"dgCMatrix")) {
    return(.adhoc_rbind_dgC_dvec(X,dvec=rep(1,ncol(X))))
  } else if (inherits(X,"Matrix")) {
    return(.adhoc_rbind_dgC_dvec(X,dvec=rep(1,ncol(X))))
  } else return(rbind(X, diag(nrow=ncol(X))))
}

.damping_to_solve <- function(X, XDtemplate=NULL, dampDpD, rhs=NULL,method="QR", .drop=TRUE) { ##  
  if (.drop) rhs <- drop(rhs) ## 1-col m/Matrix to vector ## affects indexing below but the result of the chol2inv line is always Matrix
  if (method=="QR") { ## seem always true
    if (is.null(XDtemplate)) {
      warning("Possibly inefficient call to .damping_to_solve() without precomputed XDtemplate")
      XDtemplate <- .XDtemplate(X)
    }
    nr <- nrow(XDtemplate)-length(dampDpD)
    if (inherits(XDtemplate,"Matrix")) { # both cases occur in routine use
      XD <- .Dvec_times_Matrix_lower_block(Dvec=sqrt(dampDpD),X=XDtemplate,min_row=nr)
      RP <- qr(XD) ## i.e. Matrix::qr
      RRsP <- sort.list(RP@q) 
      # solve(crossprod(XD)) = chol2inv(Matrix::qrR(RP,backPermute = FALSE))[RRsp,RRsP]
      # solve(crossprod(XD), rhs) = (Matrix::chol2inv(Matrix::qrR(RP, backPermute = FALSE)) %*% rhs[sort.list(RRsP)])[RRsP]
      if (is.null(rhs)) {
        return(list(inv=Matrix::chol2inv(qrR(RP,backPermute = FALSE)),Rperm=RP@q+1L,RRsP=RRsP))
      } else {
        ## ::chol2inv (-> essentially tcrossprod(solve()) without any particular optimization): very slow tcrossprod for large matrices,
        ## not specifically bc of the tcrossprod implementation (compared to crossprod) but because the result is *much* denser (bigranefs example) 
        ## It's then much better to avoid chol2inv(), at least when there is a RHS! 
        if (ncol(XD)<1000L) { # F I X M E ad hoc (fast tests vs 'bigranefs' test)
          if (is.matrix(rhs)) {
            resu <- (Matrix::chol2inv(qrR(RP,backPermute = FALSE)) %*% rhs[RP@q+1L,])[RRsP,]
          } else { ## if rhs is one col 
            resu <- (Matrix::chol2inv(qrR(RP,backPermute = FALSE)) %*% rhs[RP@q+1L])[RRsP]
          }
        } else {
          if (TRUE) {
            qrr <- qrR(RP,backPermute = FALSE)
            if (is.matrix(rhs)) {
              resu <- backsolve(qrr, backsolve(qrr,rhs[RP@q+1L,],transpose=TRUE))[RRsP,]
            } else { ## if rhs is one col 
              resu <- backsolve(qrr, backsolve(qrr,rhs[RP@q+1L],transpose=TRUE))[RRsP]
            }
          } else {
            solveR <- solve(qrR(RP,backPermute = FALSE))
            if (is.matrix(rhs)) {
              resu <- (solveR %*% .crossprod(solveR,rhs[RP@q+1L,]))[RRsP,]
            } else { ## if rhs is one col 
              resu <- (solveR %*% .crossprod(solveR,rhs[RP@q+1L]))[RRsP]
            }
          }
        }
      }
      ## assuming rhs is 'slim', this minimizes permutations. For computation of dv_h, it is square rather than slim...
    } else {
      XD <- .Dvec_times_matrix(Dvec=c(rep(1,nr),sqrt(dampDpD)), X=XDtemplate)
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
      #if ( ! is.null( processed$residProcessed$port_env$port_fit_values$corrPars)) browser()
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
  return(structure(lambda_est,init.lambda=init.lambda)) ## keep init.lambda provisionally bc SEMwrap expects it (FIXME)
  # F I X M E: I should add a lambda_par attribute
}

.calc_initial_init_lambda <- function(lambda.Fix, nrand, processed, ranCoefs_blob, init.HLfit, fixed) {
  init.lambda <- lambda.Fix ## already the right size 'nrand' with NA's or non-fixed ones
  lambdaType <- rep("",nrand) 
  lambdaType[processed$ranCoefs_blob$is_set] <- "fix_ranCoefs" 
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
    pos <- unlist(hyper_info$ranges[hy_names])  
    lambdaType[pos] <- hy_lam_types[hyper_info$map[pos]]
  }
  lambdaType[lambdaType=="" & ranCoefs_blob$is_set] <- "outer_ranCoefs" 
  ## "fixed" (not "fix" -- confusing) is tested eg to compute dfs p_lambda, in summary, in calc_logdisp_cov.R, in .get_logdispObject.
  lambdaType[lambdaType=="" & ! is.na(processed$lambda.Fix)] <- "fixed" 
  # F I X M E it seems that lambdaType evaluation up to this point can be preprocessed EXCEPT "outer_ranCoefs" (ranCoefs_blob$is_set !=processed$ranCoefs_blob$is_set)
  lambdaType[(lambdaType=="" & ! (is.na(lambda.Fix)))] <- "outer" # outer charac. by difference between processed$lambda.Fix and lambda.Fix
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


.update_port_fit_values <- function(old_obj, new_obj, port_fit_values, models, processed, control.HLfit,
                                    lambda_est, #we only need lambda_est to initiate the next inner optim of the mean fit
                                    phifit # we need more info to initiate both inner and outer optim of the next phifit
                                    ) {
  if ( new_obj> old_obj) {
    ## Use FALSE to inhibit all port_env usage:
    if ( ! identical(control.HLfit$write_port_env,FALSE)) assign("objective",new_obj,envir=processed$port_env) 
    if (abs(old_obj-new_obj)<5) { # from tests on fit_SPDE_spaMM <- fitme(.... Loaloa ...)
      # Update port_fit_values and assign it: 
      if (models[["eta"]]=="etaHGLM") {
        port_fit_values$lambda_est <- lambda_est ## mean fit lambda; to be used by .HLfit_finalise_init_lambda
        # 
      }
      processed$port_env$port_fit_values <- port_fit_values # NOT only for phiHGLM
      #assign("best_fit_values",port_fit_values,envir=processed$port_env) ## keep best values in all cases
      #
      # Updates for a phifit:
      if (models[["phi"]]=="phiHGLM") {
        # phifit (fitme -> HLfit) has already written in processed$residProcessed$port_env$port_fit_values (and this may be erased by the phifitarglist code): 
        # but the condition was on its own APHLs and it may actually have erased the info 
        ## To be used only in iter=0; for iter>0 see parallel code in HLfit_body()
        phifit_init_HLfit <- list(v_h=phifit$v_h)
        if (is.null(processed$residModel$fixed$fixef)) phifit_init_HLfit$fixef <- fixef(phifit) ## : Use an init if the parameters are not fixed. 
        processed$residProcessed$port_env$port_fit_values[["init_HLfit"]] <- phifit_init_HLfit
        processed$residProcessed$port_env$port_fit_values$corrPars <- .new_phifit_init_corrPars(phifit,innershift=0) 
        ## $corrPars : canonical and with full info about types including "outer"
        processed$residProcessed$port_env$port_fit_values$lambda_object <- phifit$lambda.object ## disp fit lambda info.
      }          
    }
  } else { ## keep objective value; two cases for init values:
    if (abs(old_obj-new_obj)<1) {
      ## small decrease: Keep old port_env_values
    } else {
      processed$port_env$port_fit_values <- NULL ## remove starting value, no useful for large variations 
      processed$residProcessed$port_env$port_fit_values <- NULL
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

# .backsolve() and .chol2solve are not used in operational code. They are best used once at a time, with variable thr_backsolve values.
.backsolve <- function(r, x=NULL, upper.tri = TRUE, transpose = FALSE) {
  anyMatrix <- (inherits(r,"Matrix") || inherits(x,"Matrix"))
  if (anyMatrix) warning("*M*atrix argument in call to .backsolve()")
  if ( anyMatrix || ncol(r)< .spaMM.data$options$thr_backsolve) { 
    if (is.null(x)) x <- diag(nrow=ncol(r))
    return(backsolve(r=r, x=x, upper.tri = upper.tri, transpose=transpose))
  } else {
    return(.Rcpp_backsolve(r=r, x=x, upper_tri = upper.tri, transpose=transpose))
  }
}

.chol2solve <- function(r, x=NULL) {
  anyMatrix <- (inherits(r,"Matrix") || inherits(x,"Matrix"))
  if (anyMatrix) warning("*M*atrix argument in call to .chol2solve()")
  if ( anyMatrix || ncol(r)< .spaMM.data$options$thr_backsolve) { 
    if (is.null(x)) x <- diag(nrow=ncol(r))
    return(backsolve(r, backsolve(r,x,transpose=TRUE)))
  } else {
    return(.Rcpp_chol2solve(r=r, x=x))
  }
}
