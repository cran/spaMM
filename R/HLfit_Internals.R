
## declared "arglist" to print a clean summary instead of very long list
## the print(summary()) is required if the "arglist" is a member (eg as.list(<call>));
## summary alone would print nothing
print.arglist <- function(x,...) {print(summary(x,...))}
##

.Dvec_times_m_Matrix <- function(Dvec, X) {
  if (inherits(X,"Matrix")) {
    return(.Dvec_times_Matrix(Dvec=Dvec, X=X))
  } else return(.Dvec_times_matrix(Dvec=Dvec, X=X))
}

.Dvec_times_matrix <- function(Dvec, X) { ## for *m*atrix input
  if (nrow(X)!=length(Dvec)) {
    stop("nrow(X)!=length(Dvec) ") ## fatal error for eigen code...
  } else if (ncol(X)==0L) {
    return(X)
  } else .Rcpp_sweepZ1W(X,Dvec) ## Rcpp (fast !) version of sweep ( MARGIN=1L ) which is also X * Dvec
}

.m_Matrix_times_Dvec <- function(X, Dvec) {
  if (inherits(X,"Matrix")) {
    return(.Matrix_times_Dvec(X=X, Dvec=Dvec))
  } else return(sweep(X, MARGIN=2L, Dvec, `*`))
}

.Matrix_times_Dvec <- function(X,Dvec) {
  if (inherits(X,"ddiMatrix")) {
    if (X@diag=="U") { ## diag + unitary => identity matrix
      X <- Diagonal(x=Dvec)
    } else X@x <- X@x * Dvec ## raw diag matrix
  } else if (inherits(X,c("dgCMatrix","dtCMatrix","dsCMatrix"))) {
    col_indices <- rep(1L:(ncol(X)),diff(X@p))
    X@x <- X@x * Dvec[col_indices]    
    ## a triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
    if ( methods::.hasSlot(X, "diag") && X@diag=="U") Matrix::diag(X) <- Dvec
  } else {
    warning("inefficient code in .make_Xscal or .Matrix_times_Dvec") ## eg dgeMatrix: dense matrix in the S4 Matrix representation
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
    } else if (inherits(X,c("dgCMatrix","dtCMatrix","dsCMatrix"))) {
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
  } else if (inherits(X,c("dgCMatrix","dtCMatrix","dsCMatrix"))) {
    which_i_affected_rows <- X@i>(n_u_h-1L)
    col_indices <- rep(1L:(ncol(X)),diff(X@p)) ##• from 1!
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
  } else if (inherits(X,c("dgCMatrix","dtCMatrix","dsCMatrix"))) {
    X@x <- X@x * Dvec[X@i+1L] 
    ## a triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
    if ( methods::.hasSlot(X, "diag") && X@diag=="U") Matrix::diag(X) <- Dvec
  } else if (inherits(X,c("dgeMatrix"))) { ## dense Matrix
    X <- X * Dvec 
  } else {
    # dgeMatrix occurs in Matern fitted by sparse coor methods (occurs if other ranef sparsifies stuff !?)
    warning("Possibly inefficient code in .Dvec_times_Matrix") ## eg dgeMatrix: dense matrix in the S4 Matrix representation
    ## warning if a dgeMatrix has been created = possibly inefficient code
    ## Other matrix formats are possibly not excluded
    ## if it is really a dgeMatrix then Dvec * X appears correct. (FIXME double check and implement ?)
    X <- Diagonal(x=Dvec) %*% X
  } 
  return(X)
}

.Dvec_times_Matrix_lower_block <- function(Dvec,X,min_row) {
  if (inherits(X,"ddiMatrix")) {
    if (X@diag=="U") { ## diag + unitary => identity matrix
      X <- Diagonal(x=c(rep(1,min_row),Dvec))
    } else X@x <- X@x * c(rep(1,min_row),Dvec) ## raw diag matrix
  } else if (inherits(X,c("dgCMatrix","dtCMatrix"))) {
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



## the following fns try to keep the input class in output, but are called with dense matrices (except irst tested case).
# les Matrix::(t)crossprod  paraissent (parfois au moins) remarquablement inefficaces !!
# idem pour Diagonal()
.ZWZtwrapper <- function(ZAL,w) { ## in .AUGI0_ZX_sparsePrecision; in calc_asDmLR_invV_from_fitobject...
  if (inherits(ZAL,"Matrix")) {
    if (inherits(ZAL,"ddiMatrix")) { ## if diagonal
      if ( ZAL@diag=="U") {
        return(Diagonal(x=w)) ## diagonal and unitary => identity
      } else return(Diagonal(x = w * diag(ZAL)^2)) ## this case may never occur
    } else {
      ZW <- .Matrix_times_Dvec(ZAL,w)
      return(Matrix::tcrossprod( ZW ,ZAL)) ## (Z W) %*% Zt
    }
  } else return(.ZWZt(ZAL,w))
}

.ZtWZwrapper <- function(ZAL,w) { ## used in several contexts
  # if (ncol(ZAL)==0L) {
  #   stop(".ZtWZwrapper called with ncol(ZAL)=0") ## temporary devel code since all calls are in principle protected 
  # } else 
  if (inherits(ZAL,"Matrix")) {
    if (inherits(ZAL,"ddiMatrix")) { ## if diagonal
      if ( ZAL@diag=="U") {
        return(Diagonal(x=w)) ## diagonal and unitary => identity
      } else return(Diagonal(x = w * diag(ZAL)^2)) ## this case may never occur
    } else  {
      DZAL <- ZAL
      DZAL@x <- DZAL@x * w[DZAL@i+1L] ## W Z
      ## a triangular matrix with unitary diagonal may be stored as @diag=="U" and only non-diag elements specified...
      if (  methods::.hasSlot(DZAL, "diag") &&  DZAL@diag=="U") Matrix::diag(DZAL) <- w
      return(Matrix::crossprod(ZAL, DZAL)) ## t(Z) %*% (W Z)
    }
  } else return(.ZtWZ(ZAL,w))
}

.process_ranCoefs <- function(ranCoefs_blob=NULL,ZAlist, ranCoefs, cum_n_u_h) {
  if (is.null(ranCoefs_blob)) { ## in call from .preprocess(), allows simplified user input as numeric vector
    Xi_cols <- attr(ZAlist, "Xi_cols")
    isRandomSlope <- Xi_cols>1L ## FIXME seems oK for later code but semantically sloppy, cf (X-1|id) terms
    hasRandomSlope <- any(isRandomSlope)
    if (hasRandomSlope && ! is.null(ranCoefs)) {
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
    } else return(list(isRandomSlope=isRandomSlope, is_set=(isRandomSlope & FALSE)))
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
      ##  => final return() 
    } else return(ranCoefs_blob)
  }
  ## newly_set contains some TRUE
  nrand <- length(ZAlist)
  if (is.null(LMatrices <- ranCoefs_blob$LMatrices)) LMatrices <- vector("list", nrand)
  if (is.null(lambda_est <- ranCoefs_blob$lambda_est)) lambda_est <- numeric(cum_n_u_h[nrand+1L])
  tol_ranCoefs <- 1e-8
  for (rt in which(isRandomSlope & newly_set)) {
    Xi_ncol <- Xi_cols[rt]
    # assume input is marginal variances + correlation as in user input, but ordered as in lower.tri
    diagpos <- cumsum(seq(Xi_ncol))
    ## fills first by the Corrrrr mat
    compactcovmat <- diag(nrow=Xi_ncol)
    compactcovmat[lower.tri(compactcovmat)] <- ranCoefs[[rt]][-diagpos]
    compactcovmat <- as.matrix(forceSymmetric(compactcovmat,uplo="L")) ## upper.tri does not fill symmetrically to lower.tri !
    marg_sd <- diag(x=sqrt(ranCoefs[[rt]][diagpos]))
    ## from Corrrr to cov: 
    compactcovmat <- marg_sd %*% compactcovmat %*% marg_sd
    latentL_blob <- .calc_latentL(compactcovmat,dtC=TRUE)
    ## we have a repres in terms of ZAL and of a diag matrix of variances; only the latter affects hlik computation
    longLv <- .makelong(latentL_blob$u,longsize=ncol(ZAlist[[rt]])) ## the variances are taken out in $d
    attr(longLv,"latentL_blob") <- latentL_blob ## kept for updating in next iteration, strucList, predVar, and output
    attr(longLv,"par") <- ranCoefs[[rt]]
    attr(longLv,"ranefs") <-  structure(attr(ZAlist,"exp_ranef_strings")[rt], 
                                              type= attr(ZAlist,"exp_ranef_types")[rt] )
    attr(longLv, "corr.model") <- "random-coef"
    LMatrices[[rt]] <- longLv
    blocksize <- ncol(ZAlist[[rt]])/Xi_ncol 
    u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
    lambda_est[u.range] <- rep(latentL_blob$d,rep(blocksize,Xi_ncol)) 
    lambda_est[lambda_est<tol_ranCoefs] <- tol_ranCoefs ## arbitrarily small eigenvalue is possible for corr=+/-1 even for 'large' parvec
  }
  return(list(isRandomSlope=isRandomSlope, is_set=newly_set, LMatrices=LMatrices, lambda_est=lambda_est))
}

.calcRanefPars <- function(corrEstList=NULL, ## potentially a list...
                          lev_lambda,
                          ranefEstargs,
                          ranCoefs_blob,
                          lambda.Fix, ## distinctly used for CARdispGammaGLM
                          rand.families,
                          psi_M,
                          verbose,
                          control,
                          iter ## ajustement gamma(identity...)
) {
  ## Build pseudo response for lambda GLM/HGLM
  glm_lambda <- NULL
  next_LMatrices <- prev_LMatrices <- ranefEstargs$prev_LMatrices
  if ( ! is.null(prev_LMatrices) && ! is.list(prev_LMatrices)) prev_LMatrices <- list(dummyid=prev_LMatrices)
  if ( ! is.null(corrEstList) && ! is.list(corrEstList)) corrEstList <- list(corr_est=corrEstList)
  next_corrEstList <- list()
  ranefs <- attr(ranefEstargs$ZAlist,"exp_ranef_strings")
  nrand <- length(ranefEstargs$ZAlist) ## Cf notes du 3/6/2015
  done <- rep(FALSE,nrand)
  u_h <- ranefEstargs$u_h
  cum_n_u_h <- ranefEstargs$cum_n_u_h
  resp_lambda <- matrix(0,cum_n_u_h[nrand+1L],1L)
  next_lambda_est <- numeric(length(u_h)) 
  #########################
  isRandomSlope <- ranCoefs_blob$isRandomSlope 
  if (any(isRandomSlope)) {
    is_set <- ranCoefs_blob$is_set
    #next_LMatrices[is_set] <- ranCoefs_blob$LMatrices[is_set] ## not necessary: prefilled by prev_LMatrices
    done[is_set] <- TRUE
    var_ranCoefs <- ( isRandomSlope & ! is_set ) 
    if (any(var_ranCoefs)) {
      ## handling correlation in random slope models # slmt pr gaussian ranefs, verif dans preprocess
      ranefEstargs$var_ranCoefs <- var_ranCoefs 
      LMatricesBlob <- do.call(.spaMM.data$options$covEstmethod, ranefEstargs)  
      next_corrEstList$random_coeff <- LMatricesBlob$optr_par
      next_LMatrices[var_ranCoefs] <- LMatricesBlob$updated_LMatrices[var_ranCoefs] 
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
    next_corrEstList$random_coeff <- NULL
  }
  ### next the other LMatrix models
  for (it in seq_along(prev_LMatrices) ) { ## this loop will ignore ranefs not affected by any lmatrix
    attr_lmatrix <- attributes(prev_LMatrices[[it]])
    ## find ZAlist elements affected by LMatrix element
    affected <- which(ranefs %in% attr_lmatrix[["ranefs"]] & ! done)
    ## then for each L matrix we need to select the relevant blocks of random effects
    if (length(affected)>1L) {
      stop("code needed for length(affected)>1")
    } else if (length(affected)==1L) {
      if ( ! is.null(corrEstList$corr_est) && attr_lmatrix[["corr.model"]] %in% c("adjacency")) { 
        ## the following conforms to an interface where 
        ##  next_lambda_est is a vector (heterosc) of what should go in the sigma_aug matrix
        ## consequently, the LMatrix for this model should be decomp$u, constant 
        adj_symSVD <- ranefEstargs$processed$AUGI0_ZX$envir$adj_symSVD ## may be NULL
        if (is.null(adj_symSVD)) {
          stop("is.null(adj_symSVD)")
          adj_symSVD <- attr_lmatrix[[ attr_lmatrix[["type"]] ]]  ## older conception
        }
        adjd <- adj_symSVD$adjd
        locdf <- data.frame(adjd=adjd) ## $adjd, not $d which is (1/(1-rho * $adjd)): adj, not corr
        u.range <- (cum_n_u_h[affected]+1L):cum_n_u_h[affected+1L]
        locdf$resp <- resp_lambda[u.range] <- u_h[u.range]^2
        ## here CAR allows REML contrary to the SEM CAR, hence leverages
        glm_lambda <- .calc_CARdispGammaGLM(data=locdf, lambda.Fix=lambda.Fix[affected], lev=lev_lambda[u.range],control=control)
        attr(glm_lambda,"whichrand") <- affected
        next_lambda_est[u.range] <- fitted(glm_lambda) ## prediction of heteroscedastic variances
        coeffs <- coefficients(glm_lambda)
        if (is.na(lambda.Fix[affected])) {  
          next_corrEstList$adjacency <- list(rho = - coeffs[["adjd"]]/ coeffs[["(Intercept)"]]) ## FR->FR different corr_est s'écrasent les uns les autres
        } else {
          next_corrEstList$adjacency <- list(rho = - coeffs[1]*lambda.Fix)  
        }
        done[affected] <- TRUE
      } ## Matern remains undone
    } ## else (length(affected)=0), for previously processed random slope models
  }        
  ### next the (no L matrices) or (Matern model and other fixed L matrix cases)
  for (it in which( ! done )) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    if (is.na(unique.lambda <- lambda.Fix[it])) {
      resp_lambda[u.range] <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=1) ## must give d1 in table p 989 de LeeN01
      unique.lambda <- sum(resp_lambda[u.range])/sum(1-lev_lambda[u.range]) ## NOT in linkscale 
      unique.lambda <- max(unique.lambda,1e-8) # FR->FR still corrected
      unique.lambda <- min(unique.lambda,.spaMM.data$options$maxLambda)  
      if (ranefEstargs$lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { ## Gamma(identity)
        unique.lambda <- pmin(unique.lambda,1-1/(2^(iter+1)))  ## impose lambda<1 dans ce cas 
      }
    } 
    next_lambda_est[u.range] <- rep(unique.lambda,length(u.range))
  }
  if (verbose["trace"]) { 
    ## this is cryptic but a better output would require a lot of reformatting as at the end of HLfit_body. 
    print(paste("unique(next_lambda_est)=",paste(signif(unique(next_lambda_est),4),collapse=" ")),quote=FALSE)
    print_corr_est <- unlist(next_corrEstList)
    if ( ! is.null(print_corr_est)) print(paste("corr_est=",paste(signif(print_corr_est,4),collapse=" ")),quote=FALSE)
  }
  return(list(next_LMatrices=next_LMatrices,
              resp_lambda=resp_lambda, ## for final glm...
              next_corrEstList=next_corrEstList,
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
    if (class(glm_lambda$data)=="environment") {
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
  w.resid <- structure(as.vector(GLMweights/phi_est),unique= (attr(GLMweights,"unique") && length(phi_est)==1L))
  if (ilg) {
    res[["w_resid"]] <- w.resid
    return(res)
  } else {
    return(w.resid)
  }
}

.calc_H_global_scale <- function(w.resid) {
  if (is.list(w.resid)) {
    return(exp( - mean(log(w.resid$w_resid))))
  } else return(exp( - mean(log(w.resid)))) ## generalization of exp(mean(log(phi_est)))
}

.calc_weight_X <- function(w.resid, H_global_scale) {
  if (is.list(w.resid)) {
    return(sqrt(H_global_scale*w.resid$w_resid))
  } else return(sqrt(H_global_scale*w.resid))
}

## spaMM_Gamma() fixes Gamma()$dev.resids(1e10+2,1e10,1) is < 0
# dev.resids() must be >0 for computation deviance_residual in fitting Gamma GLMM, and alos for $aic() computation.
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
  variance <- function(mu) mu^2
  validmu <- function(mu) all(mu > 0)
  dev.resids <- function(y, mu, wt) {
    if (any(mu < 0)) return(Inf) ## 2015/04/27; maybe not useful
    ## otherwise see deviance.gamma function locally defined in statmod::glmgam.fit
    dev_res <- -2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
    dev_res[dev_res<.Machine$double.eps] <- .Machine$double.eps ##  ## FR: added this
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
    shape <- gamma.shape(object)$alpha * wts 
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
           ## uses C code:
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
           -family$aic(y=y, mu=NB_shape/(exp(-theta)-1), wt=1)/2 ## handles truncation from given untruncated mu
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

.muetafn <- function(eta,BinomialDen,processed) { ## note outer var BinomialDen 
  family <- processed$family
  ## patches for borderline eta's
  if (family$link =="log") {
    eta[eta>30] <-30 ## 100 -> mu = 2.688117e+43 ; 30 -> 1.068647e+13
  } else if (family$family == "COMPoisson" && family$link =="loglambda") {
    etamax <- 30*environment(family$aic)$nu
    eta[eta>etamax] <- etamax ## using log(mu) ~ eta/nu for large nu
  } else if (family$link=="inverse" && family$family=="Gamma") {
    etamax <- sqrt(.Machine$double.eps)
    eta[eta>etamax] <- etamax ## both eta and mu must be >0
  }
  mu <- family$linkinv(eta) ## linkinv(eta) is FREQS for binomial, COUNTS for poisson...
  if (family$link %in% c("logit","probit","cloglog","cauchit")) {
    mu[mu > (1-1e-12)] <- (1-1e-12)
    mu[mu < (1e-12)] <- (1e-12)
  }
  dmudeta <- family$mu.eta(eta) ## aberrant at hoc code for cloglog 'elsewhere'...
  Vmu <- family$variance(mu) 
  if (family$family=="binomial") {
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
    GLMweights <- eval(processed$prior.weights) * dmudeta^2 /Vmu ## must be O(n) in binomial cases
    attr(GLMweights,"unique") <- FALSE ## might actually be true sometimes
  }
  if (identical(family$zero_truncated,TRUE)) { 
    if (family$family=="poisson") { 
      ## D[Log[1 - E^-E^theta], {theta, 2}] /. {theta -> Log[mu]} // Simplify
      expmu <- exp(mu)
      p0 <- 1/expmu
      dlogMthdth <- -mu/(1-expmu) ## useful to correct rhs
      d2logMthdth2 <- -(1 + expmu * (mu-1))* mu/((expmu-1)^2)
      exp2mu <- expmu^2
      d3logMthdth3 <- mu * (1 + 3 * mu *(expmu - exp2mu) + (mu^2) *(expmu + exp2mu) + exp2mu - 2 * expmu )/(expmu-1)^3
    } else if (family$family=="negbin") { 
      ## D[Log[1 - E^-E^theta], {theta, 2}] /. {theta -> Log[mu]} // Simplify
      shape <- .get_family_par(family,"shape")
      p0 <- (shape/(mu+shape))^shape ## (1-p)^r
      dlogMthdth <- (mu * p0)/(1-p0)
      d2logMthdth2 <- -((mu * p0 *(shape *(-1 + p0) + mu *(-1 + shape + p0)))/(shape *(-1 + p0)^2))
      d3logMthdth3 <- -((mu * p0 * (shape^2 *(-1 + p0)^2 + 3 * mu * shape *(-1 + p0) * (-1 + shape + p0) + 
                                      mu^2 *(3 *shape *(-1 + p0) + 2 *(-1 + p0)^2 + shape^2 *(1 + p0))))/(
                                        shape^2 *(-1 + p0)^3))
    } 
    truncGLMweights <- GLMweights*(1+d2logMthdth2/Vmu) 
    WU_WT <- GLMweights/truncGLMweights 
    return(list(mu=mu, dmudeta=dmudeta, p0=p0,
                GLMweights=list(truncGLMweights=truncGLMweights,WU_WT=WU_WT,dlogMthdth=dlogMthdth, 
                                d2logMthdth2=d2logMthdth2, d3logMthdth3=d3logMthdth3)))
  } else return(list(mu=mu,dmudeta=dmudeta,GLMweights=GLMweights))
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
  blob <- vector("list",nrand)
  for (it in seq_len(nrand)) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    blob[[it]] <- .updateWranef(rand.family=rand.families[[it]],lambda[u.range],u_h[u.range],v_h[u.range])
  }
  w.ranef <- unlist(lapply(blob,getElement, name="w.ranef"))
  w.ranef[w.ranef>1e10] <- 1e10 ## patch useful to avoid singular d2hdv2 in PLoG model
  dlogWran_dv_h <- unlist(lapply(blob,getElement, name="dlogWran_dv_h"))
  dvdu <- unlist(lapply(blob, getElement, name="dvdu"))
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
         probit = -eta*dnorm(eta) * BinomialDen,
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
      ## numerator is D[D[1/(1 + E^-\[Eta]), \[Eta]] /. {E^-\[Eta]->(1-\[Mu])/\[Mu]} ,\[Mu]]=1-2 mu 
      if (calcCoef1) coef1 <-(1-2*mu/BinomialDen)/dmudeta  
      dlW_deta <-(1-2*mu/BinomialDen)  
    } else if (family$family=="Gamma") { ## link= "inverse" !
      ## numerator is D[D[-1/\[Eta], \[Eta]] /. {\[Eta] -> -1/\[Mu]}, \[Mu]] =2 mu 
      if (calcCoef1) coef1 <- 2*mu /dmudeta
      dlW_deta <- 2*mu
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
    muFREQS <- mu/BinomialDen
    dlW_deta <- -2*eta - dnorm(eta)*(1-2*muFREQS)/(muFREQS*(1-muFREQS))
    if (calcCoef1) {
      coef1 <- dlW_deta *(muFREQS*(1-muFREQS))/ (BinomialDen * dnorm(eta)^2) 
      coef1[coef1>1e100] <- 1e100
      coef1[coef1< -1e100] <- -1e100
    }
  } else if (family$family=="Gamma" && family$link=="log") { ## ad hoc non canonical case 
    if (calcCoef1) coef1 <- rep(0L,length(mu))
    dlW_deta <- rep(0L,length(mu)) ## because they both involve dW.resid/dmu= 0
  } else { ## "negbin" only called with non-canonical links always end here 
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

.safesolve_qr_matrix <- function(qr.a,B,silent=TRUE,stop.on.error=TRUE) { ## solve.qr with fall-back; qr.a should be a qr object, B a matrix
  ## there was a 'Matrix' subcode prior to 10/03/2013
  res <- try(solve.qr(qr.a,B),silent=silent)
  if (inherits(res,"try-error")) { ##FR->FR sb systematique qd phi -> 0; slow step.
    # pivI <- sort.list(qr.a$pivot)  ## inverse perm such as pivI[$pivot]=$pivot[pivI]= identity
    #solveA <- try(solve(qr.R(qr.a)[,pivI])%*%t(qr.Q(qr.a)),silent=silent)  
    ## solveA is solve(<original 'a' matrix>) using the QR decomp... but this may work when solve.qr fails !
    solveA <- try(backsolve(qr.R(qr.a),t(qr.Q(qr.a))[qr.a$pivot,]))
    if (inherits(solveA,"try-error")) {
      if (stop.on.error) {
        stop("backsolve() failed in .safesolve_qr_matrix().") ## perhaps recover A by qr.X and solve(A) ?
      } else return(solveA) ## passes control to calling function
    } else res <- solveA %*% B  
  }
  return(res)
}

.safesolve_qr_vector <- function(qr.a,b,silent=TRUE,stop.on.error=TRUE) { ## solve.qr with fall-back; qr.a should be a qr object, b must be a vector
  if (class(qr.a)=="sparseQR") { ## pas de 'essai' en var locale !
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
  if (inherits(x,"Matrix")) {
    if (inherits(x,"ddiMatrix") ) return(x@diag=="U")
    if (! isDiagonal( x ) ) return( FALSE )
    ## hence now diagonal:
    return(max(abs(range(diag(x))-1))<tol)
  } else {
    if (matrixcheck) {
      return(ncol(x)==nrow(x) && max(abs(x- diag(ncol(x))))<tol)
    } else return(FALSE)  
  }
}

.solveWrap_vector <- function(A,b,...) { ## chol versions not used...
  if (inherits(A,"diagonalMatrix")) return(b/diag(A)) 
  if (inherits(A,"RcppChol")) { ## no pivoting; A$L is .L.ower tri
    return(forwardsolve(A$L,forwardsolve(A$L,b),transpose=TRUE)) ## inv(LLt).b=inv(Lt).(invL.b)
  }
  if (inherits(A,"cholL_LLt")) { ## as results from t(base::chol) with default pivot=FALSE; A is .L.ower tri 
    return(forwardsolve(A,forwardsolve(A,b),transpose=TRUE)) 
  }
  if (inherits(A,"Rcpp_sparseQR")) { ## PIVOTING
    if (length(b)==0L) {return(A$Q_ap)} ## appropriate dimensions
    dim(b) <- c(1,length(b)) ## conversion to "1-rwo matrix" without copy contrary to t(b)
    b <- b %*% (A$Q_ap)
    dim(b) <- c(length(b),1) ## transposition without copy    
    solved <- try(solve(A$R_ap,b)[A$pivI,],silent=TRUE)
    return(solved) ## gives control to calling function 
  } 
  ## next line should become obsolete ?
  if (inherits(A,"sparseQR")) return(Matrix::solve(A,b)) ## as produced by Matrix::qr; return value not documented, but sparse storage is used 
  ## all other cases   
  .safesolve_qr_vector(A,b,...)
}

.solveWrap_matrix <- function(A,B,...) { ## chol versions not used...
  if (inherits(A,"diagonalMatrix")) return(B/diag(A)) ## works if A is the matrix, not its diagonal...
  if (inherits(A,"RcppChol")) { ## no pivoting; A$L is .L.ower tri
    return(forwardsolve(A$L,forwardsolve(A$L,B),transpose=TRUE)) 
  }
  ## *!* FR->FR confusions ahead
  # class is never "cholL_LLt"; (identical(attr(A,"type"),"cholL_LLt"))  would be a more meaningful test
  # But the code would be wrong anyway, whe A is an LMatrix, not a list with an $L element.
  # => Distinction type de decomp / classe de matrice à garder
  if (inherits(A,"cholL_LLt")) { ## as results from t(base::chol) with default pivot=FALSE; A is .L.ower tri 
    return(forwardsolve(A,forwardsolve(A$L,B),transpose=TRUE)) 
  }
  if (inherits(A,"Rcpp_sparseQR")) { ## PIVOTING
    ## Q and R need not be sparse (even if stored as sparse matrices), can still be sparse in simple aplications
    if (.is_identity(B,matrixcheck=FALSE)) {
      solved <- try(solve(A$R_ap,t(A$Q_ap))[A$pivI,],silent=TRUE)
#    } else if (inherits(B,"sparseMatrix")) { 
#      solved <- try(solve(as(A$R_ap,"dtCMatrix"),t(A$Q_ap) %*% B,sparse=TRUE)[A$pivI,],silent=TRUE)
    } else solved <- try(backsolve(A$R_ap,crossprod(A$Q_ap, B))[A$pivI,],silent=TRUE)    
    return(solved) ## gives control to calling function 
  }
  if (inherits(A,"sparseQR")) { ## as produced by Matrix::qr; return value not documented (!), but sparse storage is used
    return(suppressMessages(solve(A,B,sparse=inherits(B,"sparseMatrix")))) 
  }
  ## all other cases
  .safesolve_qr_matrix(A,B,...)
}

# returns QR or diagonalMatrix, so the calling code must handle that
.QRwrap_but_diag <- function(mat, ## now M or m
                   useEigen=TRUE ## TRUE seems important for large square matrices such as d2hdv2
                   ## FALSE is faster for wAugX for linear solving
                   ) {
  if (inherits(mat,"diagonalMatrix")) { 
    return(mat)
  } else if (inherits(mat,"Matrix") && ncol(mat)<=nrow(mat)) { ## ncol may be > nrow in get_predVar -> ...
    # ... -> .get_logdispObject -> .QRwrap_but_diag(ZAL [cf example cov1 <- get_predVar(fitobject,newdata=moregroups,...) ]
    QR <- qr(mat) ## i.e. Matrix::qr
  } else {
    QR <- qr(mat)
  }
  return(QR)
} 

sym_eigen <- function(X) {
  if (inherits(X,"sparseMatrix")) {
    X <- as.matrix(X) ## dumb, but this is what RSpectra:::eigs_real_sym() does when full eigen is required. 
  }
  if (is.integer(X)) X <- 1.0*X
  return(.selfAdjointSolverCpp(X))
}

.Cholwrap <- function(mat) {
  if (inherits(mat,"diagonalMatrix")) { 
    chol <- diag(ncol(mat))
    class(chol) <- c("RcppChol",class(chol)) ## FR->FR alternatively return a diagonalMatrix ?  
  } else if (inherits(mat,"Matrix")) {
    chol <- t(Matrix::chol(mat))
  } else if (.spaMM.data$options$USEEIGEN) {
    chol <- .RcppChol(mat) ##
    if ( chol$Status==1L) { 
      return(chol$L) 
    } else stop("chol$Status !=1L") ## best used in combination with try()
  } else chol <- t(chol(mat)) ## !! this is the R matrix; pivot=FALSE by default
} 

.get_qr <- function(mat,provide=TRUE) {
  envir <- attr(mat,"envir") ## environment
  envir$callcount <- envir$callcount+1L
  if (is.null(envir$qr)) { 
    if (provide) envir$qr <- .QRwrap_but_diag(mat)
  } 
  return(envir$qr)
}


#LogAbsDetWrap <- function(...) .LogAbsDetWrap(...) ## 

.LogAbsDetWrap <- function(mat,logfac=0,provide.qr=FALSE) { ## M or m
  if (ncol(mat)==0) return(0) ## GLM fitted by ML: d2hdbv2 is 0 X 0 matrix 
  # un piege est que mat/(2*pi) conserve les attributes de mat (telle qu'une décomp QR de mat...)
  # il nefaut  dont pas demander LogAbsDetWrap(mat/(2*pi))
  if ( ! is.null(envir <- attr(mat,"envir"))) {
    qrmat <- .get_qr(mat,provide=provide.qr) ##  LogAbsDetWrap' own provide.qr=FALSE
  } else qrmat <- NULL
  if ( ! is.null(qrmat)) { ## 
    if (inherits(qrmat,"qr")) {
      lad <- sum(log(abs(diag(qr.R(qrmat)))))
    } else if (inherits(qrmat,"sparseQR")) { ## if d2hdv2 is a Matrix
      lad <- sum(log(abs(diag(qrmat@R)))) # _@_ ... Matrix::qr.R() serait surement plus recommande 
    } else if (inherits(qrmat,"diagonalMatrix")) { 
      lad <- sum(log(abs(diag(qrmat))))  
    }
  } else if (inherits(mat,"Matrix")) {
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

.tcrossprod <-  function(x,y=NULL) {
  if (is.null(x)) return(NULL) ## allows lapply(,.tcrossprod) on a listof (matrix or NULL)
  if (inherits(x,"Matrix") || inherits(y,"Matrix")) {
    if (is.null(y)) {
      return(Matrix::tcrossprod(x))
    } else return(Matrix::tcrossprod(x,y))
  } else {
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

.crossprod <- function(x,y=NULL) {
  if (is.null(x)) return(NULL) ## allows lapply(,.tcrossprod) on a listof (matrix or NULL)
  if (inherits(x,"Matrix") || inherits(y,"Matrix")) {
    if (is.null(y)) {
      return(Matrix::crossprod(x))
    } else return(suppressMessages(Matrix::crossprod(x,y))) ## suppressMessages for case (Matrix,matrix) or (Matrix,dgeMatrix)
  } else {
    resu <- .crossprodCpp(x,y)
    if (is.null(y)) {
      colnames(resu) <- rownames(resu) <- colnames(x)
    } else {
      rownames(resu) <- colnames(x)
      colnames(resu) <- colnames(y)
    }
    return(resu)
  }
}

.get_beta_w_cov <- function(res) {
  if (is.null(beta_w_cov <- res$envir$beta_w_cov)) { 
    beta_cov <- .get_beta_cov(res) ## beta_v_cov needed
    beta_w_cov <- attr(beta_cov,"beta_v_cov")
    invL <- .get_invL(res) ## correlation matrix of ranefs is solve((t(invL)%*%(invL)))
    # invL is currently a single matrix for allranefs. de facto a block matrix when several ranefs
    if ( ! is.null(invL)) {
      pforpv <- ncol(beta_cov)
      v.range <- pforpv+seq(ncol(invL))
      ## A function for symmetric sandwich product is missing. 
      beta_w_cov[v.range,] <- .crossprod(invL, beta_w_cov[v.range,])
      beta_w_cov[,v.range] <- beta_w_cov[,v.range] %*% invL # implies invL' %*% . %*% invL on the v.range,v.range block
      beta_w_cov <- Matrix::symmpart(beta_w_cov)
      eigenvalues <- eigen(beta_w_cov,only.values = TRUE)$values
      # beta_w_cov[v.range,-(v.range)] <- .crossprod(invL, beta_w_cov[v.range,-(v.range)])
      # beta_w_cov[-v.range,v.range] <- t(beta_w_cov[v.range,-(v.range)]) ## usually this block is small, so little speed gain
      # beta_w_cov[v.range,(v.range)] <- .crossprod(invL,beta_w_cov[v.range,(v.range)] %*% invL) ## still not enforcing symmetry
      # eigenvalues <- eigen(beta_w_cov,only.values = TRUE)$values
      # if ( inherits(eigenvalues,"complex")) { ## matrix perceived as asymmetric
      #   beta_w_cov <- Matrix::symmpart(beta_w_cov)
      #   eigenvalues <- eigen(beta_w_cov,only.values = TRUE,symmetric = TRUE)$values
      # }
      attr(beta_w_cov,"min_eigen") <- min(eigenvalues)
    }
    res$envir$beta_w_cov <- beta_w_cov
  } 
  return(beta_w_cov)
}

get_ZALMatrix <- function(object,as_matrix) {
  if ( ! missing(as_matrix)) stop("'as_matrix' is deprecated")
  if (length(ZAlist <- object$ZAlist)) { ## ou tester if (object$models[["eta"]]=="etaGLM")
    if (is.null(object$envir$ZALMatrix)) {
      object$envir$ZALMatrix <- .compute_ZAL(XMatrix=object$strucList, ZAlist=ZAlist,as_matrix=FALSE) 
    }
    return(object$envir$ZALMatrix)
  } else return(NULL) 
}


.eval_as_mat_arg <- function(object) { 
  (
    object$HL[1L]=="SEM" || ## SEM code does not yet handle sparse as it uses a dense Sig matrix
    ! identical(object$QRmethod,"sparse") ## => conversion to matrix if object$QRmethod is NULL
  )
}

.get_invColdoldList <- function(res,regul.threshold=1e-7) {
  ## the validity of this fn is tested by checking that Evar (in predVar) is null i nthe case where newdata=ori data
  # one needs to force the computation of Evar for that test (newdata && [Neednewnew= variances$cov])
  if (is.null(invColdoldList <- res$envir$invColdoldList)) { 
    ## returns a list of inv(Corr) from the LMatrix
    strucList <- res$strucList
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
          affecteds <- which(ranefs %in% attr(lmatrix,"ranefs"))
          ## end of designL.from.corr implies either type is cholL_LLt or we have decomp $u and $d 
          type <-  attr(lmatrix,"type")
          invCoo <- NULL
          if ( ! is.null(latentL_blob <- attr(lmatrix,"latentL_blob"))) { ## from .process_ranCoefs
            design_u <- latentL_blob$u 
            if (is.null(design_u)) {
              invCoo <- stop("code missing here (1)") ## there may be a $compactchol_Q in spprec case.
            } else invCoo <- .makelong(solve(.tcrossprod(design_u)),longsize=ncol(lmatrix))
          } else if (type == "from_Q_CHMfactor")  {
            invCoo <- tcrossprod(as(attr(lmatrix,"Q_CHMfactor"),"sparseMatrix")) ## FIXME more direct syntax to recover the orig matrix ?
          } else if (type == "cholL_LLt")  {
            Rmatrix <- t(lmatrix)
          } else { ## Rcpp's symSVD, or R's eigen() => LDL (also possible bad use of R's svd, not normally used)
            condnum <- kappa(lmatrix,norm="1")
            if (condnum<1/regul.threshold) {
              decomp <- attr(lmatrix,attr(lmatrix,"type")) ## of corr matrix !
              invCoo <-  try(.ZWZt(decomp$u,1/decomp$d),silent=TRUE)
              if (inherits(invCoo,"try-error")) invCoo <- NULL
            }
            if (is.null(invCoo)) Rmatrix <- qr.R(qr(t(lmatrix))) 
          }
          if (is.null(invCoo)){ ## see comments on chol2inv in .get_invL()
            singular <- which(abs(diag(Rmatrix))<regul.threshold) 
            if (length(singular)) {
              if (spaMM.getOption("wRegularization")) warning("regularization required.")
              nc <- ncol(Rmatrix)
              diagPos <- seq.int(1L,nc^2,nc+1L)[singular]
              Rmatrix[diagPos] <- sign(Rmatrix[diagPos])* regul.threshold
            }
            invCoo <- chol2inv(Rmatrix) ## 
          }
          for (aff in affecteds) invColdoldList[[aff]] <- invCoo
        } 
      }
      res$envir$invColdoldList <- invColdoldList
    } else return(NULL)
  } 
  return(invColdoldList)
}

.get_info_crits <- function(object) {
  if (is.null(info_crits <- object$envir$info_crits)) { 
    pforpv <- object$dfs[["pforpv"]]
    p_phi <- object$dfs[["p_phi"]]
    p_lambda <- object$dfs[["p_lambda"]]
    APHLs <- object$APHLs
    w.resid <- object$w.resid
    predictor <- object$predictor
    info_crits <- list()
    if  ( ! is.null(resid_fit <- object$resid_fit)) { ## indicates a phiHGLM: we need to get info from it
      # input p_phi (above) is typically set to NA, and will be ignored
      resid_fit <- object$resid_fit
      info_crits_phi <- .get_info_crits(resid_fit)
      phi_pd <- length(resid_fit$y)-info_crits_phi$GoFdf
      p_phi <- phi_pd+sum(resid_fit$dfs) ## all df's absorbed by the phi model
    }
    names_est_ranefPars <- (unlist(.get_methods_disp(object)))  
    p_GLM_family <- length(intersect(names_est_ranefPars,c("NB_shape","NU_COMP")))
    p_phi <- p_phi+p_GLM_family ## effectively part of the model for residual error structure
    # poisson-Gamma and negbin should have similar similar mAIC => NB_shape as one df or lambda as one df   
    forAIC <- APHLs
    if (object$models[[1]]=="etaHGLM") {
      if (object$HL[1]=="SEM") {
        forAIC <- list(p_v=APHLs$logLapp,p_bv=APHLs$logLapp,clik=APHLs$clik)
      } 
      ZAL <- .compute_ZAL(XMatrix=object$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg(object)) 
      d2hdv2 <- .calcD2hDv2(ZAL,w.resid,object$w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
      # if standard ML: there is an REMLformula ~ 0 (or with ranefs ?); processed$X.Re is 0-col matrix
      # if standard REML: REMLformula is NULL: $X.Re is X.pv, processed$X.Re is NULL
      # non standard REML: other REMLformula: $X.Re and processed$X.Re identical, and may take essentially any value
      # if (identical(attr(object$REMLformula,"isML"),TRUE)) {
      #   Md2hdbv2 <- - d2hdv2 
      #   Md2clikdbv2 <-  as.matrix(.ZtWZwrapper(ZAL,w.resid))
      # } else {
      #   ## REML standard || REML non standard
      #   X.Re <- object$distinctX.Re ## null if not distinct from X.pv
      #   if (is.null(X.Re)) X.Re <- object$X.pv ## standard REML
      #   ## diff de d2hdbv2 slmt dans dernier bloc (-> computation pd)
      #   hessnondiag <- .crossprod(ZAL, sweep(X.Re, MARGIN = 1, w.resid, `*`))
      #   Md2hdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
      #                                cbind2(hessnondiag, - d2hdv2))) 
      #   Md2clikdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.Re,w.resid), t(hessnondiag)),
      #                                   cbind2(hessnondiag, .ZtWZwrapper(ZAL,w.resid))))            
      # }
      X.pv <- object$X.pv
      if ( ncol(X.pv)>0 ) { ## the projection matrix for the response always includes X even for REML!
        hessnondiag <- .crossprod(ZAL, sweep(X.pv, MARGIN = 1, w.resid, `*`))
        Md2hdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.pv,w.resid), t(hessnondiag)),
                                     cbind2(hessnondiag, - d2hdv2))) 
        Md2clikdbv2 <- as.matrix(rbind2(cbind2(.ZtWZwrapper(X.pv,w.resid), t(hessnondiag)),
                                        cbind2(hessnondiag, .ZtWZwrapper(ZAL,w.resid))))            
      } else {
        Md2hdbv2 <- - d2hdv2 
        Md2clikdbv2 <-  as.matrix(.ZtWZwrapper(ZAL,w.resid))
      }
      p_corrPars <- length(intersect(names_est_ranefPars,names(object$corrPars)))
      info_crits$mAIC <- -2*forAIC$p_v + 2 *(pforpv+p_lambda+p_corrPars+p_phi)
      info_crits$dAIC <- -2*forAIC$p_bv + 2 * (p_lambda+p_phi+p_corrPars) ## HaLM07 (eq 10) focussed for dispersion params
      #                                                                             including the rho param of an AR model
      eigvals <- eigen(Md2hdbv2/(2*pi),symmetric=TRUE,only.values = TRUE)$values
      eigvals[eigvals<1e-12] <- 1e-12
      if (min(eigvals)>1e-11) {
        qr.Md2hdbv2 <- .QRwrap_but_diag(Md2hdbv2)
        ## dans un LMM avec estimation ML, pd = sum(lev_phi), mais pas de simplif plus generale 
        pd <- sum(diag(.solveWrap_matrix(qr.Md2hdbv2,Md2clikdbv2,stop.on.error=FALSE)))
        if (inherits(pd,"try-error")) {
          warning("Computation of cAIC/GoF df's failed because the 'd2hdbv2' matrix appears singular")
          pd <- NA
        }
      } else pd <- Inf
      info_crits$GoFdf <- length(object$y) - pd ## <- nobs minus # df absorbed in inference of ranefs
      ## eqs 4,7 in HaLM07
      info_crits$cAIC <- -2*forAIC$clik + 2*(pd+p_phi) ## no p_lambda !
      # print(c(pd,p_phi))
      # but Yu and Yau then suggest caicc <- -2*clik + ... where ... involves d2h/db d[disp params] and d2h/d[disp params]2
    } else { ## fixed effect model
      info_crits$mAIC <- -2*forAIC$p_v+2*(pforpv+p_phi) 
    }
    object$envir$info_crits <- info_crits
  } 
  return(info_crits)
}

.calcD2hDv2 <- function(ZAL,w.resid,w.ranef) { ## FR->FR review usagesof this function
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
    d2hdv2 <- Matrix::crossprod(x=ZAL,y= Diagonal(x= - w.resid) %*% ZAL)    
    nc <- ncol(d2hdv2)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    d2hdv2[diagPos] <- d2hdv2[diagPos] - w.ranef 
  } else d2hdv2 <- .Rcpp_d2hdv2(ZAL,w.resid,w.ranef) 
  #cat("\n new d2hdv2")
  structure(d2hdv2,envir=list2env(list(tag="d2hdv2",callcount=0L),parent=environment(HLfit_body)))
}


.get_logdispObject <- function(res) {
  if (is.null(res$envir$logdispObject) && res$models[["eta"]]=="etaHGLM" ) { 
    asDmLR_invV <- .calc_asDmLR_invV_from_fitobject(res)
    dvdloglamMat <- res$envir$dvdloglamMat
    dvdloglamMat_needed <- ( is.null(dvdloglamMat) && 
# (comment this => allows random slope)  all(unlist(attr(res$ZAlist,"namesTerms"))=="(Intercept)") && ## (1|.) or CAR or Matern
                               any(res$lambda.object$type!="fixed") ) ## some lambda params were estimated
    dvdlogphiMat <- res$envir$dvdlogphiMat
    dvdlogphiMat_needed <- (is.null(dvdlogphiMat) && 
                              res$models[["phi"]]=="phiScal") ## cf comment in calc_logdisp_cov
    if (dvdloglamMat_needed || dvdlogphiMat_needed) {
      ZAL <- get_ZALMatrix(res) ## should later simplify as =(res$QRmethod=="dense")) FIXME       
      d2hdv2 <- .calcD2hDv2(ZAL,res$w.resid,res$w.ranef) 
    }
    if (dvdloglamMat_needed) { 
      cum_n_u_h <- attr(res$ranef,"cum_n_u_h")
      psi_M <- rep(attr(res$rand.families,"unique.psi_M"),diff(cum_n_u_h))
      dlogfthdth <- (psi_M - res$ranef)/res$lambda.object$lambda_est ## the d log density of th(u)
      dvdloglamMat <- .calc_dvdloglamMat_new(dlogfthdth=dlogfthdth,
                                             cum_n_u_h=cum_n_u_h,
                                             lcrandfamfam=attr(res$rand.families,"lcrandfamfam"),
                                             rand.families=res$rand.families,
                                             u_h=res$ranef,d2hdv2=d2hdv2,stop.on.error=TRUE)
    }
    if (dvdlogphiMat_needed) {
      muetablob <- res$muetablob
      dh0deta <- ( res$w.resid *(res$y-muetablob$mu)/muetablob$dmudeta ) ## (soit Bin -> phi fixe=1, soit BinomialDen=1)
      dvdlogphiMat  <- .calc_dvdlogphiMat_new(dh0deta=dh0deta, ZAL=ZAL,
                                              d2hdv2=d2hdv2, stop.on.error=TRUE)
    }
    ## This call gets args (except res) from the envir of locfn def'd below:
    res$envir$logdispObject <- .calc_logdisp_cov(res, dvdloglamMat=dvdloglamMat, 
                                       dvdlogphiMat=dvdlogphiMat, asDmLR_invV=asDmLR_invV)
  } 
  return(res$envir$logdispObject)
}

## slow computation the one time .get_logdispObject() is called, for variances$disp (no need to store the result in an $envir)
.calc_asDmLR_invV_from_fitobject <- function(object) { ## used by .get_logdispObject
  ## Store Compute inv[{precmat=inv(L invWranef Lt)} +ZtWZ] as two matrix nXr and rXn rather than their nXn product
  if ("AUGI0_ZX_sparsePrecision" %in% object$MME_method) {
    RES <- list(ZAfix = .post_process_ZALlist(object$ZAlist, as_matrix=.eval_as_mat_arg(object))  )
    ## code clearly related to .Sigsolve_sparsePrecision() algorithm:
    RES$n_x_r <- t(object$envir$ZtW)
    RES$r_x_n <- solve(object$envir$G_CHMfactor, object$envir$ZtW)
    return(RES)
  } else {
    ZAL <- .compute_ZAL(XMatrix=object$strucList, ZAlist=object$ZAlist,as_matrix=.eval_as_mat_arg(object)) 
    if (inherits(ZAL,"diagonalMatrix")) {
      # direct but ad hoc ; commented code below is equivalent with probably little overhead (but less readable)
      w.resid <- object$w.resid
      return(list(n_x_r=Diagonal(x = w.resid),
                  r_x_n=Diagonal(x = w.resid/(w.resid + object$w.ranef/diag(ZAL)^2)))) 
    } else {
      ZAfix <- .post_process_ZALlist(object$ZAlist,as_matrix=FALSE)
      wrZ <- .Dvec_times_m_Matrix(object$w.resid, ZAfix) # suppressMessages(sweep(t(ZA), 2L, w.resid,`*`)) 
      #if (is.null(object$envir$invG)) { # it's one-time computation so no use to store invG
        invL <- .get_invL(object)
        precmat <- .ZtWZwrapper(invL,object$w.ranef)
        ZtwrZ <- .crossprod(ZAfix, wrZ) 
        inv2 <- suppressMessages(precmat+ZtwrZ) ## suppress signature message
        invG <- try(solve(inv2),silent=TRUE)
        if (inherits(invG,"try-error") || anyNA(invG)) {
          invG <- ginv(as.matrix(inv2)) ## FIXME quick patch at least
        } else invG <- invG
      #}
      ## avoid formation of a large nxn matrix:
      return(list(n_x_r=wrZ, r_x_n=.tcrossprod(invG, wrZ)))
    } ## removed primitive QR version at version 2.2.138
  }
}

.get_glm_phi <- function(fitobject) {
  if (is.null(fitobject$envir$glm_phi)) { 
    glm_phi_args <- c(fitobject$phi.object$glm_phi_args, 
                      list(formula=attr(fitobject$resid.predictor,"oriFormula"),
                           lev=fitobject$lev_phi, data=fitobject$data,  
                           family= fitobject$resid.family)
    )
    fitobject$envir$glm_phi <-  do.call(".calc_dispGammaGLM", glm_phi_args)
  } 
  return(fitobject$envir$glm_phi)
}



.calc_wAugX <- function(augX,sqrt.ww) {
  wAugX <- .Dvec_times_m_Matrix(sqrt.ww, augX) # sweep(TT,MARGIN=1,sqrt.ww,`*`) # rWW %*%TT
  return(  structure(wAugX,envir=list2env(list(tag="wAugX",callcount=0L),parent=environment(HLfit_body))) )
}

.calc_TT <- function(AUGI0_ZX,ZAL) {
  if (inherits(ZAL,"Matrix")) {
    TT <- suppressMessages(cbind2(
      rbind2(AUGI0_ZX$X.pv,AUGI0_ZX$ZeroBlock), 
      rbind2(ZAL, AUGI0_ZX$I)
    )) 
  } else {
    TT <- cbind(
      rbind(AUGI0_ZX$X.pv,AUGI0_ZX$ZeroBlock), 
      rbind(ZAL, AUGI0_ZX$I)
    ) 
  }
  return(TT) ## aug design matrix 
}

.calc_beta_cov_spprec <- function(AUGI0_ZX,BLOB,w.resid) { ## actually full beta_v_cov is computed !
  if (is.null(chol_Md2hdv2 <- BLOB$chol_Md2hdv2)) {
    factor_Md2hdv2 <- as.matrix(Matrix::solve(BLOB$chol_Q,as(BLOB$G_CHMfactor,"sparseMatrix")))
    BLOB$chol_Md2hdv2 <- chol_Md2hdv2 <- chol(tcrossprod(factor_Md2hdv2)) #Md2hdv2 = U'.U 
  }
  r12 <- Matrix::solve(BLOB$G_CHMfactor, BLOB$ZtW %*% AUGI0_ZX$X.pv,system="L") ## solve(as(BLOB$G_CHMfactor,"sparseMatrix"), BLOB$ZtW %*% AUGI0_ZX$X.pv)
  r22 <- chol(.ZtWZwrapper(AUGI0_ZX$X.pv,w.resid)-crossprod(r12)) ## both lines as explained in working doc
  n_u_h <- ncol(chol_Md2hdv2)
  pforpv <- ncol(AUGI0_ZX$X.pv)
  Rmat_crossT <- rbind(cbind(chol_Md2hdv2,as.matrix(r12)),
                       cbind(matrix(0,nrow=pforpv,ncol=n_u_h),r22)) ## R_a in the working doc
  v_beta_cov <- chol2inv(Rmat_crossT) ## it happens that the beta,beta block is solve(crossprod(r22)), which is used .calc_inv_beta_cov() but not here.
  seqp <- seq_len(pforpv)
  perm <- c(n_u_h+seqp,seq_len(n_u_h))
  beta_v_cov <- v_beta_cov[perm,perm,drop=FALSE]
  beta_cov <- beta_v_cov[seqp,seqp,drop=FALSE]
  colnames(beta_cov) <- rownames(beta_cov) <- colnames(AUGI0_ZX$X.pv)
  attr(beta_cov,"beta_v_cov") <- beta_v_cov ## with a L that differs from other MME_methods hance beta_v_cov too (!) 
  # ./. but we can check that is is = solve(crossprod(wAugX)) by reconstructing wAugX with the matching L as in .get_LSmatrix()
  return(beta_cov)
}

.calc_beta_cov_others <- function(wAugX=NULL, 
                          AUGI0_ZX, 
                          ZAL, ww) {
  if (is.null(wAugX)) {
    if (is.null(ZAL)) {
      wAugX <- .calc_wAugX(augX=AUGI0_ZX$X.pv,sqrt.ww=sqrt(ww))
    } else {
      TT <- .calc_TT(AUGI0_ZX=AUGI0_ZX,ZAL)
      wAugX <- .calc_wAugX(augX=TT,sqrt.ww=sqrt(ww))
    }
  } ## wAugX is in XZ_OI order 
  if (inherits(wAugX,"Matrix")) {
    mMatrix_method <- .spaMM.data$options$Matrix_method 
  } else mMatrix_method <- .spaMM.data$options$matrix_method
  # hack to recycle sXaug code; 
  wAugX <- do.call(mMatrix_method,list(Xaug=wAugX, weight_X=rep(1,nrow(AUGI0_ZX$X.pv)), 
                                       w.ranef=rep(1,ncol(AUGI0_ZX$I)), ## we need at least its length for get_from Matrix methods
                                       H_global_scale=1))
  beta_v_cov <- get_from_MME(wAugX,"beta_v_cov_from_wAugX") ## = solve(crossprod(wAugX))
  pforpv <- ncol(AUGI0_ZX$X.pv)
  beta_cov <- beta_v_cov[seq_len(pforpv),seq_len(pforpv),drop=FALSE]
  colnames(beta_cov) <- rownames(beta_cov) <- colnames(AUGI0_ZX$X.pv)
  attr(beta_cov,"beta_v_cov") <- beta_v_cov
  return(beta_cov)
}

.get_beta_cov <- function(res) { ## this is messy: F I X M E
  # (1) res$beta_cov is provided by the fit, without a beta_v_cov attribute
  # (2) What this function provides is a full beta_v_cov, althouth as an attribute of another copy of beta_cov, res$envir$beta_cov.
  if (is.null(res$envir$beta_cov)) { 
    if ( "AUGI0_ZX_sparsePrecision" %in% res$MME_method) {
      ## use .calc_beta_cov_spprec special code NOT calling get_from_MME()
      nrd <- length(res$w.ranef)
      pforpv <- ncol(res$X.pv)
      AUGI0_ZX <- list(I=suppressWarnings(as(Diagonal(n=nrd),"CsparseMatrix")),
                       ZeroBlock=Matrix(0,nrow=nrd,ncol=pforpv),X.pv=res$X.pv)
      res$envir$beta_cov <- .calc_beta_cov_spprec(AUGI0_ZX=AUGI0_ZX, BLOB=res$envir, w.resid=res$w.resid)
    } else {
      ## uses .calc_beta_cov_others -> get_from_MME()
      ZAL <- get_ZALMatrix(res) ## should later simplify as =(res$QRmethod=="dense")) FIXME 
      nrd <- length(res$w.ranef)
      pforpv <- ncol(res$X.pv)
      if (inherits(ZAL,"Matrix")) {
        AUGI0_ZX <- list(I=suppressWarnings(as(Diagonal(n=nrd),"CsparseMatrix")),
                         ZeroBlock=Matrix(0,nrow=nrd,ncol=pforpv),X.pv=res$X.pv)
      } else {
        AUGI0_ZX <- list(I=diag(nrow=nrd),ZeroBlock=matrix(0,nrow=nrd,ncol=pforpv),X.pv=res$X.pv)
      }
      res$envir$beta_cov <- .calc_beta_cov_others(AUGI0_ZX=AUGI0_ZX,ZAL=ZAL,ww=c(res$w.resid,res$w.ranef))
    }
  } 
  return(res$envir$beta_cov)
}


# not doc'ed (no mention of augmented model in doc)
.get_LSmatrix <- function(object,augmented=FALSE) { 
  ## gets inv(tX_a invSig_a X_a).tX_a invSig_a that gives hat(beta,v_h)
  ww <- c(object$w.resid, object$w.ranef)
  sqrt.ww <- sqrt(ww)
  pforpv <- ncol(object$X.pv)
  nrd <- length(object$w.ranef)
  nobs <- nrow(object$X.pv)
  ZAL <- get_ZALMatrix(object)  
  augX <- cbind2(
    rbind2(object$X.pv, matrix(0,nrow=nrd,ncol=pforpv)), 
    rbind2(ZAL, diag(nrow=nrd))
  ) ## template with ZAL block to be filled later
  wAugX <- .calc_wAugX(augX=augX,sqrt.ww=sqrt.ww)
  beta_cov <- .get_beta_cov(object)
  beta_v_cov <- attr(beta_cov,"beta_v_cov")
  augXWXXW <- beta_v_cov %*% crossprod(wAugX, diag(x=sqrt.ww)) ## hmmm solve(crossprod(wAugX)) %*% crossprod(wAugX, diag(x=sqrt.ww))
  if (augmented) {
    return(augXWXXW)
  } else {
    return(augXWXXW[seq_len(pforpv),seq_len(nobs)])
  }
}


.calc_newZACvar <- function(newZAlist,cov_newLv_oldv_list) {
  newZACvarlist <- vector("list",length(newZAlist))
  for (new_rd in seq_along(newZAlist)) {
    terme <- newZAlist[[new_rd]] %id*id% (cov_newLv_oldv_list[[new_rd]])[]
    terme <- as.matrix(terme)
    newZACvarlist[[new_rd]] <- as.matrix(terme)
  }
  return(do.call(cbind,newZACvarlist))
}

## returns a list !!
## input XMatrix is either a single LMatrix which is assumed to be the spatial one, or a list of matrices 
.compute_ZAXlist <- function(XMatrix,ZAlist) {
  ## ZAL is nobs * (# levels ranef) and ZA too
  ## XMatrix is (# levels ranef) * (# levels ranef) [! or more generally a list of matrices!]
  ## the levels of the ranef must match each other in multiplied matrices
  ## the only way to check this is to have the levels as rownames and colnames and to check these
  if (is.null(ZAlist)) return(list())
  ## ELSE
  ZAX <- ZAlist
  if ( ! is.null(XMatrix) && length(ZAlist)>0 ) {
    if ( ! inherits(XMatrix,"list")) XMatrix <- list(dummyid=XMatrix)
    LMlen <- length(XMatrix)
    for (ii in seq_len(LMlen)) {
      xmatrix <- XMatrix[[ii]]
      ## find ZAlist elements affected by LMatrix element
      affecteds <- which(attr(ZAlist,"exp_ranef_strings") %in% attr(xmatrix,"ranefs"))
      for (it in affecteds) {
        ZA <- ZAlist[[it]]
        if (.is_identity(ZA)) {
          ZAX[[it]] <- xmatrix          
        } else {
          locnc <- ncol(ZA)
          locnr <- nrow(xmatrix)
          if ( locnc %% locnr !=0) {
            mess <- paste("The number of levels of the grouping variable in random term ", attr(ZAlist,"exp_ranef_strings")[it],sep="")
            mess <- paste(mess,"\n  is not a multiple of the dimension of the correlation matrix.") ## by distMatrix checking in corrHLfit or no.info check somewhere...
            stop(mess)
          }         
          nblocks <- locnc %/% locnr 
          if (nblocks>1) {
            locZA <- ZA
            for (bt in 1:nblocks) 
              locZA[,locnr*(bt-1)+(1:locnr)] <- locZA[,locnr*(bt-1)+(1:locnr)] %*% xmatrix[] ## [] to handle ff_matrix
            ZAX[[it]] <- locZA
          } else {
            ## With a proxy::dist or 'crossdist' matrix, it is likely that ZA was = I and we don't reach this code;
            ## However, exceptions can occur: cf Infusion with CIpoint = MLE => replicate in points where MSEs are to be estimated
            ## Then the xmatrix must have been converted from proxy style to a matrix.
            # In random slope models, xmatrix can be a Matrix
            zax <- ZA %*% xmatrix[]
            attr(zax,"corr.model") <- attr(xmatrix,"corr.model")
            if (identical(attr(zax,"corr.model"),"random-coef")) colnames(zax) <- colnames(ZA)
            #if (named) colnames(ZAX[[it]) <- colnames(ZA) ## names needed for match_old_new_levels()
            ZAX[[it]] <- zax
          }
        }
      }
    }
  }
  return(ZAX)
}

# even though the Z's were sparse postmultplication by LMatrix leads some of the ZAL's to dgeMatrix (dense)
.post_process_ZALlist <- function(ZALlist, as_matrix ) {
  nrand <- length(ZALlist)
  if ( as_matrix ) {
    ZALlist <- lapply(ZALlist,as.matrix) 
    ZAL <- do.call(cbind,ZALlist)
  } else {
    for (it in seq_len(length(ZALlist))) if (inherits(ZALlist[[it]],"dgeMatrix")) ZALlist[[it]] <- as(ZALlist[[it]],"dgCMatrix")
    ## but leave diagonal matrix types unchanged 
    ZAL <- suppressMessages(do.call(cbind,ZALlist)) ## suppress signature messages cbind2(Matrix,Matrix) ## (fixme: hides suboptimal code) 
  } 
  return(ZAL)
}

.compute_ZAL <- function(XMatrix,ZAlist, as_matrix) {
  ZALlist <- .compute_ZAXlist(XMatrix,ZAlist)
  ZAL <- .post_process_ZALlist(ZALlist, as_matrix )
  return(ZAL)
}



## cette fonction marche que si on a fixed effect + un terme aleatoire....
.eval_corrEst_args <- function(family,rand.families,predictor,data,X.Re,
                              REMLformula,ranFix,
                              term=NULL,
                              Optimizer) {
  ## ici on veut une procedure iterative sur les params de covariance
  #  HLCor.args$processed <- processed ## FR->FR dangerous in early development
  corrEst.args <- list(family=family,rand.family=rand.families) ## but rand.families must only involve a single spatial effect 
  loc.oriform <- attr(predictor,"oriFormula")
  loc.lhs <- paste(loc.oriform)[[2]]
  ## build formula, by default with only spatial effects
  if (is.null(term)) term <- .findSpatial(loc.oriform)
  corrEst.form <-  as.formula(paste(loc.lhs," ~ ",paste(term)))
  corrEst.args$data <- data ## FR->FR way to use preprocess ???                    
  # if standard ML: there is an REMLformula ~ 0; ____processed$X.Re is 0-col matrix, not NULL____
  # if standard REML: REMLformula is NULL: processed$X.Re is NULL
  # non standard REML: other REMLformula: processed$X.Re may take essentially any value
  if (is.null(X.Re) ) { ## processed$X.Re should be NULL => actual X.Re=X.pv => standard REML 
    corrEst.args$REMLformula <- predictor ## standard REML 
  } else corrEst.args$REMLformula <- REMLformula ## _ML_ _or_ non-standard REML
  if (NCOL(X.Re)>0) { ## some REML correction (ie not ML)
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

## Aggregate info on corrpars, inner-estimated and inner-fixed.
## $CorrPars is only for info in messages() and return value, 
.get_CorrEst_and_RanFix <- function(ranFix, ## has "fix", "outer", and also "var" values ! code corrently assumes "var" <=> corr_est
                             corr_est ## correlation parameters estimated within HLfit_body
                             ) {
  CorrEst_and_RanFix <- ranFix
  corrNames_in_corr_est <- names(corr_est) # rho,nu,  pas trRho, trNu 
  CorrEst_and_RanFix[corrNames_in_corr_est] <- corr_est
  if (is.null(alltypelist <- attr(ranFix,"type"))) {
    alltypelist <- list()
    alltypelist[names(ranFix)] <- "fix"
  }
  alltypelist[corrNames_in_corr_est] <- "var" ## eg rho from HLCor(adjacency) ## maybe not useful as alredy set to "var"
  attr(CorrEst_and_RanFix,"type") <- alltypelist
  #
  corrNames_in_both <- intersect(names(CorrEst_and_RanFix),c("nu","rho","Nugget","ARphi"))
  corrPars <- CorrEst_and_RanFix[corrNames_in_both] ## as for functions in .corrMM_LRT that always look in phi, lambda, rather than .Fix  
  attr(corrPars,"type") <- attr(CorrEst_and_RanFix,"type")[corrNames_in_both]
  return(list(corrPars=corrPars, ## correlation parameters with the appropriate types "fix" or "var"
              CorrEst_and_RanFix=CorrEst_and_RanFix) ## correlation params + whatever else was in ranFix
         ) 
}

.make_lambda_object <- function(nrand, lambda_models, cum_n_u_h, lambda_est, init.lambda,                        
                               process_resglm_blob, ZAlist, next_LMatrices) {
  ## redefine names(namesTerms) and namesTerms elements
  print_namesTerms <- attr(ZAlist,"namesTerms") ## a list, which names correspond to the grouping variable, and elements are the names of the coefficients fitted
  namesnames <- names(print_namesTerms)
  for (it in seq_len(length(namesnames))) if (nchar(namesnames[it])>10) namesnames[it] <- paste(substr(namesnames[it],0,9),".",sep="")
  names(print_namesTerms) <- make.unique(namesnames,sep=".") ## makes group identifiers unique (names of coeffs are unchanged); using base::make.unique
  print_lambda <- vector("list",nrand) 
  for (it in seq_len(nrand)) {
    plam <- process_resglm_blob$print_lambdas[[it]]
    if (anyNA(plam)) { ## those for which no glm was available, such as fixed lambdas...
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      plam <- unique(lambda_est[u.range])
    }
    print_lambda[[it]] <- structure(plam, names=names(print_namesTerms)[it])
  }
  attr(print_lambda,"cum_n_u_h") <- cum_n_u_h
  lambda.object <- list(lambda_est = lambda_est,  ## full vector for simulate() calc_logdisp_cov()
                        lambda=print_lambda)  ## nrand-elements list, in output -> used by simulate, useful for init another fit, may substitute to the coefficients_lambdaS when the latter have not bee computed from  glm, etc.
  lambda.object$type <- type <- attr(init.lambda,"type")
  if (any(type=="inner")) { ## modifies default namesTerms
    coefficients_lambdaS <- process_resglm_blob$coefficients_lambdaS
    for (it in seq_len(length(coefficients_lambdaS))) { ## detect exceptions
      if (type[it]=="inner") {
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

.make_strucList <- function(ZAlist, next_LMatrices, coefficients_lambdaS, ranCoefs_blob) {
  strucList <- vector("list", length(ZAlist)) 
  for (it in seq_len(length(ZAlist))) {
    lmatrix <- next_LMatrices[[it]]
    if ( ! is.null(lmatrix)) {
      corr.model <- attr(lmatrix, "corr.model") 
      if (is.null(corr.model)) warning('attr(next_LMatrix, "corr.model") is NULL')
      if (corr.model=="random-coef") {
        ## keep in mind that str(S4...) does not show extra attributes
        ## lmatrix is an S4 object...
        ## ranefs itself has attribute type="(.|.)"
        if (ranCoefs_blob$is_set[it]) {
          attr(lmatrix,"cov.mat") <- attr(lmatrix,"latentL_blob")$compactcovmat 
        } else {
          attr(lmatrix,"cov.mat") <- .ZWZt(attr(lmatrix,"latentL_blob")$u,
                                           exp(coefficients_lambdaS[[it]])) 
        }
        longLmatrix <- .makelong(attr(lmatrix,"latentL_blob")$u,longsize=ncol(ZAlist[[it]]))
        lmatrix <- do.call(structure,c(list(.Data=longLmatrix),
                                       attributes(lmatrix)[c("latentL_blob","par","ranefs","corr.model","cov.mat")]))
      } 
      strucList[[it]] <- lmatrix  
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
    if (models[["phi"]]=="phiScal") {
      H_global_scale <- phi_est[1L]
    } else H_global_scale <- exp(mean(log(phi_est))) ## fixem a more trivial scaling (no scaling?) would be sufficient ?
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

.damping_to_solve <- function(X, dampDpD, rhs=NULL,method="QR", .drop=TRUE) { ##  
  if (.drop) rhs <- drop(rhs) ## 1-col m/Matrix to vector ## affects indexing below but the result of the chol2inv line is always Matrix
  if (method=="QR") {
    if (inherits(X,"Matrix")) {
      XD <- rbind2(X, Diagonal(x=sqrt(dampDpD))) ## fixme might even be possible to store an rbind result...
      RP <- qr(XD) ## i.e. Matrix::qr
      RRsP <- sort.list(RP@q) 
      # solve(crossprod(XD)) = chol2inv(Matrix::qrR(RP,backPermute = FALSE))[RRsp,RRsP]
      # solve(crossprod(XD), rhs) = (Matrix::chol2inv(Matrix::qrR(RP, backPermute = FALSE)) %*% rhs[sort.list(RRsP)])[RRsP]
      if (is.null(rhs)) {
        return(list(inv=Matrix::chol2inv(Matrix::qrR(RP,backPermute = FALSE)),Rperm=RP@q+1L,RRsP=RRsP))
      } else {
        if (is.matrix(rhs)) {
          resu <- (Matrix::chol2inv(Matrix::qrR(RP,backPermute = FALSE)) %*% rhs[RP@q+1L,])[RRsP,]
        } else { ## if rhs is one col 
          resu <- (Matrix::chol2inv(Matrix::qrR(RP,backPermute = FALSE)) %*% rhs[RP@q+1L])[RRsP]
        }
      }
      ## assuming rhs is 'slim', this minimizes permutations. For computation of dv_h, it is square rather than slim...
    } else {
      XD <- rbind(X, diag(x=sqrt(dampDpD),ncol=length(dampDpD)))
      RP <- .lmwithQRP(XD,yy=NULL,returntQ=FALSE,returnR=TRUE)
      RRsP <- sort.list(RP$perm)
      # dVscaled <- chol2inv(RP$R_scaled)[RRsP,RRsP] %*% rhs
      if (is.null(rhs)) {
        return(list(inv=chol2inv(RP$R_scaled),Rperm=RP$perm+1L,RRsP=RRsP))
      } else {
        if (is.matrix(rhs)) {
          resu <- (chol2inv(RP$R_scaled) %*% rhs[RP$perm+1L,])[RRsP,]  ## assuming rhs is 'slim', this minimizes permutations        
        } else { ## if rhs is one col 
          resu <- (chol2inv(RP$R_scaled) %*% rhs[RP$perm+1L])[RRsP]
        }
      }
    }
    if (.drop) resu <- drop(resu)
    return(resu)
  } else {
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
