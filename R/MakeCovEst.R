# HLfit(Reaction ~ 0 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ 1 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy,HLmethod="ML")

.calc_latentL <- function(compactcovmat,dtC=TRUE) {
  ## the Cholesky factorization is not unique, so there's no reason for the factorization of compactcovmat and compactprecmat to match
  #  unless one is deduced from the other. L and L_Q are distinct matrices used jointly in a fit (L in ZAL...).
  #  Hence it's pathetic to recompute a Chol from the cov somewhere and another from the prec elsewhere.
  #  In sparse precision code We want chol_Q to be lower triangular (dtCMatrix can be Up or Lo) 
  #    to obtain a dtCMatrix by bdiag(list of lower tri dtCMatrices)  
  if (dtC) {
    # compactprecmat <- solve(compactcovmat)
    blob <- sym_eigen(compactcovmat) ## sym_eigen better handles nearly singular matrices but then the sparse precision code is not numerically robust;
    blob$d <- pmax(blob$d,1e-8) ## FIXME ad hoc value
    solveu <- solve(blob$u)
    compactprecmat <- t(solveu) %*% diag(1/blob$d) %*% solveu
    ## lower tri dtC chol_Q:
    chol_prec <- as(Cholesky(as(compactprecmat,"sparseMatrix"),LDL=FALSE,perm=FALSE), "sparseMatrix")
    dcp <- diag(chol_prec)
    chol_Q <- chol_prec %*% Diagonal(x=1/dcp)
    design_u <- t(solve(as.matrix(chol_Q))) ## upper tri (while sym_eigen provides a "full" U matrix)
    ## design_u -> as.matrix bc of the code .ZWZt(attr(lmatrix, "latentL_blob")$u, exp(coefficients_lambdaS[[it]]))
    blob <- list(u=design_u, # design matrix for iterative algo
                 d=1/dcp^2, # given variances of latent independent ranefs in outer optim algo
                 compactchol_Q=chol_Q, # precision factor for sparse_precision algo
                 compactprecmat=compactprecmat, compactcovmat=compactcovmat) 
    # COV= blob$u %*% diag(blob$d) %*% t(blob$u); prec=  blob$compactchol_Q %*% diag(1/blob$d) %*% t(blob$compactchol_Q)
  } else { ## older code sufficient for iterative algo
    blob <- sym_eigen(compactcovmat) ## COV= blob$u %*% diag(blob$d) %*% t(blob$u) ## 
    ## the following reordering affects the order of <fit>$lambda elements hence the summary.HLfit(.,details=TRUE). It should nto affect the quality of the fit. 
    blib <- Matrix::expand(Matrix::lu(t(blob$u))) ## use pivoting in lu as a useful permutation... (t(blob$u)= PLU)
    blob <- list(u=t(as.matrix(with(blib,L %*% U))),
                 d=blob$d[blib$P@perm], compactcovmat=compactcovmat) ## COV= blob$u %*% diag(blob$d) %*% t(blob$u) again
  }
  return(blob)
  ## for sparse precision we want chol_Q to be (dtCMatrix: Csparse triangular) so that efficient solve methods can be used.
  ## This is not provided by this function
}

## Build Lmatrix between all pairs of u_h (nr*nr) from parameter estimates (2*2) for the design matrix, longsize is final dim of matrix
.makelong <- function(Lcompact,longsize,as_matrix=FALSE) { ## always returns a Matrix unless explicitly as_matrix
  blocksize <- longsize/ncol(Lcompact)
  if (longsize>180L) {
    longLv <- Diagonal(n=longsize) ## declaration
  } else longLv <- diag(nrow=longsize) ## declaration
  for (it in seq_len(ncol(Lcompact))) {
    urange1 <- (it-1)*blocksize + seq(blocksize)
    longLv[cbind(urange1,urange1)] <- Lcompact[it,it]
    for (jt in seq_len(it-1)) {
      urange2 <- (jt-1)*blocksize + seq(blocksize)
      longLv[cbind(urange1,urange2)] <- Lcompact[it,jt]
      longLv[cbind(urange2,urange1)] <- Lcompact[jt,it]
    }
  }
  if (as_matrix) {
    longLv <- as.matrix(longLv)
  } else {
    ## returns a *M*atrix in all cases: 
    if (inherits(Lcompact,"dtCMatrix")) {
      longLv <- as(longLv,"dtCMatrix")
    } else if ( ! inherits(longLv,"Matrix")) {
      longLv <- Matrix(longLv)
    }
  }
  return(longLv) 
} ## end def makelong

.ad_hoc_calc_inits_ranCoef <- function(Xi_ncol,tol_ranCoefs=1e-08,init) { ## FIXME not yer user control
  lowerb <- upperb <- matrix(NA,nrow=Xi_ncol,ncol=Xi_ncol)
  lowerbloc <- lower.tri(lowerb,diag=TRUE) ## a matrix of T/F !
  diag(lowerb) <- log((tol_ranCoefs))/2
  diag(upperb) <- log(1/(tol_ranCoefs))/2
  lowerb[lower.tri(lowerb)] <-   -(1-tol_ranCoefs)
  upperb[lower.tri(upperb)] <-   (1-tol_ranCoefs)
  if (is.null(init)) {
    init <- (upperb+lowerb)/2
    diag(init) <- 0
    init <- init[lowerbloc]        
  }
  upperb <- upperb[lowerbloc]
  lowerb <- lowerb[lowerbloc]
  init <- pmin(upperb,init) ## possible 8th-decimal issues...
  init <- pmax(lowerb,init)
  return(list(lowerb=lowerb,upperb=upperb,init=init))
}

.makeLcovLt <- function(parvec, Xi_ncol=NULL) {
  if (is.null(Xi_ncol)) Xi_ncol <- floor(sqrt(2*length(parvec)))
  compactLv <- matrix(0,nrow=Xi_ncol,ncol=Xi_ncol)
  lowerbloc <- lower.tri(compactLv,diag=TRUE) ## a matrix of T/F !
  compactLv[lowerbloc] <- parvec
  compactLv[t(lowerbloc)] <- parvec
  sigmas <- diag(exp(diag(compactLv))) 
  diag(compactLv) <- 1
  resu <- sigmas %*% compactLv %*%sigmas
  resu
}


.makeCovEst1 <- function(u_h,ZAlist,cum_n_u_h,prev_LMatrices,
                         var_ranCoefs,w.resid,processed,phi_est,lcrandfamfam,family,
                        as_matrix,v_h, MakeCovEst_pars_not_ZAL_or_lambda
) {
  nrand <- length(ZAlist)
  locX.Re <- processed$X.Re ## may be NULL 
  if (is.null(locX.Re)) locX.Re <- processed$AUGI0_ZX$X.pv
  locpredictor <- processed$predictor
  updated_LMatrices <- working_LMatrices <- prev_LMatrices
  #if (is.null(next_LMatrices)) next_LMatrices <- list() ## NULL wrong for next_LMatrices[[rt]] <- <*M*atrix>
  Xi_cols <- attr(ZAlist,"Xi_cols")
  loc_lambda_est <- numeric(length(u_h))
  for (rt in seq_len(nrand)) {
    if ( var_ranCoefs[rt]) { ## inner estimation of cov mat of u_h 
      Xi_ncol <- Xi_cols[rt]
      blocksize <- ncol(ZAlist[[rt]])/Xi_ncol 
      COVpredUnderHuncorr <- matrix(0,ncol=Xi_ncol,nrow=Xi_ncol) ## var on diag, corr outside diag
      ##prevL <- attr(prev_LMatrices[[rt]],"latentL_blob")$u
      u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
      ########## brute force optimization
      objfn <- function(parvec) {
        compactcovmat <- .makeLcovLt(parvec, Xi_ncol)
        ## cosmetic / interpretative permutation
        ## assignments as design matrix and lambda values:
        blob <- .calc_latentL(compactcovmat,dtC=FALSE)
        loc_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) 
        loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 ## arbitrarily small eigenvalue is possible for corr=+/-1 even for 'large' parvec
        ## we have a repres in terms of ZAL and of a diag matrix of variances; only the latter affects hlik computation
        longLv <- .makelong(blob$u,longsize=ncol(ZAlist[[rt]])) ## the variances are taken out in $d
        working_LMatrices[[rt]] <- longLv
        attr(working_LMatrices[[rt]],"ranefs") <- attr(ZAlist,"exp_ranef_strings")[[rt]] ## FR->FR  revoir pour matrices affectant +s termes ?
        locZAL <- .compute_ZAL(XMatrix=working_LMatrices,ZAlist=ZAlist,as_matrix=as_matrix) 
        # processed$AUGI0_ZX$Xaug <- NULL ## if .make_Xscal( . , use_Xaug=TRUE)
        locw.ranefSblob <- .updateW_ranefS(cum_n_u_h,processed$rand.families,lambda=loc_lambda_est,u_h,v_h) 
        locarglist <- c(MakeCovEst_pars_not_ZAL_or_lambda, list(ZAL=locZAL, lambda_est=loc_lambda_est, wranefblob=locw.ranefSblob))
        #locarglist$eta <- NULL
        #locarglist$muetablob <- NULL
        #locarglist$w.resid <- NULL
        auglinmodblob <- do.call(".solve_IRLS_as_ZX",locarglist)
        if ( ncol(locX.Re) ) { lik_obj="p_bv"} else lik_obj="p_v" ## FIXME could be taken out of locfn ?
        ## FR->FR but its not clear that non-standard REML is handled by calc_APHLs_from_ZX !!
        REMLcrit <- .calc_APHLs_from_ZX(auglinmodblob,which=lik_obj,processed=processed)[[lik_obj]]
        return(REMLcrit)
      } ## currently this refits the fixed effects together with the other params... probably not optimal...
      ####  
      init_blob <- .ad_hoc_calc_inits_ranCoef(Xi_ncol,init=attr(prev_LMatrices[[rt]],"par"))
      lowerb <- init_blob$lowerb
      upperb <- init_blob$upperb
      init <- init_blob$init
      parscale <- (upperb-lowerb)
      if (TRUE) {
        objfn_nloptr <- function(x) { return( - objfn(x)) }
        nloptr_controls <- spaMM.getOption("nloptr") 
        optr <- nloptr::nloptr(x0=init,eval_f=objfn_nloptr,lb=lowerb,ub=upperb,
                       opts=nloptr_controls)
        while (optr$status==5L) { ## 5 => termination bc maxeval has been reached
          prevlik <- optr$objective
          reinit <- pmax(lowerb,pmin(upperb,optr$solution))
          optr <- nloptr::nloptr(x0=reinit,eval_f=objfn_nloptr,lb=lowerb,ub=upperb,
                                 opts=nloptr_controls)
          loc_ftol <- max(1e-8, optr$options$ftol_abs)
          if (- optr$objective < - prevlik+loc_ftol) break ## no progress in <= maxeval iterations
        }
        optr$par <- optr$solution
      } else { ## deprecated:
        ################# OPTIM
        optr <- optim(init,objfn,lower=lowerb,upper=upperb,method="L-BFGS-B", ## deprecated
                      control=list(parscale=parscale,fnscale=-1))
      }
      ################# 
      ## reproduces representation in objfn
      COVcorr <- .makeLcovLt(optr$par, Xi_ncol)
      blob <- .calc_latentL(COVcorr,dtC=FALSE)
      loc_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) 
      loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 
      next_LMatrix <- .makelong(blob$u,longsize=ncol(ZAlist[[rt]])) ## il faut updater pour estimer les ranef correctement...
      attr(next_LMatrix,"latentL_blob") <- blob ## kept for updating in next iteration and for output
      attr(next_LMatrix,"par") <- optr$par ## kept for updating in next iteration and for output
      attr(next_LMatrix,"ranefs") <-  structure(attr(ZAlist,"exp_ranef_strings")[rt], 
                                                ## type is "(.|.)" if LMatrix is for random slope ## ajout 2015/06
                                                type= attr(ZAlist,"exp_ranef_types")[rt] )
      attr(next_LMatrix, "corr.model") <- "random-coef"
      updated_LMatrices[[rt]] <- next_LMatrix # (fixme ??) working_LMatrices[[rt]] <- 
    } else updated_LMatrices[rt] <- list(NULL) ## this erases Matern, AR1, and other fixed LMatrices, so updated_LMatrices is not to be used in objfn()   
  } ## loop on rt = ranefs
  return(list(updated_LMatrices=updated_LMatrices, 
              next_lambda_est=loc_lambda_est,
              optr_par=optr$par))
} ## end def makeCovEst1

## parait lent à converger vers -1; c'est là que tester d'abord les bornes 0 et 1 pour choisir init peut être bien... 
## experimental version of MakeCovEst stored in MakeCovEst.R.txt
.makeCovEst2 <- function(u_h,ZAlist,cum_n_u_h,prev_LMatrices,
                       var_ranCoefs,w.resid,processed,prevZAL,clik) {
  nrand <- length(ZAlist)
  locX.Re <- processed$X.Re ## may be NULL
  if (is.null(locX.Re)) locX.Re <- processed$AUGI0_ZX$X.pv
  based2hdv2 <- structure( - crossprod(prevZAL, diag(w.resid)) %*% prevZAL, 
                           envir=list2env(list(tag="d2hdv2",callcount=0L),parent=environment(HLfit_body)))
  next_LMatrices <- prev_LMatrices
  Xi_cols <- attr(ZAlist,"Xi_cols")
  Lu <- u_h
  loc_lambda_est <- numeric(length(u_h))
  for (rt in seq_len(length(ZAlist))) {
    if (var_ranCoefs[rt]) {
      Xi_ncol <- Xi_cols[rt]
      blocksize <- ncol(ZAlist[[rt]])/Xi_ncol 
      COVpredUnderHuncorr <- matrix(0,ncol=Xi_ncol,nrow=Xi_ncol) ## var on diag, corr outside diag
      compactLv <- matrix(0,nrow=Xi_ncol,ncol=Xi_ncol)
      lowerbloc <- lower.tri(compactLv,diag=TRUE) ## a matrix of T/F !
      u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
      ########## brute force optimization
      makeLcovLt <- function(parvec) {
        compactLv[lowerbloc] <- parvec
        compactLv[t(lowerbloc)] <- parvec
        sigmas <- diag(exp(diag(compactLv))) 
        diag(compactLv) <- 1
        resu <- sigmas %*% compactLv %*%sigmas
        resu
      }
      ####
      # pour une repres non diagonal je devrais reconstruire une Lmatrix... (je l'ai  !)
      if (ncol(locX.Re)) hessnondiag <- crossprod(prevZAL,sweep(locX.Re,MARGIN=1,w.resid,`*`))  
      objfn <- function(parvec) {
        compactcovmat <- makeLcovLt(parvec)
        blob <- .calc_latentL(compactcovmat,dtC=FALSE) 
        prevL <- attr(prev_LMatrices[[rt]],"latentL_blob")$u
        corrSig <- compactcovmat
        if ( ! is.null(prevL)) corrSig <- crossprod(prevL, corrSig) %*% prevL ## corrects by what is absorbed in constant ZAL
        longcorrSig <- .makelong(corrSig,longsize=ncol(ZAlist[[rt]])) ## car t(prevL) =solve(prevL)
        solvecompact <- blob$u %*% diag(1/blob$d) %*% t(blob$u) ## plus stable
        corrinvSig <- solvecompact
        if ( ! is.null(prevL)) corrinvSig <- crossprod(prevL, corrinvSig) %*% prevL
        longcorrinvSig <- .makelong(corrinvSig,longsize=ncol(ZAlist[[rt]])) ## car t(prevL) =solve(prevL)
        # equiv entre ligne suivante et une repres diagonale car identique a pre/ post par matrice unitaire (! important que unitaire !!!! )
        locd2hdv2 <- structure( based2hdv2 - longcorrinvSig, 
                                envir=list2env(list(tag="d2hdv2",callcount=0L),parent=environment(HLfit_body)))  
        # pour une repres non diagonal je devrais reconstruire une Lmatrix... (je l'ai  !)
        if (ncol(locX.Re)==0L) { ## fit ML: p_bv=p_v hence d2hdpbv reduces to d2hdv2
          lad <- - .LogAbsDetWrap( - locd2hdv2,logfac=-log(2*pi))
        } else { 
          Md2hdbv2 <- rbind(cbind(.ZtWZwrapper(locX.Re,w.resid), t(hessnondiag)),
                            cbind(hessnondiag, - locd2hdv2)) 
          lad <- - .LogAbsDetWrap(Md2hdbv2,logfac=-log(2*pi))
        } ## le lad est OK = ladbv de la version standard
        likranu <- ( blocksize* .LogAbsDetWrap(corrSig,logfac=log(2*pi)) + u_h %*% longcorrinvSig %*% u_h)/2
        # pas passer par logdet prevL car ?? LogAbsDetWrap(prevL non sym) incorrect ?? 
        obj <- clik - likranu + lad/2 ## different from REMLcrit (normal) but still appears to be maximized
        ## ca marche presque mais l'approx du gradient est très mauvaise pour les corr extremes...
        ######################################
        return(obj)
        ## a comparer a:
        #return(REMLcrit)
      } 
      
      ####  
      lowerb <- upperb <- matrix(NA,nrow=Xi_ncol,ncol=Xi_ncol)
      diag(lowerb) <- log(sqrt(1e-08))
      diag(upperb) <- log(sqrt(1e08))
      init <- attr(prev_LMatrices[[rt]],"par")
      if (is.null(init)) {
        lowerb[lower.tri(lowerb)] <-   -0.99
        upperb[lower.tri(upperb)] <-   0.99
        init <- (upperb+lowerb)/2
        diag(init) <- 0
        init <- init[lowerbloc]        
      } else {
        lowerb[lower.tri(lowerb)] <-  (-99+init[2])/100 
        upperb[lower.tri(upperb)] <-   (99+init[2])/100          
      }
      upperb <- upperb[lowerbloc]
      lowerb <- lowerb[lowerbloc]
      parscale <- (upperb-lowerb)        
      ################# OPTIM
      optr <- optim(init,objfn,lower=lowerb,upper=upperb,method="L-BFGS-B", 
                    control=list(parscale=parscale,fnscale=-1))
      # print(optr$par)
      ################# 
      ## standardized representation
      COVcorr <- makeLcovLt(optr$par)
      blob <- .calc_latentL(COVcorr,dtC=FALSE)
      loc_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) 
      loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 
      next_LMatrix <- .makelong(blob$u,longsize=ncol(ZAlist[[rt]])) ## il faut updater pour estimer les ranef correctement...
      attr(next_LMatrix,"latentL_blob") <- blob ## kept for updating in next iteration and for output
      attr(next_LMatrix,"par") <- optr$par ## kept for updating in next iteration and for output
      attr(next_LMatrix,"ranefs") <-  structure(attr(ZAlist,"exp_ranef_strings")[rt], 
                                                ## type is "(.|.)" if LMatrix is for random slope ## ajout 2015/06
                                                type= attr(ZAlist,"exp_ranef_types")[rt] )
    } else next_LMatrix <- NULL
    next_LMatrices[[rt]] <- next_LMatrix
  } ## loop on rt = ranefs
  return(list(next_LMatrices=next_LMatrices,next_lambda_est=loc_lambda_est,
              optr_par=optr$par))
} ## end def MakeCovEst2
