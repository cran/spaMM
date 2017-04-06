# HLfit(Reaction ~ 0 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ 1 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy,HLmethod="ML")

# for fixed u_h, numerically maximize p_bv (p_v) wrt correlation params; atroce car pour chaque va de param -> objfn -> auglinmodfit 
# par contre aucune tentative de corriger les corr mat. la prevL n'est pas utilisée, elle impacte seulement u_h en input
.calc_latentL <- function(compactcovmat) {
  blob <- sym_eigen(compactcovmat) ## COVcorr= blob$u %*% diag(blob$d) %*% t(blob$u) ## 
  blib <- Matrix::expand(Matrix::lu(t(blob$u))) ## use pivoting in lu as a useful permutation...
  blob <- list(u=t(as.matrix(with(blib,L %*% U))),d=blob$d[blib$P@perm])
  return(blob)
}

## Build Lmatrix between all pairs of u_h (nr*nr) from parameter estimates (2*2) for the design matrix, longsize is final dim of matrix
.makelong <- function(Lcompact,longsize) {
  blocksize <- longsize/ncol(Lcompact)
  longLv <- diag(longsize) ## declaration
  for (it in seq_len(ncol(Lcompact))) {
    urange1 <- (it-1)*blocksize + seq(blocksize)
    longLv[cbind(urange1,urange1)] <- Lcompact[it,it]
    for (jt in seq_len(it-1)) {
      urange2 <- (jt-1)*blocksize + seq(blocksize)
      longLv[cbind(urange1,urange2)] <- Lcompact[it,jt]
      longLv[cbind(urange2,urange1)] <- Lcompact[jt,it]
    }
  }
  Matrix(longLv) ## Matrix: 02/2016
} ## end def makelong

makeCovEst1 <- function(u_h,ZAlist,cum_n_u_h,prev_LMatrices,
                        userLfixeds,w.resid,processed,phi_est,lcrandfamfam,family,
                        as_matrix,v_h,auglinfixedpars
) {
  nrand <- length(ZAlist)
  locX.Re <- processed$X.Re ## may be NULL 
  if (is.null(locX.Re)) locX.Re <- processed$X.pv
  locpredictor <- processed$predictor
  next_LMatrices <- prev_LMatrices
  if (is.null(next_LMatrices)) next_LMatrices <- list() ## NULL wrong for next_LMatrices[[rt]] <- <*M*atrix>
  Xi_cols <- attr(ZAlist,"Xi_cols")
  Lu <- u_h
  loc_lambda_est <- numeric(length(u_h))
  for (rt in seq_len(nrand)) {
    ## estimate correlation matrix 
    Xi_ncol <- Xi_cols[rt]
    ## cov mat of u_h if not fixed by user ## standard REML method 
    if ( Xi_ncol>1 && ! userLfixeds[rt]) {
      blocksize <- ncol(ZAlist[[rt]])/Xi_ncol 
      COVpredUnderHuncorr <- matrix(0,ncol=Xi_ncol,nrow=Xi_ncol) ## var on diag, corr outside diag
      ##prevL <- attr(prev_LMatrices[[rt]],"Lcompact")
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
      objfn <- function(parvec) {
        compactcovmat <- makeLcovLt(parvec)
        ## cosmetic / interpretative permutation
        ## assignments as design matrix and lambda values:
        blob <- .calc_latentL(compactcovmat)
        loc_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) 
        loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 ## arbitrarily small eigenvalue is possible for corr=+/-1 even for 'large' parvec
        ## we have a repres in terms of ZAL and of a diag matrix of variances; only the latter affects hlik computation
        longLv <- .makelong(blob$u,longsize=ncol(ZAlist[[rt]])) ## the variances are taken out in $d
        next_LMatrices[[rt]] <- longLv
          attr(next_LMatrices[[rt]],"ranefs") <- attr(ZAlist,"ranefs")[[rt]] ## FR->FR  revoir pour matrices affectant +s termes ?
        ZALlist <- .compute_ZAXlist(XMatrix=next_LMatrices,ZAlist=ZAlist)
        locZAL <- .post_process_ZALlist(ZALlist,as_matrix=as_matrix) 
        if (inherits(locZAL,"Matrix")) {
          ddi_or_matrix_locZAL <- as.matrix(locZAL)
        } else ddi_or_matrix_locZAL <- locZAL
        locw.ranefSblob <- .updateW_ranefS(cum_n_u_h,processed$rand.families,lambda=loc_lambda_est,u_h,v_h) 
        ## FR->FR auglinfixedpars is an ambiguous name since this contains $w.resid which is updated within auglinmodfit
        locarglist <- c(auglinfixedpars, list(ZAL=locZAL, lambda_est=loc_lambda_est, wranefblob=locw.ranefSblob))
        auglinmodblob <- do.call("fit_as_ZX",locarglist)
        if ( ncol(locX.Re)>0L ) { lik_obj="p_bv"} else lik_obj="p_v" ## could be taken out of locfn ?
        ## FR->FR but its not clear that non-standard REML is handled by calc_APHLs_from_ZX !!
        REMLcrit <- calc_APHLs_from_ZX(auglinmodblob,which=lik_obj,processed=processed)[[lik_obj]]
        return(REMLcrit)
      } ## currently this refits the fixed effects together with the other params... probably not optimal
      ####  
      lowerb <- upperb <- matrix(NA,nrow=Xi_ncol,ncol=Xi_ncol)
      diag(lowerb) <- log(sqrt(1e-08))
      diag(upperb) <- log(sqrt(1e08))
      lowerb[lower.tri(lowerb)] <-   -(1-1e-08)
      upperb[lower.tri(upperb)] <-   (1-1e-08)
      init <- attr(prev_LMatrices[[rt]],"par")
      if (is.null(init)) {
        init <- (upperb+lowerb)/2
        diag(init) <- 0
        init <- init[lowerbloc]        
      }
      upperb <- upperb[lowerbloc]
      lowerb <- lowerb[lowerbloc]
      init <- pmin(upperb,init) ## possible 8th-decimal issues...
      init <- pmax(lowerb,init)
      parscale <- (upperb-lowerb)
      if (TRUE) {
        objfn_nloptr <- function(x) { return( - objfn(x)) }
        nloptr_controls <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1.0e-4,maxeval=-1,print_level=0) ## DEFAULT
        optr <- nloptr(x0=init,eval_f=objfn_nloptr,lb=lowerb,ub=upperb,
                       opts=nloptr_controls)
        optr$par <- optr$solution
      } else {
        ################# OPTIM
        optr <- optim(init,objfn,lower=lowerb,upper=upperb,method="L-BFGS-B",
                      control=list(parscale=parscale,fnscale=-1))
      }
      ################# 
      ## reproduces representation in objfn
      COVcorr <- makeLcovLt(optr$par)
      blob <- .calc_latentL(COVcorr)
      loc_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) 
      loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 
      next_LMatrix <- .makelong(blob$u,longsize=ncol(ZAlist[[rt]])) ## il faut updater pour estimer les ranef correctement...
      attr(next_LMatrix,"Lcompact") <- blob$u ## kept for updating in next iteration and for output
      attr(next_LMatrix,"par") <- optr$par ## kept for updating in next iteration and for output
      thisranef <- attr(ZAlist,"ranefs")[rt]
      attr(thisranef,"type") <- attr(attr(ZAlist,"ranefs"),"type")[rt] ## ie simply "(.|.)": LMatrix with such type is for random slope ## ajout 2015/06
      attr(next_LMatrix,"ranefs") <- thisranef
      attr(next_LMatrix, "corr.model") <- "random-coef"
      next_LMatrices[[rt]] <- next_LMatrix
    } else next_LMatrices[rt] <- list(NULL)
  } ## loop on rt = ranefs
  return(list(next_LMatrices=next_LMatrices, 
              ## : a list of nrand elements, each NULL except for random slope terms
              next_lambda_est=loc_lambda_est,
              optr_par=optr$par))
} ## end def makeCovEst1

## parait lent à converger vers -1; c'est là que tester d'abord les bornes 0 et 1 pour choisir init peut être bien... 
## experimental version of MakeCovEst stored in MakeCovEst.R.txt
makeCovEst2 <- function(u_h,ZAlist,cum_n_u_h,prev_LMatrices,
                       userLfixeds,w.resid,processed,prevZAL,clik) {
  nrand <- length(ZAlist)
  locX.Re <- processed$X.Re ## may be NULL
  if (is.null(locX.Re)) locX.Re <- processed$X.pv
  based2hdv2 <- structure( - crossprod(prevZAL, diag(w.resid)) %*% prevZAL, 
                           envir=list2env(list(tag="d2hdv2",callcount=0L),parent=environment(HLfit_body)))
  next_LMatrices <- prev_LMatrices
  Xi_cols <- attr(ZAlist,"Xi_cols")
  Lu <- u_h
  loc_lambda_est <- numeric(length(u_h))
  for (rt in seq_len(length(ZAlist))) {
    ## estimate correlation matrix 
    Xi_ncol <- Xi_cols[rt]
    ## cov mat of u_h if not fixed by user ## standard REML method 
    if ( Xi_ncol>1 && ! userLfixeds[rt]) {
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
      if (ncol(locX.Re)>0L) hessnondiag <- crossprod(prevZAL,sweep(locX.Re,MARGIN=1,w.resid,`*`))  
      objfn <- function(parvec) {
        compactcovmat <- makeLcovLt(parvec)
        blob <- sym_eigen(compactcovmat) ## COVcorr= blob$u %*% diag(blob$d) %*% t(blob$u) ## u is unitary !!
        ## cosmetic / interpretative permutation
        blib <- Matrix::expand(Matrix::lu(t(blob$u))) ## use pivoting in lu as a useful permutation...
        blob <- list(u=t(as.matrix(with(blib,L %*% U))),d=blob$d[blib$P@perm]) 
        prevL <- attr(prev_LMatrices[[rt]],"Lcompact")
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
          lad <- - LogAbsDetWrap( - locd2hdv2,logfac=-log(2*pi))
        } else { 
          Md2hdbv2 <- rbind(cbind(.ZtWZwrapper(locX.Re,w.resid), t(hessnondiag)),
                            cbind(hessnondiag, - locd2hdv2)) 
          lad <- - LogAbsDetWrap(Md2hdbv2,logfac=-log(2*pi))
        } ## le lad est OK = ladbv de la version standard
        #likranu <- ( LogAbsDetWrap(longcorrSig,logfac=log(2*pi)) + u_h %*% corrinvSig %*% u_h)/2
        likranu <- ( blocksize* LogAbsDetWrap(corrSig,logfac=log(2*pi)) + u_h %*% longcorrinvSig %*% u_h)/2
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
      blob <- .calc_latentL(COVcorr)
      loc_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) 
      loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 
      next_LMatrix <- .makelong(blob$u,longsize=ncol(ZAlist[[rt]])) ## il faut updater pour estimer les ranef correctement...
      attr(next_LMatrix,"Lcompact") <- blob$u ## kept for updating in next iteration and for output
      attr(next_LMatrix,"par") <- optr$par ## kept for updating in next iteration and for output
      thisranef <- attr(ZAlist,"ranefs")[rt]
      attr(thisranef,"type") <- attr(attr(ZAlist,"ranefs"),"type")[rt] ## ie simply "(.|.)": LMatrix with such type is for random slope ## ajout 2015/06
      attr(next_LMatrix,"ranefs") <- thisranef
    } else next_LMatrix <- NULL
    next_LMatrices[[rt]] <- next_LMatrix
  } ## loop on rt = ranefs
  return(list(next_LMatrices=next_LMatrices,next_lambda_est=loc_lambda_est,
              optr_par=optr$par))
} ## end def MakeCovEst2
