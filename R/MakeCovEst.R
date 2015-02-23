# HLfit(Reaction ~ 0 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ 1 + (Days|Subject), data = sleepstudy,HLmethod="ML")
# HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy,HLmethod="ML")


## experimental version of MakeCovEst stored in MakeCovEst.R.txt
makeCovEst2 <- function(u_h,ZAlist,cum_n_u_h,X_lamres,prev_LMatrices,
                       userLfixeds,hessUL,hessFac,w.resid,processed,prevZAL,clik) {
  nrand <- length(ZAlist)
  X.Re <- processed$X.Re
  based2hdv2 <- - t(prevZAL) %*% diag(w.resid) %*% prevZAL
  next_LMatrices <- prev_LMatrices
  Xi_cols <- attr(X_lamres,"Xi_cols")
  cum_Xi_cols <- attr(X_lamres,"cum_Xi_cols") 
  Lu <- u_h
  loc_lambda_est <- numeric(length(u_h))
  for (rt in seq_len(length(ZAlist))) {
    ## estimate correlation matrix 
    Xi_ncol <- Xi_cols[rt]
    blocksize <- ncol(ZAlist[[rt]])/Xi_ncol 
    ## cov mat of u_h if not fixed by user ## standard REML method 
    if ( Xi_ncol>1 && ! userLfixeds[rt]) {
      COVpredUnderHuncorr <- matrix(0,ncol=Xi_ncol,nrow=Xi_ncol) ## var on diag, corr outside diag
      compactLv <- matrix(0,nrow=Xi_ncol,ncol=Xi_ncol)
      lowerbloc <- lower.tri(compactLv,diag=TRUE) ## a matrix of T/F !
      ##
      ## Build Lmatrix between all pairs of u_h (nr*nr) from parameter estimates (2*2) for the design matrix,
      makelong <- function(Lcompact) {
        longLv <- diag(ncol(ZAlist[[rt]])) ## declaration
        for (it in seq_len(Xi_ncol)) {
          urange1 <- (it-1)*blocksize + seq(blocksize)
          diag(longLv)[urange1] <- Lcompact[it,it]
          for (jt in seq_len(it-1)) {
            urange2 <- (jt-1)*blocksize + seq(blocksize)
            diag(longLv[urange1,urange2]) <- Lcompact[it,jt]
            diag(longLv[urange2,urange1]) <- Lcompact[jt,it]
          }
        }
        longLv
      } ## end def makelong
      ##
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
      if (ncol(X.Re)>0L) hessnondiag <- crossprod(prevZAL,sweep(X.Re,MARGIN=1,w.resid,`*`))  
      objfn <- function(parvec) {
        compactcovmat <- makeLcovLt(parvec)
        blob <- selfAdjointSolverCpp(compactcovmat) ## COVcorr= blob$u %*% diag(blob$d) %*% t(blob$u) ## u is unitary !!
        ## cosmetic / interpretative permutation
        blib <- Matrix::expand(Matrix::lu(t(blob$u))) ## use pivoting in lu as a useful permutation...
        blob <- list(u=t(as.matrix(with(blib,L %*% U))),d=blob$d[blib$P@perm]) 
        prevL <- attr(prev_LMatrices[[rt]],"Lcompact")
        corrSig <- compactcovmat
        if ( ! is.null(prevL)) corrSig <- t(prevL) %*% corrSig %*% prevL ## corrects by what is absorbed in constant ZAL
        longcorrSig <- makelong(corrSig) ## car t(prevL) =solve(prevL)
        solvecompact <- blob$u %*% diag(1/blob$d) %*% t(blob$u) ## plus stable
        corrinvSig <- solvecompact
        if ( ! is.null(prevL)) corrinvSig <- t(prevL) %*% corrinvSig %*% prevL
        longcorrinvSig <- makelong(corrinvSig) ## car t(prevL) =solve(prevL)
        # equiv ligne suivante et une repres diagonale car identique a pre/ post par matrice unitaire (! important que unitaire !!!! )
        locd2hdv2 <- based2hdv2 - longcorrinvSig ## 
        # pour une repres non diagonal je devrais reconstruire une Lmatrix... (je l'ai  !)
        if (ncol(X.Re)==0L) { ## fit ML: p_bv=p_v hence d2hdpbv reduces to d2hdv2
          lad <- - LogAbsDetWrap( - locd2hdv2,logfac=-log(2*pi))
        } else { 
          Md2hdbv2 <- rbind(cbind(ZtWZ(X.Re,w.resid), t(hessnondiag)),
                            cbind(hessnondiag, - locd2hdv2)) 
          lad <- - LogAbsDetWrap(Md2hdbv2,logfac=-log(2*pi))
        } ## le lad est OK = ladbv de la version standard
        #likranu <- ( LogAbsDetWrap(longcorrSig,logfac=log(2*pi)) + u_h %*% corrinvSig %*% u_h)/2
        likranu <- ( blocksize* LogAbsDetWrap(corrSig,logfac=log(2*pi)) + u_h %*% longcorrinvSig %*% u_h)/2
        if (is.nan(lad) || is.infinite(lad)){## because of determinant of nearly singular matrix
          zut <- abs(eigen(longcorrSig*(2*pi),only.values = T)$values) 
          zut[zut<1e-12] <- 1e-12
          likranu <- sum(log(zut)) ## L-BFGS-B requires a non trivial value
        }
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
        lowerb[2,1] <-   -0.99
        upperb[2,1] <-   0.99
        init <- (upperb+lowerb)/2
        diag(init) <- 0
        init <- init[lowerbloc]        
      } else {
        lowerb[2,1] <-  (-99+init[2])/100 
        upperb[2,1] <-   (99+init[2])/100          
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
      blob <- selfAdjointSolverCpp(COVcorr) ## COVcorr= blob$u %*% diag(blob$d) %*% t(blob$u) ## 
      blib <- Matrix::expand(Matrix::lu(t(blob$u))) ## use pivoting in lu as a useful permutation...
      blob <- list(u=t(as.matrix(with(blib,L %*% U))),d=blob$d[blib$P@perm]) ## + jolies façon de permuter $d ?
      loc_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) 
      loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 
      Lcompact <- blob$u  #
      next_LMatrix <- makelong(Lcompact) ## il faut updater pour estimer les ranef correctement...
      attr(next_LMatrix,"Lcompact") <- Lcompact ## kept for updating in next iteration and for output
      attr(next_LMatrix,"par") <- optr$par ## kept for updating in next iteration and for output
      attr(next_LMatrix,"ranefs") <- attr(ZAlist,"ranefs")[rt]
    } else next_LMatrix <- NULL
    next_LMatrices[[rt]] <- next_LMatrix
  } ## loop on rt = ranefs
  return(list(next_LMatrices=next_LMatrices,next_lambda_est=loc_lambda_est,
              latest.unique.cov=optr$par[2]))
} ## end def MakeCovEst2
