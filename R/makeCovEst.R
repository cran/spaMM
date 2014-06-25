makeCovEst <-
function(u_h,ZAlist,cum_n_u_h,X_lamres,prev_LMatrices,userLfixeds,lev_lambda) {
  nrand <- length(ZAlist)
  next_LMatrices <- list()
  Xi_cols <- attr(X_lamres,"Xi_cols")
  cum_Xi_cols <- attr(X_lamres,"cum_Xi_cols") 
  Lu <- u_h
  next_lambda_est <- numeric(length(u_h))
  for (rt in seq_len(length(ZAlist))) {
    ## estimate correlation matrix
    Xi_ncol <- Xi_cols[rt]
    blocksize <- ncol(ZAlist[[rt]])/Xi_ncol 
    ## cov mat of u_h if not fixed by user ## standard REML method
    if ( Xi_ncol>1 && ! userLfixeds[rt]) {
      COVpredUnderHuncorr <- matrix(0,ncol=Xi_ncol,nrow=Xi_ncol) ## var on diag, corr outside diag
      for (xit in seq_len(Xi_ncol)) {
        col1 <- cum_Xi_cols[rt]+xit
        urange1 <- which(X_lamres[,col1]==1L)
        unique.lambda <- sum(Lu[urange1]^2)/sum(1-lev_lambda[urange1]) ## NOT in linkscale 
        unique.lambda <- pmax(unique.lambda,1e-8) # as for phi
        unique.lambda <- pmin(unique.lambda,.SpaMM$maxLambda)  
        # if (verbose["trace"]) {print(paste("lambda=",signif(unique.lambda,4)),quote=F)}
        COVpredUnderHuncorr[xit,xit] <- unique.lambda
        for(xjt in seq_len(xit-1L)) {
          col2 <- cum_Xi_cols[rt]+xjt
          urange2 <- which(X_lamres[,col2]==1L)
          ## not fully intuitve but cov leverages were tested:
          unique.cov <- sum(Lu[urange1]*Lu[urange2])/sqrt(sum(1-lev_lambda[urange1])*sum(1-lev_lambda[urange2])) 
          ## makes sure it is feasible despite numerical inacc 
          maxcov <- sqrt(prod(diag(COVpredUnderHuncorr)[c(col1,col2)]))
          unique.cov <- sign(unique.cov) * min(abs(unique.cov),maxcov)
          COVpredUnderHuncorr[xit,xjt] <- COVpredUnderHuncorr[xjt,xit] <- unique.cov 
        }
      } 
      ## deduce the correlation matrix of L u_h 
      if (!is.null(prev_LMatrices)) {
        prevL <- attr(prev_LMatrices[[rt]],"Lcompact")
        COVcorr <-  prevL %*% COVpredUnderHuncorr %*% t(prevL)
      } else COVcorr <- COVpredUnderHuncorr
      ## we need torepresent this as a design matrix times the variances of uncorrelated ranef. The neat solution is:
      blob <- selfAdjointSolverCpp(COVcorr) ## COVcorr= blob$u %*% diag(blob$d) %*% t(blob$u) ## 
      u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
      next_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) ## doivent être ceux du modèle fitté par modèle augmenté => diag compact avant correction
      Lcompact <- blob$u  
      ## Build Lmatrix between all pairs of u_h (nr*nr) from parameter estimates (2*2) for the design matrix,
      ## and keep the compact form as attribute for updating this compact form in the next iteration, and for output
      next_LMatrix <- diag(ncol(ZAlist[[rt]]))
      for (it in seq_len(Xi_ncol)) {
        urange1 <- (it-1)*blocksize + seq(blocksize)
        diag(next_LMatrix)[urange1] <- Lcompact[it,it]
        for (jt in seq_len(it-1)) {
          urange2 <- (jt-1)*blocksize + seq(blocksize)
          diag(next_LMatrix[urange1,urange2]) <- Lcompact[it,jt]
          diag(next_LMatrix[urange2,urange1]) <- Lcompact[jt,it]
        }
      }
      attr(next_LMatrix,"ranefs") <- attr(ZAlist,"ranefs")[rt]
      attr(next_LMatrix,"Lcompact") <- Lcompact ## version compacte de next_LMatrix
    } else next_LMatrix <- NULL
    next_LMatrices[[rt]] <- next_LMatrix
  } ## loop on rt = ranefs
  return(list(next_LMatrices=next_LMatrices,next_lambda_est=next_lambda_est))
}
