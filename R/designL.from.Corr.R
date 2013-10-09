designL.from.Corr <-
function(m,try.chol=TRUE,try.eigen=FALSE,threshold=1e-06,debug=FALSE) {
  if (try.chol) {
     ## warning("(only debugging message): try.chol is true in designL.from.Corr")
     L <- try(chol(m),silent=TRUE) ## dim -> unique geo coordinates
     if (class(L)=="try-error") { ## attempt at regularization 050213
       mreg <- m *(1-1e-8)
       diag(mreg) <- diag(mreg) + 1e-8 
       L <- try(chol(mreg),silent=TRUE) 
     }
  }
  if ( (! try.chol) || class(L)=="try-error" ) { ## default : eigen
    ## slower by more stable. Not triangular but should not be a pb
    if (try.eigen) {
      LDL <- try(eigen(m,symmetric=TRUE),silent=TRUE) ## may _hang_ in R2.15.2 on nearly-I matrices
      if (class(LDL) != "try-error") e.val <- LDL$values 
    }
    if ( (! try.eigen) || class(LDL)=="try-error" || any(e.val < -1e-08)) {
      SVD <- try(svd(m)) 
      if (class(SVD)=="try-error") {
        print("Singular value decomposition failed.") ## can still try other package like the svd package...
        return(try(stop(),silent=T)) ## passes control to calling function
      } else {
        e.val <- SVD$d
        if (any(SVD$d< -1e-08)) {
          mess <- pastefrom("correlation matrix has suspiciously large negative eigenvalue(s).")
          print(mess)
          return(try(stop(),silent=T)) ## passes control to calling function
        } else { ## we have a not-too-suspect SVD decomp
          e.val[e.val< threshold]<- threshold
          ## must be valid for sym (semi) PD matrices using U, V being eigenvectors of m %*% t(m)
          L <- sweep(SVD$u,2,sqrt(SVD$d),`*`); L <- L %*% t(SVD$v)  ##L <- SVD$u %*% sqrt(diag(SVD$d)) %*% t(SVD$v)
        }      
      } 
    } else { ## we have an LDL decomp
      e.val[e.val< threshold]<- threshold
      L <- sweep(LDL$vectors,2,sqrt(e.val),`*`); L <- L %*% t(LDL$vectors) ## L <- LDL$vectors %*% diag(sqrt(e.val)) %*% t(LDL$vectors)    ## symmetric
    }
    if (debug) {
       cholL <- t(chol(m))
       err <- (L %*% t(L) - cholL %*% t(cholL))
       if (any(is.na(err))) stop("NA in designL...")
       if (max(abs(err))> 1e-04 ) stop("max(abs(err))> 1e-04 in designL...")
       print(c(min(L),max(L),min(cholL),max(cholL)))
    } 
  } else { ## exists L from chol
    L<-t(L) ## IMPORTANT ! 
  } ## such that L %*% t(L) = the correlation matrix in both cases
  colnames(L) <- colnames(m)
  rownames(L) <- colnames(m) ## for checks in HLfit
  return(L)
}
