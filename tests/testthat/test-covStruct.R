cat(cli::col_yellow("\ntest-covStruct: two precision matrices\n"))

{ # two corrMatrix'es, both precision: permutation check (T. Tsang, 08/2024)
  # derived from pedigree example
  A <- matrix(NA, ncol=6,nrow=6)
  A[lower.tri(A,diag=TRUE)] <- c(8,0,4,4,4,2, 8,4,0,2,5, 8,2,5,4.5, 8,5,2.5, 9,5.5, 9)/8
  A <- Matrix::forceSymmetric(A,uplo = "L")
  colnames(A) <- rownames(A) <- 1:6
  
  ## data simulation
  cholA <- chol(A)
  varU <- 0.4; varE <- 0.6; rep <- 20
  n <- rep*6
  set.seed(108)
  bStar <- rnorm(6, sd=sqrt(varU))
  b <- crossprod(as.matrix(cholA),bStar)
  ID <- rep(1:6, each=rep)
  e0 <- rnorm(n, sd=sqrt(varE))
  y <- b[ID]+e0
  obs <- data.frame(y=y,IDgen=ID,IDenv=ID) ## two copies of ID for readability of GLMM results
  
  obs$y01 <- ifelse(y<1.3,0,1)
  prec_mat <- solve(A)
  colnames(prec_mat) <- rownames(prec_mat) <- rownames(A) # important
  B <- A*0.05  #generate another covariance matrix
  diag(B) <- diag(as.matrix(A)) #only the off-diagonals of B are different with A
  
  prec_mat_B <- solve(B) #name it again
  colnames(prec_mat_B) <- rownames(prec_mat_B) <- rownames(A) # important
  
  (fitBA <- fitme(y01 ~ 1+corrMatrix(1|IDenv)+ corrMatrix(1|IDgen) , 
                  covStruct=list(precision=prec_mat_B,precision=prec_mat),
                  data=obs, family=binomial(), method="REML"))
  
  (fitAB <- fitme(y01 ~ 1+corrMatrix(1|IDgen)+corrMatrix(1|IDenv), 
                  covStruct=list(precision=prec_mat,precision=prec_mat_B),
                  data=obs, family=binomial(), method="REML"))
  
  testthat::expect_true(diff(range(logLik(fitAB),logLik(fitBA)))<1e-08) # post v4.5.0
  
}

