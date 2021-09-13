cat(crayon::yellow(" -> test-composite-nested: "))

{
  data("blackcap")
  # MLdistMat2 <- as.matrix(proxy::dist(blackcap[,c("latitude","longitude")]))
  MLcorMat2 <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                          nu=0.6285603,rho=0.0544659)
  cap_mv <- blackcap
  cap_mv$name <- as.factor(rownames(blackcap))                
  cap_mv$grp <- 1L+(blackcap$migStatus>1)   
  ch <- chol(proxy::as.matrix(MLcorMat2,diag=1))
  set.seed(123)
  v1 <- tcrossprod(ch,t(rnorm(14,sd=1)))
  v2 <- tcrossprod(ch,t(rnorm(14,sd=1)))
  cap_mv$status <- 2*v1 + v2
  cap_mv$status2 <- 2*v1 - v2
  # plot(cap_mv)
}

if (TRUE) {
  cat(crayon::yellow("check predict mv() not composite first (problem pre-v3.8.34)"))
  (basic_rC1 <- fitmv(submodels=list(mod1=list(status ~ 1+ (mv(1,2)|name), fixed=list(phi=0.1)),
                                    mod2=list(status2 ~ 1+ (mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                     data=cap_mv))
  (basic_rC2 <- fitmv(submodels=list(mod1=list(status ~ 1+ (0+mv(1,2)|name), fixed=list(phi=0.1)),
                                    mod2=list(status2 ~ 1+ (0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                     data=cap_mv))
  (p11 <- predict(basic_rC1)[,1])
  (p12 <- predict(basic_rC1, newdata=basic_rC1$data)[,1] )
  (p21 <- predict(basic_rC2)[,1])
  (p22 <- predict(basic_rC2, newdata=basic_rC2$data)[,1] )
  crit <- diff(range(p11-p12,p11-p21,p21-p22))
  testthat::test_that(paste0("predict mv(): criterion was ",signif(crit,4)," >1e-7"),
                      testthat::expect_true(crit<1e-7) )
  
  # Preliminary independent fit test, pure ranCoefs
  zuta <- fitme(status ~ 1+ (1|name), #verbose=c(TRACE=TRUE),
                data=cap_mv, fixed=list(lambda=1,phi=0.1),
                control.HLfit=list(sparse_precision=FALSE))
  zutb <- fitme(status2 ~ 1+ (1|name), fixed=list(lambda=1,phi=0.1), #verbose=c(TRACE=TRUE),
                data=cap_mv, 
                control.HLfit=list(sparse_precision=FALSE))
  (zut0 <- fitmv(submodels=list(mod1=list(status ~ 1+ (0+mv(1,2)|name), fixed=list(phi=0.1)),
                                mod2=list(status2 ~ 1+ (0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                 data=cap_mv, 
                 fixed=list(ranCoefs=list("1"=c(1,0,1))),
                 control.HLfit=list(sparse_precision=FALSE)))
  logLik(zut0)-(logLik(zuta)+logLik(zutb))
  
}



cat(crayon::yellow("; checks corrMatrix composite"))

#  OK independent fit test, trivial ranCoefs + corrMatrix
zuta <- fitme(status ~ 1+ corrMatrix(1|name), #verbose=c(TRACE=TRUE),
              data=cap_mv, corrMatrix=MLcorMat2, fixed=list(lambda=1,phi=0.1),
              control.HLfit=list(sparse_precision=FALSE))
zutb <- fitme(status2 ~ 1+ corrMatrix(1|name), fixed=list(lambda=1,phi=0.1), #verbose=c(TRACE=TRUE),
              data=cap_mv, corrMatrix=MLcorMat2, 
              control.HLfit=list(sparse_precision=FALSE))
(zut0d <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                              mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
               data=cap_mv, corrMatrix=MLcorMat2, 
               fixed=list(ranCoefs=list("1"=c(1,0,1))),
               control.HLfit=list(sparse_precision=FALSE)))
logLik(zut0d)-(logLik(zuta)+logLik(zutb))
(zut0s <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                              mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
               data=cap_mv, corrMatrix=MLcorMat2, 
               fixed=list(ranCoefs=list("1"=c(1,0,1))),
               control.HLfit=list(sparse_precision=TRUE)))
logLik(zut0s)-(logLik(zuta)+logLik(zutb)) 

# and now nontrivial ranCoefs:
(zut1d <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                              mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
               data=cap_mv, corrMatrix=MLcorMat2, 
               fixed=list(ranCoefs=list("1"=c(1,0.5,1))),
               control.HLfit=list(sparse_precision=FALSE)))
(zut1s <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                              mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
               data=cap_mv, corrMatrix=MLcorMat2, 
               fixed=list(ranCoefs=list("1"=c(1,0.5,1))),
               control.HLfit=list(sparse_precision=TRUE))) 
predict(zut1d)[,1] - predict(zut1s)[,1]
logLik(zut1s)-logLik(zut1d) 

ranef(zut1s, type="uncorrelated") # distinct chol factors, distinct u_h 
ranef(zut1d, type="uncorrelated")
r1 <- ranef(zut1s)[[1]] # but the correlated ranefs are i"dentical ! up to differences in the latentL_blob$d !" (0+mv()) important here !
r2 <- ranef(zut1d)[[1]]
crit <- diff(range(r1-r2))
testthat::test_that(paste0("ranef corrMatrix(mv()...): criterion was ",signif(crit,4)," >1e-7"),
                    testthat::expect_true(crit<1e-7) )

# and now variable ranCoefs:
(zut1d <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                               mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                data=cap_mv, corrMatrix=MLcorMat2, 
                control.HLfit=list(sparse_precision=FALSE)))
(zut1s <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                               mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                data=cap_mv, corrMatrix=MLcorMat2, 
                control.HLfit=list(sparse_precision=TRUE))) 
logLik(zut1s)-logLik(zut1d) 
(p1d <- predict(zut1d))
(p2d <- predict(zut1d, newdata=zut1d$data))
(p1s <- predict(zut1s))
(p2s <- predict(zut1s, newdata=zut1s$data)) 
(crit <- diff(range(p1d-p2d,p1d-p1s,p1s-p2s)))
testthat::test_that(paste0("ranef corrMatrix(mv()...): criterion was ",signif(crit,4)," >1e-7"),
                    testthat::expect_true(crit<1e-7) )

####################################################################"


if (FALSE) {
  ZUT1 <- function(fit1, caveat="assuming lambda=1") {
    Lmat <- t(chol(proxy::as.matrix(MLcorMat2,diag=1))) # 
    
    tcrossprod(get_ZALMatrix(fit1,force_bind = TRUE)) - tcrossprod(Lmat) 
    
    I0ZX <- rbind(cbind(diag(x=14L),0*fit1$X.pv),
                  cbind(Lmat,fit1$X.pv) )
    w <- diag(sqrt(c(rep(1,14),rep(1/fit$phi[[1]],14))))
    
    wI0ZX <-w %*% I0ZX
    
    w_y <- w %*% c(rep(0,14),cap_mv$status)            
    
    Matrix::solve(crossprod(wI0ZX), crossprod(wI0ZX,w_y))[,1]
    
  }
  
  ZUT <- function(fit, caveat="assuming lambda=1") {
    if (how(fit)$MME_method=="AUGI0_ZX_sparsePrecision") {
      if (FALSE) {
        Lmat <- t(Matrix::kronecker(chol(attr(fit$strucList[[1]],"latentL_blob")$compactcovmat),
                                    t(as(Cholesky(as(proxy::as.matrix(MLcorMat2,diag=1),"sparseMatrix")),"sparseMatrix")))) # 
        Lmat <- t(Matrix::kronecker(chol(attr(fit$strucList[[1]],"latentL_blob")$compactcovmat),
                                    (solve(as(Cholesky(as(solve(proxy::as.matrix(MLcorMat2,diag=1)),"sparseMatrix")),"sparseMatrix"))))) # 
      }
      Lmat <- t(solve(BLOB$chol_Q)) # 
      if (FALSE) {
        # mais alors Md2hdv2 shoudl be 
        crossprod(Lmat)/fit$phi[[1]] +diag(28)
        # la version "chol_Q" de Lmat marche
        Lmat <- t(solve(BLOB$chol_Q)) # 
        crossprod(Lmat)/fit$phi[[1]] +diag(28) - solve(BLOB$chol_Q) %*% (BLOB$Gmat) %*% solve(t(BLOB$chol_Q))
        get_ZALMatrix(fit,force_bind = TRUE) - Lmat 
        solve(crossprod(BLOB$factor_inv_Md2hdv2)) - crossprod(Lmat)/fit$phi[[1]] +diag(28)
      }
    } else {
      Lmat <- t(kronecker(chol(attr(fit$strucList[[1]],"latentL_blob")$compactcovmat),chol(proxy::as.matrix(MLcorMat2,diag=1)))) # 
    }
    if (FALSE) {
      # BLOB$logdet_sqrt_d2hdv2 is 22.35205
      (logdet <- Matrix::determinant(crossprod(Lmat)/fit$phi[[1]] +diag(28), logarithm = TRUE)$modulus/2)
      # - get_from_MME(sXaug,"logdet_sqrt_d2hdv2") + n_u_h*log(2*pi)/2 est 3.378226
      # bref la logLik est
      fit$APHLs$clik + sum(dnorm(fit$v_h,log = TRUE))
      
      
      as.matrix(tcrossprod(get_ZALMatrix(fit,force_bind = TRUE))) - tcrossprod(Lmat) 
      
      
      # sauf que les crossprod ne sont pas unique même si les tcrossprod sont identiques
      solve(BLOB$chol_Q) %*% BLOB$chol_Q %*% (crossprod(Lmat)/fit$phi[[1]] +diag(28)) %*% t(BLOB$chol_Q) %*% solve(t(BLOB$chol_Q))
      solve(BLOB$chol_Q) %*% ( BLOB$chol_Q %*% (crossprod(Lmat)/fit$phi[[1]]) %*% t(BLOB$chol_Q)+
                                 BLOB$chol_Q %*% (diag(28)) %*% t(BLOB$chol_Q)) %*% solve(t(BLOB$chol_Q))
      solve(BLOB$chol_Q) %*% ( BLOB$chol_Q %*% (crossprod(Lmat)/fit$phi[[1]]) %*% t(BLOB$chol_Q)+
                                 precisionMatrix) %*% solve(t(BLOB$chol_Q))
      
      
      solve(BLOB$chol_Q) %*% (BLOB$Gmat) %*% solve(t(BLOB$chol_Q))
    }
    
    I0ZX <- rbind(cbind(diag(x=28L),0*fit$X.pv),
                  cbind(Lmat,fit$X.pv) )
    w <- diag(sqrt(c(rep(1,28),rep(1/fit$phi[[1]],28))))
    
    wI0ZX <-w %*% I0ZX
    
    w_y <- w %*% c(rep(0,28),cap_mv$status,cap_mv$status2)            
    
    Matrix::solve(crossprod(wI0ZX), crossprod(wI0ZX,w_y))[,1]
    
  }
  ZUT(zut0d)
  ZUT(zut1s)
  
}

