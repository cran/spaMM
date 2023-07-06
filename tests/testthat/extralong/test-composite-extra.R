cat(crayon::yellow("\ntest-composite-extra:\n"))

cat(crayon::yellow("checks AR1 composite"))

{ # test different algebras
  ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package
  
  (compAR1fitde <-  fitme(lh ~ 1 + AR1(time|time), data=ts, # control.HLfit=list(algebra="spcorr"), 
                          lower=list(phi=1e-4),
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  if ( ! compAR1fitde$how$MME_method[1]=="sXaug_EigenDense_QRP_Chol_scaled") { # yes, dense...
    stop("default MME_method has changed...")
  }
  (compAR1fit <-  fitme(lh ~ 1 + AR1(time|time), data=ts, control.HLfit=list(algebra="spcorr"), 
                        lower=list(phi=1e-4),
                        fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (compAR1fitsp <-  fitme(lh ~ 1 + AR1(time|time), data=ts, control.HLfit=list(algebra="spprec"), 
                          #lower=list(phi=1e-4, ranCoefs=list("1"=c(0.001,0.0002))), # with ad hoc lower ranCoefs... 
                          # The ad-hoc 'fixed' and/or 'lower' to avoid singularities (ARphi->1, phi->0) and associated numerical imprecisions
                          lower=list(phi=1e-4),
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (crit <- diff(range(c(logLik(compAR1fitde),logLik(compAR1fit),logLik(compAR1fitsp),-28.37433))))
  FIXME <- testthat::test_that(paste0("Whether the three algebras give consistent results for AR1(time|time): crit= ",signif(crit,4)," >1e-05"),
                               testthat::expect_true(crit<1e-05) )
  if ( ! FIXME) doSeeMe("Do see me!") 
}

cat(crayon::yellow("; checks ARp composite"))

{ # test different algebras
  ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package
  
  (compARpfitsp <-  fitme(lh ~ 1 + ARp(time|time, p=1), data=ts,  
                          lower=list(phi=1e-4),
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  if ( ! compARpfitsp$how$MME_method[1]=="AUGI0_ZX_spprec") { 
    stop("default MME_method has changed...")
  }
  (compARpfitde <-  fitme(lh ~ 1 + ARp(time|time), data=ts, control.HLfit=list(algebra="decorr"), 
                          lower=list(phi=1e-4),
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (compARpfit <-  fitme(lh ~ 1 + ARp(time|time), data=ts, control.HLfit=list(algebra="spcorr"), 
                        lower=list(phi=1e-4),
                        fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (crit <- diff(range(c(logLik(compARpfitde),logLik(compARpfit),logLik(compARpfitsp),-28.37433))))
  FIXME <- testthat::test_that(paste0("Whether the three algebras give consistent results for AR1(time|time): crit= ",signif(crit,4)," >1e-05"),
                               testthat::expect_true(crit<1e-05) )
  if ( ! FIXME) doSeeMe("Do see me!") 
}

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
  cat(crayon::yellow("; check predict mv() not composite first (problem pre-v3.8.34)"))
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

{ # test different algebras
  ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package
  
  (compMatfitde <-  fitme(lh ~ 1 + Matern(time|time), data=ts, # control.HLfit=list(algebra="spcorr"), 
                          lower=list(phi=1e-4),
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  if ( ! compMatfitde$how$MME_method[1]=="sXaug_EigenDense_QRP_Chol_scaled") { # yes, dense...
    stop("default MME_method has changed...")
  }
  # unscutable numeric changes with Matrix devel v1.6.0 => premature stop, prevented by setting xtol_rel=1e-07 or 4e-6 but not 5e-6.
  (compMatfit <-  fitme(lh ~ 1 + Matern(time|time), data=ts, control.HLfit=list(algebra="spcorr"), 
                        lower=list(phi=1e-4),
                        fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (compMatfitsp <-  fitme(lh ~ 1 + Matern(time|time), data=ts, control.HLfit=list(algebra="spprec"), 
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (crit <- diff(range(c(logLik(compMatfitde),logLik(compMatfit),logLik(compMatfitsp),-26.43478))))
  FIXME <- testthat::test_that(paste0("Whether the three algebras give consistent results for Matern(time|time): crit= ",signif(crit,4)," >1e-05"),
                               testthat::expect_true(crit<1e-05) )
  if ( ! FIXME) doSeeMe("Do see me!") 
}

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
               fixed=list(ranCoefs=list("1"=c(1,0.5,1))), # meaning dependent on using (0+mv()).
               control.HLfit=list(sparse_precision=TRUE))) 
predict(zut1d)[,1] - predict(zut1s)[,1]
logLik(zut1s)-logLik(zut1d) 

ranef(zut1s, type="uncorrelated") # distinct chol factors, distinct u_h ...
ranef(zut1d, type="uncorrelated")
r1 <- ranef(zut1s)[[1]] # ... but the correlated ranefs are identical.
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
testthat::test_that(paste0("ranef corrMatrix(mv()...): criterion was ",signif(crit,4)," >5e-7"),
                    testthat::expect_true(crit<5e-7) )

{  cat(crayon::yellow( "Matern(LHS |<nested RHS>) and more corrMatrix(LHS |<nested RHS>)... " ))
  {
    data("blackcap")
    toy <- blackcap
    toy$ID <- gl(7,2)
    grp <- rep(1:2,7)
    toy$migStatus <- toy$migStatus +(grp==2)
    toy$loc <- rownames(toy) # to use as levels matching the corrMatrix dimnames
    
    toy$grp <- factor(grp)
    toy$bool <- toy$grp==1L
    toy$boolfac <- factor(toy$bool)
    toy$num <- seq(from=1, to=2, length.out=14)
    
    ## Build a toy corrMatrix as perturbation of identity matrix:
    n_rhs <- 14L
    eps <- 0.1
    set.seed(123)
    rcov <- ((1-eps)*diag(n_rhs)+eps*rWishart(1,n_rhs,diag(n_rhs)/n_rhs)[,,1])
    # eigen(rcov)$values
    colnames(rcov) <- rownames(rcov) <- toy$loc # DON'T FORGET NAMES
  }
  
  {
    ## The unconstrained fit fits to correlation = -1 and vanishing phi by augmented-y algo, both of which complicates comparisons
    # (fit1yaug <- fitme(migStatus ~ grp + Matern(grp|latitude %in% grp), data=toy, fixed=list(ranCoefs=list("1"=c(NA,-0.5,NA)))))
    ## => Force non y-augmented algo and then consistently fitted phi across calls:
    (fit1 <- fitme(migStatus ~ grp + Matern(grp|latitude %in% grp), data=toy, init=list(phi=1e-5), fixed=list(ranCoefs=list("1"=c(NA,-0.5,NA)))))
    
    { # to test for handling of levels found only once in the data (and next, predict() checks found more issues):
      ttoy <- toy[c(1,2*seq(7)),]
      fit1t <- fitme(migStatus ~ grp + Matern(grp|latitude %in% grp), data=ttoy, init=list(phi=1e-5), fixed=list())
      predict(fit1t) 
      predict(fit1t, newdata=ttoy) 
      predict(fit1t, newdata=ttoy[rev(seq(8)),])[rev(seq(8)),] 
    }
    
    rC <- get_ranPars(fit1)$ranCoefs[[1]]
    corMat <- Corr(fit1)[[1]][1:14,1:14]/rC[1] 
    colnames(corMat) <- rownames(corMat) <- toy$loc 
    # 
    (fit1fix_nest <- fitme(migStatus ~ grp + corrMatrix(grp|loc %in% grp), data=toy, init=list(phi=1e-5), corrMatrix=corMat, fixed=list(ranCoefs=list("1"=rC))))
    (fit1fix <- fitme(migStatus ~ grp + corrMatrix(grp|loc), data=toy, init=list(phi=1e-5), corrMatrix=corMat, fixed=list(ranCoefs=list("1"=rC)))) 
    
    testthat::expect_true(diff(range(logLik(fit1),logLik(fit1fix_nest),logLik(fit1fix)))<3e-09) 
    
  }
  
}


####################################################################"


if (FALSE) { # that was devel scratch, not tidy tests
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
    if (.is_spprec_fit(fit)) {
      if (FALSE) {
        Lmat <- t(Matrix::kronecker(chol(attr(fit$strucList[[1]],"latentL_blob")$compactcovmat),
                                    t(as(Cholesky(as(proxy::as.matrix(MLcorMat2,diag=1),"sparseMatrix")),"CsparseMatrix")))) # 
        Lmat <- t(Matrix::kronecker(chol(attr(fit$strucList[[1]],"latentL_blob")$compactcovmat),
                                    (solve(as(Cholesky(as(solve(proxy::as.matrix(MLcorMat2,diag=1)),"sparseMatrix")),"CsparseMatrix"))))) # 
      }
      Lmat <- t(solve(BLOB$chol_Q)) # 
      if (FALSE) {
        # mais alors Md2hdv2 should be 
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

