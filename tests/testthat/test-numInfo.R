cat(crayon::yellow("\ntest computations of numerical Information matrix:")) 

data("wafers")
lmmfit <- fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),data=wafers)
(numinfo <- numInfo(lmmfit,FALSE))
crit <- diff(range(numinfo- numInfo(lmmfit,TRUE)))
testthat::test_that("numInfo(.,FALSE)= numInfo(.,TRUE)",
                    testthat::expect_true(crit<1e-10))
crit <- max(abs(sqrt(diag(solve(numinfo))[3:9]) - summary(lmmfit,verbose=FALSE)$beta_table[,"Cond. SE"]))
testthat::test_that("numInfo() consistent with cond.SEs",
                    testthat::expect_true(crit<2e-8))

lmmfit <- fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),data=wafers, method="REML") # variances are huge
numinfo <- numInfo(lmmfit,transf=FALSE)
crit <- diff(range(numinfo- numInfo(lmmfit,transf=TRUE)))
testthat::test_that("numInfo(.,FALSE)= numInfo(.,TRUE) (REML)",
                    testthat::expect_true(crit<1e-10)) # works maybe by 'chance', the gradient of marginal lik wrt all params still being ~0 
crit <- max(abs(sqrt(diag(solve(numInfo(lmmfit,transf=FALSE,which="beta")))) - summary(lmmfit,verbose=FALSE)$beta_table[,"Cond. SE"]))
testthat::test_that("numInfo() consistent with cond.SEs (REML)",
                    testthat::expect_true(crit<1e-8))

if (TRUE) { # Largely undocumented variants  for fixing arguments
  data("blackcap")
  if (FALSE) { # This worked for dispFn(NULL) not failing (returning numeric(0), presumably ignored)
    (mrf <- fitme(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=1),data=blackcap,
                  method="ML",fixed=list(phi=1,lambda=1,hyper=list("1"=list(hy_kap=1)))) )
    numInfo(mrf)
  }
  
  (mrf <- fitme(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=1),data=blackcap,
                method="ML",fixed=list(phi=1,hyper=list("1"=list(hy_lam=1,hy_kap=1)))) )
  numInfo(mrf)
  
  (mrf <- fitme(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=1),data=blackcap,
                method="ML",fixed=list(phi=1,hyper=list("1"=list(hy_trL=.dispFn(1),hy_trK=-4.59512)))) ) # .dispFn() is exported
  numInfo(mrf)

  if (FALSE) { # But this one does not work whatever the dispFn:
    (mrf <- fitme(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=1),data=blackcap,
                  method="ML",fixed=list(phi=1,lambda=1,corrPars=list("1"=list(kappa=1)))) )
  }
  
  # HLCor: this happens to work: 
  (mrf <- HLCor(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=1),data=blackcap,
                HLmethod="ML",ranPars=list(phi=1,lambda=1,corrPars=list("1"=list(kappa=1)))) )
  numInfo(mrf)

}

