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

if (FALSE) { 
  (lm_mix <- fitme(formula=y ~ 1, family=gaussian(), 
                   resid.model= list(formula= ~  1+(1|batch)),
                   data=wafers))
  numInfo(lm_mix) # seems OK, but *the same model* with a different parametrisation of the intercept:
  
  # logically expected and documented issues with mixed-effect residual-dispersion model
  (loglm_mix <- fitme(formula=y ~ 1, family=gaussian(log), 
                   resid.model= list(formula= ~  1+(1|batch)),
                   data=wafers))
  numInfo(loglm_mix)
  # numInfo(lm_mix, refit_hacks = list(verbose=c(TRACE=TRUE))) shows the same higher-lik as in numInfo(loglm_mix) 
  # but the gradient is larger on the log scale => detection by numInfo only in th latter case.
  
  # numInfo fixes beta values (indeed, these are the only ones considered here)
  # Using outer-beta estimation then exhibits an underlying issue that happens with fixed beta's (here simplified again as gaussian() and fixing other parameters): 
  (outer_lm_mix <- fitme(formula=y ~ 1, family=gaussian(), 
                   resid.model= list(formula= ~  1+(1|batch), fixed=list(lambda=VarCorr(lm_mix$resid_fit)[1,3],etaFix=list(beta=fixef(lm_mix$resid_fit)))),
                   init=list(beta=fixef(lm_mix)),
                   data=wafers, verbose=c(phifit=2)))
  # => symptom: outer optim leads to 1220.8641 and final refit to -1220.873.
  # Only the v_h's of the resid model are inferred, and the beta and resid-v_h are not jointly consistent with each other (since outer beta-> no new beta from new resid-v_h)
  # The  initial likelihood = the final, and logL rises above in between (only the final refit brings it back. But this final refit does not use outer-beta)
  
  # Note also that may not work on outer-beta fits: numInfo(outer_lm_mix)  (__F I X M E___)
  
  # outer-beta results are insensitive to the family link... which makes sense.
  (outer_loglm_mix <- fitme(formula=y ~ 1, family=gaussian(log), 
                         resid.model= list(formula= ~  1+(1|batch), fixed=list(lambda=VarCorr(loglm_mix$resid_fit)[1,3],etaFix=list(beta=fixef(lm_mix$resid_fit)))),
                         init=list(beta=fixef(loglm_mix)),
                         data=wafers, verbose=c(phifit=2)))
}


if (TRUE) { # Largely undocumented variants  for fixing arguments
  data("blackcap")
  if (FALSE) { # This worked for .dispFn(NULL) not failing (returning numeric(0), presumably ignored)
    # But now .dispFn(NULL) fails (intentionally) so this previously undocumented syntax is just no longer feasible. 
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

if (spaMM.getOption("example_maxtime")>0.7) { # not an actual timing
  if(requireNamespace("lme4", quietly = TRUE)) {     # check that partially fixing of ranCoefs is taken into account
    data("sleepstudy", package="lme4")
    spfit <- fitme(Reaction ~ Days + (Days|Subject), sleepstudy, fixed=list(ranCoefs=list("1"=c(NA,0.1,NA))))
    numInfo(spfit, check_deriv = TRUE) # the warning is not strictly necessary in this case, so _slightly_ confusing (despite being formally correct). 
    numInfo(spfit, check_deriv = FALSE)
  }
}



