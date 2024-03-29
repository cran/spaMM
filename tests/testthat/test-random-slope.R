cat(crayon::yellow("\ntest of random-slope model:"))
if(requireNamespace("lme4", quietly = TRUE)) {
  data("sleepstudy",package = "lme4")
  (res <- HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy)) 
  testthat::expect_equal(logLik(res),c(p_bv=-871.8141),tolerance=2e-4) 
  # lambda from lambda.object$lambda_list potentially dependent on internal representation of the cov mat
  # which itself depends on a try_chol and whether chol facto succeeds... but the following should be consistent:
  verif_lam <- testthat::expect_equal(VarCorr(res)[[2,3]],35.07163,tolerance=1e-5) # 35.07155 by spprec
  try(testthat::expect_equal(predict(res)[1:6,],predict(res,newdata=sleepstudy[1:6,],re.form= ~ (Days|Subject))[,1])) 
  ## tests of predVar with re.form and ranCoefs in test-devel-predVar-ranCoefs.R
  ## a bit slow but detects many problem: (+ effect of refit$lambda)
  if (FALSE) { # \label{ares_bobyqa_spprec} # this case produced a bug before found in .calc_cov_from_ranCoef (+ bobyqa margin is used in the bobyqa_wrapper + ).
    #spaMM.options(allow_augZXy=TRUE) # make sure of this; even if it's the default
    oldopt <- spaMM.options(optimizer="bobyqa", sparse_precision=TRUE)
    (ares <- fitme(Reaction ~ Days + AR1(1|Days) + (Days|Subject), data = sleepstudy, verbose=c(TRACE=TRUE)))
    spaMM.options(oldopt)
  }
  (ares <- fitme(Reaction ~ Days + AR1(1|Days) + (Days|Subject), data = sleepstudy, verbose=c(TRACE=FALSE)))  ## AR-lambda is ~0 hence lik is flat  wrt ARphi
  testthat::expect_equal(logLik(ares),c(p_v=-875.969672803))
  if (spaMM.getOption("example_maxtime")>1.5) { ## approx time v2.4.129 ten times faster than v2.3.33
    sm <- fitme(Reaction ~ Days + (Days|Subject) + Matern(1|Days), fixed=list(nu=0.5),data = sleepstudy)
    testthat::expect_true(diff(range((c(logLik(sm),logLik(ares)))))<1e-5)
    # some syntax checks. The fitme() fits have an addition AR1 term.
    spaMM.options(sparse_precision = FALSE)
    HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy,ranFix =list(ranCoefs=list("1"=c(612.1,0.06555,35.07)),phi=654.9))
    fitme(Reaction ~ Days + (Days|Subject) + AR1(1|Days), method="REML",
          fixed=list(lambda=c(NA,10),ranCoefs=list("1"=c(612.1,0.06555,35.07)),phi=654.9),
          data = sleepstudy,verbose=c(TRACE=FALSE))
    chkfx <- fitme(Reaction ~ Days + (Days|Subject), method="REML",
          fixed=list(ranCoefs=list("1"=c(612.666,NA,NA)),phi=654.9),
          data = sleepstudy)
    testthat::expect_equal(get_ranPars(chkfx)$ranCoefs[[1]][1],612.666)
    chkfx <- fitme(Reaction ~ Days + (Days|Subject), method="REML",
          fixed=list(ranCoefs=list("1"=c(612.666,0.06666,NA)),phi=654.9),
          data = sleepstudy,verbose=c(TRACE=FALSE))
    testthat::expect_equal(get_ranPars(chkfx)$ranCoefs[[1]][2],0.06666)
    fitme(Reaction ~ Days + (Days|Subject) + AR1(1|Days), method="REML",
          fixed=list(lambda=c(NA,10),ranCoefs=list("1"=c(NA,0,NA)),phi=654.9),
          data = sleepstudy,verbose=c(TRACE=FALSE))
    
    spaMM.options(sparse_precision = TRUE) 
    HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy,ranFix =list(ranCoefs=list("1"=c(612.1,0.06555,35.07)),phi=654.9))
    fitme(Reaction ~ Days + (Days|Subject) + AR1(1|Days), method="REML",
          fixed=list(lambda=c(NA,10),ranCoefs=list("1"=c(612.1,0.06555,35.07)),phi=654.9),
          data = sleepstudy,verbose=c(TRACE=FALSE)) 
    fitme(Reaction ~ Days + (Days|Subject) + AR1(1|Days), method="REML",
          fixed=list(lambda=c(NA,10),ranCoefs=list("1"=c(612.1,0.06555,35.07)),phi=654.9),
          data = sleepstudy,verbose=c(TRACE=FALSE))
    chkfx <- fitme(Reaction ~ Days + (Days|Subject), method="REML",
                   fixed=list(ranCoefs=list("1"=c(612.666,NA,NA)),phi=654.9),
                   data = sleepstudy)
    testthat::expect_equal(get_ranPars(chkfx)$ranCoefs[[1]][1],612.666)
    chkfx <- fitme(Reaction ~ Days + (Days|Subject), method="REML",
                   fixed=list(ranCoefs=list("1"=c(612.666,0.06666,NA)),phi=654.9),
                   data = sleepstudy,verbose=c(TRACE=FALSE))
    testthat::expect_equal(get_ranPars(chkfx)$ranCoefs[[1]][2],0.06666)
    spaMM.options(sparse_precision = NULL)
  }
} else {
  cat( "package 'lme4' not available, cannot run random-slope test.\n" )
}

if (FALSE) {
  spaMM.options(sparse_precision = TRUE)
  fitme(Reaction ~ Days + (Days|Subject), data = sleepstudy,
        fixed=list(phi=700,ranCoefs=list("1"=c(600,0.25,30)))) ## FIXME does not work without list()
  spaMM.options(sparse_precision = NULL)
}