cat("\ntest of random-slope model:")
if(require("lme4", quietly = TRUE)) {
  res <- HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy)
  testthat::expect_equal(res$APHLs$p_bv,-871.8141,tolerance=2e-4)
  # next value depends on .calc_latentL() implementation
  testthat::expect_equal(res$lambda[[2]],34.91168,tolerance=2e-4) 
  testthat::expect_equal(predict(res)[1:6,],predict(res,newdata=sleepstudy[1:6,],re.form= ~ (Days|Subject))[,1])
  if (spaMM.getOption("example_maxtime")>14) { ## approx time v2.2.80
    as <- fitme(Reaction ~ Days + AR1(1|Days) + (Days|Subject), data = sleepstudy) 
    sm <- fitme(Reaction ~ Days + (Days|Subject) + Matern(1|Days), fixed=list(nu=0.5),data = sleepstudy)
    testthat::expect_true(diff(range((c(logLik(sm),logLik(as)))))<1e-6)
    spaMM.options(sparse_precision = FALSE)
    HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy,ranFix =list(ranCoefs=list("1"=c(612.1,0.06555,35.07)),phi=654.9))
    fitme(Reaction ~ Days + (Days|Subject) + AR1(1|Days), method="REML",
          fixed=list(lambda=c(NA,10),ranCoefs=list("1"=c(612.1,0.06555,35.07)),phi=654.9),
          data = sleepstudy,verbose=c(TRACE=TRUE))
    spaMM.options(sparse_precision = TRUE) 
    HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy,ranFix =list(ranCoefs=list("1"=c(612.1,0.06555,35.07)),phi=654.9))
    fitme(Reaction ~ Days + (Days|Subject) + AR1(1|Days), method="REML",
          fixed=list(lambda=c(NA,10),ranCoefs=list("1"=c(612.1,0.06555,35.07)),phi=654.9),
          data = sleepstudy,verbose=c(TRACE=TRUE)) 
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