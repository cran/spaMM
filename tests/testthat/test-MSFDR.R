cat(crayon::yellow("\ntest of MSFDR:"))

if (spaMM.getOption("example_maxtime")>1.5) { ## fast but not routinely useful
  data("wafers")
  nullfit <- fitme(y~1+(1|batch), data=wafers,family=Gamma(log))
  fullfit <- fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), data=wafers, family=Gamma(log))
  fdr <- MSFDR(nullfit=nullfit,fullfit=fullfit,verbose=FALSE)
  if (how(fdr,verbose=FALSE)$obsInfo) {
    testthat::expect_equal(logLik(fdr), c(P_v=-1198.4051377))
  } else testthat::expect_equal(logLik(fdr), c(p_v=-1198.411))
}
