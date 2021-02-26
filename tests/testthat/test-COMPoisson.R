cat(crayon::yellow("\ntest of COMPoisson family:"))
## particularly useful to test permuted sparse algos

data("freight")
# estimation of nu in COMPoisson GLM (Sellers & Shmueli 2010):
(fitfit <- fitme(broken ~ transfers, data=freight, family = COMPoisson(),method="ML"))
testthat::expect_equal(environment(fitfit$family$aic)$nu,5.781809,tolerance=1e-6) ## FIXME need a better extractor ## decimals depend on optimzer
# COMPoissonReg::glm.cmp(broken ~ transfers, data=freight)
## Crude way to the same result:
if (spaMM.getOption("example_maxtime")>2.23) {
  objfn <- function(nu) {  
    #cat(nu," ")
    fit <- HLfit(broken ~ transfers, data=freight, family = COMPoisson(nu=nu),HLmethod="ML")
    logLik(fit)
  }
  optr <- optim(1,objfn,lower=0.05,upper=10,method="L-BFGS-B",control=list(fnscale=-1))
  testthat::expect_equal(optr$par,5.781816,tolerance=1e-6)
}
# GLMM with under-dispersed conditional response
compmm <- HLfit(broken ~ transfers+(1|id), data=freight, family = COMPoisson(nu=10),HLmethod="ML")
crit <- diff(range(compmm$lambda[[1]],0.4573305))
if (spaMM.getOption("fpot_tol")>0) {
  testthat::test_that(paste0("criterion was ",signif(crit,6)," from 0.4573305"), testthat::expect_true(crit<1e-6))
} else testthat::expect_true(crit<1e-6)
# Also test simulate() through test-DHARMa.R
