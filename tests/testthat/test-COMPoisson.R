cat(crayon::yellow("\ntest of COMPoisson family:"))
## particularly useful to test permuted sparse algos

data("freight")
# estimation of nu in COMPoisson GLM (Sellers & Shmueli 2010):
(fitfit <- fitme(broken ~ transfers, data=freight, family = COMPoisson(),method="ML"))
crit <- diff(range(residVar(fitfit,"fam_parm"),5.781809)) # more public extractor: get_inits_from_fit(fitfit)$init$COMP_nu
testthat::test_that(paste0("criterion was ",signif(crit,4)," from 5.781809"), testthat::expect_true(crit<1e-6)) # decimals depend on optimizer|COMPoisson approxs

(fitfitlog <- fitme(broken ~ transfers, data=freight, family = COMPoisson(link="log"),method="ML"))
crit <- diff(range(residVar(fitfitlog,"fam_parm"),5.733456)) 
testthat::test_that(paste0("criterion was ",signif(crit,4)," from 5.781809"), testthat::expect_true(crit<1e-6)) # decimals depend COMPoisson approxs again

# COMPoissonReg::glm.cmp(broken ~ transfers, data=freight)
## Crude way to the same result:
if (spaMM.getOption("example_maxtime")>0.9) {  
  objfn <- function(nu) {  
    #cat(nu," ")
    fit <- HLfit(broken ~ transfers, data=freight, family = COMPoisson(nu=nu),HLmethod="ML")
    logLik(fit)
  }
  optr <- optim(1,objfn,lower=0.05,upper=10,method="L-BFGS-B",control=list(fnscale=-1))
  crit <- diff(range(optr$par,5.781816)) 
  testthat::test_that(paste0("criterion was ",signif(crit,4)," from 5.781809"), testthat::expect_true(crit<1e-6)) # decimals depend COMPoisson approxs again
}
# GLMM with under-dispersed conditional response
compmm <- HLfit(broken ~ transfers+(1|id), data=freight, family = COMPoisson(nu=10),HLmethod="ML")
crit <- diff(range(compmm$lambda[[1]],0.4573305))
if (spaMM.getOption("fpot_tol")>0) {
  testthat::test_that(paste0("criterion was ",signif(crit,6)," from 0.4573305"), testthat::expect_true(crit<1e-6))
} else testthat::expect_true(crit<1e-6)
# Also test simulate() through test-DHARMa.R
