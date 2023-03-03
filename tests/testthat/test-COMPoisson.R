cat(crayon::yellow("\ntest of COMPoisson family:"))
## particularly useful to test permuted sparse algos

data("freight")
# estimation of nu in COMPoisson GLM (Sellers & Shmueli 2010):
(fitfit <- fitme(broken ~ transfers, data=freight, family = COMPoisson(),method="ML"))
crit <- diff(range(residVar(fitfit,"fam_parm"),5.781798)) 
try(testthat::test_that(paste0("criterion was ",signif(crit,4)," from 5.781798"), testthat::expect_true(crit<1e-6))) # decimals depend on optimizer|COMPoisson approxs

(fitfitlog <- fitme(broken ~ transfers, data=freight, family = COMPoisson(link="log"),method="ML"))
crit <- diff(range(residVar(fitfitlog,"fam_parm"),5.733443)) 
try(testthat::test_that(paste0("criterion was ",signif(crit,4)," from 5.733443"), testthat::expect_true(crit<1e-6))) # decimals depend COMPoisson approxs again

# COMPoissonReg::glm.cmp(broken ~ transfers, data=freight)
## Crude way to the same result:
if (spaMM.getOption("example_maxtime")>0.9) {  
  objfn <- function(nu) {  
    #cat(nu," ")
    fit <- HLfit(broken ~ transfers, data=freight, family = COMPoisson(nu=nu), method="ML")
    logLik(fit)
  }
  optr <- optim(1,objfn,lower=0.05,upper=10,method="L-BFGS-B",control=list(fnscale=-1))
  crit <- diff(range(optr$par,5.781804 )) 
  try(testthat::test_that(paste0("criterion was ",signif(crit,4)," from 5.781804"), testthat::expect_true(crit<1e-6))) # decimals depend COMPoisson approxs again
}
# GLMM with under-dispersed conditional response
(compmm <- HLfit(broken ~ transfers+(1|id), data=freight, family = COMPoisson(nu=10), method="ML"))
crit <- diff(range(compmm$lambda[[1]],0.4573305))
try(  testthat::test_that(paste0("criterion was ",signif(crit,6)," from 0.4573305"), testthat::expect_true(crit<1e-6)) )
# Also test simulate() through test-DHARMa.R

# test of COMPoisson obsInfo; fixing lambda to nonzero value bc fitted zero value makes test weak.
(compmmobs <- fitme(broken ~ transfers+(1|id), data=freight, family = COMPoisson(link="log", nu=1), fixed=list(lambda=0.1), method=c("ML","obs"))) # using obsInfo algos
# that's nu=1 ~poisson so result should be equiv to poisson(log) whether exp or obs Info is used for either family
(pmm <- fitme(broken ~ transfers+(1|id), data=freight, family = poisson(), fixed=list(lambda=0.1), method=c("ML","obs"))) 
(compmmexp <- fitme(broken ~ transfers+(1|id), data=freight, family = COMPoisson(link="log", nu=1), fixed=list(lambda=0.1), method=c("ML","exp"))) # using obsInfo algos
testthat::test_that(paste0("whether obsInfo COMPoisson(link=log, nu=1) ~ poisson()"), 
                    testthat::expect_true(diff(range(logLik(pmm),logLik(compmmobs),logLik(compmmexp)))<1e-6)) 

                    