cat("\ntest spaMM::corMatern for lme():")

data("blackcap")
blackcapD <-cbind(blackcap,dummy=1) ## obscure, isn't it? 
require(nlme)
## With method= 'ML' in lme, The correlated random effect is described 
##  as a correlated residual error and no extra residual variance is fitted:
bf <- lme(fixed = migStatus ~ means, data = blackcapD, random = ~ 1 | dummy, 
    correlation = corMatern(form = ~ longitude+latitude  | dummy), 
    method = "ML")

testthat::expect_equal(logLik(bf)[[1]],-7.941674,tolerance=1e-6)
testthat::expect_equal(exp((bf$modelStruct$corStruct)[[1]]),18.3595917,tolerance=1e-5) ## range =1/rho 
testthat::expect_equal(exp((bf$modelStruct$corStruct)[[2]]),0.62857441,tolerance=1e-5) ## nu (not range as comment said for a long time )
# parametrisation was modified in bf$modelStruct$corStruct from log(range =1/rho) to log(rho) 
