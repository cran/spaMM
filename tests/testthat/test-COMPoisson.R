cat("\ntest of COMPoisson family:")
## particularly useful to test permuted sparse algos

data(freight)
# greedy estimation of nu (solution nu=5.78...,cf Sellers & Shmueli 2010):
objfn <- function(nu) {  
  fit <- HLfit(broken ~ transfers, data=freight, family = COMPoisson(nu=nu),HLmethod="ML")
  logLik(fit)
}
optr <- optim(1,objfn,lower=0.05,upper=10,method="L-BFGS-B",control=list(fnscale=-1))
testthat::expect_equal(optr$par,5.781816,tolerance=2e-4)
# GLMM with under-dispersed conditional response
hlfit <- HLfit(broken ~ transfers+(1|id), data=freight, family = COMPoisson(nu=10),HLmethod="ML")
testthat::expect_equal(hlfit$lambda[[1]],0.4573136,tolerance=2e-4)
