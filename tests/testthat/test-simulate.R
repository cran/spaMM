cat(crayon::yellow("\ntest-simulate\n"))
data("Loaloa")
HLC <- HLCor(cbind(npos,ntot-npos)~Matern(1|longitude+latitude),
             data=Loaloa,family=binomial(),
             ranPars=list(lambda=1,nu=0.5,rho=1/0.7)) 
testthat::expect_equal(dim(simulate(HLC, nsim=2)),c(197,2)) ## matrix 2 col -> OK
testthat::expect_equal(dim(simulate(HLC, nsim=1)),NULL) ## vector length N -> OK
testthat::expect_equal(dim(simulate(HLC, type = "residual", nsim=2)),c(197,2)) ## matrix 2 col -> OK
testthat::expect_equal(dim(simulate(HLC, type = "residual", nsim=1)),NULL) ## vector length N -> OK
testthat::expect_equal(dim(simulate(HLC, type = "predVar", nsim=2,
                                    variances=list(linPred=TRUE, disp=FALSE,cov=TRUE))),c(197,2)) ## matrix 2 col -> OK
testthat::expect_equal(dim(simulate(HLC, type = "predVar", nsim=1,
                                    variances=list(linPred=TRUE, disp=FALSE,cov=TRUE))),NULL) ## vector length N -> OK


cat("simulate on trivial toy examples:\n") 
x0 <- c(-1, 1)
var(x0)
fit0 <- HLfit(x0~1,data=data.frame(x0=x0)) 
vcov(fit0)
sim0 <- simulate(fit0, nsim=10000, seed=1) # ignores uncertainty
var(t(sim0))
sim0 <- simulate(fit0, nsim=10000, seed=1, type="predVar", variances=list(predVar=TRUE)) # Accounting for uncertainty in fixed effects (cf Spencer Graves, R-devel, 2019/12/28)
var(t(sim0))
x1 <- 1 
fit1 <- HLfit(x1~1,data=data.frame(x1=x1), family=poisson)
fixef(fit1)
exp(fixef(fit1))
vcov(fit1)
sim1 <- simulate(fit1, 10000, 1)  # ignores uncertainty
var(t(sim1))
sim1 <- simulate(fit1, nsim=10000, seed=1, type="predVar", variances=list(predVar=TRUE))
var(t(sim1)) # \approx 6 

