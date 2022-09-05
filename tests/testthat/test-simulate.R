cat(crayon::yellow("\ntest-simulate:\n"))
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
#
fit1 <- HLfit(x1~1+(1|x2),data=data.frame(x1=c(1,1),x2=c(1,2)), ranFix=list(lambda=1),family=poisson)
sim1 <- simulate(fit1, nsim=3, seed=1, type="predVar", variances=list(predVar=TRUE)) # test of .calc_invV_factors() for ZAfix=ZAL= Identity

cat("simulate for Gamma response:\n") 
set.seed(123)
gr <- data.frame(y=rgamma(1000,shape=9/2,scale=2/3)) # mean mu=3=exp(1.0986), variance=2 => phi=2/9
# Here fitme uses HLfit methods which provide cond. SE for phi by default:
(gamfit <- fitme(y~1,data=gr,family=Gamma(log)))
var(simulate(gamfit,type="residual")) ## must approach 2

{ # Without newdata, etaFix and offset() should ive equivalent results
  set.seed(1)
  data_resp <- data.frame(y = rnorm(100), x = rnorm(100), ID = gl(100, 10))
  fake_fit_etaFix <- fitme(y ~ 0+offset(2+0.5*x) + (1|ID), data = data_resp,
                    fixed = list(lambda = 5, phi = 10))
  fake_fit_offset <- fitme(y ~ x + (1|ID), data = data_resp,
                    etaFix = list(beta = c("(Intercept)" = 2, x = 0.5)),
                    fixed = list(lambda = 5, phi = 10))
  set.seed(123)
  s1 <- simulate(fake_fit_etaFix)
  set.seed(123)
  s2 <- simulate(fake_fit_offset)
  diff(range(s2-s1))
  diff(range(get_predVar(fake_fit_etaFix)-get_predVar(fake_fit_offset)))
  # LM:
  fake_fit_etaFix <- fitme(y ~ 0+offset(2+0.5*x), data = data_resp,
                           fixed = list(lambda = 5, phi = 10))
  fake_fit_offset <- fitme(y ~ x, data = data_resp,
                           etaFix = list(beta = c("(Intercept)" = 2, x = 0.5)),
                           fixed = list(lambda = 5, phi = 10))
  set.seed(123)
  s1 <- simulate(fake_fit_etaFix)
  set.seed(123)
  s2 <- simulate(fake_fit_offset)
  diff(range(s2-s1))
  diff(range(get_predVar(fake_fit_etaFix)-get_predVar(fake_fit_offset)))
  
}


