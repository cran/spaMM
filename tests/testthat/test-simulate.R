cat(cli::col_yellow("\ntest-simulate:\n"))

{
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
  
  # check that ZAL is used in marginal simulation (here, further, with new locations):
  set.seed(123)
  (check_marg_ZAL <- simulate(HLC, newdata=Loaloa[c(1:3,1:3)+0.1,], sizes=HLC$BinomialDen[1:6]))
  crit <- unique(check_marg_ZAL - c(53L, 73L, 68L, 23L, 52L, 55L)) # a different result
  # has pointed to a problem in devel code wrt to selection of column 
  # at point controlled by 'cols_from_RHS' arg of .compute_ZAXlist() in .wrap_compute_ZALlist4simulate(). 
  testthat::expect_equal(crit,0) 
  
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
}

{ # Without newdata, etaFix and offset() should give equivalent results
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

if (spaMM.getOption("example_maxtime")>2) { ## NA handling
  {
    set.seed(123)
    ssize <- 10L
    shape <- 0.35
    toydata <- data.frame(
      yellow=rbinom(ssize, 16, prob=rbeta(ssize,shape,shape)),
      purple=rbinom(ssize, 16, prob=rbeta(ssize,shape,shape)), # (purple ignored below)
      blue=rbinom(ssize, 16, prob=rbeta(ssize,shape,shape)),
      phenotype=rnorm(ssize),
      grp=seq(ssize)
    )
    
    toydataNA <- toydata
    toydataNA$pheno1 <- toydata$phenotype
    toydataNA$pheno2 <- toydataNA$pheno1 + rnorm(10L,sd = 0.3)
    toydataNA$pheno3 <- toydataNA$pheno1 + rnorm(10L,sd = 0.3)
    toydataNA$yellow[1] <- NA 
    toydataNA$pheno1[2] <- toydataNA$pheno2[2] <- NA 
  }
  
  { # fitme
    foo <- simulate(mini <- fitme(yellow ~ pheno3, family = poisson(), data=toydataNA)) # 9 info lines
    ## NOT mv: the object's $data no longer contain lines with missing response values. 
    ## simulate() produces 9 values matching the 9 lines in these $data (OK for update_resp)
    crit <- length(foo)==9L
    testthat::test_that("simulate() with NA OK for simple GLM", testthat::expect_true(crit))
    # mininull <- fitme(yellow ~ 1, family = poisson(), data=toydataNA) # 
    # LRT(mini,mininull, boot.repl=3L) 
    #
    # same update_resp()-OK logic although many details differ (a pheno1 value is also missing)
    foo <- simulate(sglmm <- fitme(yellow ~ pheno1+(1|grp), family = poisson(), data=toydataNA)) # nrow(ZAL) *=* nrow(new_X_ZACblob$newX.pv)
    crit <- length(foo)==8L
    testthat::test_that("simulate() with NA OK for simple GLMM", testthat::expect_true(crit))
    # simulate(sglmm, newdata=sglmm$data) # 8 values since 2 lines removed from sglmm$data
    simulate(sglmm, newdata=toydataNA) # 9 values
    
    # By contrast, in mv fits, the $data retain all lines, so simulate() should fill with NA values
    foo <- simulate(
      fitmv(submodels = list(
        list(yellow ~ pheno1, family = poisson()),
        list(blue ~ pheno2, family = poisson()),
        list(purple ~ pheno3, family = poisson())),
        data = toydataNA)) # note that nobs attribute does not match (it should match when newdata present) 
    crit <- length(foo)==30L
    testthat::test_that("simulate() with NA OK for simple mv-GLM", testthat::expect_true(crit))
    crit <- length(na.omit(foo))==27L
    testthat::test_that("simulate() with NA OK for simple mv-GLM", testthat::expect_true(crit))
  }
  
  { # fitmv
    P3 <- fitmv(submodels = list(
      list(yellow ~ pheno1+(1|grp), family = poisson()),
      list(blue ~ pheno2, family = poisson()),
      list(purple ~ pheno3, family = poisson())),
      data = toydataNA) 
    foo <- simulate(P3) 
    crit <- length(foo)==30L
    testthat::test_that("simulate() with NA OK for mv-GLMM", testthat::expect_true(crit))
    crit <- length(na.omit(foo))==27L
    testthat::test_that("simulate() with NA OK for mv-GLMM", testthat::expect_true(crit))
    
    foo <- simulate(P3, newdata=P3$data) 
    crit <- length(foo)==30L
    testthat::test_that("simulate() with NA OK for mv-GLMM", testthat::expect_true(crit))
    # no predicted u for grp=1, but marginal simulation must be possible for grp=1:
    crit <- length(na.omit(foo))==28L
    testthat::test_that("simulate() with NA OK for mv-GLMM", testthat::expect_true(crit))
    # LRT(P3, P3, boot.repl = 3L)$bootInfo$bootreps
  }
  
  { # pois4mlogit
    foo <- simulate(byP3f <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + pheno1, family = poisson()),
      list(blue ~ offset(.dynoffset) + pheno2, family = poisson()),
      list(purple ~ offset(.dynoffset) + 0, family = poisson())),
      data = toydataNA, types=c("yellow","blue","purple"))) 
    crit <- length(foo)==30L
    testthat::test_that("simulate() with NA OK for p4m", testthat::expect_true(crit))
    # there are NA's for some simulated responses and this is OK for bootstraps:
    # LRT(byP3f, byP3f, boot.repl = 3L)$bootInfo$bootreps
    
    byP3m <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + pheno1+(1|grp), family = poisson()),
      list(blue ~ offset(.dynoffset) + pheno2+(1|grp), family = poisson()),
      list(purple ~ offset(.dynoffset) + 0, family = poisson())),
      progress=2*interactive(),               
      fixed=list(lambda=0.01), # makes it faster for the test (fitted value would be 2.475)
      data = toydataNA, types=c("yellow","blue","purple"))
    
    length(pp1 <- predict(byP3m, verbose=c(na=FALSE))) # 26 (default na.action=na.omit) 
    length(pp1 <- predict(byP3m, verbose=c(na=FALSE), na.action=na.exclude)) # 30
    rowSums(m1 <- matrix(pp1, ncol=3), na.rm=TRUE) # 0 for the row of NA's
    length(pp2 <- predict(byP3m, newdata=byP3m$data,
                          verbose=c(na=FALSE))) # 27 bc only second draw cannot be predicted
    # max(abs(range(matrix(pp2, ncol=3)-m1))) # it's not lear what I expected here. The first two rows differ.
    plot_effects(byP3m, "pheno1", submodel=1)
    
    foo <- simulate(byP3m, newdata=byP3m$data, 
                    sizes=get_drawSizes(byP3m,p4m="M")) 
    crit <- length(foo)==30L
    testthat::test_that("simulate(., newdata) with NA OK for p4m", testthat::expect_true(crit))
    crit <- length(na.omit(foo))==27L
    testthat::test_that("NA OK in simulate(., newdata) with NA for p4m", testthat::expect_true(crit))
    
    foo <- simulate(byP3m) 
    crit <- length(foo)==30L
    testthat::test_that("simulate() with NA OK for p4m", testthat::expect_true(crit))
    crit <- length(na.omit(foo))==26L
    testthat::test_that("NA OK in simulate() with NA for p4m", testthat::expect_true(crit))
    # there are NA's for some simulated responses and this is OK for bootstraps:
    # LRT(byP3m, byP3m, boot.repl = 3L)$bootInfo$bootreps
  }

}
