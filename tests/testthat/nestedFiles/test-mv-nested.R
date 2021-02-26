cat(crayon::yellow(" -> test-mv-nested:"))

library(spaMM)
options(error=recover)
{
  data("wafers")
  me <- fitme(y ~ 1+(1|batch), family=Gamma(log), data=wafers)
  set.seed(123)
  y2 <- simulate(me, type="residual")
  wafmv <- wafers
  wafmv$batch2 <- wafmv$batch
  wafmv$y2 <- y2
  wafmv$ly <- log(wafmv$y)
  wafmv$y3 <- log(y2)
  (fitme(y3 ~ 1+(1|batch), family=gaussian(log), data=wafmv))
  
  (zut1 <- fitmv(submodels=list(mod1=list(formula=ly ~ 1+(1|batch), family=gaussian()),
                         mod2=list(formula=y3 ~ 1+(1|batch2), family=gaussian())), 
                 data=wafmv))
  testthat::expect_true(diff(range( predict(zut1, newdata=zut1$data)-predict(zut1)))<1e-14)
  testthat::expect_true(diff(range( get_predVar(zut1, newdata=zut1$data)-get_predVar(zut1)))<1e-14) ## there a resid.model so nothing is done with phi
  simulate(zut1)

  cat(crayon::yellow("[ upper, get_ranPars() (-> VarCorr()) ] ; "))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=ly ~ 1+(1|batch), family=gaussian()),
                         mod2=list(formula=y3 ~ 1+(1|batch), family=gaussian())), 
                 data=wafmv, init=list(lambda=1), upper=list(lambda=0.02)))
  testthat::expect_true(diff(range(get_ranPars(zut1, which="lambda"),0.02))<1e-12) 
  
  cat(crayon::yellow("(mv()|.) and (0+mv()|.); "))
  (zut0 <- fitmv(submodels=list(mod1=list(formula=ly~X1+(0+mv(1,2)|batch)),
                         mod2=list(formula=y3~X1+(0+mv(1,2)|batch), family=gaussian())), 
                 data=wafmv))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=ly~X1+(mv(1,2)|batch)),
                         mod2=list(formula=y3~X1+(mv(1,2)|batch), family=gaussian())), 
                 data=wafmv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=y3~X1+(mv(1,2)|batch)),
                         mod1=list(formula=ly~X1+(mv(1,2)|batch))), 
                 data=wafmv))
  testthat::expect_true(diff(range(logLik(zut0),logLik(zut1),logLik(zut2)))<1e-05) 
  #  # test equivalence of the two paremtrization, and that they are indeed distinguished by the code. 
  testthat::expect_true(diff(range(ranef(zut0)[[1]][,2]-rowSums(ranef(zut1)[[1]]))/(max(abs(ranef(zut0)[[1]]))))<1e-4)
  
  { # missing data
    afers <- wafmv
    afers$y3[71:140] <- NA
    afers$ly[1:70] <- NA
    (zut0 <- fitmv(submodels=list(mod1=list(formula=ly~X1+(0+mv(1,2)|batch)),
                           mod2=list(formula=y3~X1+(0+mv(1,2)|batch), family=gaussian())), 
                   data=afers))
    (zut1 <- fitmv(submodels=list(mod1=list(formula=ly~X1+(mv(1,2)|batch)),
                           mod2=list(formula=y3~X1+(mv(1,2)|batch), family=gaussian())), 
                   data=afers))
    (zut2 <- fitmv(submodels=list(mod2=list(formula=y3~X1+(mv(1,2)|batch)),
                           mod1=list(formula=ly~X1+(mv(1,2)|batch))), 
                   data=afers))
    testthat::expect_true(diff(range(logLik(zut0),logLik(zut1),logLik(zut2)))<1e-05) 
    
  }
  
  cat(crayon::yellow("ranCoefs; "))
  (mod1 <- fitme(ly~X1, data=wafmv))
  (mod2 <- fitme(y3~X1+(X2|batch), data=wafmv))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=ly~X1),
                         mod2=list(formula=y3~X1+(X2|batch), family=gaussian())), 
                 data=wafmv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=y3~X1+(X2|batch), family=gaussian()),
                         mod1=list(formula=ly~X1)), 
                 data=wafmv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2), logLik(mod1)+logLik(mod2)))<1e-08) 
  
  (zut1 <- fitmv(submodels=list(mod1=list(formula=ly~X1+(X2|batch)),
                         mod2=list(formula=y3~X1+(X2|batch), family=gaussian())), 
                 data=wafmv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=y3~X1+(X2|batch), family=gaussian()),
                         mod1=list(formula=ly~X1+(X2|batch))), 
                 data=wafmv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2)))<1e-07) 
  
  cat(crayon::yellow("fixing ranCoefs in two ways.; "))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=ly~X1+(X2|batch), fixed=list(ranCoefs=list("1"=c(0.02, -0.1, 0.005)))),
                         mod2=list(formula=y3~X1+(X2|batch), family=gaussian())), 
                 data=wafmv))
  (zut2 <- fitmv(submodels=list(mod1=list(formula=ly~X1+(X2|batch)),
                         mod2=list(formula=y3~X1+(X2|batch), family=gaussian())), 
                 data=wafmv, fixed=list(ranCoefs=list("1"=c(0.02, -0.1, 0.005)))))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2)))<1e-08) 
  
  cat(crayon::yellow("with resid.model; "))
  # independent-fit test
  (mod1 <- fitme(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~ X3+I(X3^2), data=wafmv))
  (mod2 <- fitme(formula=y2 ~ 1+(1|batch2), family=Gamma(log), data=wafmv))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~ X3+I(X3^2)),
                                  mod2=list(formula=y2 ~ 1+(1|batch2), family=Gamma(log))), 
                          data=wafmv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=y2 ~ 1+(1|batch2), family=Gamma(log)),
                                  mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~ X3+I(X3^2))), 
                          data=wafmv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2), logLik(mod1)+logLik(mod2)))<1e-7)
  # permutation test 
  (zut1 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~ X3+I(X3^2)),
                                  mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log))), 
                          data=wafmv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log)),
                                  mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~ X3+I(X3^2))), 
                          data=wafmv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2)))<1e-10)
  
  cat(crayon::yellow("full dhglm; "))
  (mod1 <- fitme(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~ 1+(1|batch), data=wafmv))
  (mod2 <- fitme(formula=y2 ~ 1+(1|batch2), family=Gamma(log), data=wafmv))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~ 1+(1|batch)),
                                  mod2=list(formula=y2 ~ 1+(1|batch2), family=Gamma(log))), 
                          data=wafmv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=y2 ~ 1+(1|batch2), family=Gamma(log)),
                                  mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~ 1+(1|batch))), 
                          data=wafmv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2), logLik(mod1)+logLik(mod2)))<1e-4)
  # permutation test 
  (zut1 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~  1+(1|batch)),
                                  mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log))), 
                          data=wafmv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log)),
                                  mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log),resid.model= ~  1+(1|batch))), 
                          data=wafmv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2)))<1e-4)
  testthat::expect_true(diff(range( predict(zut1, newdata=zut1$data)-predict(zut1)))<1e-14)
  testthat::expect_true(diff(range( get_predVar(zut1, newdata=zut1$data)-get_predVar(zut1)))<1e-14) ## there a resid.model so nothin is done with phi
  spaMM_boot(zut1, function(v) var(v), nsim=3L, type ="marginal")$bootreps
  confint(zut1,"(Intercept)_1") # __FIXME__ messy display

  cat(crayon::yellow("fixing phi: four different ways; "))# 
  (zut1 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log), fixed=list(phi=0.001)),
                         mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log), fixed=list(phi=0.002))), 
                 data=wafmv, fixed=list(lambda=0.1))) 
  (zut2 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log)),
                         mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log))), 
                 data=wafmv, fixed=list(phi=list("1"=0.001,"2"=0.002), lambda=0.1)))
  (zut3 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log)),
                         mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log))), 
                 data=wafmv, fixed=list(phi=list(0.001,0.002), lambda=0.1)))
  (zut4 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log), fixed=list(phi=0.001)),
                         mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log))), 
                 data=wafmv, fixed=list(phi=list("2"=0.002), lambda=0.1)))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2),logLik(zut3),logLik(zut4)))<1e-14)
  
  cat(crayon::yellow("visual checks; "))# and low fixed phis are useful for the following visual checks:
  set.seed(123)
  ressim1 <- simulate(zut1, type="residual")
  ressim2 <- simulate(zut1, newdata=zut1$data, type="residual")
  margsim1 <- simulate(zut1)
  margsim2 <- simulate(zut1, newdata=zut1$data)
  plot(ressim1,ressim2);abline(0,1) # clusters on the diagonal
  plot(margsim1,margsim2);abline(0,1) # TWO indep draws of the 11 levels => clusters are not on the diagonal
  # the two submodels share the ranef values so the clusters are always on the diagonal, even for margsim:
  plot(ressim2[1:198], ressim2[1L+(1:198)]);abline(0,1) # on diagonal
  plot(margsim2[1:198], margsim2[1L+(1:198)]);abline(0,1) # on diagonal
  
  cat(crayon::yellow("confint with 'fixed'; "))# Checking that confint() obeys fixed values fixed in different ways
  (zute <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log)),
                         mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log))), 
                 data=wafmv, fixed=list(lambda=0.1))) 
  confint(zute,"(Intercept)_1")$lowerfit # estimated phis
  (zutf <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log), fixed=list(phi=0.2224)),
                         mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log), fixed=list(phi=0.2103))), 
                 data=wafmv, fixed=list(lambda=0.1))) 
  confint(zutf,"(Intercept)_1")$lowerfit # fixed phis
  (zutfm <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log), fixed=list(phi=0.2224)),
                          mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log))), 
                  data=wafmv, fixed=list(phi=list("2"=0.2103), lambda=0.1))) 
  confint(zutfm,"(Intercept)_1")$lowerfit # fixed phis
  (zutf2 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log)),
                          mod2=list(formula=y2 ~ 1+(1|batch), family=Gamma(log))), 
                  data=wafmv, fixed=list(phi=list("2"=0.2103), lambda=0.1))) 
  confint(zutf2,"(Intercept)_1")$lowerfit # fixed 2nd phi 
  
  
  cat(crayon::yellow("some of the early tests; "))# permutation test 
  (zut1 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log)),
                                  mod2=list(formula=y3 ~ 1+(1|batch), family=gaussian(log))), 
                          data=wafmv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=y3 ~ 1+(1|batch), family=gaussian(log)),
                                  mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log))), 
                          data=wafmv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2)))<1e-8)
  
  # permutation test. The init from hazards of development
  (zut1 <- fitmv(submodels=list(mod2=list(formula=y3 ~ 1+(1|batch), family=gaussian()),
                                  mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log))), 
                          control.HLfit=list(LevenbergM=TRUE),
                          data=wafmv,init=list(phi=list("1"=1 )))) 
  (zut2 <- fitmv(submodels=list(mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log)),
                                  mod2=list(formula=y3 ~ 1+(1|batch), family=gaussian())), 
                          control.HLfit=list(LevenbergM=TRUE),
                          data=wafmv,init=list(phi=list("2"=1 )))) 
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2)))<1e-10)
  (zut <- fitmv(submodels=list(mod2=list(formula=y3 ~ 1+(1|batch), family=gaussian()),
                                 mod1=list(formula=y ~ 1+(1|batch), family=Gamma(log))), 
                         data=wafmv,init=list(phi=list("1"=1,"2"=0.2 )))) 
  attr(zut,"optimInfo")$LUarglist$canon.init
  
}
  
{
  npos <- c(11,16,14,2,6,1,1,4,10,22,7,1,0,0,1,6)
  ntot <- c(36,20,19,16,17,11,5,6,37,32,19,17,12,10,9,7)
  treatment <- c(rep(1,8),rep(0,8))
  clinic <-c(seq(8),seq(8))
  clinics <- data.frame(npos=npos,nneg=ntot-npos,treatment=treatment,clinic=clinic)
  #
  (fitClinics <- HLfit(cbind(npos,nneg)~treatment+(1|clinic),family=binomial(),data=clinics))
  set.seed(123)
  y2 <- simulate(fitClinics, type="residual")
  climv <- clinics
  climv$np2 <- y2
  climv$nn2 <- climv$npos + climv$nneg - y2
  climv$clinic2 <- climv$clinic
  
  cat(crayon::yellow("binomial-poisson; "))# permutation test
  (zut1 <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                                  mod2=list(formula=np2~treatment+(1|clinic),family=poisson())), 
                          data=climv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=np2~treatment+(1|clinic),family=poisson()),
                                  mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial())), 
                          data=climv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2)))<1e-10)

  cat(crayon::yellow("predVar cov gaussian-poisson; "))
  # independent-fits test
  (fg <- fitme(formula=npos~treatment+(1|clinic),family=gaussian(),data=climv))
  (fp <- fitme(formula=np2~treatment+(1|clinic2),family=poisson(),data=climv))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=npos~treatment+(1|clinic),family=gaussian()),
                                mod2=list(formula=np2~treatment+(1|clinic2),family=poisson())), 
                 data=climv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=np2~treatment+(1|clinic2),family=poisson()),
                                mod1=list(formula=npos~treatment+(1|clinic),family=gaussian())), 
                 data=climv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2),logLik(fg)+logLik(fp)))<1e-10)
  testthat::expect_true(diff(range(get_predVar(zut1)-c(get_predVar(fg),get_predVar(fp))))<1e-5)
  testthat::expect_true(diff(range(get_predVar(zut2)-c(get_predVar(fp),get_predVar(fg))))<1e-5)
  testthat::expect_true(diff(range( get_predVar(zut1, newdata=zut1$data)-get_predVar(zut1)))<1e-14) 
  testthat::expect_true(diff(range( get_predVar(zut1, newdata=zut1$data,variances=list(cov=TRUE)) - 
                                      get_predVar(zut1,variances=list(cov=TRUE))))<1e-14) 
  testthat::expect_true(diff(range( get_predVar(zut1, newdata=zut1$data,variances=list(cov=TRUE)) - 
                                      Matrix::bdiag(get_predVar(fg,variances=list(cov=TRUE)),
                                            get_predVar(fp,variances=list(cov=TRUE)))))<1e-5) 
  get_intervals(zut1)
  get_intervals(zut1, intervals = "respVar")
  # permutation test
  (zut1 <- fitmv(submodels=list(mod1=list(formula=npos~treatment+(1|clinic),family=gaussian()),
                                mod2=list(formula=np2~treatment+(1|clinic),family=poisson())), 
                 data=climv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=np2~treatment+(1|clinic),family=poisson()),
                                mod1=list(formula=npos~treatment+(1|clinic),family=gaussian())), 
                 data=climv))
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2)))<1e-10)
  testthat::expect_true(diff(range(get_predVar(zut1)-get_predVar(zut2)[outer((1:16),c(16,0),"+")]))<1e-10)
  testthat::expect_true(diff(range( get_predVar(zut1, newdata=zut1$data)-get_predVar(zut1)))<1e-14) 
  testthat::expect_true(diff(range( get_predVar(zut1, newdata=zut1$data,variances=list(cov=TRUE)) - 
                                      get_predVar(zut1,variances=list(cov=TRUE))))<1e-14) 
  testthat::expect_true(diff(range( get_predVar(zut2, newdata=zut2$data)-get_predVar(zut2)))<1e-14) 
}


{
  meth <- "ML(1,1)" # REML seems to work but the tests are OK only with lower accuracy
  # independent-fit test
  (fb <- fitme(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial(), data=climv, method=meth))
  (fp <- fitme(formula=np2~treatment+(1|clinic2),family=poisson(), data=climv, method=meth))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                                 mod2=list(formula=np2~treatment+(1|clinic2),family=poisson())), 
                         data=climv, method=meth))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=np2~treatment+(1|clinic2),family=poisson()),
                                 mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial())), 
                         data=climv, method=meth))
  
  # test simple init
  (zut3 <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                                 mod2=list(formula=np2~treatment+(1|clinic2),family=poisson())), 
                         init=list(lambda=c("1"=1.1,"2"=2.2)), data=climv, method=meth)) 
  testthat::expect_true(identical(attr(zut3,"optimInfo")$LUarglist$canon.init, list(lambda=c("1"=1.1,"2"=2.2)))) # inits heeded
  # test fancy names in init
  zut4 <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                                 mod2=list(formula=np2~treatment+(1|clinic2),family=poisson())), 
                         init=list(lambda=c("clinic2"=2.2,"clinic"=1.1)), data=climv, method=meth) 
  testthat::expect_true(identical(attr(zut4,"optimInfo")$LUarglist$canon.init, list(lambda=c("1"=1.1,"2"=2.2)))) # inits heeded 
  
  testthat::expect_true(diff(range(logLik(zut1),logLik(zut2), logLik(zut3), logLik(zut4), logLik(fb)+logLik(fp)))<1e-10) # 1e-5 for REML
  pfb <- get_predVar(fb)
  pfp <- get_predVar(fp)
  pzut1 <- get_predVar(zut1)
  pzut1n <- get_predVar(zut1, newdata = zut1$data)
  testthat::expect_true(diff(range(c(pfb,pfp)-pzut1, pzut1n-pzut1))<1e-5) 
  pzut2 <- get_predVar(zut2)
  pzut2n <- get_predVar(zut2, newdata = zut2$data)
  testthat::expect_true(diff(range(c(pfp,pfb)-pzut2, pzut2n-pzut2))<1e-5) 
  
  cat(crayon::yellow("deliberate warnings; ")) # (Deliberately generating warnings:)
  oldopt <- options(warn=0L)
  (zut5 <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial(), init=list(lambda=1.1)),
                                 mod2=list(formula=np2~treatment+(1|clinic2),family=poisson(), init=list(lambda=2.2))), 
                         data=climv, method=meth))
  options(oldopt)
  testthat::expect_true(identical(attr(zut5,"optimInfo")$LUarglist$canon.init, NULL)) # inits NOT heeded as init is not preprocessed

  cat(crayon::yellow("rand family and many post-fit fns; "))# rand.family independent-fit test
  (zut2 <- fitmv(submodels=list(mod2=list(formula=np2~treatment+(1|clinic2),family=poisson(), rand.family=Gamma(log)), 
                                 mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial())), 
                         data=climv))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                                 mod2=list(formula=np2~treatment+(1|clinic2),family=poisson(), rand.family=Gamma(log))), 
                         data=climv))
  (fb <- fitme(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial(), data=climv))
  (fg <- fitme(formula=np2~treatment+(1|clinic2),family=poisson(), data=climv, rand.family=Gamma(log)))
  logLik(zut)-logLik(fb)-logLik(fg)
  testthat::expect_true(diff(range(logLik(zut1), logLik(fb)+logLik(fg)))<1e-9)
  
  cat(crayon::yellow("some extractors; ")) 
  testthat::expect_true(diff(range(residuals(zut1)-c(residuals(fb),residuals(fg))))<1e-5)
  formula(zut1) # list
  terms(zut1) # list 
  nobs(zut1) # single value # decision made

  # rand.family permutation test   (rand.family need to be specified in each submodel)
  (zut1 <- fitmv(submodels=list(mod1=list(cbind(npos,nneg)~treatment+(1|clinic),family=binomial(), rand.family=Gamma(log)),
                                 mod2=list(np2~treatment+(1|clinic),family=poisson(), rand.family=Gamma(log))), 
                         data=climv))
  (zut2 <- fitmv(submodels=list(mod2=list(formula=np2~treatment+(1|clinic),family=poisson(), rand.family=Gamma(log)), 
                                mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial(), rand.family=Gamma(log))), 
                 data=climv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-9)
  update_formulas(zut1,formula(zut1)) 
  update(zut1, formula.=formula(zut1))
  update(zut1, formula.=list(cbind(npos,nneg)~1+(1|clinic),
                             np2~1+(1|clinic)))
  
  # step(zut1) dos not stop() but does not steps...
  
  simulate(zut1,nsim=3) # checks that mv simulate Tpoisson works (possibly depending on same attributes as Tnegbin)
  update_resp(zut1, newresp=simulate(zut1))
  # CI for the variance of the random effect:          
  ( ci <- confint(zut1,parm=function(fit){VarCorr(fit)[1,"Variance"]}, 
                  boot_args=list(nb_cores=9, nsim=9, seed=123)) )
  # The distribution of bootstrap replicates:
  plot(ecdf(ci$call$t))
  zutnull <- update(zut1, submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial(), rand.family=Gamma(log)),
                                         mod2=list(formula=np2~1+(1|clinic),family=poisson(), rand.family=Gamma(log))))
  anova(zut1,zutnull)
  anova(zut1,zutnull, boot.repl=99, nb_cores=9)
  anova(zut1,zutnull, boot.repl=99, nb_cores=9)
  # MSFDR(zutnull,zut1) # documented problem
  
  
  cat(crayon::yellow("Tpoisson; "))# Tpoisson independent-fit test 
  (zut <- fitmv(submodels=list(mod2=list(formula=I(1L+np2)~treatment+(1|clinic2),family=Tpoisson()), 
                                 mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial())), 
                         data=climv))
  (zut <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                                 mod2=list(formula=I(1L+np2)~treatment+(1|clinic2),family=Tpoisson())), 
                         data=climv))
  (fb <- fitme(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial(), data=climv))
  (fTp <- fitme(formula=I(1L+np2)~treatment+(1|clinic2),family=Tpoisson(), data=climv))
  testthat::expect_true(diff(range(logLik(zut), logLik(fb)+logLik(fTp)))<1e-10)
  
  ## negbin(): outer-optimized dispersion parameters. 
  meth <- "ML(1,1)" 
  { # quite slow 2*40s
    # independent-fit test without fixed effects
    (zut1 <- fitmv(submodels=list(mod2=list(formula=I(20*np2)~0+(1|clinic2),family=negbin()), 
                                   mod1=list(formula=cbind(npos,nneg)~0+(1|clinic),family=binomial())), 
                           method=meth,
                           data=climv))
    (zut2 <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~0+(1|clinic),family=binomial()),
                                   mod2=list(formula=I(20*np2)~0+(1|clinic2),family=negbin())), 
                           method=meth,
                           data=climv))
    (fb <- fitme(formula=cbind(npos,nneg)~0+(1|clinic),family=binomial(),  
                 method=meth,
                 data=climv))
    (fn <- fitme(formula=I(20*np2)~0+(1|clinic2),family=negbin(),  
                 method=meth,
                 data=climv))
    testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(fb)+logLik(fn)))<1e-05)
  }
  
  cat(crayon::yellow("Tnegbin; "))## independent-fit test Tnegbin; outer-optimized dispersion parameters
  (zut1 <- fitmv(submodels=list(mod2=list(formula=I(1L+20*np2)~treatment+(1|clinic2),family=Tnegbin()), 
                                 mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial())), 
                         data=climv))
  (zut2 <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                                 mod2=list(formula=I(1L+20*np2)~treatment+(1|clinic2),family=Tnegbin())), 
                         data=climv))
  (fb <- fitme(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial(), data=climv))
  (fTn <- fitme(formula=I(1L+20*np2)~treatment+(1|clinic2),family=Tnegbin(), data=climv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(fb)+logLik(fTn)))<6e-04) 

    
  ## permutation test Tnegbin; outer-optimized dispersion parameters
  (zut1 <- fitmv(submodels=list(mod2=list(formula=I(1L+20*np2)~treatment+(1|clinic),family=Tnegbin()), 
                                 mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial())), 
                         data=climv))
  (zut2 <- fitmv(submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                                 mod2=list(formula=I(1L+20*np2)~treatment+(1|clinic),family=Tnegbin())), 
                         data=climv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<2e-5)
  simulate(zut1,nsim=3) # checks that mv simulate Tnegbin keeps required attributes 
  
  { 
    cat(crayon::yellow("COMPoisson; "))## independent-fit test; outer-optimized dispersion parameters
    data("freight") ## example from Sellers & Shmueli, Ann. Appl. Stat. 4: 943–961 (2010)
    (mod1 <- fitme(broken ~ transfers+(1|id), data=freight, family = COMPoisson()))
    freimv <- freight
    set.seed(123)
    freimv$brok2 <- simulate(mod1)
    (mod2 <- fitme(brok2 ~ transfers,family=poisson(), data=freimv))
    (zut1 <- fitmv(submodels=list(mod1=list(broken ~ transfers+(1|id),family=COMPoisson()), 
                                  mod2=list(brok2 ~ transfers,family=poisson())), 
                   data=freimv))
    (zut2 <- fitmv(submodels=list(mod2=list(brok2 ~ transfers,family=poisson()),
                                  mod1=list(broken ~ transfers+(1|id),family=COMPoisson())), 
                   data=freimv))
    testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(mod1)+logLik(mod2)))<1e-06) 
    
    
    ## permutation test COMPoisson; outer-optimized dispersion parameters
    (zut1 <- fitmv(submodels=list(mod1=list(broken ~ transfers+(1|id),family=COMPoisson()), 
                                  mod2=list(brok2 ~ transfers+(1|id),family=poisson())), 
                   data=freimv))
    (zut2 <- fitmv(submodels=list(mod2=list(brok2 ~ transfers+(1|id),family=poisson()),
                                  mod1=list(broken ~ transfers+(1|id),family=COMPoisson())), 
                   data=freimv))
    testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-4)
    simulate(zut1,nsim=3) # checks that mv simulate COMPoisson 
    
    (zutxx <- fitmv(submodels=list(mod2=list(brok2 ~ transfers+(1|id),family=COMPoisson()),
                                   mod1=list(broken ~ transfers+(1|id),family=COMPoisson())), 
                    data=freimv, lower=list(COMP_nu=c("1"=0.5,"2"=1.1)))) # names needed !
    attr(zutxx,"optimInfo")$LUarglist$init.optim$COMP_nu
    (zutxx <- fitmv(submodels=list(mod2=list(brok2 ~ transfers+(1|id),family=COMPoisson()),
                                   mod1=list(broken ~ transfers+(1|id),family=COMPoisson())), 
                    data=freimv, lower=list(COMP_nu=c("2"=1.1)))) 
    attr(zutxx,"optimInfo")$LUarglist$init.optim$COMP_nu
    (zutxx <- fitmv(submodels=list(mod2=list(brok2 ~ transfers+(1|id),family=poisson()),
                                   mod1=list(broken ~ transfers+(1|id),family=COMPoisson())), 
                    data=freimv, lower=list(COMP_nu=c("2"=1.1)))) 
    attr(zutxx,"optimInfo")$LUarglist$init.optim$COMP_nu
    (zutxx <- fitmv(submodels=list(mod2=list(brok2 ~ transfers+(1|id),family=poisson()),
                                   mod1=list(broken ~ transfers+(1|id),family=COMPoisson(1.1))), 
                    data=freimv)) 
    (zutxx <- fitmv(submodels=list(mod2=list(brok2 ~ transfers+(1|id),family=poisson()),
                                   mod1=list(broken ~ transfers+(1|id),family=COMPoisson(1.1))), 
                    data=freimv, lower=list(COMP_nu=c("2"=1.1)))) 
    (zutxx <- fitmv(submodels=list(mod2=list(brok2 ~ transfers+(1|id),family=COMPoisson(1.1)),
                                   mod1=list(broken ~ transfers+(1|id),family=COMPoisson())), 
                    data=freimv, lower=list(COMP_nu=c("2"=1.1))))  
    ## assign("last.warning", NULL, envir = baseenv()) # flush the .COMP_maxn() warnings()
  }  
}

data("Loaloa")
lll <- Loaloa
lll$ID <- lll$ID2 <- 1+(seq(nrow(lll)) %% 2)
lll$resp <- 1+lll$npos + 10*lll$ID
(tnb <- fitme(resp~1+(1|ID), data=lll,family=Tnegbin()))
set.seed(124)
lll$r1 <- simulate(tnb)
lll$r2 <- simulate(tnb)
lll$long2 <- lll$longitude 

## independent-fit test 2 Tnegbin
(tnb1 <- fitme(r1~1+(1|ID), data=lll,family=Tnegbin()))
(tnb2 <- fitme(r2~1+(1|ID), data=lll,family=Tnegbin()))
(zut1 <- fitmv(submodels=list(mod1=list(r1~1+(1|ID),family=Tnegbin()),
                               mod2=list(r2~1+(1|ID2),family=Tnegbin())), 
                       data=lll))
(zut2 <- fitmv(submodels=list(mod2=list(r2~1+(1|ID2),family=Tnegbin()),
                               mod1=list(r1~1+(1|ID),family=Tnegbin())), 
                       data=lll))
testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(tnb1)+logLik(tnb2)))<1e-06)

# permutation test
(zut1 <- fitmv(submodels=list(mod1=list(r1~1+(1|ID),family=Tnegbin()),
                               mod2=list(r2~1+(1|ID),family=Tnegbin())), 
                       data=lll))
(zut2 <- fitmv(submodels=list(mod2=list(r2~1+(1|ID),family=Tnegbin()),
                               mod1=list(r1~1+(1|ID),family=Tnegbin())), 
                       data=lll))
testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-8)

if (FALSE) { # slowish ~4s 
  cat(crayon::yellow("Matern; "))## independent-fit test with a Matern
  (tp1 <- fitme(r1~1+Matern(1|longitude+latitude), data=lll,family=poisson()))
  (tp2 <- fitme(r2~1+(1|ID), data=lll,family=poisson()))
  (zut1 <- fitmv(submodels=list(mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson()),
                                  mod2=list(r2~1+(1|ID2),family=poisson())), 
                          data=lll, verbose=c(TRACE=TRUE)))
  (zut2 <- fitmv(submodels=list(mod2=list(r2~1+(1|ID2),family=poisson()),
                                  mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson())), 
                          data=lll, verbose=c(TRACE=TRUE)))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(tp1)+logLik(tp2)))<1e-08)
}

if (FALSE) { # slow! > 4mn total
  ## independent-fit test with two Matern
  (tp1 <- fitme(r1~1+Matern(1|longitude+latitude), data=lll,family=poisson()))
  (tp2 <- fitme(r2~1+Matern(1|long2+latitude), data=lll,family=poisson()))
  (zut1 <- fitmv(submodels=list(mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson()),
                                 mod2=list(r2~1+Matern(1|long2+latitude),family=poisson())), 
                         data=lll, verbose=c(TRACE=TRUE)))
  (zut2 <- fitmv(submodels=list(mod2=list(r2~1+Matern(1|long2+latitude),family=poisson()),
                                 mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson())), 
                         data=lll, verbose=c(TRACE=TRUE)))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(tp1)+logLik(tp2)))<1e-08)
}

if (FALSE) { # slow! > 150s total
  ## independent-fit test with two Matern and fixed nu's
  (tp1 <- fitme(r1~1+Matern(1|longitude+latitude), data=lll,family=poisson(), fixed=list(nu=0.5)))
  (tp2 <- fitme(r2~1+Matern(1|long2+latitude), data=lll,family=poisson(), fixed=list(nu=1)))
  (zut1 <- fitmv(submodels=list(mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson(), fixed=list(nu=0.5)),
                                 mod2=list(r2~1+Matern(1|long2+latitude),family=poisson(), fixed=list(nu=1))), 
                         data=lll, verbose=c(TRACE=TRUE)))
  (zut2 <- fitmv(submodels=list(mod2=list(r2~1+Matern(1|long2+latitude),family=poisson(), fixed=list(nu=1)),
                                 mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson(), fixed=list(nu=0.5))), 
                         data=lll, verbose=c(TRACE=TRUE)))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(tp1)+logLik(tp2)))<1e-08)
}

if (FALSE) { # slow 22s
  # permutation test
  (zut1 <- fitmv(submodels=list(mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson(), fixed=list(nu=1)),
                                 mod2=list(r2~1+Matern(1|longitude+latitude),family=poisson())), 
                         data=lll, verbose=c(TRACE=1)))
  (zut2 <- fitmv(submodels=list(mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson()),
                                 mod2=list(r2~1+Matern(1|longitude+latitude),family=poisson(), fixed=list(nu=1))), 
                         data=lll, verbose=c(TRACE=1)))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-08)
  predict(zut1)
  simulate(zut1)
}

# default-name test 
cat(crayon::yellow("'fixed' both in submodel and global call; "))
(zut1 <- fitmv(submodels=list(mod2=list(r2~1+Matern(1|longitude+latitude),family=poisson()),
                               mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson(), fixed=list(rho=1,nu=1))), 
                       data=lll, fixed=list(lambda=c("1"=666)), verbose=c(TRACE=1)))

(zut2 <- fitmv(submodels=list(mod1=list(r1~1+Matern(1|longitude+latitude),family=poisson(), fixed=list(rho=1,nu=1)),
                       mod2=list(r2~1+Matern(1|longitude+latitude),family=poisson())), 
                       data=lll, fixed=list(lambda=c(666)), verbose=c(TRACE=1)))
testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-08)

{
  data("blackcap")
  MLdistMat2 <- as.matrix(proxy::dist(blackcap[,c("latitude","longitude")]))
  MLcorMat2 <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                         nu=0.6285603,rho=0.0544659)
  cap_mv <- blackcap
  cap_mv$name <- as.factor(rownames(blackcap))                
  cap_mv$grp <- 1L+(blackcap$migStatus>1)   
  set.seed(123)
  cap_mv$status2 <- blackcap$migStatus+ rnorm(14,sd=0.001)
  
  cat(crayon::yellow("corrMatrix vs Matern; (map_ranef too)"))# corrMatrix
  
  ## independent-fit test with compar to equivalent Matern:
  # need to fix this phi to avoid logLik uncertainty for low phi:
  (mod1 <- fitme(migStatus ~ 1+ corrMatrix(1|name),data=cap_mv, corrMatrix=MLcorMat2, fixed=list(phi=0.1)))
  (mod2 <- fitme(status2 ~ 1+ (1|grp),data=cap_mv))
  (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ corrMatrix(1|name), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ (1|grp))), 
                 data=cap_mv, corrMatrix=MLcorMat2))
  (zut1b <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ corrMatrix(1|name), fixed=list(phi=0.1), corrMatrix=MLcorMat2),
                          mod2=list(status2 ~ 1+ (1|grp))), 
                  data=cap_mv))
  (zut1c <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ corrMatrix(1|name), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ (1|grp))), 
                 data=cap_mv, covStruct=list(corrMatrix=MLcorMat2,NULL,corrMatrix=NULL))) # a bit of abuse...
  (zut2 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ (1|grp)),
                         mod1=list(migStatus ~ 1+ corrMatrix(1|name),fixed=list(phi=0.1))), 
                 data=cap_mv, corrMatrix=MLcorMat2))
  (zut3 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ Matern(1|longitude+latitude), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ (1|grp))), 
                 data=cap_mv, fixed=list(rho=0.0544659,nu=0.6285603)))
  (zut4 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ Matern(1|longitude+latitude), fixed=list(phi=0.1,rho=0.0544659,nu=0.6285603)),
                         mod2=list(status2 ~ 1+ (1|grp))), 
                 data=cap_mv))
  map_ranef(zut4)
  get_predVar(zut4, variances=list(cov=TRUE)) 
  if (FALSE) {
    get_predVar(zut4, variances=list(cov=TRUE)) -> bla
    bli <- bla
    dim(bli) <- c(14,2,14,2)
    bli[,1,,2]
    bli[,1,,1]
    bli[,2,,2]
  }
  
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut1b), logLik(zut1c), 
                                   logLik(zut2), logLik(zut3), logLik(zut4), logLik(mod1)+logLik(mod2)))<1e-10)
  testthat::expect_true(diff(range( predict(zut1, newdata=zut1$data)-predict(zut1)))<1e-14)
  
  ## permutation test with compar to equivalent Matern:
  (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+corrMatrix(1|name), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ corrMatrix(1|name))), 
                 data=cap_mv, corrMatrix=MLcorMat2))
  (zut2 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ corrMatrix(1|name)),
                         mod1=list(migStatus ~ 1+corrMatrix(1|name), fixed=list(phi=0.1))), 
                 data=cap_mv, corrMatrix=MLcorMat2))
  (zut3 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ Matern(1|longitude+latitude), fixed=list(phi=0.1,rho=0.0544659,nu=0.6285603)),
                         mod2=list(status2 ~ 1+ Matern(1|longitude+latitude))), 
                 data=cap_mv))
  (zut4 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ Matern(1|longitude+latitude)),
                         mod1=list(migStatus ~ 1+ Matern(1|longitude+latitude),fixed=list(phi=0.1,rho=0.0544659,nu=0.6285603))), 
                 data=cap_mv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(zut3), logLik(zut4)))<1e-9)
  testthat::expect_true(diff(range( predict(zut1, newdata=zut1$data)-predict(zut1)))<1e-14)
  get_predVar(zut4, variances=list(cov=TRUE)) 
  if (FALSE) { # may be used to check the mapping, but to have asym off-diag blocks, use migStatus ~ 1+ Matern(1|longitude+latitude)+(1|grp), the latter with fixed large variance
    get_predVar(zut4, variances=list(cov=TRUE)) -> bli
    nresp <- length(attr(bli, "respnames"))
    nobs <- nrow(bli)/nresp
    dim(bli) <- c(nobs,nresp,nobs,nresp)
    dimnames(bli)[[2]] <- dimnames(bli)[[4]] <- attr(bli, "respnames")
    bli[,1,,2]
    bli[,1,,1]
    bli[,2,,2]
    bli[1,,1,]
  }
  
  cat(crayon::yellow("distMatrix; "))
  (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+Matern(1|name), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ Matern(1|name))),
                 distMatrix=MLdistMat2, fixed=list(rho=0.0544659,nu=0.6285603), 
                 data=cap_mv))
  (zut2 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+Matern(1|name), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ Matern(1|name))),
                 distMatrix=2*MLdistMat2, fixed=list(rho=0.0544659/2,nu=0.6285603), 
                 data=cap_mv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-9)
  
  if (FALSE) { # Amusingly one can fit the same data by two sub-models (even tow diferent ones). 
    # we have twice the conditional likelihood, and once the ranef lik, so (as in any mv model with shared ranefs) marginal likelihoods do not add up.
    # The phi's are low and they lower value must be controlled in order to make detailed numerical comparisons
    (mod1 <- fitme(migStatus ~ 1+ corrMatrix(1|name),data=cap_mv, corrMatrix=MLcorMat))
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ corrMatrix(1|name)),
                           mod2=list(migStatus ~ 1+ corrMatrix(1|name))), 
                   data=cap_mv, corrMatrix=MLcorMat2))
    vcov(zut1) # with phi->0 and a common ranef, this may be logical
  } 
  
  
  cat(crayon::yellow("corrMatrix vs Cauchy; ")) 
  MLcorMat3 <- CauchyCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                         shape=1,longdep=0.5) # and default rho=1!
  
  ## independent-fit test with compar to equivalent Matern:
  # need to fix this phi to avoid logLik uncertainty for low phi:
  (mod1 <- fitme(migStatus ~ 1+ corrMatrix(1|name),data=cap_mv, corrMatrix=MLcorMat3, fixed=list(phi=0.1)))
  (mod2 <- fitme(status2 ~ 1+ (1|grp),data=cap_mv))
  (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ corrMatrix(1|name), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ (1|grp))), 
                 data=cap_mv, corrMatrix=MLcorMat3))
  (zut1b <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ corrMatrix(1|name), fixed=list(phi=0.1), corrMatrix=MLcorMat3),
                          mod2=list(status2 ~ 1+ (1|grp))), 
                  data=cap_mv))
  (zut1c <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ corrMatrix(1|name), fixed=list(phi=0.1)),
                          mod2=list(status2 ~ 1+ (1|grp))), 
                  data=cap_mv, covStruct=list(corrMatrix=MLcorMat3,NULL,corrMatrix=NULL))) # a bit of abuse...
  (zut2 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ (1|grp)),
                         mod1=list(migStatus ~ 1+ corrMatrix(1|name),fixed=list(phi=0.1))), 
                 data=cap_mv, corrMatrix=MLcorMat3))
  (zut3 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ Cauchy(1|longitude+latitude), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ (1|grp))), 
                 data=cap_mv, fixed=list(shape=1,rho=1,longdep=0.5)))
  (zut4 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ Cauchy(1|longitude+latitude), fixed=list(phi=0.1,rho=1,shape=1,longdep=0.5)),
                         mod2=list(status2 ~ 1+ (1|grp))), 
                 data=cap_mv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut1b), logLik(zut1c), 
                                   logLik(zut2), logLik(zut3), logLik(zut4), logLik(mod1)+logLik(mod2)))<1e-10)
  testthat::expect_true(diff(range( predict(zut1, newdata=zut1$data)-predict(zut1)))<1e-14)
  
  ## permutation test with compar to equivalent Cauchy:
  (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+corrMatrix(1|name), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ corrMatrix(1|name))), 
                 data=cap_mv, corrMatrix=MLcorMat3))
  (zut2 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ corrMatrix(1|name)),
                         mod1=list(migStatus ~ 1+corrMatrix(1|name), fixed=list(phi=0.1))), 
                 data=cap_mv, corrMatrix=MLcorMat3))
  (zut3 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+ Cauchy(1|longitude+latitude), fixed=list(phi=0.1,rho=1,shape=1,longdep=0.5)),
                         mod2=list(status2 ~ 1+ Cauchy(1|longitude+latitude))), 
                 data=cap_mv))
  (zut4 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ Cauchy(1|longitude+latitude)),
                         mod1=list(migStatus ~ 1+ Cauchy(1|longitude+latitude),fixed=list(phi=0.1,rho=1,shape=1,longdep=0.5))), 
                 data=cap_mv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(zut3), logLik(zut4)))<1e-9)
  testthat::expect_true(diff(range( predict(zut1, newdata=zut1$data)-predict(zut1)))<1e-14)
  
  cat(crayon::yellow("distMatrix with Cauchy; "))
  (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+Cauchy(1|name), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ Cauchy(1|name))),
                 distMatrix=MLdistMat2, fixed=list(rho=1,shape=1,longdep=0.5), 
                 data=cap_mv))
  (zut2 <- fitmv(submodels=list(mod1=list(migStatus ~ 1+Cauchy(1|name), fixed=list(phi=0.1)),
                         mod2=list(status2 ~ 1+ Cauchy(1|name))),
                 distMatrix=2*MLdistMat2, fixed=list(rho=1/2,shape=1,longdep=0.5), 
                 data=cap_mv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-9)
}

{cat(crayon::yellow("IMRF; "))# fit IMRF 
  { # create IMRF model
    spd <- sp::SpatialPointsDataFrame(coords = blackcap[, c("longitude", "latitude")],
                                      data = blackcap)
    ## Creating the mesh 
    mesh <- INLA::inla.mesh.2d(loc = INLA::inla.mesh.map(sp::coordinates(spd)), 
                               cutoff=30,
                               max.edge = c(3, 20)) 
    mesh$n ## 40
    matern <- INLA::inla.spde2.matern(mesh)
  }
  
  { # Aim was first to check mv IMRF predictions, but putting a (G)LM as first model caught many problems
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1),
                           mod2=list(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern))), 
                   fixed=list(phi=c(0.02,0.02)),
                   data=cap_mv)) 
    p1 <- predict(zut1) 
    p2 <- predict(zut1, newdata=zut1$data) # note difference in frame attribute. What do we want? _FIXME_
    testthat::expect_true(diff(range(p2-p1))<1e-14)
  }
  
  
  {   # independent-fit test
    # a bit tricky bc nloptr does not find the univariate optimum without some help (... tiny blackcap data)
    # and in mv case, it finds it by default or not depending on order...
    (mod1 <- fitme(migStatus ~ means + IMRF(1|longitude+latitude, model=matern), init=list(lambda=0.02),
                   data=cap_mv, verbose=c(TRACE=FALSE)))
    (mod2 <- fitme(status2 ~ 1+ Matern(1|longitude+latitude), fixed=list(phi=0.1,rho=0.0544659,nu=0.6285603), data=cap_mv,
                   verbose=c(TRACE=FALSE)))
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ means + IMRF(1|longitude+latitude, model=matern)),
                           mod2=list(status2 ~ 1+ Matern(1|longitude+latitude), fixed=list(phi=0.1))), 
                   data=cap_mv, fixed=list(rho=0.0544659,nu=0.6285603), init=list(lambda=c(0.02,NA)))) # Matern in spprec...
    (zut2 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ Matern(1|longitude+latitude), fixed=list(phi=0.1)),
                           mod1=list(migStatus ~ means + IMRF(1|longitude+latitude, model=matern))), 
                   data=cap_mv, fixed=list(rho=0.0544659,nu=0.6285603), init=list(lambda=c(NA,0.02)))) # Matern in spprec...
    testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(mod1)+logLik(mod2)))<1e-08)
  }
  {   # permutation test
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ means + IMRF(1|longitude+latitude, model=matern)),
                                  mod2=list(status2 ~ means+ IMRF(1|longitude+latitude, model=matern))), 
                   data=cap_mv)) 
    (zut2 <- fitmv(submodels=list(mod2=list(status2 ~ means+ IMRF(1|longitude+latitude, model=matern)),
                                          mod1=list(migStatus ~ means + IMRF(1|longitude+latitude, model=matern))), 
                           data=cap_mv)) # sensitive to .solve_crossr22( ., use_crossr22=TRUE)
    testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-07)
    
    # The following with status2 ~ 1+ ... rather than status2 ~ means+ ... makes lambda diverge (jointly with kappa: does not occur for fixed kappa). 
    # Then all comparisons are sensitive to floating point precision (and no stringent permutation test is passed).
    # but let us keep the zut2_testr22 computation as it has been useful in the past to 'torture-test' .calc_r22() algorithms
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ means + IMRF(1|longitude+latitude, model=matern)),
                           mod2=list(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern))), 
                   init=list(lambda=c(0.02),phi=c(0.02,0.02)),
                   data=cap_mv)) 
    crit <- diff(range(logLik(zut1), -6.4882739678))
    testthat::test_that(paste0("criterion was ",signif(crit,6)," from -6.4882739678"), # affected by use_ZA_L or .calc_r22()
                        testthat::expect_true(crit<1e-08))
    (zut2_testr22 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern)),
                                          mod1=list(migStatus ~ means + IMRF(1|longitude+latitude, model=matern))), 
                           init=list(lambda=c(0.02),phi=c(0.02,0.02)),
                           data=cap_mv)) # has been sensitive to .solve_crossr22( ., use_crossr22=TRUE)
    crit <- diff(range(-6.48830317593 , logLik(zut2_testr22)))
    testthat::test_that(paste0("criterion was ",signif(crit,6)," from -6.48830317593"), # affected by use_ZA_L or .calc_r22()
                        testthat::expect_true(crit<1e-08))
  }
  if (spaMM.getOption("example_maxtime")>119) { cat(crayon::yellow("multIMRF indep-fit tests; "))
    (mrf1fixx <- fitme(migStatus ~ 1 + (1|pos) + 
                    multIMRF(1|longitude+latitude,margin=5,levels=2), 
                  data=blackcap, fixed=list(phi=1,lambda=c("1"=0.5),
                                            hyper=list("1"=list(hy_kap=0.1,hy_lam=1)))) )
    (mrf2 <- fitme(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern), 
                   data=cap_mv)) 
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1 + (1|pos) + multIMRF(1|longitude+latitude,margin=5,levels=2)), 
                                  mod2=list(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern))), 
                   fixed=list(phi=1,lambda=c("1"=0.5), 
                              hyper=list("1"=list(hy_kap=0.1,hy_lam=1))), 
                   data=cap_mv))
    testthat::expect_true(diff(range(logLik(zut1), logLik(mrf1fixx)+logLik(mrf2)))<1e-10)
    #
    (mrf1fix <- fitme(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2), 
                      data=blackcap, fixed=list(hyper=list("1"=list(hy_kap=0.1,hy_lam=1)))) )
    (mrf2 <- fitme(status2 ~ 1, data=cap_mv)) 
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2)), 
                                  mod2=list(status2 ~ 1)), 
                   fixed=list(hyper=list("1"=list(hy_kap=0.1,hy_lam=1))), 
                   data=cap_mv))
    testthat::expect_true(diff(range(logLik(zut1), logLik(mrf1fix)+logLik(mrf2)))<1e-10)

    (mrf1fix <- fitme(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2), 
                      data=blackcap, fixed=list(hyper=list("1"=list(hy_kap=0.1)))) )
    (mrf2 <- fitme(status2 ~ 1, data=cap_mv)) 
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2)), 
                                  mod2=list(status2 ~ 1)), 
                   fixed=list(hyper=list("1"=list(hy_kap=0.1))), 
                   data=cap_mv))
    testthat::expect_true(diff(range(logLik(zut1), logLik(mrf1fix)+logLik(mrf2)))<1e-10)
    
    (mrf1fix <- fitme(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2), 
                      data=blackcap, fixed=list(hyper=list("1"=list(hy_lam=1)))) )
    (mrf2 <- fitme(status2 ~ 1, data=cap_mv)) 
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2)), 
                                  mod2=list(status2 ~ 1)), 
                   fixed=list(hyper=list("1"=list(hy_lam=1))),
                   data=cap_mv))
    testthat::expect_true(diff(range(logLik(zut1), logLik(mrf1fix)+logLik(mrf2)))<1e-10)
    
    {  # particularly slow
      {
        # now there is a  global maximum at phi=0 which may be missed when phi outer optim is used,
        # But not missed when default augZXy is used. augZXy is used by default for mrf1, but not in mv  => discrepancy
        spaMM.options(allow_augZXy=FALSE)
        logLik(mrf1F <- fitme(migStatus ~ 1 +multIMRF(1|longitude+latitude,margin=5,levels=2),data=cap_mv)) # bad
        logLik(fitme(migStatus ~ 1 +multIMRF(1|longitude+latitude,margin=5,levels=2),data=cap_mv, init=list(phi=(1e-4)))) # good
        spaMM.options(allow_augZXy=NULL)  # 
        # This may lose accuracy when .calc_r22() is modified! :
        logLik(mrf1T <- fitme(migStatus ~ 1 +multIMRF(1|longitude+latitude,margin=5,levels=2),data=cap_mv)) # good; use y-augmented matrix 
        logLik(mrf1I <- fitme(migStatus ~ 1 +multIMRF(1|longitude+latitude,margin=5,levels=2),data=cap_mv, init.HLfit=list(phi=1))) # good; use y-augmented matrix !
        # so to compare with mv fit we must provide a low init phi to the submodel  
      }
      
      (mrf2 <- fitme(status2 ~ 1, data=cap_mv)) 
      (zut0 <- fitmv(submodels=list(mod1=list(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2)), 
                                    mod2=list(status2 ~ 1)), 
                     data=cap_mv, init=list(phi=list("1"=1e-4))))
      testthat::expect_true(diff(range(logLik(zut0), logLik(mrf1T)+logLik(mrf2)))<3e-6)
      
      (mrf3 <- fitme(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern), data=cap_mv)) 
      # inits to try to speed the fit () (and lambda=c("3"=.) works as it should)
      (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2)), 
                                    mod2=list(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern))), 
                     data=cap_mv, init=list(lambda=c("3"=5),phi=list("1"=1e-4,"2"=0.05)),verbose=c(TRACE=TRUE)))
      (zut2 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern)),
                                    mod1=list(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2))), 
                     data=cap_mv, init=list(lambda=c("1"=5),phi=list("1"=0.05,"2"=1e-4)))) 
      testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(mrf1T)+logLik(mrf3)))<3e-6)
    }
    
  } else cat(crayon::bgGreen("\n multIMRF indep-fit tests are slow (~119s). Run them once in a while; "))
  
  { cat(crayon::yellow("multIMRF permutation tests; ")) # ~ 13s
    # reason for init phi as above: to avoid a local maximum
    (zut1 <- fitmv(submodels=list(mod1=list(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2)), 
                                  mod2=list(status2 ~ 1+ multIMRF(1|longitude+latitude,margin=5,levels=2))), 
                   data=cap_mv, init=list(phi=list("1"=1e-4,"2"=1e-4))))
    (zut2 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ multIMRF(1|longitude+latitude,margin=5,levels=2)),
                                  mod1=list(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2))), 
                   data=cap_mv, init=list(phi=list("1"=1e-4,"2"=1e-4)))) 
    testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-08)
    
    if (FALSE) {
      # The following code fits a single (repeated) kappa value, and only one lambda hyperparam. 
      # So the single IMRF from multIMRF(. ,levels=1) is possibly
      # recognized as first level of multIMRF(. ,levels=2). Which may be nice, or may not be the intent. (__FIXME__).
      # (What does .calc_normalized_ZAlist() do ?)
      # Adding a fictitious argument bla=666 has no effect as it is not retained in the expanded formula.
      (zut2 <- fitmv(submodels=list(mod2=list(status2 ~ 1+ multIMRF(1|longitude+latitude,margin=5,levels=1)),
                                    mod1=list(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2))), 
                     data=cap_mv, init=list(phi=list("1"=1e-4,"2"=1e-4)))) 
    }
  }
}

{ cat(crayon::yellow("adjacency; "))
  data("scotlip")
  (mod1 <- fitme(cases ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)),
        adjMatrix=Nmatrix, family=poisson(), data=scotlip) )
  scotmv <- scotlip
  set.seed(123)
  scotmv$cases2 <- simulate(mod1)
  scotmv$code2 <- scotmv$gridcode
  (mod2 <- fitme(cases2 ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)),
                adjMatrix=Nmatrix, family=poisson(), data=scotmv) )
  
  {   # independent-fit test
    # a bit tricky bc nloptr does not find the univariate optimum without some help (... tiny blackcap data)
    # and in mv case, it finds it by default or not depending on order...
    (zut1 <- fitmv(submodels=list(mod1=list(cases ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)), family=poisson()),
                           mod2=list(cases2 ~ I(prop.ag/10)+adjacency(1|code2)+offset(log(expec)), family=poisson())), 
                   data=scotmv,covStruct=list(adjMatrix=Nmatrix,adjMatrix=Nmatrix))) # Matern in spprec...
    (zut2 <- fitmv(submodels=list(mod1=list(cases2 ~ I(prop.ag/10)+adjacency(1|code2)+offset(log(expec)), family=poisson()),
                           mod2=list(cases ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)), family=poisson())), 
                   data=scotmv,covStruct=list(adjMatrix=Nmatrix,adjMatrix=Nmatrix))) # Matern in spprec...
    testthat::expect_true(diff(range(logLik(zut1), logLik(zut2), logLik(mod1)+logLik(mod2)))<1e-08)
  }
  {   # permutation test
    # here to the order affects nloptr... tiny blackcap data again. Note lambda divergence.
    # we can avoid that by either providing fixed or init phi (the latter slow)
    (zut1 <- fitmv(submodels=list(mod1=list(cases ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)), family=poisson()),
                           mod2=list(cases2 ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)), family=poisson())), 
                   data=scotmv,adjMatrix=Nmatrix)) 
    (zut2 <- fitmv(submodels=list(mod2=list(cases2 ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)), family=poisson()),
                           mod1=list(cases ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)), family=poisson())), 
                   data=scotmv,adjMatrix=Nmatrix)) 
    testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-08)
    # 
    (zut1 <- fitmv(submodels=list(mod1=list(cases ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec))),
                           mod2=list(cases2 ~ I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)))), 
                   data=scotmv,adjMatrix=Nmatrix)) 
    confint(zut1, parm="(Intercept)_1")
  }
  
}

{ 
  data("sleepstudy",package = "lme4")
  (mod0 <- fitme(Reaction ~ Days + AR1(1|Days), data = sleepstudy, fixed=list(ARphi=0.5,lambda=1000)))
  # bc the unconstrained fit has essentially lambda=0, which cause later ambiguities
  # Don't give to much importance to the result: the fixed Days capture the trend and the 'correlated' rnaefs are not really correlated
  # plot(ranef(mod0)[[1]])
  sleepmv <- sleepstudy
  set.seed(124)
  sleepmv$Reaction <- simulate(mod0)
  sleepmv$reac2 <- simulate(mod0)
  sleepmv$days2 <- sleepmv$Days
  
  cat(crayon::yellow("non-standard REML (both ways); "))
  (re1 <- fitme(Reaction ~ 1 + (1|Days), REMLformula=~Days, data = sleepmv, method="REML"))
  (re2 <- fitme(reac2 ~ Days + (1|days2), REMLformula=~1, data = sleepmv, method="REML"))
  (rezut1 <- fitmv(submodels=list(mod1=list(Reaction ~ 1 + (1|Days), REMLformula=~Days),
                                  mod2=list(reac2 ~ Days + (1|days2), REMLformula=~1)), 
                   data=sleepmv, method="REML"))
  testthat::expect_true(diff(range(logLik(rezut1), logLik(re1)+logLik(re2)))<1e-08)
  
  cat(crayon::yellow("AR1; "))
  (mod1 <- fitme(Reaction ~ Days + AR1(1|Days), data = sleepmv))
  (mod2 <- fitme(reac2 ~ Days + AR1(1|Days), data = sleepmv))
  (zut1 <- fitmv(submodels=list(mod1=list(Reaction ~ Days + AR1(1|Days)),
                                mod2=list(reac2 ~ Days + AR1(1|days2))), 
                 data=sleepmv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(mod1)+logLik(mod2)))<1e-08)
  
  (zut1 <- fitmv(submodels=list(mod1=list(Reaction ~ Days + AR1(1|Days)),
                                mod2=list(reac2 ~ Days + AR1(1|Days))), 
                 data=sleepmv))
  (zut2 <- fitmv(submodels=list(mod2=list(reac2 ~ Days + AR1(1|Days)), 
                                mod1=list(Reaction ~ Days + AR1(1|Days))),
                 data=sleepmv))
  testthat::expect_true(diff(range(logLik(zut1), logLik(zut2)))<1e-08)
  
}

if(FALSE) { cat(crayon::yellow("simulation study; "))
  
  
  {
    set.seed(123)
    replic <- function(nind=1000L, cor=-0.5, lambda=c(1,0.2), binsize=1L, method=if (binsize<4L) {"PQL/L"} else {"ML"}, return.fit=FALSE, ...) {
      cat(".")
      u1 <- rnorm(nind)
      u2 <- rnorm(nind)
      lam1 <- lambda[1]
      lam2 <- lambda[2]
      L <- chol(matrix(c(lam1,sqrt(lam1*lam2)*cor,sqrt(lam1*lam2)*cor,lam2),ncol=2))
      v <- t(L) %*% rbind(u1,u2)
      eta1 <- v[1,]
      eta2 <- 1+v[2,]
      surv <- rbinom(length(eta1),binsize,inv.logit(eta1))
      feco <- rpois(length(eta2), lambda = exp(eta2))
      lfh <- data.frame(id=seq_along(eta1), surv=surv, feco=feco, binsize=binsize)
      (fitlfh <- fitmv(submodels=list(list(cbind(surv,binsize-surv) ~ 1+(mv(1,2)|id), family=binomial()),
                                      list(feco ~ 1+(mv(1,2)|id), family=poisson())),
                       data=lfh, method=method, ...))
      if (return.fit) {
        return(fitlfh)
      } else {
        corr <- VarCorr(fitlfh)[2,"Corr."]
        cat(corr)
        corr
      }
    }
    replicate(3L, replic(nind=1000L))
    replicate(10L, replic(nind=3000L, binsize=10)) # much better...
    replicate(10L, replic(nind=300L, binsize=10)) # much better...
    PQLLdist <- replicate(50L, replic(nind=300L, binsize=2)) # passable
    PQLdist <- replicate(50L, replic(nind=300L, binsize=2, method="PQL")) # passable
    # one can guess a *small* advantage of REML:
    plot(ecdf(PQLdist))
    plot(ecdf(PQLLdist), add=TRUE, col="red")
    #
    PQLdist <- replicate(10L, replic(nind=1000L, binsize=1, method="PQL")) # binary...
    PQLdist <- replicate(10L, replic(nind=1000L, cor=0, binsize=1, method="PQL")) # binary...
    PQLdist <- replicate(10L, replic(nind=300L, cor=-0.5, binsize=1, method="PQL")) # binary & small sample size -> extreme negative
    PQLdist <- replicate(10L, replic(nind=10000L, cor=-0.5, binsize=1, method="PQL")) # large sample size does not improve this
    
    
    (zut <- replic(nind=300L, cor=-0.5, binsize=1, method="PQL", return.fit=TRUE, verbose=c(TRACE=TRUE)))
    
    
    replic(nind=1000L, binsize=1, return.fit=TRUE)
    
    plot(eta1,eta2)
    cov(t(v))
    plot(eta1,surv)
    plot(eta2,feco)
    plot(eta1,feco)
    cov(cbind(surv,feco))
    
  }
  
  {
    set.seed(123)
    replic2 <- function(nind=1000L, lambda=0.2, binsize=1L, method=if (binsize<4L) {"PQL/L"} else {"ML"}, return.fit=FALSE, ...) {
      cat(".")
      u1 <- sqrt(lambda)*rnorm(nind)
      eta_f <- 1+u1
      feco <- rpois(length(eta_f), lambda = exp(eta_f))
      surv <- rbinom(length(eta_f),binsize,inv.logit(-(log(feco)-1)))
      lfh <- data.frame(id=seq_along(eta_f), id2=seq_along(eta_f), surv=surv, feco=feco, binsize=binsize)
      (fitlfh <- fitmv(submodels=list(list(cbind(surv,binsize-surv) ~ 1+(mv(1,2)|id), family=binomial()),
                                      list(feco ~ 1+(mv(1,2)|id), family=poisson())),
                       data=lfh, method=method, ...))
      if (return.fit) {
        return(fitlfh)
      } else {
        corr <- VarCorr(fitlfh)[2,"Corr."]
        cat(corr)
        corr
      }
    }
    replicate(10L, replic2(nind=1000L))
    replicate(10L, replic2(nind=1000L, binsize=10))
    replicate(10L, replic2(nind=3000L))
    replicate(10L, replic2(nind=3000L, binsize=10)) # 
    
    (zut <- replic2(nind=3000L, binsize=10, return.fit=TRUE))
    with(zut$data, plot(surv,feco))
    
    fitme(cbind(surv,binsize-surv) ~ 1+(1|id), family=binomial(),
          data=zut$data)
    fitmv(submodels=list(list(cbind(surv,binsize-surv) ~ 1+(1|id), family=binomial()),
                         list(feco ~ 1+(1|id2), family=poisson())),
          data=zut$data)
    fitmv(submodels=list(list(cbind(surv,binsize-surv) ~ 1+(mv(1,2)|id), family=binomial()),
                         list(feco ~ 1+(mv(1,2)|id), family=poisson())),
          data=zut$data)
    
    (debu <- fitmv(submodels=list(list(cbind(surv,binsize-surv) ~ 1+(mv(1,2)|id), family=binomial()),
                                  list(feco ~ 1+(mv(1,2)|id2), family=poisson())),
                   data=zut$data[1:20,]))
    
    
    
    plot(eta1,eta2)
    cov(t(v))
    plot(eta1,surv)
    plot(eta2,feco)
    plot(eta1,feco)
    cov(cbind(surv,feco))
    
  }
  
  {
    
    simfun <- function(nind=570L, cor=-0.5, lambda=c(0.2,0.1), return.fit=FALSE, ...) {
      cat(".")
      u <- rnorm(2*nind)
      lam1 <- lambda[1]
      lam2 <- lambda[2]
      L <- t(chol(matrix(c(lam1,sqrt(lam1*lam2)*cor,sqrt(lam1*lam2)*cor,lam2),ncol=2)))
      v <- L %*% matrix(u,nrow=2)
      if (return.fit) plot(t(v))
      lfh <- data.frame(id=seq_len(nind), id2=seq_len(nind), feco= rpois(nind, lambda = exp(1+v[1,])), 
                        growth=rgamma(nind,shape=1/0.2, scale=0.2*exp(1+v[2,]))) # mean=exp(1+v[2,]), var= 0.2*mean^2
      (fitlfh <- fitmv(submodels=list(list(feco ~ 1+(0+mv(1,2)|id), family=poisson()),
                                      list(growth ~ 1+(0+mv(1,2)|id), family=Gamma(log))),
                       data=lfh, method=method, ...))
      if (return.fit) {
        return(fitlfh)
      } else {
        corr <- VarCorr(fitlfh)[2,"Corr."]
        cat(corr)
        corr
      }
    }
    set.seed(123)
    (zut <- simfun(cor=-0.5,return.fit=TRUE))
    plot(zut$data)
    replicate(3L, simfun(3000))
    replicate(3L, simfun(10000))
  }
  
}

if (FALSE) {
  library(aster)
  data(echinacea)
  asNA.acyclic <- function(data, pred) {
    order_pred <- order(pred)
    predvars <- names(pred)
    ord_desc_vars <- names(pred[order_pred])
    for (ordered_it in order_pred) {
      descvar <- ord_desc_vars[ordered_it]
      predecessor <- pred[descvar]
      if (predecessor>0L) {
        predvar <- predvars[predecessor]
        # Predecessor variable being 0 (no survival, or no flowers) means
        #    that descendant variable cannot be observed:
        data[is.na(data[[predvar]]) | data[[predvar]]==0L, descvar] <- NA
      }
    }
    return(data)
  }
  # Note the names, important here:
  varpred <- c(ld02=0, ld03=1, ld04=2, fl02=1, fl03=2, fl04=3,
               hdct02=4, hdct03=5, hdct04=6)
  NAechin <- asNA.acyclic(echinacea, pred=varpred)
  
  yearvars <- c("ld02", "ld03", "ld04")
  byyear <- reshape(NAechin, varying = list(yearvars), direction = "long",
                    timevar = "varld",times = as.factor(yearvars), v.names = "respld")
  yearvars <- c("fl02", "fl03", "fl04")
  flbyyear <- reshape(NAechin, varying = list(yearvars), direction = "long",
                      timevar = "varfl",times = as.factor(yearvars), v.names = "respfl")
  yearvars <- c("hdct02", "hdct03", "hdct04")
  hdctbyyear <- reshape(NAechin, varying = list(yearvars), direction = "long",
                        timevar = "varhdct",times = as.factor(yearvars), v.names = "resphdct")
  byyear$varfl <- flbyyear$varfl
  byyear$varhdct <- hdctbyyear$varhdct
  byyear$respfl <- flbyyear$respfl
  byyear$resphdct <- hdctbyyear$resphdct
  byyear <- byyear[ ,c(12,4:6,10,13,14,11,15,16)]
  head(byyear)
  
  
  ## ----asMM, size="small"---------------------------------------------------------------------------------------------------------------------------------------
  asMM <- fitmv(submodels=list(
    list(respld ~ varld + nsloc + ewloc+(1|id), family=binomial()),
    list(respfl ~ varfl + nsloc + ewloc+(1|id), family=binomial()),
    list(resphdct ~ varhdct + nsloc + ewloc+pop, family=Tpoisson())),
    data=byyear)
  

  if (requireNamespace("multcomp", quietly = TRUE)) {
    library(multcomp)
    #summary(glht(asMM,mcp("varld" = "Tukey"), coef.=fixef.HLfit)) # documented limitation
  }
}

summary(warnings())
  