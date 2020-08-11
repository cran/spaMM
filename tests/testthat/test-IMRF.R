cat(crayon::yellow("\ntest IMRF:")) 

data("blackcap")
(mrf <- HLCor(migStatus ~ 1 + multIMRF(1|latitude+longitude,margin=5,levels=1),data=blackcap,
              HLmethod="ML",ranPars=list(phi=1,lambda=1,corrPars=list("1"=list(kappa=1)))) )
(p1 <- predict(mrf)[2:3,])
(p2 <- predict(mrf, newdata=mrf$data[2:3,]))
testthat::expect_true(diff(range(p1-p2))<1e-12)
(p1 <- get_predVar(mrf)[2:3])
(p2 <- get_predVar(mrf, newdata=mrf$data[2:3,]))
testthat::expect_true(diff(range(p1-p2))<1e-12)

if (spaMM.getOption("example_maxtime")>5) {
  # Using 'hyper' to control fixed hyper-parameters
  (mrf <- fitme(migStatus ~ 1 + (1|pos) + multIMRF(1|latitude+longitude,margin=5,levels=2),data=blackcap,
                method="ML",fixed =list(phi=1,lambda=c("1"=0.666),hyper=list("1"=list(hy_kap=0.1,hy_lam=1)))) )
  (mrf <- fitme(migStatus ~ 1 + (1|pos) + multIMRF(1|latitude+longitude,margin=5,levels=2),data=blackcap,
                method="ML",fixed =list(phi=1,hyper=list("1"=list(hy_kap=0.1,hy_lam=1)))) )
  # Using 'hyper' to control initial hyper-parameters
  (mrf <- fitme(migStatus ~ 1 + multIMRF(1|latitude+longitude,margin=5,levels=2),data=blackcap,
                method="ML",fixed =list(phi=1),init=list(hyper=list("1"=list(hy_kap=0.1,hy_lam=1)))) )
  # *Independent* IMRF terms (often giving dubious results)
  (mrf <- HLCor(migStatus ~ 1 + IMRF(1|latitude+longitude,margin=5, nd=4L)
                + IMRF(1|latitude+longitude,margin=5, nd=7L),
                data=blackcap, HLmethod="ML",
                ranPars=list(phi=1,lambda=c(1/4,1/16),
                             corrPars=list("1"=list(kappa=0.1),"2"=list(kappa=0.1)))))
}
if (spaMM.getOption("example_maxtime")>1) {
  # trying combinations of fixed/ estimated parameters for predVar:
  # but most important for the predVar comparison with the (1|pos) moved
  (mrf <- fitme(migStatus ~ 1 + (1|pos) + multIMRF(1|latitude+longitude,margin=2,levels=2, coarse=4),data=blackcap,
                method="ML",fixed=list(lambda=1,hyper=list("1"=list(hy_kap=1)))) )
  p1 <- get_predVar(mrf)
  p2 <- get_predVar(mrf, re.form= ~multIMRF(1|latitude+longitude,margin=2,levels=2, coarse=4) + (1|pos)) # sinversion of ranefs
  testthat::expect_true(diff(range(p1-p2))<1e-12)
  get_predVar(mrf, re.form= ~multIMRF(1|latitude+longitude,margin=2,levels=2, coarse=4))
  simulate(mrf, re.form= ~multIMRF(1|latitude+longitude,margin=2,levels=2, coarse=4))
  #
  (mrf <- fitme(migStatus ~ 1 + (1|pos) + multIMRF(1|latitude+longitude,margin=2,levels=2, coarse=4),data=blackcap,
                method="ML",fixed=list(phi=1,hyper=list("1"=list(hy_kap=1)))) )
  p1 <- get_predVar(mrf)
  p2 <- get_predVar(mrf, re.form= ~multIMRF(1|latitude+longitude,margin=2,levels=2, coarse=4) + (1|pos)) # sinversion of ranefs
  testthat::expect_true(diff(range(p1-p2))<1e-12)
  get_predVar(mrf, re.form= ~multIMRF(1|latitude+longitude,margin=2,levels=2, coarse=4))
  simulate(mrf, re.form= ~multIMRF(1|latitude+longitude,margin=2,levels=2, coarse=4))
}
if (spaMM.getOption("example_maxtime")>6) {
  (mrf <- fitme(migStatus ~ 1 + (1|pos) + multIMRF(1|latitude+longitude,margin=2,levels=3, coarse=7),data=blackcap,
                method="ML",fixed=list(phi=0.1,lambda=0.1,hyper=list("1"=list(hy_kap=1)))) )
  set.seed(123)
  plot(predict(mrf),rowMeans(simulate(mrf,nsim=100,type="predVar", newdata=mrf$data, # using 'newdata' only to check more code
                                      variances=list(linPred=TRUE,disp=FALSE)))) 
  # (type="predVar" is not statistically meaningful here, but this makes a test of its code)
}

if (file.exists((privtest <- "C:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_other_pack/test-LatticeKrig.R"))) {
  source(privtest) 
}

