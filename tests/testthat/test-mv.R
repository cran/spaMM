cat(crayon::yellow("\ntest-mv:"))

if (FALSE) {
  source(paste0(projpath(),"/package/tests/testthat/nestedFiles/test-mv-nested.R"))
} else {
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
  update_resp(zut1,newresp = simulate(zut1))
}

