cat(crayon::yellow("\ntest obsInfo:"))

# data("Salamanders", package = "glmmTMB") 
# (foo <- fitme(count  ~  spp  *  mined  +  (1  |site), data=Salamanders, family=negbin(link=log), method=c("ML","obs"), verbose=c(TRACE=F)))
# testthat::expect_equal(logLik(foo), c(p_v=-815.678491998 ))

data(scotlip)
(foo <- fitme(cases ~ I(prop.ag/10)+(1|gridcode),
      family=negbin(link=log), data=scotlip, method=c("ML","obs")))
testthat::expect_equal(logLik(foo), c(P_v=-181.60802361 ))


data(wafers)
(foo <- fitme(y ~1+(1|batch),family=Gamma(log),  #fixed=list(lambda=c("1"=0.01213306), phi=0.2255785), 
       data=wafers,method=c("ML","obs"), verbose=c(TRACE=F)))
testthat::expect_equal(logLik(foo), c(P_v=-1224.65219293 ))

{
  data("clinics")

  (foo <- fitme(cbind(npos,nneg)~1+(1|clinic), method=c("ML","obs"),
                family=binomial(cloglog),data=clinics))
  testthat::expect_equal(logLik(foo), c(P_v=-40.4244068868))
}

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
  (me1 <- fitme(formula=ly ~ 1+(1|batch), family=Gamma(log), method=c("ML","obs"), data=wafmv))
  (me2 <- fitme(formula=y3 ~ 1+(1|batch2), family=Gamma(log), method=c("ML","obs"), data=wafmv))
  (zut1 <- fitmv(submodels=list(mod1=list(formula=ly ~ 1+(1|batch), family=Gamma(log)),
                                mod2=list(formula=y3 ~ 1+(1|batch2), family=Gamma(log))), 
                 method=c("ML","obs"),
                 data=wafmv))
  testthat::expect_true(diff(range(logLik(me1)+logLik(me2),logLik(zut1)))<1e-5)
  testthat::expect_true(diff(range( predict(zut1, newdata=zut1$data)-predict(zut1)))<1e-14)
  testthat::expect_true(diff(range( get_predVar(zut1, newdata=zut1$data)-get_predVar(zut1)))<1e-14) 
  update_resp(zut1,newresp = simulate(zut1))
}

if (spaMM.getOption("example_maxtime")>3) {
  data("scotlip")
  adjfit <- fitmv(submodels=list(mod1=list(formula=cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),family=poisson()),
                                 mod2=list(formula=cases~I(prop.ag/10) +(1|gridcode)+offset(log(expec)),family=poisson(),rand.family=list(Gamma(log)))),
                  adjMatrix=Nmatrix,data=scotlip)
  std_lev <- hatvalues(adjfit, type="std", force=TRUE, which="both")
  testthat::test_that("Whether hatvalues() works on mv fit moreover, with obsInfo)",
                      testthat::expect_true(diff(range(c(std_lev$ranef -adjfit$lev_lam,std_lev$resid -adjfit$lev_phi)))<1e-14) )
}
