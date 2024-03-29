cat(crayon::yellow("\ntest old examples and new tests:\n"))
# spaMM

data("scotlip") ## loads 'scotlip' data frame, but also 'Nmatrix'

# including a test of mgcv::negbin handling mechanism
# and apparently there are not so many tests of .solve_v_h_IRLS() in routine tests.
(hl <- try(fitme(I(1+cases)~I(prop.ag/10)+offset(log(expec))+adjacency(1|gridcode),
           family=negbin(), adjMatrix=Nmatrix, data=scotlip), silent=TRUE)) 
if (inherits(hl,"try-error")) {
  if (length(grep("mgcv",attr(hl,"condition")$message))) {
    message(".checkRespFam() detected that spaMM::negbin is masked by mgcv::negbin. Refitting as part of tests...")
    (hl <- try(fitme(I(1+cases)~I(prop.ag/10)+offset(log(expec))+adjacency(1|gridcode),
                     family=spaMM::negbin(), adjMatrix=Nmatrix, data=scotlip))) 
  } else stop(hl) # other unhandled error
}

(hl1 <- corrHLfit(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
          data=scotlip,family=poisson(),
          adjMatrix=Nmatrix)) ## 1D optimization -> optimize
testthat::expect_equal(hl1$APHLs$p_v,-161.5140,tolerance=1e-4)

(hl2 <- HLCor(cases~I(prop.ag/10) +adjacency(1|gridcode)+offset(log(expec)),
            data=scotlip,family=poisson(), adjMatrix=Nmatrix))
testthat::expect_equal(hl2$APHLs$p_v,-161.5141,tolerance=1e-4)
testthat::expect_true(diff(range(AIC(hl2,verbose=FALSE)-AIC(hl1,verbose=FALSE)))<1e-2)

data("salamander")
hl <- HLfit(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male),family=binomial(),
      rand.family=list(gaussian(),Beta(logit)),data=salamander,HLmethod="ML",control.HLfit = list(LevenbergM=FALSE))

testthat::expect_equal(hl$APHLs$p_v,-238.715,tolerance=1e-3)

## Nested effects
# lmer syntax allowing several degrees of nesting
hl <- HLfit(cbind(Mate,1-Mate)~1+(1|Female/Male),
      family=binomial(),rand.family=Beta(logit),data=salamander,HLmethod="ML",control.HLfit = list(LevenbergM=FALSE))

testthat::expect_equal(hl$APHLs$p_v,-243.6668,tolerance=1e-4)

# A syntax described in ?formula ## removed from the example()
hl <- HLfit(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male %in% Female),
            ranFix=list(lambda=c('Female' = 0.127517,'Male %in% Female' = 4.64595e-07)),
            family=binomial(),rand.family=Beta(logit),data=salamander,HLmethod="ML",control.HLfit = list(LevenbergM=FALSE))

testthat::expect_equal(hl$APHLs$p_v,-243.6668,tolerance=1e-4)

### check NULL auglinmodblob
d <- data.frame(y = 1:10)
summary(fitme(y ~ 0, data = d))


# test of scaling of ranCoef predictor
(ranSlope1 <- fitme(I(1+cases)~I(prop.ag/10)+adjacency(0+expec|gridcode),
                    family=poisson(), adjMatrix=Nmatrix, data=scotlip))
scotlip$verif <- scotlip$expec/2
(ranSlope2 <- fitme(I(1+cases)~I(prop.ag/10)+adjacency(0+verif|gridcode),
                    family=poisson(), adjMatrix=Nmatrix, data=scotlip))
testthat::expect_true(abs(ranSlope2$lambda-4*ranSlope1$lambda)<1e-5)

# test dfs
data("wafers")
testthat::expect_true(sum(unlist(HLfit(y~X1+(X2|batch), data=wafers)$dfs))==6L) # 3 df in p_lambda for ranCoefs
testthat::expect_true(sum(unlist(fitme(y~X1+(X2|batch)+(X2|batch), # stupid formula but effective test
                                       data=wafers, fixed=list(ranCoefs=list("1"=c(NA, -0.1, NA))))$dfs))==8L) # 5 df in p_lambda for ranCoefs
# check with mixed inner and outer-estimated ranefs (=> bug before v3.12.4):
testthat::expect_true(fitme(y~X1+(X2|batch)+(1|batch), 
                            data=wafers, fixed=list(ranCoefs=list("1"=c(NA, -0.1, NA))))$dfs$p_lambda==3L) # 2 df in p_lambda for ranCoefs
testthat::expect_true(fitme(y~X1+(X2|batch)+(1|batch), 
                            data=wafers)$dfs$p_lambda==4L) # 3 df in p_lambda for ranCoefs
# testthat::expect_true(HLfit(y~X1+(X2|batch)+(1|batch), data=wafers)$dfs$p_lambda==4L) # 3 df in p_lambda for ranCoefs (slow)

if (FALSE) { # examples from update.Rd in handy test form
  ## First the fit to be updated:
  wFit <- HLfit(y ~X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
                resid.model = ~ X3+I(X3^2) ,data=wafers)
  
  newresp <- simulate(wFit)
  update_resp(wFit,newresp=newresp)
  
  # For estimates given by Lee et al., Appl. Stochastic Models Bus. Ind. (2011) 27:  315-328:
  # Refit with given beta or/and phi values:
  
  betavals <- c(5.55,0.08,-0.14,-0.21,-0.08,-0.09,-0.09)
  # reconstruct fitted phi value from predictor for log(phi)
  Xphi <- with(wafers,cbind(1,X3,X3^2)) ## design matrix
  phifit <- exp(Xphi  %*%  c(-2.90,0.1,0.95))
  upd_wafers <- wafers
  upd_wafers$off_b <- wFit$`X.pv` %*% betavals
  update(wFit,formula.= . ~ offset(off_b)+(1|batch), data=upd_wafers,
         ranFix=list(lambda=exp(-3.67),phi=phifit))
  
  ## There are subtlety in performing REML fits of constrained models,
  ##   illustrated by the fact that the following fit does not recover 
  ##   the original likelihood values, because dispersion parameters are
  ##   estimated but the REML correction changes with the formula:
  upd_wafers$off_f <- wFit$`X.pv` %*%  fixef(wFit) ## = predict(wFit,re.form=NA,type="link")
  (diffwFit <- update(wFit,formula.= . ~ offset(off_f)+(1|batch), data=upd_wafers))
  testthat::expect_equal(logLik(diffwFit),c(p_bv=-1159.594584758636))
  ## To maintain the original REML correction, Consider instead
  (rewFit <- update(wFit,formula.= . ~ offset(off_f)+(1|batch), data=upd_wafers,
         REMLformula=formula(wFit)))  ## recover original p_v and p_bv     
  testthat::expect_true(diff(range(unlist(wFit$APHLs)-unlist(rewFit$APHLs)))<1e-4)
  ## Alternatively, show original wFit as differences from betavals:  
  (rerewFit <- update(wFit,formula.= . ~ . +offset(off_f), data=upd_wafers))
  testthat::expect_true(diff(range(fixef(rerewFit)))<1e-14)
}

{ # check of REMLformula and keepInREML
  l1 <- logLik(unconstr <- fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), data=wafers, family=Gamma(log), 
        method="REML"))
  l2 <- logLik(fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), data=wafers, family=Gamma(log), 
        method="REML", etaFix=list(beta=structure(fixef(unconstr), keepInREML=TRUE))))
  testthat::test_that("whether keepInREML uses model formula as default REMLformula",testthat::expect_true(l1-l2<1e-14))
  
  l3 <- logLik(suppressMessages(fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), data=wafers, family=Gamma(log), 
        method="REML", etaFix=list(beta=structure(fixef(unconstr), keepInREML=TRUE)),
        REMLformula=y ~(1|batch))))
  l4 <- logLik(fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), data=wafers, family=Gamma(log), 
                    method="ML"))
  testthat::test_that("whether keepInREML correctly handles nondefault REMLformula",
                      testthat::expect_true(l3-l4<1e-14))
}

# Code added to check model frame issues when updating a model with variable in resid model not in mean response model (+ syntax I() ). 
data("wafers")
## Gamma GLMM with log link
m1 <- HLfit(y ~X1+X2+I(X2^2),family=Gamma(log),
            resid.model = ~ X3+I(X3^2) ,data=wafers,method="ML") 
update_resp(m1,newresp=simulate(m1))

# residual dispersion model with partial etaFix
(wfit <- fitme(y ~ X1+X2+X1*X3+X2*X3+I(X2^2), family=Gamma(log),
               rand.family=inverse.Gamma(log), 
               resid.model = list(formula= ~ X3+I(X3^2)+(1|batch), 
                                  fixed=list(lambda=0.4), # let's check that too
                                  etaFix=list(beta=c("(Intercept)"=-2.93))) ,
               data=wafers))
testthat::test_that("whether partially-fixed etaFix works in resid.model",
                    testthat::expect_equal(fixef(wfit$resid_fit)[[1L]],-2.93))
testthat::test_that("whether fixed lambda works in resid.model",
                    testthat::expect_equal((wfit$resid_fit$lambda)[[1L]],0.4))

# Handling of ~ . :
testthat::test_that(
  "Warning expected for ~ . ", 
  {
    warn <- ""
    withCallingHandlers(
      fitme(Sepal.Length ~ . -Species + (1|Species), data=iris), warning=function(w) {
        warn <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      })
    testthat::expect_true(
      warn=="It looks like there is a '.' in the RHS of the formula.\n Fitting may be successful, but post-fit functions such as predict() will fail."
    )}
)

# Handling of ~ . :
testthat::test_that(
  "Warning expected for Earth + .(.|latitude+longitude) ", 
  {
    warn <- ""
    withCallingHandlers(
      fitme(migStatus ~ 1 + Matern(1|latitude+longitude),data=blackcap,
            method="ML", fixed=list(nu=0.5,phi=1e-6),
            control.dist=list(dist.method="Earth")), warning=function(w) {
        warn <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      })
    testthat::expect_true(
      warn=="Hmm... the first coordinate should be longitude, but seems to be latitude."
    )}
)
