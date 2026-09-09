cat(cli::col_yellow("\ntest pois4mlogit:\n"))
if (FALSE) { 
  source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-data_paternity.R")
  source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-paternity-interaction.R")
  source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-script-TC-20251030.R") 
  source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-composite-antisym.R") 
  # some of the fits are slow: (~1mn total):
  source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-paternity-3ranefs.R") 
  # This one is definitely slower (2-3 mn ?):
  if (FALSE) source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-confint.R") 
} else  cat(cli::bg_green(cli::col_black("\nRun private test-pois4mlogit-... R files once in a while.")))


#### Fitting a binomial(logit) model by a bivariate poisson(log) surrogate:
## Toy data 
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
}

if (spaMM.getOption("example_maxtime")>15) { 
  # using a small tol value allows more precise compar of numInfo's
  { # without ranCoefs
    (BbyB <- fitme(cbind(yellow,blue) ~ 0+phenotype+(1|grp), family = binomial(), 
                   data=toydata[1:5,]))
    BbyB$v_h # -1.449113 -2.977750 -1.175442  1.570128 -1.192737
    
    # *is equivalent to*:   
    
    (BbyP <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
      list(blue ~ offset(.dynoffset) + 0, family = poisson())), 
      control=list(p4m="oH"), progress=1L, verbose=c(TRACE=FALSE), tol=1e-7,
      data = toydata[1:5,], types=c("yellow","blue")))
    BbyP$v_h # 
    testthat::test_that('deparse(getCall(BbyP)[[1]])=="pois4mlogit"',
                        testthat::expect_true(deparse(getCall(BbyP)[[1]])=="pois4mlogit"))

    (BbyP_sp <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
      list(blue ~ offset(.dynoffset) + 0, family = poisson())), tol=1e-7,
      control=list(p4m="oH"), progress=1L, verbose=c(TRACE=FALSE), control.HLfit=list(algebra="spprec"),
      data = toydata[1:5,], types=c("yellow","blue")))
    
    (BbyPlam5 <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
      list(blue ~ offset(.dynoffset) + 0, family = poisson())),  fixed=list(lambda=5),
      control=list(p4m="oH"), progress=1L, verbose=c(TRACE=FALSE), tol=1e-7,
      data = toydata[1:5,], types=c("yellow","blue")))
    
    
    (BbyPdyndyn <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
      list(blue ~ offset(.dynoffset) + 0, family = poisson())),  fixed=list(lambda=5),
      control=list(p4m="W"), progress=1L, verbose=c(TRACE=FALSE), tol=1e-7,
      data = toydata[1:5,], types=c("yellow","blue"), to.long=FALSE))
    (crit <- abs(diff(range(logLik(BbyPlam5),logLik(BbyPdyndyn)))))
    testthat::test_that(paste0("p4m='H' works (for moderate lambda!)",signif(crit,4)," >1e-5"),
                        testthat::expect_true(crit<1e-5) )
    
    
    # see test-pois4mlogit-confint.R for tests of confint procedure
    
    if (FALSE) { # 112s (on "occupied" CPU), spcorr by default, & logL differs at tenth decimal from BbyP
      (BbyP_LM <- pois4mlogit(submodels = list(
        list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
        list(blue ~ offset(.dynoffset) + 0, family = poisson())), 
        control=list(p4m="oH"), progress=1L, verbose=c(TRACE=TRUE), control.HLfit=list(LevenbergM=TRUE),
        data = toydata[1:5,], types=c("yellow","blue"), n_iter=100L, tol=1e-7)) 
    }
    
    (infoB <- numInfo(BbyB))
    (infoP <- numInfo(BbyP)) # accurate in absolute terms, but not so much in relative ones. 
    (crit <- max(abs((infoB-infoP)/(0.001+abs(infoB))))) # becomes relative for small values
    testthat::test_that(paste0("Whether numInfo() gives correct result for BbyP:",signif(crit,4)," >1e-3"),
                        testthat::expect_true(crit<1e-3) )
    (infoP_sp <- numInfo(BbyP_sp))
    (crit <- max(abs((infoB-infoP_sp)/(0.001+abs(infoB)))))
    testthat::test_that(paste0("Whether numInfo() gives correct result for BbyP_sp:",signif(crit,4)," >1e-3"),
                        testthat::expect_true(crit<1e-3) )
    
    simulate(BbyP, nsim=3L)
    
    if (FALSE) { # REML...
      # the p4m fit matches the *outer-optimized fitme* one (expected), which is now the default:
      (fitme(cbind(yellow,blue) ~ 0+phenotype+(1|grp), family = binomial(), 
             data=toydata[1:5,], method="REML"))
      (pois4mlogit(submodels = list(
        list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
        list(blue ~ offset(.dynoffset) + 0, family = poisson())), 
        control=list(p4m="oH"), progress=1L, verbose=c(TRACE=FALSE), method="REML",
        data = toydata[1:5,], types=c("yellow","blue"), n_iter=100L, tol=1e-7))
    }
    
  } 
  
  { # with ranCoefs
    data(clinics)
    (fitClinics <- fitme(cbind(npos,nneg)~1+treatment+(treatment|clinic),
                         family=binomial(),data=clinics))
    
    (clinicp4m <- pois4mlogit(submodels = list(
      list(npos ~ offset(.dynoffset) + 1+ treatment+(treatment|clinic), family = poisson()),
      list(nneg ~ offset(.dynoffset) + 0, family = poisson())), control=list(p4m="oH"), 
      progress=1L+interactive(), verbose=c(TRACE=FALSE),
      lower=list(ranCoefs=list("1"=c(0,-0.9,0))), 
      upper=list(ranCoefs=list("1"=c(200,0.9,200))), # control.HLfit=list(LevenbergM=FALSE),
      data = clinics, types=c("npos","nneg"), n_iter=100L, tol=c(1e-5,1e-5)))
  }
  
} else {cat(cli::bg_green(cli::col_black("Mixed-effect models not run in 'fast' tests.")))}

## Standard binomial fit (purple flowers are ignored here)
(byB <- fitme(cbind(yellow,blue) ~ phenotype, family = binomial(), 
      data=toydata))


{
  ## Surrogate fit
  (byP2 <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
    list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
    data = toydata, types=c("yellow","blue")))
  
  #### Recover summaries of the binomial fit from the surrogate fit:
  
  ## Coefficients of binomial model recovered as:
  fixef(byP2)[1:2]-fixef(byP2)[3:4]
  
  ## SEs of coefficients of binomial model recovered as:  
  {
    P2B <- rbind(c(1,0,-1,0),c(0,1,0,-1))
    P2B %*% vcov(byP2) %*% t(P2B)
  } 
  # practically equivalent to 
  vcov(byB)
  
  ## logLiks
  str(long2 <- reshape2long(toydata, c("yellow","blue")))
  
  # Fits on long data:
  (byBlong  <- fitme(cbind(yellow,blue) ~ phenotype, family = binomial(), 
                     data=long2))
  (byPlong <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
    list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
    data = toydata, to.long=TRUE, types=c("yellow","blue")))
  
  crit <- diff(range(logLik(byBlong), logLik(byPlong)))
  testthat::test_that("Consistency between to.long=FALSE and =TRUE", 
                      testthat::expect_true(crit < 1e-10))

  # This illustrates that the logLiks of different surrogate models differ 
  # only by a constant independent of the fitted model. 
  # Likelihood ratios between surrogate models are then 
  # identical to likelihood ratios between corresponding binomial models,
  # provided a single data format is used throughout the comparisons.
  
}   
     
## NA handling: cf test-simulate

if (FALSE) { ## trying to test identifiability... but some fits converge immediately...
  (unident_ranef <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
    list(blue ~ offset(.dynoffset) + 0+(1|grp), family = poisson())), 
    control=list(p4m="oH"), progress=1L, verbose=c(TRACE=FALSE), 
    data = toydata[1:5,], types=c("yellow","blue")))
  # => low fitted variance appears typical in this unidentifiable case. 
  # There must be a warning.
  
  X_4to3 <- 
    matrix(c(1,0,0,
             0,1,0,
             1,0,0,
             0,0,1), nrow=4, ncol=3, byrow=TRUE,
           dimnames=list(NULL, c("(Intercept)","phenotype_1","phenotype_2")))
  
  (unidentif <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
    list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
    data = toydata, X2X=X_4to3, types=c("yellow","blue"), progress=1L+interactive()))
  
  if (FALSE) {
    # low_pot bug... but before that, .makeMatp4m() generates a rank deficient sXaug 
    # (correct given the wrong formulas)
    # so let's say this does not count as worth debugging.
    (unidentif <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + phenotype+(1|grp), family = poisson()),
      list(blue ~ offset(.dynoffset) + phenotype+(1|grp), family = poisson())),
      data = toydata, X2X=X_4to3, types=c("yellow","blue"), progress=1L+interactive()))
  }
  if (FALSE) {
    # Ultimately same cause, different bug. But after many more iterations + diagnosis
    (unidentif <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
      list(blue ~ offset(.dynoffset) + phenotype+(1|grp), family = poisson())),
      data = toydata, X2X=X_4to3, types=c("yellow","blue"), progress=1L+interactive()))
  }
}

