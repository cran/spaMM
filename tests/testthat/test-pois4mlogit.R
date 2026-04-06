cat(cli::col_yellow("\ntest pois4mlogit:\n"))
if (FALSE) { 
  source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-data_paternity.R")
  source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-paternity-interaction.R")
  source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-script-TC-20251030.R") 
  # some of the fits are slow: (~1mn total)
  source("D:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/test-pois4mlogit-paternity-3ranefs.R") 
} else  cat(cli::bg_green(cli::col_black("\nRun private test-pois4mlogit...R once in a while.")))


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

if (spaMM.getOption("example_maxtime")>60) { # cf D:/home/francois/travail/stats/spaMMplus/devel/pois4mlogit/p4m_DEBUG.R
  { # without ranCoefs
    (BbyB <- fitme(cbind(yellow,blue) ~ 0+phenotype+(1|grp), family = binomial(), 
                   data=toydata[1:5,]))
    BbyB$v_h # -1.449113 -2.977750 -1.175442  1.570128 -1.192737
    
    # *is equivalent to*:   
    
    (BbyP <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
      list(blue ~ offset(.dynoffset) + 0, family = poisson())), 
      control=list(p4m="oH"), progress=1L, verbose=c(TRACE=FALSE), 
      data = toydata[1:5,], types=c("yellow","blue"), n_iter=100L, tol=1e-7, to.long=FALSE))
    BbyP$v_h # 
    testthat::test_that('deparse(getCall(BbyP)[[1]])=="pois4mlogit"',
                        expect_true(deparse(getCall(BbyP)[[1]])=="pois4mlogit"))
    (BbyP_sp <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
      list(blue ~ offset(.dynoffset) + 0, family = poisson())), 
      control=list(p4m="oH"), progress=1L, verbose=c(TRACE=FALSE), control.HLfit=list(algebra="spprec"),
      data = toydata[1:5,], types=c("yellow","blue"), n_iter=100L, tol=1e-7, to.long=FALSE))
    if (FALSE) { # 205s & logL differs at fifth decimals
      (BbyP_LM <- pois4mlogit(submodels = list(
        list(yellow ~ offset(.dynoffset) + 0+ phenotype+(1|grp), family = poisson()),
        list(blue ~ offset(.dynoffset) + 0, family = poisson())), 
        control=list(p4m="oH"), progress=1L, verbose=c(TRACE=TRUE), control.HLfit=list(LevenbergM=TRUE),
        data = toydata[1:5,], types=c("yellow","blue"), n_iter=100L, tol=1e-7, to.long=FALSE)) 
    }
    
  } 
  
  { # with ranCoefs
    data(clinics)
    (fitClinics <- fitme(cbind(npos,nneg)~1+treatment+(treatment|clinic),
                         family=binomial(),data=clinics))
    
    (clinicp4m <- pois4mlogit(submodels = list(
      list(npos ~ offset(.dynoffset) + 1+ treatment+(treatment|clinic), family = poisson()),
      list(nneg ~ offset(.dynoffset) + 0, family = poisson())), control=list(p4m="oH"), progress=2L, verbose=c(TRACE=FALSE),
      lower=list(ranCoefs=list("1"=c(0,-0.9,0))), 
      upper=list(ranCoefs=list("1"=c(200,0.9,200))), # control.HLfit=list(LevenbergM=FALSE),
      data = clinics, types=c("npos","nneg"), n_iter=100L, tol=c(1e-5,1e-5), to.long=FALSE))
  }
  
}

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
  # The logLiks differ between the binomial and poisson surrogate fits because the combinatorial coefficients differ.
  # The relationship between these logLiks is simple when the long form of the data is used:
  str(long2 <- reshape2long(toydata, c("yellow","blue")))
  
  # Fits on long data:
  (byBlong  <- fitme(cbind(yellow,blue) ~ phenotype, family = binomial(), 
                     data=long2))
  (byPlong <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
    list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
    data = toydata, to.long=TRUE, types=c("yellow","blue")))
  
  # The logLik of the surrogate model is
  (logL1 <- sum(byPlong$eta*byPlong$y) - sum(exp(byPlong$eta))) # = logLik(byPlong)
  # while le logLik of the binomial one can be obtained as 
  (logL2 <- sum(byPlong$eta*byPlong$y)) # = logLik(byBlong)
  # the difference is 156  __= total multinomial sample size__
  crit <- diff(range(logL1,logL2-156, logLik(byPlong)))
  testthat::test_that("Consistency between to.long=FALSE and =TRUE", 
                      testthat::expect_true(crit < 1e-10))

  # This illustrates that the logLiks of different surrogate models differ 
  # only by a constant independent of the fitted model. 
  # Likelihood ratios between surrogate models are then 
  # identical to likelihood ratios between corresponding binomial models,
  # provided a single data format is used throughout the comparisons.
  
}   
     
{ ## NA handling
  toydataNA <- toydata
  toydataNA$pheno1 <- toydataNA$pheno2 <- toydataNA$pheno3 <- toydata$phenotype
  toydataNA$yellow[1] <- NA 
  toydataNA$pheno1[2] <- toydataNA$pheno2[2] <- NA 
  (byP3 <- pois4mlogit(submodels = list(
      list(yellow ~ offset(.dynoffset) + pheno1, family = poisson()),
      list(blue ~ offset(.dynoffset) + pheno2, family = poisson()),
      list(purple ~ offset(.dynoffset) + pheno3, family = poisson())),
      data = toydataNA, types=c("yellow","blue","purple")))
  length(pp1 <-predict(byP3, verbose=c(na=FALSE))) # 27 bc only second draw cannot be predicted
  rowSums(m1 <- matrix(pp1, ncol=3))
  length(pp2 <- predict(byP3, newdata=byP3$data,
                      verbose=c(na=FALSE))) # 27 bc only second draw cannot be predicted
  max(abs(range(matrix(pp2, ncol=3)-m1)))
  plot_effects(byP3, "pheno1", submodel=1)
}  

if (FALSE) { ## trying to test identifiability... but this converges immediately...
  X_4to3 <- 
    matrix(c(1,0,0,
             0,1,0,
             1,0,0,
             0,0,1), nrow=4, ncol=3, byrow=TRUE,
           dimnames=list(NULL, c("(Intercept)","phenotype_1","phenotype_2")))
  
  (unidentif <- pois4mlogit(submodels = list(
    list(yellow ~ offset(.dynoffset) + phenotype, family = poisson()),
    list(blue ~ offset(.dynoffset) + phenotype, family = poisson())),
    data = toydata, X2X=X_4to3, types=c("yellow","blue"), progress=2))
}

