cat(crayon::yellow("\ntest LLMs:")) # there is an expanded version test-devel-LMM.R

TRACEv <- FALSE
# TRACEv <- interactive()

if (spaMM.getOption("example_maxtime")>39) {
  data(scotlip)
  
  { cat("check negbin2...")
    (nb2 <- fitme(cases~1+offset(log(expec))+(1|id),family=negbin2(), data=scotlip, method=c("ML","obs")) )# -181.3871 
    (nb2strict <- fitme(cases~1+offset(log(expec))+(1|id),family=negbin2(LLgeneric = FALSE), data=scotlip, method=c("ML","obs"))) # -181.3871 
    testthat::test_that("check negbin(2) mixed model",
                        testthat::expect_true(diff(c(range(logLik(nb2),logLik(nb2strict),-181.3870997165374)))<1e-10))
    
    if (FALSE) { # slow
      confint(nb,parm="(Intercept)") # with two fits approaching logLik(nb)-1.92073=-183.3078
      confint(nb2,parm="(Intercept)") # practically equivalent result
    }
    
    tp <- fitme(I(1+cases)~1+(1|id),family=Tpoisson(), data=scotlip, method=c("ML")) # -182.0702 # that the canonical link so "obs" is not distinct.
    testthat::test_that("check Tpoisson() mixed model",
                        testthat::expect_true(diff(c(range(logLik(tp),-182.07022358 )))<1e-5))
    # small numerical issues for Tnegbin with large shape, which should converge to Tpoisson irrespective of obsInfo since log is canonical for poisson:
    fitme(I(1+cases)~1+(1|id),family=negbin(trunc=0, shape=1e10), data=scotlip, method=c("ML")) # -182.0711 after -182.0691 and -182.0708 (Hexp must be ~Hobs as shape -> infty)
    fitme(I(1+cases)~1+(1|id),family=negbin2(trunc=0, shape=1e10), data=scotlip, method=c("ML","obs")) # -182.0691
    
    (tnb2 <- fitme(I(1+cases)~1+(1|id),family=negbin2(trunc=0), data=scotlip, method=c("ML","obs"))) # # -181.7404 
    testthat::test_that("check truncated negbin(2) mixed model",
                        testthat::expect_true(diff(c(range(logLik(tnb2),-181.7403819160504 )))<1e-10))
    
    (psqrt <- fitme(cases~1+(1|id),family=Poisson(link="sqrt"), data=scotlip, method=c("ML","obs"))) # -181.9246
    testthat::test_that("check poisson(sqrt) mixed model",
                        testthat::expect_true(diff(c(range(logLik(psqrt),-181.9246119586463  )))<1e-6))
    
    tnb2sqrt <- fitme(I(1+cases)~1+(1|id),family=negbin2(link="sqrt", trunc=0), data=scotlip, method=c("ML","obs")) #     -181.2723
    testthat::test_that("check truncated negbin(2)(sqrt) mixed model",
                        testthat::expect_true(diff(c(range(logLik(tnb2sqrt),-181.272264836057 )))<1e-10))
    tpsqrt <- fitme(I(1+cases)~1+(1|id),family=Tpoisson(link="sqrt"), data=scotlip, method=c("ML","obs")) #             -182.2697
    tpsqrtnb2 <- fitme(I(1+cases)~1+(1|id),family=negbin2(link="sqrt", trunc=0, shape=1e6), data=scotlip, method=c("ML","obs")) #     -182.2697 
    testthat::test_that("check truncated <poisson|negbin[2](shape=1e6)>(sqrt) mixed model",
                        testthat::expect_true(diff(c(range(logLik(tpsqrt),logLik(tpsqrtnb2),-182.269735148878  )))<1e-4))
  }
  
  { cat("check negbin1 (incl. test resid.model in locally available file)...")
    ### negbin 1 
    ## Fixed-effect model
    nb1llm <- fitme(cases~I(prop.ag/10)+offset(log(expec)),family=negbin1(shape=1/2.94), data=scotlip) # -176.0048
    testthat::test_that("check negbin1() LLM",
                        testthat::expect_true(diff(c(range(logLik(nb1llm),-176.004792436594  )))<1e-10))
    ## mixed model
    # Convenient example with negative weights for devel:
    nb1spc <- fitme(cases~1+(1|id),family=negbin1(), data=scotlip, verbose=c(TRACE=TRACEv)) #  -181.6869
    nb1dec <- fitme(cases~1+(1|id),family=negbin1(), data=scotlip, verbose=c(TRACE=TRACEv), control.HLfit=list(algebra="decorr")) #  -181.6869
    nb1spp <- fitme(cases~1+(1|id),family=negbin1(), data=scotlip, verbose=c(TRACE=TRACEv), control.HLfit=list(algebra="spprec")) #  -181.6869 
    spaMM.options(Hobs_Matrix_method= "def_sXaug_Matrix_QRP_CHM_scaled")
    nb1QRP <- fitme(cases~1+(1|id),family=negbin1(), data=scotlip, verbose=c(TRACE=TRACEv)) #  -181.6869
    spaMM.options(Hobs_Matrix_method= "def_sXaug_Matrix_CHM_H_scaled") # default
    testthat::test_that("check negbin1 mixed model",
                        testthat::expect_true(diff(c(range(logLik(nb1spc),logLik(nb1spp),logLik(nb1dec),logLik(nb1QRP),-181.6869459852397   )))<2e-6))
    
    logLik(glm(cases~I(prop.ag/10)+offset(log(expec)),family=poisson(), data=scotlip, method="spaMM_glm.fit", control=list(trace=FALSE)))
    logLik(glm(cases~I(prop.ag/10)+offset(log(expec)),family=negbin1(shape=10000), data=scotlip, method="llm.fit", control=list(trace=FALSE))) 
    
    if (file.exists((privtest <- "C:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_other_pack/test-resid.model.R"))) {
      source(privtest)  # negbin1 and negbin2: testing robustness to scaling of predictor variables, + local maximum issue in negbin2 fit
    } 
    
    { # diverging shape param
      # negbin1 is coherent with the poisson fit :
      onb1 <- fitme(cases~1+offset(log(expec))+(1|id),family=negbin1(), data=scotlip)  # -181.4601
      spaMM.options(Hobs_Matrix_method= "def_sXaug_Matrix_QRP_CHM_scaled")
      onb1QRP <- fitme(cases~1+offset(log(expec))+(1|id),family=negbin1(), data=scotlip) # -181.4601
      spaMM.options(Hobs_Matrix_method= "def_sXaug_Matrix_CHM_H_scaled") # default
      testthat::test_that("check negbin1 mixed model with diverging shape", # sensitive to .NB_shapeFn()
                          testthat::expect_true(diff(c(range(logLik(onb1QRP),logLik(onb1QRP),-181.4600881552111   )))<1e-6)) # (2.691147e-07)
      
    }
    
    if (FALSE) { # remarkable problem of  L-BFGS-B
      (tnb1 <- fitme(I(1+cases)~1+(1|id),family=negbin1(trunc=0), data=scotlip, verbose=c(TRACE=TRUE), control=list(optimizer="L-BFGS-B"),
                     init=list(NB_shape=0.2123))) # -181.7984   
    }
    
    (tnb1 <- fitme(I(1+cases)~1+(1|id),family=negbin1(trunc=0), data=scotlip, verbose=c(TRACE=TRACEv), 
                  init=list())) # -181.7984   
    tnb1de <- fitme(I(1+cases)~1+(1|id),family=negbin1(trunc=0), data=scotlip, verbose=c(TRACE=TRACEv), 
                    init=list(),control.HLfit=list(algebra="decorr")) # -181.7984  
    tnb1sp <- fitme(I(1+cases)~1+(1|id),family=negbin1(trunc=0), data=scotlip, verbose=c(TRACE=TRACEv), 
                    init=list(),control.HLfit=list(algebra="spprec")) # -181.7984  
    testthat::test_that("check truncated negbin1()",
                        testthat::expect_true(diff(c(range(logLik(tnb1),logLik(tnb1de),logLik(tnb1sp),-181.7984443645773  )))<2e-6))
    
    {
      # Low-shape torture test: many nonSPD 
      # the default method finds nonSPD matrices including at the attained fit. Then the Lev_M results are hard to characterize. 
      tnb1 <- fitme(I(1+cases)~1+(1|id),family=negbin1(trunc=0,shape=0.076247), data=scotlip, verbose=c(TRACE=interactive()), 
                    fixed=list(lambda=2.20016)) 
      tnb1de <- fitme(I(1+cases)~1+(1|id),family=negbin1(trunc=0,shape=0.076247), data=scotlip, verbose=c(TRACE=interactive()), 
                      fixed=list(lambda=2.20016), control.HLfit=list(algebra="decorr")) 
      tnb1sp <- fitme(I(1+cases)~1+(1|id),family=negbin1(trunc=0,shape=0.076247), data=scotlip, verbose=c(TRACE=interactive()), 
                      fixed=list(lambda=2.20016), control.HLfit=list(algebra="spprec")) 
      if (FALSE) {
        # Further, in QRP_CHM LevM does not use the same Hessian matrix. Effect was first seen in levM_v_h (iter 28 gainratio and newlik (hlik): use TRACE=2 to see it), and the final result differs.    
        spaMM.options(Hobs_Matrix_method= "def_sXaug_Matrix_QRP_CHM_scaled")
        tnb1CHM <- fitme(I(1+cases)~1+(1|id),family=negbin1(trunc=0,shape=0.076247), data=scotlip, verbose=c(TRACE=TRUE), 
                         fixed=list(lambda=2.20016)) 
        spaMM.options(Hobs_Matrix_method= "def_sXaug_Matrix_CHM_H_scaled") # default
      }
      
    }
  }
  
  { cat("check beta_resp (incl. resid.model)...")
    set.seed(123)
    beta_dat <- data.frame(y=runif(100),grp=sample(2,100,replace = TRUE), x_het=runif(100),
                           y2=runif(100))
    beta_llm <- fitme(y ~1, family=beta_resp(), data= beta_dat) # 0.08819
    testthat::test_that("check beta LLM",
                        testthat::expect_true(diff(c(range(logLik(beta_llm),0.0881926514984  )))<1e-10))
    beta_llmm <- fitme(y ~1+(1|grp), family=beta_resp(), data= beta_dat) # 0.08908       
    beta_llmm_fix <- fitme(y ~1+(1|grp), family=beta_resp(prec=residVar(beta_llmm, which="fam_parm")), data= beta_dat) # 0.08908       
    testthat::test_that("check beta LLMM",
                        testthat::expect_true(diff(c(range(logLik(beta_llmm),logLik(beta_llmm_fix),logLik(beta_llmm),0.0890820607358  )))<1e-6))
    beta_llmm_het <- fitme(y ~1+(1|grp), family=beta_resp(), data= beta_dat, resid.model= ~ x_het) # 0.08908       
    testthat::test_that("check LRT with beta_resp residual-dispersion", # test of p-value... checks that correct df number (1)
                        testthat::expect_true(abs(diff(c(LRT(beta_llmm,beta_llmm_het)$basicLRT$p_value,0.231086)))<1e-6))
    testthat::test_that("check LRT with beta_resp residual-dispersion", # test of p-value... checks that correct df number (2)
                        testthat::expect_true(abs(diff(c(LRT(beta_llmm_fix,beta_llmm_het)$basicLRT$p_value,0.4881745)))<1e-6))
    
    { ### check interactions with other packages
      # Here the one case where the family's own $simulate is called:
      beta_glm <- glm(y ~1, family=beta_resp(prec=1), data= beta_dat, method="llm.fit")
      str(simulate(beta_glm))
      beta_glmw <- glm(y ~1, family=beta_resp(prec=1), data= beta_dat, method="llm.fit", weights = rep(1:2,50))
      str(simulate(beta_glmw))
      # see also test-DHARMa
    }
    
    { # mv fit, permutation test, with one resid.model for beta response
      (zut1 <- fitmv(list(list(y ~1+(1|grp), family=beta_resp(), resid.model = ~x_het),
                 list(y2 ~1+(1|grp), family=beta_resp())), 
            data= beta_dat))
      (zut2 <- fitmv(list(list(y2 ~1+(1|grp), family=beta_resp()),
                 list(y ~1+(1|grp), family=beta_resp(), resid.model = ~x_het)), 
            data= beta_dat))
      testthat::test_that("Permutation check for mv fit with resid.model for beta response", 
                          testthat::expect_true(diff(range(logLik(zut1),logLik(zut2)))<1e-8))
      get_fittedPars(zut1)
      crit <- diff(range(numInfo(zut1)-numInfo(zut1,transf=TRUE))) # warnings OK for lambda at its lower bound
      testthat::test_that("numInfo(.,transf) for mv fit with resid.model for beta response", 
                          testthat::expect_true(crit<3e-6))
      # example continued for prediction with missing data, in fitmv.Rd
    }
  }

}



