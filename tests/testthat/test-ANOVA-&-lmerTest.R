cat(crayon::yellow("\ntest lmerTest interface and other ANOVA tables:")) 
# See test-rank for additional tests of anova()

{ 
  { # LM (more LMs with sleepstudy data below)
    fit <- lm(sr ~ ., data = LifeCycleSavings)
    spfit <- fitme(sr ~ pop15+pop75+dpi+ddpi , data = LifeCycleSavings)
    testthat::test_that("whether lm() and fitme ML yield identical anova tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(fit))-unlist(anova(spfit)))))<1e-10))
    spfitRE <- fitme(sr ~ pop15+pop75+dpi+ddpi , data = LifeCycleSavings, method="REML")
    testthat::test_that("whether fitme ML and fitme REML yield identical anova tables for LM",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(spfit))-unlist(anova(spfitRE)))))<1e-10))
  }
  
  { # GLM
    dat <- data.frame(treatment=gl(3,3), outcome= gl(3,1,9), 
                      counts=c(18,17,15,20,10,20,25,13,12)) # showing data
    glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(), data=dat)
    spglm <- fitme(counts ~ outcome + treatment, family = poisson(), data=dat)
    # no F test for poisson
    testthat::test_that("whether fitme and poisson glm yield identical anova Cp tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(glm.D93, test = "Cp"))-unlist(anova(spglm, test = "Cp")))))<1e-10))
    testthat::test_that("whether fitme and poisson glm yield identical anova Chisq tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(glm.D93, test = "Chisq"))-unlist(anova(spglm, test = "Chisq")))))<1e-10))
    testthat::test_that("whether fitme and poisson glm yield identical anova Rao [score test] tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(glm.D93, test = "Rao"))-unlist(anova(spglm, test = "Rao")))))<1e-5)) # less accurate
    
    clotting <- data.frame(
      u = c(5,10,15,20,30,40,60,80,100),
      lot1 = c(118,58,42,35,27,25,21,19,18),
      lot2 = c(69,35,26,21,18,16,13,12,12))
    fit <- glm(lot1 ~ log(u), data = clotting, family = Gamma)
    spglm <- fitme(lot1 ~ log(u), data = clotting, family = Gamma, method="REML")
    # F test returned for Gamma, but large differences are expected from different phi estimation
    #testthat::test_that("whether fitme and poisson glm yield identical anova F tables",
    #                    testthat::expect_true(diff(range(na.omit(unlist(anova(fit, test = "F"))-unlist(anova(spglm, test = "F")))))<1e-3)) 
    testthat::test_that("whether fitme and Gamma glm yield identical anova Cp tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(fit, test = "Cp"))-unlist(anova(spglm, test = "Cp")))))<1e-3)) # some difference -- expected from different phi estimation
    # spglm <- fitme(lot1 ~ log(u), data = clotting, family = Gamma, method="ML") # differs more on Cp
    testthat::test_that("whether fitme and Gamma glm yield identical anova Chisq tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(fit, test = "Chisq"))-unlist(anova(spglm, test = "Chisq")))))<1e-10))
    testthat::test_that("whether fitme and Gamma glm yield identical anova Rao [score test] tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(fit, test = "Rao"))-unlist(anova(spglm, test = "Rao")))))<1e-5)) 
    testthat::test_that("whether fitme and Gamma glm yield identical drop1 F tables (apart from AICs)",
                        testthat::expect_true(diff(range(na.omit(unlist(drop1(fit, test = "F"))-unlist(drop1(spglm, test = "F")))[-(4:5)]))<1e-3)) 
    testthat::test_that("whether fitme and Gamma glm yield identical drop1 Chisq tables (apart from AICs and scaled dev)",
                        testthat::expect_true(diff(range(na.omit(unlist(drop1(fit, test = "Chisq"))-unlist(drop1(spglm, test = "Chisq")))[-(4:6)]))<1e-3)) 
    testthat::test_that("whether fitme and Gamma glm yield identical drop1 Rao [score test] tables (apart from AICs and and scaled dev)",
                        testthat::expect_true(diff(range(na.omit(unlist(drop1(fit, test = "Rao"))-unlist(drop1(spglm, test = "Rao")))[-(4:6)]))<1e-3)) # some fifference -- expected from different phi estimation
  }

  # LM and LMM
  if (spaMM.getOption("example_maxtime")>3) { # not an actual timing
    if(requireNamespace("lme4", quietly = TRUE)) {
      data("sleepstudy", package="lme4")
      
      { # LM
        lmf <- lm(formula = Reaction ~ Days+I(Days^2), data=sleepstudy)
        glmf <- glm(formula = Reaction ~ Days+I(Days^2), data=sleepstudy)
        splm <- fitme(formula = Reaction ~ Days+I(Days^2), data=sleepstudy, method="REML")
        anova(lmf)
        ano1 <- anova(glmf, test="F")
        ano2 <- anova(splm)
        crit <- diff(range(ano1[2:3,"Pr(>F)"]-ano2[1:2,"Pr(>F)"]))
        testthat::test_that("equivalent anova glm() and fitme()",
                            testthat::expect_true(crit<1e-10))
      }
      
      if(requireNamespace("lmerTest", quietly = TRUE)) { # LMM
        
        spfit <- fitme(Reaction ~ Days + (1|Subject), sleepstudy)
        spfit_lmlt <- as_LMLT(spfit) 
        ano3 <- lmerTest::contest(spfit_lmlt, L=c(1,0))
        spfit_lmlt <- as_LMLT(spfit, transf=FALSE) 
        ano2 <- lmerTest::contest(spfit_lmlt, L=c(1,0))
        fm <- lme4::lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
        ano1 <- lmerTest::contest(fm, L=c(1,0))
        crit <- abs(1-ano1[,"F value"]/ano2[,"F value"]) # p-value too low for informative comparison; compar of relative F values
        testthat::test_that("contest() results unchanged", testthat::expect_true(crit<1e-6)) 
        
        spfit <- fitme(Reaction ~ Days + (1+Days|Subject), sleepstudy)
        spfit_lmlt <- as_LMLT(spfit) 
        ano3 <- lmerTest::contest(spfit_lmlt, L=c(1,0)) 
        spfit_lmlt <- as_LMLT(spfit, transf=FALSE) 
        ano2 <- lmerTest::contest(spfit_lmlt, L=c(1,0)) 
        fm <- lme4::lmer(Reaction ~ Days + (1+Days|Subject), sleepstudy, REML=FALSE)
        ano1 <- lmerTest::contest(fm, L=c(1,0))
        crit <- diff(range(ano3[,"F value"],ano2[,"F value"], 1436.88833634)) # but slight difference with lmer; already observed with brute refit rather than hlcorcall hack
        testthat::test_that("equivalent contest() merMod and fitme()",
                            testthat::expect_true(crit<1e-6))
        
        fm <- lme4::lmer(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject), sleepstudy)
        spfit <- fitme(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject), sleepstudy, method="REML")
        spfit_lmlt <- as_LMLT(spfit) 
        ano3 <- anova(spfit_lmlt, type="1")
        ano2 <- anova(spfit_lmlt, type="1", transf=FALSE)
        ano1 <- anova(lmerTest::as_lmerModLmerTest(fm), type="1")
        crit <- diff(range(ano1[,"Pr(>F)"]-ano2[,"Pr(>F)"]))
        testthat::test_that("equivalent type-1 anova merMod and fitme()",
                            testthat::expect_true(crit<1e-6))
        ano3 <- anova(spfit_lmlt, type="2")
        ano2 <- anova(spfit_lmlt, type="2", transf=FALSE)
        ano1 <- anova(lmerTest::as_lmerModLmerTest(fm), type="2")
        crit <- diff(range(ano1[,"Pr(>F)"]-ano2[,"Pr(>F)"]))
        testthat::test_that("equivalent type-2 anova merMod and fitme()",
                            testthat::expect_true(crit<1e-6))
        ano3 <- anova(spfit_lmlt, type="3")
        ano2 <- anova(spfit_lmlt, type="3", transf=FALSE)
        ano1 <- anova(lmerTest::as_lmerModLmerTest(fm), type="3")
        crit <- diff(range(ano1[,"Pr(>F)"]-ano2[,"Pr(>F)"]))
        testthat::test_that("equivalent type-3 anova merMod and fitme()",
                            testthat::expect_true(crit<1e-6))
      } else message("Package 'lmerTest' not available for testing as_LMLT().")
    } 
  }
}
