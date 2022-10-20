cat(crayon::yellow("\ntest lmerTest interface and other ANOVA tables:")) 

{ 
    fit <- lm(sr ~ ., data = LifeCycleSavings)
    spfit <- fitme(sr ~ pop15+pop75+dpi+ddpi , data = LifeCycleSavings)
    testthat::test_that("whether lm() and fitme ML yield identical anova tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(fit))-unlist(anova(spfit)))))<1e-10))
    spfitRE <- fitme(sr ~ pop15+pop75+dpi+ddpi , data = LifeCycleSavings, method="REML")
    testthat::test_that("whether fitme ML and fitme REML yield identical anova tables for LM",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(spfit))-unlist(anova(spfitRE)))))<1e-10))
    spfitRE <- fitme(sr ~ pop15+pop75+dpi+ddpi , data = LifeCycleSavings, method="REML")
    
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
    testthat::test_that("whether fitme and poisson glm yield identical anova Cp tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(fit, test = "Cp"))-unlist(anova(spglm, test = "Cp")))))<1e-3)) # some fifference -- expected from different phi estimation
    # spglm <- fitme(lot1 ~ log(u), data = clotting, family = Gamma, method="ML") # differs more on Cp
    testthat::test_that("whether fitme and poisson glm yield identical anova Chisq tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(fit, test = "Chisq"))-unlist(anova(spglm, test = "Chisq")))))<1e-10))
    testthat::test_that("whether fitme and poisson glm yield identical anova Rao [score test] tables",
                        testthat::expect_true(diff(range(na.omit(unlist(anova(fit, test = "Rao"))-unlist(anova(spglm, test = "Rao")))))<1e-5)) 
    
    # ___F I X M E___? Add LMM tests (there are some in private tests using sleepstudy data)
  }
