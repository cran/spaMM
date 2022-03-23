## Tidy tests for routine checks

# library(spaMM)
# options(error=recover)

if (spaMM.getOption("example_maxtime")>10) {
  cat(crayon::cyan("\ntest-corrFamilies.R"))
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
  }
  
  {
    cat(crayon::yellow("MaternIMRFa; ")) 
    { # create IMRF model
      ## Creating the mesh 
      mesh <- INLA::inla.mesh.2d(loc = blackcap[, c("longitude", "latitude")], 
                                 cutoff=30,
                                 max.edge = c(3, 20)) 
      mesh$n ## 40
      matern <- INLA::inla.spde2.matern(mesh)
    }
    { 
      (fit_IMRF <- fitme(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern), # verbose=c(TRACE=TRUE), 
                         fixed=list(phi=c(0.02)),
                         data=cap_mv))
      (fit_reg <- fitme(status2 ~ 1+ MaternIMRFa(1|longitude+latitude, mesh=mesh, fixed=c(alpha=2)), # verbose=c(TRACE=TRUE), 
                        fixed=list(phi=c(0.02)), 
                        data=cap_mv))
      # ! outside points: rowSums(fit_reg$ranef_info$sub_corr_info$AMatrices[[1]])
      # alternative syntax
      (fit_cF <- fitme(status2 ~ 1+ corrFamily(1|longitude+latitude), # verbose=c(TRACE=TRUE), 
                       fixed=list(phi=c(0.02)), covStruct=list(corrFamily=MaternIMRFa(mesh=mesh, fixed=c(alpha=2))),
                       data=cap_mv))
      testthat::test_that("fitme MaternIMRFa OK",
                          testthat::expect_true(diff(range(logLik(fit_IMRF), logLik(fit_cF), logLik(fit_reg)))<1e-5))
      p1 <- predict(fit_IMRF)[c(3,1,3),]
      p2 <- predict(fit_reg, newdata=cap_mv[c(3,1,3),])
      p3 <- predict(fit_cF, newdata=cap_mv[c(3,1,3),])
      testthat::test_that("predict MaternIMRFa OK",
                          testthat::expect_true(diff(range(p1-p2, p1-p3))<1e-6))
      p1 <- get_predVar(fit_IMRF)[c(3,1,3)]
      p2 <- get_predVar(fit_reg, newdata=cap_mv[c(3,1,3),])
      p3 <- get_predVar(fit_cF, newdata=cap_mv[c(3,1,3),])
      testthat::test_that("get_predVar MaternIMRFa OK",
                          testthat::expect_true(diff(range(p1-p2, p1-p3))<1e-6))
      cat(crayon::blue("Two expected warnings as bc fitmv with unregistered corrFamily:"))
      (zut_IMRF <- fitmv(submodels=list(mod1=list(migStatus ~ 1),
                                        mod2=list(status2 ~ 1+ IMRF(1|longitude+latitude, model=matern))), # verbose=c(TRACE=TRUE), 
                         fixed=list(phi=c(0.02,0.02)), covStruct=list(corrFamily=MaternIMRFa(mesh=mesh, fixed=c(alpha=2))),
                         data=cap_mv))
      
      # Warning at submodel preprocessing, as expected since global 'covStruct' arg not part of submodel:
      (zut_cF <- fitmv(submodels=list(mod1=list(migStatus ~ 1),
                                      mod2=list(status2 ~ 1+ corrFamily(1|longitude+latitude))), # verbose=c(TRACE=TRUE), 
                       fixed=list(phi=c(0.02,0.02)), covStruct=list(corrFamily=MaternIMRFa(mesh=mesh, fixed=c(alpha=2))),
                       data=cap_mv))
      #
      # ... but this works without warning:
      zut_reg <- try(fitmv(submodels=list(mod1=list(migStatus ~ 1),
                                          mod2=list(status2 ~ 1+ MaternIMRFa(1|longitude+latitude, mesh=mesh, fixed=c(alpha=2)))), 
                           fixed=list(phi=c(0.02,0.02)), 
                           data=cap_mv))
      
      testthat::test_that("fitmv cF OK",
                          testthat::expect_true(diff(range(logLik(zut_IMRF), logLik(zut_cF), logLik(zut_reg)))<1e-5))
      
    }
  }
  
  {
    cat(crayon::yellow("ARp; ")) 
    {
      ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package
      AR1fit <- fitme(lh ~ 1 + AR1(1|time), data=ts, method="REML")
      ARpfit <- fitme(lh ~ 1 + ARp(1|time), data=ts, method="REML")
      p1 <- predict(AR1fit, newdata=ts[2:4,])
      p2 <- predict(ARpfit, newdata=ts[2:4,])
      testthat::expect_true(diff(range(p1-p2))<1e-6)
      
      AR2fit <- fitme(lh ~ 1 + ARp(1|time, p=2), data=ts, method="REML")
      LRT(AR2fit,ARpfit)
      
      # Need to force nonzero phi (=> uncertainty in fitted values) for meaningful preVar comparison:
      AR1fit <- fitme(lh ~ 1 + AR1(1|time), data=ts, method="REML", fixed=list(phi=0.1))
      ARpfit <- fitme(lh ~ 1 + ARp(1|time), data=ts, method="REML", fixed=list(phi=0.1))
      p1 <- get_predVar(AR1fit)[2:4]
      p2 <- get_predVar(ARpfit)[2:4]
      p3 <- get_predVar(AR1fit, newdata=ts[2:4,])
      p4 <- get_predVar(ARpfit, newdata=ts[2:4,])
      testthat::expect_true(diff(range(p1-p2,p1-p3,p1-p4))<1e-8)
    }
    cat(crayon::yellow("ARMA; ")) 
    {
      ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package
      ARMAfit <- fitme(lh ~ 1 + ARMA(1|time,p=1,q=1), data=ts, method="REML")
      p1 <- predict(ARMAfit)[2:4,]
      p2 <- predict(ARMAfit, newdata=ts[2:4,])
      testthat::expect_true(diff(range(p1-p2))<1e-6)
      p1 <- get_predVar(ARMAfit)[2:4]
      p2 <- get_predVar(ARMAfit, newdata=ts[2:4,])
      testthat::expect_true(diff(range(p1-p2))<1e-6)
      
      # Need to force nonzero phi (=> uncertainty in fitted values) for meaningful predVar comparison:
      AR1fit <- fitme(lh ~ 1 + AR1(1|time), data=ts, method="REML", fixed=list(phi=0.1))
      AR10fit <- fitme(lh ~ 1 + ARMA(1|time,p=1,q=1, fixed=c(q1=0)), data=ts, method="REML", fixed=list(phi=0.1))
      p1 <- get_predVar(AR1fit, newdata=ts[2:4,])
      p2 <- get_predVar(AR10fit, newdata=ts[2:4,])
      testthat::expect_true(diff(range(p1-p2))<1e-6)
    }
    { 
      
      cat(crayon::yellow("AR1 fitmv; ")) 
      {
        # good test-data because they are not ordered by time in the data.frame (min year= 1932)
        # But the rho estimate is effectively 1, so there are numerical problems
        load(file = paste0(projpath(),"/../tests_misc/misc/Sel_T.rda"))
        
        Phen_Sel <- droplevels(subset(Sel_T, Trait_Categ == 'Phenological'))
        
        rough_means <- with(Phen_Sel, tapply(Selection_mean, id, mean))
        
      }
      {
        (sub_AR1 <- fitme(Selection_mean ~ Year+AR1(1|Year), 
                          prior.weights = 1/Selection_SE^2,
                          upper=list(corrPars=list("1"=list(ARphi=0.999))),
                          data=Phen_Sel))
        (zut_AR1 <- fitmv(submodels=list(mod1=list(Selection_mean ~ 1),
                                         mod2=list(Selection_mean ~ Year+AR1(1|Year), 
                                                   prior.weights = quote(1/Selection_SE^2))), 
                          data=Phen_Sel))
        cat(crayon::blue("One expected warning as bc fitmv with unregistered corrFamily:"))
        (zut_cF <- fitmv(submodels=list(mod1=list(Selection_mean ~ 1),
                                        mod2=list(Selection_mean ~ Year+corrFamily(1|Year), 
                                                  prior.weights = quote(1/Selection_SE^2))), # verbose=c(TRACE=TRUE), 
                         covStruct=list(corrFamily=ARp()),
                         data=Phen_Sel))
        (zut_cF_reg <- fitmv(submodels=list(mod1=list(Selection_mean ~ 1),
                                            mod2=list(Selection_mean ~ Year+ARp(1|Year), 
                                                      prior.weights = quote(1/Selection_SE^2))), # verbose=c(TRACE=TRUE), 
                             data=Phen_Sel))
        try(testthat::test_that("fitmv cF AR1 OK",
                                testthat::expect_true(diff(range(logLik(zut_AR1), logLik(zut_cF), logLik(zut_cF_reg)))<1e-5)))
      }
      
    }
  }
  { #### Dyadic data
    {#### Simulate dyadic data
      
      set.seed(123)
      nind <- 10
      x <- runif(nind^2)  # but for easy comparison of predictiosn for reciprocal pairs, this is not used 
      id12 <- expand.grid(id1=seq(nind),id2=seq(nind))
      id1 <- id12$id1
      id2 <- id12$id2
      u <-  rnorm(50,mean = 0, sd=0.5)
      
      ## Same with non-additive individual effects:
      dist.u <- abs(u[id1] -  u[id2])
      z <- 0.1 + 1*x + dist.u + rnorm(nind^2,sd=0.2)
      
      dyaddf <- data.frame(x=x, z=z, id1=id1,id2=id2)
      
      dyaddf <- dyaddf[- seq.int(1L,nind^2,nind+1L),] 
    }
    { cat(crayon::yellow("diallel; "))
      
      (diallel_fit <- fitme(z ~1 +diallel(1|id1+id2), 
                            data=dyaddf)) 
      (diallel_p <- fitme(z ~1 +diallel(1|id1+id2), 
                          data=dyaddf[sample(nrow(dyaddf)),])) 
      (diallel_cF <- fitme(z ~1 +corrFamily(1|id1+id2), 
                           covStruct=list(corrFamily=diallel(tpar=0.42)), 
                           data=dyaddf)) 
      testthat::test_that("Check that diallel OK wrt permutations and OK as covStruct",
                          testthat::expect_true(diff(range(logLik(diallel_fit),logLik(diallel_p),logLik(diallel_cF)))<1e-12)) 
      recip_1_3 <- c(2L,19L) # reciprocal pairs 1:3 and 3:1
      (p1 <- predict(diallel_fit)[recip_1_3,])
      testthat::test_that("predict diallel  OK wrt reciprocal pairs",
                          testthat::expect_true(diff(p1)<1e-14))
      p2 <- predict(diallel_fit, newdata=diallel_fit$data[recip_1_3,])
      p3 <- predict(diallel_p, newdata=diallel_fit$data[recip_1_3,])
      p4 <- predict(diallel_cF, newdata=diallel_fit$data[recip_1_3,])
      testthat::test_that("predict diallel  OK wrt permutations and OK as covStruct",
                          testthat::expect_true(diff(range(p1-p2,p1-p3,p1-p4))<1e-8))
      (p1 <- get_predVar(diallel_fit)[recip_1_3])
      testthat::test_that("get_predVar diallel  OK wrt reciprocal pairs",
                          testthat::expect_true(diff(p1)<1e-14))
      p2 <- get_predVar(diallel_fit, newdata=diallel_fit$data[recip_1_3,])
      p3 <- get_predVar(diallel_p, newdata=diallel_fit$data[recip_1_3,])
      p4 <- get_predVar(diallel_cF, newdata=diallel_fit$data[recip_1_3,])
      testthat::test_that("get_predVar diallel  OK wrt permutations and OK as covStruct",
                          testthat::expect_true(diff(range(p1-p2,p1-p3,p1-p4))<1e-8))
    }
    
    { cat(crayon::yellow("ranGCA; ")) 
      
      (ranGCA_fit <- fitme(z ~1 +ranGCA(1|id1+id2), 
                           data=dyaddf)) 
      (ranGCA_p <- fitme(z ~1 +ranGCA(1|id1+id2), 
                         data=dyaddf[sample(nrow(dyaddf)),])) 
      (ranGCA_cF <- fitme(z ~1 +corrFamily(1|id1+id2), 
                          covStruct=list(corrFamily=ranGCA()), 
                          data=dyaddf)) 
      testthat::test_that("Check that ranGCA OK wrt permutations and OK as covStruct",
                          testthat::expect_true(diff(range(logLik(ranGCA_fit),logLik(ranGCA_p),logLik(ranGCA_cF)))<1e-12)) 
      recip_1_3 <- c(2L,19L) # reciprocal pairs 1:3 and 3:1
      (p1 <- predict(ranGCA_fit)[recip_1_3,])
      testthat::test_that("predict ranGCA  OK wrt reciprocal pairs",
                          testthat::expect_true(diff(p1)<1e-14))
      p2 <- predict(ranGCA_fit, newdata=ranGCA_fit$data[recip_1_3,])
      p3 <- predict(ranGCA_p, newdata=ranGCA_fit$data[recip_1_3,])
      p4 <- predict(ranGCA_cF, newdata=ranGCA_fit$data[recip_1_3,])
      testthat::test_that("predict ranGCA  OK wrt permutations and OK as covStruct",
                          testthat::expect_true(diff(range(p1-p2,p1-p3,p1-p4))<1e-8))
      (p1 <- get_predVar(ranGCA_fit)[recip_1_3])
      testthat::test_that("get_predVar ranGCA  OK wrt reciprocal pairs",
                          testthat::expect_true(diff(p1)<1e-14))
      p2 <- get_predVar(ranGCA_fit, newdata=ranGCA_fit$data[recip_1_3,])
      p3 <- get_predVar(ranGCA_p, newdata=ranGCA_fit$data[recip_1_3,])
      p4 <- get_predVar(ranGCA_cF, newdata=ranGCA_fit$data[recip_1_3,])
      testthat::test_that("get_predVar ranGCA  OK wrt permutations and OK as covStruct",
                          testthat::expect_true(diff(range(p1-p2,p1-p3,p1-p4))<1e-8))
    }
  }
}
  
