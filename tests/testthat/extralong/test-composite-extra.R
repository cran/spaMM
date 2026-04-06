cat(cli::col_yellow("\ntest-composite-extra:\n"))

cat(cli::col_yellow("checks AR1 composite"))

if ( ! exists("doSeeMe")) doSeeMe <- spaMM.getOption("doSeeMe") # in principle provided by parent tests file

{ # compare different algebras (default choice of method tested in test-determine_spprec.R)
  ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package
  
  (compAR1fitsp <-  fitme(lh ~ 1 + AR1(time|time), data=ts, control.HLfit=list(algebra="spprec"), 
                          lower=list(phi=1e-4),
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (compAR1fitspc <-  fitme(lh ~ 1 + AR1(time|time), data=ts, control.HLfit=list(algebra="spcorr"), 
                        lower=list(phi=1e-4),
                        fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (compAR1fitdec <-  fitme(lh ~ 1 + AR1(time|time), data=ts, control.HLfit=list(algebra="decorr"), 
                          lower=list(phi=1e-4),
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (crit <- diff(range(c(logLik(compAR1fitdec),logLik(compAR1fitspc),logLik(compAR1fitsp),-28.37433))))
  FIXME <- try(testthat::test_that(paste0("Whether the three algebras give consistent results for AR1(time|time): crit= ",signif(crit,4)," >1e-05"),
                                   testthat::expect_true(crit<1e-05) ), silent=TRUE)
  doSeeMe(FIXME) 
  cat(cli::col_yellow("Warning expected here:"))
  (p1 <- predict(compAR1fitsp))
  (p2 <- predict(compAR1fitspc, newdata=compAR1fitsp$data))
  (p3 <- predict(compAR1fitdec, newdata=compAR1fitsp$data))
  cat(cli::col_yellow("Warning expected here:"))
  (p4 <- predict(compAR1fitsp, newdata=compAR1fitsp$data)) 
  (crit <- diff(range(c(p1-p2,p1-p3,p1-p4))))
  FIXME <- testthat::test_that(paste0(
    "Whether the three algebras give consistent results\n for predict(<composite AR1>, newdata): crit= ",signif(crit,4)," >1e-05"),
                               testthat::expect_true(crit<1e-05) )
}

{ # Same with ARp(., p=1)
  
  cat(cli::col_yellow("; checks ARp composite"))
  ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package
  
  (compAR1fitsp <-  fitme(lh ~ 1 + ARp(time|time, p=1), data=ts, control.HLfit=list(algebra="spprec"), 
                          lower=list(phi=1e-4),
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (compAR1fitspc <-  fitme(lh ~ 1 + ARp(time|time, p=1), data=ts, control.HLfit=list(algebra="spcorr"), 
                           lower=list(phi=1e-4),
                           fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (compAR1fitdec <-  fitme(lh ~ 1 + ARp(time|time, p=1), data=ts, control.HLfit=list(algebra="decorr"), 
                           lower=list(phi=1e-4),
                           fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (crit <- diff(range(c(logLik(compAR1fitdec),logLik(compAR1fitspc),logLik(compAR1fitsp),-28.37433))))
  FIXME <- try(testthat::test_that(paste0("Whether the three algebras give consistent results for ARp(time|time, p=1): crit= ",signif(crit,4)," >1e-05"),
                                   testthat::expect_true(crit<1e-05) ), silent=TRUE)
  doSeeMe(FIXME) 
}

{ # Same with Matern
  cat(cli::col_yellow("; checks Matern composite"))
  ts <- data.frame(lh=lh,time=seq(48)) ## using 'lh' data from 'stats' package
  
  (compMatfitde <-  fitme(lh ~ 1 + Matern(time|time), data=ts, control.HLfit=list(algebra="decorr"), 
                          lower=list(phi=1e-4),
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  # Inscrutable numeric changes with Matrix devel v1.6.0 => premature stop, prevented by setting xtol_rel=1e-07 or 4e-6 but not 5e-6.
  (compMatfit <-  fitme(lh ~ 1 + Matern(time|time), data=ts, control.HLfit=list(algebra="spcorr"), 
                        lower=list(phi=1e-4),
                        fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (compMatfitsp <-  fitme(lh ~ 1 + Matern(time|time), data=ts, control.HLfit=list(algebra="spprec"), 
                          fixed=list(ranCoefs=list("1"=c(NA,0.1,NA)))))
  (crit <- diff(range(c(logLik(compMatfitde),logLik(compMatfit),logLik(compMatfitsp),-26.43478))))
  FIXME <- try(testthat::test_that(paste0("Whether the three algebras give consistent results for Matern(time|time): crit= ",signif(crit,4)," >1e-05"),
                                   testthat::expect_true(crit<1e-05) ), silent=TRUE)
  doSeeMe(FIXME) 
  
  (p1 <- predict(compMatfitsp))
  (p2 <- predict(compMatfit))
  (p3 <- predict(compMatfit, newdata=compMatfit$data))
  (p4 <- predict(compMatfitsp, newdata=compMatfit$data))
  (crit <- diff(range(c(p1-p2,p1-p3,p1-p4))))
  FIXME <- testthat::test_that(paste0(
    "Whether spcorr and spprec give consistent results\n for predict(<composite Matern>, newdata): crit= ",signif(crit,4)," >1e-05"),
    testthat::expect_true(crit<1e-05) )
  
}

cat(cli::col_yellow("more checks of AR1 and ARp composite"))
{
  data("Orthodont",package = "nlme")
  fix_rc <- list("1"=c(0.01,0.5,0.01))
  (fit1 <- fitme(distance ~ age + AR1(age|age), 
                 data = Orthodont,method="REML", 
                 control.HLfit=list(algebra="spcorr"), 
                 fixed=list(ARphi=0.1,ranCoefs=fix_rc)))
  (fit1sp <- fitme(distance ~ age + AR1(age|age), 
                   data = Orthodont,method="REML", 
                   control.HLfit=list(algebra="spprec"), 
                   fixed=list(ARphi=0.1,ranCoefs=fix_rc)))
  (fit1de <- fitme(distance ~ age + AR1(age|age), 
                   data = Orthodont,method="REML", 
                   control.HLfit=list(algebra="decorr"), 
                   fixed=list(ARphi=0.1,ranCoefs=fix_rc)))
  (fitp <- fitme(distance ~ age + ARp(age|age,p=1, fixed=c(p1=0.1)), 
                 data = Orthodont,method="REML", 
                 control.HLfit=list(algebra="spcorr"), 
                 fixed=list(ranCoefs=fix_rc)))
  (fitpsp <- fitme(distance ~ age + ARp(age|age,p=1, fixed=c(p1=0.1)), 
                   data = Orthodont,method="REML", 
                   control.HLfit=list(algebra="spprec"), 
                   fixed=list(ranCoefs=fix_rc)))
  (fitpde <- fitme(distance ~ age + ARp(age|age,p=1, fixed=c(p1=0.1)), 
                   data = Orthodont,method="REML", 
                   control.HLfit=list(algebra="decorr"), 
                   fixed=list(ranCoefs=fix_rc)))
  (crit <- diff(range(logLik(fit1),logLik(fit1sp),logLik(fit1de),logLik(fitp),logLik(fitpsp),logLik(fitpde))))
  testthat::test_that("composite AR1 & ARp OK", 
                      testthat::expect_true(crit<1e-12))
  p1 <- predict(fit1, newdata=fit1$data)
  pp <- predict(fitp, newdata=fit1$data) 
  p1s <- predict(fit1sp, newdata=fit1$data)
  pps <- predict(fitpsp, newdata=fit1$data) 
  p1d <- predict(fit1de, newdata=fit1$data)
  ppd <- predict(fitpde, newdata=fit1$data) 
  (crit <- diff(range(p1-pp,p1-p1s,p1-p1d,p1s-pps,p1d-ppd)))
  testthat::test_that("predict() composite AR1 & ARp OK", 
                      testthat::expect_true(crit<1e-12))
  p1 <- get_predVar(fit1, newdata=fit1$data)
  pp <- get_predVar(fitp, newdata=fit1$data) 
  p1s <- get_predVar(fit1sp, newdata=fit1$data)
  pps <- get_predVar(fitpsp, newdata=fit1$data) 
  p1d <- get_predVar(fit1de, newdata=fit1$data)
  ppd <- get_predVar(fitpde, newdata=fit1$data) 
  (crit <- diff(range(p1-p1s,p1-p1d, pp-pps,pp-ppd, p1-pp)))
  testthat::test_that("get_predVar() composite AR1 & ARp OK",
                     testthat::expect_true(crit<1e-6))
  #
  
}

{
  data("blackcap")
  # MLdistMat2 <- as.matrix(proxy::dist(blackcap[,c("latitude","longitude")]))
  MLcorMat2 <- MaternCorr(proxy::dist(blackcap[,c("latitude","longitude")]),
                          nu=0.6285603,rho=0.0544659)
  cap_mv <- blackcap
  cap_mv$name <- as.factor(rownames(blackcap))                
  cap_mv$grp <- 1L+(blackcap$migStatus>1)   
  ch <- chol(proxy::as.matrix(MLcorMat2,diag=1))
  set.seed(123)
  v1 <- tcrossprod(ch,t(rnorm(14,sd=1)))
  v2 <- tcrossprod(ch,t(rnorm(14,sd=1)))
  cap_mv$status <- 2*v1 + v2
  cap_mv$status2 <- 2*v1 - v2
  # plot(cap_mv)
}

if (TRUE) {
  cat(cli::col_yellow("; check predict mv() not composite first (problem pre-v3.8.34)"))
  (basic_rC1 <- fitmv(submodels=list(mod1=list(status ~ 1+ (mv(1,2)|name), fixed=list(phi=0.1)),
                                    mod2=list(status2 ~ 1+ (mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                     data=cap_mv))
  (basic_rC2 <- fitmv(submodels=list(mod1=list(status ~ 1+ (0+mv(1,2)|name), fixed=list(phi=0.1)),
                                    mod2=list(status2 ~ 1+ (0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                     data=cap_mv))
  (p11 <- predict(basic_rC1)[,1])
  (p12 <- predict(basic_rC1, newdata=basic_rC1$data)[,1] )
  (p21 <- predict(basic_rC2)[,1])
  (p22 <- predict(basic_rC2, newdata=basic_rC2$data)[,1] )
  crit <- diff(range(p11-p12,p11-p21,p21-p22))
  testthat::test_that(paste0("predict mv(): criterion was ",signif(crit,4)," >1e-7"),
                      testthat::expect_true(crit<1e-7) )
  
  # Preliminary independent fit test, pure ranCoefs
  zuta <- fitme(status ~ 1+ (1|name), #verbose=c(TRACE=TRUE),
                data=cap_mv, fixed=list(lambda=1,phi=0.1),
                control.HLfit=list(sparse_precision=FALSE))
  zutb <- fitme(status2 ~ 1+ (1|name), fixed=list(lambda=1,phi=0.1), #verbose=c(TRACE=TRUE),
                data=cap_mv, 
                control.HLfit=list(sparse_precision=FALSE))
  (zut0 <- fitmv(submodels=list(mod1=list(status ~ 1+ (0+mv(1,2)|name), fixed=list(phi=0.1)),
                                mod2=list(status2 ~ 1+ (0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                 data=cap_mv, 
                 fixed=list(ranCoefs=list("1"=c(1,0,1))),
                 control.HLfit=list(sparse_precision=FALSE)))
  logLik(zut0)-(logLik(zuta)+logLik(zutb))
  
}





#  OK independent fit test, trivial ranCoefs + corrMatrix
zuta <- fitme(status ~ 1+ corrMatrix(1|name), #verbose=c(TRACE=TRUE),
              data=cap_mv, corrMatrix=MLcorMat2, fixed=list(lambda=1,phi=0.1),
              control.HLfit=list(sparse_precision=FALSE))
zutb <- fitme(status2 ~ 1+ corrMatrix(1|name), fixed=list(lambda=1,phi=0.1), #verbose=c(TRACE=TRUE),
              data=cap_mv, corrMatrix=MLcorMat2, 
              control.HLfit=list(sparse_precision=FALSE))
cat(cli::col_yellow("; checks corrMatrix composite"))
(zut0d <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                              mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
               data=cap_mv, corrMatrix=MLcorMat2, 
               fixed=list(ranCoefs=list("1"=c(1,0,1))),
               control.HLfit=list(sparse_precision=FALSE)))
logLik(zut0d)-(logLik(zuta)+logLik(zutb))
(zut0s <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                              mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
               data=cap_mv, corrMatrix=MLcorMat2, 
               fixed=list(ranCoefs=list("1"=c(1,0,1))),
               control.HLfit=list(sparse_precision=TRUE)))
logLik(zut0s)-(logLik(zuta)+logLik(zutb)) 

# and now nontrivial ranCoefs:
(zut1d <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                              mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
               data=cap_mv, corrMatrix=MLcorMat2, 
               fixed=list(ranCoefs=list("1"=c(1,0.5,1))),
               control.HLfit=list(sparse_precision=FALSE)))
(zut1s <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                              mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
               data=cap_mv, corrMatrix=MLcorMat2, 
               fixed=list(ranCoefs=list("1"=c(1,0.5,1))), # meaning dependent on using (0+mv()).
               control.HLfit=list(sparse_precision=TRUE))) 
crit <- max(abs(range(predict(zut1d)[,1] - predict(zut1d,newdata=zut1d$data)[,1])))
testthat::test_that(paste0(
  "Whether predict(<composite corMatrix>, newdata) is OK: crit= ",signif(crit,4)," >1e-12"),
  testthat::expect_true(crit<1e-12) )
(crit <- max(abs(range(predict(zut1d,newdata=zut1d$data)[,1] - predict(zut1s,newdata=zut1d$data)[,1]))))
testthat::test_that(paste0(
  "Whether the two algebras give consistent results\n for predict(<composite corMatrix>, newdata) is OK: crit= ",signif(crit,4)," >1e-12"),
  testthat::expect_true(crit<1e-12) )

logLik(zut1s)-logLik(zut1d) 

ranef(zut1s, type="uncorrelated") # distinct chol factors, distinct u_h ...
ranef(zut1d, type="uncorrelated")
r1 <- ranef(zut1s)[[1]] # ... but the correlated ranefs are identical.
r2 <- ranef(zut1d)[[1]]
crit <- diff(range(r1-r2))
testthat::test_that(paste0("ranef corrMatrix(mv()...): criterion was ",signif(crit,4)," >1e-7"),
                    testthat::expect_true(crit<1e-7) )

{ # and now variable ranCoefs:
  { # correlation fixed to 0: => predVar seems OK 
    (zut1d <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                                   mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                    data=cap_mv, corrMatrix=MLcorMat2, 
                    fixed=list(ranCoefs=list("1"=c(NA,0,NA))),
                    control.HLfit=list(sparse_precision=FALSE)))
    (zut1s <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                                   mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                    data=cap_mv, corrMatrix=MLcorMat2, 
                    fixed=list(ranCoefs=list("1"=c(NA,0,NA))), # meaning dependent on using (0+mv()).
                    control.HLfit=list(sparse_precision=TRUE))) 
    (pVd <- get_predVar(zut1d))
    (pVs <- get_predVar(zut1s))
    (crit <- diff(range(pVd-pVs)))
    testthat::test_that(paste0("ranef corrMatrix(mv()...): criterion was ",signif(crit,4)," >5e-7"),
                        testthat::expect_true(crit<5e-7) ) 
  }
  if (FALSE) { # correlation fixed to 0.5: => predVar not OK  
    (zut1d <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                                   mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                    data=cap_mv, corrMatrix=MLcorMat2, 
                    fixed=list(ranCoefs=list("1"=c(NA,0.5,NA))),
                    control.HLfit=list(sparse_precision=FALSE)))
    (zut1s <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                                   mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                    data=cap_mv, corrMatrix=MLcorMat2, 
                    fixed=list(ranCoefs=list("1"=c(NA,0.5,NA))), # meaning dependent on using (0+mv()).
                    control.HLfit=list(sparse_precision=TRUE))) 
    (pVd <- get_predVar(zut1d))
    (pVs <- get_predVar(zut1s))
    (crit <- diff(range(pVd-pVs)))
    # test fails:
    testthat::test_that(paste0("ranef corrMatrix(mv()...): criterion was ",signif(crit,4)," >5e-7"),
                        testthat::expect_true(crit<5e-7) ) 
    
  }
  (zut1d <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                                 mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                  data=cap_mv, corrMatrix=MLcorMat2, 
                  control.HLfit=list(sparse_precision=FALSE)))
  (zut1s <- fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1)),
                                 mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|name), fixed=list(phi=0.1))), #verbose=c(TRACE=TRUE),
                  data=cap_mv, corrMatrix=MLcorMat2, 
                  control.HLfit=list(sparse_precision=TRUE))) 
  logLik(zut1s)-logLik(zut1d) 
  (p1d <- predict(zut1d))
  (p2d <- predict(zut1d, newdata=zut1d$data))
  (pVd <- get_predVar(zut1d))
  (pVd <- get_predVar(zut1d, newdata=zut1d$data,variances=list(cov=F)))
  # cat(cli::col_yellow("Warning expected here:"))
  (p1s <- predict(zut1s))
  cat(cli::col_yellow("Warning expected here:"))
  (p2s <- predict(zut1s, newdata=zut1s$data)) 
  (pVs <- get_predVar(zut1s))
  (pVs <- get_predVar(zut1s, newdata=zut1s$data,variances=list(cov=F)))
  (crit <- diff(range(p1d-p2d,p1d-p1s,p1s-p2s)))
  testthat::test_that(paste0("ranef corrMatrix(mv()...): criterion was ",signif(crit,4)," >5e-7"),
                      testthat::expect_true(crit<5e-7) )
  if (FALSE) {
    # The test fails because I don't have code to handle uncertainty in ranCoefs' correlation parameters
    # and the returned result depend on the Cholesky representation.
    (crit <- diff(range(pVd-pVs)))
    testthat::test_that(paste0("ranef corrMatrix(mv()...): criterion was ",signif(crit,4)," >5e-7"),
                        testthat::expect_true(crit<5e-7) ) 
  }
}

{  cat(cli::col_yellow( "Matern(LHS |<nested RHS>) and more corrMatrix(LHS |<nested RHS>)... " ))
  {
    data("blackcap")
    toy <- blackcap
    toy$ID <- gl(7,2)
    grp <- rep(1:2,7)
    toy$migStatus <- toy$migStatus +(grp==2)
    toy$loc <- rownames(toy) # to use as levels matching the corrMatrix dimnames
    
    toy$grp <- factor(grp)
    toy$bool <- toy$grp==1L
    toy$boolfac <- factor(toy$bool)
    toy$num <- seq(from=1, to=2, length.out=14)
    
    ## Build a toy corrMatrix as perturbation of identity matrix:
    n_rhs <- 14L
    eps <- 0.1
    set.seed(123)
    rcov <- ((1-eps)*diag(n_rhs)+eps*rWishart(1,n_rhs,diag(n_rhs)/n_rhs)[,,1])
    # eigen(rcov)$values
    colnames(rcov) <- rownames(rcov) <- toy$loc # DON'T FORGET NAMES
  }
  
  {
    ## The unconstrained fit fits to correlation = -1 and vanishing phi by augmented-y algo, both of which complicates comparisons
    # (fit1yaug <- fitme(migStatus ~ grp + Matern(grp|latitude %in% grp), data=toy, fixed=list(ranCoefs=list("1"=c(NA,-0.5,NA)))))
    ## => Force non y-augmented algo and then consistently fitted phi across calls:
    (fit1 <- fitme(migStatus ~ grp + Matern(grp|latitude %in% grp), data=toy, init=list(phi=1e-5), fixed=list(ranCoefs=list("1"=c(NA,-0.5,NA)))))
    
    { # to test for handling of levels found only once in the data (and next, predict() checks found more issues):
      ttoy <- toy[c(1,2*seq(7)),]
      fit1t <- fitme(migStatus ~ grp + Matern(grp|latitude %in% grp), data=ttoy, init=list(phi=1e-5), fixed=list())
      predict(fit1t) 
      predict(fit1t, newdata=ttoy) 
      predict(fit1t, newdata=ttoy[rev(seq(8)),])[rev(seq(8)),] 
    }
    
    rC <- get_ranPars(fit1)$ranCoefs[[1]]
    corMat <- Corr(fit1)[[1]][1:14,1:14]/rC[1] 
    colnames(corMat) <- rownames(corMat) <- toy$loc 
    # 
    (fit1fix_nest <- fitme(migStatus ~ grp + corrMatrix(grp|loc %in% grp), data=toy, init=list(phi=1e-5), corrMatrix=corMat, fixed=list(ranCoefs=list("1"=rC))))
    (fit1fix <- fitme(migStatus ~ grp + corrMatrix(grp|loc), data=toy, init=list(phi=1e-5), corrMatrix=corMat, fixed=list(ranCoefs=list("1"=rC)))) 
    
    testthat::expect_true(diff(range(logLik(fit1),logLik(fit1fix_nest),logLik(fit1fix)))<3e-09) 
    
  }
  
}

if (TRUE) { # 
  cat(cli::col_yellow("; checks IMRF composite"))
  
  # test by compar with example in doc for IMRF
  data("blackcap") ## toy examples; but IMRF may be useful only for much larger datasets
  data("small_spde") ## load object of class 'inla.spde2', created and saved by :
  {
    (fit_SPDE <- fitme(migStatus ~ 1 + IMRF(1|longitude+latitude, model=small_spde), 
                      data=blackcap))
    # it seems important to rescale 'means'...
    (fit_SPDE_fixrc <- fitme(migStatus ~ 1 + IMRF(1+I(means/100)|longitude+latitude, model=small_spde), 
                          data=blackcap, fixed=list(ranCoefs=list("1"=c(0.08846,0,0.00001)))))
    (fit_SPDE_rc <- fitme(migStatus ~ 1 + IMRF(1+I(means/100)|longitude+latitude, model=small_spde), 
                          data=blackcap, fixed=list(ranCoefs=list("1"=c(NA,0,0.00001)))))
    
    (crit <- diff(range(c(logLik(fit_SPDE_fixrc),logLik(fit_SPDE_rc)))))
    testthat::expect_true(crit<1e-8)
    # also close to fit_SPDE since we forced the variance of the slope to 1e-5
    
    if ("INLA" %in% installed.packages()) {
      # same comparison for MaternIMRFa with fixed alpha=2 ~ default IMRF
      spd <- sp::SpatialPointsDataFrame(coords = blackcap[, c("longitude", "latitude")],
                                        data = blackcap)
      small_mesh <- fmesher::fm_mesh_2d_inla(loc = fmesher::fm_mesh_2d_map(sp::coordinates(spd)),
                                       max.n=100, # only for demonstration purposes
                                       max.edge = c(3, 20))
      (fit_SPDE_cF <- fitme(migStatus ~ 1 + MaternIMRFa(1|longitude+latitude, 
                                                        mesh=small_mesh,
                                                        fixed=c(alpha=2)), 
                            data=blackcap, 
                            fixed=list(phi=1e-6))) # note 'failed' optimisation without this
      (crit <- diff(range(c(logLik(fit_SPDE_cF),logLik(fit_SPDE)))))
      testthat::expect_true(crit<1e-5) # just
      
      (fit_SPDE_rc_cF <- fitme(migStatus ~ 1 + MaternIMRFa(1+I(means/100)|longitude+latitude, 
                                                           mesh=small_mesh,
                                                           fixed=c(alpha=2)), 
                               data=blackcap, 
                               fixed=list(ranCoefs=list("1"=c(NA,0,0.00001)),
                                          phi=1e-6))) 
      # (crit <- diff(range(c(logLik(fit_SPDE_cF),logLik(fit_SPDE_rc_cF)))))
      # testthat::expect_true(crit<3e-5) # not quite identical but does have to
      #    since forcing the variance of the slope to 1e-5 only approximates the simpler model.
    }
    
  }
}


####################################################################"


if (FALSE) { # that was devel scratch, not tidy tests
  ZUT1 <- function(fit1, caveat="assuming lambda=1") {
    Lmat <- t(chol(proxy::as.matrix(MLcorMat2,diag=1))) # 
    
    tcrossprod(get_ZALMatrix(fit1,force_bind = TRUE)) - tcrossprod(Lmat) 
    
    I0ZX <- rbind(cbind(diag(x=14L),0*fit1$X.pv),
                  cbind(Lmat,fit1$X.pv) )
    w <- diag(sqrt(c(rep(1,14),rep(1/fit$phi[[1]],14))))
    
    wI0ZX <-w %*% I0ZX
    
    w_y <- w %*% c(rep(0,14),cap_mv$status)            
    
    Matrix::solve(crossprod(wI0ZX), crossprod(wI0ZX,w_y))[,1]
    
  }
  
  ZUT <- function(fit, caveat="assuming lambda=1") {
    if (.is_spprec_fit(fit)) {
      if (FALSE) {
        Lmat <- t(Matrix::kronecker(chol(attr(fit$strucList[[1]],"latentL_blob")$compactcovmat),
                                    t(as(Cholesky(as(proxy::as.matrix(MLcorMat2,diag=1),"sparseMatrix")),"CsparseMatrix")))) # 
        Lmat <- t(Matrix::kronecker(chol(attr(fit$strucList[[1]],"latentL_blob")$compactcovmat),
                                    (solve(as(Cholesky(as(solve(proxy::as.matrix(MLcorMat2,diag=1)),"sparseMatrix")),"CsparseMatrix"))))) # 
      }
      Lmat <- t(solve(BLOB$chol_Q)) # 
      if (FALSE) {
        # mais alors Md2hdv2 should be 
        crossprod(Lmat)/fit$phi[[1]] +diag(28)
        # la version "chol_Q" de Lmat marche
        Lmat <- t(solve(BLOB$chol_Q)) # 
        crossprod(Lmat)/fit$phi[[1]] +diag(28) - solve(BLOB$chol_Q) %*% (BLOB$Gmat) %*% solve(t(BLOB$chol_Q))
        get_ZALMatrix(fit,force_bind = TRUE) - Lmat 
        solve(crossprod(BLOB$factor_inv_Md2hdv2)) - crossprod(Lmat)/fit$phi[[1]] +diag(28)
      }
    } else {
      Lmat <- t(kronecker(chol(attr(fit$strucList[[1]],"latentL_blob")$compactcovmat),chol(proxy::as.matrix(MLcorMat2,diag=1)))) # 
    }
    if (FALSE) {
      # BLOB$logdet_sqrt_d2hdv2 is 22.35205
      (logdet <- Matrix::determinant(crossprod(Lmat)/fit$phi[[1]] +diag(28), logarithm = TRUE)$modulus/2)
      # - get_from_MME(sXaug,"logdet_sqrt_d2hdv2") + n_u_h*log(2*pi)/2 est 3.378226
      # bref la logLik est
      fit$APHLs$clik + sum(dnorm(fit$v_h,log = TRUE))
      
      
      as.matrix(tcrossprod(get_ZALMatrix(fit,force_bind = TRUE))) - tcrossprod(Lmat) 
      
      
      # sauf que les crossprod ne sont pas unique même si les tcrossprod sont identiques
      solve(BLOB$chol_Q) %*% BLOB$chol_Q %*% (crossprod(Lmat)/fit$phi[[1]] +diag(28)) %*% t(BLOB$chol_Q) %*% solve(t(BLOB$chol_Q))
      solve(BLOB$chol_Q) %*% ( BLOB$chol_Q %*% (crossprod(Lmat)/fit$phi[[1]]) %*% t(BLOB$chol_Q)+
                                 BLOB$chol_Q %*% (diag(28)) %*% t(BLOB$chol_Q)) %*% solve(t(BLOB$chol_Q))
      solve(BLOB$chol_Q) %*% ( BLOB$chol_Q %*% (crossprod(Lmat)/fit$phi[[1]]) %*% t(BLOB$chol_Q)+
                                 precisionMatrix) %*% solve(t(BLOB$chol_Q))
      
      
      solve(BLOB$chol_Q) %*% (BLOB$Gmat) %*% solve(t(BLOB$chol_Q))
    }
    
    I0ZX <- rbind(cbind(diag(x=28L),0*fit$X.pv),
                  cbind(Lmat,fit$X.pv) )
    w <- diag(sqrt(c(rep(1,28),rep(1/fit$phi[[1]],28))))
    
    wI0ZX <-w %*% I0ZX
    
    w_y <- w %*% c(rep(0,28),cap_mv$status,cap_mv$status2)            
    
    Matrix::solve(crossprod(wI0ZX), crossprod(wI0ZX,w_y))[,1]
    
  }
  ZUT(zut0d)
  ZUT(zut1s)
  
}

