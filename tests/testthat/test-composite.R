cat(crayon::yellow("\ntest-composite.R: "))

# more tests of composite in 
# "tests_private/VarCorr.R" (numeric LHS); 
# "nested/test-composite-nested.R" (mv()); 

if (FALSE) { # __F I X M E__ Interesting alternative numerical setings: only small effects on numerical precision, and some effects on speed (visibly on spherical fit)
  spaMM.options(
    nloptr=list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-5, print_level=0), # distinct xtol_rel
    xtol_abs_factors=c(rcLam=5e-6,rcCor=5e-5,others=5e-11,abs=1e-7), # distinct rcLam and rcCor
    doSeeMe=warning
  )
}

doSeeMe <- spaMM.getOption("doSeeMe")

if (spaMM.getOption("example_maxtime")>8) {
  { # testthat's on examples for composite-ranef. 
    
    ## Data preparation
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
    eigen(rcov)$values
    colnames(rcov) <- rownames(rcov) <- toy$loc
    
    ##### Illustrating the different LHS types
    
    ### <LHS> is logical (TRUE/FALSE) => No induced random-coefficient C matrix; 
    #   corrMatrix affects only responses for which <LHS> is TRUE:
    #
    (fit1 <- fitme(migStatus ~ bool + corrMatrix(bool|loc), data=toy, corrMatrix=rcov))
    #
    # Matrix::image(get_ZALMatrix(fit))
    predict(fit1)
    predict(fit1, newdata=fit1$data[14:1,])[14:1,]
    get_predVar(fit1)
    get_predVar(fit1, newdata=fit1$data[14:1,])[14:1]
    
    hatvalues(fit1)
    
    ### <RHS> is a factor built from a logical => same a 'logical' case above:
    #
    (fit2 <- fitme(migStatus ~ boolfac + corrMatrix(boolfac|loc), data=toy, corrMatrix=rcov))
    (crit <- diff(range(c(logLik(fit1),logLik(fit2)))))
    testthat::test_that(paste0("bool vs boolfac LHS: criterion was ",signif(crit,4)," >1e-12"),
                        testthat::expect_true(crit<1e-12) )
    #
    # Matrix::image(get_ZALMatrix(fit))
    
    
    ### <RHS> is a factor not built from a logical: 
    # (grp|.) and (0+grp|.) lead to equivalent fits of the same composite model, 
    #   but contrasts are not used in the second case and the C matrices differ,
    #   as for standard random-coefficient models.
    #
    (fit1 <- fitme(migStatus ~ grp +  corrMatrix(grp|loc), data=toy, corrMatrix=rcov))
    (fit2 <- fitme(migStatus ~ grp +  corrMatrix(0+grp|loc), data=toy, corrMatrix=rcov))
    (crit <- diff(range(c(logLik(fit1),logLik(fit2),-17.07165))))
    testthat::test_that(paste0("grp vs 0+grp LHS: criterion was ",signif(crit,4)," >1e-4"),
                        testthat::expect_true(crit<1e-4) )
    # 
    # => same fits, but different internal structures:
    Matrix::image(fit1$ZAlist[[1]]) # (contrasts used) 
    Matrix::image(fit2$ZAlist[[1]]) # (contrasts not used)
    # Also compare ranef(fit1) versus ranef(fit2) 
    #
    #
    ## One can fix the C matrix, as for standard random-coefficient terms 
    #
    (fit1 <- fitme(migStatus ~ grp +  corrMatrix(0+grp|loc),data=toy, corrMatrix=rcov, 
                   fixed=list(ranCoefs=list("1"=c(1,0.5,1)))))
    #       
    # same result without contrasts hence different 'ranCoefs':             
    #
    (fit2 <- fitme(migStatus ~ grp +  corrMatrix(grp|loc), data=toy, corrMatrix=rcov, 
                   fixed=list(ranCoefs=list("1"=c(1,-0.5,1)))))
    (crit <- diff(range(c(logLik(fit1),logLik(fit2)))))
    testthat::test_that(paste0("grp vs 0+grp LHS, fixed ranCoefs: criterion was ",signif(crit,4)," >1e-12"),
                        testthat::expect_true(crit<1e-12) )
    
    
    ### <LHS> is numeric (but not '0+numeric'):
    # composite model with C being 2*2 for Intercept and numeric variable
    #
    (fitme(migStatus ~ num +  corrMatrix(num|loc), data=toy, corrMatrix=rcov))
    
    ### <LHS> is 0+numeric: no random-coefficient C matrix 
    #  as the Intercept is removed, but the correlated random effects 
    #  arising from the corrMatrix are multiplied by sqrt(<numeric variable>)
    #
    (fitme(migStatus ~ num +  corrMatrix(0+num|loc), data=toy, corrMatrix=rcov))
  }
  
  if (FALSE) { # next block has similar tests, and more
    data("sleepstudy", package="lme4")
    set.seed(123)
    #rcov18 <- (17*diag(18)+rWishart(1,18,diag(18)/18)[,,1])/18
    rcov18 <-  (diag(18)+rWishart(1,18,diag(18)/18)[,,1])/2
    colnames(rcov18) <- rownames(rcov18) <- levels(sleepstudy$Subject)
    #eigen(rcov)$values
    # OK even for clearly non-unit rcov:
    (res1 <- fitme(Reaction ~ Days + corrMatrix(Days|Subject), data = sleepstudy, corrMatrix=rcov18))
    (res2 <- fitme(Reaction ~ Days + corrMatrix(Days|Subject), data = sleepstudy, corrMatrix=rcov18, control.HLfit=list(sparse_precision=TRUE)))
    (res3 <- fitme(Reaction ~ Days + corrMatrix(Days|Subject), data = sleepstudy, corrMatrix=rcov18, control=list(refit=TRUE),
                   control.HLfit=list(sparse_precision=TRUE)))
    (res4 <- fitme(Reaction ~ Days + corrMatrix(Days|Subject), data = sleepstudy, corrMatrix=rcov18, control=list(refit=TRUE),
                   control.HLfit=list(sparse_precision=F)))
    (crit <- diff(range(c(logLik(res1),logLik(res2),logLik(res3),logLik(res4)))))
    FIXME <- testthat::test_that(paste0("sleepstudy spprec F/T and refit F/T: criterion was ",signif(crit,4)," >1e-09"),
                        testthat::expect_true(crit<1e-09) )
    if ( ! FIXME) doSeeMe("Do see me!") 
  }
  
  { # spprec T/F, refit T/F; two corrMatrix terms vs ranCoefs;
    requireNamespace("nlme")
    data("Orthodont",package = "nlme")
    (fit <- fitme(distance ~ age+(age|Subject), data = Orthodont))
    n_rhs <- 27L
    set.seed(123)
    rcov27 <- ((n_rhs-1L)*diag(n_rhs)+rWishart(1,n_rhs,diag(n_rhs)/n_rhs)[,,1])/n_rhs
    #eigen(rcov27)$values
    colnames(rcov27) <- rownames(rcov27) <- levels(Orthodont$Subject)
    (fit1 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=rcov27))
    (fit2 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=rcov27, control=list(refit=list(ranCoefs=TRUE))))
    #trace(spaMM:::.makeCovEst1,print=TRUE)
    (fit3 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=rcov27, control.HLfit=list(sparse_precision=TRUE)))
    (fit4 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=rcov27, control.HLfit=list(sparse_precision=TRUE),
                   control=list(refit=list(ranCoefs=TRUE))))
    (crit <- diff(range(c(logLik(fit1),logLik(fit2),logLik(fit3),logLik(fit4)))))
    FIXME <- testthat::test_that(paste0("Orthodont spprec F/T and refit F/T: criterion was ",signif(crit,4)," >1e-09"),
                        testthat::expect_true(crit<1e-09) ) # was <1e09 until I changed as_precision to use chol2inv(chol()). Back to excellent precision much later
    if ( ! FIXME) doSeeMe("Do see me!") 
    if (FALSE) {
      trace(spaMM:::.makeCovEst1,print=TRUE)
      (fit <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=rcov27, control.HLfit=list(sparse_precision=TRUE),
                    control=list(refit=list(ranCoefs=TRUE, phi=TRUE))))
      untrace(spaMM:::.makeCovEst1)
    } # : has had a periodic behavior, but appears to work now (with 27 .makeCovEst1 calls)
    
    (fit1 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=rcov27, 
                   fixed=list(phi=1.7,ranCoefs=list("1"=c(NA,-0,NA)))))
    (fit2 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=rcov27, 
                   fixed=list(phi=1.7,ranCoefs=list("1"=c(NA,-0,NA))),control.HLfit=list(sparse_precision=TRUE)))
    (fit3 <- fitme(distance ~ age+ corrMatrix(1|Subject)+corrMatrix(0+age|Subject), data = Orthodont, 
                   covStruct=list(corrMatrix=rcov27,corrMatrix=rcov27), 
                   fixed=list(phi=1.7),control.HLfit=list(sparse_precision=F))) 
    (fit4 <- fitme(distance ~ age+ corrMatrix(1|Subject)+corrMatrix(0+age|Subject), data = Orthodont, 
                   covStruct=list(corrMatrix=rcov27,corrMatrix=rcov27), 
                   fixed=list(phi=1.7),control.HLfit=list(sparse_precision=T))) 
    (crit <- diff(range(c(logLik(fit1),logLik(fit2),logLik(fit3),logLik(fit4)))))
    FIXME <- testthat::test_that(paste0("two corrMatrix terms vs ranCoefs: spprec F/T: criterion was ",signif(crit,4)," >1e-07"),
                        testthat::expect_true(crit<1e-07) )
    if ( ! FIXME) doSeeMe("Do see me!") 
    predict(fit1)[1:6]
    predict(fit1, newdata=fit1$data[1:6,])
    predict(fit1, newdata=fit1$data[6:1,])[6:1,]  ## OK
    
    predict(fit4)[1:6]
    predict(fit4, newdata=fit1$data)[1:6,]
    predict(fit4, newdata=fit1$data[6:1,])[6:1,]
    
    wh <- "predVar"
    dsp <- TRUE 
    (p1 <- get_predVar(fit1,which=wh, variances=list(disp=dsp))[1:6])
    p1n <- get_predVar(fit1,which=wh, variances=list(disp=dsp), newdata=fit1$data[1:6,])[1:6]
    p2 <- get_predVar(fit2,which=wh, variances=list(disp=dsp))[1:6] ## OK
    p2n <- get_predVar(fit2,which=wh, variances=list(disp=dsp), newdata=fit1$data[1:6,])[1:6]
    (crit <- diff(range(c(p1-p1n,p1-p2,p1-p2n))))
    FIXME <- testthat::test_that(paste0("predVar spprec F/T: criterion was ",signif(crit,4)," >1e-07"),
                        testthat::expect_true(crit<1e-07) )
    if ( ! FIXME) doSeeMe("Do see me!") 
    
    (p3 <- get_predVar(fit3,which=wh, variances=list(disp=dsp))[1:6]) ## 3 vs 4 OK
    p3n <- get_predVar(fit3,which=wh, variances=list(disp=dsp), newdata=fit1$data[1:6,])[1:6]
    p4 <- get_predVar(fit4,which=wh, variances=list(disp=dsp))[1:6]
    p4n <- get_predVar(fit4,which=wh, variances=list(disp=dsp), newdata=fit1$data[1:6,])[1:6]
    (crit <- diff(range(c(p3-p3n,p3-p4,p4-p4n))))
    FIXME <- testthat::test_that(paste0("not composite predVar spprec F/T: criterion was ",signif(crit,4)," >1e-08"),
                        testthat::expect_true(crit<1e-08) )
    if ( ! FIXME) doSeeMe("Do see me!") 

    (crit <- diff(range(c(p1-p3))))
    FIXME <- testthat::test_that(paste0("two corrMatrix terms vs ranCoefs (Important test .calc_logdisp_cov() code for ranCoefs!)",
                                        signif(crit,4)," >1e-05"),
                                 testthat::expect_true(crit<1e-05) )
    if ( ! FIXME) doSeeMe("Do see me!") 
    
    
    # clear predVar difference between fixing the corr or not in orrMatrix(age|Subject).
    get_predVar(fit1 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=rcov27, 
                              fixed=list(phi=1.7,ranCoefs=list("1"=c(NA,0.57546616,NA)))))[1:6]
    get_predVar(fit1 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=rcov27, 
                              fixed=list(phi=1.7,ranCoefs=list())))[1:6]

    if (FALSE) {
      # The difference disappears when fixing the lambdas:
      (fit1 <- fitme(distance ~ age+(age|Subject), data = Orthodont, #corrMatrix=rcov27, 
                     fixed=list(ranCoefs=list("1"=c(2,-0,2)))))
      (fit2 <- fitme(distance ~ age+(age|Subject), data = Orthodont, #corrMatrix=rcov27, 
                     fixed=list(ranCoefs=list("1"=c(2,-0,2))),control.HLfit=list(sparse_precision=TRUE)))
      (fit3 <- fitme(distance ~ age+ (1|Subject)+(0+age|Subject), data = Orthodont, 
                     #covStruct=list(corrMatrix=rcov27,corrMatrix=rcov27), 
                     fixed=list(lambda=c(2,2)),control.HLfit=list(sparse_precision=F))) 
      (fit4 <- fitme(distance ~ age+ (1|Subject)+(0+age|Subject), data = Orthodont, 
                     #covStruct=list(corrMatrix=rcov27,corrMatrix=rcov27), 
                     fixed=list( lambda=c(2,2)),control.HLfit=list(sparse_precision=T)))  
    }
    
    
    diag27 <- diag(27L)
    colnames(diag27) <- rownames(diag27) <- colnames(rcov27)
    (fit1 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=diag27, 
                   fixed=list(ranCoefs=list("1"=c(NA,0,NA)))))
    (fit2 <- fitme(distance ~ age+corrMatrix(age|Subject), data = Orthodont, corrMatrix=diag27, 
                   fixed=list(ranCoefs=list("1"=c(NA,-0,NA))),control.HLfit=list(sparse_precision=TRUE)))
    (fit3 <- fitme(distance ~ age+(age||Subject), data = Orthodont, 
                   fixed=list()))
    (fit4 <- fitme(distance ~ age+(age||Subject), data = Orthodont, 
                   fixed=list(),control.HLfit=list(sparse_precision=TRUE)))
    (crit <- diff(range(c(logLik(fit1),logLik(fit2),logLik(fit3),logLik(fit4)))))
    testthat::test_that(paste0("two uncorr corrMatrix() terms by ranCoefs: spprec F/T: criterion was ",signif(crit,4)," >5e-10"),
                        testthat::expect_true(crit<5e-10) )
    
    
  }
} 
