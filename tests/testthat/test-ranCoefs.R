# spaMM.options(example_maxtime=70) # keep it handy...
if (file.exists((privdata <- paste0(projpath(),"/../tests_misc/ranCoefs/all_fitness.txt")))) {
  cat(crayon::yellow("\ntest ranCoefs:"))
  my.data <- read.table(privdata, header = TRUE, sep = "\t",dec = ".")
  my.data$line <- factor(as.character(my.data$line))
  my.data <- na.omit(my.data)
  if (TRUE) {
    (fitme3 <- fitme(total_red ~ sex*env + (1|rep) + (0 + env|line),
                     data = my.data, method="ML"))
    how(fitme3) # 0.7s v3.0.34
    # final precision depends on two steps ! And ultimately on HLfit/hatval precision in the refit
    testval3o <- 1559.81264462379 ## sph gives practically the same result but longer 
    zut <- try(testthat::expect_equal((res3o <- attr(attr(fitme3,"optimInfo")$optim.pars,"optr")$objective),
                                      testval3o,tol=1e-8),silent=TRUE)
    if (inherits(zut,"try-error")) {
      if (res3o>testval3o) strg <- " poorer" else strg <- " better" 
      cat(paste0("objective from fitme3's optimInfo=",
                 res3o, strg, " than tested value ",testval3o," by ",res3o-testval3o ,"\n"))
    }
    testval3f <- -1559.812644545275
    zut <- try(testthat::expect_equal((res3f <- logLik(fitme3)[[1L]]),testval3f,tol=1e-8),silent=TRUE) 
    if (inherits(zut,"try-error")) {
      if (res3f<testval3f) strg <- " poorer" else strg <- " better" 
      cat(paste0("logLik(fitme3)=",
                 res3f, strg, " than tested value ",testval3f," by ",res3f-testval3f ,"\n"))
    }
    # check effect of changing contrasts for factors between fit and predict: 
    p1 <- predict(fitme3)[2]
    # options(contrasts=c("contr.treatment", 'contr.poly'))
    p2 <- predict(fitme3,newdata=my.data[2,])[1]
    oldopt <- options(contrasts=c("contr.poly", 'contr.treatment'))
    p3 <- predict(fitme3,newdata=my.data[2,])[1]
    options(oldopt)
    testthat::expect_true(diff(range(c(p1,p2,p3)))<1e-8) 
    #
  }
  if (spaMM.getOption("example_maxtime")>69) {
    if (file.exists((privtest <- paste0(projpath(),"/package/tests_private/test-rc_transf.R")))) {
      source(privtest)
    }
    (HLfit3 <- HLfit(total_red ~ sex*env + (1|rep) + (0 + env|line),
                     data = my.data, HLmethod="ML"))
    how(HLfit3) # 2.2s v3.0.34
    testval3h <- -1559.81264432437 # chol -- better than ever
    zut <- try(testthat::expect_equal((res3h <- logLik(HLfit3)[[1L]]),testval3h,tol=1e-8),silent=TRUE) 
    if (inherits(zut,"try-error")) {
      if (res3h<testval3h) strg <- " poorer" else strg <- " better" 
      cat(paste0("logLik(HLfit3)=",
                 res3h, strg, " than tested value ",testval3h," by ",res3h-testval3h ,"\n"))
    }
    if (TRUE) { # a check of get_inits_from_fit() usage
      reinit <- get_inits_from_fit(HLfit3)
      abyss <- fitme(total_red ~ sex*env + (1|rep) + (0 + env|line),
                     data = my.data, fixed=reinit$init.HLfit["ranCoefs"])
    }
    (fitme6 <- fitme(total_red ~ sex*env + (1|rep) + (0 + sex:env|line),
                     data = my.data, method="ML"))  
    how(fitme6) # 2.5s v3.0.34
    testval6o <- 1536.080805076491 # chol -- good 
    zut <- try(testthat::expect_equal((res6o <- attr(attr(fitme6,"optimInfo")$optim.pars,"optr")$objective), 
                                      testval6o,tol=1e-8),silent=TRUE)
    if (inherits(zut,"try-error")) {
      if (res6o>testval6o) strg <- " poorer" else strg <- " better" 
      cat(paste0("objective from fitme6's optimInfo=",
                 res6o, strg, " than tested value ",testval6o," by ",res6o-testval6o ,"\n"))
    }
    testval6f <- -1536.080805075333 # chol -- good 
    zut <- try(testthat::expect_equal((res6f <- logLik(fitme6)[[1L]]),testval6f,tol=1e-8),silent=TRUE) 
    if (inherits(zut,"try-error")) {
      if (res6f<testval6f) strg <- " poorer" else strg <- " better" 
      cat(paste0("logLik(fitme6)=",
                 res6f, strg, " than tested value ",testval6f," by ",res6f-testval6f ,"\n"))
    }
    (HLfit6 <- HLfit(total_red ~ sex*env + (1|rep) + (0 + sex:env|line),
                     data = my.data, HLmethod="ML"))
    how(HLfit6) # 12s v3.0.34; quite longer (35.7s) in v3.0.44 which gives more exact result.
    testval6h <- -1536.08079543165 # chol -- good 
    zut <- try(testthat::expect_equal((res6h <- logLik(HLfit6)[[1L]]),testval6h,tol=1e-8),silent=TRUE) 
    if (inherits(zut,"try-error")) { 
      if (res6h<testval6h) strg <- " poorer" else strg <- " better" 
      cat(paste0("logLik(HLfit6)=",
                 res6h, strg, " than tested value ",testval6h," by ",res6h-testval6h ,"\n"))
    }
  }
}

