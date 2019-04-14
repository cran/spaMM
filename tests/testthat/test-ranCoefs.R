if (file.exists((privdata <- "C:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_private/all_fitness.txt"))) {
  cat("\ntest ranCoefs:")
  my.data <- read.table(privdata, header = TRUE, sep = "\t",dec = ".")
  my.data$line <- factor(as.character(my.data$line))
  my.data <- na.omit(my.data)
  if (TRUE || spaMM.getOption("example_maxtime")>2.5) {
    # -1559.813 1.8s v 2.4.146
    (fitme3 <- fitme(total_red ~ sex*env + (1|rep) + (0 + env|line),
                     data = my.data, method="ML"))
    how(fitme3)
    # final precision depends on two steps ! And ultimately on HLfit/hatval precision in the refit
    testval3o <- 1559.81295627
    zut <- try(testthat::expect_equal((res3o <- attr(attr(fitme3,"optimInfo")$optim.pars,"optr")$objective),
                                      testval3o,tol=1e-8),silent=TRUE)
    if (inherits(zut,"try-error")) {
      if (res3o>testval3o) strg <- " poorer" else strg <- " better" 
      cat(paste0("objective from fitme3's optimInfo=",
                 res3o, strg, " than tested value ",testval3o," by ",res3o-testval3o ,"\n"))
    }
    testval3f <- -1559.81272901
    zut <- try(testthat::expect_equal((res3f <- logLik(fitme3)[[1L]]),testval3f,tol=1e-8),silent=TRUE) ## -1536.080818 is quite possible
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
    # -1559.815 # 3.6 s v 2.4.146
    (HLfit3 <- HLfit(total_red ~ sex*env + (1|rep) + (0 + env|line),
                     data = my.data, HLmethod="ML"))
    how(HLfit3)
    testval3h <- -1559.81431734
    zut <- try(testthat::expect_equal((res3h <- logLik(HLfit3)[[1L]]),testval3h,tol=1e-8),silent=TRUE) ## -1536.080818 is quite possible
    if (inherits(zut,"try-error")) {
      if (res3h<testval3h) strg <- " poorer" else strg <- " better" 
      cat(paste0("logLik(HLfit3)=",
                 res3h, strg, " than tested value ",testval3h," by ",res3h-testval3h ,"\n"))
    }
    #  -1536.081 24.9s v 2.4.146
    (fitme6 <- fitme(total_red ~ sex*env + (1|rep) + (0 + sex:env|line),
                     data = my.data, method="ML"))  
    how(fitme6)
    testval6o <- 1536.08100566145
    zut <- try(testthat::expect_equal((res6o <- attr(attr(fitme6,"optimInfo")$optim.pars,"optr")$objective),
                                      testval6o,tol=1e-8),silent=TRUE)
    if (inherits(zut,"try-error")) {
      if (res6o>testval6o) strg <- " poorer" else strg <- " better" 
      cat(paste0("objective from fitme6's optimInfo=",
                 res6o, strg, " than tested value ",testval6o," by ",res3o-testval3o ,"\n"))
    }
    testval6f <- -1536.08082814
    zut <- try(testthat::expect_equal((res6f <- logLik(fitme6)[[1L]]),testval6f,tol=1e-8),silent=TRUE) ## -1536.080818 is quite possible
    if (inherits(zut,"try-error")) {
      if (res6f<testval6f) strg <- " poorer" else strg <- " better" 
      cat(paste0("logLik(fitme6)=",
                 res6f, strg, " than tested value ",testval6f," by ",res6f-testval6f ,"\n"))
    }
    #  -1536.083 37.9 v 2.4.146 
    (HLfit6 <- HLfit(total_red ~ sex*env + (1|rep) + (0 + sex:env|line),
                     data = my.data, HLmethod="ML"))
    how(HLfit6)
    testval6h <- -1536.083037
    zut <- try(testthat::expect_equal((res6h <- logLik(HLfit6)[[1L]]),testval6h,tol=1e-8),silent=TRUE) ## -1536.080818 is quite possible
    if (inherits(zut,"try-error")) {
      if (res6h<testval6h) strg <- " poorer" else strg <- " better" 
      cat(paste0("logLik(HLfit6)=",
                 res6h, strg, " than tested value ",testval6h," by ",res6h-testval6h ,"\n"))
    }
  }
}

