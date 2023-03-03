cat(crayon::yellow("test LRT with bootstrap, parallel or not:\n"))
if (spaMM.getOption("example_maxtime")>(13.7+24)) { ## user time + system.time for parallel setup
  cat("test LRT()")
  data("salamander")
  fullfit <- HLfit(cbind(Mate,1-Mate)~TypeF+(1|Female)+(1|Male),family=binomial(),data=salamander,
                  HLmethod="ML")
  nullfit <- HLfit(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male),family=binomial(),data=salamander,
                  HLmethod="ML")
  set.seed(123)
  pv1 <- LRT(nullfit,fullfit,boot.repl=10)$BartBootLRT$p_value
  set.seed(123)
  pv3 <- LRT(nullfit,fullfit,boot.repl=10,nb_cores=5)$BartBootLRT$p_value
  crit <- diff(range(pv1-pv3))
  testthat::test_that(paste0("Max difference in prediction was",signif(crit,6)," >1e-6"), testthat::expect_true(crit<1e-6)) 
} else cat("increase example_maxtime (38s) to run LRT() test.\n")

if (spaMM.getOption("example_maxtime")>18) { ## user time + system.time for parallel setup
  cat("test spaMM_boot() with different backends:\n")
  data("blackcap")
  
  # Generate fits of null and full models:
  lrt <- fixedLRT(null.formula=migStatus ~ 1 + Matern(1|longitude+latitude),
                  formula=migStatus ~ means + Matern(1|longitude+latitude), 
                  HLmethod='ML',data=blackcap)
  
  myfun <- function(y, what=NULL, lrt, ...) { 
    data <- lrt$fullfit$data
    data$migStatus <- y ## replaces original response (! more complicated for binomial fits)
    full_call <- getCall(lrt$fullfit) ## call for full fit
    full_call$data <- data
    ReSu <- eval(full_call) ## fits the full model on the simulated response
    if (!is.null(what)) ReSu <- eval(what)(ReSu=ReSu)
    return(ReSu) ## the fit, or anything produced by evaluating 'what'
  }

  ## nb_cores=1, serial pbapply 
  set.seed(123)
  spaMM_boot(lrt$nullfit, simuland = myfun, nsim=4, type="marginal", 
             what=quote(function(ReSu) fixef(ReSu)[2]), lrt=lrt,nb_cores=1L)[["bootreps"]]    
  
  ## foreach+doSNOW (this was slow until the firewall settings were modified...)
  if ( (! "covr" %in% loadedNamespaces()) && 
       file.exists((privtest <- "C:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_other_pack/test-doSNOW.R"))) {
    source(privtest)
  } # i.e.: 
  # library("doSNOW")
  # set.seed(123)
  # spaMM_boot(lrt$nullfit, simuland = myfun, nsim=4,
  #            what=quote(fixef(ReSu)[2L]), lrt=lrt,nb_cores=4L)[["bootreps"]]
  # unloadNamespace("doSNOW")
  
  ## no doSNOW => parallel pbapply 
  set.seed(123)
  # if (exists("ReSu")) rm(ReSu) ## otherwise any error in the following code may result in a confusing message
  spaMM_boot(lrt$nullfit, simuland = myfun, nsim=4, control.foreach = list(.errorhandling="pass"),
             type="marginal", what=quote(function(ReSu) fixef(ReSu)[2]), lrt=lrt,nb_cores=4L)[["bootreps"]]    
  
} 

