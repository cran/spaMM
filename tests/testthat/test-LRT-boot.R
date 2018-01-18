if (spaMM.getOption("example_maxtime")>(13.7+24)) { ## user time + system.time for parallele setup
  cat("test LRT with bootstrap, parallel or not:\n")
  data("salamander")
  fullfit <-HLfit(cbind(Mate,1-Mate)~TypeF+(1|Female)+(1|Male),family=binomial(),data=salamander,
                  HLmethod="ML",control.HLfit = list(LevenbergM=FALSE))
  nullfit <-HLfit(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male),family=binomial(),data=salamander,
                  HLmethod="ML",control.HLfit = list(LevenbergM=FALSE))
  set.seed(123)
  pv1 <- LRT(nullfit,fullfit,boot.repl=10)$BartBootLRT$p_value
  if (require(doSNOW,quietly=TRUE)) {
    set.seed(123)
    pv2 <- LRT(nullfit,fullfit,boot.repl=10,nb_cores=2)$BartBootLRT$p_value
    unloadNamespace("doSNOW")
    testthat::expect_equal(pv1,pv2,tolerance=1e-6)
  }
  set.seed(123)
  pv3 <- LRT(nullfit,fullfit,boot.repl=10,nb_cores=2)$BartBootLRT$p_value
  testthat::expect_equal(pv1,pv3,tolerance=1e-6)
} else cat("test LRT with bootstrap, parallel or not: increase example_maxtime (38s) to run this test.\n")

