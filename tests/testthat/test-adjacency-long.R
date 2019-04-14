cat("\nTest of adjacency (long fits):")

if (spaMM.getOption("example_maxtime")>360) { # 
  ## fixme This is sparse_precision LevM: improve efficiency? 
  ## example suggested by Jeroen van den Ochtend jeroenvdochtend@gmail.com Jeroen.vandenochtend@business.uzh.ch
  ## it's no use to try sparse_precision=FALSE bc bc the augmented sigma matrix is huge
  ###################### spaMM.options(sparse_precision=FALSE) 
  ## make sure that the wrong value is not set:
  data("adjlg")
  oldop <- spaMM.options(sparse_precision=NULL) 
  system.time({
    IRLS.Frailty <- fitme(BUY ~ factor(month) + AGE + GENDER + X1*X2 + adjacency(1|ID),
                          data=adjlg,family = binomial(link = cloglog),method = "ML",
                          control.HLfit=list(LevenbergM=FALSE), ## inhibits default for binary data 
                          verbose=c(TRACE=interactive()), # to trace convergence 
                          adjMatrix=adjlgMat
    ) ## F I X M E refitting lambda (by request) gives a lower lik...
  }) ##  ~107 (v.2.4.102)  ## 102 in v.2.5.34 ## 99.99 in v2.6.53
  expectedMethod <- "AUGI0_ZX_sparsePrecision" 
  if (interactive()) {
    if (! (expectedMethod %in% IRLS.Frailty$MME_method)) {
      message(paste('! ("',expectedMethod,'" %in% IRLS.Frailty$MME_method): was a non-default option selected?'))
    }
  } else testthat::expect_true(expectedMethod %in% IRLS.Frailty$MME_method) 
  system.time({
    LevM.Frailty <- fitme(BUY ~ factor(month) + AGE + GENDER + X1*X2 + adjacency(1|ID),
                          data=adjlg,family = binomial(link = cloglog),method = "ML",
                          verbose=c(TRACE=interactive()), # to trace convergence 
                          #fixed=list(rho = -0.0294184,  lambda = 0.241825),
                          adjMatrix=adjlgMat
    )
  }) ## ~359 (v.2.4.102) => 250 in v.2.5.34; 242.44 in v2.6.53
  # spprec_LevM_D=="colSums" gives the highest lik, but "1" is fastest; "rowSums" may be slowest.
  spaMM.options(oldop)
}