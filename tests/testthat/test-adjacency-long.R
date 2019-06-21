cat("\nTest of adjacency (long fits):")

if (spaMM.getOption("example_maxtime")>227) { # 
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
    ) ## F I X M E refitting lambda (by request) gives a lower lik... (-1552.946 v2.7.19)
  }) ##  ~107 (v.2.4.102)  ## 102 in v.2.5.34 ## 99.99 in v2.6.53 # 64.43utilin v2.7.1 # 55.12util in v2.7.6 #53.49util in v2.7.27  
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
  }) ## ~359 (v.2.4.102) => 250 in v.2.5.34; 242.44 in v2.6.53 # 156.44util in v2.7.1 # 152.23util in v2.7.6 # 146.38util in v.2.7.27
  # spprec_LevM_D=="colSums" gives the highest lik, but "1" is fastest; "rowSums" may be slowest.
  spaMM.options(oldop)
}