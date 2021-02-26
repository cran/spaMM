cat(crayon::yellow("\nTest of adjacency (long fits):"))

if (spaMM.getOption("example_maxtime")>71) { # actually quite faster ~20s but don't run it in every 'long test'
  ## example suggested by Jeroen van den Ochtend jeroenvdochtend@gmail.com Jeroen.vandenochtend@business.uzh.ch
  data("adjlg")
  oldop <- spaMM.options(sparse_precision=NULL, warn=FALSE) 
  system.time({
    fit.Frailty <- fitme(BUY ~ factor(month) + AGE + GENDER + X1*X2 + adjacency(1|ID),
                          data=adjlg,family = binomial(link = cloglog),method = "ML",
                          #control.HLfit=list(LevenbergM=FALSE), 
                          verbose=c(TRACE=interactive()), # to trace convergence 
                          adjMatrix=adjlgMat
    ) ## _F I X M E_ refitting lambda (by request) gives a lower lik... (-1552.946 v2.7.19 & v3.1.2) (point estimates are clearly different)
  }) # now it times as LevM=FALSE bc no longer start as LM by default
  # timings:
  # v3.6.7 21.36util ; v3.5.141 (New default use_ZA_L=TRUE) 29.24user (and LevM 70.25user)
  ##### Old default TRY_ZAX=NULL (ie use_ZA_L=NULL)
  # Presumably below 36 in v3.5.138. Need to reassess on an unperturbed CPU...
  ## 48.79util in v.3.0.35 #53.49util in v2.7.27 # 55.12util in v2.7.6 # 64.43utilin v2.7.1 ## 99.99 in v2.6.53 ## 102 in v.2.5.34 ## ~107 (v.2.4.102) 
  # LevenbergM=TRUE
  # 239.57util in v3.0.42; 244.64util in v3.0.35 
  ## v3.0.25-3.0.35 redefine the LevM controls:
  ## ~359 (v.2.4.102) => 250 in v.2.5.34; 242.44 in v2.6.53 # 156.44util in v2.7.1 # 152.23util in v2.7.6 # 146.38util in v.2.7.27
  expectedMethod <- "AUGI0_ZX_sparsePrecision" 
  if (interactive()) {
    if (! (expectedMethod %in% fit.Frailty$MME_method)) {
      message(paste('! ("',expectedMethod,'" %in% IRLS.Frailty$MME_method): was a non-default option selected?'))
    }
  } else testthat::expect_true(expectedMethod %in% IRLS.Frailty$MME_method) 
  # Older comment: spprec_LevM_D=="colSums" gives the highest lik, but "1" is fastest; "rowSums" may be slowest.
  spaMM.options(oldop)
}