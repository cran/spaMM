cat(crayon::yellow("\nTest of adjacency (long fits):"))

if (spaMM.getOption("example_maxtime")>61) { # actually faster ~23s
  ## example suggested by Jeroen van den Ochtend jeroenvdochtend@gmail.com Jeroen.vandenochtend@business.uzh.ch
  data("adjlg")
  fit.Frailty <- fitme(BUY ~ factor(month) + AGE + GENDER + X1*X2 + adjacency(1|ID),
                          data=adjlg,family = binomial(link = cloglog),method = c("ML","exp"),
                          #control.HLfit=list(LevenbergM=FALSE), 
                          verbose=c(TRACE=interactive()), # to trace convergence 
                          adjMatrix=adjlgMat) ## _F I X M E_ refitting lambda (by request) gives a lower lik... (-1552.946 v2.7.19 & v3.1.2) (point estimates are clearly different)
  how(fit.Frailty) # 23s on 4.1.73; 26.4.s on 3.9.56. Longer since, possible effect of change in .dispFn...
  # timings previously using system.time()  [how() reports longer times ?! mysteries...]
  # v3.9.6   20.27        4.37       24.59
  # v3.6.7   21.36util ; 
  # v3.5.141 (New default use_ZA_L=TRUE) 29.24user (and LevM 70.25user)
  ##### Old default TRY_ZAX=NULL (ie use_ZA_L=NULL)
  # Presumably below 36 in v3.5.138. 
  ## 48.79util in v.3.0.35 #53.49util in v2.7.27 # 55.12util in v2.7.6 # 64.43utilin v2.7.1 ## 99.99 in v2.6.53 ## 102 in v.2.5.34 ## ~107 (v.2.4.102) 
  # LevenbergM=TRUE
  # 239.57util in v3.0.42; 244.64util in v3.0.35 
  ## v3.0.25-3.0.35 redefine the LevM controls:
  ## ~359 (v.2.4.102) => 250 in v.2.5.34; 242.44 in v2.6.53 # 156.44util in v2.7.1 # 152.23util in v2.7.6 # 146.38util in v.2.7.27
  if (interactive()) {
    if (! .is_spprec_fit(fit.Frailty)) {
      message(paste('! .is_spprec_fit(fit.Frailty): was a non-default option selected?'))
    }
  } else testthat::expect_true(expectedMethod %in% how(IRLS.Frailty, verbose=FALSE)$MME_method) 
  # Older comment: spprec_LevM_D=="colSums" gives the highest lik, but "1" is fastest; "rowSums" may be slowest.
} else if (spaMM.getOption("example_maxtime")>20) cat(crayon::bgGreen("\nIncrease maxtime above 61 to run the adjacency-long test !"))
