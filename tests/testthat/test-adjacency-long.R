cat(cli::col_yellow("\nTest of adjacency (long fits):"))

if (spaMM.getOption("example_maxtime")>61) { # actually faster ~19s
  ## example suggested by Jeroen van den Ochtend jeroenvdochtend@gmail.com Jeroen.vandenochtend@business.uzh.ch
  data("adjlg")
  fit.Frailty <- fitme(BUY ~ factor(month) + AGE + GENDER + X1*X2 + adjacency(1|ID),
                          data=adjlg,family = binomial(link = cloglog),method = c("ML","exp"),
                          #control.HLfit=list(LevenbergM=FALSE), 
                          verbose=c(TRACE=interactive()), # to trace convergence 
                          adjMatrix=adjlgMat) ## _F I X M E_ refitting lambda (by request) gives a lower lik... (-1552.946 v2.7.19 & v3.1.2) (point estimates are clearly different)
  how(fit.Frailty) 
  # "Model fitted by spaMM::fitme, version 4.7.7, in 18.99s using sparse-precision methods (with exp. info. matrix)."
  # 23s on 4.1.73; 26.4.s on 3.9.56. Longer since, possible effect of change in .dispFn...
  # timings previously using system.time()  [how() reports longer times ?! mysteries...]
  ################ LevenbergM=TRUE
  # "Model fitted by spaMM::fitme, version 4.7.7, in 68.62s using sparse-precision methods (with exp. info. matrix)."  # 239.57util in v3.0.42; 244.64util in v3.0.35 
  ## v3.0.25-3.0.35 redefine the LevM controls:
  ## ~359 (v.2.4.102) => 250 in v.2.5.34; 242.44 in v2.6.53 # 156.44util in v2.7.1 # 152.23util in v2.7.6 # 146.38util in v.2.7.27
  if (interactive()) {
    if (! .is_spprec_fit(fit.Frailty)) {
      message(paste('! .is_spprec_fit(fit.Frailty): was a non-default option selected?'))
    }
  } else testthat::expect_true(expectedMethod %in% how(IRLS.Frailty, verbose=FALSE)$MME_method) 
} else if (spaMM.getOption("example_maxtime")>20) cat(cli::bg_green(cli::col_black("\nIncrease maxtime above 61 to run the adjacency-long test !")))
