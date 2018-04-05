cat("\ntest of bboptim compatibility:")

if (require("blackbox",quietly=TRUE)) {
  example("bboptim") ## detected calcpredVar issue for etaFix   
}

# also feasible but slow:
if (require("Infusion",quietly=TRUE)) {
  #Infusion.options(example_maxtime=20)  ## allows basic example
  #Infusion.options(example_maxtime=120)  ## allows refine()
  example("Infusion") ## detected calc_logdisp_cov issues   
}

if (require("IsoriX",quietly=TRUE)) { ## Checks that exports/imports are OK
  IsoriX.options(spaMM.options("example_maxtime"))
  IsoriX.options(dont_ask = TRUE) ## for plots
  chk <- try(example("isofit"),silent=TRUE) ## to detect e.g. do_TRACE issues... or corrPars$rho...
  if (inherits(chk,"try-error")) {
    warning(paste('example("isofit")',"produced",chk))
  } else if ( ! is.null(chk$value)) { ## NULL value if example_maxtime too low
    testthat::expect_equal(chk$value$mean.fit$corrPars[["2"]]$rho,0.0001176491)
    testthat::expect_equal(chk$value$disp.fit$corrPars[["2"]]$rho,0.0001176491) ## lower bound for rho in both cases
  }
}