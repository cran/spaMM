cat("\ntest of bboptim compatibility:")

if (requireNamespace("blackbox",quietly=TRUE)) {
  blackbox::blackbox.options(spaMM.options("example_maxtime"))
  #blackbox.options(example_maxtime=5) ## allows basic example
  example("bboptim",package="blackbox") ## detected calcpredVar issue for etaFix   
}

# also feasible but slow:
if (requireNamespace("Infusion",quietly=TRUE)) {
  Infusion::Infusion.options(spaMM.options("example_maxtime"))
  #Infusion.options(example_maxtime=20)  ## allows basic example
  #Infusion.options(example_maxtime=120)  ## allows refine()
  example("Infusion",package="Infusion",ask=FALSE) ## detected calc_logdisp_cov issues   
}

if (requireNamespace("IsoriX",quietly=TRUE)) { ## Checks that exports/imports are OK
  IsoriX::IsoriX.options(spaMM.options("example_maxtime"))
  IsoriX::IsoriX.options(dont_ask = TRUE) ## for plots
  chk <- try(example("isofit",package="IsoriX",ask=FALSE),silent=TRUE) ## to detect e.g. do_TRACE issues... or corrPars$rho...
  if (inherits(chk,"try-error")) {
    warning(paste('example("isofit")',"produced",chk))
  } else if ( exists("GermanFit")) { ## !exists if IsoriX.getOption("example_maxtime") too low
    testthat::expect_equal(get_ranPars(GermanFit$mean.fit,which="corrPars")[["2"]]$rho,0.0001176491)
    #testthat::expect_equal(get_ranPars(GermanFit$disp.fit,which="corrPars")$rho,0.0001176491) ## lower bound for rho in both cases
  }
  chk <- try(example("isoscape",package="IsoriX",ask=FALSE),silent=TRUE) 
  if (inherits(chk,"try-error")) {
    warning(paste('example("isoscape")',"produced",chk))
  } else if ( exists("plot.mean.respVar")) { ## !exists if IsoriX.getOption("example_maxtime") too low
    testthat::expect_equal(isoscape$disp.predVar@data@values[1],0.080403352)
    testthat::expect_equal(isoscape$mean.predVar@data@values[1],148.3338)
  }
  chk <- try(example("calibfit",package="IsoriX",ask=FALSE),silent=TRUE) 
  if (inherits(chk,"try-error")) {
    warning(paste('example("calibfit")',"produced",chk))
  } else if ( exists("calib")) { ## !exists if IsoriX.getOption("example_maxtime") too low
    testthat::expect_equal(logLik(calib$calib.fit),c(p_bv=-1135.841334))
  }
  chk <- try(example("isofind",package="IsoriX",ask=FALSE),silent=TRUE) ## ~101s but the condition is 200
  if (inherits(chk,"try-error")) {
    warning(paste('example("isofind")',"produced",chk))
  } else if ( exists("assignment.GP")) { ## !exists if IsoriX.getOption("example_maxtime") too low
    testthat::expect_equal(assignment.GP$group$pv@data@values[12240],0.798571335)
  }
}