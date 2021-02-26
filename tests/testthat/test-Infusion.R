cat(crayon::yellow("\ntest of Infusion compatibility:"))

if (requireNamespace("Infusion",quietly=TRUE)) {
  Infusion::Infusion.options(spaMM.options("example_maxtime"))
  #Infusion.options(example_maxtime=20)  ## allows basic example
  #Infusion.options(example_maxtime=120)  ## allows refine()
  if (FALSE) { # for manual testing
    library(Infusion)
    Infusion.options(example_maxtime=60)
  }
  if (do.call("require",list(package="Rmixmod", quietly = TRUE))) {
    example("example_raw",package="Infusion",ask=FALSE) ## detected calc_logdisp_cov issues 
  } else {
    cat( "package 'Rmixmod' not available, cannot run Infusion test.\n" )
  }
}

