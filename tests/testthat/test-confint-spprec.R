if(spaMM.getOption("example_maxtime")>12L) {
  old_opt <- spaMM.options(sparse_precision=TRUE) # keep this one...
  source(paste0(spaMM::projpath(),"/package/tests/testthat/test-confint.R") )
  spaMM.options(old_opt)
}
