print("test of random-slope model:")
if(require("lme4", quietly = TRUE)) {
  res <- HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy)
  expect_equal(res$APHLs$p_bv,-871.8141,tolerance=2e-4)
} else {
  cat( "package 'lme4' not available, cannot run random-slope test.\n" )
}
