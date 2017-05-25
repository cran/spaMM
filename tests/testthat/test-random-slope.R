cat("\ntest of random-slope model:")
if(require("lme4", quietly = TRUE)) {
  res <- HLfit(Reaction ~ Days + (Days|Subject), data = sleepstudy)
  expect_equal(res$APHLs$p_bv,-871.8141,tolerance=2e-4)
  expect_equal(res$lambda[[2]],34.91168,tolerance=2e-4) # checks correct reporting of lambdas
  expect_equal(predict(res)[1:6,],predict(res,newdata=sleepstudy[1:6,])[,1],re.form= ~ (Days|Subject))
} else {
  cat( "package 'lme4' not available, cannot run random-slope test.\n" )
}
