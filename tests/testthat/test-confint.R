# confint
cat("\ntest confint(). for HLfit() fit (see spaMMIntro for test with corrHLfit):\n")
data(wafers)
wfit <- HLfit(y ~X1+(1|batch),family=Gamma(log),data=wafers,HLmethod="ML")
ci <- confint(wfit,"X1")

testthat::expect_equal(ci$interval[[1]],0.01828659,tolerance=1e-4)
testthat::expect_equal(ci$interval[[2]],0.17271333,tolerance=1e-4)