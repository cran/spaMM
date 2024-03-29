cat(crayon::yellow("\ntest fixedLRT:\n"))
# fixedLRT

data("blackcap")
## result comparable to the corrHLfit examples based on blackcap
set.seed(123)
nb_cores <- 1 ## don't use parallel here (slow for boot.repl=3 ) and may be other 'check' issues
# nb_cores <- parallel::detectCores(logical=FALSE)-1
(fl <- fixedLRT(null.formula=migStatus ~ 1 + Matern(1|longitude+latitude),
         formula=migStatus ~ means + Matern(1|longitude+latitude), 
         HLmethod='ML',data=blackcap,init=list(phi=1e-6),boot.repl=3,nb_cores=nb_cores)) #, control=list(optimizer="bobyqa")) # may help if .dispFn() is modified

## phi=1e-6 must be automatically converted to 1e-4 : potential test of changes in .calc_inits_dispPars()
## but the init forces inner estimation of phi, which affects the following tests.
testthat::expect_equal(fl$basicLRT$p_value,0.00806473,tolerance=2e-5)
testthat::expect_equal(fl$BartBootLRT$p_value,structure(0.008546614,boot_type="marginal"),tolerance=1e-5)
testthat::expect_true(length(fl$bootInfo$RNGstates)==626) ## check that it's there

# With prior weights

if (spaMM.getOption("example_maxtime")>3) {
  data("hills", package="MASS")
  # serial
  (fl <- fixedLRT(null.formula=time ~ 1,
                  formula=time ~ climb, data = hills, weights.form= ~ 1/dist^2, method="ML",
                  boot.repl=10,nb_cores=1L, seed=123) )
  # parallel
  (fl <- fixedLRT(null.formula=time ~ 1,
                  formula=time ~ climb, data = hills, weights.form= ~ 1/dist^2, method="ML",
                  boot.repl=10,nb_cores=2L, seed=123) )
}

