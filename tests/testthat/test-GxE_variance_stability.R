# From https://www.r-bloggers.com/2019/06/gen-experiments-fitting-a-stability-variance-model-with-r/
#library(spaMM)

cat(crayon::yellow("\ntest of yield stability analysis:"))

if (spaMM.getOption("example_maxtime")>1.5) {
  if (requireNamespace("agridat", quietly = TRUE)) {
    
    library("nlme") # more is needed than the functions imported in spaMM
    
    data("onofri.winterwheat", package="agridat")
    
    { # fit without any form of random coeff LHS
      model.mix <- lme(yield ~ gen, 
                       random=list(year = pdBlocked(list(pdIdent(~1),
                                                         pdIdent(~block - 1),
                                                         pdIdent(~gen - 1)))), # pdIdent -> effect on RHS
                       data=onofri.winterwheat)
      # same as
      spfit <- spaMM::fitme(yield ~ gen + (1|year)+(1|block %in% year)+(1|gen %in% year), data=onofri.winterwheat, method="REML")
      
      crit <- diff(range(c(logLik(model.mix),logLik(spfit))))
      testthat::test_that("whether lme() and spaMM give identical results", testthat::expect_true(crit<1e-8))
    }
    
    { # Fit  With different variances for each gen, using ranCoef syntax + activating isDiagFamily method by the fixed argument:
      model.mix <- lme(yield ~ gen, 
                       random=list(year = pdBlocked(list(pdIdent(~1),
                                                         pdIdent(~block - 1),
                                                         pdDiag(~gen - 1)))), # pdDiag -> gen both in LHS and RHS
                       data=onofri.winterwheat)
      # same as
      ranCoefs_for_diag <- function(nlevels) { # convenience function
        diagmat <- matrix(NA, ncol=nlevels,nrow=nlevels)
        diagmat[lower.tri(diagmat,diag=FALSE)] <- 0
        diagmat[lower.tri(diagmat,diag=TRUE)]
      }
      spfit <- spaMM::fitme(yield ~ gen + (1|year)+(1|block %in% year)+(0+gen|gen %in% year), data=onofri.winterwheat, method="REML",
                            fixed=list(ranCoefs=list("3"=ranCoefs_for_diag(8))))
      crit <- diff(range(c(logLik(model.mix),logLik(spfit))))
      testthat::test_that("whether lme() and spaMM give identical results (pDiag case)", testthat::expect_true(crit<1e-8))
      
    }
    
    detach("package:nlme")
    
  } else cat(crayon::bgGreen("Data 'onofri.winterwheat' not available for testing."))
}
