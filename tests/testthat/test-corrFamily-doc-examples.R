cat(crayon::yellow("\ntest from test-corrFamily-doc-examples.Rd:\n"))

if (requireNamespace("agridat", quietly = TRUE)) {
  # spaMM.options(example_maxtime=60) # keep it handy...
  
  {cat(crayon::yellow("\ntest from corrFamily-design.Rd:\n"))
    
    method <- "REML"
    
    if (exists("by_cF")) rm("by_cF")
    example("corrFamily-design", echo=FALSE)
    if (exists("by_cF")) { # if the examples have been run, check their results:
      testthat::test_that("Check that corrFamily and ranCoef fit are consistent",testthat::expect_true(diff(range(logLik(by_rC),logLik(by_cF)))<1e-08))  
      testthat::test_that("Check that logLik(Toepfit)==-152.4187",testthat::expect_true(diff(range(logLik(Toepfit),-152.4187))<1e-4))
      
      { # permutation check ; the restricted range are designed to minimize numerical discrepancies due to singularity issues...
        (Toep1 <- spaMM::fitme(
          yield ~ 1 + corrFamily(1|year) + (1|gen), data=onofri.winterwheat, method=method,
          covStruct=list(corrFamily=list(Cf=Toepfn, tpar=rep(1e-4,6))), 
          lower=list(corrPars=list("1"=rep(-0.1,6))), control.HLfit=list(algebra="decorr"),
          upper=list(corrPars=list("1"=rep(0.1,6)))))
        (Toep2 <- spaMM::fitme(
          yield ~ 1 + (1|gen) + corrFamily(1|year), data=onofri.winterwheat, method=method,
          covStruct=list("1"=NULL, corrFamily=list(Cf=Toepfn, tpar=rep(1e-4,6))), 
          lower=list(corrPars=list("2"=rep(-0.1,6))), control.HLfit=list(algebra="decorr"),
          upper=list(corrPars=list("2"=rep(0.1,6)))))
        testthat::expect_true(diff(range(logLik(Toep1),logLik(Toep2)))<1e-10)
      }
      
      { # Checks the tpar check...
        drop0Toepfn <- function(v) {
          toepmat <- Matrix::drop0(toeplitz(c(1,v)))
          toepmat <- regularize(toepmat, maxcondnum=1e12)
          rownames(toepmat) <- colnames(toepmat) <- unique(onofri.winterwheat$year)
          toepmat
        } 
        bla <- tryCatch(Toepfit <- spaMM::fitme(
          yield ~ gen + corrFamily(1|year), data=onofri.winterwheat, method=method,
          covStruct=list(corrFamily=list(Cf=drop0Toepfn, tpar=rep(0,6))), 
          lower=list(corrPars=list("1"=rep(-0.999,6))), control.HLfit=list(algebra="decorr"),
          upper=list(corrPars=list("1"=rep(0.999,6)))),
          warning = function(w){
            substr(w$message,1,8)
          })
        testthat::test_that("Check of 'tpar' check successful", # Catching "Cf(tpar) was not a least sparse matrix in the corrFamily. Check 'tpar'.",
                            testthat::expect_true(bla=="Cf(tpar)"))
      }
    }
  }
  
  {cat(crayon::yellow("\ntest from ARp.Rd:\n"))
    
    if (exists("AR3_fix")) rm("AR3_fix")
    example("ARp", echo=FALSE)
    if (exists("AR3_fix")) { # if the examples have been run, check their results:
      testthat::test_that("Check that 'fixed' and covStruct=list(corrFamily=<constructor>(., fixed=.)) give equivalent results",
                          testthat::expect_true(diff(range(logLik(AR3fix),logLik(AR3_fix)))<1e-14))
    }
  }
  
} else (message("package 'agridat' not available for test-corrFamily.")) 






