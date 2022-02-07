cat(crayon::yellow("\ntest from corrFamily.Rd:\n"))

if (spaMM.getOption("example_maxtime")>3.3) {
  if (requireNamespace("agridat", quietly = TRUE)) {
    
    data("onofri.winterwheat", package="agridat")
    
    method <- "REML"
    # method <- "ML"
    
    ## Fitting a Toeplitz correlation model for temporal correlations
    Toepfn <- function(v) {
      toepmat <- Matrix::forceSymmetric(toeplitz(c(1,v)))
      toepmat <- regularize(toepmat, maxcondnum=1e12)
      rownames(toepmat) <- unique(onofri.winterwheat$year)
      toepmat
    }

    # 
    (Toepfit <- spaMM::fitme(
      yield ~ gen + corrFamily(1|year), data=onofri.winterwheat, method=method,
      covStruct=list(corrFamily=list(f=Toepfn, tpar=rep(1e-4,6))), 
      lower=list(corrPars=list("1"=rep(-0.999,6))), control.HLfit=list(algebra="decorr"),
      upper=list(corrPars=list("1"=rep(0.999,6)))))
    
    testthat::test_that("Check that logLik(Toepfit)==-152.4187",testthat::expect_true(diff(range(logLik(Toepfit),-152.4187))<1e-4))
    
    eigen(Corr(Toepfit)[[1]])$values # last one ~ 0
    
    { # permutation check ; the restricted range are designed to minimize numerical discrepancies due to singularity issues...
      (Toep1 <- spaMM::fitme(
        yield ~ 1 + corrFamily(1|year) + (1|gen), data=onofri.winterwheat, method=method,
        covStruct=list(corrFamily=list(f=Toepfn, tpar=rep(1e-4,6))), 
        lower=list(corrPars=list("1"=rep(-0.1,6))), control.HLfit=list(algebra="decorr"),
        upper=list(corrPars=list("1"=rep(0.1,6)))))
      
      (Toep2 <- spaMM::fitme(
        yield ~ 1 + (1|gen) + corrFamily(1|year), data=onofri.winterwheat, method=method,
        covStruct=list("1"=NULL, corrFamily=list(f=Toepfn, tpar=rep(1e-4,6))), 
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
        covStruct=list(corrFamily=list(f=drop0Toepfn, tpar=rep(0,6))), 
        lower=list(corrPars=list("1"=rep(-0.999,6))), control.HLfit=list(algebra="decorr"),
        upper=list(corrPars=list("1"=rep(0.999,6)))),
        warning = function(w){
          substr(w$message,1,7)
        })
      testthat::test_that("Check of 'tpar' check successful",testthat::expect_true(bla=="f(tpar)"))
    }
    
    { # more primitive implementation with map+fixed & internal regularization with warnings
      toep_map <- Matrix::drop0(toeplitz(seq(0,6)))
      rownames(toep_map) <- colnames(toep_map) <- unique(onofri.winterwheat$year)
      
      suppressWarnings(Toepfit <- spaMM::fitme(
        yield ~ gen + corrFamily(1|year), data=onofri.winterwheat, method=method,
        covStruct=list(corrFamily=list(map=toep_map, fixed="unit diag")), 
        lower=list(corrPars=list("1"=rep(-0.999,6))), 
        upper=list(corrPars=list("1"=rep(0.999,6)))))
      
      # The condition for automatic internal regularization being that chol() failed, it is not applied exactly
      # as in the given corrFamily$f, so the result is slightly different:
      testthat::test_that("Check that logLik(Toepfit)==-152.4236",testthat::expect_true(diff(range(logLik(Toepfit),-152.4236))<1e-4))
      
    }
    
    
    # one variance among years per each of 8 genotypes. 
    # This can be fitted as a constrained random-coefficient model
    
    # Diagonal matrix of NA's, represented as vector for its lower triangle:
    ranCoefs_for_diag <- function(nlevels) { 
      vec <- rep(0,nlevels*(nlevels+1L)/2L)
      vec[cumsum(c(1L,rev(seq(nlevels-1L)+1L)))] <- NA
      vec
    } 
    
    (by_rC <- spaMM::fitme(yield ~ 1 + (0+gen|year), data=onofri.winterwheat, method="REML",
                           fixed=list(ranCoefs=list("1"=ranCoefs_for_diag(8)))))
    
    # how(by_rC)
    
    gy_levels <- paste0(gl(8,1,length =56,labels=levels(onofri.winterwheat$gen)),":",
                        gl(7,8,labels=unique(onofri.winterwheat$year)))

    # diagf <- function(logvar) {
    #   corrm <- kronecker(Matrix::.symDiagonal(n=7),diag(x=exp(logvar)))
    #   rownames(corrm) <- colnames(corrm) <- gy_levels
    #   corrm
    # } 
    
    corr_map <- Matrix::forceSymmetric(kronecker(Matrix::.symDiagonal(n=7),diag(x=seq(8))))
    rownames(corr_map) <- gy_levels
    diagf <- function(logvar) {
      corrm <- corr_map
      corrm@x <- exp(logvar)[corrm@x]
      corrm
    } 
    
    (by_cF <- spaMM::fitme(
      yield ~ 1 + corrFamily(1|gen %in% year), data=onofri.winterwheat, method="REML",
      covStruct=list(corrFamily = list(f=diagf, tpar=rep(1,8))), 
      fixed=list(lambda=1), # verbose=c(TRACE=TRUE),
      # init=list(corrPars=list("1"=rep(log(O.1),8))), # 'init' optional 
      lower=list(corrPars=list("1"=rep(log(1e-6),8))), # 'lower' and 'upper' required
      upper=list(corrPars=list("1"=rep(log(1e6),8)))))       # 0.68s
    
    # how(by_cF)
    
    testthat::test_that("Check that corrFamily and ranCoef fit are consistent",testthat::expect_true(diff(range(logLik(by_rC),logLik(by_cF)))<1e-08)) 

  } else (message("package 'agridat' not available for test-corrFamily.")) 
}
