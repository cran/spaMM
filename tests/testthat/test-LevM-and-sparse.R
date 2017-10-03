if (spaMM.getOption("example_maxtime")>12) { 
  cat("\ntest Levenberg and/or sparse algo, for HGLM with offset:\n")
  # from devel/sparsePrecision/scratch
  
  data("scotlip")
  
  precmat <- diag(56)-0.1*Nmatrix
  colnames(precmat) <- rownames(precmat) <- seq(56)
  covmat <- solve(precmat)
  colnames(covmat) <- rownames(covmat) <- seq(56)
  
  oldop <- spaMM.options(sparse_precision=FALSE)
  covfit <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                  covStruct=list(corrMatrix=covmat),
                  rand.family=list(gaussian(),Gamma(log)), verbose=c(TRACE=1L),
                  #fixed=list(lambda=c(0.1,0.05)), 
                  family=poisson(),data=scotlip)
  testthat::expect_true("dgCMatrix" %in% covfit$MME_method)
  testthat::expect_equal(logLik(covfit),c(p_v=-168.129669688))
  covfitLM <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                    covStruct=list(corrMatrix=covmat),
                    rand.family=list(gaussian(),Gamma(log)), verbose=c(TRACE=1L),
                    #fixed=list(lambda=c(0.1,0.05)), 
                    family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=c(LevenbergM=TRUE))) ## 
  testthat::expect_equal(logLik(covfit),logLik(covfitLM),tolerance=2e-4) ## diff 1.212993e-08 in v2.1.81
  
  
  spaMM.options(sparse_precision=TRUE)
  precfit <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                   covStruct=list(precision=precmat),
                   rand.family=list(gaussian(),Gamma(log)), verbose=c(TRACE=1L),
                   #fixed=list(lambda=c(0.1,0.05)), 
                   family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=FALSE))
  testthat::expect_true("AUGI0_ZX_sparsePrecision" %in% precfit$MME_method)
  testthat::expect_equal(logLik(covfit),logLik(precfit))
  precfitLM <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                     covStruct=list(precision=precmat),
                     rand.family=list(gaussian(),Gamma(log)), verbose=c(TRACE=1L),
                     #fixed=list(lambda=c(0.1,0.05)), 
                     family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=TRUE))
  testthat::expect_equal(logLik(covfit),logLik(precfitLM))
  spaMM.options(oldop)
  
}

