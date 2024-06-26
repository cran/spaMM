cat(crayon::yellow("\ntest-adjacency-corrMatrix: adjacency (dense,sparse) vs. corrMatrix() (dense,sparse * LevM or not), for HGLM with offset:\n"))

data("scotlip")

adjfitsp <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+(1|gridcode)+offset(log(expec)),
                adjMatrix=Nmatrix,
                rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                fixed=list(rho=0.1), 
                family=poisson(),data=scotlip)
#expectedMethod <- "AUGI0_ZX_spprec" ## has been dense for a long time
#actualMethods <- how(adjfit, verbose=FALSE)$MME_method
if (interactive()) {
  if ( ! .is_spprec_fit(adjfitsp)) {
    message(paste0('Actual method for adjfitsp differs from expected spprec: was a non-default option selected?'))
  }
} else testthat::expect_true(expectedMethod %in% adjfit$MME_method) 
adjfit <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+(1|gridcode)+offset(log(expec)),
                  adjMatrix=Nmatrix,
                  rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                  fixed=list(rho=0.1), 
                  family=poisson(),data=scotlip, 
                  control.HLfit=list(algebra="decorr"))
testthat::expect_true( ! .is_spprec_fit(adjfit))
if (spaMM.getOption("EigenDense_QRP_method")==".lmwithQR") {
  crit <- diff(range(logLik(adjfit),logLik(adjfitsp)))
  if (spaMM.getOption("fpot_tol")>0) {
    testthat::test_that(paste0("criterion was ",signif(crit,6)," from -168.1298"), testthat::expect_true(crit<2e-6) )
  } else testthat::expect_true(crit<2e-6)
  testthat::expect_true(max(abs(range(get_predVar(adjfit)-get_predVar(adjfitsp))))<4e-7)
} else {
  testthat::expect_true(diff(range(logLik(adjfit),logLik(adjfitsp)))<2e-8) 
  testthat::expect_true(max(abs(range(get_predVar(adjfit)-get_predVar(adjfitsp))))<8e-6)
}

testthat::expect_true(diff(range(predict(adjfit)[2:4,]-predict(adjfit,newdata=scotlip[2:4,])))<1e-12)


## same using corrMatrix()
if (spaMM.getOption("example_maxtime")>6.90) { 
  precmat <- diag(56)-0.1*Nmatrix   ## equivalent to adjacency model with rho=0.1
  colnames(precmat) <- rownames(precmat) <- seq(56)
  covmat <- solve(precmat)
  colnames(covmat) <- rownames(covmat) <- seq(56)
  precfit <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                   covStruct=list(precision=precmat),
                   rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                   #fixed=list(lambda=c(0.1,0.05)), 
                   family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=FALSE))
  testthat::expect_true(.is_spprec_fit(precfit))
  precfitLM <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                     covStruct=list(precision=precmat),
                     rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                     #fixed=list(lambda=c(0.1,0.05)), 
                     family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=TRUE))
  covfit <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                  covStruct=list(corrMatrix=covmat),
                  rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                  #fixed=list(lambda=c(0.1,0.05)), 
                  family=poisson(),data=scotlip)
  testthat::expect_true("matrix" %in% how(covfit, verbose=FALSE)$MME_method)
  testthat::expect_equal(logLik(covfit),c(p_v=-168.12966973),tolerance=5e-5) ## all methods are equally sensitive to the initial value (note that one lambda->0)
  testthat::expect_true(max(abs(range(get_predVar(covfit)-get_predVar(precfit))))<9e-6)  
  if (spaMM.getOption("EigenDense_QRP_method")==".lmwithQR") {
    testthat::expect_true(max(abs(range(get_predVar(adjfit)-get_predVar(covfit))))<9e-6) 
  } else testthat::expect_true(max(abs(range(get_predVar(adjfit)-get_predVar(covfit))))<6e-6) 
  covfitLM <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                    covStruct=list(corrMatrix=covmat),
                    rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                    #fixed=list(lambda=c(0.1,0.05)), 
                    family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=TRUE)) ## 
  if (spaMM.getOption("EigenDense_QRP_method")==".lmwithQR") {
    crit <- diff(range(logLik(covfit),logLik(covfitLM),logLik(precfitLM),logLik(adjfit),logLik(adjfitsp)))
    if (spaMM.getOption("fpot_tol")>0) {
      testthat::test_that(paste0("criterion was ",signif(crit,6)," from -168.12966973"), testthat::expect_true(crit<2e-6) )
    } else testthat::expect_true(crit<2e-6)
  } else testthat::expect_true(diff(range(logLik(covfit),logLik(covfitLM),logLik(precfitLM),logLik(adjfit),logLik(adjfitsp)))<3e-8) 
  # : correctness sensitive to w.resid <- damped_WLS_blob$w.resid.
  
  if (FALSE) { ## single IRLS fit sensitive to w.resid <- damped_WLS_blob$w.resid
    aaaa <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                  covStruct=list(corrMatrix=covmat),
                  rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=3L),
                  fixed=list(lambda=c(0.05,0.05)), 
                  family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=TRUE,sparse_precision=FALSE))
    bbbb <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                  covStruct=list(precision=precmat),
                  rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=3L),
                  fixed=list(lambda=c(0.05,0.05)), 
                  family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=TRUE,sparse_precision=FALSE))
  }
}

