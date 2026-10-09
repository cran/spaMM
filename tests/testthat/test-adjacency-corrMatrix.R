cat(cli::col_yellow("\ntest-adjacency-corrMatrix: adjacency (dense,sparse) vs. corrMatrix() (dense,sparse * LevM or not), for HGLM with offset:\n"))

data("scotlip")

adjfitsp <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+(1|gridcode)+offset(log(expec)),
                adjMatrix=Nmatrix,
                rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                fixed=list(rho=0.1), 
                family=poisson(),data=scotlip, 
                control.HLfit=list(algebra="spprec"))
adjfit <- fitme(cases~I(prop.ag/10) +adjacency(1|gridcode)+(1|gridcode)+offset(log(expec)),
                  adjMatrix=Nmatrix,
                  rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                  fixed=list(rho=0.1), 
                  family=poisson(),data=scotlip, 
                  control.HLfit=list(algebra="decorr"))
if (spaMM.getOption("EigenDense_QRP_method")==".lmwithQR") {
  crit <- diff(range(logLik(adjfit),logLik(adjfitsp)))
  if (spaMM.getOption("fpot_tol")>0) {
    testthat::test_that(paste0("criterion was ",signif(crit,6)," from -168.1298"), testthat::expect_true(crit<2e-6) )
  } else testthat::expect_true(crit<2e-6)
  # precision has apparently changed in [dependencies installed with?] R devel-to-become-v4.6.0 
  testthat::expect_true(max(abs(range(get_predVar(adjfit)-get_predVar(adjfitsp))))<3e-6) 
} else {
  testthat::expect_true(diff(range(logLik(adjfit),logLik(adjfitsp)))<2e-8) 
  testthat::expect_true(max(abs(range(get_predVar(adjfit)-get_predVar(adjfitsp))))<8e-6)
}

testthat::expect_true(diff(range(predict(adjfit)[2:4,]-predict(adjfit,newdata=scotlip[2:4,])))<1e-12)


## same using corrMatrix()
if (spaMM.getOption("example_maxtime")>6.90) { 
  precmat <- diag(56)-0.1*Nmatrix   ## equivalent to adjacency model with rho=0.1
  precmat <- as(precmat,"sparseMatrix")
  colnames(precmat) <- rownames(precmat) <- seq(56)
  covmat <- solve(precmat)
  colnames(covmat) <- rownames(covmat) <- seq(56)
  (precfit <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                   covStruct=list(precision=precmat),
                   rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                   family=poisson(),data=scotlip))
  (is_spp <- .is_spprec_fit(precfit))
  if (is_spp) {
    precfitT <- precfit
    precfitF <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                       covStruct=list(corrMatrix=covmat),
                       rand.family=list(gaussian(),Gamma(log)), 
                       control.HLfit=list(sparse_precision=FALSE),
                       family=poisson(),data=scotlip)
    if (how(precfitF, verbose = FALSE)$fit_time < 0.95*how(precfit, verbose = FALSE)$fit_time) {
      message("'spprec' was chosen, but alternative was faster")
    } else message("'spprec' was chosen... OK")
  } else {
    precfitF <- precfit
    (precfitT <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                      covStruct=list(corrMatrix=covmat),
                      rand.family=list(gaussian(),Gamma(log)), 
                      control.HLfit=list(sparse_precision=TRUE),
                      family=poisson(),data=scotlip))
    if (how(precfitT, verbose = FALSE)$fit_time < 0.95*how(precfit, verbose = FALSE)$fit_time) {
      message("'spprec' was NOT chosen, but it was faster")
    } else message("'spprec' was NOT chosen... OK")
  }
  testthat::expect_equal(logLik(precfit),c(p_v=-168.12966973),tolerance=5e-5) ## all methods are equally sensitive to the initial value (note that one lambda->0)
  testthat::expect_true(max(abs(range(get_predVar(precfitT)-get_predVar(precfitF))))<9e-6)  
  precfitLM <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                     covStruct=list(precision=precmat),
                     rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                     family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=TRUE))
  covfitLM <- fitme(cases~I(prop.ag/10) +corrMatrix(1|gridcode)+(1|gridcode)+offset(log(expec)),
                    covStruct=list(corrMatrix=covmat),
                    rand.family=list(gaussian(),Gamma(log)), #verbose=c(TRACE=1L),
                    #fixed=list(lambda=c(0.1,0.05)), 
                    family=poisson(),data=scotlip,control.HLfit=list(LevenbergM=TRUE)) ## 
  if (spaMM.getOption("EigenDense_QRP_method")==".lmwithQR") {
    crit <- diff(range(logLik(precfitT),logLik(precfitF),logLik(precfitLM),logLik(covfitLM)))
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
