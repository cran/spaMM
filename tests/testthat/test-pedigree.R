cat(crayon::yellow("\ntest pedigree:\n"))

{
  if(FALSE) {
    # if(requireNamespace("pedigreemm", quietly=TRUE)) {
    #   # derived from help("pedigreemm")
    #   p1 <- new("pedigree",
    #             sire = as.integer(c(NA,NA,1, 1,4,5)),
    #             dam  = as.integer(c(NA,NA,2,NA,3,2)),
    #             label = as.character(1:6))
    #   oldMDCopt <- options(Matrix.warnDeprecatedCoerce = 0) # pedigreemm not updated for Matrix 1.5.0
    #   A <- pedigreemm::getA(p1) ## relationship matrix
    #   options(oldMDCopt)
    # }
  }
  { ## by package 'kinship':
    # kinship2 package looks more 'active' but it seems to assume that one parent is known, both are known, 
    # which is not the case for 4th individual in the above example (=Table 2.1 or Mrode 2005). The syntax is
    # ped2 <- kinship2::pedigree(1:6, c(0,0,1, 1,4,5), c(0,0,2,2,3,2), c(1,2,2,1,1,3))
    # A_for_fixed4_th <- 2 * kinship2::kinship(ped2)
  }
  { # Manual version to avoid external dependencies:
    A <- matrix(NA, ncol=6,nrow=6)
    A[lower.tri(A,diag=TRUE)] <- c(8,0,4,4,4,2, 8,4,0,2,5, 8,2,5,4.5, 8,5,2.5, 9,5.5, 9)/8
    A <- forceSymmetric(A,uplo = "L")
    colnames(A) <- rownames(A) <- 1:6
  }
  ## data simulation
  cholA <- chol(A)  
  varU <- 0.4; varE <- 0.6; rep <- 20
  n <- rep*6
  set.seed(108)
  bStar <- rnorm(6, sd=sqrt(varU))
  b <- crossprod(as.matrix(cholA),bStar)
  ID <- rep(1:6, each=rep)
  e0 <- rnorm(n, sd=sqrt(varE))
  y <- b[ID]+e0
  obs <- data.frame(y=y,IDgen=ID,IDenv=ID) ## two copies of ID for readability of GLMM results
  ## fits
  verif1 <- fitme(y ~ 1+ corrMatrix(1|IDgen) , corrMatrix=A,data=obs,method="ML") ## tests the full augZXy method
  testthat::expect_true(diff(c(logLik(verif1),-132.037970787))<1e-6)
  obs$y01 <- ifelse(y<1.3,0,1)
  verif2 <- fitme(y01 ~ 1+ corrMatrix(1|IDgen)+(1|IDenv), corrMatrix=A,data=obs, 
        family=binomial(), method="ML")
  ## test-adjacency-corrMatrix also tests variants of corrMatrix + covStruct, but is longer
  prec_mat <- solve(A)
  colnames(prec_mat) <- rownames(prec_mat) <- rownames(A) # important
  verif3 <- fitme(y01 ~ 1+ corrMatrix(1|IDgen)+(1|IDenv) , covStruct=list(precision=prec_mat),
                  data=obs,family=binomial(), method="ML")
  # or
  verif4 <- HLCor(y01 ~ 1+ corrMatrix(1|IDgen)+(1|IDenv) , covStruct=list(precision=prec_mat),
                  data=obs,family=binomial(), HLmethod="ML")
  crit <- diff(range((c(logLik(verif2),logLik(verif3),logLik(verif4),-13.943504068))))
  testthat::test_that(paste0("criterion was ",signif(crit,6)," from -13.943504068"), testthat::expect_true(crit<2e-6))
}
