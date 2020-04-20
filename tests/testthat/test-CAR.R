cat(crayon::yellow("\ntest CAR and SEM:\n"))

data("scotlip")

# see also 'old donttest' examples

## same without optim: run in scotlip examples; cf also autoregressive.Rd for ML fits

#set.seed(124)
{
  set.seed(129) ## many samples will diverge (ML on binary response) or undemonstratively be fitted by extreme rho values
  eigenv <- eigen(Nmatrix, symmetric=TRUE) 
  Lmat <- eigenv$vectors %*% diag(sqrt(1/(1-0.17*eigenv$values)))
  lp <- 0.1 + 3* Lmat %*% rnorm(ncol(Lmat)) ## single intercept beta =0.1; lambda=3
  resp <- rbinom(ncol(Lmat),1,1/(1+exp(-lp)))
  donn <- data.frame(npos=resp,nneg=1-resp,gridcode=scotlip$gridcode)
}

# CAR by Laplace with 'inner' estimation of rho
blob1 <- HLCor(cbind(npos,nneg)~1 +adjacency(1|gridcode),
          adjMatrix=Nmatrix,family=binomial(probit),data=donn,HLmethod="ML",control.HLfit = list(LevenbergM=FALSE)) ## 2 s.
testthat::expect_equal(get_ranPars(blob1,which = "corrPars")[["1"]]$rho, 0.07961275)
#AIC(blob1)

if (FALSE) { ## HLCor/corrHlfit already compared on scotlip by test-spaMM.R
  # corrHLfit without corners was poor here
  # CAR by Laplace with 'outer' estimation of rho
  blob2 <- fitme(cbind(npos,nneg)~1 +adjacency(1|gridcode),
                 adjMatrix=Nmatrix,family=binomial(probit),data=donn,method="ML",control.HLfit = list(LevenbergM=FALSE)) 
  #AIC(blob2) 
  testthat::expect_true(diff(range(AIC(blob2,verbose=FALSE)-AIC(blob1,verbose=FALSE)))<0.1) # effective-df calculation sensitive to small difs in fit
}

if ( (! "covr" %in% loadedNamespaces()) && 
     file.exists((privtest <- "C:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_other_pack/test-probitgem.R"))) {
  source(privtest)  
} # including another AIC() check
