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
  CARSEMd <- data.frame(npos=resp,nneg=1-resp,gridcode=scotlip$gridcode)
}

# CAR by Laplace with 'inner' estimation of rho
blob1 <- HLCor(cbind(npos,nneg)~1 +adjacency(1|gridcode),
          adjMatrix=Nmatrix,family=binomial(probit),data=CARSEMd,HLmethod="ML") ## ~1.27 s.
# it is crashes the session after installing a new R version, recompile probitgem with a clean /src directory before trying anything else. 
crit <- diff(range(get_ranPars(blob1,which = "corrPars")[["1"]]$rho, 0.07961275))
#if (spaMM.getOption("fpot_tol")>0) {
  try(testthat::test_that(paste0("criterion was ",signif(crit,6)," from 0.07961275"), testthat::expect_true(crit<5e-9))) # bobyqa finds 0.04582924 ('flat' p_bv)
#} else testthat::expect_true(crit<5e-9)

#AIC(blob1)
VarCorr(blob1)

if (FALSE) { ## HLCor/corrHlfit already compared on scotlip by test-spaMM.R
  # corrHLfit without corners was poor here
  # CAR by Laplace with 'outer' estimation of rho
  blob2 <- fitme(cbind(npos,nneg)~1 +adjacency(1|gridcode),
                 adjMatrix=Nmatrix,family=binomial(probit),data=CARSEMd,method="ML",control.HLfit = list(LevenbergM=FALSE)) 
  #AIC(blob2) 
  testthat::expect_true(diff(range(AIC(blob2,verbose=FALSE)-AIC(blob1,verbose=FALSE)))<0.1) # effective-df calculation sensitive to small difs in fit
}

if ( (! "covr" %in% loadedNamespaces()) && 
     file.exists((privtest <- "C:/home/francois/travail/stats/spaMMplus/spaMM/package/tests_other_pack/test-probitgem.R"))) {
  source(privtest)  
} # including another AIC() check

## test handling of missing data
if (spaMM.getOption("example_maxtime")>0.7) {
  scotli <- scotlip
  scotli$cases[1] <- NA
  rownames(Nmatrix) <- colnames(Nmatrix) <- scotli$gridcode # needed here
  mdfit <- fitme(cases~I(prop.ag/10)+adjacency(1|gridcode)+offset(log(expec)), data=scotli, adjMatrix=Nmatrix, family=poisson)
  predict(mdfit, newdata=scotlip) # on 56 positions
} 

