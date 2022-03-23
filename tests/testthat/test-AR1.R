cat(crayon::yellow("\ntest AR1:\n"))

if (spaMM.getOption("example_maxtime")>4) { # 
  set.seed(123)
  nobs <- 500
  distm <- as.matrix(dist(1:nobs)) 
  m <- (-0.4)^distm 
  cholm <- t(chol(m))
  eta <- 1+ cholm %*% rnorm(nobs) ## hglm_lambda=1
  obs <- rpois(nobs,exp(eta))
  plot(obs)
  fake <- data.frame(obs=obs,age=1:nobs)
  fitar1 <- fitme(obs ~ 1+AR1(1|age),family=poisson(),data=fake,verbose=c(TRACE=0.5),method="REML")
  crit <- diff(range(logLik(fitar1), c(p_bv=-1269.060222885349)))
  try(testthat::test_that(paste0("criterion was ",signif(crit,4)," from 1269.060222885349"), testthat::expect_true(crit<1e-9))) 
  # There should be better tests elsewhere:
  # fitar1_cF <- fitme(obs ~ 1+corrFamily(1|age),family=poisson(),data=fake,verbose=c(TRACE=0.5),
  #                    covStruct=list("1"=ARp()), method="REML")
  # crit <- diff(range(logLik(fitar1_cF),logLik(fitar1)))
  # try(testthat::test_that(paste0("logLik(fitar1_cF) was ",signif(crit,4)," logLik(fitar1)"), testthat::expect_true(crit<1e-9))) 
}

## same with nested AR1 within individual
## Large data necess for good estimation of ARphi and other params 

# quick version for routine tests
if (TRUE) {
  rngcheck <- ("sample.kind" %in% names (formals(RNGkind)))
  if (rngcheck) suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"  )) ## necessary for the test on range(c(p1,p2,p3, 2.35717079935))
  set.seed(123)
  age    <-    (rep(c(1:((Nage <- 30))),times=(Nind <- 30)))
  ind    <- 	 rep(c(1:Nind),	 each=(Nage))
  distm <- as.matrix(dist(age)) 
  blocks <- proxy::dist(ind,`==`)
  blocks[blocks==0] <- NA
  distm <- distm * as.matrix(blocks) ## nice hack
  distm[is.na(distm)] <- 1e100 ## temporary hack
  m <- (0.4)^distm  ## intriguingly negative ARphi seems much easier to estimate than positive ones.
  cholm <- t(chol(m))
  eta <- 1+ cholm %*% rnorm(Nage*Nind) ## encore lambda=1
  obs <- rpois(Nind*Nage,exp(eta))
  plot(obs)
  fake <- data.frame(obs=obs,age=age,ind=as.factor(ind+1L), ## as.factor( [all > 1] ) to test the dark side of uniqueGeo
                     idx=seq_len(length(obs)))
  ## the sample() calls provides a check that permutations of the data have no effect
  ## checks the sparse->non-sparse case (assuming .determine_spprec() returns FALSE)
  zut <- corrHLfit(obs ~ 1+AR1(1|idx %in% ind),family=poisson(),data=fake[20+sample(20),],ranFix=list(ARphi=0.7040234,lambda=0.7308))
  rezut <- corrHLfit(obs ~ 1+AR1(1|idx %in% ind),family=poisson(),data=fake[20+sample(20),],ranFix=list(ARphi=0.7040234,lambda=0.7308), 
                     control.HLfit=list(sparse_precision=TRUE))
  rerezut <- corrHLfit(obs ~ 1+AR1(1|idx %in% ind),family=poisson(),data=fake[20+sample(20),],ranFix=list(ARphi=0.7040234,lambda=0.7308), 
                       control.HLfit=list(sparse_precision=FALSE))
  # The data are permuted between each fit, which could contributed to (in principle trivial) differences among fits
  testthat::expect_true(diff(range((c(logLik(zut),logLik(rezut),logLik(rerezut),-47.3130016607291))))<1e-8)
  ## check predict on each fit and subset of (permuted) data:
  p1 <- predict(zut,newdata=rezut$data[rownames(rezut$data)>30,])["39"] 
  p2 <- predict(rezut,newdata=rerezut$data[rownames(rerezut$data)>30,])["39"]
  p3 <- predict(rerezut,newdata=zut$data[rownames(zut$data)>30,])["39"]
  crit <- diff(range(c(p1,p2,p3, 2.35717079935)))## last decimals sensitive to d_relV_b_tol
  if (spaMM.getOption("fpot_tol")>0) {
    testthat::test_that(paste0("criterion was ",signif(crit,6)," from 2.35717079935"), testthat::expect_true(crit<1e-10)) 
  } else testthat::expect_true(crit<1e-10)
  if (rngcheck) RNGkind("Mersenne-Twister", "Inversion", "Rejection"  )
}


if (spaMM.getOption("example_maxtime")>6) {
  set.seed(123)
  age    <-    (rep(c(1:((Nage <- 30))),times=(Nind <- 30)))
  ind    <- 	 rep(c(1:Nind),	 each=(Nage))
  distm <- as.matrix(dist(age)) 
  blocks <- proxy::dist(ind,`==`)
  blocks[blocks==0] <- NA
  distm <- distm * as.matrix(blocks) ## nice hack
  distm[is.na(distm)] <- 1e100 ## temporary hack
  m <- (0.4)^distm  ## intriguingly negative ARphi seems much easier to estimate than positive ones.
  cholm <- t(chol(m))
  eta <- 1+ cholm %*% rnorm(Nage*Nind) ## encore lambda=1
  obs <- rpois(Nind*Nage,exp(eta))
  plot(obs)
  
  fake <- data.frame(obs=obs,age=age,ind=ind,idx=seq_len(length(obs)))
  ## the sample() provides a check that permutations of the data have no effect
  ## checks the sparse->non-sparse case
  (zut <- corrHLfit(obs ~ 1+AR1(1|idx %in% ind),family=poisson(),data=fake[20+sample(20),]))
  (rezut <- corrHLfit(obs ~ 1+AR1(1|idx %in% ind),family=poisson(),data=fake[20+sample(20),], 
                      control.HLfit=list(sparse_precision=TRUE)))
  rerezut <- corrHLfit(obs ~ 1+AR1(1|idx %in% ind),family=poisson(),data=fake[20+sample(20),], 
                       control.HLfit=list(sparse_precision=FALSE))
  crit <- diff(range(c(logLik(zut),logLik(rezut),logLik(rerezut))))
  if (spaMM.getOption("fpot_tol")>0) {
    testthat::test_that(paste0("criterion was ",signif(crit,6)," from  -47.31300"), testthat::expect_true(crit<1e-8) )
  } else testthat::expect_true(crit<1e-8)
  ## full data
  fit_ar1nested <- corrHLfit(obs ~ 1+AR1(1|age %in% ind),family=poisson(),data=fake,verbose=c(TRACE=interactive())) 
  testthat::expect_equal(logLik(fit_ar1nested), c(p_bv=-2295.67792783))
}

if (spaMM.getOption("example_maxtime")>0.5) {
  requireNamespace("nlme")
  data("Orthodont",package = "nlme")
  if (TRUE) { # fitme has (finally) become as fast as corrHLfit on this example
    checkinput <- fitme(distance ~ age + factor(Sex)+( 1 | Subject)+ AR1(1|age %in% Subject), fixed=list(phi=1e-6), 
                        data = Orthodont,method="REML")  
  } else { 
    checkinput <- corrHLfit(distance ~ age + factor(Sex)+( 1 | Subject)+ AR1(1|age %in% Subject), ranFix=list(phi=1e-6),
                            data = Orthodont,HLmethod="REML")
  }
  testthat::expect_equal(logLik(checkinput), c(p_bv=-218.69839984))
  # consistent with 
  # lme(distance ~ age + factor(Sex),random = ~ 1 | Subject, cor=corCAR1(form=~age|Subject),data = Orthodont)
  # which is faster (FIXME: .assign_geoinfo_and_LMatrices_but_ranCoefs() for AR1 not efficient; more work needed to handle nested AR1 efficiently)
}