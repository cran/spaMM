cat(crayon::yellow("\nAdditional tests for predict():"))

# test of bug corrected in version 2.6.0 (cf NEWS.Rd: "Bug in predict()'ing with an unusual combination of random effects...")
npos <- c(11,16,14,2,6,1,1,4,10,22,7,1,0,0,1,6)
ntot <- c(36,20,19,16,17,11,5,6,37,32,19,17,12,10,9,7)
treatment <- c(rep(1,8),rep(0,8))
clinic <-c(seq(8),seq(8))
clinics <- data.frame(npos=npos,nneg=ntot-npos,treatment=treatment,clinic=clinic)
fitobject <- HLfit(cbind(npos,nneg)~treatment+(1|clinic),family=binomial(),data=clinics,HLmethod="ML")
p1 <- predict(fitobject)[c(4,3),]
p2 <- predict(fitobject, newdata=rbind(clinics[c(4,3),3:4],c(1,9)))[1:2,] 
## The newdata have a >1-subset (not one, not all) of the original levels of the 'clinic' factor.
## For un-correlated random effects, earlier versions returned an incorrect result in this case.
crit <- diff(range(p1-p2))
if (spaMM.getOption("fpot_tol")>0) {
  testthat::test_that(paste0("test-more-predict: criterion was ",signif(crit,4)," >1e-12"), testthat::expect_true(crit<1e-12)) 
} else testthat::expect_true(crit<1e-12)

                    
