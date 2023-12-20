cat(crayon::yellow("\nAdditional tests for predict():"))

# test of bug corrected in version 2.6.0 (cf NEWS.Rd: "Bug in predict()'ing with an unusual combination of random effects...")
data("clinics")
fitobject <- HLfit(cbind(npos,nneg)~treatment+(1|clinic),family=binomial(),data=clinics,HLmethod="ML")
p1 <- predict(fitobject)[c(4,3),]
p2 <- predict(fitobject, newdata=rbind(clinics[c(4,3),3:4],c(1,9)))[1:2,] 
## The newdata have a >1-subset (not one, not all) of the original levels of the 'clinic' factor.
## For un-correlated random effects, earlier versions returned an incorrect result in this case.
crit <- diff(range(p1-p2))
if (spaMM.getOption("fpot_tol")>0) {
  testthat::test_that(paste0("test-more-predict: criterion was ",signif(crit,4)," >1e-12"), testthat::expect_true(crit<1e-12)) 
} else testthat::expect_true(crit<1e-12)

if (spaMM.getOption("example_maxtime")>5) {
  ## Parallelisation of predict.HLfit()
  library(spaMM)
  data("Loaloa")
  
  foo <- fitme(cbind(npos,ntot-npos)~1 +Matern(1|longitude+latitude),
               data=Loaloa, family=binomial(), fixed=list(lambda=2.602,nu=0.4085199,rho=0.9484547))
  str1 <- tail(capture.output(str(
    predict(foo, blockSize = 10, newdata=Loaloa))),6) # serial
  str2 <- tail(capture.output(str(
    predict(foo, blockSize = 10, newdata=Loaloa, cluster_args=list(spec=3)))),6) # PSOCK cluster
  identical(str1,str2)
  if (.Platform$OS.type !="windows") {
    str3 <- tail(capture.output(str(
      predict(foo, blockSize = 10, newdata=Loaloa, cluster_args=list(spec=3, type="FORK")))),6)
    identical(str1,str3)
  } else message("\nTest of FORK parallelisation cannot be run on Windows.")
}
