cat("\ntest of prediction variance:")
# examples from Booth & Hobert 1998 JASA

npos <- c(11,16,14,2,6,1,1,4,10,22,7,1,0,0,1,6)
ntot <- c(36,20,19,16,17,11,5,6,37,32,19,17,12,10,9,7)
treatment <- c(rep(1,8),rep(0,8))
clinic <-c(seq(8),seq(8))
clinics <- data.frame(npos=npos,nneg=ntot-npos,treatment=treatment,clinic=clinic)
fitobject <- HLfit(cbind(npos,nneg)~treatment+(1|clinic),family=binomial(),data=clinics,HLmethod="ML")
res <- sqrt(attr(predict(fitobject,variances=list(linPred=TRUE,dispVar=TRUE,fixed=TRUE)),"predVar"))
expect_equal(res[6],c("6"=0.7316108),tolerance=2e-4)

if(require("rsae", quietly = TRUE)) {
  data(landsat)
  fitobject <- HLfit(HACorn ~ PixelsCorn + PixelsSoybeans + (1|CountyName),data=landsat[-33,],HLmethod="ML")
  newXandZ <- unique(data.frame(PixelsCorn=landsat$MeanPixelsCorn,PixelsSoybeans=landsat$MeanPixelsSoybeans,CountyName=landsat$CountyName))
  res <- attr(predict(fitobject,newdata=newXandZ,variances = list(linPred=TRUE)),"predVar")
  expect_equal(res[12],c("32"=26.56215),tolerance=2e-4)
  res <- attr(predict(fitobject,newdata=newXandZ,variances = list(linPred=TRUE,dispVar=TRUE)),"predVar")
  expect_equal(res[12],c("32"=27.80684),tolerance=2e-4)
} else {
  cat( "package 'rsae' not available, cannot run prediction variance test.\n" )
}


