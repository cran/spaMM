cat(crayon::yellow("\ntest multinomial:"))

if (spaMM.getOption("example_maxtime")>1.7) {
  # extended from example in multinomial.Rd
  
  # good for detecting when one has forgotten to use setProcessed...
  
  set.seed(123)
  genecopy1 <- sample(4,size=50,prob=c(1/2,1/4,1/8,1/8),replace=TRUE)
  genecopy2 <- sample(4,size=50,prob=c(1/2,1/4,1/8,1/8),replace=TRUE)
  alleles <- c("122","124","126","128")
  
  genoInSpace <- data.frame(type1=alleles[genecopy1],type2=alleles[genecopy2],x=runif(50),y=runif(50))
  ## Fitting distinct variances of random effects for each binomial response
  method <- "PQL"
  multifit <- corrHLfit(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
                        family=multi(responses=c("type1","type2")),
                        ranFix=list(rho=1,nu=0.5), method=method)
  if (packageVersion("spaMM")>"3.5.59") {
    abyss <- capture.output(summary(multifit)) # to capture display bugs
    abyss <- capture.output(how(multifit)) # to capture how() bugs
    testthat::expect_true(length(unique(unlist(lapply(multifit, get_ranPars, which="lambda"))))>1L)
  } 
  multifit <- fitme(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
                    family=multi(responses=c("type1","type2")),
                    fixed=list(rho=1,nu=0.5), init=list(lambda=NaN), method=method)
  if (packageVersion("spaMM")>"3.5.59") {
    abyss <- capture.output(summary(multifit)) # to capture display bugs
    testthat::expect_true(length(unique(unlist(lapply(multifit, get_ranPars, which="lambda"))))>1L)
  } 
  ## Fitting the same variance for all binomial responses           
  multifit <- fitme(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
                    family=multi(responses=c("type1","type2")),
                    fixed=list(rho=1,nu=0.5), method=method)
  if (packageVersion("spaMM")>"3.5.59") {
    abyss <- capture.output(summary(multifit)) # to capture display bugs
    testthat::expect_true(length(unique(unlist(lapply(multifit, get_ranPars, which="lambda"))))==1L)
  } 
  multifit <- corrHLfit(cbind(npos,nneg)~1+Matern(1|x+y),data=genoInSpace, 
                        family=multi(responses=c("type1","type2")),
                        ranFix=list(rho=1,nu=0.5), init.corrHLfit=list(lambda=1), method=method)
  if (packageVersion("spaMM")>"3.5.59") {
    abyss <- capture.output(summary(multifit)) # to capture display bugs
    testthat::expect_true(length(unique(unlist(lapply(multifit, get_ranPars, which="lambda"))))==1L)
  } 
}
