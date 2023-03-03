cat(crayon::yellow("\ntest dyadic fixed-effect constructs:\n"))

if (spaMM.getOption("example_maxtime")>0.7) {
  
  #### Simulate dyadic data
  
  set.seed(123)
  nind <- 10       # Beware data grow as O(nind^2)
  x <- runif(nind^2) 
  id12 <- expand.grid(id1=seq(nind),id2=seq(nind))
  id1 <- id12$id1
  id2 <- id12$id2
  u <-  rnorm(nind,mean = 0, sd=0.5)
  
  ## additive individual effects:
  y <-  0.1 + 1*x + u[id1] +  u[id2] + rnorm(nind^2,sd=0.2)
  
  ## anti-smmetric individual effects:
  t <-  0.1 + 1*x + u[id1] - u[id2] + rnorm(nind^2,sd=0.2)
  
  dyaddf <- data.frame(x=x, y=y, t=t, id1=id1,id2=id2, fa1=as.factor(id1), fa2=as.factor(id2))
  # : note that this contains two rows per dyad, which avoids identifiability issues.
  
  # Enforce that interactions are between distinct individuals (not essential for the fit):
  dyaddf <- dyaddf[- seq.int(1L,nind^2,nind+1L),] 
  
  # scramble the data so that input factors are in no partiular order
  set.seed(123)
  dyaddf <- dyaddf[sample(nrow(dyaddf)),]
  
  
  # Fits:
  
  (addfiti <- fitme(y ~x +X.GCA(id1:id2), data=dyaddf))
  (addfitf <- fitme(y ~x +X.GCA(fa1:fa2), data=dyaddf))
  foo <- rev(2:4)
  p1 <- predict(addfiti)[foo]
  testthat::test_that("predict X.GCA  consistent with fitme(y ~x +GCA(fa1:,fa2), data=dyaddf) using lmDiallel:GCA():",
                      testthat::expect_true(diff(range(p1-c(1.7466509, 0.1610820, 0.3324041)))<1e-7))
  p2 <- predict(addfitf)[foo]
  (p3 <- predict(addfiti, newdata=dyaddf[foo,]))
  p4 <- predict(addfitf, newdata=dyaddf[foo,])
  testthat::test_that("predict X.GCA  OK wrt permutations and factor coding",
                      testthat::expect_true(diff(range(p1-p2,p1-p3,p1-p4))<1e-8))
  
  (mvcheck <- fitmv(list(y ~x +X.GCA(fa1:fa2),y ~x +X.GCA(fa1:fa2)), data=dyaddf))
  predict(mvcheck, newdata=dyaddf[foo,])
  
  (antifiti <- fitme(t ~x +X.antisym(id1:id2), data=dyaddf))
  (antifitf <- fitme(t ~x +X.antisym(fa1:fa2), data=dyaddf))
  (p1 <- predict(antifiti)[foo])
  p2 <- predict(antifitf)[foo]
 ( p3 <- predict(antifiti, newdata=dyaddf[foo,]))
  p4 <- predict(antifitf, newdata=dyaddf[foo,])
  testthat::test_that("predict X.antisym  OK wrt permutations and factor coding",
                      testthat::expect_true(diff(range(p1-p2,p1-p3,p1-p4))<1e-8))
  
  if (file.exists(privtest <- paste0(spaMM::projpath(),"/package/tests_other_pack/test-lmDiallel.R"))) {
    source(privtest) 
  }
  
}


