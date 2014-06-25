sampleNearby <-
function(focalPts,n=NULL,stepsizes) {
  d <- ncol(focalPts)
  nr <- nrow(focalPts)
  if (n<nr) {
    subfocal <- focalPts[sample(seq_len(nr),n),,drop=FALSE]
  } else subfocal <- focalPts
  randpts <- apply(subfocal,1,function(v) {
    rdxy <- runif(d)* sample(c(-1,1),d,replace=TRUE)* stepsizes/2
    v + rdxy 
  })
  if (d==1L) {
    randpts <- as.matrix(randpts) 
  } else randpts <- t(randpts)
  colnames(randpts) <- colnames(focalPts) ## (in 1D at least) rbind checks that the names match...
  randpts
}
