\name{pedigree}
\alias{pedigree}
\title{Fit mixed-effects models incorporating pedigrees}
\description{This illustrates how to use spaMM for quantitative genetic analyses. spaMM appears competitive in terms of speed for GLMMs with large data sets, particularly when using the PQL method, which may be a quite good approximation in such cases. For large pedigrees it may be useful to compute the inverse of the relationship matrix using some efficient ad hoc algorithm, then to provide it as argument of the fit using the \code{\link{covStruct}(list(precision=...))} syntax.}
\examples{
\dontrun{
if(require("pedigreemm", quietly=TRUE)) {
## derived from help("pedigreemm")
  p1 <- new("pedigree",
           sire = as.integer(c(NA,NA,1, 1,4,5)),
           dam  = as.integer(c(NA,NA,2,NA,3,2)),
           label = as.character(1:6))
  A <- pedigreemm::getA(p1) ## relationship matrix
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
  fitme(y ~ 1+ corrMatrix(1|IDgen) , corrMatrix=A,data=obs,method="REML")
  obs$y01 <- ifelse(y<1.3,0,1)
  fitme(y01 ~ 1+ corrMatrix(1|IDgen)+(1|IDenv), corrMatrix=A,data=obs, 
        family=binomial(), method="REML")
  prec_mat <- solve(A)
  colnames(prec_mat) <- rownames(prec_mat) <- rownames(A) # important
  fitme(y01 ~ 1+ corrMatrix(1|IDgen)+(1|IDenv) , covStruct=list(precision=prec_mat),
        data=obs,family=binomial(), method="REML")
}
}
}

