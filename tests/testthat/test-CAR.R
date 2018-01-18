cat("\ntest CAR and SEM:\n")

data("scotlip")

# see also 'old donttest' examples

## same without optim: run in scotlip examples; cf also autoregressive.Rd for ML fits

set.seed(124)
eigenv <- sym_eigen(Nmatrix) ## could use eigen(,symmetric=TRUE)
Lmat <- eigenv$u %*% diag(sqrt(1/(1-0.17*eigenv$d)))
lp <- 0.1 + 3* Lmat %*% rnorm(ncol(Lmat)) ## single intercept beta =0.1; lambda=3
resp <- rbinom(ncol(Lmat),1,1/(1+exp(-lp)))
donn <- data.frame(npos=resp,nneg=1-resp,gridcode=scotlip$gridcode)

if (FALSE) { ## HLCor/corrHlfit already compared on scotlip by test-spaMM.R
  # corrHLfit without corners was poor here
  # CAR by Laplace with 'outer' estimation of rho
  # *** fitme is not very convicing, stops early ***
  blob <- fitme(cbind(npos,nneg)~1 +adjacency(1|gridcode),
                    adjMatrix=Nmatrix,family=binomial(probit),data=donn,method="ML",control.HLfit = list(LevenbergM=FALSE)) 
  AIC(blob)
}

# CAR by Laplace with 'inner' estimation of rho
blob <- HLCor(cbind(npos,nneg)~1 +adjacency(1|gridcode),
          adjMatrix=Nmatrix,family=binomial(probit),data=donn,HLmethod="ML",control.HLfit = list(LevenbergM=FALSE)) ## 2 s.
AIC(blob)



if (require("probitgem",quietly=TRUE)) {
  # CAR by SEMs +optimsmooth ... slow
  if (spaMM.getOption("example_maxtime")>38) {
    blob <- corrHLfit(cbind(npos,nneg)~1 +adjacency(1|gridcode),
            adjMatrix=Nmatrix,family=binomial(probit),data=donn,HLmethod="SEM")
  }
  
  # CAR by single SEM
  blob <- HLCor(cbind(npos,nneg)~1 +adjacency(1|gridcode),
                adjMatrix=Nmatrix,family=binomial(probit),data=donn,HLmethod="SEM") ## -1 s, with more variance of logLik
  AIC(blob)
} else {
  cat( "package 'probitgem' not available, cannot run SEM test.\n" )
}

