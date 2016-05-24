## y* = X beta + Z b + e
## b :=autocorrelated ranefs, =Lv for iid v, = Uw for i but not id w in CAR
SEMbetalambda <- function(beta_eta,lambda,
                          corr_est=NULL,
                          lambda.Fix=NA, 
                          nSEMiter=200, # must be >9 cf check of non-default value in preprocess
                          ngibbs=20,
                          nMCint=NULL, 
                          SEMseed=NULL,SEMsample=NULL,
                          symSVD, ## with a $dim member, not simply a "dim" attribute
                          SEMlogL=if (nrow(X.pv)>250L) {"pMVN"} else {"pmvnorm"}, ## integration of dim nobs, not n_u_h
                          control_pmvnorm=list(maxpts=nMCint),  ## 'for development purposes, not documented...'
                          control_pMVN=list(nrep=nMCint),  ## 'not documented...'
                          ZA,X.pv,qr.XtX,
                          ZAL,whichy1,off,
                          stop.on.error,
                          verbose=FALSE,
                          X_lamres,
                          control.glm,
                          mc){ 
  if (FALSE) {
    adaptive <- TRUE
    nSEMiter <- 10000L
    pilotLength <- 50L
  } else {
    adaptive <- FALSE
    pilotLength <- nSEMiter
  }
  if (is.null(SEMsample)) SEMsample <- ceiling(nSEMiter/2):nSEMiter
  if(!is.null(SEMseed)) { ## explicit value of SEMseed argument
    set.seed(SEMseed) ## so that estimates of beta,lambda are repeatable ## comment ne pas avoir a retirer les memes nombres XXXX fois ?
    #      cat(paste("SEMseed=",SEMseed))
  } # else HLCorobj call
  nobs <- nrow(X.pv)
  betaMat <- matrix(0,nrow=nSEMiter,ncol=ncol(X.pv))
  colnames(betaMat) <- colnames(X.pv) 
  betaMat[1,] <- beta_eta
  lambdaVec <- numeric(nSEMiter)
  condVar <- rep(0,nSEMiter)
  condVar[1] <- lambdaVec[1] <- lambda
  corr.model <- symSVD$corr.model  
  estimCARrho <- ( corr.model=="adjacency" &&  ! is.null(corr_est$rho)) 
  if (estimCARrho) {
    rhoVec <- numeric(nSEMiter)
    rhoVec[1] <- corr_est$rho
  }
  n_u_h <- symSVD$dim[2L] ## locally true as implied by code below.
  EbGivenY <- matrix(0,nrow=nSEMiter,ncol=n_u_h)
  decomp <- symSVD$symsvd
  if ( ! is.identity(ZA)) {
    ZAisI <- FALSE
    ZA <- as.matrix(ZA) ## useless to stay Matrix
    ZAE <- ZA %id*id% decomp$u ## dense in spatial model  
  } else ZAisI <- TRUE
  whichy0 <- (! whichy1) ##FR->FR in preprocess ?
  ny1 <- sum(whichy1) ##FR->FR in preprocess ?
  ny0 <- nobs - ny1
  if (SEMlogL == "pmvnorm") {
    pmvnorm.lower <- rep(0,nobs)
    pmvnorm.lower[whichy0] <- -Inf
    pmvnorm.upper <- rep(0,nobs)
    pmvnorm.upper[whichy1] <- Inf
  }
  gibbsSample <- ceiling(ngibbs/2):ngibbs
  # prevmsglength <- 0L
  i <- 2
  while (TRUE) {
    ## whatever depends on rhoVec[i-1] (fixed in the gibbs block)
    if ( estimCARrho ) {
      decomp$d <- 1/(1-rhoVec[i-1]*decomp$adjd) ## le decomp$d originel est Corr, pas pour L
    }
    if (i==2L || estimCARrho) { 
      if (corr.model=="identity") {
        LMatrix <- invLMatrix <- decomp$u
      } else {
        LMatrix <- ZWZt(decomp$u,sqrt(decomp$d))
        invLMatrix <- ZWZt(decomp$u,1/sqrt(decomp$d)) ## no wrapper (testing identity is a loss of time)
      }
      if( ! ZAisI) {
        ZAEdEAZ <- ZWZt(ZAE,decomp$d) ## FR->FR large (nresp * nresp)
        forV <- selfAdjointSolverCpp(ZAEdEAZ) ## so that inv(V) = ZWZt(forV$u,1/(1+lambda_est * forV$d))
        ##             [ D = lambda Corr]  . Z'      . forV$u but without the lambda factor
        if (corr.model=="identity") {
          ranefCorr <- decomp$u
        } else {
          ranefCorr <- tcrossprodCpp(LMatrix) ## reconstruction compliqu?e de la corr matrix ?
        }
        LHSCorrblob <- as.matrix(ranefCorr %id*% t(ZA) %*% forV$u)
      }
    }
    ## whatever (also) depends on lambdaVec[i-1] (fixed in the gibbs block)
    if(ZAisI) {
      CondNorm <- CondNormfn(decomp,lambdaVec[i-1])
    } else { ## D - D Z' inv(V) Z D
      ##          [D = lambda *Corr] - [LHSblob= lambda Corr Z' forV$u]. 1/(1+lambda_est * forV$d) . t(LHSblob) 
      ## with lambda^2 /(1+lambda d) = lambda/(1/lambda + d)
      condCov <- lambdaVec[i-1] * ranefCorr - ZWZt(LHSCorrblob,lambdaVec[i-1]/(1/lambdaVec[i-1] + forV$d)) ## ZWZt must be a slow step
      condL <- RcppChol(as.matrix(condCov))$L ## such that only tcrossprod(condL) = tcrossprod(tcrossprod(condL)) when ZAisI
      ## not more code because I will try to perform only matrix * vector operations
    }
    # S part of the SEM algorithm
    # we use a Gibbs sampling algorithm
    if (i==2) { #!# test added 04/2016
      randbGivenObs <- sqrt(lambdaVec[i-1]) * (LMatrix %*% rnorm(n_u_h,0))
    } 
    rand_augY <- rep(0,nobs)
    fix <- X.pv %*% betaMat[i-1,] + off
    # estimates, not 'used' in the gibbs loop, but modified by it.  
    if (estimCARrho) {
      sum_condw2 <- rep(0,n_u_h)      
    } else sum_condv2 <- rep(0,n_u_h)
    sumcond_bMeans <- rep(0,n_u_h) ## initial value of sum, 
    sum_rand_augY_Means <- rep(0,nobs) ## initial value of sum,
    gibbs_subarglist <- list(fix=fix, y=as.numeric(whichy1), 
                          ZAisI=ZAisI, n_u_h=n_u_h, lambda=lambdaVec[i-1],
                          eps=.Machine$double.eps )
    if (ZAisI) {
      gibbs_subarglist <- c(gibbs_subarglist,list(CondNorm=CondNorm,ZA=NULL,forV=NULL,LHSCorrblob=NULL,condL=NULL))
    } else {
      gibbs_subarglist <- c(gibbs_subarglist,list(CondNorm=NULL,ZA=ZA,forV=forV,LHSCorrblob=LHSCorrblob,condL=condL))
    }
    gibbs_arglist <- c(gibbs_subarglist,list(ngibbs=ngibbs, 
                                             gibbsSample=gibbsSample-1, ## C iterator vs R iterator... 
                                             estimCARrho=estimCARrho,
                                          decomp=decomp))
    if (TRUE) { 
      gibbsblob <- do.call("Rcpp_gibbs",gibbs_arglist)
      sumcond_bMeans <- gibbsblob$sumcond_bMeans
      sum_rand_augY_Means <- gibbsblob$sum_rand_augY_Means
      if (estimCARrho) {
        sum_condw2 <- gibbsblob$sum_condv2_or_w2
      } else sum_condv2 <- gibbsblob$sum_condv2_or_w2
    } else { ## pure R code kept for debugging
      for (k in 1:ngibbs) { ## versions < 1.8.30 have a partial C++ implementation 
        # random generation of augY given obs: y and v (fixed beta, fixed lambda)
        ## for Matrix ZA, ZA %*% vector is a considerable waste of time, avoided by ZA %id*%
        rand_Esp_augY <- fix + ZA %id*% randbGivenObs  
        rand_augY[whichy1] <- rntpos(ny1,rand_Esp_augY[whichy1],1)
        rand_augY[whichy0] <- rntneg(ny0,rand_Esp_augY[whichy0],1)
        ## whatever depends on augY
        if(ZAisI) {
          randcond_bMean <- CondNorm$condLvReg %*% (rand_augY-fix) # (dim n_u_h) <- (dim nobs)
          randcond_brand <- CondNorm$sqrtCondCovLv %*% rnorm(n_u_h,0)
        } else { ##randcond_bMean <- lambda_est * (ranefCorr %*% (t(ZA) %*% solve(augYCov,augY-fix))) ## DZ'inv(V)(y-X beta) in Searle p. 275
          ##             [D = lambda *Corr]   .     Z'    .   inv(V).(augY-fix)   with initial lambda brought inside
          ##randcond_bMean <- ranefCorr %*% t(t(forV$u %*% t((t(augY-fix) %*% forV$u)/(1/lambdaVec[i-1] + forV$d))) %*% ZA) ## only t(vector) ## SLOW
          ##randcond_bMean <- LHSCorrblob %*% t((t(augY-fix) %*% forV$u)/(1/lambdaVec[i-1] + forV$d)) ## only t(vector)
          # faster code for the same:
          locv <- rand_augY-fix
          dim(locv) <- c(1,nrow(locv)) ## fast transposition
          locv <- (locv %*% forV$u)/(1/lambdaVec[i-1] + forV$d)
          dim(locv) <- c(ncol(locv),1) ## fast transposition
          randcond_bMean <- LHSCorrblob %*% locv
          randcond_brand <- condL %*% rnorm(n_u_h,0) ## has zero expectation
        }
        ## augY should be fix + randbGivenObs + one-epsilon-per-individual 
        # random generation of v given (y and) augmented Y 
        randbGivenObs <- randcond_bMean + randcond_brand ## to compute condw2/condv2 and to initiate next iteration
        #if (k==2) print(c(i,k,mean(randcond_bMean),mean(randbGivenObs)))
        if (k %in% gibbsSample) { ## computes sums rather than store a matrix
          sumcond_bMeans <- sumcond_bMeans + randcond_bMean ## sufficient for E[b|obs] since randcond_brand has zero expectation
          sum_rand_augY_Means <- sum_rand_augY_Means + rand_augY
          if (estimCARrho) {
            sum_condw2 <- sum_condw2 + (t(decomp$u) %*% randbGivenObs)^2 ## w indep mais pas i.d., pour estim lambda et rho
          } else sum_condv2 <-  sum_condv2 + (invLMatrix %*% randbGivenObs)^2 ## v iid
        }
      } ## end ngibbs loop
    }
    ### end of E step
    ### M step of the SEM algorithm
    # determination of beta by standard least square #betaMat[i,] <- lm((z-Lvs)~X.pv-1)$coeff
    EbGivenY[i,] <- sumcond_bMeans/length(gibbsSample) ##  E[b|...]
    EaugY <- sum_rand_augY_Means/length(gibbsSample)
    if (ZAisI) {
      betaMat[i,] <- solveWrap.vector( qr.XtX , t(X.pv) %*% (EaugY- EbGivenY[i,] -off) ,stop.on.error=stop.on.error)
    } else betaMat[i,] <- solveWrap.vector( qr.XtX , t(X.pv) %*% (EaugY-(ZA %id*% EbGivenY[i,]) -off) ,stop.on.error=stop.on.error)
    if ( estimCARrho ) {
      condw2 <- sum_condw2/length(gibbsSample)
    } else condv2 <- sum_condv2/length(gibbsSample) ## in output...
    # Estim of (co)variance parameters
    if (is.na(lambda.Fix)) {
      if ( estimCARrho ) {
        locdf <- data.frame(v2=condw2,adjd=decomp$adjd) ## FR->FR hmf. adjd as part of an X_lamres computed by preprocess ?
        ## these two functions are (equivalently) maximized by the GLM:
        #         bordel <- function(lr) {
        #           covmat <- diag(lr[1]/(1-lr[2]*decomp$adjd))
        #           dmvnorm(as.numeric(vs), mean = rep(0, nrow(covmat)), sigma = covmat, log = TRUE)
        #         }
        #         bordel2 <- function(lr) {
        #           covmat <- lr[1] * ZWZt(decomp$u,1/(1-lr[2]*decomp$adjd)) ## pour randbGivenObs
        #           dmvnorm(as.numeric(randbGivenObs), mean = rep(0, nrow(covmat)), sigma = covmat, log = TRUE)
        #         }        
        #print(fixef(HLfit(v2 ~ 1+ adjd,family=spaMM_Gamma("inverse"),prior.weights=rep(1/2,n_u_h),data=locdf)))
        resglm <- spaMM_glm(v2 ~ 1+ adjd, family=spaMM_Gamma("inverse"), strict=TRUE,
                            weights=rep(1/2,n_u_h), data=locdf,
                            start=c(1/mean(locdf$v2),0)) ## is ML (no EQL issue)
        coeffs <- coefficients(resglm) 
        lambdaVec[i] <- 1/coeffs[1]
        rhoVec[i] <- - coeffs[2]/coeffs[1] ## note that lambdas[k]= mean((1-rhos[k]*decomp$adjd)*vs^2)
        if (verbose) prevmsglength <- overcat(paste("iter", i,": lambda=",signif(lambdaVec[i],4),"; rho=",signif(rhoVec[i],4)), prevmsglength)
      } else {
        lambdaVec[i] <- mean(condv2)
        if (verbose) prevmsglength <- overcat(paste("iter",i,": lambda=",signif(lambdaVec[i],4)), prevmsglength)
      }
    } else {
      if ( estimCARrho ) {
        locdf <- data.frame(v2=condw2,adjd=decomp$adjd) ## FR->FR hmf. adjd as part of an X_lamres computed by preprocess ?
        resglm <- spaMM_glm(v2 ~ adjd -1, offset=rep(1/lambda.Fix,nrow(locdf)),
                            family=spaMM_Gamma("inverse"), strict=TRUE,
                            weights=rep(1/2,n_u_h), data=locdf, start=c(0)) ## is ML (no EQL issue)
        coeffs <- coefficients(resglm) 
        lambdaVec[i] <- lambda.Fix
        rhoVec[i] <- - coeffs[1]*lambda.Fix 
        if (verbose) prevmsglength <- overcat(paste("iter", i,": rho=",signif(rhoVec[i],4)), prevmsglength)
      } else lambdaVec[i] <- lambda.Fix
    }
    if (adaptive && i==pilotLength) {} ## see memo in adaptiveDiagnostics.R.txt
    if (nSEMiter==i) break ## either !adaptive or (adaptive && nSEMiter=pilotLength=i)
    i <- i+1
  }  ## end nSEMiter loop
  if (verbose) {
    cat("\n") 
    plot(lambdaVec[1:nSEMiter])
    if (estimCARrho) plot(rhoVec)
  }
  
  if (FALSE) {
    # Some diagnostics
    adjsample <- seq(min(SEMsample),max(SEMsample),length.out = min(200,max(SEMsample)-min(SEMsample)))
    adj <- diag(0,length(adjsample))
    diag(adj[-1,]) <- diag(adj[,-1]) <- 1
    beta1frame <- data.frame(beta1=betaMat[adjsample,1L],id=seq(length(adjsample)))
    adjfit <- HLCor(beta1 ~1+adjacency(1|id),adjMatrix=adj,data=beta1frame,HLmethod="ML")
    effdf <- var(beta1frame$beta1)/vcov(adjfit)[1] ## effective sample size
    if (effdf<30) message(paste("Suggested chain length: >",ceiling(nSEMiter*30/effdf)))
  }
  #
  if (estimCARrho) {
    corr_est$rho <- mean(rhoVec[SEMsample])
    decomp$d <- 1/(1-corr_est$rho*decomp$adjd) 
    LMatrix <- ZWZt(decomp$u,sqrt(decomp$d))
    invLMatrix <- ZWZt(decomp$u,1/sqrt(decomp$d))
    if (ZAisI) {ZAL <- LMatrix} else { ZAL <- ZA %id*id% LMatrix }
  }
  #
  logLargs <- list(betaMat=betaMat, lambdaVec=lambdaVec, 
                   SEMsample=SEMsample, ZAL=ZAL, nobs=nobs, X.pv=X.pv, off=off, 
                   SEMlogL=SEMlogL,
                   pMVN_arglist = c(control_pMVN, 
                                    list(limits = as.vector(fix),ismax=whichy1))
  )
  if (SEMlogL=="pmvnorm") {
    pmvnorm_arglist <- c(control_pmvnorm,
                                  list(lower=pmvnorm.lower,upper=pmvnorm.upper,mean=as.vector(fix),
                                       abseps=0,releps=0.5))
    # mvt.f shows that the algo continues (within some limits) 
    #  if the error is higher than the MAX of the abseps and releps-induced constraints. 
    # For large negative logL, we can only control releps, and we must set set abseps=0 to control releps.
    if (is.null(pmvnorm_arglist$maxpts)) ##  postprocessing in corrHLfit provides a maxpts
      pmvnorm_arglist$maxpts <- quote(1000L*nobs) ## "sensible ... to start with MAXPTS = 1000*N" (mvt.f)
    ## but the default in optim through smooth is 250L*nobs
    logLargs$pmvnorm_arglist <- pmvnorm_arglist
  }
  #
  if (verbose) cat("Estimating the likelihood...")
  resu <- do.call("logL_from_SEMsample",logLargs)
  lambda <- resu$lambda ## Don't forget this 2015/09/16
  
  if (verbose) cat(" SEM done.\n")
  if (estimCARrho) {
    ## to provide a conforming glm to HLfit but also coefficients to be substituted to the coefficents of this glm
    locdf$resp <- condw2
    resu$glm_lambda <- calc_CARdispGammaGLM(data=locdf, lambda.Fix=lambda.Fix,control=control.glm)
    resu$glm_lambda$coeffs_substitute <- c("(Intercept)"= 1/lambda,"adjd"= - corr_est[["rho"]]/lambda)
  } else {
    ## calc_dispGammaGLM needs 'data' even if 'data' has zero cols  
    Intcol <- which(colnames(X_lamres)=="(Intercept)")
    loclamdata <- data.frame(X_lamres[, - Intcol,drop=FALSE])
    ## FR->FR in a "later" version the formula and 'data' should be provided by processed$ as X_lamres come from there
    ## FR->FR could be tested now ?
    if (ncol(loclamdata)>0L) {
      loclamformula <- as.formula(paste("~",paste(colnames(loclamdata),collapse="+")))
    } else loclamformula <- ~1
    resu$glm_lambda <- calc_dispGammaGLM(
      formula=loclamformula,
      dev.res=as.matrix(condv2),
      lev=rep(0,n_u_h),
      data=loclamdata,
      etastart=rep(log(lambda),n_u_h), 
      control=control.glm)
    resu$glm_lambda$coeffs_substitute <- c("(Intercept)"= log(lambda))
  }
  ## output for prediction only: 
  resu$v_h <- invLMatrix %*% colMeans(EbGivenY[SEMsample,,drop=FALSE]) ## v iid from b not iid
  return(resu)
} ## end def of SEMbetalambda

