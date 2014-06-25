HLfit <-
function(formula,
                  data,family=gaussian(),rand.family=gaussian(), 
                  resid.formula = ~ 1 ,REMLformula=NULL,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),
                  HLmethod="HL(1,1)",
                  control.HLfit=list(),
                  init.HLfit = list(), 
                  ranFix=list(), ## phi, lambda, possibly nu, rho if not in init.HLfit
                  etaFix=list(), ## beta, v_h (or even u_h)
                  prior.weights= rep(1,nobs),
                  processed=NULL
) {

  #####################################################################
  # local fn defs
  # attention au piege vite oublié
  # locfn1 <- fn() {... use global, e.g. mu}
  # locfn2 <- fn() {... modif mu; locfn1()}
  # => locfn2->locfn1-> lit mu global pas local a locfn2
  #####################################################################
  multiHLfit <- function() {
    fitlist <- lapply(data,function(dt){
      locmc <- mc
      if (family$family=="multi") {
        locmc$family <- family$binfamily
      }
      locmc$data <- dt
      eval(locmc)
    })
    liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
    liks<- apply(liks,1,sum)
    attr(fitlist,"APHLs") <- as.list(liks)
    attr(fitlist,"sortedTypes") <- attr(data,"sortedTypes")
    attr(fitlist,"responses") <- attr(data,"responses")
    class(fitlist) <- c("HLfitlist",class(fitlist))     
    return(fitlist)
  }
  #####################################################################
  #####################################################################
  
  
  
  mc <- match.call() ## ## potentially used by getCall(object) in update.HL... if HLfit was called by HLCor through a do.call() this contains the body of the function 
  ## Pour resoudre le probleme de memoire (mais pas du programmeur): 
  ## In that case HLCor removes this from the HLfit object and gives its own call. Otherwise we can improve a bit by 
  ## mc[[1]] <-  call("HLfit")[[1]] ## replace the body with the call; eval(mc) will still work
  ## but all other arguments are still evaluated... cf HLCor
  if (missing(data)) data <- environment(formula)
  if (family$family=="multi") {
    if ( ! inherits(data,"list")) {
      if(family$binfamily$family=="binomial") {
        familyargs <- family
        familyargs$family <- NULL
        familyargs$binfamily <- NULL
        data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }
    }
  }    
  ################# data LIST ##############################################
  if ( inherits(data,"list")) return(multiHLfit())
  ##########################################################################
  corrNames <- intersect(c("nu","rho","Nugget","ARphi"),names(init.HLfit))
  if (length(corrNames)>0) {
    corr_est <- init.HLfit[corrNames]
  } else corr_est <- NULL
  ###################################################
  useSparseM <- FALSE ## TRUE (and require(SparseM) and LevenbergM=FALSE) 
  ##   allows a quick test of a sparse Cholesky alternative to QR computations.
  ##   Not conclusive because wAugX not very sparse
  warningList<-list()
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["warn"])) verbose["warn"] <- TRUE
  if (is.na(verbose["summary"])) verbose["summary"] <- FALSE
  ##
  if (family$family %in% c("poisson","binomial")) {
    phi.Fix<-1 
  } else {
    phi.Fix <- ranFix$phi
    if (any(phi.Fix==0)) {
      mess <- pastefrom("phi cannot be fixed to 0.",prefix="(!) From ")
      stop(mess)
    }
  } ## immediately used in preprocess call:
  if (is.null(processed)) {
    validrows <- validRows(formula=formula,resid.formula=resid.formula,data=data) ## will remove rows with NA's in required variables
    data <- data[validrows,,drop=FALSE] ## ## before Predictor is called and an LMatrix is added, etc.
    processed <- preprocess(control.HLfit=control.HLfit,HLmethod=HLmethod,predictor=formula,phi.Fix=phi.Fix,
                            resid.predictor=resid.formula,REMLformula=REMLformula,data=data,family=family,
                            BinomialDen=BinomialDen,rand.families=rand.family)
  } 
  predictor <- processed$predictor
  rand.families <- processed$rand.families
  lcrandfamfam <- processed$lcrandfamfam
  LevenbergM <- processed$LevenbergM
  betaFirst <- processed$betaFirst ## avoids explicitly solving as an aug.lin.mod.
  HL <- processed$HL
  if (HL[1]=="SEM") {
    SEMargs <- processed[names(processed) %in% c("SEMseed","nMCint","nSEMiter","ngibbs","SAEMsample")]
  }
  stop.on.error <- processed$stop.on.error ## to control issues with matrix computations; F by default
  AIC <- processed$AIC ## whether to compute any AIC stuff; F by default
  essai <- processed$essai ## to control any tested new code...
  version <- processed$version ## to replicate some selected features of old versions
  useSparseQR <- processed$useSparseQR ## to control any tested new code...
  conv.threshold <- processed$conv.threshold
  iter.mean.dispFix <- processed$iter.mean.dispFix
  iter.mean.dispVar <- processed$iter.mean.dispVar
  max.iter <- processed$max.iter
  #    ps_v.threshold <- processed$ps_v.threshold
  resid.predictor <- processed$resid.predictor 
  BinomialDen <- processed$BinomialDen
  loglfn.fix <- processed$loglfn.fix
  y <- processed$y
  REMLformula <- processed$REMLformula
  X.Re <- processed$X.Re
  X.pv <- processed$X.pv
  ### a bit of post processing
  nobs <- NROW(X.pv)
  pforpv <- ncol(X.pv)
  if ( ! is.null(REMLformula) && (ncol(X.Re) != pforpv)) { ## differences affects only REML estimation of dispersion params, ie which p_bv is computed
    distinct.X.ReML <- TRUE ## true in the ML case [ncol(X.Re)=0] if pforpv>0
  } else {
    distinct.X.ReML <- FALSE ## the X of REML is the standard one  
  }
  ###
  canonicalLink <- processed$canonicalLink
  LMMbool <- processed$LMMbool
  models <- processed$models
  #### Note that HLCor modifies the L matrix (inprocessed$predictor if required) => ZAL cannot be preprocessed by corHLfit and must be recomputed each time 
  if (models[["eta"]]=="etaHGLM") { ## Design matriCES for random effects in particular, prob only a match between the levels or the ranef and the observ. Ie Z, not ZAL 
    lambda.family <- processed$lambdaFamily
    LMatrix <- attr(predictor,"LMatrix")
    ZAlist <- processed$ZAlist ## : ZAlist is a list of design matrices 
    Groupings <- attr(ZAlist,"Groupings")
    ZAL <- attr(predictor,"ZALMatrix")
    if ( is.null(ZAL)) { ## reconstruct ZAL from Z (Z from spMMFactorList, L from user)
      ZALlist <- compute.ZALlist(LMatrix=LMatrix,ZAlist=ZAlist,Groupings=Groupings)
    } ## end case reconstruction of ZAL from ZA and L
  } else { ## models[["eta"]] = "etaGLM"
    ZALlist <- NULL
    u_h <- v_h <- lev_lambda <- numeric(0)
  } 
  ### a bit of post processing // repeat of code in preprocess...
  nrand <- length(ZALlist)
  lambda.Fix <- ranFix$lambda
  if (any(lambda.Fix==0)) {
    mess <- pastefrom("lambda cannot be fixed to 0.",prefix="(!) From ")
    stop(mess)
  }
  vec_n_u_h <- rep(0, nrand)
  for (i in 1:nrand) vec_n_u_h[i] <- ncol(ZALlist[[i]]) ## nb cols each design matrix = nb realizations each ranef
  cum_n_u_h <- cumsum(c(0, vec_n_u_h)) ## if two ranef,  with q=(3,3), this is 0,3,6. cum_n_u_h[nrand+1] is then 6, the total # of realizations
  ZAL <- post.process.ZALlist(ZALlist,predictor=predictor)
  I <- diag(rep(1,cum_n_u_h[nrand+1L]))
  ZALI <- rbind(ZAL,I) 
  ###
  X_lamres <- processed$X_lamres
  X_disp <- processed$X_disp ## may be NULL
  off <- attr(processed$predictor,"offset")
  ##################
  if (is.character(init.HLfit)) { ## at this point init.HLfit is a string or not. Elsewhere it can be a list
    spaMM.options(INIT.HLFITNAME=init.HLfit) ## if a string, copied in...
  } else {
    spaMM.options(INIT.HLFITNAME=NA)  
    # init.HLfitName <- NULL
    unknowns <- names(init.HLfit)[!names(init.HLfit) %in% c("fixef","phi","lambda","v_h","rho","nu","Nugget","ARphi")] 
    if (length(unknowns)>0) {
      mess <- pastefrom("unhandled elements in 'init.HLfit'.",prefix="(!) From ")
      message(mess)
      if ("beta" %in% unknowns) message("  Use 'fixef' rather than 'beta' in 'init.HLfit'.")
      stop()
    }
  }
  ###################
  if ( ! is.null(corr_est)) {
    ## ici on veut une procedure iterative sur les params de covariance
    #  HLCor.args$processed <- processed ## FR->FR dangerous in early development
    corrEst.args <- list(family=family,rand.family=rand.families) ## but rand.families must only involve a single spatial effect 
    loc.oriform <- attr(predictor,"oriFormula")
    loc.lhs <- paste(loc.oriform)[[2]]
    ## build formula with only spatial effects
    corrEst.form <-  as.formula(paste(loc.lhs," ~ ",paste(findSpatial(loc.oriform))))
    corrEst.args$data <- data ## FR->FR way to use preprocess ???                    
    if (ncol(X.Re)>0) { ## some REML correction
      if (distinct.X.ReML) {
        corrEst.args$REMLformula <- REMLformula
      } else corrEst.args$REMLformula <- predictor ## REML without an explicit formula
      corrEst.args$objective <- "p_bv"
    } else corrEst.args$objective <- "p_v" 
    corrEst.args$ranFix <- ranFix ## maybe not very useful
    corrEst.args$control.corrHLfit$Optimizer<- control.HLfit$Optimizer ## (may be NULL => L-BFGS-B) 
    corrEst.args$control.corrHLfit$optim$control$maxit <- 1 
    corrEst.args$control.corrHLfit$nlminb$control$iter.max <- 2 ## 1 => convergence trop lente
    corrEst.args$control.corrHLfit$optimize$tol <- 1e10 
    corrEst.args$control.corrHLfit$corners <- FALSE ## 
  }
  
  
  
  #################### MORE LOCAL FNS DEFS ###################################
  ##mais processed controle le default nMCint
  ## all per-iteration stats are taken from gibbsSample
  ## and all final stats are the means, from iterations SAEMsample, of the per-iteration stats 
  SEMbetalambda <- function(beta_eta,nSEMiter=100,ngibbs=20,nMCint=10000,SEMseed=NULL,SAEMsample=NULL){ ## beta_eta as explicit argument so that the iterative aspect is explicit
    if (nSEMiter<10) stop("(!) In 'SEMbetalambda', 'nSEMiter' should be >9")
    if (is.null(SAEMsample)) SAEMsample <- (nSEMiter/2):nSEMiter
    if(!is.null(SEMseed)) set.seed(SEMseed) ## comment ne pas avoir a retirer les memes nombres XXXX fois ?
    betaMat <- matrix(0,nrow=nSEMiter,ncol=ncol(X.pv))
    EuGivenY=matrix(0,nrow=nSEMiter,ncol=length(y))
    lambdaVec <- numeric(nSEMiter)
    betaMat[1,] <- beta_eta
    lambdaVec[1] <- unique(lambda_est)
    condVar=rep(0,nSEMiter)
    condVar[1] <- lambdaVec[1]
    ZA <- ZAlist[[1]] ## FR->FR ad hoc
    if ( ! attr(ZA,"identityMatrix")) {
      stop("! attr(ZA,'identityMatrix'): more code needed in SEM algo") ## CondNormf not adequate
      ZAisI <- FALSE
    } else ZAisI <- TRUE
    whichy1 <- (y==1) ##FR->FR in preprocess ?
    whichy0 <- (! whichy1) ##FR->FR in preprocess ?
    ny1 <- sum(whichy1) ##FR->FR in preprocess ?
    ny0 <- sum(whichy0) ##FR->FR in preprocess ?
    if (ny0+ny1 != nrow(X.pv)) {
      stop("(!) SEM procedure: the data do not seem binary; other binomial data are not handled.")
    }
    decomp <- attr(LMatrix,attr(LMatrix,"type"))
    ## whatever does not depend on lambda:
    if(ZAisI) {
      invLMatrix <- ZWZt(decomp$u,1/sqrt(decomp$d))
    } else {
      ranefCorr <- tcrossprodCpp(LMatrix) ## FR->FR not elegant
      ZAE <- ZA %*% decomp$u  
      ZAEdEAZ <- ZWZt(ZAE,decomp$d)
      forV <- selfAdjointSolverCpp(ZAEdEAZ) ## so that inv(V) = ZWZt(forV$u,1/(1+lambda_est * forV$d))
      ##             [ D = lambda Corr]  . Z'      . forV$u but without the lambda factor
      LHSCorrblob <- ranefCorr %*% t(ZA) %*% forV$u
    }
    gibbsSample <- (ngibbs/2):ngibbs 
    for (i in 2:nSEMiter)
    {
      ## whatever depends on lambdaVec[i-1] (fixed in the gibbs block)
      if(ZAisI) {
        CondNorm <- CondNormfn(LMatrix,lambdaVec[i-1])
      } else { ## D - D Z' inv(V) Z D
        ##          [D = lambda *Corr] - [LHSblob= lambda Corr Z' forV$u]. 1/(1+lambda_est * forV$d) . t(LHSblob) 
        ## with lambda^2 /(1+lambda d) = lambda/(1/lambda + d)
        condCov <- lambdaVec[i-1] * ranefCorr - ZWZt(LHSCorrblob,lambdaVec[i-1]/(1/lambdaVec[i-1] + forV$d))
        condL <- RcppChol(condCov)$L ## such that only tcrossprod(condL) = tcrossprod(tcrossprod(condL)) when ZAisI
        ## not more code because I will try to perform only matrix * vector operations
      }
      # S part of the SEM algorithm
      # random generation of z and v given y
      # we use a Gibbs sampling algorithm
      rvGivenObs <- sqrt(lambdaVec[i-1]) * (ZAL %*% rnorm(n_u_h,0))
      augY <- rep(0,nrow(ZAL))
      fix <- X.pv %*% betaMat[i-1,] + off
      lambdas <- numeric(ngibbs)
      condMeans <- matrix(0,nrow=ngibbs,ncol=n_u_h)
      for (k in 1:ngibbs) {
        # random generation of augY given obs: y and v (fixed beta, fixed lambda)
        moy.augY <- fix + rvGivenObs  
        augY[whichy1] <- rntpos(ny1,moy.augY[whichy1],1)
        augY[whichy0] <- rntneg(ny0,moy.augY[whichy0],1)
        ## whatever depends on augY
        if(ZAisI) {
          condMean <- CondNorm$condLvReg %*% (augY-fix)
          condv <- CondNorm$sqrtCondCovLv %*% rnorm(n_u_h,0)
        } else { ##condMean <- lambda_est * (ranefCorr %*% (t(ZA) %*% solve(augYCov,augY-fix))) ## DZ'inv(V)(y-X beta) in Searle p. 275
          ##             [D = lambda *Corr]   .     Z'    .   inv(V).(augY-fix)   with initial lambda brought inside
          condMean <- ranefCorr %*% t(t(forV$u %*% t((t(augY-fix) %*% forV$u)/(1/lambdaVec[i-1] + forV$d))) %*% ZA) ## only t(vector)
          condv <- condL %*% rnorm(n_u_h,0)
        }
        ## augY should be fix + rvGivenObs + one-epsilon-per-individual 
        # random generation of v given (y and) augmented Y 
        rvGivenObs <- condMean + condv
        condMeans[k,] <- condMean
        lambdas[k] <- mean((invLMatrix %*% rvGivenObs)^2)
      } ## end ngibbs loop
      EuGivenY[i,] <- colMeans(condMeans[gibbsSample,])
      # M part of the SEM algorithm
      # determination of beta by standard least square #betaMat[i,] <- lm((z-Lvs)~X.pv-1)$coeff
      betaMat[i,] <- solveWrap.vector( qr.XtX , t(X.pv) %*% (augY-EuGivenY[i,]-off) ,stop.on.error=stop.on.error) 
      if (is.null(lambda.Fix)) {
        lambdaVec[i] <- mean(lambdas[gibbsSample])
        #if(ZAisI) {
        #  CondNorm <- CondNormfn(LMatrix,lambdaVec[i])
        #} 
      } else lambdaVec[i] <- lambda.Fix
    }  ## end nSEMiter loop
    beta_eta <- colMeans(betaMat[SAEMsample,,drop=FALSE]) ## SAEMsample no longer useful ?
    lambda_est <- mean(lambdaVec[SAEMsample])
    ## simulate final likelihood(with high variance...)
    fix <- X.pv %*% beta_eta + off
    binLikcond <- numeric(nrow(ZAL))
    logtotLikcond <- numeric(nMCint)
    Lik <- 0
    for (it in 1:nMCint) {
      eta <- fix + sqrt(lambda_est) * (ZAL %*% rnorm(n_u_h,0))
      mu <- family$linkinv(eta) ## ou pnorm()
      binLikcond <- (1-mu)* (whichy0)+ mu*(whichy1)
      logtotLikcond[it] <- sum(log(binLikcond))
    }
    maxlogtotLikcond <- max(logtotLikcond)
    logLik <- log(mean(exp(logtotLikcond-maxlogtotLikcond)))+maxlogtotLikcond
    return(list(beta_eta=beta_eta,lambda_est=lambda_est,v_h=colMeans(EuGivenY[SAEMsample,]),p_v=logLik))
  } ## end local def of SEMbetalambda
  
  u_h_from_v_h <- function(v_h) {
    anyinf <- FALSE
    u_list <- lapply(seq(nrand), function(it){
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      vh <- v_h[u.range]
      uh <- rand.families[[it]]$linkinv(vh)
      if (any(is.infinite(uh))) {
        anyinf <<- TRUE
        if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="log") {
          ## Gamma(log)$linkinv is pmax(exp(eta), .Machine$double.eps), ensuring that all gamma deviates are >= .Machine$double.eps
          ## we ensure that log(u_h) has symmetric bounds on log scale (redefine Gamma()$linkfun ?)
          uh <- pmin(uh,1/.Machine$double.eps)
          vh <- rand.families[[it]]$linkfun(uh)
        } else {
          mess <- pastefrom("infinite 'u_h'.",prefix="(!) From ") 
          warning(mess)
        }
      } 
      attr(uh,"vh") <- vh 
      return(uh) ## wait for problems to happen...
    })
    u_h <- unlist(u_list)
    if (anyinf) attr(u_h,"v_h") <- unlist(lapply(u_list, function(v){attr(v,"vh")}))
    return(u_h)
  }

  `calc.p_v` <- function(returnLad=FALSE,only.h=FALSE) { 
    theta <- theta.mu.canonical(mu/BinomialDen,family$family)  
    if (family$family=="binomial") {
      clik <- sum(loglfn.fix(theta,y/BinomialDen,BinomialDen,1/(phi_est))) ## freq a revoir
    } else {
      #browser()
      clik <- sum(loglfn.fix(theta,y,prior.weights/phi_est)) ## note (prior) weights meaningful only for gauss/ Gamma 
    }
    if (models[[1]]!="etaHGLM" && models[3]!="phiHGLM") return(list(clik=clik,p_v=clik))
    ## ELSE
    likranU <- unlist(lapply(seq(nrand), function(it) {
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      loglfn.ranU(lcrandfamfam[it],u_h[u.range],1/lambda_est[u.range])    
    }))
    #if(returnLad) browser()
    log.du_dv <- - log(dvdu) 
    likranV <- sum(likranU + log.du_dv)
    ##### HLIK
    hlik <- clik+likranV      
    if (only.h) return(list(hlik=hlik))
    ##### P_V
    if (.SpaMM$USEEIGEN) {
      lad <- LogAbsDetCpp(d2hdv2/(2*pi))
    } else lad <-determinant(d2hdv2/(2*pi))$modulus[1] 
    if (is.nan(lad) || is.infinite(lad)){## because of determinant of nearly singular matrix
      zut <- abs(eigen(d2hdv2/(2*pi),only.values = T)$values) ## 05/01/13
      zut[zut<1e-12] <- 1e-12
      lad <- sum(log(zut)) ## L-BFGS-B requires a non trivial value
    }
    p_v <- hlik-lad/2
    resu <- list(clik=clik,hlik=hlik,p_v=p_v)
    if(returnLad) resu$lad <- lad
    if (HL[1]==2) { ## uses second-order correction as in LeeN01, p.996
      if (family$family!="binomial") {
        stop("HL(2,.) not implemented for non-binomial models")
      } 
      muFREQS <- as.numeric(mu/BinomialDen)
      d3bTh <- BinomialDen * muFREQS*(1-muFREQS)*(1-2*muFREQS) ## b'''(th) ## BinomialDen * muFREQS*(1-muFREQS) doit etre les blob$GLMweights
      d4bTh <- BinomialDen * muFREQS*(1-muFREQS)*(1-6*muFREQS+6*muFREQS^2) ##b''''(th)=D[Log[1 + E^th], {th, 4}] /. {E^th -> p/(1 - p), E^(k_ th) :> (p/(1 - p))^k}
      ### ## code based on the notation of Appendix B of LeeL12 (last page of supp mat)   
      ### initial version with qr. Slow. QR not useful here. 
      #    if (is.null(qr.d2hdv2)) qr.d2hdv2 <- qr(d2hdv2) ## BUT best computed before the calc.p_v call, indeed will be needed for next HL(1) in *G*LMM, therefore not locally in this function
      #    B <- try(solve.qr(qr.d2hdv2),silent=TRUE)   ## best for the sandwiches below would be a factorization CCt of this inverse... but directly!
      #    if (class(B)=="try-error") {
      #      second.corr <- - exp(700) ## drastically penalizes the likelihood when d2hdv2 is nearly singular 
      #    } else { ## code based on the notation of Appendix B of LeeL12 (last page of supp mat)
      #      B <- - B
      #      ZBZt <- ZAL %*% B %*% t(ZAL) ## we need the full matrix for the last term...
      #      diagZBZt <- diag(ZBZt)
      #      HabcdBacBbd <- - sum(d4bTh*diagZBZt^2) ## minus added 040213 ...
      #      b3.diagZBZT <- as.numeric(diagZBZt * d3bTh) ## inner product vector = vector
      #      acoefs <- as.numeric(b3.diagZBZT %*% ZAL) ## vector which 'a'th element is Sum_i d3bTh_i * (ZBZt)_ii ZAL_ia
      #      HabcHrstBarBbcBst <- acoefs %*% B %*% acoefs
      #      ZBZtcube <- ZBZt * ZBZt * ZBZt
      #      HabcHrstBarBbsBct <- d3bTh %*% ZBZtcube %*% d3bTh
      #      second.corr <- HabcdBacBbd/8 + HabcHrstBarBbcBst/8 + HabcHrstBarBbsBct/12 ## - F/24
      #    }
      L <- try(t(chol( - d2hdv2)),silent=TRUE) 
      if (class(L)=="try-error") { ## but fall back 
        L <- designL.from.Corr( - d2hdv2,try.chol=F) ## recycles code for 'square root'; but no longer triangular
        if (class(L)=="try-error") {
          second.corr <- - exp(700) ## drastically penalizes the likelihood when d2hdv2 is nearly singular 
        } else { ## code based on the notation of Appendix B of LeeL12 (last page of supp mat)
          invL.ZALt <- solve(L,t(ZAL)) ## L^{-1}.t(ZAL)
          ZBZt <- crossprod(invL.ZALt) ## ZAL.t(L)^{-1}.L^{-1}.t(ZAL)
          diagZBZt <- diag(ZBZt)
          HabcdBacBbd <- - sum(d4bTh*diagZBZt^2) ## minus sign from habcd =- b''''(th) zzzz
          b3.diagZBZT <- as.numeric(diagZBZt * d3bTh) ## inner product vector = vector
          acoefs <- as.numeric(b3.diagZBZT %*% ZAL) ## vector which 'a'th element is Sum_i d3bTh_i * (ZBZt)_ii ZAL_ia
          HabcHrstBarBbcBst <- sum(solve(L,acoefs)^2) ## strictly, sum(solve(L, - acoefs)^2) starting from  habc=-b'''(th) zzz
          ZBZtcube <- ZBZt * ZBZt * ZBZt
          HabcHrstBarBbsBct <- d3bTh %*% ZBZtcube %*% d3bTh ## again, - d3bTh %*% ZBZtcube %*% (-d3bTh)
          second.corr <-  HabcdBacBbd/8 + HabcHrstBarBbcBst/8 + HabcHrstBarBbsBct/12 ## - F/24 (ps_bv =p_bv + second.corr = p_bv -F/24)
        }    
      } else { ## clearly the fastest code
        invL.ZALt <- forwardsolve(L,t(ZAL)) ## L^{-1}.t(ZAL)
        ZBZt <- crossprod(invL.ZALt) ## ZAL.t(L)^{-1}.L^{-1}.t(ZAL)
        diagZBZt <- diag(ZBZt)
        HabcdBacBbd <- - sum(d4bTh*diagZBZt^2)
        b3.diagZBZT <- as.numeric(diagZBZt * d3bTh) ## inner product vector = vector
        acoefs <- as.numeric(b3.diagZBZT %*% ZAL) ## vector which 'a'th element is Sum_i d3bTh_i * (ZBZt)_ii ZAL_ia
        HabcHrstBarBbcBst <- sum(forwardsolve(L,acoefs)^2)
        ZBZtcube <- ZBZt * ZBZt * ZBZt  ## I tried things with sweep() but there is always a triple inner product 
        HabcHrstBarBbsBct <- d3bTh %*% ZBZtcube %*% d3bTh ## again, - d3bTh %*% ZBZtcube %*% (-d3bTh)
        second.corr <- HabcdBacBbd/8 + HabcHrstBarBbcBst/8 + HabcHrstBarBbsBct/12 ## - F/24
      }
      resu$ga <- HabcdBacBbd/8
      resu$bu <- HabcHrstBarBbcBst/8
      resu$zo <- HabcHrstBarBbsBct/12
      ## ## second.corr <- log(1+second.corr)
      resu$second.corr <- second.corr
      resu$p_v <- p_v + second.corr
      
      ## print(c(p_v,ps_v))
    }
    return(resu)
  }
  
  
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  ## syntax check on etaFix$beta (12/2013)
  if ( length(etaFix$beta)>0 ) {
    if ( length(etaFix$beta)!=ncol(X.pv) ) {
      message("(!) An incomplete etaFix$beta vector was provided.")
      message("  This is highly dubious. If you want to fix some elements and fix others")
      message("  It is recommended to use a restricted model formula plus an offset.")
      stop("    I exit.")
    } else {
      ## correct length, but this won't be taken into account if the elemnts are not named
      if (is.null(names(etaFix$beta))) {
        message("(!) The elements of etaFix$beta should be named and the names should match the column names of the design matrix.")
        stop("    I exit.")
      }
    }
  } 


  ### case where nothing to fit #############################################
  if (is.null(corr_est) && 
        length(etaFix$beta)==ncol(X.pv) &&
        !is.null(phi.Fix) &&
        (models[[1]]=="etaGLM" || (!is.null(etaFix$v_h) &&  !is.null(lambda.Fix))) 
      ) { ## nothing to fit. We just want a likelihood
    ### a bit the same as max.iter<1 ... ?
    phi_est <- phi.Fix
    if ( ! is.null(etaFix$beta) ) { ## can be false if the whole of the fixed part is in the offset 
      eta <- off + X.pv %*% etaFix$beta
    } else eta <- off
    if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
      ## we need u_h in calc.p_v() and v_h here for eta...
      v_h <- etaFix$v_h
      u_h <- etaFix$u_h
      if (is.null(u_h)) {u_h <- u_h_from_v_h(v_h)}
      lambda_est <- lambda.Fix
      eta <- eta + ZAL %*% etaFix$v_h ## updated at each iteration
    } ## FREQS
    ## conversion to mean of response variable (COUNTS for binomial)
    blob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen,priorWeights=prior.weights) 
    mu <- blob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
    w.resid<-as.vector(blob$GLMweights /phi_est) ## 'weinu', must be O(n) in all cases
    if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
      wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## no fit, likelihood computation
      dvdu <- wranefblob$dvdu
      w.ranef <- wranefblob$w.ranef
      d2hdv2 <- calcDhDv2(ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    }
    return(list(APHLs=calc.p_v())) ### RETURN !!!!!!!!! ## FR->FR but p_bv is not returned...
  } 
  
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  
  `provide.resglm` <- function() { 
    if (family$family=="binomial" && ncol(y)==1) { 
      ##  && ncol(y)==1: attempt to implement the cbind() for y itself syntax throughout. But fails later on 'y - mu'...
      begform <-"cbind(y,BinomialDen-y)~"  
    } else {begform <-"y~"}
    ###################################################if (pforpv==0) {endform <-"0"} else 
    if(pforpv>0) {
      endform <-"X.pv-1" ## pas besoin de rajouter une constante vue qu'elle est deja dans X
    } else {
      if (family$family %in% c("binomial","poisson")) {
        endform <- "1" ## no meaningful glm without fixed effect in this case !
      } else {endform <- "0"}
    }
    locform <- as.formula(paste(begform, endform))
    resglm<-glm(locform,family=family,offset=off,weights=prior.weights) ## warning "glm.fit: fitted probabilities numerically 0 or 1 occurred" => separation or large offset
    if (pforpv>0) {
      if (max(abs(c(resglm$coefficients)),na.rm=TRUE)>1e10) { ## na.rm v1.2 pour param non estimables (cas normal)
        message("(!) Apparent divergence of estimates in a *glm* analysis of the data.")
        message("    Check your data for separation or bad offset values.")
        stop("    I exit.") 
      } else names(resglm$coefficients) <- colnames(X.pv) ## because this if unfortunately not the case... 
    } 
    return(resglm)
  }

  generateInitLambda <- function() {
    if (is.null(lambda.Fix)) { 
      init.lambda <- init.HLfit$lambda
      if (is.null(init.lambda) ) {
        #### initial values for lambda
        # first rough estimate of lambda assuming a single rand.family=gaussian(identity)
        # then distribute the variation over the different rand families
        # then account for non gaussian(id) rand families
        # (1)
        if (family$family=="binomial" && family$link=="logit") {
          fv<-fitted(resglm)
          init.lambda <- sum((resid(resglm)^2)/(resglm$prior.weights*fv*(1-fv)))/resglm$df.residual
        } else {
          resdisp <-as.numeric(deviance(resglm)/resglm$df.residual) 
          if (family$family=="poisson" && family$link=="log") {
            init.lambda <- pmax(0.00001,log(resdisp))
          } else init.lambda <- resdisp/5 ## assume that most of the variance is residual
        } 
        # (2)
        init.lambda <- init.lambda/nrand        
        #
        ## allows for different rand.family
        init.lambda <- unlist(lapply(seq(nrand), function(it) {
          if(lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="log") {
            objfn <- function(lambda) {psigamma(1/lambda,1)-init.lambda}
            adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
          } else if(lcrandfamfam[it]=="beta" && rand.families[[it]]$link=="logit") {
            #ad hoc approximation which should be quite sufficient; otherwise hypergeometric fns.
            objfn <- function(lambda) {8* lambda^2+3.2898*lambda/(1+lambda)-init.lambda}
            adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
          } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") {
            ## (pi^2)/6 is upper bound for expected value
            if (init.lambda > 1.64491 ) { 
              adhoc <- 100000 ## so that psigamma(1+1/100000,1) ~  1.64491
            } else {
              objfn <- function(lambda) {psigamma(1+1/lambda,1)-init.lambda}
              adhoc <- uniroot(objfn,interval=c(1e-8,1e8))$root
            }
          } else if(lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="-1/mu") {
            adhoc <- (sqrt(1+4*init.lambda)-1)/2 # simple exact solution
          } else adhoc <- init.lambda
          adhoc
        }))
      } 
    } else init.lambda <- lambda.Fix
    return(init.lambda)
  }
  
  `compute.sscaled` <- function() { ## needs qr.d2hdv2, ZAL, stop.on.error, d2hdv2, rWW, ZALI, family, dmudeta, BinomialDen, mu, eta... ETC
    ########## HL(1,.) adjustment for mean ################## and specifically the a(1) term in LeeL 12 p. 963
    ## if LMM (ie resp gaussian, ranef gaussian), all coef<x> are 0 -> correction is 0 (but this fn must not have been called)
    ## if (gaussian, not gaussian) d3 nonzero
    ## if (non gaussian, gaussian), d3 zero (!maybe not for all possible cases) but d1,d2 nonzero 
    vecdi1 <- NA; vecdi2 <- NA; vecdi3 <-NA
    if (all(dlogWran_dv_h==0L)) vecdi3 <- 0
    # coef1 is the factor of P_ii in d1
    # coef2 is the factor between P_jj and K1 in d2
    ## We first handle the canonical link cases, where comput. of coef1 depends only on the link  
    ## here w=dmudeta; d1=dwdmu dmudeta /w^2 = dlogwdeta/w = (d2mu/deta2)/(dmu/deta) /w =
    ##      (d2mu/deta2)/(dmu/deta)^2 = (d(dmudeta)/dmu)/dmudeta where d(dmudeta)/dmu is the numerator as detailed:
    if (canonicalLink) {
      if (family$family=="gaussian") {
        vecdi1 <- 0
        vecdi2 <- 0
      } else if (family$family=="poisson") {
        ## numerator is D[D[E^\[Eta], \[Eta]] /. {E^\[Eta] -> \[Mu]}, \[Mu]] =1 
        coef1 <- 1/dmudeta
        coef2 <- rep(1,nobs)
      } else if (family$family=="binomial") {
        ## numerator is D[D[1/(1 + E^-\[Eta]), \[Eta]] /. {E^-\[Eta]->(1-\[Mu])/\[Mu]} ,\[Mu]]=1-2 mu 
        coef1 <-(1-2*mu/BinomialDen)/dmudeta  
        coef2 <-(1-2*mu/BinomialDen)  
      } else if (family$family=="Gamma") { ## link= "inverse" !
        ## numerator is D[D[-1/\[Eta], \[Eta]] /. {\[Eta] -> -1/\[Mu]}, \[Mu]] =2 mu 
        coef1 <- 2*mu /dmudeta
        coef2 <- 2*mu
      } 
    } else if (family$family=="binomial" && family$link=="probit") { ## ad hoc non canonical case 
      muFREQS <- mu/BinomialDen
      coef2 <- -2*eta - dnorm(eta)*(1-2*muFREQS)/(muFREQS*(1-muFREQS))
      dnorm.eta <- dnorm(eta) ## if eta diverges, dnorm -> 0, coef1 -> +/- Inf hence correction:
      # dnorm.eta[dnorm.eta==0] <- .Machine$double.eps
      coef1 <- coef2 *(muFREQS*(1-muFREQS))/ (BinomialDen * dnorm.eta^2) 
      coef1[coef1>1e100] <- 1e100
      coef1[coef1< -1e100] <- -1e100
    } else if (family$family=="Gamma" && family$link=="log") { ## ad hoc non canonical case 
      vecdi1 <- 0
      vecdi2 <- 0 ## because they both involve dW.resid/dmu= 0
    } else {
      ## general code handling non canonical links
      tmblob <- thetaMuDerivs(mu,BinomialDen,family$family)
      Dtheta.Dmu <- tmblob$Dtheta.Dmu # calcul co fn de muFREQS puis / BinomialDen
      D2theta.Dmu2 <- tmblob$D2theta.Dmu2 # calcul co fn de muFREQS puis / BinomialDen ^2
      d2mudeta2 <- d2mudeta2fn(link=family$link,mu=mu,eta=eta,BinomialDen=BinomialDen)
      ## if the family link eta(mu) equals the canonical link theta(mu), then theta=eta, the following line is null  
      D2theta.Deta2_Dtheta.Deta <- (D2theta.Dmu2 * dmudeta^2 + Dtheta.Dmu * d2mudeta2)/(Dtheta.Dmu * dmudeta)
      coef2 <- (d2mudeta2/dmudeta) + D2theta.Deta2_Dtheta.Deta 
      coef1 <- coef2 / (Dtheta.Dmu * dmudeta^2) ## note that coef2 is indep of the BinomialDen, but coef1 depends on it 
    }
    if (any(is.na(c(vecdi1,vecdi2,vecdi3)))) { ## then we need to compute some of them
      ## P is P in LeeL appendix p. 4 and is P_R in MolasL p. 3307 
      ## looks like leverage computation, but the ZALI columns are considered instead of the full augmented design matrix 
      wAugZALI <- sweep(ZALI,MARGIN=1,sqrt.ww,`*`) # rWW %*% ZALI ## rWW previously computed for leverage computation
      ## next two lines will be slow for very large matrices but the leverages() function using RcppEigen is even slower
      partition <- attr(ZAL,"partition")
      if ( !is.null(partition) ) { ## use block structure in ZAL;  
        nrI <- nrow(I)
        Pdiag <- rep(0,nrow(wAugZALI))
        abyss <- sapply(seq_len(length(partition)-1),function(v) {
          sequ <- (partition[v]+1):partition[v+1] 
          colonnes <- wAugZALI[,sequ,drop=F]
          qrcolonnes <- qr(colonnes)
          levs <- rowSums(qr.qy(qrcolonnes, diag(1, nrow = nrow(colonnes), ncol = ncol(colonnes)))[c(sequ,sequ+nrI),,drop=FALSE]^2)
          Pdiag[c(sequ,sequ+nrI)] <<- levs ## proceeds by this side-effect
        })
      } else { ## straightforward but does not use block structure
        qrwAugZALI <- qr(wAugZALI)
        Pdiag <- rowSums(qr.qy(qrwAugZALI, diag(1, nrow = nrow(wAugZALI), ncol = ncol(wAugZALI)))^2)
      }
      #
      seqn_u_h <- seq_len(cum_n_u_h[nrand+1L])
      Pdiagn <- Pdiag[1:nobs]
      if (is.na(vecdi1)) vecdi1 <- Pdiagn * coef1
      # K2 is K2 matrix in LeeL appendix p. 4 and is -D in MolasL p. 3307 
      # W is Sigma^-1 ; TWT = t(ZALI)%*%W%*%ZALI = ZAL'.Wresid.ZAL+Wranef = -d2hdv2 !
      K2 <- solveWrap.matrix(qr.d2hdv2,t(ZAL),stop.on.error=stop.on.error) ## t(ZAL) missing in some implementations... no effect with Gaussian ranefs...  
      if (class(K2)=="try-error") {
        mess <- pastefrom("problem in 'K2' computation.",prefix="(!) From ") ## cf BB 877
        warning(mess)
        K2 <- ginv(d2hdv2) %*% t(ZAL)            
      } ## so far we have computed (d2hdv2)^{-1}.t(Z)= -(TWT)^{-1}.t(Z)
      if (is.na(vecdi2)) {
        # K1 is K1 in LeeL appendix p. 4 and is A=-ZD in MolasL p. 3307-8 
        K1 <- ZAL %*% K2
        vecdi2<-rep(0,nobs) 
        for (i in 1:nobs) vecdi2[i] <- sum( Pdiagn * coef2 * K1[1:nobs,i] )                 
      }
      # coef3 =(1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
      ## This (dlogWran_dv_h) was computed when w.ranef was computed
      if (is.na(vecdi3)) { 
        vecdi3<-rep(0,nobs)
        for (i in 1:nobs) vecdi3[i] <- sum( Pdiag[nobs+seqn_u_h] * dlogWran_dv_h[seqn_u_h] * K2[seqn_u_h , i] ) ## d3 reste nul pour gaussian ranef
      }
      vecdi <- vecdi1+vecdi2+vecdi3 ## k_i in MolasL; le d_i de LeeL app. p. 4
      sscaled <- vecdi /2  ## sscaled := detadmu s_i= detadmu d*dmudeta/2 =d/2 in LeeL12; or dz1 = detadmu (y*-y) = detadmu m_i=0.5 k_i dmudeta = 0.5 k_i in MolasL 
    } else sscaled <-0
    if (any(is.infinite(sscaled))) { ## FR->FR debugging code
      mess <- pastefrom("infinite 'sscaled'.",prefix="(!) From ") 
      stop(mess)
    }
    if (any(is.nan(sscaled))) { ## FR->FR debugging code
      mess <- pastefrom("NaN 'sscaled'.",prefix="(!) From ") 
      stop(mess)
    }
    return(sscaled)
  }
  
  ##############################################################################################
  ######### Initial estimates for mu by GLM ####################
  ## FR->FR UGLY: to fit a GLM, glm() is called here... 
  if ( ( pforpv>0 && is.null(init.HLfit$fixef)) || is.null(phi.Fix) || is.null(init.HLfit$v_h) || is.null(lambda.Fix) ) { 
    ## all cases where an initial resglm is needed (even when pforpv=0, may be needed to provide init phi or init lambda)
    resglm <- provide.resglm()   
  }
  if (pforpv>0) { 
    beta_eta <- numeric(ncol(X.pv))
    beta_eta <- init.HLfit$fixef
    if (is.null(beta_eta) ) {
      beta_eta<-c(resglm$coefficients)[1:pforpv] ## FR->FR peut-on legitimement avoir NA ici ?
    } 
    beta_eta[names(etaFix$beta)] <- etaFix$beta
  } else {
    beta_eta <- numeric(0)
    se_beta <-numeric(0)
  }
  if(FALSE && any(is.na(beta_eta))) { # FR->FR note FALSE in the if... 
    ## FR->FR das preprocess en utilisant lm pour faire le boulot simplement ? 
    ## stop("(!) Some coefficients of the predictor are not estimable. Check that all predictor 'variables' are variable.")
    XpvOri <- X.pv
    pforpvori <- pforpv
    beta_etaOri <- beta_eta
    X.pv <- X.pv[,which(!is.na(beta_eta)),drop=F]
    pforpv <- ncol(X.pv)
    beta_eta <- beta_eta[which(!is.na(beta_eta))]
  } else beta_etaOri <- NULL
  ## Initial estimate for phi ####
  if (is.null(phi.Fix)) { ## at this point, means that not poisson nor binomial
    phi_est <- init.HLfit$phi ## must be a list of 'predictor' values not linear coefficients of predictor 
    if (is.null(phi_est) ) {
      phi_est <- as.numeric(deviance(resglm)/resglm$df.residual)
      if (models[[3]] != "phiScal") {
        phi_est <- rep(phi_est,nobs) ## moche ## why is this necess ?
      }
    } 
  } else {
    phi_est <- phi.Fix
  }
  ##
  ######## initialize v_h #############################
  if (models[[1]]=="etaHGLM") { ## the basic case (LMM, GLMM...)
    if (is.null(LevenbergM)) { ##FR->FR cf comment in preprocess()...
      if (HL[1]=="SEM") {
        LevenbergM <- FALSE 
        betaFirst <- FALSE 
      } else if (betaFirst) {
        LevenbergM <- FALSE 
      } else {
        if (LMMbool) {  
          LevenbergM <- FALSE ## because no reweighting when beta_eta changes => no IRWLS necess   
        } else LevenbergM <- TRUE
      }
    }
    psi_M <- unlist(lapply(seq(nrand), function(it) {
      lpsi <- switch(lcrandfamfam[it], 
                     gaussian = 0,
                     gamma = 1, 
                     beta = 1/2, 
                     "inverse.gamma" = 1
      )
      rep(lpsi,vec_n_u_h[it])
    }))
    v_h <- etaFix$v_h
    if (is.null(v_h) ) v_h <- init.HLfit$v_h
    if (is.null(v_h) ) {
      v_h <- unlist(lapply(seq(nrand), function(it){
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        rand.families[[it]]$linkfun(psi_M[u.range]) ## v as link(mean(u)) 
      }))
    }
    u_h <- u_h_from_v_h(v_h) 
    checkv_h <- attr(u_h,"v_h") ## verif is v_h was corrected 
    if (!is.null(checkv_h)) v_h <- checkv_h
    init.lambda <- generateInitLambda()
    ## one could imagine fixing some of then but not others...  
    if  (length(init.lambda)==nrand) {
      lambda_est <- rep(init.lambda,vec_n_u_h)
    } else if (length(init.lambda)==1) { ## typically what the current default resglm provides even for nrand>1
      lambda_est <- rep(init.lambda,cum_n_u_h[nrand+1L])
    } else if (length(init.lambda)==cum_n_u_h[nrand+1L]) {
      lambda_est <- init.lambda
    } else {stop("Initial lambda cannot be mapped to levels of the random effect(s).")}
  }
  if (models[[3]]=="phiHGLM") {
    stop("random effects in predictor or residual variance (phi) not yet implemented")
    ## there is a buggy template code with comments in version 260812 of HLfit
  }
  #print(paste("ZAL[1,1] in HLfit ",ZAL[1,1]))
  ## predictor from initial values
  if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
    eta <- off + X.pv %*% beta_eta + ZAL %*% v_h ## updated at each iteration
  } else  eta <- off + X.pv %*% beta_eta ## no iteration hence no updating  ## FREQS
  ## conversion to mean of response variable (COUNTS for binomial)
  blob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen,priorWeights=prior.weights) 
  mu <- blob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  dmudeta<- blob$dmudeta ## if Bin/Pois, must be O(n)
  Vmu <- blob$Vmu ## if Bin/Pois, O(n)
  w.resid<-as.vector(blob$GLMweights /phi_est) ## 'weinu', must be O(n) in all cases
  
  if (models[[1]]=="etaHGLM") {
    wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## initilization !
    w.ranef <- wranefblob$w.ranef
    dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
    dvdu <- wranefblob$dvdu
  }
  betaV <- c(beta_eta,v_h) 
  conv.phi <- FALSE; conv.lambda <- FALSE; conv.corr <- FALSE
  if (models[[1]]=="etaHGLM") {
    Sig <- ZWZt(ZAL,1/w.ranef) + diag(1/w.resid) ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid) 
    d2hdv2 <- calcDhDv2(ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    #Sig <- sweep(ZAL,MARGIN=2,1/w.ranef,`*`)  %*% t(ZAL) + diag(1/w.resid) ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid) 
    #d2hdv2 <- - sweep(t(ZAL),MARGIN=2,w.resid,`*`) %*% ZAL - diag(w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    qr.d2hdv2 <- NULL
    OO1 <- matrix(0,cum_n_u_h[nrand+1L],pforpv)
    TT <- rbind(cbind(X.pv,ZAL),cbind(OO1,I))  ## aug design matrix
    if ( distinct.X.ReML ) {  
      OO1leve <- matrix(0,cum_n_u_h[nrand+1L],ncol(X.Re))
      TTleve <- rbind(cbind(X.Re,ZAL),cbind(OO1leve,I))
    }
    if (length(etaFix$beta)==ncol(X.pv) && !is.null(etaFix$v_h)) {
      maxit.mean <- 0 ## used in test near the end...
    } else if ( LMMbool ) {
      #if ( LevenbergM ) { ## only for testing LevenbergM in a case it is not needed
      #  maxit.mean <- 3 
      #} else 
      maxit.mean <- 1 ## sufficient for LMM as Hessian does not vary with beta_eta  => quadratic function
    } else { ## even h maximization in *G*LMMs 
      if ( ! is.null(phi.Fix) && ! is.null(lambda.Fix)) { ## allFix hence no true outer iteration 
        maxit.mean <- iter.mean.dispFix 
      } else maxit.mean <- iter.mean.dispVar # If phi.Fix and lambda.Fix, the only way to have 'outer' convergence is to have 'inner' convergence
    } 
  } else if (models[[1]]=="etaGLM") {
    Sig <- diag(1/w.resid)  
    TT <- X.pv
    if ( ! is.null(phi.Fix)) { ## 
      maxit.mean <- iter.mean.dispFix 
    } else maxit.mean <- iter.mean.dispVar # If phi.Fix and lambda.Fix, the only way to have 'outer' convergence is to have 'inner' convergence
  }
  iter<-0
  next_LMatrices <- NULL
  ########################################
  ######### Main loop ####################
  ########################################
  if (HL[1]=="SEM") { ## specif probit
    n_u_h <- vec_n_u_h[1] ## ugly but coherent handling of info # levels ranef
    qr.XtX <- QRwrap(crossprod(X.pv)) ##qr(t(X.pv)%*%X.pv) ## FR->FR facilement recodable Rcpp
    SEMargs$beta_eta <- beta_eta
    betalambda <- do.call(SEMbetalambda,SEMargs)
    beta_eta <- betalambda$beta_eta
    lambda_est <- betalambda$lambda_est
    u_h <- v_h <- betalambda$vs
    p_v <- betalambda$p_v
    tXinvS <- NULL
  } else while ( TRUE ) { ##the main loop with steps: new linear predictor, new leverages, new phi, new w.resid, new lambda, new fn(lambda)
    if (models[[1]]=="etaHGLM") {
      if (.SpaMM$USEEIGEN) levQ <- NULL ## allows testing if it is computed before leverage computation 
      if (LevenbergM) {
        damping <- 1/1000000 ## as suggested by Madsen-Nielsen-Tingleff... ## mauvais resultats si on part + haut
        dampingfactor <- 2
      }
      for (innerj in 1:maxit.mean) {  ## breaks when conv.threshold is reached
        ##################
        ### Inner Loop ### IF random effect *in mean predictor*: estim beta,v [only v if pforpv=0] for given phi,lambda
        ################## 
        sqrt.w1 <- sqrt(w.resid) ## if maxit.mean>1 GLM weights have been changed and the following must be recomputed
        sqrt.w2 <- sqrt(w.ranef) ##         
        sqrt.ww <- c(sqrt.w1,sqrt.w2) 
        wAugX <- sweepZ1W(TT,sqrt.ww) # sweep(TT,MARGIN=1,sqrt.ww,`*`) # rWW %*%TT
        if (useSparseM) {
          ## do nothing, no SQR computation => implemntation limited, cf comments on useSparseM
        } else {
          if ( ! .SpaMM$USEEIGEN)  {
            if ( useSparseQR ) wAugX <- Matrix(wAugX,sparse=TRUE) 
            SQR <- qr(wAugX)  ## FR->FR SLOW
          }
        }
        old_betaV <- betaV
        ######## According to 'theorem 1' of LeeL12
        ## new beta estimate from z1-a(i)
        ## where z1 is
        z1 <- eta+(y-mu)/dmudeta-off ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
        ## and a(i) (for HL(i,1)) is a(0) or a(0)+ something
        ## and a(0) depends on z2, as follows :
        if (all(lcrandfamfam=="gaussian")) {
          z2 <- rep(0,length(w.ranef)) 
          a <- rep(0,nobs)
        } else { ## HGLM: nonzero z2, nonzero a(0) ## this could perhaps make a separate block, but nevertheless sometimes qr.d2hdv2 computation...
          psi_corr <- unlist(lapply(seq(nrand), function(it) {
            u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
            if (lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") { 
              return(2*u_h[u.range]-(1+lambda_est[u.range])*u_h[u.range]^2) ## LeeL01 p.1003; to cast the analysis into the form of z2  
            } else {   
              return(psi_M[u.range])  
            } 
          }))
          z2 <- v_h + (psi_corr-u_h)*dvdu ## update since u_h,v_h updated (yes)
          #        nXn  .   nXn      nX'r'    'r'X'r'       'r'X'r'    'r'
          # a <- Sig %*% Wresid %*% ZAL %*% solve(-d2hdv2) %*% Wranef %*% z2 ## p. 963 l. 1-2; a(0) supp mat p. 6 
          aa <- w.ranef * z2
          if (is.null(qr.d2hdv2)) { 
            a <- try(solve(d2hdv2, - aa),silent=TRUE)
            if (class(a)=="try-error") {
              qr.d2hdv2 <- QRwrap(d2hdv2) 
              a <- solveWrap.vector(qr.d2hdv2,  -aa,stop.on.error=stop.on.error)
              ## FR->FR patch: 
              if (class(a)=="try-error") {
                  mess <- pastefrom("the Hessian matrix appears singular. Extreme lambda/phi value and/or extremely correlated random effects?",prefix="(!) From ")
                  message(mess)
                  cat(paste("max(lambda estimates)=",max(lambda_est)))
                  if (length(ranFix)>0) {
                    cat("; correlation parameters=")
                    cat(paste(names(ranFix),"=",ranFix))
                  }
                  largeLambdaMessages()
                  #      if (is.null(options()$error)) { ## default if not error=recover or error=stop
                  #        return(list(error="Singular augmented 'Sigma' matrix")) ## HLCor has code to handle return(list(error=...))
                  #      } else stop() ## will call options()$error i.e. (ideally) recover
                  stop()
              }
            }  
          } else { ## we already have a qr, we use it
            a <- solveWrap.vector(qr.d2hdv2, -aa,stop.on.error=stop.on.error)
          }    
          a <- Sig %*% ( w.resid * (ZAL %*% a) ) ## a(0) in LeeL12
          a <- as.numeric(a) ## incase it became a Matrix, which oddly does not fit with z1-a below...
        }         
        ## and the 'something' for a(1) is computed as follows
        if (HL[1]>0 && (! LMMbool )) { #pforpv>0 && removed since sscaled used for u_h estimation too...
          if (is.null(qr.d2hdv2)) qr.d2hdv2 <- try(QRwrap(d2hdv2),silent=TRUE) 
          if (class(qr.d2hdv2)=="try-error") {
            mess <- pastefrom("problem in 'qr.d2hdv2' computation.",prefix="(!) From ")
            stop(mess)
          }
          sscaled <- compute.sscaled() ## s detadmu
        } else sscaled <- 0L 
        ######## new estimates (tentative if LevenbergMM) 
        if (HL[1]=="SEM") {
          tXinvS <- NULL
          v_h <- solve(-d2hdv2) %*% (t(ZAL)%*% ((z1 - X.pv %*% beta_eta) * w.resid)+ z2*w.ranef)          
          betaV <- c(beta_eta,v_h)
          dbetaV <- betaV - old_betaV
        } else if (betaFirst)  { ### LeeL12 top p. 6 Appendix (code non optimise, useful for checking other algorithms) 
          ## c'est bien equivalent au calcul de Augz en fonction de sscaled essaye ci dessous selon LeeL12
          tXinvS <- calc.tXinvS(Sig,X.pv,stop.on.error,lambda_est,ranFix) ## note dependence v_h -> eta -> Sig...
          ## from a(0) to a(1) LeeL12 p. 963 col 1 l 2
          a <- as.numeric(a + Sig%*% (w.resid * sscaled)) ## in case it became a Matrix...
          rhs <-  tXinvS %*% (z1-a) ## already correct in 1.0
          qr.XtinvSX <- QRwrap(tXinvS%*%X.pv) ## looks contrived but avoids computing sqrt(Sig) (! not diag !); and XinvS%*%X.pv is small...
          beta_eta <- solveWrap.vector( qr.XtinvSX , rhs ,stop.on.error=stop.on.error) 
          v_h <- solve(-d2hdv2) %*% (t(ZAL)%*% ((z1 - X.pv %*% beta_eta) * w.resid)+ z2*w.ranef)  ## FR->FR pbb nonoptimal solve here        
          betaV <- c(beta_eta,v_h)
          dbetaV <- betaV - old_betaV
        } else { ### true augmented model, whether LevenbergM or not;
          # browser()
          tXinvS <- NULL
          if (LMMbool) {
            Augz <- c(z1,z2) ## sscaled=0L (la formule generale s'applique mais perd du temps)
          } else if(version=="1.0") { ## using control.hlfit$version... 
            Augz <- c(z1-a,z2) ## better lik for corrHLfit scotlip than the LeeL12 algo... innocent car petit de toute façon
          } else { ## solution of augmented system
            ## (1)the first operation in LevenbergMstep is to substract (eta^0,v^0): LM_wAugz <- wAugz - wAugX %*% betaV
            ##    so we keep (eta^0,v^0) here
            ## (2) what's needed here is the factor of T w on the RHS, not the factor of XinvS in the direct eqns above
            ##    As proven this gives z1-a in one algo and z1- sscaled in the other (!= version 1.0)
            Augz <- c(z1- sscaled,
                      z2+ ((1/w.ranef) * t(ZAL)) %*% (sscaled * w.resid ))  ## 
            ## z2 correction  constructed from eqs (3)(4) of LeeL12 App. p5 and consistent with eq B1 MolasL:
            ## so that row 2 of wAugX.(z1-sscaled,z2+...) = Z'W1 z1 + W2 z2 => Z W1 sscaled - W2 (...) =0 => (...)=
          }
          wAugz <- Augz*sqrt.ww  
          #if (any(is.infinite(wAugz))) browser()
          #if (any(is.nan(wAugz))) browser()
          if ( maxit.mean > 1 && LevenbergM) { ## default LevenbergM
            LM_wAugz <- wAugz - wAugX %*% betaV
            ## verif correct rhs: verif_dbetaV <- safesolve.qr.vector(SQR, LM_wAugz,stop.on.error=stop.on.error)            
            if (.SpaMM$USEEIGEN) {
              bla <- LevenbergMstepCallingCpp(wAugX=wAugX,LM_wAugz=LM_wAugz,damping=damping)
            } else bla <- LevenbergMstep(wAugX=wAugX,LM_wAugz=LM_wAugz,damping=damping,stop.on.error=stop.on.error)
            dbetaV <- bla$dbetaV
            betaV <- betaV + dbetaV
            #if (any(is.nan(betaV))) browser()
            denomGainratio <- bla$denomGainratio
          } else { ## simple aug lin or basic IRLS depending on maxit.mean
            ## QR appears faster than alternatives with crossprod(wAugX); see version 040213
            #if (useSparseM) {
            #### but slm appars to have been removed from Matrix, without notice... (30/12/2013)
            #  slmfit <- Matrix::slm(as.matrix(wAugz) ~ wAugX -1,tmpmax=100000) ## ! tmpmax... 
            #  betaV <- slmfit$coefficients
            #  ## les leverages sont diag(crossprod(forwardsolve(slmfit$chol, t(wAugX))))
            #} else {
            if (.SpaMM$USEEIGEN) {
              betaVQ <- lmwithQ(wAugX,wAugz)
              betaV <- betaVQ$coef   
              #if (any(is.nan(betaV))) browser()
              levQ <- betaVQ$Q[,seq_len(ncol(wAugX))] ## Eigen's HouseholderQ returns a square matrix...
            } else {
              betaV <- solveWrap.vector(SQR, wAugz,stop.on.error=stop.on.error) ## qr.coef(SQR, wAugz) ## vector
            }
            #}
            if (class(betaV)=="try-error") betaV <- ginv(wAugX)%*% wAugz ## occurred with large lambda either as 'init.HLfit', or by the iterative algo
            betaV <- as.numeric(betaV) #! utile ? LevenbergM produces a matrix anyway
            if (maxit.mean>1) dbetaV <- betaV - old_betaV
          } ## endif LevenbergM else...
          if (pforpv>0) {
            beta_eta <- betaV[seq_len(pforpv)]
            names(beta_eta) <- colnames(X.pv)
            beta_eta[names(etaFix$beta)] <- etaFix$beta
            if (is.null(etaFix$v_h)) v_h <- betaV[-seq_len(pforpv)] 
          } else {if (is.null(etaFix$v_h)) v_h <- betaV}
        } ## end true augmented model       
        u_h <- u_h_from_v_h(v_h)
        checkv_h <- attr(u_h,"v_h") ## verif is v_h was corrected 
        if (!is.null(checkv_h)) v_h <- checkv_h
        # print(paste(innerj," ",paste(beta_eta,collapse=" ")),quote=F)
        ####### new values of everything, only tentative if LevenbergM  
        eta <- off + X.pv %*% beta_eta + ZAL %*% v_h ## updated at each inner iteration
        ## update functions u_h,v_h
        keep_wranefblob <- wranefblob #######################
        wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h)
        w.ranef <- wranefblob$w.ranef
        dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
        dvdu <- wranefblob$dvdu
        keep_blob <- blob ############################
        blob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen,priorWeights=prior.weights) 
        mu <- blob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
        dmudeta<-blob$dmudeta
        Vmu <- blob$Vmu ## if Bin/Pois, O(n)
        ## update functions of v_h -> blob
        w.resid<-as.vector(blob$GLMweights /phi_est) ## 'weinu', must be O(n) in all cases
        #### update fns of v_h -> blob -> w.resid
        if ( LevenbergM ) { ## for LMM, d2hdv2 is constant // mu, hence it is constant over LevenbergMarquardt iterations
          if (pforpv>0) keep_Sig <- Sig
          keep_d2hdv2 <- d2hdv2
        }
        if (pforpv>0) {
          Sig <- ZWZt(ZAL,1/w.ranef) # + diag(1/w.resid) ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta estimate
          diag(Sig) <- diag(Sig) + 1/w.resid ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta estimate
        }
        d2hdv2 <- calcDhDv2(ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
        qr.d2hdv2 <- NULL
        ######### for LevenbergM, we check there was a progress, and restore everything otherwise
        if (LevenbergM) {      
          if (HL[1]==0) {
            currentlik <- calc.p_v(only.h=TRUE)$hlik ##  use h in PQL/L (> v1.0)  
          } else currentlik <- calc.p_v()$p_v
          if (innerj==1) { ## not a good idea to use the p_v computed for all u_h initially set to zero
            gainratio <- 1
          } else {
            if (currentlik==-Inf) { ## obs in binary probit with extreme eta... 
              gainratio <- -1
            } else {
              gainratio <- 2*(currentlik-oldlik)/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
            }
          }
          if (class(try(if (gainratio<0) {}))=="try-error") { ## nnote that NaN *is* numeric!
            mess <- pastefrom("'try(if (gainratio<0)...) failed.",prefix="(!) From ") 
            stop(mess)
          }
          
          if (gainratio<0) { ## restore everything
            currentlik <- oldlik ## important to restore this for the next test !
            betaV <- old_betaV
            if (pforpv>0) {
              beta_eta <- betaV[seq_len(pforpv)]
              names(beta_eta) <- colnames(X.pv)
              v_h <- betaV[-seq_len(pforpv)] 
            } else v_h <- betaV
            u_h <- u_h_from_v_h(v_h) 
            checkv_h <- attr(u_h,"v_h") ## verif is v_h was corrected 
            if (!is.null(checkv_h)) v_h <- checkv_h
            eta <- off + X.pv %*% beta_eta + ZAL %*% v_h ## updated at each inner iteration
            wranefblob <- keep_wranefblob ########################## restoration
            w.ranef <- wranefblob$w.ranef ; dlogWran_dv_h <- wranefblob$dlogWran_dv_h ; dvdu <- wranefblob$dvdu
            blob <- keep_blob ######################################### restoration
            mu <- blob$mu ; dmudeta<-blob$dmudeta ; Vmu <- blob$Vmu ; w.resid<-as.vector(blob$GLMweights /phi_est) 
            if (pforpv>0) Sig <- keep_Sig ###### more restoration
            d2hdv2 <- keep_d2hdv2 ###### more restoration
            #            qr.d2hdv2 <- keep_qr.d2hdv2 ###### more restoration
            damping <- damping*dampingfactor
            dampingfactor <- dampingfactor*2 
          } else { ## cf Madsen-Nielsen-Tingleff again, and as in levmar library by Lourakis
            oldlik <- currentlik
            damping <- damping * max(1/3,1-(2*gainratio-1)^3)  
            dampingfactor <- 2
          }
          #print(paste("damping: ",damping))
        }
        #if(length(beta_eta)>0) browser()
        ########## nearly done with one inner iteration
        if (verbose["trace"]) {
          cat(paste("Inner iteration ",innerj,sep=""))
          if (LevenbergM) cat(paste(", likelihood= ",currentlik,sep=""))
          cat("\n")
          if (innerj>1) cat("norm.dbetaV=",sqrt(sum(dbetaV^2)))
          cat(paste("; beta_eta=",paste(beta_eta,collapse=", ")))
          cat("\n")
          print("================================================")
        } 
        if (maxit.mean>1) {
          ## the convergence on v_h^2 must be relative to lambda; this raises questions about the lowest meaningful lambda values.
          relvariation <- dbetaV*(c(rep(1,pforpv),w.ranef)) ## 21/01/2013 
          if (mean(abs(relvariation)) < conv.threshold) break; ## FR->FR mean(abs) is not standard ?  
        }
      } ## end for (innerj in 1:maxit.mean)
      ######################
      ### end Inner Loop ### HL(.,.) estim of beta, v for given lambda,phi
      ######################
    } else if (models[[1]]=="etaGLM") {
      if (pforpv>0) {
        for (innerj in 1:maxit.mean) {  ## breaks when conv.threshold is reached
          old_beta_eta <- beta_eta
          z1 <- eta+(y-mu)/dmudeta-off ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
          tXinvS <- calc.tXinvS(Sig,X.pv,stop.on.error,lambda_est,ranFix)
          rhs <-  tXinvS %*% z1
          qr.XtinvSX <- QRwrap(tXinvS%*%X.pv)
          beta_eta <- solveWrap.vector( qr.XtinvSX , rhs ,stop.on.error=stop.on.error) 
          names(beta_eta) <- colnames(X.pv)
          beta_eta[names(etaFix$beta)] <- etaFix$beta ## added 03/2014
          dbetaV <- beta_eta - old_beta_eta
          eta <- off + X.pv %*% beta_eta ## updated at each inner iteration
          blob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen,priorWeights=prior.weights) 
          mu <- blob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
          dmudeta<-blob$dmudeta
          Vmu <- blob$Vmu ## if Bin/Pois, O(n)
          ## update functions of v_h -> blob
          w.resid <- as.vector(blob$GLMweights /phi_est) ## 'weinu', must be O(n) in all cases
          Sig <- diag(1/w.resid) ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta estimate
          ########## nearly done with one inner iteration
          if (verbose["trace"]) {
            print(paste("Inner iteration ",innerj,sep=""))
            print_err <- c(beta_eta=beta_eta)
            if (innerj>1) print_err <- c(norm.dbetaV=sqrt(sum(dbetaV^2)),print_err)
            print(print_err)
            print("================================================")
          } 
          if (maxit.mean>1) {
            if (mean(abs(dbetaV)) < conv.threshold) break; ## FR->FR mean(abs) is not standard ?  
          }
        } ## end for (innerj in 1:maxit.mean)
      }
    }
    ##########
    if (models[[1]]=="etaHGLM") {
      if (is.null(lambda.Fix) || is.null(phi.Fix)) {
        if (maxit.mean==0) {
          stop("(!) Computation of leverages with maxit.mean=0: check that this is meaningful.")
        } # ELSE rWW was updated in the inner loop for betaV
        if ( distinct.X.ReML) { 
          #          wAugXleve <- sweep(TTleve,MARGIN=1,sqrt.ww,`*`) # rWW %*% TTleve ## note that TTleve should have been updated if ZAL has been so.
          #          SQRleve <- qr(wAugXleve)
          #          hatval <- rowSums(qr.qy(SQRleve, diag(1, nrow = nrow(wAugXleve), ncol = ncol(wAugXleve)))^2) ## even for ML, nontrivial "partition" between phi and lambda
          wAugXleve <- sweepZ1W(TTleve,sqrt.ww)
          hatval <- leverages(wAugXleve)
        } else { ## basic REML, leverages from the same matrix used for estimation of betaV (even simply V)
          #        if (pforpv==0) { ## SQR not previously computed
          #          wAugX <- rWW%*%TT 
          #          SQR <- qr(wAugX)
          #        } ## else this is already ~ up to date from the inner loop 
          if (.SpaMM$USEEIGEN) {
            #            if (is.null(levQ)) levQ <- Rcpp_qr_Q(wAugX)[,seq_len(ncol(wAugX))]
            #            hatval <- rowSums(levQ^2)
            if (is.null(levQ)) {
              if (FALSE)  {
                ## wAugX updated not only by change in lambda, phi but also GLM weights -> leverage comput difficult to optimize  
                ## the following could be useful if the GLM wights are unity, phiScal et lamScal...
                # mais il manque bouts pour pour produire <u> et <d> | unperturbed RpR = u d u'                
                #                hatval <- LevPerturbedQCpp(perturbedwAugX=wAugX,pforREML=ncol(X.Re),
                #                                           RpRu = <u>,RpRd=<d>,lamOverLam0=lambda/lambda0,phiOverPhi0=phi/phi0)
              } else hatval <- leverages(wAugX)
            } else hatval <- rowSums(levQ^2)   
          } else hatval <- rowSums(qr.qy(SQR, diag(1, nrow = nrow(wAugX), ncol = ncol(wAugX)))^2) ## == but faster than rowSums(qr.Q(SQR)^2) !
          # slmfit from old versions of SparseM has been tried:         hatval <- diag(crossprod(SparseM::forwardsolve(slmfit$chol, t(wAugX)))) 
        }
        if (any(abs(hatval) > 1 - 1e-8)) {
          hatval <- ifelse(abs(hatval) > 1 - 1e-8, 1 - 1e-8,hatval)
          warningList$leveLam1 <-TRUE
        }
        lev_phi <- hatval[1:nobs] ## for the error residuals (phi)
        lev_lambda <- hatval[(nobs+1L):(nobs+cum_n_u_h[nrand+1L])]  ## for the ranef residuals (lambda)
      }
    } else { ## GLM
      # rWW <- diag(sqrt(w.resid))  
      if ( distinct.X.ReML ) { 
        wAugXleve <- sweep(X.Re,MARGIN=1,sqrt(w.resid),`*`) # rWW%*%X.Re
    #    SQRleve <- qr(wAugXleve)
    #    lev_phi <- rowSums(qr.qy(SQRleve, diag(1, nrow = nrow(wAugXleve), ncol = ncol(wAugXleve)))^2) ## even for ML, nontrivial "partition" between phi and lambda
        lev_phi <- leverages(wAugXleve)
      } else { ## basic REML, leverages from the same matrix used for estimation of beta
        wAugX <- sweep(X.pv,MARGIN=1,sqrt(w.resid),`*`) # rWW %*% X.pv 
    #    SQR <- qr(wAugX)
    #    lev_phi <- rowSums(qr.qy(SQR, diag(1, nrow = nrow(wAugX), ncol = ncol(wAugX)))^2)
        lev_phi <- leverages(wAugX)
      }
      wrongphi <- which(lev_phi > 1 - 1e-8) 
      if (length(wrongphi)>0) {
        lev_phi[wrongphi] <- 1 - 1e-8
        #        warningList$leveLam1 <-T ## Lam ?
      }
    }
    d2mudeta2 <- d2mudeta2fn(link=family$link,mu=mu,eta=eta,BinomialDen=BinomialDen)
    if (verbose["trace"]) {print(paste("beta=",paste(signif(beta_eta,4),collapse=", ")),quote=F)}
    if (HL[2]>0) { ## LeeN01.... ie the + in 'EQL+'
      #### Then we CORRECT the leverages previously computed from  the hat matrix
      ## (0): hat matrix -> p, notEQL -> tilde(p), (1): full correction -> q 
      ## first the d log hessian / d log lambda or phi corrections then the notEQL correction
      ## For the d log hessian first the derivatives of GLM weights wrt eta 
      ##################### noter que c'est le coef2 de HL(1,.), but mu,eta may have been updated since coef2 was computed
      if (canonicalLink) {
        dlW_deta <- d2mudeta2 / dmudeta
      } else if (family$family=="binomial" && family$link=="probit") { ## ad hoc non canonical case 
        muFREQS <- mu/BinomialDen
        dlW_deta <- -2*eta - dnorm(eta)*(1-2*muFREQS)/(muFREQS*(1-muFREQS))
      } else if (family$family=="Gamma" && family$link=="log") { ## ad hoc non canonical case 
        dlW_deta <- rep(0L,length(eta)) ## because the GLM weight is 1 ## correct rep v1.2
      } else {
        ## we need to update more functions of mu...
        tmblob <- thetaMuDerivs(mu,BinomialDen,family$family)
        Dtheta.Dmu <- tmblob$Dtheta.Dmu
        D2theta.Dmu2 <- tmblob$D2theta.Dmu2
        ## ... to compute this:
        D2theta.Deta2_Dtheta.Deta <- (D2theta.Dmu2 * dmudeta^2 + Dtheta.Dmu * d2mudeta2)/(Dtheta.Dmu * dmudeta)
        dlW_deta <- d2mudeta2 / dmudeta + D2theta.Deta2_Dtheta.Deta
      }
      ## we join this with the deriv of log w.ranef wrt v_h
      if (models[[1]]=="etaHGLM") {
        dlW_deta_or_v <- c(dlW_deta, dlogWran_dv_h)  ## vector with n+'r' elements
        # dlogWran_dv_h is 0 gaussian ranef; d2mudeta2 is 0 for identity link => vector is 0 for LMM
        ## else we continue the computation of the d log hessian term
        if (any(dlW_deta_or_v!=0L)) {
          ## 
          if(models[[1]]=="etaHGLM" && is.null(lambda.Fix)) {
            ############### all random effect models are canonical conjugate except the inverse.Gamma(log) ############### 
            dlogfthdth <- (psi_M - u_h)/lambda_est ## the d log density of th(u)
            neg.d2f_dv_dloglam <- unlist(lapply(seq(nrand), function(it) {
              u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
              if (lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") { ## g(u) differs from theta(u) : cf oklink dans preprocess pour detection des cas
                ## same computation as alternative case, except that first we consider dlogfthdv=dlogfthdth * dth/dv=1/u for th(u)=-1/u, v=log(u)
                return(dlogfthdth[u.range] / u_h[u.range])  
              } else { ## v=g(u) = th(u)  
                return(dlogfthdth[u.range]) ## (neg => -) (-)(psi_M-u)/lambda^2    *    lambda.... 
              } 
            }))
            if (is.null(qr.d2hdv2)) {
              if (.SpaMM$USEEIGEN) {
                dvdloglamMat <- pseudoSolvediag(d2hdv2,as.vector(neg.d2f_dv_dloglam)) ## FR->FR dangereux car contient solve(R)
              } else {
                qr.d2hdv2 <- QRwrap(d2hdv2)
                ## mff solve(A,diag(b)) est pareil que solve(A,diag(1)) * b ('*')  
                dvdloglamMat <- solveWrap.matrix(qr.d2hdv2, diag( as.vector(neg.d2f_dv_dloglam) ), stop.on.error=stop.on.error) # rXr       
              }
            } else {
              dvdloglamMat <- solveWrap.matrix(qr.d2hdv2, diag( as.vector(neg.d2f_dv_dloglam) ), stop.on.error=stop.on.error) # rXr                     
            }
            if (is.null(qr.d2hdv2))  
            if (class(dvdloglamMat)=="try-error") {
              mess <- pastefrom("problem in dvdloglamMat computation.",prefix="(!) From ")
              warning(mess)
              dvdloglamMat <- sweep(ginv(d2hdv2),MARGIN=2,as.vector(neg.d2f_dv_dloglam),`*`) ## ginv(d2hdv2) %*% diag( as.vector(neg.d2f_dv_dloglam))
            }
            # next line uses only vector X matrix :
            dleve <- ((hatval * dlW_deta_or_v) %*% ZALI ) %*% dvdloglamMat # (r+n) . (r+n)Xr . rXr = r (each element is a sum over r+n terms= a trace)
            lev_lambda <- lev_lambda - as.vector(dleve)  
          }
          ## 
          if(is.null(phi.Fix)) {
            dh0deta<-( w.resid *(y-mu)/dmudeta ) ## 12/2013 supp BinomialDen (soit Bin -> phi fixe=1, soit BinomialDen=1)
            ## cf calcul dhdv, but here we want to keep each d/d phi_i distinct hence not sum over observations i 
            ## code corrected here 12/2013; this is dh0dv = neg.d2h0_dv_dlogphi (eta always linear in v :-) and w.resid always propto 1/phi)
            neg.d2h0_dv_dlogphi <- sweep(t(ZAL),MARGIN=2,as.vector(dh0deta),`*`) ## dh0dv <- t(ZAL) %*% diag(as.vector(dh0deta)) ## nXr each ith column is a vector of derivatives wrt v_k
            if (is.null(qr.d2hdv2)) qr.d2hdv2 <- QRwrap(d2hdv2) 
            dvdlogphiMat <- solveWrap.matrix(qr.d2hdv2, neg.d2h0_dv_dlogphi , stop.on.error=stop.on.error)  # rXn       
            if (class(dvdlogphiMat)=="try-error") {
              mess <- pastefrom("problem in dvdlogphiMat computation.",prefix="(!) From ")
              stop(mess) ## warning + ginv for phi... !
            }
            dleve <- ((hatval * dlW_deta_or_v) %*% ZALI) %*% dvdlogphiMat # (r+n) . (r+n)Xr . rXn = n (each element is a sum over r+n terms= a trace)
            lev_phi <- lev_phi - as.vector(dleve)  
          } 
        }
      } 
    }
    if (HL[2]>1) {
      stop("Need a_i correction in Table 7 of NohL07 ie derivatives of second order correction wrt dips param.")
    }
    ## see instances of digamma in Genstat code
    if (models[[1]]=="etaHGLM") {
      if (is.null(lambda.Fix)) {
        if (HL[3]!=0 ) { ## ie , p_bv(h), not EQL p_bv(q+), LeeNP p89; distinction does not arise for PQL <=> Gaussian ranefs...  
          ## d h/ d !log! lambda correction (nul for gaussian ranef)
          ## ! correction for not using the deviance residuals as approx for the distribution of the random effects. It's not specifically ReML !
          ## this is a trick for still using deviances residuals in the Gamma GLM
          notEQL <- unlist(lapply(seq(nrand), function(it) {
            u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
            loclambda <- lambda_est[u.range]
            blob <- switch(lcrandfamfam[it], 
                           gaussian=rep(0,length(u.range)),
                           gamma=1+2*(log(loclambda)+digamma(1/loclambda))/loclambda,## cf notes on p. 89 of the book
                           "inverse.gamma"=1+2*(log(loclambda)-loclambda+digamma(1+(1/loclambda)) )/loclambda, ## appears to be the same as for the gamma case [digamma(1+x)=digamma(x)+1/x]... 
                           beta=1-2*(digamma(1/loclambda)/loclambda)+2*(digamma(1/(2*loclambda))/loclambda)+log(4)/loclambda
            ) ## as in HGLMMM PoissonEst.r
            blob    
          }))
          lev_lambda <- lev_lambda + notEQL   
        }
        wronglambda <- which(lev_lambda > 1 - 1e-8) 
        if (length(wronglambda)>0) {
          lev_lambda[wronglambda] <- 1 - 1e-8
          #        warningList$leveLam1 <-T ## Lam ?
        }
      }
    }
    if (HL[3]!=0 && family$family=="Gamma" && is.null(phi.Fix) ) { ## Gaussian: no lev correction; poisson and binom: no phi    
      ## d h/ d !log! phi correction (zero for gaussian residual error, no phi for other distribs)
      ## exact correction for not using the deviance residuals as approx for the likelihood of the residual error. It's not specifically ReML ! 
      ## this is a trick for still using deviances residuals in the Gamma GLM
      notEQL <- 1+2*(log(phi_est)+digamma(1/phi_est))/phi_est ## LNP p. 89 and as in HGLMMM IWLS_Gamma; en jouant sur HL[3] on voit que ca ameliore les estim
      lev_phi <- lev_phi + notEQL   
    }    
    ## compute deviance_residuals 
    ## updated residuals from updated mu must be used (LeeNP p.161) [not so in dhglmfit !!]
    deviance_residual <- family$dev.resids(y,mu,wt=1) 
    ######### Dispersion Estimates for phi #####################
    RespLink_disp<-"log"  ## currently not under user control
    if (is.null(phi.Fix)) { ## if phi is estimated (phi.Fix set to 1 for Poisson, Binomial)
      ## leverages have been computed before the  inner loop, which did not change the design matrices 
      Offset_disp <- attr(resid.predictor,"offset") 
      if (models[[3]]=="phiScal") {
        next_phi_est <- sum(deviance_residual*prior.weights)/sum((1-lev_phi)) ## NOT in linkscale
        if (next_phi_est<1e-8) next_phi_est <- 1e-8 # e-10 not high enough (d2hdv2 -> lad -> p_v unstable)
        if (verbose["trace"]) {print(paste("phi_est=",signif(next_phi_est,4)),quote=F)}
      } else if (models[[3]]=="phiGLM") { ## there is a phi predictor to estimate but no ranef in this predictor
        resglm_phi <- dispGammaGLM(dev.res=deviance_residual*prior.weights,lev=lev_phi,X=X_disp,offset=Offset_disp)
        # glm(I(deviance_residual/(1-lev_phi))~X_disp-1,weights=(1-lev_phi)/2,family=Gamma(log))
        if (verbose["trace"]) {print(paste("phi coefficients=",paste(signif(resglm_phi$coefficients,4),collapse=", ")),quote=F)}
        next_phi_est <- resglm_phi$fitted.values
        lowphi <- which(next_phi_est < 1e-08)
        next_phi_est[lowphi] <- 1e-08 ## to avoid problems with nearly singular matrices
      } else if (models[[3]]=="phiHGLM") { ## random effect(s) in predictor for phi
        stop("random effects in predictor or residual variance (phi) not yet implemented")
        ## there is a template code with comments in version 260812 of HLfit
        reshglm_phi <- list()
      } 
      if (all(abs(next_phi_est-phi_est) < conv.threshold* (phi_est+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.phi <- TRUE ## 'weak convergence'... 
      } else conv.phi <- FALSE
    } else {conv.phi <- TRUE} ## there is a phi.Fix, already -> phi_est
    ######### Dispersion Estimates for lambda #####################
    if (models[[1]]=="etaHGLM" && is.null(lambda.Fix)) { ## lambda must be estimated     
      if (any(abs(lev_lambda) > 1 - 1e-8)) {
        lev_lambda<- ifelse(abs(lev_lambda) > 1 - 1e-8, 1 - 1e-8,lev_lambda)
        warningList$leveLam1 <- TRUE
      }
      ## Build pseudo response for lambda GLM/HGLM
      resp_lambda <- matrix(0,cum_n_u_h[nrand+1L],1L)
      #########################
      for (it in 1:nrand) {
        u.range <- (cum_n_u_h[it]+1L):cum_n_u_h[it+1L]
        resp_lambda[u.range] <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=1) ## must give d1 in table p 989 de LeeN01...
      }
      ## then analyse pseudoresponse
      if (all(models[[2]]=="lamScal")) { ## simplified estimate
        if (any(attr(X_lamres,"Xi_cols")>1)) {
          ## handling correlation in random slope models ### valid slmt pr gaussian ranefs, verif dans preprocess
          LMatricesBlob <- makeCovEst(u_h,ZAlist=ZAlist,cum_n_u_h=cum_n_u_h
                                      ,X_lamres=X_lamres,prev_LMatrices=next_LMatrices,
                                      userLfixeds=attr(ZALlist,"userLfixeds"),
                                      lev_lambda=lev_lambda)
          next_LMatrices <- LMatricesBlob$next_LMatrices ## a list of matrices with NULL elements for non-random-slope terms
          next_lambda_est <- LMatricesBlob$next_lambda_est ## a full-length vector with values only in the appropriate u ranges 
          ## FR->FR code a nettoyer en appliquant CovEst a chaque u.range...
        } else next_lambda_est <- numeric(length(u_h)) ## next_LMatrices remains an empty list()
        ## then fill all missing values 
        for (it in seq_len(nrand)) {
          if (attr(X_lamres,"Xi_cols")[it]==1) {
            u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
            unique.lambda <- sum(resp_lambda[u.range])/sum(1-lev_lambda[u.range]) ## NOT in linkscale 
            unique.lambda <- pmax(unique.lambda,1e-8) # as for phi
            unique.lambda <- pmin(unique.lambda,.SpaMM$maxLambda)  
            if (verbose["trace"]) {print(paste("lambda=",signif(unique.lambda,4)),quote=F)}
            next_lambda_est[u.range] <- rep(unique.lambda,length(u.range))
          }
        }
        ##########################################################
      } else if (any(models[[2]]=="lamHGLM")) { ## if ranef in predictor lambda...
        stop("random effects in predictor or ranef variance (lambda) not yet implemented")
        ## there is a template code with comments in version 260812 of HLfit
        reshglm_lambda <- list()
      } else { ## any mixture of lamScal and lamGLM and we can use a single X_lamres for all of them
        stop("Linear predictor for ranef variance (lambda) not yet implemented")
        ## FR->FR code non teste
        resglm_lambda <- dispGammaGLM(dev.res=resp_lambda,lev=lev_lambda,X=X_lamres)
        next_lambda_est <- resglm_lambda$fitted.values ## $fitted.value is NOT in linkscale, contrairement a $coefficients
        lowlambda <- which(next_lambda_est < 1e-08)
        next_lambda_est[lowlambda] <- 1e-08 ## to avoid problems with nearly singular matrices
      }
      if ( 
          ## for low values, precision on lambda must be O(v_h^2) ... need precision in relative terms:
          all(abs(log(pmax(next_lambda_est,1e-06)/pmax(lambda_est,1e-06))) < conv.threshold)  
        && 
          all(abs(next_lambda_est-lambda_est) < conv.threshold* (lambda_est+0.1)) ## ie 1e-6 ~ 1e-5*(1e-6+0.1) 
      ) { 
        conv.lambda <- TRUE ## 'weak convergence'... 
      } else conv.lambda <- FALSE
    } else { 
      # lambda_est remains = lambda.Fix
      conv.lambda <- TRUE
    } ## end if null lambda.Fix else ...
    if (! is.null(corr_est)) {
      corrEst.args$formula <- Predictor(formula=corrEst.form,offset=off + X.pv %*% beta_eta) ## FR->FR ugly input for offset...
      corrEst.args$init.corrHLfit <- corr_est ## this entails use of optim() (or another Optimizer) on these parameters
      if (nrand>1) stop("code needed for corr Estim within HLfit with multiple lambda parameters") ## FR->FR
      corrEst.args$ranFix$lambda <- unique(lambda_est)
      corrEst.args$ranFix$phi <- phi_est 
      corrEst.args$init.HLfit$v_h <- v_h ## substantial gain of time (no need for inner call to provide.resglm which takes time) 
      ## corrEst.args$HLmethod <- .... ## default REML  ~ ML here
      #       if (FALSE) { ## seems to work...
      #       locprocessed <- preprocess(control.HLfit=control.HLfit,HLmethod=HLmethod,
      #                                  predictor=Predictor(formula=corrEst.form,offset=off + X.pv %*% beta_eta),phi.Fix=phi_est,                 
      #                                  resid.predictor=resid.formula, ## must be ignored, but no default... =>preprocess could be improved
      #                                  REMLformula=corrEst.args$REMLformula,data=data,
      #                                  family=family,BinomialDen=BinomialDen,rand.family=rand.family)
      #       corrEst.args$processed <- locprocessed ## risky
      #       }
      pff <- do.call("corrHLfit",corrEst.args)
      next_corr_est <- pff$corrPars[corrNames] ## rho,nu,  pas trRho, trNu
      #FR->FR maybe conv_threshold a bit strict here...
      if (all(abs(log(unlist(next_corr_est)/unlist(corr_est))) < conv.threshold) ) { ## 
        conv.corr <- TRUE ## this is the simplest, best case. ## but if slow geometric decrease to 0, this is never true 
      } else 
      if (all(abs(unlist(next_corr_est)-unlist(corr_est)) < conv.threshold* (unlist(corr_est)+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.corr <- TRUE ## 'weak convergence'... 
      } else conv.corr <- FALSE
    } else {
      conv.corr <- TRUE
    }
    iter<-iter+1 ## here first from 0 to 1
    ## We need to make sure either that convergence of lambda occurred on a relative log scale ( loop not stopping at max.iter !) so that the v_h are very accurate on same scale
    ## or that the v_h's are computed with the very latest lambda, otherwise a call with ranFix$lambda does not yield the same result as estimated lambda
    if ( conv.phi && conv.lambda && conv.corr) {
      ## do not update phi and lambda so that the v_h where computed from the latest lambda_est in particular
      break 
    } else if (iter>=max.iter) { ## did not converge...
      break 
    } else { ## update and continue
      #if(length(beta_eta)>0) browser()
      if ( is.null(phi.Fix)) {
        phi_est <- next_phi_est
        w.resid <- as.vector(blob$GLMweights /phi_est) ## 'weinu', must be O(n) in all cases; blob was updated when eta was
      }
      if (length(next_LMatrices)>0) {
        ZALlist <- compute.ZALlist(LMatrix=next_LMatrices,ZAlist=ZAlist,Groupings=Groupings)
        ## meme block que pour corr.est
        ZAL <- post.process.ZALlist(ZALlist,predictor=predictor)
        TT <- rbind(cbind(X.pv,ZAL),cbind(OO1,I))  ## aug design matrix
        ZALI <- rbind(ZAL,I) 
        if (distinct.X.ReML) {
          TTleve<-rbind(cbind(X.Re,ZAL),cbind(OO1leve,I))
        }
      }
      if ( ! is.null(corr_est)) {
        corr_est <- next_corr_est 
        LMatrix <- attr(pff$predictor,"LMatrix")
        ZALlist <- compute.ZALlist(LMatrix=LMatrix,ZAlist=ZAlist,Groupings=Groupings)
        ZAL <- post.process.ZALlist(ZALlist,predictor=predictor)
        TT <- rbind(cbind(X.pv,ZAL),cbind(OO1,I))  ## aug design matrix
        ZALI <- rbind(ZAL,I) 
        if (distinct.X.ReML) {
          TTleve<-rbind(cbind(X.Re,ZAL),cbind(OO1leve,I))
        }
      }
      if (models[[1]]=="etaHGLM" && is.null(lambda.Fix)) { ## lambda was modified
        lambda_est <- next_lambda_est
      }
      if (models[[1]]=="etaHGLM" && (is.null(lambda.Fix) || ! is.null(corr_est))) { ## lambda or u_h were modified
        wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h)
        w.ranef <- wranefblob$w.ranef
        dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
        dvdu <- wranefblob$dvdu
      } 
      if (models[[1]]=="etaHGLM") {
        if (is.null(lambda.Fix) || is.null(phi.Fix) || ! is.null(corr_est)) { ## w.ranef or w.resid or ZAL were modified 
          Sig <- ZWZt(ZAL,1/w.ranef) # + diag(1/w.resid) ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid)
          diag(Sig) <- diag(Sig) + 1/w.resid ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta estimate
          d2hdv2 <- calcDhDv2(ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
          qr.d2hdv2 <- NULL
        }
      } else { ## no random effect
        if ( is.null(phi.Fix)) Sig <- diag(1/w.resid) 
      }
      if (verbose["trace"]) {
        print(paste("iteration ",iter,sep=""))
        ## inappropriately large output
        #if ( is.null(phi.Fix)) {print.arg <- c(`next_phi_est`=next_phi_est)} else {print.arg <- c(`phi.Fix`=phi.Fix)} 
        #if ( is.null(lambda.Fix)) {print.arg <- c(print.arg,`next_lambda_est`=next_lambda_est)} else {print.arg <- c(print.arg,`lambda.Fix`=lambda.Fix)} 
        #print(print.arg)
        print("================================================")
      } 
    } 
  } ## end main loop while ( TRUE )
  ########################################
  ######### END main loop ################
  ########################################
  if (verbose["trace"]) {
    if (iter==max.iter) {
      mess <- paste("(beta,v)/lambda/phi iterations failed to converge in",max.iter,"iterations")
      mess <- pastefrom(mess,prefix="(!) From ")
      message(mess)
    } else {
      message(paste("(beta,v)/lambda/phi iterations in HLfit() converged in",iter,"iterations"))
    }
  }
  #if (family$family %in% c("gaussian","Gamma")) {
  #  mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-lev_phi))
  #}
  if (pforpv>0 && max.iter >0) { ## condition on max.iter <=> some params have been fitted
    if (is.null(tXinvS)) tXinvS <- calc.tXinvS(Sig,X.pv,stop.on.error,lambda_est,ranFix) ## slow when phi ->0 ...     
    if (class(tXinvS)=="try-error") {
      beta_se <- rep(Inf,pforpv) ## maybe...
    } else {
      beta_cov <- try(solve(tXinvS%*%X.pv),silent=TRUE) ## solve(small matrix !)
      if (class(beta_cov)=="try-error") {
        beta_se <- rep(Inf,pforpv) ## maybe...
      } else {
        beta_se <- diag(beta_cov)
        if (any(beta_se<0)) { ## divergence of tXinvS%*%X.pv leads to negative variance estimates
          beta_se <- rep(Inf,pforpv) ## maybe... 
        } else beta_se <- sqrt(beta_se)
      }
    }
    #    if (any(is.infinite(beta_se))) {
    #      message("Suspected divergence of lambda estimates. Check model formula (wrong offset for example),")
    #      message(" otherwise try increasing control.HLfit$iter.mean.dispVar")
    #    }
  } else beta_se <- NULL 
  ######################
  ######################
  ######################
  ##### LAMBDA
  if ((HL[1]!="SEM") && models[[1]]=="etaHGLM" && is.null(lambda.Fix)) {
    if (all(models[[2]]=="lamScal")) { ## there is a single X_lamres for all lambda's
      ## to compute the se we need the GLM residuals etc. So if the GLM has not been previously used it's better to use it here
      resglm_lambda <- dispGammaGLM(dev.res=resp_lambda,lev=lev_lambda,X=X_lamres)
      summ_glm_lambda <- summary(resglm_lambda,dispersion=2) ## 2 dispersion=for the gamma GLM, probably 
      p_lambda <- ncol(X_lamres) 
      linkscale.lambda <- summ_glm_lambda$coefficients[1:p_lambda]
      lambda_se <- summ_glm_lambda$coefficients[(p_lambda+1L):(2L*p_lambda)] ## ici probleme de transfo log      
      if (any(attr(X_lamres,"Xi_cols")>1)) { ## 
        ## here dispGammaGLM has reestimated the (log of) the lambda_est given as $d by random slope procedure
        ## so we can use this; but we need to display the corrmat...
      } 
    } else {
      stop("From HLfit: 'lamHGLM' and 'lamGLM' not fully implemented.")
      ## there is a template code with comments in version 260812 of HLfit
    }
  }       
  ##### PHI
  if ( is.null(phi.Fix)) {
    if (models[[3]]=="phiHGLM") {
      ## there is a template code with comments in version 260812 of HLfit
    } else {
      #          if (models[[3]]=="phiScal") {
      ## We need the $coefficients. So if the GLM has not been previously called it's better to call it here
      resglm_phi <- dispGammaGLM(dev.res=deviance_residual*prior.weights,lev=lev_phi,X=X_disp,offset=Offset_disp)
      res2<-summary(resglm_phi)
      #          }
      p_disp<-ncol(X_disp) 
      linkscale.phi <- res2$coefficients[1:p_disp]
      phi_se <- res2$coefficients[(p_disp+1L):(2L*p_disp)]
    }
  } 
  ########## LIKELIHOODS
  #    theta<-theta.mu.canonical(mu/BinomialDen,family$family)  
  if (HL[1]=="SEM") {
    APHLs <- list(p_v=p_v,p_bv=p_v) ## FR->FR attributes ?
  } else {
    if (models[[1]]=="etaHGLM" && pforpv==0) { 
      d2hdv2 <- calcDhDv2(ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    }
    calcpv <- calc.p_v(returnLad=TRUE)
    if (models[[1]] != "etaHGLM" && models[3] != "phiHGLM") { ## ie GLM, not HGLM
      ml <- calcpv$clik ## vanilla likelihood
      d2hdx2 <- - ZtWZ(X.Re,w.resid)  ## t(X.Re)%*%Wresid%*%X.Re ## X should be the one for leverages
      if (.SpaMM$USEEIGEN) {
        lad <- LogAbsDetCpp(d2hdx2/(2*pi))
      } else lad <- determinant(d2hdx2/(2*pi))$modulus[1] ## REML for estim phi
      rl <- ml - lad/2
      cAIC<- -2*ml+2*pforpv
      d2hdbv2 <- - d2hdx2 ## FR->FR util de deux notations ?
      hlik <-ml 
      ladbv <- 0
    } else { ## add likelihood of ranef
      if (models[[1]]=="etaHGLM") {
        clik <- calcpv$clik
        hlik <- calcpv$hlik
        p_v <- calcpv$p_v 
        ## see readable account of aic in HaLM07
        if (ncol(X.Re)>0) {
          d2hdbv2 <-rbind(cbind(ZtWZ(X.Re,w.resid), crossprod(X.Re,sweep(ZAL,MARGIN=1,w.resid,`*`))),
                          cbind(crossprod(ZAL,sweep(X.Re,MARGIN=1,w.resid,`*`)), - d2hdv2))
          #         d2hdbv2 <-rbind(cbind((t(X.Re)%*%Wresid%*%X.Re),(t(X.Re)%*%Wresid%*%ZAL)),
          #                         cbind((t(ZAL)%*%Wresid%*%X.Re),(-1*d2hdv2)))
          if (.SpaMM$USEEIGEN) {
            ladbv <- LogAbsDetCpp(d2hdbv2/(2*pi))
          } else ladbv <- determinant(d2hdbv2/(2*pi))$modulus[1]
        } else {
          d2hdbv2 <- d2hdv2
          ladbv <- calcpv$lad
        }
        if (AIC) { ## diff de d2hdbv2 slmt dans dernier bloc
          # d2clikdbv2 <-rbind(cbind((t(X.Re)%*%Wresid%*%X.Re),(t(X.Re)%*%Wresid%*%ZAL)),
          #                  cbind((t(ZAL)%*%Wresid%*%X.Re),(t(ZAL)%*%Wresid%*%ZAL)))
          d2clikdbv2 <-rbind(cbind(ZtWZ(X.Re,w.resid), crossprod(X.Re,sweep(ZAL,MARGIN=1,w.resid,`*`))),
                                cbind(crossprod(ZAL,sweep(X.Re,MARGIN=1,w.resid,`*`)), ZtWZ(ZAL,w.resid)))
        }
      } 
    }
    if (models[[3]]=="phiHGLM") {
      mess <- pastefrom("correction needed for p_bv for phi DHGLMs.")
      stop(mess)
    } else hv10<-0 ## code cleanup 20/01/13
    if (models[[3]]=="lamHGLM") {
      mess <- pastefrom("correction needed for p_bv for lambda DHGLMs.")
      stop(mess)
    } else hv20<-0 ## idem
    #### distinct handling of AIC and p_bv (L-BFGS-B requires a non trivial value):
    if (is.nan(ladbv) || AIC ) { 
      eigvals <- eigen(d2hdbv2/(2*pi),only.values = T)$values
    }
    if (is.nan(ladbv)) { 
      zut <- abs(eigvals) ## 05/01/13
      zut[zut<1e-12] <- 1e-12
      ladbv <- sum(log(zut)) ## 
    }
    p_bv <- hlik-(hv10+hv20+ladbv/2)  
    if ( ! is.null(calcpv$second.corr)) p_bv <- p_bv + calcpv$second.corr
    if ( AIC ) {
      # a debugging issue is that options(error=recover) acts before tryCatch gets the return value
      # from its first argument. So a tryCatch on solve is not a good idea.
      if (min(eigvals)>1e-12) {
        qr.d2hdbv2 <- QRwrap(d2hdbv2)
        pd <- sum(diag(solveWrap.matrix(qr.d2hdbv2,d2clikdbv2,stop.on.error=stop.on.error)))
        if (class(pd)=="try-error") {
          warning("Computation of cAIC failed because the 'd2hdbv2' matrix appears singular")
          pd <- NA
        }
      } else pd <- Inf
      ## eqs 4,7 in HaLM07
      cAIC <- -2*clik + 2*pd ## 
      # there is also a "focussed" AIC in HaLM07 that would be   
      # - 2 p_bv + 2 * <number of dispersion parameters> (ie lambda,phi,nu,rho...)
      ## that would be used to select among different dispersion models 
      # discussion in section 7 of the paper suggests using an AIC based on p_v for selection among different fixed effect component models
      # but Yu and Yau then suggest caicc <- -2*clik + ... where ... involves d2h/db d[disp params] and d2h/d[disp params]2
    } else {cAIC <-NULL}
    if (models[[1]] != "etaHGLM") {
      APHLs <- list(p_v=ml,p_bv=p_bv) ## FR->FR rename ?
    } else APHLs <- c(calcpv,list(p_bv=p_bv))
    APHLs$cAIC <- cAIC
  }
  ###################
  ###################
  ## BUILD RETURN VALUE
  ###################
  ###################
  #
  ###################
  ## LIKELIHOODS
  ###################
  res<-list(APHLs=APHLs)
  ###################
  ## DATA
  ###################
  res$data <- data ## very useful for simulate...
  if (family$family=="binomial") {
    res$weights <- BinomialDen
  }
  res$y <- y ## counts for Pois/bin
  res$prior.weights <- prior.weights ## see Gamma()$simulate
  ###################
  ## MODEL info
  ###################
  res$family <- family
  res$X <- X.pv
  res$ranFix <- ranFix ## currently as a uniform template consistent with projected changes ; excpt that lamFix, phiFix info is now in lambda.object, etc
  corrFixNames <- intersect(names(ranFix),c("nu","rho","Nugget","ARphi"))
  corrPars <- ranFix[corrFixNames] ## idem. cf functions in corrMM.LRT that always look in phi, lambda, rather than .Fix. 
  typelist <- list()
  typelist[corrFixNames] <- "fix"
  corrPars[corrNames] <- corr_est[corrNames]
  typelist[corrNames] <- "var"
  attr(corrPars,"type") <- typelist
  res$corrPars <- corrPars 
  ## FR->FR il serait logique ? de regrouper $ranFix et $corrPars dans la sortie ?
  res$models <- models
  res$predictor <- predictor ##  all post fitting functions expect PROCESSED predictor 
  res$ZALMatrix <- ZAL ## ZAL used by simulate.HL (the $LMatrix is in the $predictor)...
  if (models[[1]] == "etaHGLM") res$ZAlist <- ZAlist ## ... but we more generally need ZA for prediction variance
  res$REMLformula <- REMLformula
  ###################
  ## ALGORITHM
  ###################
  res$HL <- HL ## info on fitting method
  ###################
  ## FITTED VALUES
  ###################
  if (family$family=="binomial") {
    res$fv <- mu/BinomialDen ## cf glm(binomial): fitted values are frequencies 
  } else {res$fv <- mu} ## fitted values may be counts (cf poisson), or reals
  ###################
  ## FIXEF, ETA
  ###################
  if (!is.null(beta_etaOri)) {
    names(beta_etaOri) <- colnames(XpvOri) 
    beta_etaOri[names(beta_eta)] <- beta_eta ## keeps the original NA's
    res$fixef <- beta_etaOri
    se <- beta_etaOri ## also put NA's
    se[names(beta_se)] <- beta_se 
    res$fixef_se <- se
  } else {
    names(beta_eta) <- colnames(X.pv)
    res$fixef <- beta_eta
    res$fixef_se <- beta_se
  }
  res$eta <- eta ## convenient for defining starting values...
  ###################
  ## LEVERAGES and REML (ie either phi OR lambda was estimated)
  ###################
  if (HL[1]!="SEM") { ## both lev_phi and deviance_residual missing otherwise
    if (is.null(phi.Fix) || is.null(lambda.Fix)) { ## in either case all leverages are computed and it makes sense to consider the residuals
      res$lev_phi <- lev_phi 
      res$std_dev_res <- sign(y-mu) * deviance_residual*prior.weights/(phi_est*(1-lev_phi)) ## should all have variance 1
    }
    if (is.null(lambda.Fix)) res$lev_lambda <- lev_lambda
  }  
  if (HL[1]!="SEM") {
    res$RespLink_disp <- RespLink_disp
  }
  if ( distinct.X.ReML ) res$X.Re <- X.Re
  ###################
  ## ALL other LAMBDA returns
  ###################
  res$rand.families <- rand.families 
  ##
  res$ranef <- u_h
  res$v_h <- v_h
  if (nrand>0) {
    print_lambda <- lapply(seq(nrand), function(it) {
      if (models[[2]][it]=="lamScal") {
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        unique(lambda_est[u.range])        
      } else {
        print_lambda <- lambda_est ## pseudocode
      }
    })
    if (all(models[[2]]=="lamScal")) print_lambda <- unlist(print_lambda) ## format for easy display... but also used by simulate...
    attr(print_lambda,"n_u_h") <- vec_n_u_h
  } else {
    print_lambda <- lambda_est <- NULL
  }  
  res$lambda <- print_lambda
  res$fittedLambda <- lambda_est
  if (models[[1]]=="etaHGLM") {
    if (is.null(lambda.Fix)) {
      if (HL[1]=="SEM") {
        lambda.object <- list(linkscale.lambda=log(lambda_est),lambda_se=NA,namesX_lamres="(Intercept)",namesTerms=NULL)
      } else {
        namesTerms <- attr(ZAlist,"namesTerms") ## for each random term, the names of the coefficients fitted
        namesnames <- unlist(lapply(names(namesTerms),function(st) {
          if (nchar(st)>10) st <- paste(substr(st,0,9),".",sep="")
          st
        }))
        names(namesTerms) <- make.unique(namesnames,sep=".") ## makes group identifiers unique (names of coeffs are unchanged)
        lambda.object <- list(linkscale.lambda=linkscale.lambda,lambda_se=lambda_se,
                              namesX_lamres=colnames(X_lamres),namesTerms = namesTerms)
        attr(lambda.object,"warning") <- resglm_lambda$warning$message ## may be NULL
        if (! is.null(next_LMatrices)) {
          res$cov.mats <- lapply(next_LMatrices,function(mat) {
            ZWZt(attr(mat,"Lcompact"),exp(linkscale.lambda))
          })
        }
      }
    } else {
      lambda.object <- list(lambda.fix=lambda.Fix)
    }
    res$lambda.object <- lambda.object
  }
  ###################
  ## ALL other PHI returns
  ###################
  if ( is.null(phi.Fix) ) {
    res$resid.predictor <- resid.predictor
  } else {
    res$resid.predictor <- NULL
  }
  res$phi <- phi_est
  if (is.null(phi.Fix)) {
    res$phi.object <- list(linkscale.phi=linkscale.phi,phi_se=phi_se,namesX_disp=colnames(X_disp))
  } else {
    res$phi.object <- list(phi.Fix=phi.Fix)
  }
  ###################
  ## something to be checked
  ###################
  if (max.iter<1) return(res) ## FR->FR  !! NOT class==HLfit !! FR->FR 06/2014: no longer clear what for 
  ##The idea was to put below this line everything that requires the computation of a fit as defined by max.iter
  ## now, either only a likelihood is computed, or the first iteration of the main loop *before* the test on max.iter (late modif of code ?) may provide much of what follows...
  ###################
  ## private hack
  ###################
  #    if ( ! is.null(init.HLfitName)) {
  if ( ! is.na(spaMM.getOption("INIT.HLFITNAME"))) {
    nextinit.HLfit <- list()
    nextinit.HLfit$fixef <- beta_eta
    nextinit.HLfit$v_h <- v_h
    if (is.null(lambda.Fix)) nextinit.HLfit$lambda <- lambda_est
    spaMM.options(INIT.HLFITNAME=nextinit.HLfit)
    ##assign(init.HLfitName, nextinit.HLfit,pos=".GlobalEnv")
  }  
  ###################
  ## WARNINGS
  ###################
  ## translation of warnings in user-more friendly form ##FR -> FR  a revoir
  if ( ! is.null(warningList$resLam0) && warningList$resLam0) { 
    warningList$resLam0 <- "lambda residuals numerically 0 were replaced by 1e-6"
  }
  if ( ! is.null(warningList$resLamInf) && warningList$resLamInf) { 
    warningList$resLamInf <- "lambda residuals numerically >1e10 were replaced by 1e10"
  }
  if (! is.null(warningList$leveLam1) && warningList$leveLam1) {
    warningList$leveLam1 <- "lambda leverages numerically 1 were replaced by 1 - 1e-8"
  }
  if ( ! is.null(warningList$resPhi0) && warningList$resPhi0) { 
    warningList$resPhi0 <- "phi residuals numerically 0 were replaced by 1e-6"
  }
  if ( ! is.null(warningList$resPhiInf) && warningList$resPhiInf) { 
    warningList$resPhiInf <- "phi residuals numerically >1e10 were replaced by 1e10"
  }
  if (! is.null(warningList$levePhi1) && warningList$levePhi1) {
    warningList$levePhi1 <- "phi leverages numerically 1 were replaced by 1 - 1e-8"
  }
  if (! is.null(warningList$negLevLam) && warningList$negLevLam) {
    warningList$negLevLam <- "Negative leverages for lambda were replaced by 1e-8"
  }
  if ((HL[1]!="SEM") && maxit.mean>1 && 
        ((models[[1]]=="etaHGLM") || ((models[[1]]=="etaGLM") && pforpv>0)) ## cases where iterations are needed
      && innerj==maxit.mean) {
    warningList$innerNotConv <- paste("linear predictor estimation did not converge. Try increasing 'max.iter.mean' above ",maxit.mean,sep="")
  }
  if (iter==max.iter) {
    warningList$mainNotCov <- paste("Joint estimation of fixed effects and dispersion parameters \n  did not converge. Try increasing 'max.iter' above ",max.iter,sep="")
  }
  res$warnings <- warningList
  ###################
  ## SUMMARY, RETURN
  ###################
  class(res) <- c("HLfit",class(res)) 
  if (verbose["summary"]) {
    summary(res) 
  }
  if (verbose["warn"]) {
    seriousWarnings <- warningList[intersect(c("innerNotConv","mainNotCov"),names(warningList))]
    if (length(seriousWarnings)>0 ) { 
      abyss <- sapply(length(seriousWarnings),function(i) {warning(paste("In HLfit :\n",seriousWarnings[[i]],sep=""),call.=FALSE)}) 
      warningList[setdiff(names(warningList),c("innerNotConv","mainNotCov"))] <- NULL
    }
  } 
  if (verbose["trace"]) {
    if (length(warningList)>0 ) {
      abyss <- sapply(length(warningList),function(i) {cat(warningList[[i]]);cat("\n")}) 
    }
  }
  res$call <- mc
  return(res)
}
