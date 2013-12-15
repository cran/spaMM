HLfit <-
function(predictor,
                  data,family=gaussian(),rand.family=gaussian(), 
                  resid.predictor = ~ 1 ,REMLformula=NULL,
                  verbose=c(F,F),
                  HLmethod="HL(1,1)",
                  control.HLfit=list(),
                  init.HLfit = list(), 
                  ranFix=list(), ## phi, lambda
                  etaFix=list(), ## beta, v_h (or even u_h)
                  processed=NULL
                 ) {
    ###################################################
    mc <- match.call()
    cstForBinary <- 0.34584297349917959320 ## ((16 \[Sqrt]3)/(15 \[Pi]))^2 ## not for normal use
    warningList<-list()
    if (length(verbose)==1) verbose <- c(verbose,F)
    lambda.Fix <- ranFix$lambda
    if (any(lambda.Fix==0)) {
      mess <- pastefrom("lambda cannot be fixed to 0.",prefix="(!) From ")
      stop(mess)
    }
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
      processed <- preprocess(control.HLfit=control.HLfit,HLmethod=HLmethod,predictor=predictor,phi.Fix=phi.Fix,
                                                 resid.predictor=resid.predictor,REMLformula=REMLformula,data=data,family=family,BinomialDen=BinomialDen,rand.family=rand.family)
    } 
    vUpdating <- processed$vUpdating
    LevenbergM <- processed$LevenbergM
    spam <- processed$spam 
    stop.on.error <- processed$stop.on.error ## to control issues with matrix computations; F by default
    AIC <- processed$AIC ## whether to compute any AIC stuff; F by default
    essai <- processed$essai ## to control any tested new code...
    conv.threshold <- processed$conv.threshold
    iter.mean.dispFix <- processed$iter.mean.dispFix
    iter.mean.dispVar <- processed$iter.mean.dispVar
    max.iter <- processed$max.iter
#    ps_v.threshold <- processed$ps_v.threshold
    resid.predictor <- processed$resid.predictor 
    REMLformula <- processed$REMLformula
    X.Re <- processed$X.Re
    ### a bit of post processing
    if ( ! is.null(REMLformula) ) { ## differences affects only REML estimation of dispersion params, ie which p_bv is computed
       adjREML <- T
    } else {
       adjREML <- F
    }
    ###
    HL <- processed$HL
    BinomialDen <- processed$BinomialDen
    loglfn.fix <- processed$loglfn.fix
    X.pv <- processed$X.pv
    ### a bit of post processing
    n <- nrow(X.pv)
    pforpv <- ncol(X.pv)
    ###
    y <- processed$y
    canonicalLink <- processed$canonicalLink
    ### a bit of post processing
    if (family$family=="gaussian" && rand.family$family=="gaussian" && canonicalLink) {
      LMMbool <- TRUE
    } else LMMbool <- FALSE
    ###
    models <- processed$models
    #### here problem that  HLCor modifies the L matrix => it cannot be preprocessed 
    if (models[["eta"]]=="etaHGLM") { ## Design matriCES for random effects in particular, prob only a match between the levels or the ranef and the observ. Ie Z, not ZL 
      namesRE <- processed$namesRE
      Zlist <- processed$Zlist ## : Zlist is a list of design matrices 
      if ( ! is.null(predictor$ZLMatrix)) { ## if ZL directly given
         ZL <- predictor$ZLMatrix
      } else { ## reconstruct ZL from Z (Z from spMMFactorList, L from user)
        ## we need to check for user's confusion when we multiply Z by LMatrix
        ## ZL is nobs * (# levels ranef) and Z too
        ## LMatrix is (# levels ranef) * (# levels ranef)
        ## the levels of the ranef must match eachother in the two matrices
        ## the only way to check this is to have the levels as rownames and colnames and to check these
        ZL <- list()
        LMatrix<-predictor$LMatrix
        if ( is.null(LMatrix)) {
            ZL[[1]] <- Zlist[[1]]
        } else if ( ( ! is.null(LMatrix)) && length(Zlist)>0) {
          if (class(LMatrix)=="matrix") {
            ## Either the user provided a distm. No unique geographic location has been determined. nrow(LMatrix) should be nrow(data)
            ## or distm and LMatrix are constructed from unique geographic locations. Fewer possibilites for errors
            ## nrow(LMatrix) should be ncol(Zlist) in both cases
            if (ncol(Zlist[[1]])!=nrow(LMatrix)) {
               warning(paste("The number of levels of the grouping variable in random term (...|",namesRE[1],")",sep=""))
               warning("  is not the dimension of the correlation matrix (This error shoud have been caught before).") ## by distm checking in corrHLfit or no.info check somewhere...
               stop("I exit.")
            } 
            ### it's difficult to make checks on names at this step
            ## rownames(LMatrix) are the names of first occurrences of unique geographic locations, or (if user provided distm) whatever was in this distm.
            ## Zlist inherits anything from the spMMFactorList call that is hard to manipulate 
            ZL[[1]] <- Zlist[[1]] %*% LMatrix
          } else if (class(LMatrix)=="list") {
            LMlen <- length(LMatrix)
            for (ii in min(c(length(Zlist),LMlen))) {
               if ( ! all(colnames(Zlist[[ii]])==rownames(LMatrix[[ii]]))) {
                  stop("The colnames of the design matrix Z in eta=...+Zv should be the rownames of the design matrix L  in v=Lu")
               }
               ZL[[ii]] <- Zlist[[ii]] %*% LMatrix[[ii]]
            }
          } 
        }
      } ## end case reconstruction of ZL from Z and L
    } else { ## models[["eta"]] = "etaGLM"
      ZL <- NULL
      ## there are really to cases here (1) true GLM with constant phi, eta (2) phi(H)GLM. In both case we could need to bypass the augmented GLM fit
      mess <- pastefrom("no random effect in model formula. More code needed to bypass the augmented linear model step.",prefix="(!) From ")
      stop(mess)
    } ## 
    ####
    ### a bit of post processing // repeat of code in translate...
    nrand <- length(ZL)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- ncol(ZL[[i]]) ## nb cols each design matrix = nb realizations each ranef
    qcum <- cumsum(c(0, q)) ## if two ranef,  with q=(3,3), this is 0,3,6. qcum[nrand+1] is then 6, the total # of realizations
    if (nrand==1) ZL <- ZL[[1]] ## temporary hack...
    I <- diag(rep(1,qcum[nrand+1L]))
    ZLI <- rbind(ZL,I) ## used where ?
    ###
    X_lambda <- processed$X_lambda
    X_disp <- processed$X_disp ## may be NULL
    off <- processed$off
    RandDist <- processed$RandDist
    ##################
    RespLink_lambda <- "log"
    if (is.character(init.HLfit)) { ## at this point init.HLfit is a string or not. Elsewhere it can be a list
       spaMM.options(INIT.HLFITNAME=init.HLfit) ## if a string, copied in...
       # init.HLfitName <- init.HLfit ## if a string, copied in 
       ### init.HLfit <- get(init.HLfitName) ## not a $-member !
       #print(init.HLfit)
       #init.HLfit <- init.HLfit$latestVals
    } else {
      spaMM.options(INIT.HLFITNAME=NA)  
      # init.HLfitName <- NULL
      unknowns <- names(init.HLfit)[!names(init.HLfit) %in% c("fixef","phi","lambda","v_h")] 
      if (length(unknowns)>0) {
        mess <- pastefrom("unhandled elements in 'init.HLfit'.",prefix="(!) From ")
        message(mess)
        if ("beta" %in% unknowns) message("  Use 'fixef' rather than 'beta' in 'init.HLfit'.")
        stop()
      }
    }
    
#################### LOCAL FNS DEFS ###################################
`calc.p_v` <- function() { 
  theta <- theta.mu.canonical(mu/BinomialDen,family$family)  
  if (family$family=="binomial") {
    clik <- sum(loglfn.fix(theta,y/BinomialDen,BinomialDen,1/phi_est)) ## freq a revoir
  } else {
    clik <- sum(loglfn.fix(theta,y,1/phi_est)) 
  }
  if (models[1]!="etaHGLM" && models[3]!="phiHGLM") return(list(clik=clik))
  ## ELSE
  likranU <- loglfn.ranU(RandDist,u_h,1/lambda_est) 
  log.du_dv <- - log(dvdu) 
  likranV <- sum(likranU + log.du_dv)
  ##### HLIK
  hlik<-clik+likranV      
  ##### P_V
  lad <-determinant(d2hdv2/(2*pi))$modulus[1] 
  if (is.nan(lad) || is.infinite(lad)){## because of determinant of nearly singular matrix
    zut <- abs(eigen(d2hdv2/(2*pi),only.values = T)$values) ## 05/01/13
    zut[zut<1e-12] <- 1e-12
    lad <- sum(log(zut)) ## L-BFGS-B requires a non trivial value
  }
  p_v <- hlik-lad/2
  resu <- list(clik=clik,hlik=hlik,p_v=p_v)
  if (family$family=="binomial" && max(BinomialDen)==1 && FALSE) { ## may efficiently restrict rho but not lambda...
    ##je ne peux pas estimer les v par unefonction qui ne depend pas des valeurs de v... Doc cane peut pas marcher...
     resu$p_vOri <- p_v
     ## compute approximation for independent binary samples
     logitp <- (off + X.pv %*% beta_eta)/sqrt(1+cstForBinary * lambda_est*colSums(ZL^2)) ## approx ... one could even consider an nintegrate here...
     logitp[logitp > 700] <- 700 ## because -> exp(700) below... 
     loglikapp <- sum(log((1 + y * (exp(logitp) - 1))/(exp(logitp) + 1)))
     loglikapp <- max(loglikapp,-exp(700)) ## avoids -Inf
     zut <- crossprod(ZL) ## takes out lambda from d2hdv2...
     hum <- eigen(zut,only.values=T)$values 
     indepness <- mean(sort(hum,decreasing=T)[-1]) ## [0,1] will reach maxi 1 for zut = I.  Uses the fact that the sum of eigenvalues remains 1
     ## [0,1] [0,1007] [exp(-1000),exp(7)=1096...]
     weiindep <- exp(1007*indepness-1000) ## exp(0)=1 if all subdom eigen are 1000/1007
     if (weiindep>0.001 && ranFix$rho<30 && lambda_est[1]>50) {browser()}
     if (weiindep>1) {cat("@")} else {cat("#")}
     resu$p_v <- (p_v+weiindep*loglikapp)/(1+weiindep)
  }
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
#      ZBZt <- ZL %*% B %*% t(ZL) ## we need the full matrix for the last term...
#      diagZBZt <- diag(ZBZt)
#      HabcdBacBbd <- - sum(d4bTh*diagZBZt^2) ## minus added 040213 ...
#      b3.diagZBZT <- as.numeric(diagZBZt * d3bTh) ## inner product vector = vector
#      acoefs <- as.numeric(b3.diagZBZT %*% ZL) ## vector which 'a'th element is Sum_i d3bTh_i * (ZBZt)_ii ZL_ia
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
        invL.ZLt <- solve(L,t(ZL)) ## L^{-1}.t(ZL)
        ZBZt <- crossprod(invL.ZLt) ## ZL.t(L)^{-1}.L^{-1}.t(ZL)
        diagZBZt <- diag(ZBZt)
        HabcdBacBbd <- - sum(d4bTh*diagZBZt^2) ## minus sign from habcd =- b''''(th) zzzz
        b3.diagZBZT <- as.numeric(diagZBZt * d3bTh) ## inner product vector = vector
        acoefs <- as.numeric(b3.diagZBZT %*% ZL) ## vector which 'a'th element is Sum_i d3bTh_i * (ZBZt)_ii ZL_ia
        HabcHrstBarBbcBst <- sum(solve(L,acoefs)^2) ## strictly, sum(solve(L, - acoefs)^2) starting from  habc=-b'''(th) zzz
        ZBZtcube <- ZBZt * ZBZt * ZBZt
        HabcHrstBarBbsBct <- d3bTh %*% ZBZtcube %*% d3bTh ## again, - d3bTh %*% ZBZtcube %*% (-d3bTh)
        second.corr <-  HabcdBacBbd/8 + HabcHrstBarBbcBst/8 + HabcHrstBarBbsBct/12 ## - F/24 (ps_bv =p_bv + second.corr = p_bv -F/24)
      }    
    } else { ## clearly the fastest code
      invL.ZLt <- forwardsolve(L,t(ZL)) ## L^{-1}.t(ZL)
      ZBZt <- crossprod(invL.ZLt) ## ZL.t(L)^{-1}.L^{-1}.t(ZL)
      diagZBZt <- diag(ZBZt)
      HabcdBacBbd <- - sum(d4bTh*diagZBZt^2)
      b3.diagZBZT <- as.numeric(diagZBZt * d3bTh) ## inner product vector = vector
      acoefs <- as.numeric(b3.diagZBZT %*% ZL) ## vector which 'a'th element is Sum_i d3bTh_i * (ZBZt)_ii ZL_ia
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

if ( length(etaFix$beta)==ncol(X.pv) && !is.null(etaFix$v_h) && !is.null(phi.Fix) && !is.null(lambda.Fix) ) { ## nothing to fit. We just want a likelihood
  ### a bit the same as max.iter<1 ... ?
  ## we need u_h in calc.p_v() and v_h here for eta...
  v_h <- etaFix$v_h
  u_h <- etaFix$u_h
  if (is.null(u_h)) { u_h <- rand.family$linkinv(v_h)}
  lambda_est <- lambda.Fix
  phi_est <- phi.Fix
  if (models[1]=="etaHGLM") { ## linear predictor for mean with ranef
    eta <- off + X.pv %*% etaFix$beta + ZL %*% etaFix$v_h ## updated at each iteration
  } else  eta <- off + X.pv %*% etaFix$beta ## no iteration hence no updating  ## FREQS
  ## conversion to mean of response variable (COUNTS for binomial)
  blob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen) 
  mu <- blob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  w.resid<-as.vector(blob$GLMweights /phi_est) ## 'weinu', must be O(n) in all cases
  wranefblob <- updateWranef(rand.family,lambda_est,u_h,v_h)
  dvdu <- wranefblob$dvdu
  w.ranef <- wranefblob$w.ranef
  d2hdv2 <- - sweep(t(ZL),MARGIN=2,w.resid,`*`) %*% ZL - diag(w.ranef) ##  - t(ZL) %*% diag(w.resid) %*% ZL - diag(w.ranef)
  return(list(APHLs=calc.p_v())) ### RETURN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      } else endform <- "0"
      locform <- as.formula(paste(begform, endform))
      resglm<-glm(locform,family=family,offset=off) ## warning "glm.fit: fitted probabilities numerically 0 or 1 occurred" => separation or large offset
      if (pforpv>0) {
        if (max(abs(c(resglm$coefficients)))>1e10) {
           message("(!) Apparent divergence of estimates in a *glm* analysis of the data.")
           message("    Check your data for separation or bad offset values.")
           stop("    I exit.") 
        } else names(resglm$coefficients) <- colnames(X.pv) ## because this if unfortunately not the case... 
      } 
      return(resglm)
}

`compute.a1.adj` <- function() { ## needs qr.d2hdv2, ZL, stop.on.error, d2hdv2, rWW, ZLI, family, dmudeta, BinomialDen, mu, eta... ETC
          ########## HL(1,.) adjustment for mean ##################
          ## if LMM (ie resp gaussian, ranef gaussian), all coef<x> are 0 -> correction is 0
          ## if (gaussian, not gaussian) d3 nonzero
          ## if (non gaussian, gaussian), d3 zero but d1,d2 nonzero 
          # the aim of this is to compute the HL(1,.) estimate of beta 
          # and specifically the a(1) term in LeeL 12 p. 963
          # W is Sigma^-1 ; TWT = t(ZLI)%*%W%*%ZLI = ZL'.Wresid.ZL+Wranef = -d2hdv2 !
          # K2 is K2 matrix in LeeL appendix p. 4 and is D in MolasL p. 3307 
          K2 <- safesolve.qr.matrix(qr.d2hdv2,t(ZL),stop.on.error=stop.on.error) ## t(ZL) missing in some implementations... no effect with Gaussian ranefs...  
          if (class(K2)=="try-error") {
            mess <- pastefrom("problem in 'K2' computation.",prefix="(!) From ") ## cf BB 877
            warning(mess)
            K2 <- ginv(d2hdv2) %*% t(ZL)            
          }
          K2 <- - K2 ## '- try-error' causes a really silly error...
          # K1 is K1 in LeeL appendix p. 4 and is A (actually A11 for the first ranef) in MolasL p. 3307-8 
          K1 <- - ZL%*%K2
          ## Wresid in W is always O(n)
          ## P is P in LeeL appendix p. 4 and is P_R in MolasL p. 3307 
          ## looks like leverage computation, but the ZLI columns are considered instead of the full augmented design matrix 
          wAugZLI <- rWW%*%ZLI ## rWW previously computed for leverage computation
          qrwAugZLI <- qr(wAugZLI)
          Pdiag <- rowSums(qr.qy(qrwAugZLI, diag(1, nrow = nrow(wAugZLI), ncol = ncol(wAugZLI)))^2)
          ## K2 is not needed for gaussian ranef... but K1 is still required...
          # coef1 is the factor of P_ii in d1
          # coef2 is the factor between P_jj and K1 in d2
          ## We first handle the canonical link cases, where comput. of coef1 depends only on the link  
          ## here (d2mu/deta2)/(dmu/deta)^2 = (d(dmudeta)/dmu)/dmudeta where d(dmudeta)/dmu is as detailed:
          if (canonicalLink) {
             if (family$family=="gaussian") {
                   coef1 <- rep(0,n)
                   coef2 <- rep(0,n)
            } else if (family$family=="poisson") {
                   ## numerator is D[D[E^\[Eta], \[Eta]] /. {E^\[Eta] -> \[Mu]}, \[Mu]] =1 
                   coef1 <- 1/dmudeta
                   coef2 <- rep(1,n)
            } else if (family$family=="binomial") {
                  ## numerator is D[D[1/(1 + E^-\[Eta]), \[Eta]] /. {E^-\[Eta]->(1-\[Mu])/\[Mu]} ,\[Mu]]=1-2 mu 
                  coef1 <-(1-2*mu/BinomialDen)/dmudeta  
                  coef2 <-(1-2*mu/BinomialDen)  
            } else if (family$family=="Gamma") {
                   ## numerator is D[D[-1/\[Eta], \[Eta]] /. {\[Eta] -> -1/\[Mu]}, \[Mu]] =2 mu 
                   coef1 <- 2*mu /dmudeta
                   coef2 <- 2*mu
            } 
          } else { ## general code handling non canonical links
            tmblob <- thetaMuDerivs(mu,BinomialDen,family$family)
            Dtheta.Dmu <- tmblob$Dtheta.Dmu
            D2theta.Dmu2 <- tmblob$D2theta.Dmu2
            ## if the family link eta(mu) equals the canonical link theta(mu), then theta=eta, the following line is null  
            d2mudeta2 <- d2mudeta2fn(link=family$link,mu=mu,eta=eta,BinomialDen=BinomialDen)
            D2theta.Deta2_Dtheta.Deta <- (D2theta.Dmu2 * dmudeta^2 + Dtheta.Dmu * d2mudeta2)/(Dtheta.Dmu * dmudeta)
            D2mu.Deta2_Dmu.Deta <- d2mudeta2/dmudeta
            coef2 <- D2mu.Deta2_Dmu.Deta + D2theta.Deta2_Dtheta.Deta 
            coef1 <- coef2 / (Dtheta.Dmu * dmudeta^2) ## note that coef2 is indep of the BinomialDen, but coef1 depends on it 
          }
          n_u_h <- ncol(ZL)
          seqn_u_h <- (1:n_u_h)
          # coef3 =(1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
          ## This (dlogWran_dv_h) was computed when w.ranef was computed
          if (dlogWran_dv_h!=0L) coef3 <- rep(dlogWran_dv_h,n_u_h)
          Pdiagn <- Pdiag[1:n]
          d1 <- Pdiagn * coef1
          d2<-rep(0,n)
          d3<-rep(0,n)
          for (i in 1:n) {
            #####  d1[i]<-Pdiag[i]*coef1[i]
            #####  for (qq in 1:n) { d2[i]<-d2[i]+Pdiag[qq]*coef2[qq]*K1[qq,i] }
            d2[i] <- sum( Pdiagn * coef2 * K1[1:n,i] )                 
            if (dlogWran_dv_h!=0L) { ## to save some time...
              ##### for (qq in 1:n_u_h) { d3[i]<-d3[i]+Pdiag[n+qq]*coef3[qq]*K2[qq,i]} ## d3 reste nul pour gaussian ranef
              d3[i] <- sum( Pdiag[n+seqn_u_h] * coef3[seqn_u_h] * K2[seqn_u_h , i] ) ## d3 reste nul pour gaussian ranef
            }
          }
          d<-d1+d2+d3 ## cf d_i p.4 supp mat
          sscaled<-d/2  ## we had s<-d*dmudeta/2 then Sig%*%Wresid%*%s*detadmu ... so sscaled := s*detadmu
          ## up to here, a is zero if ranef are gaussian... au it should remain so if response is Gaussian ??
          da <- Sig%*% (w.resid * sscaled) ## from a(0) to a(1) LeeL12 p. 963 col 1 l 2
          return(da)
}


LevenbergMstep <- function() {
  LM_wAugz <- wAugz - wAugX %*% betaV
  ## verif correct rhs: verif_dbetaV <- safesolve.qr.vector(SQR, LM_wAugz,stop.on.error=stop.on.error)
  rhs <- crossprod(wAugX,LM_wAugz) ## t(wAugX) %*% LM_wAugz ## backward: A'.z* ## rhs needed further below to assess gainratio
  ApAdDpD <- crossprod(wAugX) ## t(wAugX) %*% wAugX ## A'.A=R'.R 
  dampDpD <- damping * diag(ApAdDpD) ## faster than computing qr.R(SQR ) and using it specially for this computation... ## vector
  diag(ApAdDpD) <- diag(ApAdDpD) + dampDpD  ## A'.A + damping*diag(A'.A)
  ## ## code that oddly solved B y = rhs= A'z* for B:=(A'.A + damping*diag(A'.A)) by solving B'B y =B'(A'z*) is in version 040213. For damping=0, this solved A.A'.A'.A y= A.A'.A'z*...
  ### attempt to solve 'normal equations' directly ## LMcorrected useful only here or for the ginv() call
  dbetaV <- try(solve(ApAdDpD,rhs),silent=TRUE)
  if (class(dbetaV)=="try-error") {
    ### then QR, using a standard trick: the solution of (A'A+damping D'D)y = -A'f is the solution of extended system (A // sqrt(damping) D) y = -( f // 0)
    corrD <- sqrt(dampDpD) ##  D'.D = damping* diag(A'.A) (computing diag A'A through qrR(A) is slower) ## vector
    trick <- rbind(wAugX,diag(corrD)) ## matrix
    dbetaV <- safesolve.qr.vector(qr(trick),c(LM_wAugz,rep(0,length(corrD))),stop.on.error=stop.on.error) ## 
    ## and there is nmore to this as explained by Mor'e and illustrated by the following code
#set.seed(123)
#m <- matrix(runif(50),ncol=5)
#J <- Matern.corr(m,nu=1);qrJ<- qr(J);qrJ$pivot
#b <-rnorm(10)
#H <- t(J) %*% J
#D <- sqrt(diag(H))
#qrP <- qr(rbind(qr.R(qrJ),diag(sqrt(0.001)*D))) ## that's the RHS of (3.4)
#tW <- qr.Q(qrP) # 3.5 de MOr\'e means tW.(R_lam //0) = (R//D_lam)
#u <- t(tW) %*% c(t(qr.Q(qrJ)) %*% b,rep(0,5)) ## Mor\'e's Q ist(Q)
#solve(qr.R(qrP),u[1:5])
# = :
#solve(H+0.001*diag(D^2) ,t(J) %*% b)
    ## but it misses the permutations to be generally valid !!!!!!!!!!!!!!!!!!!!!!          
    }
    if (class(dbetaV)=="try-error") {
      dbetaV <- ginv(ApAdDpD) %*% rhs
    }
    return(list(dbetaV=dbetaV,denomGainratio=sum(dbetaV*(rhs+ dampDpD * dbetaV))))
}

##############################################################################################
    ######### Initial estimates for mu by GLM ####################
    if ( ( pforpv>0 && is.null(init.HLfit$fixef)) || is.null(init.HLfit$v_h) || is.null(phi.Fix) || is.null(lambda.Fix)) { 
       ## all cases where an initial resglm is needed (even when pforpv=0, may be needed to provide init phi or init lambda)
       resglm <- provide.resglm()   
    }
    if (pforpv>0) { 
      beta_eta <- numeric(ncol(X.pv))
      beta_eta <- init.HLfit$fixef
      if (is.null(beta_eta) ) {
         beta_eta<-c(resglm$coefficients)[1:pforpv] 
      } 
      beta_eta[names(etaFix$beta)] <- etaFix$beta
    } else {
      beta_eta <- numeric(0)
      se_beta <-numeric(0)
    }
    ## Initial estimate for phi ####
    if (is.null(phi.Fix)) { ## at this point, means that not poisson nor binomial
      phi_est <- init.HLfit$phi
      if (is.null(phi_est) ) {
        phi_est <- as.numeric(deviance(resglm)/resglm$df.residual)
      } 
    } else {
      phi_est <- phi.Fix
    }
    ##
    ######## initialize v_h #############################
    if (models[1]=="etaHGLM") { ## the basic case (LMM, GLMM...)
       if (is.null(vUpdating)) {vUpdating <- FALSE}
       if (is.null(LevenbergM)) {
         if (vUpdating) {
           LevenbergM <- FALSE 
         } else {
           if (LMMbool) {  
             LevenbergM <- FALSE
           } else LevenbergM <- TRUE
         }
       }
       if (is.null(spam)) {
         spam <- FALSE
       }
       psi_M <- switch(RandDist, 
          gaussian = 0,
          gamma = 1, 
          beta = 1/2, 
          "inverse.gamma" = 1
       )
       v_h <- etaFix$v_h
       if (is.null(v_h) ) v_h <- init.HLfit$v_h
       if (is.null(v_h) ) {
         ## rand.family$family gives the distrib of u; the mean of u is psi_M
         vv <- rand.family$linkfun(psi_M) ## v as link(mean(u))
         v_h<-rep(vv,qcum[nrand+1L])
       } 
       u_h <- rand.family$linkinv(v_h)
       lambda.family <-Gamma(link=RespLink_lambda)
       ## init values for lambda
       if (is.null(lambda.Fix)) { ## uses a default value inspired from the 'hglm' package
          init.lambda <- init.HLfit$lambda
          if (is.null(init.lambda) ) {
             resdisp <-as.numeric(deviance(resglm)/resglm$df.residual) ## FR->FR cf overdispersion... better initial lambda values ?
             #if (family$family=="gaussian" && ! is.null(phi.Fix)) { ## a revoir...
             #  init.lambda <- resdisp - phi.Fix
             #  if (init.lambda<1e-3) init.lambda <- resdisp/5 
             #} else {
               init.lambda <- resdisp/5 
             #}
          } 
       } else init.lambda <- lambda.Fix
       if  (length(init.lambda)==1) {
         lambda_est<-rep(init.lambda,qcum[nrand+1L]) ## lambda was originally a matrix with ncol= nb real ranef
       } else lambda_est <- init.lambda
    }
    if (models[3]=="phiHGLM") {
       stop("random effects in predictor or residual variance (phi) not yet implemented")
       ## there is a buggy template code with comments in version 260812 of HLfit
    }
#print(paste("ZL[1,1] in HLfit ",ZL[1,1]))
## predictor from initial values
  if (models[1]=="etaHGLM") { ## linear predictor for mean with ranef
    eta <- off + X.pv %*% beta_eta + ZL %*% v_h ## updated at each iteration
  } else  eta <- off + X.pv %*% beta_eta ## no iteration hence no updating  ## FREQS
  ## conversion to mean of response variable (COUNTS for binomial)
  blob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen) 
  mu <- blob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  dmudeta<- blob$dmudeta ## if Bin/Pois, must be O(n)
  Vmu <- blob$Vmu ## if Bin/Pois, O(n)
  w.resid<-as.vector(blob$GLMweights /phi_est) ## 'weinu', must be O(n) in all cases
  if (models[1]=="etaHGLM") {
    wranefblob <- updateWranef(rand.family,lambda_est,u_h,v_h)
    w.ranef <- wranefblob$w.ranef
    dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
    dvdu <- wranefblob$dvdu
  }
  betaV <- c(beta_eta,v_h) 
  conv.phi <- FALSE; conv.lambda <- FALSE
  if (models[1]=="etaHGLM") {
    Sig<- sweep(ZL,MARGIN=2,1/w.ranef,`*`)  %*% t(ZL) + diag(1/w.resid) ## ZL %*% diag(1/w.ranef) %*% t(ZL) + diag(1/w.resid) 
    d2hdv2 <- - sweep(t(ZL),MARGIN=2,w.resid,`*`) %*% ZL - diag(w.ranef) ##  - t(ZL) %*% diag(w.resid) %*% ZL - diag(w.ranef)
    qr.d2hdv2 <- NULL
    OO1<-matrix(0,qcum[nrand+1L],pforpv)
    TT <- rbind(cbind(X.pv,ZL),cbind(OO1,I))  ## aug design matrix
    if (length(etaFix$beta)==ncol(X.pv) && !is.null(etaFix$v_h)) {
      maxit.mean <- 0 ## used in test near the end...
    } else if ( LMMbool ) {
      if ( LevenbergM ) { ## only for testing purposes
        maxit.mean <- 3 
      } else maxit.mean <- 1 ## sufficient for LMM 
    } else { ## even h maximization in *G*LMMs 
      if ( ! is.null(phi.Fix) && ! is.null(lambda.Fix)) { ## allFix hence no true outer iteration 
        maxit.mean <- iter.mean.dispFix 
      } else maxit.mean <- iter.mean.dispVar # If phi.Fix and lambda.Fix, the only way to have 'outer' convergence is to have 'inner' convergence
    } 
  }
  iter<-0
  ########################################
  ######### Main loop ####################
  ########################################
  while ( TRUE ) { ##the main loop with steps: new linear predictor, new leverages, new phi, new w.resid, new lambda, new fn(lambda)
    if (models[1]=="etaHGLM") {
      if (LevenbergM) {
        damping <- 1/1000000 ## as suggested by Madsen-Nielsen-Tingleff... ## mauvais resultats si on part + haut
        dampingfactor <- 2
      }
      ##################
      ### Inner Loop ### IF random effect *in mean predictor*: estim beta,v [only v if pforpv=0] for given phi,lambda
      ################## 
      for (innerj in 1:maxit.mean) {  ## breaks when conv.threshold is reached
        w <- sqrt(c(w.resid,w.ranef)) ## if maxit.mean>1 GLM weights have been changed and the following must be recomputed
        rWW <- diag(w) ## needed early, for a1.adj... 
        wAugX <- rWW%*%TT 
        SQR <- qr(wAugX)
        old_betaV <- betaV
        if (vUpdating) { ######## new v_h estimates in the vUpdating case
          ############### all random effect models are canonical conjugate except the inverse.Gamma(log) ############### 
          dlogfthdth <- (psi_M - u_h)/lambda_est ## the d log density of th(u)
          if (RandDist=="inverse.gamma" && rand.family$link=="log") { ## g(u) differs from theta(u)
            ## general: dlog f(v)/dv = d[log (f(th)*dth/dv)]/d v = d[log f(th)]/d v + d[log (dth/dv)]/d v 
            ## =    dlogfthdth * (dth/dv) + d[log (dth/dv)]/d v
            dfdv <- dlogfthdth / u_h - 1
            ## other non conjugate canonical links not handled anywhere in this code
          } else { ## canonical conjugate 
            dfdv <- dlogfthdth 
          } ## dfdv is f(v) PART of dhdv for link v=g(u) equal to the conjugate canonical link theta(u); used explicity only for vUpdating, but implicitly otherwise
          ##### (y|v) PART of dhdv
          ## mu is predictor of counts ie same scale as y
          dcdv <- t(t(w.resid * (y - mu)/dmudeta) %*% ZL) ## same as t(ZL)%*% (w.resid*(y-mu)/dmudeta) but only vectors are transposed ## NB w.resid/dmudeta = 1/Phi=1 in poisson case
          dhdv <- dcdv + dfdv
          if (is.null(qr.d2hdv2)) qr.d2hdv2 <- try(qr(d2hdv2)) ## FR->FR: no code if the try fails
          dv_h <- safesolve.qr.vector(qr.d2hdv2,dhdv,stop.on.error=stop.on.error)
          if (class(dv_h)=="try-error") {
            mess <- pastefrom("problem in dv_h computation. Try option 'control.HLfit=list(vUpdating=FALSE)'",prefix="(!) From ")
            stop(mess)
          }
          v_h <- v_h-dv_h ## new estimate of random effect
          u_h <- rand.family$linkinv(v_h)
        } 
        ######## According to 'theorem 1' of LeeL12
        ## new beta estimate from z1-a(i)
        ## where z1 is
        z1 <- eta+(y-mu)/dmudeta-off ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
        ## and a(i) (for HL(i,1)) is a(0) or a(0)+ something
        ## and a(0) depends on z2, as follows :
        if (RandDist=="gaussian") {
             z2 <- rep(0,length(w.ranef)) 
             a <- rep(0,n)
        } else { ## HGLM: nonzero z2, nonzero a(0) ## this could perhaps make a separate block, but nevertheless sometimes qr.d2hdv2 computation...
             z2 <-v_h+(psi_M-u_h)*dvdu ## update since u_h,v_h updated (yes)
             #        nXn  .   nXn      nX'r'    'r'X'r'       'r'X'r'    'r'
             # a <- Sig %*% Wresid %*% ZL %*% solve(-d2hdv2) %*% Wranef %*% z2 ## p. 963 l. 1-2; a(0) supp mat p. 6 
             a <-  w.ranef * z2
             if (is.null(qr.d2hdv2)) { 
               a <- try(solve(d2hdv2, - a))
               if (class(a)=="try-error") {
                 qr.d2hdv2 <- qr(d2hdv2) 
                 a <- safesolve.qr.vector(qr.d2hdv2, - a,stop.on.error=stop.on.error)
               }  
             } else { ## we already have a qr, we use it
                 a <- safesolve.qr.vector(qr.d2hdv2, - a,stop.on.error=stop.on.error)
             }    
             if (class(a)=="try-error") {
               mess <- pastefrom("problem in 'a' correction computation.",prefix="(!) From ")
               stop(mess)
             }
             a <- Sig %*% ( w.resid * (ZL %*% a) ) ## a(0) in LeeL12
             a <- as.numeric(a) ## incase it became a Matrix, which oddly does not fit with z1-a below...
        }         
        ## and the 'something' for a(1) is computed as follows
        if (pforpv>0 && HL[1]>0 && (! LMMbool )) {
          if (is.null(qr.d2hdv2)) qr.d2hdv2 <- try(qr(d2hdv2)) 
          da <- compute.a1.adj()
          a <- as.numeric(a + da) ## in case it became a Matrix, which oddly does not fit with z1-a below...
        } ############ end computation a(1) end HL(1,.) specific code
        ######## new estimates (tentative if LM), with three subcases (only betaV estimates in the vUpdating case)
        if ( vUpdating) {
          XinvS <- calc.XinvS(Sig,X.pv,stop.on.error)
          rhs <-  XinvS %*% (z1-a)
          qr.XtinvSX <- qr(XinvS%*%X.pv)
          beta_eta <- safesolve.qr.vector( qr.XtinvSX , rhs ,stop.on.error=stop.on.error) 
          if (class(beta_eta)=="try-error") beta_eta <- ginv( XinvS%*%X.pv)%*% rhs 
          betaV <- c(beta_eta,v_h)
          dbetaV <- betaV - old_betaV
        } else { ### true augmented model, whether LevenbergM or not
          XinvS <- NULL
          if (RandDist=="inverse.gamma" && rand.family$link=="log") {
            Augz <- c(z1-a, z2 - 1/w.ranef) ## far fetched version of scaled gradient
          } else Augz <- c(z1-a,z2) ## would seem to be the more general code
          wAugz <- Augz*w
          if ( maxit.mean > 1 && LevenbergM) {
            bla <- LevenbergMstep()
            dbetaV <- bla$dbetaV
            betaV <- betaV + dbetaV
            denomGainratio <- bla$denomGainratio
          } else { ## simple aug lin or basic IRLS depending on maxit.mean
            ## QR appears faster than alternatives with crossprod(wAugX); see version 040213
            betaV <- safesolve.qr.vector(SQR, wAugz,stop.on.error=stop.on.error) ## qr.coef(SQR, wAugz) ## vector
            if (class(betaV)=="try-error") betaV <- ginv(wAugX)%*% wAugz ## occurred with large lambda either as 'init.HLfit', or by the iterative algo
            betaV <- as.numeric(betaV) #! utile ? LevenbergM produces a matrix anyway
            if (maxit.mean>1) dbetaV <- betaV - old_betaV
          } ## endif LevenbergM else...
          if (pforpv>0) {
            beta_eta <- betaV[seq(pforpv)]
            names(beta_eta) <- colnames(X.pv)
            beta_eta[names(etaFix$beta)] <- etaFix$beta
            if (is.null(etaFix$v_h)) v_h <- betaV[-seq(pforpv)] 
          } else {if (is.null(etaFix$v_h)) v_h <- betaV}
          u_h <- rand.family$linkinv(v_h)
        } ## end true augmented model       
#       print(paste(innerj," ",paste(beta_eta,collapse=" ")),quote=F)
        ####### new values of everything, only tentative if LevenbergM  
        eta <- off + X.pv %*% beta_eta + ZL %*% v_h ## updated at each inner iteration
        ## update functions u_h,v_h
        keep_wranefblob <- wranefblob #######################
        wranefblob <- updateWranef(rand.family,lambda_est,u_h,v_h)
        w.ranef <- wranefblob$w.ranef
        dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
        dvdu <- wranefblob$dvdu
        keep_blob <- blob ############################
        blob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen) 
        mu <- blob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
        dmudeta<-blob$dmudeta
        Vmu <- blob$Vmu ## if Bin/Pois, O(n)
        ## update functions of v_h -> blob
        w.resid<-as.vector(blob$GLMweights /phi_est) ## 'weinu', must be O(n) in all cases
        #### update fns of v_h -> blob -> w.resid
        if ( LevenbergM ) { ## for LMM, d2hdv2 is constant // mu, hence it is constant over LevenbergMarquardt iterations
          if (pforpv>0) keep_Sig <- Sig
          keep_d2hdv2 <- d2hdv2
#          keep_qr.d2hdv2 <- qr.d2hdv2
        }
        if (pforpv>0) Sig <- sweep(ZL,MARGIN=2,1/w.ranef,`*`)  %*% t(ZL) + diag(1/w.resid) ## ZL %*% diag(1/w.ranef) %*% t(ZL) + diag(1/w.resid) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta estimate
        d2hdv2 <- - sweep(t(ZL),MARGIN=2,w.resid,`*`) %*% ZL - diag(w.ranef) ##  - t(ZL) %*% diag(w.resid) %*% ZL - diag(w.ranef)
        qr.d2hdv2 <- NULL
        ######### for LevenbergM, we check there was a progress, and restore everything otherwise
        if (LevenbergM) {      
          currentlik <- calc.p_v()$p_v             
          if (innerj==1) { ## not a good idea to use the p_v computed for all u_h are initially set to zero
            gainratio <- 1
          } else {
            gainratio <- 2*(currentlik-oldlik)/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
          }
          if (gainratio<0) { ## restore everything
            currentlik <- oldlik ## important to restore this for the next test !
            betaV <- old_betaV
            if (pforpv>0) {
              beta_eta <- betaV[seq(pforpv)]
              names(beta_eta) <- colnames(X.pv)
              v_h <- betaV[-seq(pforpv)] 
            } else v_h <- betaV
            u_h <- rand.family$linkinv(v_h) ###############  restore u_h !
            eta <- off + X.pv %*% beta_eta + ZL %*% v_h ## updated at each inner iteration
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
        ########## nearly done with one inner iteration
        if (verbose[2]) {
          print(paste("Inner iteration ",innerj,sep=""))
          print_err<-c(lambda_est=lambda_est)
          if (innerj>1) print_err <- c(norm.dbetaV=sqrt(sum(dbetaV^2)),print_err)
          print(print_err)
          print(c(beta_eta=beta_eta))
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
    } else {
      mess <- pastefrom("eta GLM meaningful but not yet handled.",prefix="(!) From ") ##FR->FR si etaGLM, phiGLM 
      stop(mess)
    }
    ##########
    if (models[1]=="etaHGLM") {
      if ( adjREML ) { ## bidouille pour LRT on REML fits
         OO1leve<-matrix(0,qcum[nrand+1L],ncol(X.Re))
         TTleve<-rbind(cbind(X.Re,ZL),cbind(OO1leve,I))
         wAugXleve <- rWW%*%TTleve
         SQRleve <- qr(wAugXleve)
         hatval <- rowSums(qr.qy(SQRleve, diag(1, nrow = nrow(wAugXleve), ncol = ncol(wAugXleve)))^2)
      } else {
         if (pforpv==0) { ## SQR not previously computed
           w <- sqrt(c(w.resid,w.ranef)) ## if maxit.mean>1 GLM weights have been changed and the following must be recomputed
           rWW <- diag(w) ## needed early, for a1.adj... 
           wAugX <- rWW%*%TT 
           SQR <- qr(wAugX)
         }
         hatval <- rowSums(qr.qy(SQR, diag(1, nrow = nrow(wAugX), ncol = ncol(wAugX)))^2)
      }
      if (any(abs(hatval) > 1 - 1e-8)) {
		hatval <- ifelse(abs(hatval) > 1 - 1e-8, 1 - 1e-8,hatval)
        warningList$leveLam1 <-T
      }
      lev_phi <- hatval[1:n] ## for the error residuals (phi)
      lev_lambda <- hatval[(n+1L):(n+qcum[nrand+1L])]  ## for the ranef residuals (lambda)
    } else {
      diag<-glm.diag(resglm) ## function from boot package
      lev_phi<-diag$h
    }
    d2mudeta2 <- d2mudeta2fn(link=family$link,mu=mu,eta=eta,BinomialDen=BinomialDen)
    if (verbose[2]) {print(paste("beta=",paste(signif(beta_eta,4),collapse=", ")),quote=F)}
    if (HL[2]>0) { ## f_i correction in NohL07 Table 7
      #### Then we CORRECT the leverages previously computed from  the hat matrix
      ## first the d log hessian / d log lambda corrections
      if (canonicalLink) {
        dlW_deta <- d2mudeta2 / dmudeta
      } else {
        ## we need to update more functions of mu...
        tmblob <- thetaMuDerivs(mu,BinomialDen,family$family)
        Dtheta.Dmu <- tmblob$Dtheta.Dmu
        D2theta.Dmu2 <- tmblob$D2theta.Dmu2
        ## ... to compute this:
        D2theta.Deta2_Dtheta.Deta <- (D2theta.Dmu2 * dmudeta^2 + Dtheta.Dmu * d2mudeta2)/(Dtheta.Dmu * dmudeta)
        dlW_deta <- d2mudeta2 / dmudeta + D2theta.Deta2_Dtheta.Deta
      }
      dlW_deta_or_v <- c(dlW_deta, rep(dlogWran_dv_h,length(u_h)) ) ## vector with n+'r' elements 
      # dlogWran_dv_h is 0 gaussian ranef; d2mudeta2 is 0 for identity link => vector is 0 for LMM
      if (any(dlW_deta_or_v!=0L)) {
        if(is.null(lambda.Fix)) {
          ############### all random effect models are canonical conjugate except the inverse.Gamma(log) ############### 
          dlogfthdth <- (psi_M - u_h)/lambda_est ## the d log density of th(u)
          if (RandDist=="inverse.gamma" && rand.family$link=="log") { ## g(u) differs from theta(u)
            neg.d2f_dv_dloglam <- dlogfthdth / u_h ## (-) (-)(psi_M-u)/lambda^2    *    lambda....  ## needed later
            ## other non conjugate canonical links not handled anywhere in this code
          } else { ## canonical conjugate 
            neg.d2f_dv_dloglam <- dlogfthdth 
          } 
          if (is.null(qr.d2hdv2)) qr.d2hdv2 <- qr(d2hdv2) 
          dvdloglamMat <- safesolve.qr.matrix(qr.d2hdv2, diag( as.vector(neg.d2f_dv_dloglam) ), stop.on.error=stop.on.error) # rXr       
          if (class(dvdloglamMat)=="try-error") {
            mess <- pastefrom("problem in dvdloglamMat computation.",prefix="(!) From ")
            warning(mess)
            dvdloglamMat <- sweep(ginv(d2hdv2),MARGIN=2,as.vector(neg.d2f_dv_dloglam),`*`) ## ginv(d2hdv2) %*% diag( as.vector(neg.d2f_dv_dloglam))
          }
          Zdv_dloglamMat <- ZLI %*% dvdloglamMat # (r+n)Xr . rXr = (r+n)Xr
          dleve <- (hatval * dlW_deta_or_v) %*% Zdv_dloglamMat # (r+n) . (r+n)Xr = r (each element is a sum over r+n terms= a trace)
          lev_lambda <- lev_lambda - as.vector(dleve)  
        }
        if(is.null(phi.Fix)) {
          dh0deta<-( w.resid *(y/BinomialDen-mu)/dmudeta )
          ## cf calcul dhdv, but here we want to keep each d/d phi_i distinct hence not sum over observations i 
          dh0dv <- sweep(t(ZL),MARGIN=2,as.vector(dh0deta),`*`) ## dh0dv <- t(ZL) %*% diag(as.vector(dh0deta)) ## nXr each ith column is a vector of derivatives wrt v_k
          if (is.null(qr.d2hdv2)) qr.d2hdv2 <- qr(d2hdv2) 
          dvdlogphiMat <- safesolve.qr.matrix(qr.d2hdv2, dh0dv , stop.on.error=stop.on.error)  # rXn       
          if (class(dvdlogphiMat)=="try-error") {
            mess <- pastefrom("problem in dvdlogphiMat computation.",prefix="(!) From ")
            stop(mess) ## warning + ginv for lambda... !
          }
          Zdv_dlogphiMat <- ZLI %*% dvdlogphiMat # (r+n)Xr . rXn = (r+n)Xn
          dleve <- (hatval * dlW_deta_or_v) %*% Zdv_dlogphiMat # (r+n) . (r+n)Xn = n (each element is a sum over r+n terms= a trace)
          lev_phi <- lev_phi - as.vector(dleve)  
        } 
      }
    }
    if (HL[2]>1) {
      stop("Need a_i correction in Table 7 of NohL07 ie derivatives of second order correction wrt dips param.")
    }
    if (HL[3]!=0 && is.null(lambda.Fix) ) { ## ie , p_bv(h), not EQL p_bv(q+), LeeNP p89  ## PQL ? no formal def for non Gaussian ranefs...  
      ## d h/ d !log! lambda correction (nul for gaussian ranef)
      notEQL <- switch(RandDist, 
                        gaussian=0,
                        gamma=1+2*(log(lambda_est)+digamma(1/lambda_est))/lambda_est,## cf notes on p. 89 of the book
                        "inverse.gamma"=1+2*(log(lambda_est)-lambda_est+digamma(1+(1/lambda_est)) )/lambda_est, ## appears to be the same as for the gamma case... 
                        beta=1-2*(digamma(1/lambda_est)/lambda_est)+2*(digamma(1/(2*lambda_est))/lambda_est)+log(4)/lambda_est
                      ) ## as in HGLMMM PoissonEst.r
      lev_lambda <- lev_lambda + notEQL   
    }
    if (HL[3]!=0 && family$family=="Gamma" &&is.null(phi.Fix) ) { ##     
      ## d h/ d !log! phi correction (zero for gaussian residual error, no phi for other distribs)
      notEQL <- 1+2*(log(phi_est)+digamma(1/phi_est))/phi_est ## as in HGLMMM IWLS_Gamma; en jouant sur HL[3] on voit que ca ameliore les estim
      lev_phi <- lev_phi + notEQL   
    }    
    if (any(lev_lambda<0)) {
      warningList$negLevLam <- TRUE
      lev_lambda[lev_lambda<0] <- 1e-8 ##FR->FR quick patch
    }
    ## compute deviance_residuals 
    ## updated residuals from updated mu must be used (LeeNP p.161) [not so in dhglmfit !!]
    deviance_residual <- family$dev.resids(y,mu,wt=1)
    ######### Dispersion Estimates for phi #####################
    RespLink_disp<-"log"  ## currently not under user control
    if (is.null(phi.Fix)) { ## if phi is estimated (phi.Fix set to 1 for Poisson, Binomial)
      ## leverages have been computed before the  inner loop, which did not change the design matrices 
      Offset_disp <- resid.predictor$offset
      if (models[3]=="phiScal") {
         next_phi_est <- sum(deviance_residual)/sum(1-lev_phi) ## NOT in linkscale
         if (next_phi_est<1e-8) next_phi_est <- 1e-8 # e-10 not high enough (d2hdv2 -> lad -> p_v unstable)
         if (verbose[2]) {print(paste("phi_est=",signif(next_phi_est,4)),quote=F)}
      } else if (models[3]=="phiGLM") { ## there is a phi predictor to estimate but no ranef in this predictor
        resglm_phi <- dispGammaGLM(dev.res=deviance_residual,lev=lev_phi,X=X_disp,offset=Offset_disp)
        if (verbose[2]) {print(paste("phi coefficients=",paste(signif(resglm_phi$coefficients,4),collapse=", ")),quote=F)}
        next_phi_est <- resglm_phi$fitted.values
        lowphi <- which(next_phi_est < 1e-08)
        next_phi_est[lowphi] <- 1e-08 ## to avoid problems with nearly singular matrices
      } else if (models[3]=="phiHGLM") { ## random effect(s) in predictor for phi
        stop("random effects in predictor or residual variance (phi) not yet implemented")
        ## there is a template code with comments in version 260812 of HLfit
        reshglm_phi <- list()
      } 
      if (all(abs(next_phi_est-phi_est) < conv.threshold* (phi_est+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.phi <- TRUE ## 'weak convergence'... 
      } else conv.phi <- FALSE
    } else {conv.phi <- TRUE} ## there is a phi.Fix, already -> phi_est
    ######### Dispersion Estimates for lambda #####################
    if ( is.null(lambda.Fix)) { ## lambda must be estimated     
      ## first builds pseudo response for lambda GLM/HGLM
      z_dimension<-rep(0,nrand)
      for (i in 1:nrand) z_dimension[i]<-qcum[i+1L]-qcum[i]  ## number of levels for each ranef
      psi<-matrix(psi_M,qcum[nrand+1L],1L) ## a column matrix of length the sum of numbers of levels of ranefs...
      resp_lambda<-matrix(0,qcum[nrand+1L],1L)
#########################
    for (i in 1:nrand) {
          u.range <- (qcum[i]+1L):qcum[i+1L]
          pseudoDev <-rand.family$dev.resids(u_h,psi_M,wt=1) ## must give d1 in table p 989 de LeeN01...
          ## si plusieurs ranef, le truc est de calculer trop de valeurs puis de prendre un sous ensemble
          ## a revoir lejour ou +s ranefs; idealement il faudrait un boucle externe sur nrand pour pouvoir faire lamScal sur certains ranef et GLM sur d'autres
          resp_lambda[u.range]<-pseudoDev[u.range]
    }
    if (any(abs(lev_lambda) > 1 - 1e-8)) {
		 lev_lambda<- ifelse(abs(lev_lambda) > 1 - 1e-8, 1 - 1e-8,lev_lambda)
         warningList$leveLam1 <-T
	}
    ## then analyse pseudoresponse
    ## 'lambda' better defined as a vector which size is the number of realizations; 
    ## ** this vector should be updated immediately when the lambda estimates are updated **
    ## lambda_est may be a vector of predicted values, or a scalar, depending how it is computed
    ## The code in logLik.hglm, for example, handles both cases in a minimal way
    if (models[2]=="lamScal") {
         unique.lambda <- sum(resp_lambda)/sum(1-lev_lambda) ## NOT in linkscale ## assumes a single ranef...
         if (unique.lambda<1e-8) unique.lambda <- 1e-8 # as for phi
         if (unique.lambda>1e10) unique.lambda <- 1e10 # 
         if (verbose[2]) {print(paste("lambda=",signif(unique.lambda,4)),quote=F)}
         next_lambda_est <- rep(unique.lambda,nrow(resp_lambda))
      } else if (models[2]=="lamGLM") { ## there is a lambda to estimate but no ranef in its linear predictor
         resglm_lambda <- dispGammaGLM(dev.res=resp_lambda,lev=lev_lambda,X=X_lambda)
         next_lambda_est <- resglm_lambda$fitted.values ## $fitted.value is NOT in linkscale, contrairement a $coefficients
         lowlambda <- which(next_lambda_est < 1e-08)
         next_lambda_est[lowlambda] <- 1e-08 ## to avoid problems with nearly singular matrices
      } else  if (models[2]=="lamHGLM") { ## if ranef in predictor lambda...
           stop("random effects in predictor or ranef variance (lambda) not yet implemented")
           ## there is a template code with comments in version 260812 of HLfit
           reshglm_lambda <- list()
      } 
      if (all(abs(log(next_lambda_est/lambda_est)) < conv.threshold) ) { ## for low values, precision on lambda must be O(v_h^2) ... need precision in relative terms
        conv.lambda <- TRUE ## this is the simplest, best case. ## but if slow geometric decrease of lambda to 0, this is never true 
      } else if (all(abs(next_lambda_est-lambda_est) < conv.threshold* (lambda_est+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.lambda <- TRUE ## 'weak convergence'... 
      } else conv.lambda <- FALSE
    } else { 
       lambda_est <-lambda.Fix
       conv.lambda <- TRUE
    } ## end if null lambda.Fix else ...
    iter<-iter+1 ## here first from 0 to 1
    ## We need to make sure either that convergence of lambda occurred on a relative log scale ( loop not stopping at max.iter !) so that the v_h are very accurate on same scale
    ## or that the v_h's are computed with the very latest lambda, otherwise a call with ranFix$lambda does not yield the same result as estimated lambda
    if ( conv.phi && conv.lambda) {
      ## do not update phi and lambda so that the v_h where computed from the latest lambda_est in particular
      break 
    } else if (iter>=max.iter) { ## did not converge...
      break 
    } else { ## update and continue
      if ( is.null(phi.Fix)) {
        phi_est <- next_phi_est
        w.resid<-as.vector(blob$GLMweights /phi_est) ## 'weinu', must be O(n) in all cases
      }
      if ( is.null(lambda.Fix)) {
        lambda_est<-next_lambda_est
        wranefblob <- updateWranef(rand.family,lambda_est,u_h,v_h)
        w.ranef <- wranefblob$w.ranef
        dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
        dvdu <- wranefblob$dvdu
      }
      if (verbose[2]) {
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
    if (verbose[1]) {
       if (iter==max.iter) {
         mess <- paste("(beta,v)/lambda/phi iterations failed to converge in",max.iter,"iterations")
         mess <- pastefrom(mess,prefix="(!) From ")
         message(mess)
       } else {
          message(paste("(beta,v)/lambda/phi iterations in HLfit() converged in",iter,"iterations"))
       }
    }
    if (family$family %in% c("gaussian","Gamma")) {
       mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-lev_phi))
    }
    if (pforpv>0 && models[1]=="etaHGLM" && max.iter >0) { ## condition on max.iter <=> some params have been fitted
      if (is.null(XinvS)) XinvS <- calc.XinvS(Sig,X.pv,stop.on.error)      
      beta_cov <- try(solve(XinvS%*%X.pv),silent=T) ## solve(small matrix !)
      if (class(beta_cov)=="try-error") {
        beta_se <- rep(Inf,pforpv) ## maybe...
      } else {
        beta_se <- diag(beta_cov)
        if (any(beta_se<0)) { ## divergence of XinvS%*%X.pv leads to negative variance estimates
           beta_se <- rep(Inf,pforpv) ## maybe... 
        } else beta_se <- sqrt(beta_se)
      }
    } else beta_se <- NULL 
    ######################
    ######################
    ######################
    ##### LAMBDA
    if (is.null(lambda.Fix)) {
      if (models[2]=="lamHGLM") { 
       ## there is a template code with comments in version 260812 of HLfit
      } else {
        if (models[2]=="lamScal") {
          ## to compute the se we need the GLM residuals etc. So if the GLM has not been previously used it's better to use it here
          resglm_lambda <- dispGammaGLM(dev.res=resp_lambda,lev=lev_lambda,X=X_lambda)
          res3<-summary(resglm_lambda,dispersion=2) ## 2 dispersion=for the gamma GLM, probably 
        }
        p_lambda<-ncol(X_lambda) ## # of parameters... except if two ranefs...
        linkscale.lambda<-res3$coefficients[1:p_lambda]
        lambda_se<-res3$coefficients[(p_lambda+1L):(2L*p_lambda)] ## ici probleme de transfo log
      }
    }       
    ##### PHI
    if ( is.null(phi.Fix)) {
       if (models[3]=="phiHGLM") {
       ## there is a template code with comments in version 260812 of HLfit
       } else {
#          if (models[3]=="phiScal") {
            ## We need the $coefficients. So if the GLM has not been previously called it's better to call it here
            resglm_phi <- dispGammaGLM(dev.res=deviance_residual,lev=lev_phi,X=X_disp,offset=Offset_disp)
            res2<-summary(resglm_phi)
#          }
          p_disp<-ncol(X_disp) 
          linkscale.phi <- res2$coefficients[1:p_disp]
          phi_se <- res2$coefficients[(p_disp+1L):(2L*p_disp)]
       }
    } 
    ########## LIKELIHOODS
    theta<-theta.mu.canonical(mu/BinomialDen,family$family)  
    Wresid <- diag(w.resid)
    if (models[1]=="etaHGLM" && pforpv==0) d2hdv2 <- - sweep(t(ZL),MARGIN=2,w.resid,`*`) %*% ZL - diag(w.ranef) ##  - t(ZL) %*% diag(w.resid) %*% ZL - diag(w.ranef)
    calcpv <- calc.p_v()
    if (models[1]!="etaHGLM" && models[3]!="phiHGLM") { ## ie GLM, not HGLM
      ml <- calcpv$clik ## vanilla likelihood
      d2hdx2<- -t(X.Re)%*%Wresid%*%X.Re ## X should be the one for leverages
      lad<- determinant(d2hdx2/(2*pi))$modulus[1] ## REML for estim phi
      rl <- ml - lad/2
      cAIC<- -2*ml+2*pforpv
    } else { ## add likelihood of ranef
      if (models[1]=="etaHGLM") {
        clik <- calcpv$clik
        hlik <- calcpv$hlik
        p_v <- calcpv$p_v 
      }
      ## see readable account of aic in HaLM07
      d2hdbv2<-rbind(cbind((t(X.Re)%*%Wresid%*%X.Re),(t(X.Re)%*%Wresid%*%ZL)),
                                    cbind((t(ZL)%*%Wresid%*%X.Re),(-1*d2hdv2)))
      d2clikdbv2<-rbind(cbind((t(X.Re)%*%Wresid%*%X.Re),(t(X.Re)%*%Wresid%*%ZL)),
                             cbind((t(ZL)%*%Wresid%*%X.Re),(t(ZL)%*%Wresid%*%ZL)))
      if (models[3]=="phiHGLM") {
        mess <- pastefrom("correction needed for p_bv for phi DHGLMs.")
        stop(mess)
      } else hv10<-0 ## code cleanup 20/01/13
      if (models[3]=="lamHGLM") {
        mess <- pastefrom("correction needed for p_bv for lambda DHGLMs.")
        stop(mess)
      } else hv20<-0 ## idem
      ##### P_BV
      lad <- determinant(d2hdbv2/(2*pi))$modulus[1]
      #### distinct handling of AIC and p_bv (L-BFGS-B requires a non trivial value):
      if (is.nan(lad) || AIC ) { 
          eigvals <- eigen(d2hdbv2/(2*pi),only.values = T)$values
      }
      if (is.nan(lad)) { 
          zut <- abs(eigvals) ## 05/01/13
          zut[zut<1e-12] <- 1e-12
          lad <- sum(log(zut)) ## 
      }
      p_bv<- hlik-(hv10+hv20+lad/2)  
      if ( ! is.null(calcpv$second.corr)) p_bv <- p_bv + calcpv$second.corr
      if ( AIC ) {
# a debugging issue is that options(error=recover) acts before tryCatch gets the return value
# from its first argument. So a tryCatch on solve is not a good idea.
          if (min(eigvals)>1e-12) {
             qr.d2hdbv2 <- qr(d2hdbv2)
             pd <- sum(diag(safesolve.qr.matrix(qr.d2hdbv2,d2clikdbv2,stop.on.error=stop.on.error)))
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
    }
    if (models[1]!="etaHGLM") {
       APHLs<-list(c.lik=clik)
    } else APHLs<-c(calcpv,list(p_bv=p_bv))
    APHLs$cAIC <- cAIC
    ## Building the return object
    if (models[1]!="etaHGLM") { ## quick patch...
       resglm$APHLs <- APHLs 
       if(verbose[1]) summary(resglm)
       return(resglm) ## NOT of HLfit class...
    } ################################ ELSE
    res<-list(APHLs=APHLs)
    res$data <- data ## very useful for simulate...
    res$family <- family
    res$HL <- HL ## info on fitting method
    res$rand.family <- rand.family
    ##
    res$ranef <- u_h
    res$v_h <- v_h
    if (family$family=="binomial") {
       res$weights <- BinomialDen
    }
    res$y <- y ## counts for Pois/bin
    if ( is.null(phi.Fix) ) {
      res$resid.predictor <- resid.predictor
    } else {
      res$resid.predictor <- NULL
    }
#    res$nu <- 1/phi_est
    res$phi <- phi_est
    res$X <- X.pv
    res$ranFix <- ranFix ## currently as a uniform template consistent with projected changes ; excpt that lamFix, phiFix info is now in lambda.object, etc
    res$corrPars <- ranFix[names(ranFix) %in% c("rho","nu","Nugget")] ## idem. cf functions in corrMM.LRT thta always look in phi, lambda, rather than .Fix. 
    if (family$family=="binomial") {
       res$fv <- mu/BinomialDen ## cf glm(binomial): fitted values are frequencies 
    } else {res$fv <- mu} ## fitted values may be counts (cf poisson), or reals
    res$models <- models
    res$predictor <- predictor ##  all functions must expect predictor, not formula 
    res$ZLMatrix <- ZL ## used by simulate.HL
    res$REMLformula <- REMLformula
    res$fixef <- beta_eta
    res$eta <- eta ## convenient for defining starting values...
    names(res$fixef) <- colnames(X.pv) 
    if (models[2]=="lamScal") { 
      print_lambda <- lambda_est[1]
    } else print_lambda <- lambda_est
    res$lambda <- print_lambda
    if (models[1]=="etaHGLM") {
      res$lev_phi <- lev_phi ## exists even if no phi to be estimated
      res$std_dev_res <- sign(y-mu) * deviance_residual/(phi_est*(1-lev_phi))
    } 
    res$lev_lambda <- lev_lambda
    if (max.iter<1) return(res) ## FR->FR  !! NOT class==HLfit !!
    ## FR->FR the idea was to put below this line everything that requires the computation of a fit as defined by max.iter
    ## now, either only a likelihood is computed, or the firts iteration of the main loop *before* the test on max.iter (late modif of code ?) may provide much of what follows...
    res$RespLink_disp <- RespLink_disp
    if ( adjREML ) res$X.Re <- X.Re
#    if ( ! is.null(init.HLfitName)) {
    if ( ! is.na(spaMM.getOption("INIT.HLFITNAME"))) {
       nextinit.HLfit <- list()
       nextinit.HLfit$fixef <- beta_eta
       nextinit.HLfit$v_h <- v_h
       if (is.null(lambda.Fix)) nextinit.HLfit$lambda <- lambda_est
       spaMM.options(INIT.HLFITNAME=nextinit.HLfit)
       ##assign(init.HLfitName, nextinit.HLfit,pos=".GlobalEnv")
    }  
    res$fixef_se <- beta_se
    if (is.null(lambda.Fix)) {
      res$lambda.object <- list(linkscale.lambda=linkscale.lambda,lambda_se=lambda_se,namesX_lambda=colnames(X_lambda),namesRE=namesRE)
    } else {
      res$lambda.object <- list(lambda.fix=lambda.Fix)
	}
    if (is.null(phi.Fix)) {
       res$phi.object <- list(linkscale.phi=linkscale.phi,phi_se=phi_se,namesX_disp=colnames(X_disp))
    } else {
       res$phi.object <- list(phi.Fix=phi.Fix)
	}
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
    if (maxit.mean>1 && innerj==maxit.mean) {
        warningList$innerNotConv <- paste("linear predictor estimation did not converge. Try increasing 'max.iter.mean' above ",maxit.mean,sep="")
    }
    if (iter==max.iter) {
        warningList$mainNotCov <- paste("Joint estimation did not converge. Try increasing 'max.iter' above ",max.iter,sep="")
    }
    res$warnings <- warningList
    class(res) <- c("HLfit",class(res)) 
    if (verbose[1]) summary(res) 
    if (verbose[1] && length(warningList)>0 ) { ## snif...
       silent<-sapply(length(warningList),function(i) {cat(warningList[[i]]);cat("\n")}) 
    }
    res$call <- mc
    return(res)
}
