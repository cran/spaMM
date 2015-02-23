


auglinmodfit <- function(TT,ZAL,lambda_est,wranefblob,d2hdv2,w.resid,beta_eta,
                         maxit.mean,eta,u_h,v_h,Sig,
                         control.HLfit,
                         X.pv,
                         etaFix,
                         cum_n_u_h,psi_M,
                         muetablob, 
                         family,prior.weights,phi_est,verbose,
                         ranFix, ## only for error message in calc.tXinvS
                         corrPars, ## only for error message in intervalStep
                         processed
) {
  if (inherits(eta,"Matrix")) eta <- as.matrix(eta) ## eh oui...
  ### 
  BinomialDen <- processed$BinomialDen
  conv.threshold <- processed$conv.threshold
  stop.on.error <- processed$stop.on.error
  LevenbergM <- processed$LevenbergM
  betaFirst <- processed$betaFirst ## avoids explicitly solving as an aug.lin.mod.
  LMMbool <- processed$LMMbool
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  y <- processed$y
  off <- attr(processed$predictor,"offsetObj")$vector
  HL <- processed$HL
  mu <- muetablob$mu
  dmudeta <- muetablob$dmudeta
  pforpv <- ncol(X.pv)
  nobs <- NROW(X.pv)
  betaV <- c(beta_eta,v_h)   
  ## t(ZAL) used repeatedly, for constant ZAL, within this function 
  #if ( ! inherits(ZAL,"identityMatrix")) {
    tZAL <- as.matrix(t(ZAL))
  #} else tZAL <- ZAL
  `compute.sscaled` <- function(sqrt.ww,qr.d2hdv2) { ## needs qr.d2hdv2, ZAL, stop.on.error, d2hdv2, rWW, ZALI, family, dmudeta, BinomialDen, mu, eta... ETC
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
    if (processed$canonicalLink) {
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
      Pdiag <- calc.Pdiag(ZAL=ZAL,sqrt.ww=sqrt.ww)
      #
      nrand <- length(lcrandfamfam)
      seqn_u_h <- seq_len(cum_n_u_h[nrand+1L])
      Pdiagn <- Pdiag[1:nobs]
      if (is.na(vecdi1)) vecdi1 <- Pdiagn * coef1
      # K2 is K2 matrix in LeeL appendix p. 4 and is -D in MolasL p. 3307 
      # W is Sigma^-1 ; TWT = t(ZALI)%*%W%*%ZALI = ZAL'.Wresid.ZAL+Wranef = -d2hdv2 !
      K2 <- solveWrap.matrix(qr.d2hdv2,tZAL,stop.on.error=stop.on.error) ## t(ZAL) missing in some implementations... no effect with Gaussian ranefs...  
      if (class(K2)=="try-error") {
        mess <- pastefrom("problem in 'K2' computation.",prefix="(!) From ") ## cf BB 877
        warning(mess)
        K2 <- ginv(d2hdv2) %*% tZAL            
      } ## so far we have computed (d2hdv2)^{-1}.t(Z)= -(TWT)^{-1}.t(Z)
      if (is.na(vecdi2)) {
        # K1 is K1 in LeeL appendix p. 4 and is A=-ZD in MolasL p. 3307-8 
        K1 <- ZAL %id*id% K2
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
  
  
  
  
  
  if (LevenbergM) {
    damping <- 1/1000000 ## as suggested by Madsen-Nielsen-Tingleff... ## mauvais resultats si on part + haut
    dampingfactor <- 2
  }
  levQ <- NULL  
  SQR <- NULL ## only useful for testing not-useeigen code
  w.ranef <- wranefblob$w.ranef ; dlogWran_dv_h <- wranefblob$dlogWran_dv_h ; dvdu <- wranefblob$dvdu
  for (innerj in 1:maxit.mean) {
    ## breaks when conv.threshold is reached
    ##################
    ### Inner Loop ### IF random effect *in mean predictor*: estim beta,v [only v if pforpv=0] for given phi,lambda
    ################## 
    sqrt.w1 <- sqrt(w.resid) ## if maxit.mean>1 GLM weights have been changed and the following must be recomputed
    sqrt.w2 <- sqrt(w.ranef) ##         
    sqrt.ww <- c(sqrt.w1,sqrt.w2) 
    if (inherits(TT,"Matrix")) {
      wAugX <- Diagonal(x=sqrt.ww) %*% TT ## apparemment sweep vriament pas efficace sur Matrix
    } else wAugX <- sweepZ1W(TT,sqrt.ww) # sweep(TT,MARGIN=1,sqrt.ww,`*`) # rWW %*%TT
    old_betaV <- betaV
    ######## According to 'theorem 1' of LeeL12
    ## new beta estimate from z1-a(i)
    ## where z1 is
    z1 <- eta+(y-mu)/dmudeta-off ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
    if (inherits(z1,"Matrix")) z1 <- as.numeric(z1) ## conversion from 1-col dgCMatrix because c(z1 dgCMatrix,z2 numeric) does not work
    ## and a(i) (for HL(i,1)) is a(0) or a(0)+ something
    ## and a(0) depends on z2, as follows :
    if (all(lcrandfamfam=="gaussian")) {
      z2 <- rep(0,length(w.ranef)) 
      a <- rep(0,nobs)
    } else { ## HGLM: nonzero z2, nonzero a(0) ## this could perhaps make a separate block, but nevertheless sometimes qr.d2hdv2 computation...
      nrand <- length( rand.families)
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
      if (is.null(attr(d2hdv2,"qr"))) { 
        a <- try(solve(d2hdv2, - aa),silent=TRUE)
        if (class(a)=="try-error") {
          attr(d2hdv2,"qr") <- QRwrap(d2hdv2)
          a <- solveWrap.vector(attr(d2hdv2,"qr"),  -aa,stop.on.error=stop.on.error)
          ## FR->FR patch: 
          if (class(a)=="try-error") {
            mess <- pastefrom("the Hessian matrix appears singular. Extreme lambda/phi value and/or extremely correlated random effects?",
                              prefix="(!) From ")
            message(mess)
            cat(paste("max(lambda estimates)=",max(lambda_est)))
            if (length(corrPars)>0) {
              cat("; correlation parameters=")
              cat(paste(names(corrPars),"=",corrPars))
            }
            largeLambdaMessages()
            stop()
          }
        }  
      } else { ## we already have a qr, we use it
        a <- solveWrap.vector(attr(d2hdv2,"qr"), -aa,stop.on.error=stop.on.error)
      }    
      a <- Sig %*% ( w.resid * (ZAL %id*id% a) ) ## a(0) in LeeL12
      a <- as.numeric(a) ## incase it became a Matrix, which oddly does not fit with z1-a below...
    }         
    ## and the 'something' for a(1) is computed as follows
    if (HL[1]>0 && (! LMMbool )) { #pforpv>0 && removed since sscaled used for u_h estimation too...
      if (is.null(attr(d2hdv2,"qr"))) {attr(d2hdv2,"qr") <- QRwrap(d2hdv2)} ## Cholwrap tested
      sscaled <- compute.sscaled(sqrt.ww=sqrt.ww,qr.d2hdv2=attr(d2hdv2,"qr")) ## s detadmu
    } else sscaled <- 0L 
    ######## new estimates (tentative if LevenbergMM) 
    if (HL[1]=="SEM") {
      tXinvS <- NULL
      v_h <- solve(-d2hdv2, (tZAL %*% ((z1 - X.pv %*% beta_eta) * w.resid)+ z2*w.ranef))          
      betaV <- c(beta_eta,v_h)
      dbetaV <- betaV - old_betaV
    } else if (betaFirst)  { ### LeeL12 top p. 6 Appendix (code non optimise, useful for checking other algorithms) 
      ## c'est bien equivalent au calcul de Augz en fonction de sscaled essaye ci dessous selon LeeL12
      tXinvS <- calc.tXinvS(Sig,X.pv,stop.on.error,lambda_est,ranFix) ## note dependence v_h -> eta -> Sig...
      ## from a(0) to a(1) LeeL12 p. 963 col 1 l 2
      a <- as.numeric(a + Sig%*% (w.resid * sscaled)) ## in case it became a Matrix...
      rhs <-  tXinvS %*% (z1-a) ## already correct in 1.0
      qr.XtinvSX <- QRwrap(tXinvS%*%X.pv) ## looks contrived but avoids computing sqrt(Sig) (! not diag !); and XinvS%*%X.pv is small
      beta_eta <- solveWrap.vector( qr.XtinvSX , rhs ,stop.on.error=stop.on.error) 
      v_h <- solve(-d2hdv2, (tZAL %*% ((z1 - X.pv %*% beta_eta) * w.resid)+ z2*w.ranef))        
      betaV <- c(beta_eta,v_h)
      dbetaV <- betaV - old_betaV
    } else { ### true augmented model, whether LevenbergM or not;
      tXinvS <- NULL
      if (LMMbool) {
        Augz <- c(z1,z2) ## sscaled=0L (la formule generale s'applique mais perd du temps)
      ## ici il avant code "version 1.0"  
      } else { ## solution of augmented system
        ## (1) what's needed here is the factor of T w on the RHS, not the factor of XinvS in the direct eqns above
        ##    As proven this gives z1-a in one algo and z1- sscaled in the other (!= version 1.0)
        ## (2) the first operation in LevenbergMstep is to substract (eta^0,v^0): LM_wAugz <- wAugz - wAugX %*% betaV
        ##    so we keep (eta^0,v^0) here
        Augz <- c(z1- sscaled,
                  z2+ ((1/w.ranef) * tZAL) %*% (sscaled * w.resid ))  ## 
        ## z2 correction  constructed from eqs (3)(4) of LeeL12 App. p5 and consistent with eq B1 MolasL:
        ## so that row 2 of wAugX.(z1-sscaled,z2+...) = Z'W1 z1 + W2 z2 => Z W1 sscaled - W2 (...) =0 => (...)=
      }
      wAugz <- Augz*sqrt.ww  
      
      if ( maxit.mean > 1 && LevenbergM) { ## default(for *G*LMM) LevenbergM
        LM_wAugz <- wAugz - as.matrix(wAugX %*% betaV)
        ## verif correct rhs: verif_dbetaV <- safesolve.qr.vector(SQR, LM_wAugz,stop.on.error=stop.on.error)            
        if (.spaMM.data$options$USEEIGEN) {
          bla <- LevenbergMstepCallingCpp(wAugX=as.matrix(wAugX),LM_wAugz=LM_wAugz,damping=damping)
        } else bla <- LevenbergMstep(wAugX=wAugX,LM_wAugz=LM_wAugz,damping=damping,stop.on.error=stop.on.error)
        dbetaV <- bla$dbetaV
        betaV <- betaV + dbetaV
        denomGainratio <- bla$denomGainratio
      } else { ## simple aug lin or basic IRLS depending on maxit.mean
        ## QR appears faster than alternatives with crossprod(wAugX); see version 040213
        if (.spaMM.data$options$USEEIGEN) {
          if ( !is.null(control.HLfit$intervalInfo)) {
            currentp_v <- calc.p_v(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,
                                   d2hdv2=d2hdv2,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,
                                   processed=processed,family=family,prior.weights=prior.weights)$p_v
            intervalBlob <- intervalStep(old_betaV=old_betaV,wAugX=wAugX,wAugz=wAugz,
                                         currentp_v=currentp_v,
                                         intervalInfo=control.HLfit$intervalInfo,corrPars=corrPars)
            levQ <- intervalBlob$levQ ; betaV <- intervalBlob$betaV
          } else { ## default for LMMs
            if (inherits(wAugX,"Matrix")) {
              ## DEFAULT COMES LAST
              if (.spaMM.data$options$processedQRmethod == "Matrix::qr") {
                  #xpx <- crossprod(wAugX)
                #xpy <- crossprod(wAugX, wAugz)
                #betaV <- solve(xpx,xpy) 
                #levQ <- (wAugX %*% (chol2inv(chol(xpx)))%*% t(wAugX)) ## suremet améliorable
                ## see http://cran.r-project.org/web/packages/Matrix/vignettes/Comparisons.pdf
                qrwAugX <- Matrix::qr(wAugX)
                betaV <- Matrix::qr.coef(qrwAugX,wAugz)
                levQ <- as.matrix(Matrix::qr.Q(qrwAugX))
              } else if (.spaMM.data$options$processedQRmethod == "lmwithSparseQ") {
                message("lmwithSparseQ called -- should be devel code only") ## protection...
                betaVQ <- lmwithSparseQ(wAugX,wAugz) ## tragically slow,; cf system.time's in commented code below
                betaV <- betaVQ$coef
                pivI <- sort.list(betaVQ$P) ## no pivoting with lmwithQ, pivoting with sparse
                levQ <- as.matrix(betaVQ$Q_ap[,pivI][,seq_len(ncol(wAugX))]) # using Q_ap, not Q
                #print(system.time(Rcpp_sparseQR(wAugX)))
                #print(system.time(Rcpp_QR(as.matrix(wAugX))))
                if (FALSE) { ## confinement de code de debugage 
                  qrwAugX <- Matrix::qr(wAugX)
                  essai <- Matrix::qr.Q(qrwAugX)
                  print(max(abs(diag(tcrossprod(essai))-diag(tcrossprod(levQ)))))
                  essai <- Matrix::qr.coef(qrwAugX,wAugz)
                  print(max(abs(essai-betaV)))
                }
              } else if (.spaMM.data$options$processedQRmethod == "lmwithQ_sparseZAL") {     
                betaVQ <- lmwithQ(as.matrix(wAugX),wAugz) ## slow step
                betaV <- betaVQ$coef
                levQ <- betaVQ$Q[,seq_len(ncol(wAugX))] ## Eigen's HouseholderQ returns a square matrix...
              } else { 
                #stop("Unknown (processed) QRmethod")
                ## lmwithQ_denseZAL but ZAL still Matrix: 
                ## only valid case would be when ZAL was Identity. (in which cas matrix::qr may be useful *if the matrix is large*)
                wAugX <- as.matrix(wAugX) ## because both TT and wAugX kept being Matrix              }
              }
            } ## end all Matrix cases 
            ## matrix case not exclusive to Matrix case because of the latest subcase above 
            if (is.matrix(wAugX)) { 
              ## wAugX is matrix not Matrix (lmwithQ_denseZAL), with useEigen
              betaVQ <- lmwithQ(wAugX,wAugz) 
              betaV <- betaVQ$coef
              levQ <- betaVQ$Q[,seq_len(ncol(wAugX))] ## Eigen's HouseholderQ returns a square matrix...
            }
          }
          #if (any(is.nan(betaV))) browser()
        } else { ## basic IRLS without use eigen (lmwithQ_denseZAL)
          SQR <- QRwrap(wAugX)
          betaV <- solveWrap.vector(SQR, wAugz,stop.on.error=stop.on.error) ## qr.coef(SQR, wAugz) ## vector
        }
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
    u_h <- u_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam) 
    checkv_h <- attr(u_h,"v_h") ## verif is v_h was corrected 
    if (!is.null(checkv_h)) v_h <- checkv_h
    # print(paste(innerj," ",paste(beta_eta,collapse=" ")),quote=F)
    ####### new values of everything, only tentative if LevenbergM  
    eta <- as.matrix(off + X.pv %*% beta_eta + ZAL %id*id% v_h) ## updated at each inner iteration
    ## update functions u_h,v_h
    keep_wranefblob <- wranefblob ####################### 
    wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h)
    w.ranef <- wranefblob$w.ranef
    dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
    dvdu <- wranefblob$dvdu
    keep_muetablob <- muetablob ############################ straight from function argument
    muetablob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen,priorWeights=prior.weights) 
    mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
    dmudeta <- muetablob$dmudeta
    ## update functions of v_h -> blob
    w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    #### update fns of v_h -> blob -> w.resid
    if ( LevenbergM ) { ## for LMM, d2hdv2 is constant // mu, hence it is constant over LevenbergMarquardt iterations
      if (pforpv>0) keep_Sig <- Sig
      keep_d2hdv2 <- d2hdv2
    }
    if (pforpv>0) {
      Sig <- Sigwrapper(ZAL,1/w.ranef,1/w.resid) ## update v_h -> blob$GLMweights -> w.resid -> Sig -> next beta estimate
    }
    d2hdv2 <- calcDhDv2(ZAL,w.resid,w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    ######### for LevenbergM, we check there was a progress, and restore everything otherwise
    if (LevenbergM) {      
      if (HL[1]==0L) {
        currentlik <- calc.p_v(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,
                               d2hdv2=d2hdv2,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,
                               processed=processed,family=family,prior.weights=prior.weights,only.h=TRUE)$hlik ##  use h in PQL/L (> v1.0)  
      } else currentlik <- calc.p_v(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,
                                    d2hdv2=d2hdv2,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,
                                    processed=processed,family=family,prior.weights=prior.weights)$p_v
      if (innerj==1L) { ## not a good idea to use the p_v computed for all u_h initially set to zero
        gainratio <- 1
      } else {
        if (currentlik==-Inf) { ## obs in binary probit with extreme eta... 
          gainratio <- -1
        } else {
          gainratio <- 2*(currentlik-oldlik)/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
        }
      }
      ## amusing debugging code 
      #if (class(try(if (gainratio<0) {}))=="try-error") { ## nnote that NaN *is* numeric!
      #  mess <- pastefrom("'try(if (gainratio<0)...) failed.",prefix="(!) From ") 
      #  stop(mess)
      #}
      if (gainratio<0) { ## restore everything
        currentlik <- oldlik ## important to restore this for the next test !
        betaV <- old_betaV
        if (pforpv>0) {
          beta_eta <- betaV[seq_len(pforpv)]
          names(beta_eta) <- colnames(X.pv)
          v_h <- betaV[-seq_len(pforpv)] 
        } else v_h <- betaV
        u_h <- u_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam) 
        checkv_h <- attr(u_h,"v_h") ## verif is v_h was corrected 
        if (!is.null(checkv_h)) v_h <- checkv_h
        eta <- as.matrix(off + X.pv %*% beta_eta + ZAL %id*id% v_h) ## updated at each inner iteration
        wranefblob <- keep_wranefblob ########################## restoration
        w.ranef <- wranefblob$w.ranef ; dlogWran_dv_h <- wranefblob$dlogWran_dv_h ; dvdu <- wranefblob$dvdu
        muetablob <- keep_muetablob ######################################### restoration
        mu <- muetablob$mu ; dmudeta <- muetablob$dmudeta ; w.resid <-as.vector(muetablob$GLMweights /phi_est) 
        if (pforpv>0) Sig <- keep_Sig ###### more restoration
        d2hdv2 <- keep_d2hdv2 ###### more restoration
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
      relvariation <- dbetaV*(c(rep(1,pforpv),w.ranef)) 
      if (mean(abs(relvariation)) < conv.threshold) break; ## FR->FR mean(abs) is not standard ?  
    }
  } ## end for (innerj in 1:maxit.mean)
  return(list(beta_eta=beta_eta,v_h=v_h,u_h=u_h,eta=eta,wranefblob=wranefblob,
              w.resid=w.resid,Sig=Sig,d2hdv2=d2hdv2,wAugX=wAugX,tXinvS=tXinvS,
              sqrt.ww=sqrt.ww,innerj=innerj,levQ=levQ,SQR=SQR,muetablob=muetablob)
  )
} ## end auglinmodfit
