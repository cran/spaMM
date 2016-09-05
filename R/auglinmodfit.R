####### new values of everything, only local to this function
eval_gain_LevM_HGLM <- function(betaV, conv_dbetaV, pforpv, X.pv, etaFix,
                                u_h_v_h_from_v_h_arglist,
                                updateW_ranefS_arglist, wranefblob, ## need one of these two
                                maxit.mean, processed, lambda_est, off, phi_est, ddi_or_matrix_ZAL, 
                                calc_gainratio, LevenbergMstep_result, oldlik,provide.qr) {
  betaV <- betaV + conv_dbetaV ## betaV may be corrected below
  beta_eta <- betaV[seq_len(pforpv)] ## provide it for pforpv=0 too
  if (pforpv>0) {
    if (is.null(etaFix$v_h)) v_h <- betaV[-seq_len(pforpv)] 
  } else {if (is.null(etaFix$v_h)) v_h <- betaV}
  u_h <- do.call("u_h_v_h_from_v_h",c(u_h_v_h_from_v_h_arglist,list(v_h=v_h)))
  if (maxit.mean > 1L && !is.null(attr(u_h,"v_h"))) { ## second test = if upper.v_h or lower.v_h non NULL
    v_h <- attr(u_h,"v_h")
    if (is.null(etaFix$v_h)) betaV[(pforpv+1L):length(betaV)] <- v_h ## to be copied in old_betaV, in valid space
    ## but conv_dbetaV unchanged for assessing convergence
  }
  ## update functions u_h,v_h
  if ( ! processed$GLMMbool) { 
    wranefblob <- do.call("updateW_ranefS",c(updateW_ranefS_arglist,list(u_h=u_h,v_h=v_h)))
  } ## else keep input wranefblob since lambda_est not changed
  eta <- off + drop(X.pv %*% beta_eta) + drop(ddi_or_matrix_ZAL %*% v_h) 
  muetablob <- muetafn(eta=eta,BinomialDen=processed$BinomialDen,processed=processed) 
  mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  ## update functions of v_h -> blob
  w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
  d2hdv2 <- calcD2hDv2(ddi_or_matrix_ZAL,w.resid,wranefblob$w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
  locarglist <- list(mu=mu,u_h=u_h,dvdu=wranefblob$dvdu,lambda_est=lambda_est,phi_est=phi_est,
                     d2hdv2=d2hdv2, processed=processed,provide.qr=provide.qr)
  if (processed$HL[1]==0L) {
    newAPHLs <- do.call(calc.p_v,c(locarglist,list(only.h=TRUE)))
    currentlik <- unlist(newAPHLs["hlik"]) ##  use h in PQL/L (> v1.0)  ## keep name
  } else {
    newAPHLs <- do.call("calc.p_v",locarglist)
    currentlik <- unlist(newAPHLs["p_v"]) ## keep name
  }

  if (calc_gainratio) {
    if (currentlik==-Inf) { ## obs in binary probit with extreme eta... 
      gainratio <- -1
    } else {
      summand <- conv_dbetaV*(LevenbergMstep_result$rhs+ LevenbergMstep_result$dampDpD * conv_dbetaV) 
      ## In the summand, all terms should be positive. conv_dbetaV*rhs should be positive. However, numerical error may lead to <0 or even -Inf
      ## Further, if there are both -Inf and +Inf elements the sum is NaN and HLfit fails.
      summand[summand<0] <- 0
      denomGainratio <- sum(summand)
      gainratio <- 2*(currentlik-oldlik)/denomGainratio ## cf computation in MadsenNT04 below 3.14, but generalized for D' != I ; numerator is reversed for maximization
    }
  } else {gainratio <- 1}
  attr(d2hdv2,"APHLs") <- newAPHLs
  return(list(gainratio=gainratio,currentlik=currentlik,betaV=betaV,beta_eta=beta_eta,
              v_h=v_h,u_h=u_h,eta=eta,wranefblob=wranefblob,muetablob=muetablob,w.resid=w.resid,
              d2hdv2=d2hdv2))
}


simple_WLS_with_Eigen <- function(wAugX,wAugz,QRmethod) { ## in auglinmodfit(): no need for LevenbergM when LMM 
  qr_wAugX <- betaVQ <- betaV <- NULL
  if (inherits(wAugX,"Matrix")) {
    ## DEFAULT COMES LAST
    if (QRmethod == "Matrix::qr") {
      ## see http://cran.r-project.org/web/packages/Matrix/vignettes/Comparisons.pdf
      qr_wAugX <- attr(wAugX,"get_qr")(wAugX)
      betaV <- Matrix::qr.coef(qr_wAugX,wAugz)
    } else { 
      wAugX <- as.matrix(wAugX) ## so that the next block for 'm'atrix could be executed
      betaV <- .lm.fit(x=wAugX,y=wAugz)$coefficients ## SEE BELOW
    }
  } ## end all Matrix cases 
  ## matrix case not exclusive to Matrix case because of a possible intervening wAugX <- as.matrix(wAugX)
  if (inherits(wAugX,"matrix")) { 
    # if ( ! is.null(get_qr <- attr(wAugX,"get_qr")) {
    #   qr_wAugX <- get_qr(wAugX) 
    #   betaV <- solveWrap.vector(qr_wAugX,wAugz)
    # }
    if (FALSE) {
      ###### fastLmPure
      ## 0 for the column-pivoted QR decomposition, 
      ## 1 for the unpivoted QR decomposition, 
      ## 2 for the LLT Cholesky, 3 for the LDLT Cholesky, ...................
      ## bechmarks: http://dirk.eddelbuettel.com/blog/2011/07/05/
      ##            http://stackoverflow.com/questions/30420185/fastlm-is-much-slower-than-lm
      ## In my experience (denser matrices ?) .lm.fit remains faster
      # betaV <- RcppEigen::fastLmPure(X=wAugX,y=wAugz,method=1)$coefficients
      ######
    } else betaV <- .lm.fit(x=wAugX,y=wAugz)$coefficients
  }
  return(list(betaV=betaV,betaVQ=betaVQ))
}

calc_sscaled <- function(vecdisneeded, dlogWran_dv_h, coef12, Pdiag, n_u_h, nobs, K2, ZAL) {
  if  (any(vecdisneeded[-3L])) {
    coef12 <- coef12 ## eval promise
    coef1 <- coef12$coef1 # coef1 is the factor of P_ii in d1
    coef2 <- coef12$dlW_deta # coef2 is the factor between P_jj and K1 in d2
  }
  vecdis <- vecdisneeded
  vecdis[vecdisneeded] <- NA
  vecdis[!vecdisneeded] <- 0
  vecdi1 <- vecdis[1L]
  vecdi2 <- vecdis[2L]
  vecdi3 <- vecdis[3L]
  if (any(vecdisneeded)) { ## buut its stupid to call the function otherwise
    ## here version 1.5.3 has an interesting signed.wAugX concept
    ## P is P in LeeL appendix p. 4 and is P_R in MolasL p. 3307 
    Pdiag <- Pdiag ## eval promise
    Pdiagn <- Pdiag[1:nobs]
    seqn_u_h <- seq_len(n_u_h)
    if (vecdisneeded[1L]) vecdi1 <- Pdiagn * coef1
    # K2 is K2 matrix in LeeL appendix p. 4 and is -D in MolasL p. 3307 
    # W is Sigma^-1 ; TWT = t(ZALI)%*%W%*%ZALI = ZAL'.Wresid.ZAL+Wranef = -d2hdv2 !
    if (vecdisneeded[2L]) { # ( ZAL %*% K2 ) is K1 in LeeL appendix p. 4 and is A=-ZD in MolasL p. 3307-8 
      vecdi2 <- drop( ((Pdiagn * coef2) %*id% ZAL) %*% K2)
    }
    # coef3 =(1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
    if (vecdisneeded[3L]) {  ## d3 reste nul pour gaussian ranef
      vecdi3 <- drop( (Pdiag[nobs+seqn_u_h] * dlogWran_dv_h[seqn_u_h]) %*% K2)
    }
    vecdi <- vecdi1+vecdi2+vecdi3 ## k_i in MolasL; le d_i de LeeL app. p. 4
    sscaled <- vecdi /2  ## sscaled := detadmu s_i= detadmu d*dmudeta/2 =d/2 in LeeL12; or dz1 = detadmu (y*-y) = detadmu m_i=0.5 k_i dmudeta = 0.5 k_i in MolasL 
  } else sscaled <- 0
  return(sscaled)
}

Sig_times_b <- function(Sig0,ZAL,w.ranef,w.resid,b) { # Sig= [Sig0=Z.(1/w.ranef).Z^t+1/w.resid]
  if (is.null(Sig0)) { ## w.ranef is variable
    v1 <- ZAL %*% (t(t(b) %*% ZAL)/w.ranef)
  } else {
    v1 <- Sig0 %*% b
  }
  v2 <- b/w.resid
  return(as.numeric(v1+v2))
}

auglinmodfit <- function(TT, ## used to update wAugX
                         ZAL, ## constant within the function
                         tZAL, ## idem
                         lambda_est,wranefblob,d2hdv2,w.resid,beta_eta,
                         maxit.mean,eta,u_h,v_h,
                         control.HLfit,
                         X.pv,off,
                         etaFix, ## FR->FR still uses $v_h, put perhaps reconsider
                         cum_n_u_h,psi_M,
                         muetablob, 
                         phi_est, ## for GLM weights; constant within...
                         verbose, ## [["trace"]] used toward the end
                         ranFix, corr_est, ## both only for error message in calc_tXinvS and message in intervalStep
                         processed,
                         ZALtZAL=NULL,
                         ddi_or_matrix_ZAL ## for d2hdv2 and Sig and calc_Pdiag
) {
  BinomialDen <- processed$BinomialDen
  conv.threshold <- processed$conv.threshold
  stop.on.error <- processed$stop.on.error
  LevenbergM <- processed$LevenbergM
  betaFirst <- processed$betaFirst ## avoids explicitly solving as an aug.lin.mod.
  LMMbool <- processed$LMMbool
  family <- processed$family
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  y <- processed$y
  off <- attr(processed$predictor,"offsetObj")$total 
  HL <- processed$HL
  mu <- muetablob$mu
  dmudeta <- muetablob$dmudeta
  pforpv <- ncol(X.pv)
  nobs <- NROW(X.pv)
  betaV <- c(beta_eta,v_h)   
  nrand <- length(lcrandfamfam)
  u_h_info <- processed$u_h_info 
  lower.v_h <- u_h_info$lower.v_h
  upper.v_h <- u_h_info$upper.v_h
  boxConstraintsBool <- u_h_info$boxConstraintsBool
  #
  if (LevenbergM) {
    damping <- 1e-7 ## as suggested by Madsen-Nielsen-Tingleff... ## mauvais resultats si on part + haut
    dampingfactor <- 2
  }
  w.ranef <- wranefblob$w.ranef 
  dlogWran_dv_h <- wranefblob$dlogWran_dv_h 
  dvdu <- wranefblob$dvdu
  Sig0 <- NULL
  updates <- TRUE
  # Tried to implement conditional LevM, but this is slower ! 
  NewtonNext <- FALSE ## FALSE in first iteration at least
  if (HL[1]==0L) { likfn <- "hlik" } else {likfn <- "p_v"} ##  use h in PQL/L (> v1.5.59) 
  for (innerj in 1:maxit.mean) {
    if (updates) { ## new IRLS response and hessian; TRUE except if last step (LevM or basic IRLS) failed
      if (HL[1]>0 && (! LMMbool )) { #pforpv>0 && removed since sscaled used for u_h estimation too...
        ########## HL(1,.) adjustment for mean ################## and specifically the a(1) term in LeeL 12 p. 963
        ## if LMM (ie resp gaussian, ranef gaussian), all coef<x> are 0 -> correction is 0 (but this fn must not have been called)
        ## if (gaussian, not gaussian) d3 nonzero
        ## if (non gaussian, gaussian), d3 zero (!maybe not for all possible cases) but d1,d2 nonzero 
        coef12needed <- ! ((family$family=="gaussian" && family$link=="identity")
                           || (family$family=="Gamma" && family$link=="log")  ) ## two ad hoc cases
        vecdisneeded <- c( coef12needed, coef12needed, any(dlogWran_dv_h!=0L) )
        K2needed <- ( HL[1]!=0L && any(vecdisneeded[-1L]) )
      } else K2needed <- FALSE
      sqrt.w1 <- sqrt(w.resid) ## if maxit.mean>1 GLM weights have been changed and the following must be recomputed
      ## here version 1.5.3 has an interesting signed.wAugX concept
      if (anyNA(sqrt.w1)) { 
        stop("NA/NAN in sqrt.w1; some w.ranef<0 ?")
      } else {
        sqrt.w2 <- sqrt(w.ranef) ##         
        sqrt.ww <- c(sqrt.w1,sqrt.w2) 
      }  
      wAugX <- calc_wAugX(augX=TT,sqrt.ww=sqrt.ww)
      old_betaV <- betaV
      ######## According to 'theorem 1' of LeeL12, new beta estimate from z1-a(i), where z1 is
      z1 <- eta+(y-mu)/dmudeta-off ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
      if (inherits(z1,"Matrix")) z1 <- as.numeric(z1) ## conversion from 1-col dgCMatrix because c(z1 dgCMatrix,z2 numeric) does not work
      ## and a(i) (for HL(i,1)) is a(0) or a(0)+ something
      ## and a(0) depends on z2, as follows :
      if (all(lcrandfamfam=="gaussian")) {
        z2 <- rep(0,length(w.ranef)) 
        a <- rep(0,nobs)
      } else { ## HGLM: nonzero z2, nonzero a(0) ## this could perhaps make a separate block, but nevertheless sometimes qr.d2hdv2 computation...
        psi_corr <- unlist(lapply(seq(nrand), function(it) {
          u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
          if (lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") { 
            return(2*u_h[u.range]- (u_h[u.range]^2)*(1+lambda_est[u.range])) ## LeeL01 p.1003; to cast the analysis into the form of z2  
          } else if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { ## gamma(identity)
            return(2*u_h[u.range] - (u_h[u.range]^2)/(1-lambda_est[u.range])) ## interesting singularity 
            ## moreover pb: u_h=1, lambda =1/2 -> psi=0 -> z2=0 -> negative u_h
          } else {   
            return(psi_M[u.range])  
          } 
        }))
        z2 <- v_h + (psi_corr-u_h)*dvdu ## update since u_h,v_h updated (yes)
        #        nXn  .   nXn      nX'r'    'r'X'r'       'r'X'r'    'r'
        # a <- Sig %*% Wresid %*% ZAL %*% solve(-d2hdv2) %*% Wranef %*% z2 ## p. 963 l. 1-2; a(0) supp mat p. 6 
        aa <- w.ranef * z2
        # qr_d2hdv2 <- attr(d2hdv2,"get_qr")(d2hdv2,provide=K2needed)
        qr_d2hdv2 <- attr(d2hdv2,"get_qr")(d2hdv2,provide=TRUE) ## remarquably, TRUE seems faster...
        if (is.null(qr_d2hdv2)) { 
          a <- try(solve(d2hdv2, - aa),silent=TRUE)
          if (inherits(a,"try-error")) {
            qr_d2hdv2 <- attr(d2hdv2,"get_qr")(d2hdv2)
            a <- solveWrap.vector(qr_d2hdv2,  -aa,stop.on.error=stop.on.error)
            ## FR->FR patch: 
            if (inherits(a,"try-error")) {
              mess <- pastefrom("the Hessian matrix appears singular. Extreme lambda/phi value and/or extremely correlated random effects?",
                                prefix="(!) From ")
              message(mess)
              cat(paste("max(lambda estimates)=",max(lambda_est)))
              dispCorrPars <- get_DispCorrPars(ranFix, corr_est)
              if (length(dispCorrPars)>0) {
                cat("; dispersion and correlation parameters=")
                cat(paste(names(dispCorrPars),"=",dispCorrPars))
              }
              largeLambdaMessages()
              stop()
            }
          }  
        } else { ## we already have a qr, we use it
          a <- solveWrap.vector(qr_d2hdv2, -aa,stop.on.error=stop.on.error)
        }    
        a <- Sig_times_b(Sig0=NULL, ZAL=ZAL, w.ranef=w.ranef,w.resid=w.resid,b= w.resid * (ZAL %id*% a) )
        # a <- Sig %*% ( w.resid * (ZAL %id*% a) ) ## a(0) in LeeL12
        # a <- as.numeric(a) ## incase it became a _M_atrix, which oddly does not fit with z1-a below...
      }         
      ## and the 'something' for a(1) is computed as follows
      if (HL[1]>0 && (! LMMbool )) { #pforpv>0 && removed since sscaled used for u_h estimation too...
        ########## HL(1,.) adjustment for mean ################## and specifically the a(1) term in LeeL 12 p. 963
        ## if LMM (ie resp gaussian, ranef gaussian), all coef<x> are 0 -> correction is 0 (but this fn must not have been called)
        ## if (gaussian, not gaussian) d3 nonzero
        ## if (non gaussian, gaussian), d3 zero (!maybe not for all possible cases) but d1,d2 nonzero 
        if (K2needed) { ## do this outside calc_sscaled() to manage qr.d2hdv2 globally
          qr_d2hdv2 <- attr(d2hdv2,"get_qr")(d2hdv2)
          K2 <- solveWrap.matrix(qr_d2hdv2,tZAL,stop.on.error=stop.on.error) ## t(ZAL) missing elsewhere; for non-Gaussian ranefs.
          if (inherits(K2,"try-error")) {
            mess <- pastefrom("problem in 'K2' computation.",prefix="(!) From ") ## cf BB 877
            warning(mess)
            K2 <- ginv(d2hdv2) %*% tZAL            
          } ## so far we have computed (d2hdv2)^{-1}.t(Z)= -(TWT)^{-1}.t(Z)
        } else K2 <- NULL
        sscaled <- calc_sscaled(vecdisneeded=vecdisneeded,
                     dlogWran_dv_h=dlogWran_dv_h, ## This (dlogWran_dv_h) was computed when w.ranef was computed
                     coef12= calc.dlW_deta(dmudeta=dmudeta, family=family, 
                                           mu=mu, eta=eta, BinomialDen=BinomialDen, 
                                           canonicalLink=processed$canonicalLink,
                                           calcCoef1=TRUE), ## evaluated if any vecdisneeded[-3]
                     ## under ML, Pdiag are the leverages given the current c(w.resid,w.ranef)
                     Pdiag=calc.Pdiag(ddi_or_matrix_ZAL, c(w.resid,w.ranef),wAugZALI=wAugX[,-seq_len(pforpv)]), 
                     n_u_h=cum_n_u_h[nrand+1L], nobs=nobs, K2= K2,
                     ZAL=ZAL # vecdi2
        )
      } else sscaled <- 0L 
      ######## new estimates (tentative if LevenbergMM) 
      ## auglinmodfit not called with SEM...
      #     if (HL[1]=="SEM") { 
      #       tXinvS <- NULL
      #       v_h <- solve(-d2hdv2, (tZAL %*% ((z1 - X.pv %*% beta_eta) * w.resid)+ z2*w.ranef))          
      #       betaV <- c(beta_eta,v_h)
      #       conv_dbetaV <- betaV - old_betaV
      #     } else 
      if (betaFirst)  { ### LeeL12 top p. 6 Appendix (code non optimise, useful for checking other algorithms) 
        ## c'est bien equivalent au calcul de Augz en fonction de sscaled essaye ci dessous selon LeeL12
        Sig <- Sigwrapper(ddi_or_matrix_ZAL,1/w.ranef,1/w.resid,ZALtZAL=ZALtZAL)
        tXinvS <- calc_tXinvS(Sig,X.pv,stop.on.error) ## note dependence v_h -> eta -> Sig...
        if (inherits(tXinvS,"try-error")) {
          dispCorrPars <- get_DispCorrPars(ranFix,corr_est)
          singularSigmaMessagesStop(lambda_est=lambda_est,phi_est=phi_est,corrPars=dispCorrPars)
        }
        ## from a(0) to a(1) LeeL12 p. 963 col 1 l 2
        a <- a + Sig_times_b(Sig0=NULL, ZAL=ZAL, w.ranef=w.ranef, w.resid=w.resid,  b= (w.resid * sscaled) )
        # a <- as.numeric(a + Sig%*% (w.resid * sscaled)) ## in case it became a Matrix...
        rhs <-  tXinvS %*% (z1-a) ## already correct in 1.0
        qr.XtinvSX <- QRwrap(tXinvS%*%X.pv) ## looks contrived but avoids computing sqrt(Sig) (! not diag !); and XinvS%*%X.pv is small
        beta_eta <- solveWrap.vector( qr.XtinvSX , rhs ,stop.on.error=stop.on.error) 
        v_h <- solve(-d2hdv2, ( (tZAL %*% ((z1 - X.pv %*% beta_eta) * w.resid))+ z2*w.ranef))        
        betaV <- c(beta_eta,v_h)
        conv_dbetaV <- betaV - old_betaV
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
      } 
    } ## end if (updates) => new IRLS response and hessian
    
    if ( maxit.mean > 1L && ## FR->FR useful ? prevents testing LevenbergM on LMM
         LevenbergM && ! NewtonNext) { ## default(for *G*LMM) LevenbergM
      LM_wAugz <- wAugz - drop(wAugX %*% betaV)
      ## verif correct rhs: verif_dbetaV <- safesolve.qr.vector(qr_wAugX, LM_wAugz,stop.on.error=stop.on.error)
      ## verif correct v_h probleme specifique v ~ (+/-) Gamma
      #
      ## here version 1.5.3 has an interesting signed.wAugX concept
      if (.spaMM.data$options$USEEIGEN) {
        LevenbergMstep_result <- LevenbergMstepCallingCpp(wAugX=as.matrix(wAugX),
                                                          LM_wAugz=LM_wAugz,damping=damping)
      } else LevenbergMstep_result <- LevenbergMstep(wAugX=wAugX,LM_wAugz=LM_wAugz,damping=damping,stop.on.error=stop.on.error)
      conv_dbetaV <- LevenbergMstep_result$dbetaV## the one that will be used for assessing convergence, 
      levMblob <- eval_gain_LevM_HGLM(betaV, conv_dbetaV, pforpv, X.pv, etaFix,
                                      u_h_v_h_from_v_h_arglist=list(rand.families=rand.families,cum_n_u_h=cum_n_u_h,
                                                                    lcrandfamfam=lcrandfamfam,lower.v_h=lower.v_h,upper.v_h=upper.v_h),
                                      updateW_ranefS_arglist=list(cum_n_u_h=cum_n_u_h,rand.families=rand.families,lambda=lambda_est),
                                      wranefblob=wranefblob,
                                      maxit.mean, processed, lambda_est, off, phi_est, ddi_or_matrix_ZAL,
                                      calc_gainratio = (innerj>1L), LevenbergMstep_result, oldlik,provide.qr=K2needed)
      if (levMblob$gainratio<0) { ## unsuccesful step: do not update anything 
        damping <- damping*dampingfactor
        dampingfactor <- dampingfactor*2 
        updates <- FALSE
      } else { ## update everything (always TRUE when innerj=1)
        updates <- TRUE
        ## NewtonNext <- TRUE ## slow
        oldlik <- currentlik <- levMblob$currentlik
        betaV <- levMblob$betaV
        beta_eta <- levMblob$beta_eta
        v_h <- levMblob$v_h
        u_h <- levMblob$u_h
        eta <- levMblob$eta 
        ### here any updates except those of the IRLS eq themsleves (but eg for lkelihood computation
        if (!processed$GLMMbool) {
          wranefblob <- levMblob$wranefblob
          w.ranef <- wranefblob$w.ranef
          dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
          dvdu <- wranefblob$dvdu
        } ## else nothing changed since lambda_est not changed
        muetablob <- levMblob$muetablob 
        mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
        dmudeta <- muetablob$dmudeta
        ## update functions of v_h -> blob
        w.resid <- levMblob$w.resid
        d2hdv2 <- levMblob$d2hdv2
        ## cf Madsen-Nielsen-Tingleff again, and as in levmar library by Lourakis
        damping <- damping * max(1/3,1-(2*levMblob$gainratio-1)^3)  
        dampingfactor <- 2
      }
    } else { ## simple aug lin (default if LMMbool), or basic IRLS (still if LevenbergM but innerj=1)
      # Would be possible if LevenbergM IF NewtonNext COULD BE TRUE 
      # Also for *Intervals*  
      # 'updates' remains always TRUE in this case 
      ## QR appears faster than alternatives with crossprod(wAugX); see version 040213
      if (.spaMM.data$options$USEEIGEN) {
        if ( !is.null(control.HLfit$intervalInfo)) {
          calcLikArglist <- list(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,
                                 d2hdv2=d2hdv2, processed=processed)
          currentlik <- do.call(calc.p_v,calcLikArglist)[likfn]
          currentlik <- unlist(currentlik) ## a way of keeping the name of the lik
          #print(c(lambda_est[1],betaV[1:4],currentlik))  ## to locate convergence issue 
          intervalBlob <- intervalStep(old_betaV=old_betaV,wAugX=wAugX,wAugz=wAugz,
                                       currentlik=currentlik,
                                       intervalInfo=control.HLfit$intervalInfo,
                                       ranFix=ranFix,corr_est=corr_est, ## both only for message()
                                       likfn=likfn, QRmethod=processed$QRmethod)
          betaV <- intervalBlob$betaV
        } else { ## (default if LMMbool), or basic IRLS (still if LevenbergM but innerj=1)
          WLS_blob <- simple_WLS_with_Eigen(wAugX,wAugz,QRmethod=processed$QRmethod)
          betaV <- WLS_blob$betaV
          betaVQ <- WLS_blob$betaVQ
        }
      } else { ## basic IRLS without use eigen (lmwithQ_denseZAL)
        qr_wAugX <- attr(wAugX,"get_qr")(wAugX)
        betaV <- solveWrap.vector(qr_wAugX, wAugz,stop.on.error=stop.on.error) ## qr.coef(qr_wAugX, wAugz) ## vector
      }
      if (inherits(betaV,"try-error")) betaV <- ginv(wAugX)%*% wAugz ## occurred with large lambda either as 'init.HLfit', or by the iterative algo
      ## updating
      betaV <- drop(betaV) 
      if (maxit.mean > 1L) conv_dbetaV <- betaV - old_betaV
      if (pforpv>0) {
        beta_eta <- betaV[seq_len(pforpv)]
        names(beta_eta) <- colnames(X.pv)
        if (is.null(etaFix$v_h)) v_h <- betaV[-seq_len(pforpv)] 
      } else {if (is.null(etaFix$v_h)) v_h <- betaV}
      u_h <- u_h_v_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,lower.v_h,upper.v_h) 
      #checkv_h <- attr(u_h,"v_h") ## verif is v_h was corrected 
      #if (!is.null(checkv_h)) v_h <- checkv_h
      if (maxit.mean > 1L && !is.null(attr(u_h,"v_h"))) { ## second test = if upper.v_h or lower.v_h non NULL
        oldv_h <- v_h
        v_h <- attr(u_h,"v_h")
        if (is.null(etaFix$v_h)) betaV[(pforpv+1L):length(betaV)] <- v_h ## to be copied in old_betaV, in valid space
        ## but conv_dbetaV unchanged for assessing convergence
      }
      eta <- off + drop(X.pv %*% beta_eta) + drop(ZAL %id*% v_h) 
      ### here any updates except those of the IRLS eq themsleves (but eg for lkelihood computation
      if (!processed$GLMMbool) {
        wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h)
        w.ranef <- wranefblob$w.ranef
        dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
        dvdu <- wranefblob$dvdu
      } ## else nothing changed since lambda_est not changed
      muetablob <- muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed) 
      mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
      dmudeta <- muetablob$dmudeta
      ## update functions of v_h -> blob
      w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
      #### update fns of v_h -> blob -> w.resid
      d2hdv2 <- calcD2hDv2(ddi_or_matrix_ZAL,w.resid,w.ranef) ## update d2hdv2= - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    } ## end (if LevenbergM) ELSE IRLS...
    if ( LevenbergM && NewtonNext) { ## If conditional LevM was allowed, this would operate
      newLikArglist <- list(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,
                            d2hdv2=d2hdv2,processed=processed)
      newlik <- do.call(calc.p_v,newLikArglist)[[likfn]]
      ## The first oldlik was computed by previous LevM step that set NewtonNext to TRUE 
      if (newlik< oldlik) { # failure of Newton step => try again with LevM
        updates <- FALSE 
        NewtonNext <- FALSE
      } else {
        updates <- TRUE ## it must already be, but this anticipates future problems
        ## NewtonNext remains true
        oldlik <- newlik
      } 
    }
    ########## nearly done with one inner iteration
    if (verbose[["trace"]]) {
      cat(paste("Inner iteration ",innerj,sep=""))
      if (LevenbergM) cat(paste(", ",names(currentlik),"= ",currentlik,sep=""))
      cat("\n")
      if (innerj>1) cat("norm.dbetaV=",sqrt(sum(conv_dbetaV^2)))
      cat(paste("; beta_eta=",paste(beta_eta,collapse=", ")))
      cat("\n")
      print("================================================")
    } 
    if (maxit.mean>1) {
      ## the convergence on v_h^2 must be relative to lambda; this raises questions about the lowest meaningful lambda values.
      relvariation <- conv_dbetaV*(c(rep(1,pforpv),w.ranef)) 
      if (mean(abs(relvariation)) < conv.threshold) break; ## FR->FR mean(abs) is not standard ?  
    }
    
  } ## end for (innerj in 1:maxit.mean)
  if (.spaMM.data$options$USEEIGEN) {
    if ( !is.null(control.HLfit$intervalInfo)) {
      levQ <- intervalBlob$levQ ## maybe not optimized, but gets the Q of the *local wAugX* used by intervalStep
    } else levQ <- NULL # not useful otherwise 
  } ## else basic IRLS without use eigen 
  return(list(beta_eta=beta_eta,v_h=v_h,u_h=u_h,eta=eta,wranefblob=wranefblob,
              w.resid=w.resid,d2hdv2=d2hdv2,wAugX=wAugX,tXinvS=tXinvS,
              sqrt.ww=sqrt.ww,innerj=innerj,levQ=levQ,muetablob=muetablob)
  )
} ## end auglinmodfit
