
## more convenient public interface with a generic and methods
"get_from_MME" <- function(sXaug,which="",szAug=NULL,B=NULL,...) UseMethod("get_from_MME") 

## pure solve, not returning decomp
"get_from_MME_default" <- function(sXaug,which="",szAug=NULL,B=NULL,...) UseMethod("get_from_MME_default") 

# get_from -> sparseMatrix and default methods

get_from_MME.default <- function(sXaug,which="",szAug=NULL,B=NULL,...) {
  method <- attr(sXaug,"get_from")
  if (length(method)==0L) {
    method <- "'sXaug' has no 'get_from' attribute."
    ## : useful for trace(get_from.sparseMatrix,exit=quote(print(method)))
    get_from_MME_default.matrix(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  } else {
    ## using match.call() is terribly slow! => passing ... without match.call
    #  warning("method without ad hoc code in get_from.sparseMatrix")
    do.call(what=method,
            args=c(list(sXaug=sXaug,which=which,szAug=szAug,B=B),list(...)))
  }
}

get_from_MME.sparseMatrix <- function(sXaug,which="",szAug=NULL,B=NULL,...) {
  method <- attr(sXaug,"get_from")
  if (length(method)==0L) {
    method <- "'sXaug' has no 'get_from' attribute."
    ## : useful for trace(get_from.sparseMatrix,exit=quote(print(method)))
    get_from_MME_default.Matrix(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  # } else if ( method=="sXaug_Matrix_QRP_scaled") {
  #   get_from.sXaug_Matrix_QRP_scaled(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  # } else if ( method=="sXaug_Matrix_cholP_scaled") {
  #   get_from.sXaug_Matrix_cholP_scaled(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  # } else if ( method=="sXaug_EigenSparse_QR_scaled") {
  #   get_from.sXaug_EigenSparse_QR_scaled(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  # } else if ( method=="sXaug_EigenSparse_QRP_scaled") {
  #   get_from.sXaug_EigenSparse_QRP_scaled(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  # } else if ( method=="sXaug_EigenSparse_LDLP_scaled") {
  #   get_from.sXaug_EigenSparse_LDLP_scaled(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  } else {
    ## using match.call() is terribly slow! => passing ... without match.call
  #  warning("method without ad hoc code in get_from.sparseMatrix")
    do.call(what=method,
            args=c(list(sXaug=sXaug,which=which,szAug=szAug,B=B),list(...)))
  }
  ## direct calls of the function may be faster but require ad hoc programming...
}

## pure solve, not returning decomp
get_from_MME_default.Matrix <- function(sXaug,which="",szAug=NULL,B=NULL,...) {
  if (which=="" && ! is.null(szAug)) {
    #   ## see http://cran.r-project.org/web/packages/Matrix/vignettes/Comparisons.pdf
    #   return(Matrix::qr.coef(Matrix::qr(sXaug),szAug)) ## Vscaled.beta
    ## FR->FR needs a fully automated selection of methods:
    if (length(grep("QR", .spaMM.data$options$Matrix_method))) {
      solve_method <- ".lmwith_sparse_QRp"
    } else if (length(grep("LDL", .spaMM.data$options$Matrix_method))) {
      solve_method <- ".lmwith_sparse_LDLp"
    } else solve_method <- ".lmwith_sparse_LLp"
    sol <- do.call(solve_method,list(XX=sXaug,yy=szAug,returntQ=FALSE,returnR=FALSE,pivot=TRUE))
    return(sol$coef)
  } else stop("Unhandled arguments in get_from_MME_default.Matrix (missing method for get_from_MME() ?)")
}

## pure solve, not returning decomp
get_from_MME_default.matrix <- function(sXaug,which="",szAug=NULL,B=NULL,...) {
  if (which=="" && ! is.null(szAug)) {
    if (FALSE) {
      ###### fastLmPure
      ## 0 for the column-pivoted QR decomposition, 
      ## 1 for the unpivoted QR decomposition, 
      ## 2 for the LLT Cholesky, 3 for the LDLT Cholesky, ...................
      ## benchmarks: http://dirk.eddelbuettel.com/blog/2011/07/05/
      ##            http://stackoverflow.com/questions/30420185/fastlm-is-much-slower-than-lm
      ## In my experience (denser matrices ?) .lm.fit remains faster
      # betaV <- RcppEigen::fastLmPure(X=sXaug,y=szAug,method=1)$coefficients
      ######
    #} else return(.lm.fit(x=sXaug,y=szAug)$coefficients) ## 
    } else return(.lmwithQR(sXaug,szAug,returntQ=FALSE,returnR=FALSE)$coef) 
  } else stop("Unhandled arguments in get_from_default.matrix")
}

.calc_sXaug_Re <- function(locsXaug, ## conforming template
                          X.Re,weight_X) {
  distinct.X.ReML <- attr(X.Re,"distinct.X.ReML")
  n_u_h <- attr(locsXaug,"n_u_h")
  if ( distinct.X.ReML[1L] ) {
    locsXaug <- locsXaug[,-(n_u_h+attr(X.Re,"unrestricting_cols"))]
  } 
  extra_vars <- attr(X.Re,"extra_vars") ## may be NULL
  if (inherits(locsXaug,"Matrix")) {
    suppl_cols <- Matrix(0,ncol=length(extra_vars),nrow=nrow(locsXaug))
    suppl_cols[n_u_h+seq(nrow(X.Re)),] <- Diagonal(x=weight_X) %*% X.Re[,extra_vars]
  } else {
    suppl_cols <- matrix(0,ncol=length(extra_vars),nrow=nrow(locsXaug))
    suppl_cols[n_u_h+seq(nrow(X.Re)),] <- diag(x=weight_X) %*% X.Re[,extra_vars]
  }
  locsXaug <- cbind(locsXaug,suppl_cols)
  return(locsXaug)
}


# function to get the hatvalues (only: not the other similar computations on t_Q_scaled)
# no permutation issues for Q => a single get_hatvalues function should handle all sXaug classes
.get_hatvalues <- function(sXaug, X.Re, weight_X) {
  if ( ! is.null(X.Re) ) { ## not the standard REML
    distinct.X.ReML <- attr(X.Re,"distinct.X.ReML")
    if ( distinct.X.ReML[2L] ) { ## test FALSE for standard ML, TRUE only for some non-standard REML cases
      locsXaug <- .calc_sXaug_Re(locsXaug=sXaug,X.Re,weight_X)
      hatval <- .leveragesWrap(locsXaug) ## Rcpp version of computation through computation of Q
    } else if ( distinct.X.ReML[1L] ) { ## includes ML standard
      whichcols <- attr(X.Re,"unrestricting_cols")
      if (length(whichcols)==attr(sXaug,"pforpv")) { ## should be ML standard
        hatval <- get_from_MME(sXaug,which="hatval_Z")
      } else { ## non-standard case
        t_Q_scaled <- get_from_MME(sXaug,which="t_Q_scaled")
        n_u_h <- attr(sXaug,"n_u_h")
        ## substract cols directly from Q ! -- FR->FR working on t_Q for non-pivoted case only !
        t_Q_scaled <- t_Q_scaled[-(n_u_h+whichcols),] ## test TRUE for standard ML 
        ## [, -integer(0)] would empty the matrix...
        hatval <- colSums(t_Q_scaled*t_Q_scaled)
      }
    } else { ## REML standard
      stop("Ideally this case is not reached") ## REML: ideally X.Re is null. But...  
    }
  } else hatval <- get_from_MME(sXaug,which="hatval") # colSums(t_Q_scaled*t_Q_scaled) ## basic REML, leverages from the same matrix used for estimation of betaV 
  if (is.list(hatval)) hatval <- unlist(hatval) ## assuming order lev_lambda,lev_phi
  return(hatval)
}

.leveragesWrap <- function(X) {
  if (inherits(X,"sparseMatrix")) {
    return(rowSums(qr.Q(qr(X))^2)) ## perhaps not optimal and certainly not optimized as first met in LM with sparse X.pv 01/2018
  } else .leverages(X) 
}


.calc_dvdloglamMat_new <- function(dlogfthdth,cum_n_u_h,lcrandfamfam,rand.families,u_h,
                                   sXaug, d2hdv2=NULL, ## use either one
                                   stop.on.error) {
  neg.d2f_dv_dloglam <- vector("list",length(lcrandfamfam))
  for (it in seq_len(length(lcrandfamfam)) ) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    ## First the cases where g(u) differs from theta(u) : cf oklink dans preprocess pour detection des cas
    ## same computation as canonical case, except that first we consider dlogfthdv=dlogfthdth * [dth/dv]
    if (lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") { 
      neg.d2f_dv_dloglam[[it]] <- (dlogfthdth[u.range] / u_h[u.range])  ## [dth/dv=1/u] for th(u)=-1/u, v=log(u)
    } else if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { 
      neg.d2f_dv_dloglam[[it]] <- (dlogfthdth[u.range] / u_h[u.range]) ## gamma(identity)  ## [dth/dv=1/u] for th(u)=log(u), v=u
    } else { ## v=g(u) = th(u) : random effect model is canonical conjugate
      neg.d2f_dv_dloglam[[it]] <- (dlogfthdth[u.range]) ## (neg => -) (-)(psi_M-u)/lambda^2    *    lambda.... 
    } 
  }
  neg.d2f_dv_dloglam <- unlist(neg.d2f_dv_dloglam)
  neg.d2f_dv_dloglam <- as.vector(neg.d2f_dv_dloglam)
  if(is.null(d2hdv2)) {
    dvdloglamMat <- get_from_MME(sXaug,"solve_d2hdv2",B=diag( neg.d2f_dv_dloglam)) 
  } else {
    qr_d2hdv2 <- .get_qr(d2hdv2,provide=FALSE)
    if (is.null(qr_d2hdv2)) {
      dvdloglamMat <- try(solve(d2hdv2,diag( neg.d2f_dv_dloglam ))) 
    } else dvdloglamMat <- .solveWrap_matrix(qr_d2hdv2, diag( neg.d2f_dv_dloglam ), stop.on.error=stop.on.error) # rXr
    if (inherits(dvdloglamMat,"try-error")) {
      warning("Using ginv() in dvdloglamMat computation to overcome numerical problem.")
      dvdloglamMat <- sweep(ginv(d2hdv2),MARGIN=2,as.vector(neg.d2f_dv_dloglam),`*`) ## ginv(d2hdv2) %*% diag( as.vector(neg.d2f_dv_dloglam))
    }
  } 
  return(dvdloglamMat)
}

.calc_dvdlogphiMat_new <- function(dh0deta,ZAL,
                                   sXaug,d2hdv2=NULL, ## either one
                                   stop.on.error) {
  ## cf calcul dhdv, but here we want to keep each d/d phi_i distinct hence not sum over observations i 
  ## code corrected here 12/2013; this is dh0dv = neg.d2h0_dv_dlogphi (eta always linear in v :-) and w.resid always propto 1/phi)
  neg.d2h0_dv_dlogphi <- .m_Matrix_times_Dvec(t(ZAL), drop(dh0deta)) # dh0dv <- t(ZAL) %*% diag(as.vector(dh0deta)) ## nXr each ith column is a vector of derivatives wrt v_k
  if (is.null(d2hdv2)) {
    dvdlogphiMat <- get_from_MME(sXaug,"solve_d2hdv2",B=neg.d2h0_dv_dlogphi) 
  } else {
    qr_d2hdv2 <- .get_qr(d2hdv2)
    dvdlogphiMat <- .solveWrap_matrix(qr_d2hdv2, neg.d2h0_dv_dlogphi , stop.on.error=stop.on.error)  # rXn       
    if (inherits(dvdlogphiMat,"try-error")) stop("problem in dvdlogphiMat computation.") ## ou warning + ginv  comme dans .calc_dvdloglamMat
  }
  return(dvdlogphiMat)
}


.calc_sscaled_new <- function(vecdisneeded, dlogWran_dv_h, coef12, 
                              n_u_h, nobs, sXaug, ZAL,WU_WT) {
  if  (any(vecdisneeded[-3L])) {
    #coef12 <- coef12 ## eval promise
    coef1 <- coef12$coef1 # coef1 is the factor of P_ii in d1
  }
  vecdis <- vecdisneeded
  vecdis[vecdisneeded] <- NA
  vecdis[!vecdisneeded] <- 0
  vecdi1 <- vecdis[1L]
  vecdi2 <- vecdis[2L]
  vecdi3 <- vecdis[3L]
  if (any(vecdisneeded)) { ## but the call to .calc_sscaled_new is conditional to the same condition 
    ## here version 1.5.3 has an interesting signed.wAugX concept
    ## P is P in LeeL appendix p. 4 and is P_R in MolasL p. 3307 
    Pdiag <- get_from_MME(sXaug,"hatval_Z") 
    seqn_u_h <- seq_len(n_u_h)
    if (vecdisneeded[1L]) vecdi1 <- Pdiag$lev_phi * coef1
    # K2 = solve(d2hdv2,tZAL) is K2 matrix in LeeL appendix p. 4 and is -D in MolasL p. 3307 
    # W is Sigma^-1 ; TWT = t(ZALI)%*%W%*%ZALI = ZAL'.Wresid.ZAL+Wranef = -d2hdv2 !
    if (vecdisneeded[2L]) { # ( ZAL %*% K2 ) is K1 in LeeL appendix p. 4 and is A=-ZD in MolasL p. 3307-8 
      # vecdi2 <- as.vector( ((Pdiag$lev_phi * coef2) %*id% ZAL) %*% K2)
      coef2 <- coef12$dlW_deta # coef2 is the factor between P_jj and K1 in d2
      vecdi2 <- get_from_MME(sXaug,"solve_d2hdv2",B=as.vector((Pdiag$lev_phi * coef2) %*id% ZAL))
      if ( ! is.null(WU_WT)) vecdi2 <- vecdi2/WU_WT # zero-truncated model
      vecdi2 <- as.vector(ZAL %*% vecdi2)
    }
    # coef3 =(1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
    if (vecdisneeded[3L]) {  ## d3 reste nul pour gaussian ranef
      # vecdi3 <- as.vector( (Pdiag$lev_lambda * dlogWran_dv_h[seqn_u_h]) %*% K2)
       vecdi3 <- get_from_MME(sXaug,"solve_d2hdv2",B=as.vector(Pdiag$lev_lambda * dlogWran_dv_h[seqn_u_h]))
       if ( ! is.null(WU_WT)) vecdi3 <- vecdi3/WU_WT # zero-truncated model
       vecdi3 <- as.vector(ZAL %*% vecdi3)
    }
    vecdi <- vecdi1+vecdi2+vecdi3 ## k_i in MolasL; le d_i de LeeL app. p. 4
    sscaled <- vecdi /2  ## sscaled := detadmu s_i= detadmu d*dmudeta/2 =d/2 in LeeL12; or dz1 = detadmu (y*-y) = detadmu m_i=0.5 k_i dmudeta = 0.5 k_i in MolasL 
  } else sscaled <- 0
  return(sscaled)
}

.init_resp_z_corrections_new <- function(lcrandfamfam, w.ranef, nobs, nrand, cum_n_u_h, rand.families, u_h, lambda_est, psi_M, v_h, dvdu, sXaug, stop.on.error, ranFix, ZAL, w.resid) {
  if (all(lcrandfamfam=="gaussian")) {
    z2 <- rep(0,length(w.ranef)) 
    a <- rep(0,nobs)
  } else { ## HGLM: nonzero z2, nonzero a(0) 
    psi_corr <- vector("list",nrand)
    for (it in seq_len(nrand)) {
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      if (lcrandfamfam[it]=="inverse.gamma" && rand.families[[it]]$link=="log") { 
        psi_corr[[it]] <- (2*u_h[u.range]- (u_h[u.range]^2)*(1+lambda_est[u.range])) ## LeeL01 p.1003; to cast the analysis into the form of z2  
      } else if (lcrandfamfam[it]=="gamma" && rand.families[[it]]$link=="identity") { ## gamma(identity)
        psi_corr[[it]] <- (2*u_h[u.range] - (u_h[u.range]^2)/(1-lambda_est[u.range])) ## interesting singularity 
        ## moreover pb: u_h=1, lambda =1/2 -> psi=0 -> z2=0 -> negative u_h
      } else {   
        psi_corr[[it]] <- (psi_M[u.range])  
      } 
    }
    psi_corr <- unlist(psi_corr)
    # w.ranef v^0 + dlogfv_dv ('dlogfvdv' elsewhere) is represented as w.ranef (z2:= v_h + (psi_corr-u_h)*dvdu) 
    #    as detailed in 'Adjustments of the score equations for different random effect ($v$) distributions'
    z2 <- v_h + (psi_corr-u_h)*dvdu ## update since u_h,v_h updated (yes)
    #        nXn  .   nXn      nX'r'    'r'X'r'       'r'X'r'    'r'
    # a <- Sig %*% Wresid %*% ZAL %*% solve(-d2hdv2) %*% Wranef %*% z2 ## p. 963 l. 1-2; a(0) supp mat p. 6 
    aa <- w.ranef * z2
    a <- - get_from_MME(sXaug,"solve_d2hdv2",B= aa )
    a <- .Sig_times_b(Sig0=NULL, ZAL=ZAL, w.ranef=w.ranef,w.resid=w.resid,b= w.resid * (ZAL %id*% a) )
    # a <- Sig %*% ( w.resid * (ZAL %id*% a) ) ## a(0) in LeeL12
  }         
  return(list(a0=a,z20=z2))
}

.Sig_times_b <- function(Sig0,ZAL,w.ranef,w.resid,b) { # Sig= [Sig0=Z.(1/w.ranef).Z^t+1/w.resid]
  if (is.null(Sig0)) { ## w.ranef is variable
    v1 <- drop(t(b) %*% ZAL)
    v1 <- ZAL %*% (v1/w.ranef)
  } else {
    v1 <- Sig0 %*% b
  }
  v2 <- b/w.resid
  return(as.numeric(v1+v2))
}

.calc_zAug_not_LMM <- function(n_u_h, nobs, pforpv, y, off, ZAL, 
                      # variable within fit_as_ZX:
                      eta, muetablob, dlogWran_dv_h, sXaug, w.resid, w.ranef, 
                      init_z_args, 
                      #
                      processed) {
  GLMMbool <- processed$GLMMbool
  coef12needed <- processed$coef12needed
  
  mu <- muetablob$mu
  dmudeta <- muetablob$dmudeta
  ######## According to 'theorem 1' of LeeL12, new beta estimate from z1-a(i), where z1 is
  if (is.list(w.resid)) {
    z1 <- as.vector(eta+w.resid$WU_WT*(y-mu-w.resid$dlogMthdth)/dmudeta-off) ## MolasL10
  } else z1 <- as.vector(eta+(y-mu)/dmudeta-off) ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
  ## and a(i) (for HL(i,1)) is a(0) or a(0)+ something
  ## and a(0) depends on z2, as follows :
  if ( ! GLMMbool) {
    z2 <- do.call(".init_resp_z_corrections_new",init_z_args)$z20 
  } else z2 <- rep(0,n_u_h)
  if (processed$HL[1L]) { 
    ########## HL(1,.) adjustment for mean ################## and specifically the a(1) term in LeeL 12 p. 963
    ## if LMM (ie resp gaussian, ranef gaussian), all coef<x> are 0
    ## if (gaussian, not gaussian) d3 nonzero
    ## if (non gaussian, gaussian), d3 zero (!maybe not for all possible cases) but d1,d2 nonzero 
    vecdisneeded <- c( coef12needed, coef12needed, any(dlogWran_dv_h!=0L) )
    if (any(vecdisneeded)) {
      if (is.list(w.resid)) {
        WU_WT <- w.resid$WU_WT 
      } else WU_WT <- NULL
      sscaled <- .calc_sscaled_new(vecdisneeded=vecdisneeded,
                              dlogWran_dv_h=dlogWran_dv_h, ## dlogWran_dv_h was computed when w.ranef was computed
                              coef12= .calc_dlW_deta(dmudeta=drop(dmudeta), mu=drop(mu), eta=drop(eta), 
                                                    family=processed$family, 
                                                    BinomialDen=processed$BinomialDen, 
                                                    canonicalLink=processed$canonicalLink,
                                                    calcCoef1=TRUE,
                                                    w.resid=w.resid), ## promise evaluated if any vecdisneeded[-3]
                              n_u_h=n_u_h, nobs=nobs, 
                              sXaug=sXaug,
                              ZAL=ZAL, # vecdi2
                              WU_WT=WU_WT ## NULL except for truncated model
      )
    } else sscaled <- 0
    ## sscaled if an 
    if (is.list(w.resid)) {
      sscaled <- sscaled * w.resid$WU_WT
      y2_sscaled <- z2+ as.vector((sscaled * w.resid$w_resid ) %*% ZAL )/w.ranef ## that's the y_2 in "Methods of solution based on the augmented matrix"
    } else y2_sscaled <- z2+ as.vector((sscaled * w.resid ) %*% ZAL )/w.ranef
    zInfo <- list(sscaled=sscaled, z1=z1, z2=z2, z1_sscaled=z1-sscaled, y2_sscaled=y2_sscaled)
  } else zInfo <- list(sscaled=0, z1=z1, z2=z2, z1_sscaled=z1, y2_sscaled=z2) 
  return(zInfo)
}



.warn_intervalStep <- function(currentlik,for_intervals) {
  locmess <- paste("A higher",for_intervals$likfn,"was found than for the original fit.",
                   "\nThis suggests the original fit did not fully maximize",for_intervals$likfn,"\n or numerical accuracy issues.")
  message(locmess)
  dispCorrPars <- .get_CorrEst_and_RanFix(for_intervals$ranFix,for_intervals$corr_est)$corrPars
  if (length(dispCorrPars)>0) message(paste("Current dispersion and correlation parameters are ",
                                            paste(names(dispCorrPars),"=",signif(unlist(dispCorrPars),6),collapse=", ")))
  message("Current log likelihood is =",currentlik)                    
  message("logLik of the fit=",for_intervals$fitlik)    
}


.make_Xscal <- function(ZAL, ZAL_scaling=NULL, AUGI0_ZX, n_u_h, use_Xaug=FALSE) {
  if (use_Xaug) { ## Write and read AUGI0_ZX$Xaug; but 
    # (1) the alternative with the bind's remains faster, despite being a slow step of a fit;
    # this requires reinitializing AUGI0_ZX$Xaug when ZAL is modified, which exposes to bugs.
    if (is.null(AUGI0_ZX$Xaug)) {
      AUGI0_ZX$Xaug <- suppressMessages(cbind2(
        suppressMessages(rbind2(AUGI0_ZX$I, ZAL)), ## suppress signature message generated by eg rbind2(Diagonal(x=runif(5)), Diagonal(n=5)) 
        rbind2(AUGI0_ZX$ZeroBlock, AUGI0_ZX$X.pv)
      )) 
    }
    if (is.null(ZAL_scaling)) {
      return(AUGI0_ZX$Xaug)
    } else return(.m_Matrix_llblock_times_Dvec(AUGI0_ZX$Xaug, ZAL_scaling))
  } else { 
    if (!is.null(ZAL_scaling)) ZAL <- .m_Matrix_times_Dvec(ZAL,ZAL_scaling)
    Xscal <- suppressMessages(cbind2(
      suppressMessages(rbind2(AUGI0_ZX$I, ZAL)), ## suppress signature message generated by eg rbind2(Diagonal(x=runif(5)), Diagonal(n=5)) 
      rbind2(AUGI0_ZX$ZeroBlock, AUGI0_ZX$X.pv)
    )) 
    return(Xscal)
  }
}

## y=u_h in all cases
## for gamma ranef y = u_h and theta = -1 the function reduces to 
## -nu*y+nu*(log(nu*y))-lgamma(nu)-log(y) as it should, LeeNP p. 180
## for beta ranef y = u_h and theta = 1/2 this is also OK
## for inv gamma cf Log[PDF[InverseGammaDistribution[1 + \[Nu], \[Nu]], uh]] + theta heuristically added to fit p. 181...
## To merge this with .get_clik_fn, relationship between theta and psi_M sould be clarified...
.loglfn_ranU <- function(RandDist,y,nu) { ## functions with standardized mean and only a dispersion param
  switch(RandDist,
         gaussian = {- ((y^2)*nu+log(2*pi/nu))/2}, 
         gamma = {-nu*y+nu*(log(nu*y))-lgamma(nu)-log(y)}, ## p. 180 with psi=1 gives log pdf ranV assuming V=logU
         beta = {(nu/2-1)*log(y*(1-y))-lbeta(nu/2,nu/2)}, ## version explained p. 181 LeeNP
         ## Log[PDF[InverseGammaDistribution[1 + \[Nu], \[Nu] \[Mu]], uh]] with Mu=1 + |du/dv|
         "inverse.gamma" = {-nu/y - (2+nu)* log(y) + (1+nu)*log(nu) - lgamma(1+nu)} ## p. 181 with psi=1 gives log pdf ranV assuming V=-1/U, not log pdf ranU
  )
}

.calc_APHLs_from_ZX <- function(auglinmodblob=NULL,processed, which="p_v",
                               ## alternative to auglinmodblob, insuff pour REML non standard:
                               sXaug, phi_est, lambda_est, dvdu, u_h, mu 
                               ) {
  resu <- list()
  if (FALSE  &&  ## in particular pb with p_v
      processed$LMMbool && ! all(which=="hlik") && ! is.null(auglinmodblob)) {
    n_u_h <- length(auglinmodblob$u_h)
    weight_Xaug <- c(rep(1,n_u_h),auglinmodblob$weight_X)
    augy <- c(rep(0,n_u_h),processed$y) ## contains the offset, cntrary to augz in fit_as_ZX !
    SSE <- sum((weight_Xaug*(augy-auglinmodblob$fitted))^2) 
    ## SSE [sum of nobs+nr terms]/nobs provides an estimate of a scaling factor 
    ## not of phi (which could be  sum((y-fitted)[ypos])^2)/sum(1-lev_phi)
    nobs <- auglinmodblob$nobs
    if ("p_v" %in% which) { ## I N V A L I D  formula in REML fits
      resu$p_v <- sum(log(auglinmodblob$weight_X)) - get_from_MME(auglinmodblob$sXaug,"logdet_R_scaled_v") - nobs*(1+log(2*pi*SSE/nobs))/2
    }
    if ("p_bv" %in% which) {
      resdf <- nobs - auglinmodblob$pforpv
      resu$p_bv <- sum(log(auglinmodblob$weight_X)) - get_from_MME(auglinmodblob$sXaug,"logdet_R_scaled_b_v") - resdf*(1+log(2*pi*SSE/resdf))/2
    }
  } else { ## general code
    if ( ! is.null(auglinmodblob)) {
      sXaug <- auglinmodblob$sXaug 
      mu <- auglinmodblob$muetablob$mu
      phi_est <- auglinmodblob$phi_est
      u_h <- auglinmodblob$u_h
      lambda_est <- auglinmodblob$lambda_est
      dvdu <- auglinmodblob$wranefblob$dvdu
      pforpv <- auglinmodblob$pforpv
    }
    #
    family <- processed$family
    famfam <- family$family
    clik_fn <- processed$clik_fn
    y <- processed$y
    BinomialDen <- processed$BinomialDen
    theta <- .theta.mu.canonical(mu/BinomialDen,family)  
    if (famfam=="binomial") {
      resu$clik <- sum(clik_fn(theta,y/BinomialDen,BinomialDen,1/(phi_est))) ## freq a revoir
    } else {
      phi_est[phi_est<1e-12] <- 1e-10 ## 2014/09/04 local correction, has to be finer than any test for convergence 
      ## creates upper bias on clik but should be more than compensated by the lad
      ## correcting the lad makes an overall upper bias for small (y-theta) at "constant" corrected phi 
      ## this can be compensated by correcting the lad LESS.
      resu$clik <- sum(clik_fn(theta,y,eval(processed$prior.weights)/phi_est)) ## note (prior) weights meaningful only for gauss/ Gamma 
    }
    if (processed$models[["eta"]]=="etaGLM") {
      resu$p_v <- resu$clik
      return(resu)
    } # E L S E 
    cum_n_u_h <- processed$cum_n_u_h
    lcrandfamfam <-  attr(processed$rand.families,"lcrandfamfam")
    likranU <- vector("list",length(lcrandfamfam))
    for (it in seq_len(length(lcrandfamfam))) {
      u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
      likranU[[it]] <- .loglfn_ranU(lcrandfamfam[it],u_h[u.range],1/lambda_est[u.range])
    }
    likranU <- unlist(likranU)
    log.du_dv <- - log(dvdu) 
    likranV <- sum(likranU + log.du_dv)
    resu$hlik <- resu$clik + likranV
    #
    n_u_h <- length(lambda_est)
    # beware that computation of logdet_sqrt_d2hdv2 depends on w.ranef
    if ("p_v" %in% which || "p_bv" %in% which) {
      resu$p_v <- resu$hlik - get_from_MME(sXaug,"logdet_sqrt_d2hdv2") + n_u_h*log(2*pi)/2
    }
    if ("p_bv" %in% which) {
      X.Re <- processed$X.Re
      if ( is.null(X.Re)) {## REML standard
        # beware that computation of logdet_r22 depends on H_global_scale
        resu$p_bv <- resu$p_v - get_from_MME(sXaug,"logdet_r22") + pforpv*log(2*pi)/2
        #browser()
      } else if ( ncol(X.Re)==0L) {## ML standard
        resu$p_bv <- resu$p_v
      } else {
        #resu$p_bv <- NA ## FR->FR non-standard REML not yet handled
        locXscal <- auglinmodblob$sXaug   
        weight_X <- auglinmodblob$weight_X
        nobs <- auglinmodblob$nobs
        H_global_scale <- attr(auglinmodblob$sXaug,"H_global_scale")
        w.ranef <- attr(auglinmodblob$sXaug,"w.ranef")
        if (inherits(locXscal,"Matrix")) {
          locXscal <- .Dvec_times_Matrix_lower_block(1/weight_X,locXscal,n_u_h)
          mMatrix_method <- .spaMM.data$options$Matrix_method
        } else {
          Xrows <- n_u_h+seq(nobs)
          locXscal[Xrows,] <- diag(x=1/weight_X) %*% locXscal[Xrows,] ## get back to unweighted scaled matrix
          mMatrix_method <- .spaMM.data$options$matrix_method
        }
        locXscal <- .calc_sXaug_Re(locXscal,X.Re,rep(1,nobs))      ## or some cbind ?  
        locsXaug <- do.call(mMatrix_method,
                         list(Xaug=locXscal, weight_X=weight_X, w.ranef=w.ranef, H_global_scale=H_global_scale))
        
        resu$p_bv <- resu$p_v - get_from_MME(locsXaug,"logdet_r22") + ncol(X.Re)*log(2*pi)/2
      }
    }
  }
  # if (! any(which=="hlik")) resu$hlik <- NA
  return(resu)
}

.calc_APHLs_from_auglinmodblob <- function(auglinmodblob,processed, which, phi_est, lambda_est) {
  APHLs_args <- list(processed=processed, which=which, phi_est=phi_est, lambda_est=lambda_est)
  APHLs_args$sXaug <- auglinmodblob$sXaug
  APHLs_args$dvdu <- auglinmodblob$wranefblob$dvdu
  APHLs_args$u_h <- auglinmodblob$u_h 
  APHLs_args$mu <- auglinmodblob$muetablob$mu
  do.call(".calc_APHLs_from_ZX", APHLs_args)[[which]]
}




