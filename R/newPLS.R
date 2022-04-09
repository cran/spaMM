
## more convenient public interface with a generic and methods
"get_from_MME" <- function(sXaug,which="",szAug=NULL,B=NULL,...) UseMethod("get_from_MME") 

## pure solve, not returning decomp
"get_from_MME_default" <- function(sXaug,which="",szAug=NULL,B=NULL,...) UseMethod("get_from_MME_default") 

# get_from -> sparseMatrix and default methods

get_from_MME.default <- function(sXaug,which="",szAug=NULL,B=NULL,...) {
  method <- attr(sXaug,"get_from")
  if (length(method)==0L) {
    method <- "'sXaug' has no 'get_from' attribute." # local copy useful for tracing with exit=quote(print(method))
    get_from_MME_default.matrix(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  } else {
    ## using match.call() is terribly slow! => passing ... without match.call
    do.call(what=method,
            args=c(list(sXaug=sXaug,which=which,szAug=szAug,B=B),list(...)))
  }
}

get_from_MME.sparseMatrix <- function(sXaug,which="",szAug=NULL,B=NULL,...) {
  method <- attr(sXaug,"get_from")
  if (is.null(method)) {
    method <- "'sXaug' has no 'get_from' attribute."  # local copy useful for tracing with exit=quote(print(method))
    get_from_MME_default.Matrix(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  # } else if ( method=="sXaug_Matrix_QRP_scaled") {
  #   get_from.sXaug_Matrix_QRP_scaled(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  # } else if ( method=="sXaug_EigenSparse_QR_scaled") {
  #   get_from.sXaug_EigenSparse_QR_scaled(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  # } else if ( method=="sXaug_EigenSparse_QRP_scaled") {
  #   get_from.sXaug_EigenSparse_QRP_scaled(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  # } else if ( method=="sXaug_EigenSparse_LDLP_scaled") {
  #   get_from.sXaug_EigenSparse_LDLP_scaled(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
  } else {
    ## using match.call() is terribly slow! => passing ... without match.call
    locfn <- get(method,asNamespace("spaMM"), inherits=FALSE)
    locfn(sXaug=sXaug,which=which,szAug=szAug,B=B,...)
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
      sol <- .lmwith_sparse_QRp(XX=sXaug,yy=szAug,returntQ=FALSE,returnR=FALSE) # no pivot argument
    } else if (length(grep("LDL", .spaMM.data$options$Matrix_method))) {
      sol <- .lmwith_sparse_LDLp(XX=sXaug,yy=szAug,returntQ=FALSE,returnR=FALSE,pivot=TRUE)
    } else sol <- .lmwith_sparse_LLp(XX=sXaug,yy=szAug,returntQ=FALSE,returnR=FALSE,pivot=TRUE)
    return(sol$coef)
  } else stop("Unhandled arguments in get_from_MME_default.Matrix (missing method for get_from_MME() ?)")
}

## pure solve, not returning decomp
get_from_MME_default.matrix <- function(sXaug,which="",szAug=NULL,B=NULL,...) {
  if (which=="" && ! is.null(szAug)) {
    if (FALSE) {
      if (FALSE) {
        ###### fastLmPure
        ## 0 for the column-pivoted QR decomposition, 
        ## 1 for the unpivoted QR decomposition, 
        ## 2 for the LLT Cholesky, 3 for the LDLT Cholesky, ...................
        ## benchmarks: http://dirk.eddelbuettel.com/blog/2011/07/05/
        ##            http://stackoverflow.com/questions/30420185/fastlm-is-much-slower-than-lm
        ## In my experience (denser matrices ?) .lm.fit remains faster
        # betaV <- RcppEigen::fastLmPure(X=sXaug,y=szAug,method=1)$coefficients
        # return(betaV)
        ######
      } else return(.lm.fit(x=sXaug,y=szAug)$coefficients) ## 
    } else { 
      ## according to https://eigen.tuxfamily.org/dox/group__DenseDecompositionBenchmark.html
      ## HouseholderQR is faster thanColPivHouseholderQR. There are precision trade-offs.
      return(.lmwithQR(sXaug,szAug,returntQ=FALSE,returnR=FALSE)$coef)
    } ## Eigen QR is OK bc we don't request Q
  } else stop("Unhandled arguments in get_from_MME_default.matrix")
}

.calc_sXaug_Re <- function(locsXaug, ## conforming template
                          X.Re,weight_X) {
  distinct.X.ReML <- attr(X.Re,"distinct.X.ReML") # This fn should be called only for non-stdd REML
  n_u_h <- attr(locsXaug,"n_u_h")
  if ( distinct.X.ReML[1L] ) {
    locsXaug <- locsXaug[,-(n_u_h+attr(X.Re,"unrestricting_cols"))]
  } 
  extra_vars <- attr(X.Re,"extra_vars") ## may be NULL, in which case the following block works with dense matrices but not sparse ones...
  if (nc <- length(extra_vars)) {
    if (inherits(locsXaug,"Matrix")) {
      suppl_cols <- Matrix(0,ncol=nc,nrow=nrow(locsXaug))
    } else {
      suppl_cols <- matrix(0,ncol=nc,nrow=nrow(locsXaug))
    }
    suppl_cols[n_u_h+seq(nrow(X.Re)),] <- .Dvec_times_m_Matrix(weight_X,X.Re[,extra_vars, drop=FALSE])#  Diagonal(x=weight_X) %*% X.Re[,extra_vars]
    locsXaug <- cbind(locsXaug,suppl_cols)
  }
  return(locsXaug)
}

.calc_sXaug_Re_spprec <- function(locsXaug, ## conforming template
                                  X.Re) {
  distinct.X.ReML <- attr(X.Re,"distinct.X.ReML") # This fn should be called only for non-stdd REML
  AUGI0_ZX <- as.list(locsXaug$AUGI0_ZX)
  locX <- AUGI0_ZX$X.pv
  if ( distinct.X.ReML[1L] ) {
    locX <- locX[,-(attr(X.Re,"unrestricting_cols"))]
  } 
  if ( distinct.X.ReML[2L] ) {
    extra_vars <- attr(X.Re,"extra_vars") ## may be NULL
    suppl_cols <- X.Re[,extra_vars, drop=FALSE]
    locX <- cbind(locX,suppl_cols)
  }
  AUGI0_ZX$X.pv <- locX
  locsXaug <- def_AUGI0_ZX_sparsePrecision(AUGI0_ZX = list2env(AUGI0_ZX),
                                           w.ranef=attr(locsXaug,"w.ranef"),
                                           cum_n_u_h=attr(locsXaug,"cum_n_u_h"),
                                           w.resid=attr(locsXaug,"w.resid"),
                                           corrPars=attr(locsXaug,"corrPars") )
  return(locsXaug)
}


# function to get the hatvalues (only: not the other similar computations on t_Q_scaled)
# no permutation issues for Q => a single get_hatvalues function should handle all sXaug classes
.get_hatvalues_MM <- function(sXaug, X.Re, weight_X, B=c("phi", "lambda")) {
  if ( is.null(X.Re)) { #<REML standard>
    hatval <- get_from_MME(sXaug,which="hatval") # colSums(t_Q_scaled*t_Q_scaled) ## basic REML, leverages from the same matrix used for estimation of betaV 
  } else if ( ncol(X.Re)==0L) { #<ML standard>
    hatval <- get_from_MME(sXaug,which="hatval_Z", B=B)
  } else {#<non-standard REML>
    distinct.X.ReML <- attr(X.Re,"distinct.X.ReML")
    if (inherits(sXaug,"AUGI0_ZX_sparsePrecision")) {
      if (any(distinct.X.ReML)) {
        locsXaug <- .calc_sXaug_Re_spprec(locsXaug=sXaug,X.Re) 
        hatval <- get_from_MME(locsXaug,which="hatval") 
      } else hatval <- get_from_MME(sXaug,which="hatval") # etaFix with formula=REMLformula: REML de factor standard by non-standard REML syntax
    } else { ## not spprec: code with shortcuts and checks but less clear
      if ( distinct.X.ReML[2L] ) { 
        locsXaug <- .calc_sXaug_Re(locsXaug=sXaug,X.Re,weight_X) 
        hatval <- .leveragesWrap(locsXaug) ## Rcpp version of computation through computation of Q
      } else if ( distinct.X.ReML[1L] ) { 
        whichcols <- attr(X.Re,"unrestricting_cols")
        if (length(whichcols)==attr(sXaug,"pforpv")) { ## should be ML standard
          stop("Ideally this case is not reached") # ("hatval_Z",B=B would be the arguments for ML)
        } else { ## non-standard case
          t_Q_scaled <- get_from_MME(sXaug,which="t_Q_scaled")
          n_u_h <- attr(sXaug,"n_u_h")
          ## substract cols directly from Q ! => t_Q_scaled must have cols in the order that give the "hatval"by get_from_MME(.,which="hatval")
          t_Q_scaled <- t_Q_scaled[-(n_u_h+whichcols),] ## test TRUE for standard ML 
          ## [, -integer(0)] would empty the matrix...
          hatval <- colSums(t_Q_scaled*t_Q_scaled)
        }
      } else { # etaFix with formula=REMLformula: REML de factor standard by non-standard REML syntax
        hatval <- get_from_MME(sXaug,which="hatval") 
      }      
    }
  }
  if (is.list(hatval)) hatval <- .unlist(hatval) ## assuming order lev_lambda,lev_phi
  return(hatval)
}

.get_hatvalues_FM <- function(X.Re, augX, w.resid) { ## for (G)LM
  if ( is.null(X.Re) ) { ## basic REML, leverages from the same matrix used for estimation of beta
    if (ncol(augX)) { 
      wAugX <- .calc_wAugX(XZ_0I=augX,sqrt.ww=sqrt(w.resid)) # rWW %*% X.pv 
      lev_phi <- .leveragesWrap(wAugX)
    } else { 
      lev_phi <- rep(0,nrow(augX)) ## leveragesWrap() -> .leverages() would fail on 0-col matrix? or it is large nrow the problem ?
    }
  } else if (ncol(X.Re)) { ## non standard REML
    wAugXleve <- .calc_wAugX(XZ_0I=X.Re,sqrt.ww=sqrt(w.resid)) # rWW%*%X.Re 
    lev_phi <- .leveragesWrap(wAugXleve)
  } else { # ML: X.Re non NULL mais ncol(X.Re)=0
    lev_phi <- rep(0,nrow(X.Re)) ## leveragesWrap() -> .leverages() would fail on 0-col matrix? or it is large nrow the problem ?
  }
}

.leveragesWrap <- function(X) { 
  # Matrix::qr.Q fails is X had zero columns. The call to this function assumes ncol>0
  if (inherits(X,"sparseMatrix")) {
    return(rowSums(qr.Q(qr(X))^2)) ## perhaps not optimal. Use hlfit <- HLfit(y ~ x, data=data.test) to profile it.
  } else .leverages(X) ## requests thinQ from Eigen QR
}

.calc_neg_d2f_dv_dloglam <- function(dlogfthdth, cum_n_u_h, lcrandfamfam, rand.families, u_h) {
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
  neg.d2f_dv_dloglam <- .unlist(neg.d2f_dv_dloglam)
  return(as.vector(neg.d2f_dv_dloglam))
}

.calc_dvdloglamMat_new <- function(neg.d2f_dv_dloglam,
                                   sXaug, d2hdv2_info=NULL ## use either one
                                   ) {
  if(is.null(d2hdv2_info)) {
    if (TRUE) {
      if (is.matrix(sXaug)) { # sXaug_EigenDense_QRP_Chol_scaled case
        dvdloglamMat <- get_from_MME(sXaug,"solve_d2hdv2",B=diag( neg.d2f_dv_dloglam))
      } else  dvdloglamMat <- get_from_MME(sXaug,"solve_d2hdv2",B=Diagonal(x= neg.d2f_dv_dloglam))
    } else {
      # Avoids the inelegant test is.matrix(sXaug), but... less accurate!
      inv_d2hdv2 <- get_from_MME(sXaug,"solve_d2hdv2") # slow step ("test negbin1" with large nobs is good example)
      dvdloglamMat <- .m_Matrix_times_Dvec(inv_d2hdv2, neg.d2f_dv_dloglam)# get_from_MME(sXaug,"solve_d2hdv2",B=diag( neg.d2f_dv_dloglam)) ## square matrix, by  the formulation of the algo 
    }
  } else if (inherits(d2hdv2_info,"CHMfactor")) {# CHM of ***-*** d2hdv2
    dvdloglamMat <- solve(d2hdv2_info, Diagonal(x= - neg.d2f_dv_dloglam ))  # rXr !       # .symDiagonal() is equivalent here
  } else if (inherits(d2hdv2_info,"qr") || inherits(d2hdv2_info,"sparseQR") ) { ## much slower than using CHMfactor
    if (length(neg.d2f_dv_dloglam)>5000L) message("[one-time solve()ing of large matrix, which may be slow]") 
    dvdloglamMat <- solve(d2hdv2_info, diag( neg.d2f_dv_dloglam ))  # rXr !       
  } else if (is.environment(d2hdv2_info)) {
    # dvdloglamMat <- solve(d2hdv2_info, diag( neg.d2f_dv_dloglam ))  # rXr !       
    rhs <- .Matrix_times_Dvec(d2hdv2_info$chol_Q, neg.d2f_dv_dloglam )
    rhs <- solve(d2hdv2_info$G_CHMfactor, rhs)
    dvdloglamMat <- - .crossprod(d2hdv2_info$chol_Q, rhs) # don't forget '-'
  } else { ## then d2hdv2_info is ginv(d2hdv2) or some other form of inverse
    dvdloglamMat <- .m_Matrix_times_Dvec(as(d2hdv2_info, "dgCMatrix"), # otherwise Matrix_times_Dvec() with dsC defaults detect a problem
                                         neg.d2f_dv_dloglam) ## sweep(d2hdv2_info,MARGIN=2L,neg.d2f_dv_dloglam,`*`) ## ginv(d2hdv2) %*% diag( as.vector(neg.d2f_dv_dloglam))      
  }
  #return(as.matrix(dvdloglamMat)) ## quite dense even if many small values and we subset it;
  return(dvdloglamMat) ## ./. but as.matrix() is terribly inefficient in bigranefs case (_F I X M E_ more adaptive code?)
} ## square matrix, by  the formulation of the algo 


.calc_dvdlogphiMat_new <- function(dh0deta,ZAL,
                                   sXaug,d2hdv2_info=NULL ## either one
                                   ) {
  ## cf calcul dhdv, but here we want to keep each d/d phi_i distinct hence not sum over observations i 
  if (inherits(ZAL,"ZAXlist")) { # appeared in first try gaussian("logit")... with ARp()...
    neg.d2h0_dv_dlogphi <- .crossprod(ZAL,Diagonal(x=drop(dh0deta)))
  } else neg.d2h0_dv_dlogphi <- .m_Matrix_times_Dvec(t(ZAL), drop(dh0deta)) # n_u_h*nobs: each ith column is a vector of derivatives wrt v_k# dh0dv <- t(ZAL) %*% diag(as.vector(dh0deta)) 
  if (is.null(d2hdv2_info)) { # call by HLfit_body
    dvdlogphiMat <- get_from_MME(sXaug,"solve_d2hdv2",B=neg.d2h0_dv_dlogphi) 
  } else if (inherits(d2hdv2_info,"CHMfactor")) { # CHM of ***-*** d2hdv2
    dvdlogphiMat <- solve(d2hdv2_info, - neg.d2h0_dv_dlogphi) # efficient without as.matrix()!        
  } else if (inherits(d2hdv2_info,"qr") || inherits(d2hdv2_info,"sparseQR") ) {
    dvdlogphiMat <- solve(d2hdv2_info, as.matrix(neg.d2h0_dv_dlogphi))  # rXn       
  } else if (is.environment(d2hdv2_info)) {
    # dvdlogphiMat <- d2hdv2_info %*% neg.d2h0_dv_dlogphi # rXn     
    rhs <- d2hdv2_info$chol_Q %*% neg.d2h0_dv_dlogphi
    rhs <- solve(d2hdv2_info$G_CHMfactor,rhs)
    dvdlogphiMat <- - .crossprod(d2hdv2_info$chol_Q,rhs) # don't forget '-'
  } else { ## then d2hdv2_info is ginv(d2hdv2) or a sparse matrix inverse of (d2hdv2) (spprec code will provide a dsCMatrix)
    if (inherits(d2hdv2_info,"dsCMatrix")) {
      d2hdv2_info <- as(d2hdv2_info,"dgeMatrix") ## more efficient if inv_d2hdv2 is math-dense
      # It would be nice to store only the half matrix but then as( - d2hdv2_info, "dpoMatrix") and reversing sign afterwards. 
    }
    dvdlogphiMat <- d2hdv2_info %*% neg.d2h0_dv_dlogphi # rXn       
  }
  return(dvdlogphiMat)
}


.calc_sscaled_new <- function(vecdisneeded, dlogWran_dv_h, coef12, 
                              n_u_h, sXaug, ZAL,WU_WT) {
  if (any(vecdisneeded)) { ## but the call to .calc_sscaled_new is conditional to the same condition 
    ## here version 1.5.3 had an interesting signed.wAugX concept
    vecdi1 <- vecdi2 <- vecdi3 <- 0
    ## P is P in LeeL appendix p. 4 and is P_R in MolasL p. 3307; X cols are excluded.
    Pdiag <- get_from_MME(sXaug, which="hatval_Z", B=unique(c("phi","phi","lambda")[which(vecdisneeded)])) # currently only spprec takes advantage of this.
    if (vecdisneeded[1L]) vecdi1 <- Pdiag$lev_phi * coef12$coef1 # coef1 is the factor of P_ii in d1
    # K2 = solve(d2hdv2,tZAL) is K2 matrix in LeeL appendix p. 4 and is -D in MolasL p. 3307 
    # W is Sigma^-1 ; TWT = t(ZALI)%*%W%*%ZALI = ZAL'.Wresid.ZAL+Wranef = -d2hdv2 !
    if (vecdisneeded[2L]) { # ( ZAL %*% K2 ) is K1 in LeeL appendix p. 4 and is A=-ZD in MolasL p. 3307-8 
      # vecdi2 <- as.vector( ((Pdiag$lev_phi * coef2) %*% ZAL) %*% K2)
      coef2 <- coef12$dlW_deta # coef2 is the factor between P_jj and K1 in d2
      vecdi2 <- get_from_MME(sXaug,"solve_d2hdv2",B=as.vector((Pdiag$lev_phi * coef2) %*id% ZAL))
      vecdi2 <- as.vector(ZAL %*% vecdi2) ## equiv  post-multiplication by Z^T in the expression for D p.3307 bottom.
      
      if ( ! is.null(WU_WT)) vecdi2 <- vecdi2/WU_WT # zero-truncated model: final factor in A11 in B4
    }
    # coef3 =(1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
    if (vecdisneeded[3L]) {  ## d3 reste nul pour gaussian ranef
      # vecdi3 <- as.vector( (Pdiag$lev_lambda * dlogWran_dv_h[seqn_u_h]) %*% K2)
      seqn_u_h <- seq_len(n_u_h)
      vecdi3 <- get_from_MME(sXaug,"solve_d2hdv2",B=as.vector(Pdiag$lev_lambda * dlogWran_dv_h[seqn_u_h]))
      vecdi3 <- as.vector(ZAL %*% vecdi3) ## equiv  post-multiplication by Z^T in the expression for D p.3307 bottom.
      if ( ! is.null(WU_WT)) vecdi3 <- vecdi3/WU_WT # zero-truncated model: final factor in A22 in B4
    }
    vecdi <- vecdi1+vecdi2+vecdi3 ## k_i in MolasL; le d_i de LeeL app. p. 4
    sscaled <- vecdi /2  ## sscaled := detadmu s_i= detadmu d*dmudeta/2 =d/2 in LeeL12; or dz1 = detadmu (y*-y) = detadmu m_i=0.5 k_i dmudeta = 0.5 k_i in MolasL 
    if ( ! is.null(WU_WT)) sscaled <- sscaled * WU_WT
  } else sscaled <- 0
  return(sscaled)
}

.init_resp_z_corrections_new <- function(lcrandfamfam, w.ranef, nobs, nrand, cum_n_u_h, rand.families, u_h, lambda_est, psi_M, v_h, dvdu, 
                                         sXaug, ZAL, w.resid) {
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
    psi_corr <- .unlist(psi_corr)
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

# b is a vector !
.Sig_times_b <- function(Sig0,ZAL,w.ranef,w.resid,b) { # Sig= [Sig0=Z.(1/w.ranef).Z^t+1/w.resid]
  if (is.null(Sig0)) { ## w.ranef is variable
    v1 <- .crossprod(ZAL, b) # drop(t(b) %*% ZAL) # drop() or .crossprod OK if b is effectively a vector
    v1 <- ZAL %*% ( v1 /w.ranef)
  } else {
    v1 <- Sig0 %*% b
  }
  v2 <- b/w.resid
  return(as.numeric(v1+v2))
}

# derived from .calc_dvdlogphiMat_new() and NOT USED but handy. (but does this work when ZAL is ZAX_list?)
.calc_lhs_InvSigma_rhs <- function(lhs, rhs=t(lhs), invV_factors, w.resid) {
  ## next lines use invV= w.resid- n_x_r %*% r_x_n
  resu <- lhs %*% .Dvec_times_m_Matrix(w.resid, rhs)
  resu <- resu - (lhs %*% invV_factors$n_x_r) %*% (invV_factors$r_x_n %*% rhs)
  return(resu)
}


.calc_z1 <- function(muetablob, w.resid, y, off, cum_nobs) { # (__FIXME__) if y and off were lists, I would not need resp_range etc.
  if (is.list(w.resid)) {
    if (is.null(mvlist <- w.resid$mvlist)) {
      z1 <- as.vector(muetablob$sane_eta+w.resid$WU_WT*(y-muetablob$mu-w.resid$dlogMthdth)/muetablob$dmudeta-off) ## MolasL10
    } else {
      z1s <- vector("list",length(mvlist)) 
      for (mv_it in seq_along(mvlist)) {
        resp_range <- .subrange(cumul=cum_nobs, it=mv_it)
        z1s[[mv_it]] <- .calc_z1(muetablob=muetablob$mv[[mv_it]], w.resid=w.resid$mvlist[[mv_it]], y=y[resp_range], off=off[resp_range])
      }
      z1 <- .unlist(z1s)
    }
  } else z1 <- as.vector(muetablob$sane_eta+(y-muetablob$mu)/muetablob$dmudeta-off) ## LeeNP 182 bas. GLM-adjusted response variable; O(n)*O(1/n)
  return(z1)  
}

.calc_zAug_not_LMM <- function(n_u_h, nobs, pforpv, y, off, ZAL, 
                      # variable within fit_as_ZX:
                      muetablob, dlogWran_dv_h, sXaug, w.resid, w.ranef, 
                      ########################### ZAL_scaling,
                      init_z_args, 
                      #
                      processed) {
  GLMMbool <- attr(processed[["models"]],"GLMMbool") 
  ######## According to 'theorem 1' of LeeL12, new beta estimate from z1-a(i), where z1 is
  z1 <- .calc_z1(muetablob, w.resid, y, off, cum_nobs=attr(processed$families,"cum_nobs"))
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
    ######################## ZAL <- .m_Matrix_times_Dvec(ZAL, ZAL_scaling)
    vecdisneeded <- processed$vecdisneeded # vecdisneeded <- c( coef12needed, coef12needed, any(dlogWran_dv_h!=0L) )
    if (any(vecdisneeded)) {
      if (is.list(w.resid)) {
        WU_WT <- w.resid$WU_WT 
      } else WU_WT <- NULL
      sscaled <- .calc_sscaled_new(vecdisneeded=vecdisneeded,
                              dlogWran_dv_h=dlogWran_dv_h, ## dlogWran_dv_h was computed when w.ranef was computed
                              coef12= .calc_dlW_deta(muetablob,
                                                    family=processed$family, 
                                                    BinomialDen=processed$BinomialDen, 
                                                    calcCoef1=TRUE,
                                                    w.resid=w.resid, families=processed$families), ## promise evaluated if any vecdisneeded[-3]
                              n_u_h=n_u_h, 
                              sXaug=sXaug,
                              ZAL=ZAL, # vecdi2
                              WU_WT=WU_WT ## NULL except for truncated model
      )
      if (is.list(w.resid)) { # both the truncated and the mv cases
        y2_sscaled <- z2+ as.vector((sscaled * w.resid$w_resid ) %*% ZAL )/w.ranef ## that's the y_2 in "Methods of solution based on the augmented matrix"
        # it is unaffected by the matrix rescaling bc it is a fn of z1 and z2. But rescaled is alays taken into account bc we uuse y2_sscaled only 
        #      in the context wzAug <- c(zInfo$y2_sscaled/ZAL_scaling, (zInfo$z1_sscaled)*weight_X)
      } else y2_sscaled <- z2+ as.vector((sscaled * w.resid ) %*% ZAL )/w.ranef
    } else { # notably after observing that general code with sscaled=0 and large ZAL is slow!
      sscaled <- 0
      y2_sscaled <- z2
    }
    zInfo <- list(sscaled=sscaled, z1=z1, z2=z2, z1_sscaled=z1-sscaled, y2_sscaled=y2_sscaled)
  } else zInfo <- list(sscaled=0, z1=z1, z2=z2, z1_sscaled=z1, y2_sscaled=z2) 
  return(zInfo)
}

.oldcbind_dgC_dgC <- function(leftcols, rightcols) { # expects @x,i,p => dgCMatrix
  leftcols@p <- c(leftcols@p, leftcols@p[length(leftcols@p)] + rightcols@p[-1L])
  leftcols@i <- c(leftcols@i, rightcols@i) 
  leftcols@x <- c(leftcols@x, rightcols@x)
  if (is.null(leftcols@Dimnames[[2L]])) {
    if ( ! is.null(rightcols@Dimnames[[2L]])) {
      leftcols@Dimnames[[2L]] <- c(rep("",leftcols@Dim[2L]),rightcols@Dimnames[[2L]])
    } ## else all colnames are NULL
  } else {
    if ( is.null(rightcols@Dimnames[[2L]])) {
      leftcols@Dimnames[[2L]] <- c(leftcols@Dimnames[[2L]],rep("",rightcols@Dim[2L]))
    } else leftcols@Dimnames[[2L]] <- c(leftcols@Dimnames[[2L]],rightcols@Dimnames[[2L]])
  } 
  leftcols@Dim[2L] <- leftcols@Dim[2L]+rightcols@Dim[2L]
  # the old code should have included:
  leftcols@factors <- list()
  attr(leftcols,"is_incid") <- FALSE
  attr(leftcols,"namesTerm") <- NULL
  return(leftcols)
}

# Sligthly faster than the old, pure R version:
.cbind_dgC_dgC <- function(leftcols, rightcols) {
  res <- .RcppMatrixCb2(leftcols, rightcols)
  if (is.null(leftcols@Dimnames[[2L]])) {
    if ( ! is.null(rightcols@Dimnames[[2L]])) {
      res@Dimnames[[2L]] <- c(rep("",leftcols@Dim[2L]),rightcols@Dimnames[[2L]])
    } ## else all colnames are NULL
  } else {
    if ( is.null(rightcols@Dimnames[[2L]])) {
      res@Dimnames[[2L]] <- c(leftcols@Dimnames[[2L]],rep("",rightcols@Dim[2L]))
    } else res@Dimnames[[2L]] <- c(leftcols@Dimnames[[2L]],rightcols@Dimnames[[2L]])
  } 
  return(res)
}


.adhoc_cbind_dgC_0 <- function(leftcols, newcoln) { # expects @x,i,p => dgCMatrix
  leftcols@p <- c(leftcols@p, rep(leftcols@p[length(leftcols@p)],newcoln))
  if ( ! is.null(leftcols@Dimnames[[2L]])) leftcols@Dimnames[[2L]] <- c(leftcols@Dimnames[[2L]],rep("",newcoln)) 
  leftcols@Dim[2L] <- leftcols@Dim[2L]+newcoln
  return(leftcols)
}

.adhoc_rbind_I_dgC <- function(Ilen, ZAL) {
  newlen <- Ilen+length(ZAL@x)
  Iseq <- seq_len(Ilen)
  Ip <- c(0L,Iseq)
  newp <- Ip+ZAL@p
  Ipos <- newp[-length(newp)]+1L
  #
  newx <- numeric(newlen)
  newx[Ipos] <- 1 # "I@x"
  newx[-Ipos] <- ZAL@x
  newi <- integer(newlen)
  newi[Ipos] <- Iseq-1L
  newi[-Ipos] <- ZAL@i+Ilen
  #
  ZAL@i <- newi
  ZAL@x <- newx
  ZAL@p <- newp
  ZAL@Dim[1L] <- Ilen+ZAL@Dim[1L]
  if ( ! is.null(ZAL@Dimnames[[1L]])) ZAL@Dimnames[[1L]] <- c(rep("",Ilen),ZAL@Dimnames[[1L]]) 
  return(ZAL)
}

.make_Xscal <- function(ZAL, ZAL_scaling=NULL, processed, as_matrix) {
  if (inherits(ZAL,"ZAXlist")) ZAL <- .ad_hoc_cbind(ZAL@LIST, as_matrix=as_matrix )
  # capture programming error for ZAL_scaling:
  if (length(ZAL_scaling)==1L && ncol(ZAL)!=1L) stop("ZAL_scaling should be a full-length vector, or NULL. Contact the maintainer.")
  # ncol(ZAL)=1L could occur in 'legal' (albeit dubious) use. The total number of levels of random effects has been checked in preprocessing.
  if ( ! is.null(ZAL_scaling)) ZAL <- .m_Matrix_times_Dvec(ZAL,ZAL_scaling)
  AUGI0_ZX <- processed$AUGI0_ZX
  if (is.null(Zero_sparseX <- AUGI0_ZX$Zero_sparseX)) Zero_sparseX <- rbind2(AUGI0_ZX$ZeroBlock, AUGI0_ZX$X.pv)
  if (inherits(ZAL,"dgCMatrix")) {
    I_ZAL <- .adhoc_rbind_I_dgC(nrow(AUGI0_ZX$I), ZAL) ## this is faster...
  } else I_ZAL <- rbind2(AUGI0_ZX$I, ZAL)
  if (inherits(I_ZAL,"dgCMatrix") &&  inherits(Zero_sparseX,"dgCMatrix") ) {
    Xscal <- .cbind_dgC_dgC(I_ZAL, Zero_sparseX) # substantially faster than the general alternative 
  } else Xscal <- cbind2(I_ZAL, Zero_sparseX)
  attr(Xscal,"AUGI0_ZX") <- AUGI0_ZX # environment => cheap access to its 'envir$updateable' variable or anything else 
  return(Xscal)
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

.adhoc_cbind_dgC_sXaug_pxy_o <- function(sXaug, pwy_o, n_u_h) {
  I00_ZXy <- sXaug
  I00_ZXy@p <- c(I00_ZXy@p, I00_ZXy@p[length(I00_ZXy@p)]+length(pwy_o))
  I00_ZXy@i <- c(I00_ZXy@i, n_u_h-1L+seq_len(length(pwy_o))) ## fails if n_u_h is not integer
  I00_ZXy@x <- c(I00_ZXy@x, pwy_o)
  I00_ZXy@Dim[2L] <- I00_ZXy@Dim[2L]+1L
  I00_ZXy@Dimnames[[2L]] <- c(I00_ZXy@Dimnames[[2L]],"") ## otherwise try(chol()=> error)  (which makes a test of the rescue code...)
  return(I00_ZXy)
}

.get_R_aug_ZXy <- function(aug_ZXy, augZXy_solver, return_tri) {
  # Currently using only the diagonal (though not simply the logdet) => tri is important, lower or upper OK.  BUT .../...
  # .../... actually it's not true: I use its t(solve(.)) in a subcase
  nc <- ncol(aug_ZXy)
  solver <- augZXy_solver[1L]
  if (solver =="chol") {
    R_aug_ZXy <- try(chol(.crossprod(aug_ZXy)), silent=TRUE)
    if ( ! inherits(R_aug_ZXy,"try-error")) return(R_aug_ZXy)
    solver <- augZXy_solver[2L]
    if (is.na(solver)) solver <- "EigenQR"
  } else if (solver=="QR") solver <- "EigenQR" ## explicitation of current default meaning of "QR"
  if (solver =="EigenQR") {
    if (inherits(aug_ZXy,"Matrix")) { 
      # If ZXy is 'tall' then $R will have the correct size, but if it is 'wide', Eigen  returns a square matrix with the wide dimension, 
      # .lmwithQR's last rows contains noise (variable between different calls!) that impacts the crossprod (see example in devel/Eigen), so these rows should be removed.
      # .lmwith_sparse_QRp may be less likely to generate such noise but let's be consistent
      qrblob <- .lmwith_sparse_QRp(aug_ZXy,yy=NULL,returntQ=FALSE,returnR=TRUE)
      R_aug_ZXy <- qrblob$R
      if (nrow(aug_ZXy)<ncol(aug_ZXy)) R_aug_ZXy <- R_aug_ZXy[seq_len(nrow(aug_ZXy)),]
      if ( ! all(unique(diff(qrblob$P))==1L)) {
        R_aug_ZXy <- R_aug_ZXy[,sort.list(qrblob$P)] ## not triangular
        if (return_tri) { # eval an unpermuted triangular R
          R_aug_ZXy <- .lmwithQR(as.matrix(R_aug_ZXy) ,yy=NULL,returntQ=FALSE,returnR=TRUE)$R
        }
      }
    } else R_aug_ZXy <- .lmwithQR(aug_ZXy,yy=NULL,returntQ=FALSE,returnR=TRUE)$R
  } else if (solver =="qr") { ## tries base qr but checks pivoting, with fallback
    if (inherits(aug_ZXy,"Matrix")) {
      qrblob <- qr(aug_ZXy)
      R_aug_ZXy <- qrR(qrblob,backPermute=TRUE) ## not triangular
      if ( return_tri && ! all(unique(diff(qrblob@q))==1L)) { # eval an unpermuted triangular R
        R_aug_ZXy <- .lmwithQR(as.matrix(R_aug_ZXy) ,yy=NULL,returntQ=FALSE,returnR=TRUE)$R ## upper tri
      }
    } else {
      qrblob <- qr(aug_ZXy)
      R_aug_ZXy <- qr.R(qrblob)
      if ( return_tri && ! all(unique(diff(qrblob$pivot))==1L)) { # eval an unpermuted triangular R
        R_aug_ZXy <- .lmwithQR(R_aug_ZXy[, sort.list(qrblob$pivot)] ,yy=NULL,returntQ=FALSE,returnR=TRUE)$R ## upper tri
      } 
    }
  } else stop("unknown 'augZXy_solver' requested.")
  return(R_aug_ZXy)
}

# This function belongs to SPCORR methods (inherits(ZAL,"sparseMatrix") -> class(sXaug) <- c(class(sXaug),"sXaug_blocks")). 
# It gives (mostly) logdet terms of blocks of the chol facto of the crossproduct of the y-augm-augm weighted design matrix.
# The algo avoids forming the chol of the crossproduct and even the crossproduct itself), instead starting from 
# the more easily available Cholesky facto 'CHM_ZZ' of its upper left block Z'WZ+diag(). To obtain the logdet of the R_XX block
# of the chol factor, the crossproduct 'cross_Rxx' of this block is formed (not to be confused with the XX block of the crossproduct 
# of the y-augm-augm thing) as a difference of crossproducts of small|slim-as-X terms, and determinant(cross_Rxx) 
# is computed (no explicit chol facto; this is a small block). The 'r_yy' block is trivially 1*1 and it's its square 'ryy2' 
# which is reported (rather that its trivial logdet), computed also as (essentially) a difference of crossproducts of 
# 1-col terms. Forming these differences is however a bit clumsy.
.get_absdiagR_blocks <- function(sXaug_blocks, pwy_o, n_u_h, processed, augZXy_solver,update_info) { # for SPCORR !
  seq_n_u_h <- seq(n_u_h)
  tZW <- t(sXaug_blocks$ZW) # actually a ZL rather than a Z.
  if (is.null(template <- processed$AUGI0_ZX$template_CHM_ZZ_blocks)) { 
    cross_Z <- .tcrossprod(tZW) 
    if (inherits(cross_Z,"dsyMatrix")) { ## Matrix considered the matrix as effectively dense
      message(paste("Possibly poor selection of methods: Z stored as sparse, but Z'Z assessed as dense by Matrix's as(., 'symmetricMatrix').",
                    "control.HLfit(algebra='decorr') may be used to control this on a one-time, ad-hoc basis.",
                    ## see comments in .choose_QRmethod(). We may reach this block whenever the correlation matrix is dense.
                    collapse="\n"))
      cross_Z <- as(cross_Z,"sparseMatrix")
    } 
    CHM_ZZ <- Cholesky(cross_Z, perm=TRUE, LDL=FALSE, Imult=1) # Imult !
    if (update_info$allow) { 
      processed$AUGI0_ZX$template_CHM_ZZ_blocks <- CHM_ZZ
    }
  } else CHM_ZZ <- Matrix::.updateCHMfactor(template, tZW, mult=1) # no need to compute the crossprod for updating: the tcrossfac is enough.
  # perm <- as(CHM_ZZ,"pMatrix")@perm # remarkably slow... and using   perm <- CHM_ZZ@perm+1L + [perm,] is much slower than:
  ZtWy <- tZW %*% pwy_o
  r_Zy <- solve(CHM_ZZ, solve(CHM_ZZ,ZtWy,system="P"), system="L") # solve(CHM_ZZ, ZtWy[perm], system="L")  # 
  #
  #cross_Rxx <- .crossprod(XW,as_mat=TRUE)-.crossprod(Rzx,as_mat=TRUE) # as(,"dpoMatrix) involves a chol() factorization...
  # Calling directly .Rcpp_crossprod avoids some bureaucratic overhead (irrespective of keep_names which could rather affect later computations)
  XW <- sXaug_blocks$XW
  if (ncol(XW)) { # there are fixed effects in X
    ZtWX <- tZW %*% XW
    Rzx <- solve(CHM_ZZ, solve(CHM_ZZ,ZtWX,system="P"), system="L") # (maybe) dense but dge... # solve(CHM_ZZ, ZtWX[perm,], system="L") # 
    cross_Rxx <- .Rcpp_crossprod(XW,BB=NULL, keep_names=FALSE,as_mat=TRUE) -
      .Rcpp_crossprod(Rzx,BB=NULL, keep_names=FALSE,as_mat=TRUE) 
    XtWy <- .Rcpp_crossprod(XW, pwy_o, keep_names=FALSE,as_mat=TRUE)
    u_of_quadratic_utAu <- XtWy-.Rcpp_crossprod(Rzx, r_Zy, keep_names=FALSE,as_mat=TRUE)
    if (TRUE) { # not clear why solve(cross_Rxx,.) would work and not chol() 
      chol_XX <- chol(cross_Rxx)
      r_Xy <- backsolve(chol_XX, u_of_quadratic_utAu, transpose=TRUE) ## I wrote "transpose since chol() provides a tcrossfac". ?? 
      r_Zy_x <- r_Zy@x
      ryy2 <- sum(pwy_o*pwy_o) - sum(r_Zy_x*r_Zy_x) - sum(r_Xy*r_Xy)
      absdiagR_terms <- list(logdet_v=determinant(CHM_ZZ)$modulus[1], 
                             logdet_b=sum(log(abs(.diagfast(chol_XX)))), ryy2=ryy2)
    } else if (use_crossr22 <- TRUE) { # a bit slower (even using .Rcpp_crossprod)
      # No need for the complex Utri_chol computation here, as sum(r_Xy^2) is easy to compute without it.
      # Another place where one can avoid it is also labelled 'use_crossr22' in .solve_crossr22()
      sum_r_Ry_2 <- .crossprod(u_of_quadratic_utAu, solve(cross_Rxx, u_of_quadratic_utAu))
      ryy2 <- sum(pwy_o^2) - sum(r_Zy^2) - sum_r_Ry_2
      absdiagR_terms <- list(logdet_v=determinant(CHM_ZZ)$modulus[1], 
                             logdet_b=determinant(cross_Rxx)$modulus[1]/2, ryy2=ryy2)
    } else {
      chol_XX <- .Utri_chol_by_qr(cross_Rxx) # chol(cross_Rxx) # chol_XX matrix is quite small -> same algos as in .calc_r22()
      # ## test-poly test-random-slope test-ranCoefs; 
      # test-random-slope  is slower by .Rcpp_backsolve() but only because of more precise, but longer, outer optim in (ares <- ...)
      # this better result is by accumulated effects on the optimization path rather than by functional improvement.
      # also .Rcpp_backsolve() visibly increases range(get_predVar(twolambda)[1:5]-get_predVar(onelambda)[1:5]) 
      r_Xy <- backsolve(chol_XX, u_of_quadratic_utAu, transpose=TRUE) ## transpose since chol() provides a tcrossfac 
      # .crossprod(Rzx, r_Zy) appears to be .crossprod(ZtWX, solve(CHM_ZZ, ZtWy, system = "A")) 
      # but we still need Rzx and r_Zy 
      r_Zy_x <- r_Zy@x
      ryy2 <- sum(pwy_o*pwy_o) - sum(r_Zy_x*r_Zy_x) - sum(r_Xy*r_Xy)
      absdiagR_terms <- list(logdet_v=determinant(CHM_ZZ)$modulus[1], 
                             logdet_b=sum(log(abs(.diagfast(chol_XX)))), ryy2=ryy2)
    }
  } else { # no fixed effects... 
    r_Zy_x <- r_Zy@x
    ryy2 <- sum(pwy_o*pwy_o) - sum(r_Zy_x*r_Zy_x)
    absdiagR_terms <- list(logdet_v=determinant(CHM_ZZ)$modulus[1], 
                           logdet_b=0, ryy2=ryy2)
  }
  return(absdiagR_terms)
}

.get_absdiagR <- function(aug_ZXy, augZXy_solver) {
  R_aug_ZXy <- .get_R_aug_ZXy(aug_ZXy, augZXy_solver,return_tri=TRUE)
  nc <- ncol(aug_ZXy)
  diagPos <- seq.int(1L,nc^2,nc+1L)
  return(abs(R_aug_ZXy[diagPos]))
}

.sum_pwt_Q_y_o_2 <- function(sXaug,pwy_o) {
  if (inherits(sXaug,"AUGI0_ZX_sparsePrecision")) {
    sum_pwt_Q_y_o_2 <- .calc_sum_pwt_Q_y_o_2(sXaug,pwy_o)
  } else {
    #pwt_Q_y_o <- get_from_MME(sXaug,"t_Q_scaled")%*% c(rep(0,n_u_h),pwy_o) 
    pwt_Q_y_o <- get_from_MME(sXaug,"Qt_leftcols*B", B=pwy_o)
    sum_pwt_Q_y_o_2 <- sum(pwt_Q_y_o^2)
  }
  sum_pwt_Q_y_o_2
}

.calc_APHLs_by_augZXy_or_sXaug <- function(processed, auglinmodblob=NULL, 
                                     sXaug, W00_R_qr_ZXy=NULL, which, phi_est,
                                     update_info) { # either auglinmodblob or (sXaug|W00_R_qr_ZXy) and (possibly NULL) phi_est
  resu <- list()
  if ( ! is.null(auglinmodblob)) {
    sXaug <- auglinmodblob$sXaug
    phi_est <- auglinmodblob$phi_est
  } 
  if (!is.null(W00_R_qr_ZXy)) {
    locattr <- attributes(W00_R_qr_ZXy)
  } else locattr <- attributes(sXaug)
  #weight_X <- locattr$weight_X
  H_global_scale <- locattr$H_global_scale 
  extranorm <- locattr$extranorm 
  #if (is.null(W00_R_qr_ZXy) && inherits(sXaug,"AUGI0_ZX_sparsePrecision")) { 
#    if (is.null(weight_X)) weight_X <- 1 ## spprec case
  weight_X <- 1 ## 05/12/2019 using weight_X in this fn is actually a 'bug' (adding constant term to objective, but not affecting optimization)
  if (is.null(H_global_scale)) H_global_scale <- 1 ## spprec case
  #}
  if (is.null(extranorm)) extranorm <- H_global_scale
  n_u_h <- locattr$n_u_h
  nobs <- length(processed$y)
  pforpv <- locattr$pforpv
  resdf <- nobs - pforpv
  #
  prior_weights <- eval(processed$prior.weights)
  if (is.null(phi_est)) { ## then we estimate a factor 'lamphifac" common to lambda and phi
    # in effect we fit for phi=1 then estimate lamphifac from a sum of squares for all augmented residuals.
    augZXy_solver <- .spaMM.data$options$augZXy_solver ## ie "chol", "EigenQR", etc.
    if ( ! is.null(augZXy_solver) && ! inherits(sXaug,"AUGI0_ZX_sparsePrecision")) { # use augZXy_solver
      if (! is.null(W00_R_qr_ZXy)) { # y-augmented factor available
        absdiagR <- .get_absdiagR(W00_R_qr_ZXy, augZXy_solver)
        absdiagR[seq(n_u_h)] <- absdiagR[seq(n_u_h)] /attr(W00_R_qr_ZXy,"eigen_s_invL") # equivalent to the |Omega| term in BatesD04
        nc <- length(absdiagR)
        pwSSE <- (absdiagR[nc]^2)/extranorm
        logdet_R_scaled_b_v <- sum(log(absdiagR[-nc]))
        X_scaled_H_unscaled_logdet_r22 <- sum(log(absdiagR)[-c(seq(n_u_h),nc)]) -pforpv*log(H_global_scale)/2 
        ## -pforpv*log(H_global_scale)/2 for consistency with get_from_MME(sXaug,"logdet_r22") assuming the latter is correct
      } else if (inherits(sXaug,"sXaug_blocks")) { # SPCORR !
        pwphi <- 1/(prior_weights) ## vector
        y_o <- drop(processed$y-processed$off)
        pwy_o <- y_o*sqrt(extranorm/pwphi)
        # .spaMM.data$options$ATUER <- FALSE
        # absdiagR_terms1 <- .get_absdiagR_new(sXaug, pwy_o, n_u_h, processed, 
        #                                      augZXy_solver=augZXy_solver,
        #                                      update_info=update_info) 
        # .spaMM.data$options$ATUER <- TRUE
        # absdiagR_terms <- .get_absdiagR_new(sXaug, pwy_o, n_u_h, processed, 
        #                                     augZXy_solver=augZXy_solver,
        #                                     update_info=update_info) 
        absdiagR_terms <- .get_absdiagR_blocks(sXaug_blocks=sXaug, pwy_o, n_u_h, processed, 
                                            augZXy_solver=augZXy_solver,
                                            update_info=update_info) 
        pwSSE <- absdiagR_terms$ryy2/extranorm
        logdet_R_scaled_b_v <- absdiagR_terms$logdet_v+absdiagR_terms$logdet_b
        X_scaled_H_unscaled_logdet_r22 <- absdiagR_terms$logdet_b -pforpv*log(H_global_scale)/2 
      } else { # y-augmented factor to be constructed from sXaug: .HLfit_body_augZXy, or check_augZXy code
        pwphi <- 1/(prior_weights) ## vector
        y_o <- (processed$y-processed$off)
        pwy_o <- y_o*sqrt(extranorm/pwphi)
        if (inherits(sXaug,"dgCMatrix")) {
          I00_ZXy <- .adhoc_cbind_dgC_sXaug_pxy_o(sXaug, pwy_o, n_u_h) ## distinctly faster
        } else if (is.numeric(sXaug)) { ## not in routine tests but in CAR_timings
          I00_ZXy <- .Rcpp_dense_cbind_mat_vec(sXaug, c(rep(0,n_u_h),pwy_o)) # typically costly ./.
          # as a big matrix must be allocated each time .calc_APHLs_by_augZXy_or_sXaug) is called.
          # This is where assignment in place in a stored template would be useful, but pure R will not avoid local copies. I tried
          # I00_ZXy <- .update_I00_ZXy(sXaug, pwy_o, n_u_h)
          # but this was slow.
        } else I00_ZXy <- cbind(sXaug,c(rep(0,n_u_h),pwy_o)) ## this cbind takes time...
        # Rcpp version of cbind for sparse matrices : https://stackoverflow.com/questions/45875668/rcpp-eigen-sparse-matrix-cbind#
        # but the gain is small...
        absdiagR <- .get_absdiagR(I00_ZXy, augZXy_solver)
        nc <- length(absdiagR)
        pwSSE <- (absdiagR[nc]^2)/extranorm
        logdet_R_scaled_b_v <- sum(log(absdiagR[-nc]))
        X_scaled_H_unscaled_logdet_r22 <- sum(log(absdiagR)[-c(seq(n_u_h),nc)]) -pforpv*log(H_global_scale)/2 
        ## -pforpv*log(H_global_scale)/2 for consistency with get_from_MME(sXaug,"logdet_r22") assuming the latter is correct
      }
    } else { ## other sXaug methods not using y-augmented factor: AUGI0_ZX_sparsePrecision or devel(?) code
      pwphi <- 1/(prior_weights) ## vector
      y_o <- (processed$y-processed$off)
      pwy_o <- y_o*sqrt(extranorm/pwphi) # extranorm is for better accuracy of next step
      sum_pwt_Q_y_o_2 <- .sum_pwt_Q_y_o_2(sXaug,pwy_o)
      pwSSE <- (sum(pwy_o^2)-sum_pwt_Q_y_o_2)/extranorm ## sum() : vectors of different lengths !
      logdet_R_scaled_b_v <- get_from_MME(sXaug,"logdet_R_scaled_b_v") # logdet_R_scaled_v+X_scaled_H_unscaled_logdet_r22 ## p_bv substract all of this and p_v cancels the r22 part 
      X_scaled_H_unscaled_logdet_r22 <- get_from_MME(sXaug,"logdet_r22") # if spprec: already available in BLOB from logdet_R_scaled_b_v computation...
    }
    # We obtain phi_est IN ANOTHER MODEL than in the general formulation as this phi also impacts the ranef variances
    ## SSE [sum of nobs+nr terms]/nobs provides an estimate of a scaling factor
    ## not of phi (which could be  sum((y-fitted)[ypos])^2)/sum(1-lev_phi)
    #we have fitted for the model (lambda, 1/prior_weights) and deduce the optimal (lamphifac*lambda, lamphifac/prior_weights)
    #The hatval are thus those both for phi and lambda whose sum is the #df
    # hatval <- .get_hatvalues_MM(sXaug,X.Re=processed$X.Re, weight_X) ## in case we need them...
    # devel code for prior weights removed from [v2.7.11
    #
    p_base <- sum(log(weight_X)) - logdet_R_scaled_b_v + pforpv*log(H_global_scale)/2 ## keep H_global_scale here even when it differs from extranorm
    if (is.null(processed$X.Re)) { # canonical REML
      resu$phi_est <- lamphifac_REML <- max(pwSSE/(resdf), 1e-6) ## check with pw ## remind We obtain phi_est IN ANOTHER MODEL than in the general formulation
      X_scaled_p_bv <- p_base - resdf * (1+log(2 * pi * lamphifac_REML))/2 
    } else {
      resu$phi_est <- lamphifac_ML <- max(pwSSE/(nobs), 1e-6) 
      # X_scaled_H_unscaled_logdet_r22 must have been previously computed  in all subcases where it is needed
      resu$p_v <- p_base + X_scaled_H_unscaled_logdet_r22 - nobs * (1+log(2 * pi * lamphifac_ML))/2 
    }
  } else { ## phi_est available; no lamphifac estimation; in particular for .makeCovEst1
    pwphi <- phi_est/prior_weights ## vectorize phi if not already vector
    pwy_o <- (processed$y-processed$off)/sqrt(pwphi/extranorm) # extranorm is for better accuracy of next step
    sum_pwt_Q_y_o_2 <- .sum_pwt_Q_y_o_2(sXaug,pwy_o)
    pwSSE <- (sum(pwy_o^2)-sum_pwt_Q_y_o_2)/extranorm ## vectors of different lengths !
    logdet_R_scaled_b_v <- get_from_MME(sXaug,"logdet_R_scaled_b_v")
    # we don't assume here that phi_est is at its MLE (in contrast to null-phi_est case => Bates's formulas)
    cliklike <- (pwSSE+sum(log(2*pi*pwphi)))/2
    if (FALSE) {
      p_base <- sum(log(weight_X)) - logdet_R_scaled_b_v + pforpv*log(2*pi*H_global_scale)/2 - cliklike ## keep  H_global_scale here even when it differs from extranorm
      if (is.null(processed$X.Re)) {
        X_scaled_p_bv <- p_base 
      } else { # we don't assume here that phi_est is at its MLE (in commparison to Bates's formulas)
        if ( ! inherits(sXaug,"AUGI0_ZX_sparsePrecision")) X_scaled_H_unscaled_logdet_r22 <- get_from_MME(sXaug,"logdet_r22") 
        resu$p_v <- p_base + X_scaled_H_unscaled_logdet_r22 - pforpv*log(2*pi)/2 
      }
      old_p_v <- resu$p_v
    } 
    if (FALSE) { ## FALSE TRUE TRUE => .816
      p_base <- - cliklike + sum(log(weight_X)) - logdet_R_scaled_b_v + pforpv*log(2*pi*H_global_scale)/2 ## keep  H_global_scale here even when it differs from extranorm
      if (is.null(processed$X.Re)) { # canonical REML
        X_scaled_p_bv <- p_base 
      } else {
        if ( ! inherits(sXaug,"AUGI0_ZX_sparsePrecision")) X_scaled_H_unscaled_logdet_r22 <- get_from_MME(sXaug,"logdet_r22") 
        resu$p_v <- p_base + X_scaled_H_unscaled_logdet_r22 - pforpv*log(2*pi)/2 
      }
    }
    if (TRUE) { ## optimization fitme6 etc. is sensitive to the smallest numerical errors... even affected by order of additions and subtrations 
      p_base <- - cliklike + sum(log(weight_X)) - logdet_R_scaled_b_v + pforpv*log(H_global_scale)/2 ## keep  H_global_scale here even when it differs from extranorm
      if (is.null(processed$X.Re)) { # canonical REML
        X_scaled_p_bv <- p_base + pforpv*log(2*pi)/2
      } else {
        resu$p_v <- p_base + get_from_MME(sXaug,"logdet_r22") 
      }
    }
  }

  if ("p_bv" %in% which) {
    if (is.null(processed$X.Re)) { # canonical REML
      if ( ! is.null(X_scale <- attr(processed$AUGI0_ZX$X.pv,"scaled:scale"))) {
        resu$p_bv <- X_scaled_p_bv - sum(log(X_scale))
      } else resu$p_bv <- X_scaled_p_bv
    } else resu$p_bv <- resu$p_v
  }
  return(resu)
}

.test_augZXy <- function(augZXy_resu, augZX_resu,phi.Fix, phi_est) {
  if (!is.null(augZXy_phi <- augZXy_resu$phi_est)) { ## ie was estimted by the augZXy method
    cat("dphi:", augZXy_phi-phi_est)
  }
  if ( ! is.null(p_v <- augZXy_resu$p_v)) {
    cat("p_v:", p_v)
    zut <- abs(p_v-augZX_resu$p_v)
    if (zut>1e-6) browser()
  }
  if ( ! is.null(p_bv <- augZXy_resu$p_bv)) {
    cat("p_bv:", p_bv)
    zut <- abs(p_bv-augZX_resu$p_bv)
    if (zut>1e-6) browser()
  }
}

.calc_APHLs_from_ZX <- function(auglinmodblob=NULL,processed, which="p_v",
                               ## alternative to auglinmodblob, insuff pour REML non standard:
                               sXaug, phi_est, lambda_est, dvdu, u_h, muetablob
                               ) {
  augZX_resu <- list()
  if ( ! is.null(auglinmodblob)) {
    sXaug <- auglinmodblob$sXaug 
    muetablob <- auglinmodblob$muetablob
    phi_est <- auglinmodblob$phi_est
    u_h <- auglinmodblob$u_h
    lambda_est <- auglinmodblob$lambda_est
    dvdu <- auglinmodblob$wranefblob$dvdu
  }
  mu <- muetablob$mu
  #
  augZX_resu$clik <- .calc_clik(mu,phi_est,processed, 
                                muetaenv=muetablob) # muetaenv used in COMPoisson case
  if (all(which =="clik")) return(augZX_resu)
  if (processed$models[["eta"]]=="etaGLM") {
    augZX_resu$p_v <- augZX_resu$clik
    return(augZX_resu)
  } # E L S E 
  cum_n_u_h <- processed$cum_n_u_h
  lcrandfamfam <-  attr(processed$rand.families,"lcrandfamfam")
  likranU <- vector("list",length(lcrandfamfam))
  for (it in seq_len(length(lcrandfamfam))) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    likranU[[it]] <- .loglfn_ranU(lcrandfamfam[it],u_h[u.range],1/lambda_est[u.range])
  }
  likranU <- .unlist(likranU)
  log.du_dv <- - log(dvdu) 
  likranV <- sum(likranU + log.du_dv)
  augZX_resu$hlik <- augZX_resu$clik + likranV
  #
  n_u_h <- length(lambda_est)
  # beware that computation of logdet_sqrt_d2hdv2 depends on w.ranef
  if ("p_v" %in% which || "p_bv" %in% which) {
    augZX_resu$p_v <- augZX_resu$hlik - get_from_MME(sXaug,"logdet_sqrt_d2hdv2") + n_u_h*log(2*pi)/2
  }
  if ("p_bv" %in% which) {
    X.Re <- processed$X.Re
    if ( is.null(X.Re)) {## REML standard
      pforpv <- attr(sXaug,"pforpv")
      X_scaled_H_unscaled_logdet_r22 <- get_from_MME(sXaug,"logdet_r22") 
      augZX_resu$p_bv <- augZX_resu$p_v - X_scaled_H_unscaled_logdet_r22 + pforpv*log(2*pi)/2 
      if ( ! is.null(X_scale <- attr(processed$AUGI0_ZX$X.pv,"scaled:scale"))) {
        augZX_resu$p_bv <- augZX_resu$p_bv -sum(log(X_scale))
      } 
    } else if ( ncol(X.Re)==0L) {## ML standard
      augZX_resu$p_bv <- augZX_resu$p_v
    } else {## non-standard REML: => no X-scaling
      locXscal <- auglinmodblob$sXaug   
      if (inherits(locXscal, "AUGI0_ZX_sparsePrecision")) {
        locsXaug <- .calc_sXaug_Re_spprec(locXscal,X.Re)   
      } else {
        weight_X <- auglinmodblob$weight_X 
        nobs <- auglinmodblob$nobs
        H_global_scale <- attr(auglinmodblob$sXaug,"H_global_scale")
        w.ranef <- attr(auglinmodblob$sXaug,"w.ranef")
        if (inherits(locXscal,"Matrix")) {
          locXscal <- .Dvec_times_Matrix_lower_block(1/weight_X,locXscal,n_u_h)
          mMatrix_method <- .spaMM.data$options$Matrix_method
        } else {
          Xrows <- n_u_h+seq(nobs)
          locXscal[Xrows,] <- .Dvec_times_matrix(1/weight_X,locXscal[Xrows,]) ## get back to unweighted scaled matrix
          mMatrix_method <- .spaMM.data$options$matrix_method
        }
        locXscal <- .calc_sXaug_Re(locXscal,X.Re,rep(1,nobs))   ## non-standard REML: => no X-scaling
        locsXaug <- do.call(mMatrix_method,
                            list(Xaug=locXscal, weight_X=weight_X, w.ranef=w.ranef, H_global_scale=H_global_scale))
      }
      loc_unscaled_logdet_r22 <- get_from_MME(locsXaug,"logdet_r22") 
      augZX_resu$p_bv <- augZX_resu$p_v - loc_unscaled_logdet_r22 + ncol(X.Re)*log(2*pi)/2
    }
  }
  return(augZX_resu)
}

# .calc_APHLs_from_auglinmodblob <- function(auglinmodblob,processed, which, phi_est, lambda_est) { 
#   APHLs_args <- list(processed=processed, which=which, phi_est=phi_est, lambda_est=lambda_est)
#   APHLs_args$sXaug <- auglinmodblob$sXaug
#   APHLs_args$dvdu <- auglinmodblob$wranefblob$dvdu
#   APHLs_args$u_h <- auglinmodblob$u_h 
#   APHLs_args$mu <- auglinmodblob$muetablob$mu
#   do.call(".calc_APHLs_from_ZX", APHLs_args)[[which]]
# } 

.dsCsum <- function(A, B, keep_names=FALSE) {
  if ( any(dim(A)!=dim(B))) stop("Dimensions of the two matrices are not identical") # if unprotected, causing hard crash
  B <- .Rcpp_Csum(A,B)
  B <- forceSymmetric(B)
  if (keep_names) dimnames(B) <- dimnames(A)
  return(B)
}


