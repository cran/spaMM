## declared "arglist" to print a clean summary instead of very long list
## the print(summary) is required if the "arglist" is a member (eg as.list(<call>));
## summary alone would print nothing
print.arglist <-function(x,...) {print(summary(x,...))}
##
## GammaForDispGammaGLM function to correct a problem with the dev.resids function
#####
# in particular e.g Gamma()$dev.resids(1e10+2,1e10,1) is < 0
# then (if all y=mu, sum(dev.resids) is <0) 
# and sum(dev.resids) is used for the aic(call) in glm.fit
# where it causes a warnings for NaN in dgamma.
# We want to suppress these warnings. 
# We're not interested in this AIC anyway.
# cf also MASS p. 174...
#####

calcDhDv2 <- function(ZAL,w.resid,w.ranef) {
  d2hdv2 <- - ZtWZ(ZAL,w.resid) 
  diag(d2hdv2) <- diag(d2hdv2) - w.ranef ## small benefit that diag(w.ranef) not called on length-1 w.ranef which may occasionally be meaningful
  d2hdv2
}

`calc.w.resid` <- function(GLMweights,phi_est) {
  phi_est[phi_est<1e-12] <- 1e-11 ## 2014/09/04 local correction, cf comment in calc.p_v
  as.vector(GLMweights/phi_est)
}


GammaForDispGammaGLM <- function (link = "inverse") {
    linktemp <- substitute(link)
    if (!is.character(linktemp)) 
        linktemp <- deparse(linktemp)
    okLinks <- c("inverse", "log", "identity")
    if (linktemp %in% okLinks) 
        stats <- make.link(linktemp)
    else if (is.character(link)) 
        stats <- make.link(link)
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name)) 
                linktemp <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for gamma family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")), 
                domain = NA)
        }
    }
    variance <- function(mu) mu^2
    validmu <- function(mu) all(mu > 0)
    dev.resids <- function(y, mu, wt) {
      pmax(-2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu),1e-16)  ## FR: added the pmax
    }
    aic <- function(y, n, mu, wt, dev) {
        n <- sum(wt)
        disp <- dev/n
        -2 * sum(dgamma(y, 1/disp, scale = mu * disp, log = TRUE) * 
            wt) + 2
    }
    initialize <- expression({
        if (any(y <= 0)) stop("non-positive values not allowed for the gamma family") 
        n <- rep.int(1, nobs)
        mustart <- y
    })
    simfun <- function(object, nsim) {
        wts <- object$prior.weights
        if (any(wts != 1)) 
            message("using weights as shape parameters")
        ftd <- fitted(object)
        shape <- MASS::gamma.shape(object)$alpha * wts 
        rgamma(nsim * length(ftd), shape = shape, rate = shape/ftd)
    }
    structure(list(family = "Gamma", link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
        validmu = validmu, valideta = stats$valideta, simulate = simfun), 
        class = "family")
}

dispGammaGLM <- function(dev.res,lev,X,offset=NULL,family=GammaForDispGammaGLM(link=log),method="glm") {
  ## do not 'filter' the dev.res and lev (in any way different from the lamScal one) here otherwise the lambda estimate may be inconsistent with the v_h
  ## except that fatally, a 0L response value occurs in normal use; eg Poisson, y ~ 1+ ranef(lambda->very small), mean(y) happens to be equal to some y value -> u_h =0L
  resp <- dev.res/((1-lev)) 
  resp[resp==0L] <- 1e-100
  resp[resp>1e150] <- 1e150 ## v1.2 fix for glm -> glm.fit -> .Call(C_Cdqrls, x[good, , drop = FALSE]*w problem
  weight <- (1-lev)/2 
  etastart <- rep(family$linkfun(mean(resp)),nrow(resp))   ## glm needs a bit help...
  if (is.null(offset)) offset <- rep.int(0, nrow(resp))
  etastart <- etastart - offset
  if (method=="glm"){
    #if (interactive()) {
    #  resglm <- glm(resp~X-1,family=family,weights = weight,etastart=etastart,,offset=offset)
    # warnmess <- NULL
    #} else {
      resglm_wrap <- tryCatch.W.E(glm(resp~X-1,family=family,weights = weight,etastart=etastart,offset=offset))    
      resglm <- resglm_wrap$value
      warnmess <- resglm_wrap$warning$message ## may be NULL
    #}
    # if (verbose["trace"]) {print(paste("phi coefficients=",paste(signif(resglm$coefficients,4),collapse=", ")),quote=F)}
    Qr <- resglm$qr  
    beta_disp <- resglm$coefficients[Qr$pivot[1L:ncol(X)]] ## As in summary.glm.
  } else {
    ## this will remain similar to glm aslong as HLM -> provide.resglm -> glm.fit 
    stop("need to put back HLM.R into the sources")
    #resglm <- HLM(resp~X-1,family=family,prior.weights = weight,offset=offset) 
    beta_disp <- resglm$fixef
    summ <- NULL
    warnmess <- NULL
  }
  return(list(beta_disp=beta_disp,next_disp_est=fitted(resglm),resglm=resglm,warnmess=warnmess))
}



`selectLoglfn` <- function(family) {
   family <- tolower(family)
   switch(family,
      gaussian = function(theta,y,nu) {nu*(theta*y-(theta^2)/2)- ((y^2)*nu+log(2*pi/nu))/2}, 
      poisson = function(theta,y,nu) {
        res <- nu*(theta*y-exp(theta))   -  lfactorial(y)
        res[theta== -Inf & y==0] <- 1
        res
        },
      binomial = function(theta,freqs,sizes,nu) {nu*sizes*(freqs*theta-log(1+exp(theta))) +lchoose(sizes,round(sizes*freqs))},
      # gamma = function(theta,y,nu) {nu*(y*theta+log(-theta))+nu*(log(nu*y))-lgamma(nu)-log(y)} ## mean mu=-1/th, **** var = mu^2 / vu ****
      # same bu using ad hoc C library...
      gamma = function(theta,y,nu) {
        disp <- 1/nu
        mu <- -1/theta
        dgamma(y, shape=nu , scale = -1/(nu*theta), log = TRUE) ## from Gamma(log)$aic
      }
    )
}


`theta.mu.canonical` <-function(mu,family) { 
   ## the (fixed) canonical link between theta and mu, not the family link between eta and mu 
   switch(tolower(family),
      gaussian = mu ,
      poisson = log(mu) ,
      binomial = make.link("logit")$linkfun(mu),  # correct syntax, does not use 'non-public API' such as .Call to access code from dlls from the base packages...
## if this does no work, use 
#                 { 
#                    theta <- logit(mu)
#                    theta[theta>27.6310] <- 27.6310 ## mu>1-1e-12
#                    theta[theta < -27.6310] <- -27.6310 
#                    theta
#                 },
      gamma = -1/mu, ## "-": McC&N p. 290
   )
} ## returns values for given mu

`theta.mu.conjugate` <-function(mu,family) { 
   ## theta(u) in LeeN01... this is more pedagogy than efficient code
   switch(tolower(family),
      gaussian = theta.mu.canonical(mu,"gaussian") , ## mu 
      gamma = theta.mu.canonical(mu,"poisson"), ## log(mu)
      beta = theta.mu.canonical(mu,"binomial"), ## improved logit(mu)      
      "inverse.gamma" = theta.mu.canonical(mu,"gamma") ## -1/mu
   )
} ## returns values for given mu

## another approach is
#theta.mu.call <- switch(tolower(family),
#      poisson =  call("log",quote(mu)) ,
#   )
## followed by eval(theta.mu.call)

### there is a nice syntax
#eta.mu.expr <- parse(text=paste(family$link,"(mu)",sep="")) # "log" -> expression(log(mu))
#etadmu<-D(theta.mu.expr,"mu") ## a call; eval(dthetadmu) then suffices if mu has a value
### except that D(identity... does not work)

thetaMuDerivs <-function(mu,BinomialDen,familyfam) {
  if (familyfam=="binomial") muFREQS <- mu/BinomialDen
  ## these definitions depend only on the canonical link
  Dtheta.Dmu <- switch(tolower(familyfam),
    gaussian = rep(1,length(mu)) ,
    poisson = 1/mu ,
    binomial = 1/(muFREQS*(1-muFREQS)),
    gamma = 1/mu^2
  ) ## values for given mu
  if (familyfam=="binomial") Dtheta.Dmu <- Dtheta.Dmu/BinomialDen
  D2theta.Dmu2 <- switch(tolower(familyfam),
    gaussian = rep(0,length(mu)) ,
    poisson = -1/mu^2 ,
    binomial = -(1-2*muFREQS)/(muFREQS*(1-muFREQS))^2,
    gamma = -2/mu^3
  ) ## values for given mu
  if (familyfam=="binomial") D2theta.Dmu2 <- D2theta.Dmu2/(BinomialDen^2)
  return(list(Dtheta.Dmu=Dtheta.Dmu,D2theta.Dmu2=D2theta.Dmu2))
}

muetafn <- function(family,eta,BinomialDen,priorWeights=1) { ## note outer var BinomialDen 
  ## a patch for large eta in poisson case
  if (family$family=="poisson" && family$link =="log") {
    eta <- pmin(30,eta) ## 100 -> mu = 2.688117e+43 ; 30 -> 1.068647e+13
  }
  mu <- family$linkinv(eta) ## linkinv(eta) is FREQS for binomial, COUNTS for poisson...
  if (family$link %in% c("logit","probit","cloglog")) {
        mu[mu > (1-1e-12)] <- (1-1e-12)
        mu[mu < (1e-12)] <- (1e-12)
  }
  dmudeta <- family$mu.eta(eta) ## aberrant at hoc code for cloglog 'elsewhere'...
  Vmu <- family$variance(mu)
  if (family$family=="binomial") {
      Vmu <- Vmu * BinomialDen ## not checked for probit et cloglog
      mu <- mu *BinomialDen
      dmudeta <- dmudeta * BinomialDen
  } 
  GLMweights <- priorWeights * dmudeta^2 /Vmu ## must be O(n) in binomial cases
  if (any(is.nan(GLMweights))) {
#    calls <- sys.calls()
#    ncalls <- length(calls)
#    HLfitcall <- calls[[ncalls-1]]
#    save(HLfitcall,file=generateFileName("HLfitcall")) ## debug  code, not for package
    stop("NaN GLMweights generated in 'muetafn'")
  }
return(list(mu=mu,dmudeta=dmudeta,Vmu=Vmu,GLMweights=GLMweights))
} ## end local def muetafn

updateWranef <- function(rand.family,lambda,u_h,v_h) {
  dudv <- rand.family$mu.eta(v_h) ## general cf InvGamma with log link rand.family$mu.eta(v_h) = exp(v) =u is du/d(log u)   
#  dudv <- pmin(.Machine$double.xmax,dudv) ## patch necess for poisson gamma
  ## compute w.ranef, the scaled d^2f(v)/dv^2. See Appendix 3 of LeeN01 + my notes for what it implies...
  if (rand.family$family=="gaussian") {
    if (rand.family$link=="identity") {
      V_M <- rand.family$variance(u_h) ##rep(1,length(u_h)) ## GLMMs in general
      dlogWran_dv_h <- rep(0L,length(u_h))
    }
  } else if (rand.family$family=="Gamma") { 
    if (rand.family$link=="log") {
      V_M <- u_h ## V(u), canonical conjugate Gamma as in canonical Poisson Gamma HGLM
      dlogWran_dv_h <- rep(1L,length(u_h))
    } 
  } else if (rand.family$family=="inverse.Gamma") { ## for Gamma HGLM 
    ## the canonical form gives the density of theta(u)
    if (rand.family$link=="log") {
      w.ranef <- as.numeric(1/(u_h * lambda)) ## W1/lambda, W1 computation shown in appendix 3 of LeeN01; also in Noh and Lee's code.
      dlogWran_dv_h <- rep(-1L,length(u_h)) ## v=log u, dlogW/dv= dlog(1/u)/dv=-1
      return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))  ###### return here !
    } else if (rand.family$link=="-1/mu") {
      ## D[\[Nu] (\[Theta][u] - (-Log[-\[Theta][u]])), {\[Theta][u], 2}]
      V_M <- rand.family$variance(u_h) ## u_h^2 ## V(u), canonical conjugate HGLM 
      dlogWran_dv_h <- 2 * u_h ## no independent check 
    }
  } else if (rand.family$family=="Beta") {
    if (rand.family$link=="logit") {
      V_M <- rand.family$variance(u_h) ##  u_h*(1-u_h) ## canonical conjugate HGLM
      dlogWran_dv_h <- 1 - 2 * u_h ## D[Log[u (1 - u)] /. u -> 1/(1 + E^-v), v] /. v -> Log[u/(1 - u)] ; no independent check
    }
  }
#  V_M <- pmin(.Machine$double.xmax,V_M) ## patch necess for poisson gamma
  ## dudv/V_M may be 1 as both diverge: 
  w.ranef <- as.numeric((dudv/V_M)*(dudv/lambda)) ## semble valide quand v=g(u) = th(u): not previous return()
#if (any(is.nan(w.ranef))) browser()
#  w.ranef <- pmin(.Machine$double.xmax,w.ranef) ## patch necess for poisson gamma
  return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=1/dudv))
}

updateW_ranefS <- function(cum_n_u_h,rand.families,lambda,u_h,v_h) {
  nrand <- length(rand.families)
  blob <- lapply(seq(nrand), function(it) {
    u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
    updateWranef(rand.family=rand.families[[it]],lambda[u.range],u_h[u.range],v_h[u.range])
  })
  w.ranef <- unlist(lapply(blob,function(b) {b$w.ranef}))
  w.ranef <- pmin(w.ranef,1e10) ## patch useful to avoid singular d2hdv2 in PLoG model
  dlogWran_dv_h <- unlist(lapply(blob,function(b) {b$dlogWran_dv_h}))
  dvdu <- unlist(lapply(blob,function(b) {b$dvdu}))
  return(list(w.ranef=w.ranef,dlogWran_dv_h=dlogWran_dv_h,dvdu=dvdu))
}

d2mudeta2fn <- function(link,mu=NULL,eta=NULL,BinomialDen=NULL) { ## d2 MuCOUNTS d etaFREQS^2
  switch(link,
      identity = 0,
      log = mu, 
      inverse = 2 * mu^3 , ## canonical for Gamma()
      ## next three make sense for Binomial response data
      logit = {muFREQS <- mu/BinomialDen;
               d2muFREQS <- muFREQS*(1-muFREQS)*(1-2*muFREQS);
               d2muFREQS * BinomialDen
               },
      probit = -eta*dnorm(eta) * BinomialDen,
      cloglog = exp(eta-exp(eta))*(1-exp(eta)) * BinomialDen ## D[1 - E^-E^\[Eta], {\[Eta], 2}]
  )
} 

`safesolve.qr.matrix` <- function(qr.a,B,silent=TRUE,stop.on.error=TRUE) { ## solve.qr with fall-back; qr.a should be a qr object, B a matrix
  ## there was a 'Matrix' subcode prior to 10/03/2013
    res <- try(solve.qr(qr.a,B),silent=silent)
    if (class(res)=="try-error") { ##FR->FR sb systematique qd phi -> 0; slow step.
      # pivI <- sort.list(qr.a$pivot)  ## inverse perm such as pivI[$pivot]=$pivot[pivI]= identity
      #solveA <- try(solve(qr.R(qr.a)[,pivI])%*%t(qr.Q(qr.a)),silent=silent)  
      ## solveA is solve(<original 'a' matrix>) using the QR decomp... but this may work when solve.qr fails !
      solveA <- try(backsolve(qr.R(qr.a),t(qr.Q(qr.a))[qr.a$pivot,]))
      if (class(solveA)=="try-error") {
        if (stop.on.error) {
          mess <- pastefrom("class(solveA)='try-error'.",prefix="(!) From ")
          message(mess)
          stop("More code is needed to handle this... I exit.") ## perhaps recover A by qr.X and solve(A) ?
        } else return(solveA) ## passes control to calling function
      } else res <- solveA %*% B  
    }
  return(res)
}

`safesolve.qr.vector` <- function(qr.a,b,silent=TRUE,stop.on.error=TRUE) { ## solve.qr with fall-back; qr.a should be a qr object, b must be a vector
  if (class(qr.a)=="sparseQR") { ## pas de 'essai' en var locale !
    ## there was a 'Matrix' subcode prior to 10/03/2013; another try on 11/2013
    res <- qr.coef(qr.a,b)
  } else {
    res <- try(solve.qr(qr.a,b),silent=silent)
    if (class(res)=="try-error") {   ## then some weird code, but...
      ## we try to solve(<original 'a' matrix>) using the QR decomp... this may work when solve.qr fails !
      ## The following is equivalent to solveA <- try(solve(qr.R(qr.a)[,pivI])%*%t(qr.Q(qr.a)),silent=silent) then return solveA %*% b 
      ## but uses only backward mutliplication of vectors and transpose of vectors  
      res <- t(t(b) %*% (qr.Q(qr.a))) ## not yet res
      #pivI <- sort.list(qr.a$pivot) ## inverse perm such as pivI[$pivot]=$pivot[pivI]= identity
      #res <- try(solve(qr.R(qr.a)[, pivI]) %*% res, silent=silent)
      res <- try(backsolve(qr.R(qr.a),res[qr.a$pivot]))
      if (class(res)=="try-error") {
        if (stop.on.error) {
          mess <- pastefrom("class(res)='try-error'.",prefix="(!) From ")
          message(mess)
          stop("More code is needed to handle this... I exit.") ## perhaps recover A by qr.X and solve(A) ?
        } else return(res) ## passes control to calling function
      } 
    }  
  }
  return(res)
}

solveWrap.vector <- function(A,b,...) {
  if (inherits(A,"diagonalMatrix")) return(b/diag(A))
  if (inherits(A,"Rcpp-QR")) { ## no pivoting, $R is upper triangular
    b <- t(t(b) %*% (A$Q))
    solved <- try(backsolve(A$R,b),silent=TRUE)
    return(solved) ## gives control to calling function 
  }
  ## all other cases
  safesolve.qr.vector(A,b,...)
}

solveWrap.matrix <- function(A,B,...) {
  if (inherits(A,"diagonalMatrix")) return(B/diag(A)) ## works if A is the matrix, not its diagonal...
  if (inherits(A,"Rcpp-QR")) { ## no pivoting, $R is upper triangular
    solved <- try(backsolve(A$R,t(A$Q) %*%B),silent=TRUE)
    return(solved) ## gives control to calling function 
  }
  ## all other cases
  safesolve.qr.matrix(A,B,...)
}

QRwrap <- function(mat) {
  if (.spaMM.data$options$USEEIGEN) {
    QR <- Rcpp_QR(mat) ## the main benefit is that the solve method for Rcpp-QR objects is numerically more stable 
  } else {
    QR <-qr(mat)
  }
  #   ## ready debugging code: try(...,silent=TRUE) + 
  #   if (class(QR)=="try-error") {
  #     mess <- pastefrom("problem in QR computation.",prefix="(!) From ")
  #     stop(mess)
  #   }
  return(QR)
} 

LogAbsDetWrap <- function(mat) {
  if (.spaMM.data$options$USEEIGEN) {
    lad <- LogAbsDetCpp(mat)
  } else lad <-determinant(mat)$modulus[1] 
  return(lad)
}

## this serves mainly for the SE of beta (but also in betaFirst)
## One needs Xt.InvS.X, a small matrix but InvS is the slow step.  
`calc.tXinvS` <- function(Sig,X.pv,stop.on.error,lambda_est,ranFix) { ## slow... 
  #  qr.Sig <- qr(Sig) ## SLOW
  #  XinvS <- safesolve.qr.matrix(qr.Sig,X.pv,stop.on.error=stop.on.error) ## invSig %*% X.pv
  qr.Sig <- Rcpp_QR(Sig) ## FR->FR still SLOW
  XinvS <- solveWrap.matrix(qr.Sig,X.pv,stop.on.error=stop.on.error) ## invSig %*% X.pv
  if (class(XinvS)=="try-error") {
    if (stop.on.error) {
      mess <- pastefrom("the augmented 'Sigma' matrix appears singular. Extreme lambda/phi value and/or extremely correlated random effects?",prefix="(!) From ")
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
    } else {
      return(XinvS) ## returns a try-error
    }
  } else tXinvS <- t(XinvS) ## we have to transpose either this one or X.pv, which are of the same size 
  return(tXinvS)
}

## for SEM

rntneg <- function(n,mu,sigma2)
{
  # produce n samples from the
  # specified rigth-truncated to 0 gaussian
  # distribution
  pn <- runif(n)*pnorm(0,mu,sqrt(sigma2))
  pn[pn==0] <-  .Machine$double.eps ## because if mu is large -> qnorm(0) is -Inf which later cause NaN's
  pn[pn==1] <-  1-.Machine$double.eps 
  qnorm(pn,mu,sqrt(sigma2))
  ## alternatively use ... qn[pn==0] <- ... 
}

################################################################################

rntpos <- function(n,mu,sigma2)
{
  # produce n samples from the
  # specified left-truncated to 0 gaussian
  # distribution
  -rntneg(n,-mu,sigma2)
}



#rntpos <- function(n,mu,sigma2) {rtnorm.copy(n, mean = mu, sd = sqrt(sigma2), lower = 0, upper = Inf)}
#rntneg <- function(n,mu,sigma2) {rtnorm.copy(n, mean = mu, sd = sqrt(sigma2), lower = -Inf, upper = 0)}


## returns a list !!
## input LMatrix is either a single LMatrix whcih is assumed to be the spatial one, or a list of matrices 
`compute.ZALlist` <-function(LMatrix=NULL,CMatrix=NULL,ZAlist,Groupings) {## we need to check for user's confusion when we multiply Z by LMatrix
  ## ZAL is nobs * (# levels ranef) and ZA too
  ## LMatrix is (# levels ranef) * (# levels ranef) [! or more generally a list of matrices!]
  ## CMatrix is an alternative to LMatrix that requires that the return value is used together with coefficients [t(L_ori)]^{-1} v_ori
  ## the levels of the ranef must match each other in the two matrices
  ## the only way to check this is to have the levels as rownames and colnames and to check these
  if (is.null(ZAlist)) return(list())
  ## ELSE
  ZAL <- ZAlist
  if (is.null(LMatrix)) LMatrix <- CMatrix ## the two inputs are not further distinguished below (but imply different meanings for the results)
  if ( ! is.null(LMatrix) && length(ZAlist)>0 ) {
    if (inherits(LMatrix,"blockDiag")) {
      stop("compute.ZALlist code should be revised to handle blockDiag objects")
    } ## ELSE
    if ( ! is.list(LMatrix)) LMatrix <- list(LMatrix)
    LMlen <- length(LMatrix)
    for (ii in seq_len(LMlen)) {
      lmatrix <- LMatrix[[ii]]
      ## find ZAlist elements affected by LMatrix element
      affecteds <- which(attr(ZAlist,"ranefs") %in% attr(lmatrix,"ranefs"))
      for (it in affecteds) {
        ZA <- ZAlist[[it]]
        if (inherits(ZA,"identityMatrix")) {
          ZAL[[it]] <- lmatrix          
        } else {
          locnc <- ncol(ZA)
          locnr <- nrow(lmatrix)
          if ( locnc %% locnr !=0) {
            mess <- paste("The number of levels of the grouping variable in random term (...|",Groupings[it],")",sep="")
            mess <- paste(mess,"\n  is not the dimension of the correlation matrix.") ## by distMatrix checking in corrHLfit or no.info check somewhere...
            stop(paste(mess," I exit."))
          }         
          nblocks <- locnc %/% locnr 
          if (nblocks>1) {
            locZA <- ZA
            for (bt in 1:nblocks) 
              locZA[,locnr*(bt-1)+(1:locnr)] <- locZA[,locnr*(bt-1)+(1:locnr)] %*% lmatrix
            ZAL[[it]] <- locZA
          } else {
            ## rownames(LMatrix) are the names of first occurrences of unique geographic locations, 
            ## or (if user provided distMatrix) whatever was in this distMatrix. But with a distMatrix, it is likely that ZA was = I and we don't reach this code
            ### it's difficult to make checks on names at this step
            ## LMatrix inherits its names from thos of uniqueGeo. These are the names of first occurrences of unique geographic locations, or (if user provided distMatrix) whatever was in this distMatrix
            ## ZAlist inherits anything from the spMMFactorList call which input does not include info about rownames of data
            #             if ( ! all(attr(ZA,"colnames")==rownames(lmatrix))) {
            #               stop("The colnames of the design matrix Z in eta=...+Zv should be the rownames of the design matrix L  in v=Lu")
            #             }
            ZAL[[it]] <- ZA %*% lmatrix
          }
        }
        attr(ZAL[[it]],"userLfixed") <- attr(lmatrix,"userLfixed") ## TRUE or NULL
      }
    }
  }
  attr(ZAL,"userLfixeds") <- unlist(lapply(ZAL,function(mat) { 
    att <- attr(mat,"userLfixed") ## TRUE or NULL   
    if (is.null(att)) att <- FALSE
    att
  })) ## vector of TRUE or FALSE
  return(ZAL)
}

post.process.ZALlist <- function(ZAL,predictor) {
  nrand <- length(ZAL)
  if (nrand>1) {
    ZAL <- do.call(cbind,ZAL)
    # FR->FR ici ça serait bien de pouvoir identifier un diagonal matrix... (? pq ?)
  } else if (nrand==1) {
    ZAL <- ZAL[[1]] 
    if ((! is.null(attr(predictor,"%in%"))) && attr(predictor,"%in%") && ncol(ZAL)==nrow(ZAL)) { 
      ## test of the attribute is a heuristic way of detecting when using the block structure will lead to faster analysis
      partition <- findblocks(ZAL) 
      if ( length(partition)>1 ) {
        partition <- cumsum(c(0,partition))
        attr(ZAL,"partition") <- partition
      }
    }
  }
  I <- diag(rep(1,ncol(ZAL)))
  attr(ZAL,"ZALI") <- rbind(ZAL,I) 
  return(ZAL)
}


intervalStep <- function(old_betaV,wAugX,wAugz,currentp_v,intervalInfo,corrPars) {
  parmcol <- attr(intervalInfo$parm,"col")
  #print((control.HLfit$intervalInfo$fitp_v-currentp_v)/(control.HLfit$intervalInfo$MLparm-old_betaV[parmcol]))
  ## voir code avant 18/10/2014 pour une implem rustique de VenzonM pour debugage  
  ## somewhat more robust algo (FR->FR: still improvable ?), updates according to a quadratic form of lik near max
  ## then target.dX = (current.dX)*sqrt(target.dY/current.dY) where dX,dY are relative to the ML x and y 
  ## A nice thing of this conception is that if the target lik cannot be approached, 
  ##   the inferred x converges to the ML x => this x won't be recognized as a CI bound (important for locoptim) 
  currentDx <- (old_betaV[parmcol]-intervalInfo$MLparm)
  targetDy <- (intervalInfo$fitp_v-intervalInfo$targetp_v)
  currentDy <- (intervalInfo$fitp_v-currentp_v)
  if (currentDy <0) {
    message("An higher likelihood was found than for the original fit.\nThis suggests the original fit did not fully maximize the likelihood.")
    if (length(corrPars)>0) message(paste("Current correlation parameters are ",
                                          paste(names(corrPars),"=",signif(unlist(corrPars),6),collapse=", ")))
    message("Current likelihood is p_v=",currentp_v)                    
  } else {
    betaV <- rep(NA,length(old_betaV))
    betaV[parmcol] <- intervalInfo$MLparm + currentDx*sqrt(targetDy/currentDy)
  }
  betaVQ <- lmwithQ(wAugX[,-parmcol,drop=FALSE],wAugz-wAugX[,parmcol]*betaV[parmcol])
  betaV[-parmcol] <- betaVQ$coef
  return(list(betaVQ=betaVQ,betaV=betaV))
}

## cette fonction marche que si on a fixed effect + un terme aleatoire....
eval.corrEst.args <- function(family,rand.families,predictor,data,X.Re,
                              distinct.X.ReML,REMLformula,ranFix,
                              term=NULL,
                              Optimizer) {
  ## ici on veut une procedure iterative sur les params de covariance
  #  HLCor.args$processed <- processed ## FR->FR dangerous in early development
  corrEst.args <- list(family=family,rand.family=rand.families) ## but rand.families must only involve a single spatial effect 
  loc.oriform <- attr(predictor,"oriFormula")
  loc.lhs <- paste(loc.oriform)[[2]]
  ## build formula, by default with only spatial effects
  if (is.null(term)) term <- findSpatial(loc.oriform)
  corrEst.form <-  as.formula(paste(loc.lhs," ~ ",paste(term)))
  corrEst.args$data <- data ## FR->FR way to use preprocess ???                    
  if (ncol(X.Re)>0) { ## some REML correction
    if (distinct.X.ReML) {
      corrEst.args$REMLformula <- REMLformula
    } else corrEst.args$REMLformula <- predictor ## REML without an explicit formula
    corrEst.args$objective <- "p_bv"
  } else corrEst.args$objective <- "p_v" 
  corrEst.args$ranFix <- ranFix ## maybe not very useful
  corrEst.args$control.corrHLfit$Optimizer<- Optimizer ## (may be NULL => L-BFGS-B) 
  corrEst.args$control.corrHLfit$optim$control$maxit <- 1 
  corrEst.args$control.corrHLfit$nlminb$control$iter.max <- 2 ## 1 => convergence trop lente
  corrEst.args$control.corrHLfit$optimize$tol <- 1e10 
  corrEst.args$control.corrHLfit$corners <- FALSE ## 
  return(list(corrEst.args=corrEst.args,corrEst.form=corrEst.form))
}




corr.notEQL.lambda <- function(nrand,cum_n_u_h,lambda_est,lcrandfamfam) {  
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
  return(notEQL)
}

u_h_from_v_h <- function(v_h,rand.families,cum_n_u_h,lcrandfamfam) {
  anyinf <- FALSE
  nrand <- length(rand.families)
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





## glm convention in binomial models : eta, fitted values describes FREQUENCIES
##                                     linkinv(eta) describes frequencies, but we need mu to scale as y in the code...
## but the input response ar COUNTS
HLfit <- function(formula,
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
  family <- checkRespFam(family)
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
  ################# family multi  #########################################
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
  #corrFixNames <- intersect(names(ranFix),c("nu","rho","Nugget","ARphi"))
  #corrPars <- ranFix[corrFixNames] ## as for functions in corrMM.LRT that always look in phi, lambda, rather than .Fix. 
  # corrNames <- intersect(c("nu","rho","Nugget","ARphi"),names(init.HLfit)) ## the ones optimized within HLfit (confusing name)
  #if (length(corrNames)>0) {
  #  corr_est <- init.HLfit[corrNames]
  #} else corr_est <- NULL
  corrNames_in_ranFix <- intersect(names(ranFix),c("nu","rho","Nugget","ARphi"))
  corrNames_in_init_HLfit <- intersect(c("nu","rho","Nugget","ARphi"),names(init.HLfit)) ## the ones optimized within HLfit 
  if (length(corrNames_in_init_HLfit)>0) {
    corr_est <- init.HLfit[corrNames_in_init_HLfit]
  } else corr_est <- NULL
  ## corrPars is only for info in messages() and return value, 
  corrPars <- ranFix[corrNames_in_ranFix] ## as for functions in corrMM.LRT that always look in phi, lambda, rather than .Fix. 
  corrPars[corrNames_in_init_HLfit] <- NA ## will be filled at the end of the fit
  typelist <- list()
  typelist[corrNames_in_ranFix] <- "fix"
  if (!is.null(rFtype <- attr(ranFix,"type"))) { 
    corrNames_in_ranFix_type <- intersect(corrNames_in_ranFix,names(rFtype))
    typelist[corrNames_in_ranFix_type] <- rFtype[corrNames_in_ranFix_type]
  }
  typelist[corrNames_in_init_HLfit] <- "var"
  attr(corrPars,"type") <- typelist
  ###################################################
  warningList<-list()
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["warn"])) verbose["warn"] <- TRUE
  if (is.na(verbose["summary"])) verbose["summary"] <- FALSE
  ##
  if (family$family %in% c("poisson","binomial")) {
    phi.Fix<-1 
  } else {
    phi.Fix <- getPar(ranFix,"phi")
    if (any(phi.Fix==0)) {
      mess <- pastefrom("phi cannot be fixed to 0.",prefix="(!) From ")
      stop(mess)
    }
  } ## immediately used in preprocess call:
  if (is.null(processed)) {
    validdata <- validData(formula=formula,resid.formula=resid.formula,data=data) ## will remove rows with NA's in required variables
    if (!inherits(data,"environment")) {
      data <- data[rownames(validdata),,drop=FALSE] ##     before Predictor is called and an LMatrix is added, etc. 
    } else data <- validdata
    loclist <- list(control.HLfit=control.HLfit,HLmethod=HLmethod,predictor=formula,phi.Fix=phi.Fix,
                    resid.predictor=resid.formula,REMLformula=REMLformula,data=data,family=family,
                    rand.families=rand.family) ## BinomialDen always missing here
    processed <- do.call("preprocess",loclist)
  } 
  predictor <- processed$predictor
  rand.families <- processed$rand.families
  lcrandfamfam <- attr(rand.families,"lcrandfamfam")
  HL <- processed$HL
  if (HL[1]=="SEM") {
    SEMargs <- processed[names(processed) %in% c("SEMseed","nMCint","nSEMiter","ngibbs","SAEMsample")]
  }
  stop.on.error <- processed$stop.on.error ## to control issues with matrix computations; F by default
  AIC <- processed$AIC ## whether to compute any AIC stuff; F by default
  essai <- processed$essai ## to control any tested new code...
  conv.threshold <- processed$conv.threshold
  iter.mean.dispFix <- processed$iter.mean.dispFix
  iter.mean.dispVar <- processed$iter.mean.dispVar
  max.iter <- processed$max.iter
  #    ps_v.threshold <- processed$ps_v.threshold
  resid.predictor <- processed$resid.predictor 
  BinomialDen <- processed$BinomialDen
  y <- processed$y
  REMLformula <- processed$REMLformula
  X.Re <- processed$`X.Re`
  X.pv <- processed$`X.pv`
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
    } else {
      ZALlist <- list(`1`=ZAL) ## 12/10/2014
      attr(ZALlist,"userLfixeds") <- TRUE 
    }
    ZAL <- post.process.ZALlist(ZALlist,predictor=predictor)
  } else { ## models[["eta"]] = "etaGLM"
    ZALlist <- NULL
    ZAL <- NULL
    u_h <- v_h <- lev_lambda <- numeric(0)
  } 
  ### a bit of post processing // repeat of code in preprocess...
  nrand <- length(ZALlist)
  lambda.Fix <- getPar(ranFix,"lambda")
  if (any(lambda.Fix==0)) {
    mess <- pastefrom("lambda cannot be fixed to 0.",prefix="(!) From ")
    stop(mess)
  }
  vec_n_u_h <- rep(0, nrand)
  for (i in 1:nrand) vec_n_u_h[i] <- ncol(ZALlist[[i]]) ## nb cols each design matrix = nb realizations each ranef
  cum_n_u_h <- cumsum(c(0, vec_n_u_h)) ## if two ranef,  with q=(3,3), this is 0,3,6. cum_n_u_h[nrand+1] is then 6, the total # of realizations
  ###
  X_lamres <- processed$X_lamres
  next_cov12_est <- NULL ## will be tested
  X_disp <- processed$X_disp ## may be NULL
  p_phi <- ncol(X_disp) 
  off <- attr(processed$predictor,"offsetObj")$vector
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
    corrEstBlob <- eval.corrEst.args(family=family,rand.families=rand.families,predictor=predictor,data=data,X.Re=X.Re,
                                      distinct.X.ReML=distinct.X.ReML,REMLformula=REMLformula,ranFix=ranFix,
                                      Optimizer=control.HLfit$Optimizer)
    corrEst.args <- corrEstBlob$corrEst.args ## but corrEstBlob also has $corrEst.form which will stay there for later use
  }
  #################### MORE LOCAL FNS DEFS ###################################
  ##mais processed controle le default nMCint
  ## all per-iteration stats are taken from gibbsSample
  ## and all final stats are the means, from iterations SAEMsample, of the per-iteration stats 
  resize.lambda <- function(lambda) {
    if  (length(lambda)==nrand) {
      lambda_est <- rep(lambda,vec_n_u_h)
    } else if (length(lambda)==1) { ## typically what the current default resglm provides even for nrand>1
      lambda_est <- rep(lambda,cum_n_u_h[nrand+1L])
    } else if (length(lambda)==cum_n_u_h[nrand+1L]) {
      lambda_est <- lambda
    } else {stop("Initial lambda cannot be mapped to levels of the random effect(s).")}
    lambda_est
  }
  
  SEMbetalambda <- function(beta_eta,nSEMiter=100,ngibbs=20,nMCint=10000,SEMseed=NULL,SAEMsample=NULL){ ## beta_eta as explicit argument so that the iterative aspect is explicit
    if (nSEMiter<10) stop("(!) In 'SEMbetalambda', 'nSEMiter' should be >9")
    if (is.null(SAEMsample)) SAEMsample <- (nSEMiter/2):nSEMiter
    if(!is.null(SEMseed)) {
      set.seed(SEMseed) ## so that estimates of beta,lambda are repeatable ## comment ne pas avoir a retirer les memes nombres XXXX fois ?
    #      cat("SEMseed not null")
    } #else {cat("NULL SEMseed")}
    betaMat <- matrix(0,nrow=nSEMiter,ncol=ncol(X.pv))
    colnames(betaMat) <- colnames(X.pv) 
    EuGivenY=matrix(0,nrow=nSEMiter,ncol=length(y))
    lambdaVec <- numeric(nSEMiter)
    betaMat[1,] <- beta_eta
    lambdaVec[1] <- unique(lambda_est)
    condVar <- rep(0,nSEMiter)
    condVar[1] <- lambdaVec[1]
    ZA <- ZAlist[[1]] ## FR->FR ad hoc
    if ( ! inherits(ZA,"identityMatrix")) {
      stop("! inherits(ZA,'identityMatrix'): more code needed in SEM algo") ## CondNormf not adequate
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
      ranefCorr <- tcrossprodCpp(LMatrix) 
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
      EuGivenY[i,] <- colMeans(condMeans[gibbsSample,,drop=FALSE])
      # M part of the SEM algorithm
      # determination of beta by standard least square #betaMat[i,] <- lm((z-Lvs)~X.pv-1)$coeff
      betaMat[i,] <- solveWrap.vector( qr.XtX , t(X.pv) %*% (augY-EuGivenY[i,]-off) ,stop.on.error=stop.on.error) 
      if (is.null(lambda.Fix)) {
        lambdaVec[i] <- mean(lambdas[gibbsSample])
      } else lambdaVec[i] <- lambda.Fix
    }  ## end nSEMiter loop
    beta_eta <- colMeans(betaMat[SAEMsample,,drop=FALSE]) ## SAEMsample no longer useful ?
    lambda <- mean(lambdaVec[SAEMsample]) ## lambda_est will be given a different length
    v_h <- colMeans(EuGivenY[SAEMsample,,drop=FALSE])
    #browser()
    ## simulate final likelihood(with high variance...)
    fix <- X.pv %*% beta_eta + off
    binLikcond <- numeric(nrow(ZAL))
    logtotLikcond <- numeric(nMCint)
    if (nMCint==0) { ## use standard Laplace approx for estimating the likelihood
      ## HACK
      arglist <- as.list(mc) ## of which [[1]] is "HLfit"
      proc <- arglist$processed 
      if (! is.null(proc)) arglist$processed <- eval.update.call(proc$callargs,HLmethod="PQL/L") 
      arglist$HLmethod <- NULL ## for clarity; should be ignored anyway
      arglist$etaFix$beta <- beta_eta
      #      arglist$etaFix$v_h <- v_h ## interestingly, disastrous.
      arglist$ranFix$lambda <- lambda
      logLapp <- eval(as.call(arglist))$APHLs$p_v
      attr(logLapp,"method") <- "p_v(h) (marginal L):"
    } else {
      Lik <- 0
      for (it in 1:nMCint) {
        eta <- fix + sqrt(lambda) * (ZAL %*% rnorm(n_u_h,0))
        mu <- family$linkinv(eta) ## ou pnorm()
        binLikcond <- (1-mu)* (whichy0)+ mu*(whichy1)
        logtotLikcond[it] <- sum(log(binLikcond))
      }
      maxlogtotLikcond <- max(logtotLikcond)
      rel <- exp(logtotLikcond-maxlogtotLikcond)
      logLapp <- log(mean(rel))+maxlogtotLikcond
      seInt <- sqrt(var(log(colMeans(matrix(rel,ncol=50)))/50))
      attr(logLapp,"method") <- "  logL (MC estimate)" ## directly usable for screen output
    } 
    # if (lambda>20) {
    #       cat("\a\a\a")
    #       browser()
    # }
    return(list(beta_eta=beta_eta,lambda=lambda,v_h=v_h,logLapp=logLapp,seInt=seInt))
  } ## end local def of SEMbetalambda
      


  

  makeCovEst <- function(u_h,ZAlist,cum_n_u_h,X_lamres,prev_LMatrices,
                         userLfixeds,hessUL,hessFac,w.resid,processed) {
    nrand <- length(ZAlist)
    locprocessed <- processed
    locprocessed$ZAlist <- NULL
    locprocessed$X_lamres <- NULL
    next_LMatrices <- prev_LMatrices
    Xi_cols <- attr(X_lamres,"Xi_cols")
    cum_Xi_cols <- attr(X_lamres,"cum_Xi_cols") 
    Lu <- u_h
    loc_lambda_est <- numeric(length(u_h))
    for (rt in seq_len(length(ZAlist))) {
      ## estimate correlation matrix 
      Xi_ncol <- Xi_cols[rt]
      blocksize <- ncol(ZAlist[[rt]])/Xi_ncol 
      ## cov mat of u_h if not fixed by user ## standard REML method 
      if ( Xi_ncol>1 && ! userLfixeds[rt]) {
        COVpredUnderHuncorr <- matrix(0,ncol=Xi_ncol,nrow=Xi_ncol) ## var on diag, corr outside diag
        prevL <- attr(prev_LMatrices[[rt]],"Lcompact")
        compactLv <- matrix(0,nrow=Xi_ncol,ncol=Xi_ncol)
        lowerbloc <- lower.tri(compactLv,diag=TRUE) ## a matrix of T/F !
        ##
        ## Build Lmatrix between all pairs of u_h (nr*nr) from parameter estimates (2*2) for the design matrix,
        makelong <- function(Lcompact) {
          longLv <- diag(ncol(ZAlist[[rt]])) ## declaration
          for (it in seq_len(Xi_ncol)) {
            urange1 <- (it-1)*blocksize + seq(blocksize)
            diag(longLv)[urange1] <- Lcompact[it,it]
            for (jt in seq_len(it-1)) {
              urange2 <- (jt-1)*blocksize + seq(blocksize)
              diag(longLv[urange1,urange2]) <- Lcompact[it,jt]
              diag(longLv[urange2,urange1]) <- Lcompact[jt,it]
            }
          }
          longLv
        } ## end def makelong
        ##
        u.range <- (cum_n_u_h[rt]+1L):(cum_n_u_h[rt+1L])
        ########## brute force optimization
        makeLcovLt <- function(parvec) {
          compactLv[lowerbloc] <- parvec
          compactLv[t(lowerbloc)] <- parvec
          sigmas <- diag(exp(diag(compactLv))) 
          diag(compactLv) <- 1
          resu <- sigmas %*% compactLv %*%sigmas
          resu
        }
        ####
        objfn <- function(parvec) {
          compactcovmat <- makeLcovLt(parvec)
          ## cosmetic / interpretative permutation
          blob <- selfAdjointSolverCpp(compactcovmat) ## COVcorr= blob$u %*% diag(blob$d) %*% t(blob$u) ## 
          blib <- Matrix::expand(Matrix::lu(t(blob$u))) ## use pivoting in lu as a useful permutation...
          blob <- list(u=t(as.matrix(with(blib,L %*% U))),d=blob$d[blib$P@perm]) 
          ## assignments as design matrix and lambda values:
          loc_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) 
          loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 ## arbitrarily small eigenvalue is possible for corr=+/-1 even for 'large' parvec
          Lcompact <- blob$u  ## the variances are taken out in $d
          ## we have a repres in terms of ZAL and of a diag matrix of variances; only the latter affects hlik computation
          longLv <- makelong(Lcompact)
          next_LMatrices[[rt]] <- longLv
          attr(next_LMatrices[[rt]],"ranefs") <- attr(ZAlist,"ranefs")[[rt]] ## FR->FR  revoir pour matrices affectant +s termes ?
          ZALlist <- 
            compute.ZALlist(LMatrix=next_LMatrices,ZAlist=ZAlist,Groupings=Groupings)
          ZAL <- post.process.ZALlist(ZALlist,predictor=locprocessed$predictor)
          attr(locprocessed$predictor,"ZALMatrix") <- ZAL
          locTT <- cbind(rbind(X.pv,OO1),attr(ZAL,"ZALI"))
          locw.ranefSblob <- 
            updateW_ranefS(cum_n_u_h,rand.families,lambda=loc_lambda_est,u_h,v_h) 
          auglinmodblob <- 
            auglinmodfit(TT=locTT,ZAL=ZAL,lambda_est=loc_lambda_est,
                         wranefblob=locw.ranefSblob,
                         d2hdv2=d2hdv2,qr.d2hdv2=qr.d2hdv2,w.resid=w.resid,beta_eta=beta_eta,
                         maxit.mean=maxit.mean,eta=eta,u_h=u_h,v_h=v_h,Sig=Sig,
                         control.HLfit=control.HLfit,
                         X.pv=X.pv,etaFix=etaFix,
                         cum_n_u_h=cum_n_u_h,psi_M=psi_M,
                         muetablob=muetablob,family=family,prior.weights=prior.weights,
                         phi_est=phi_est,verbose=verbose,
                         ranFix=ranFix,
                         corrPars=corrPars, 
                         processed=processed
            )
          if ( distinct.X.ReML ) {  
            stop("'code missing for non standard REMLformula in random slope estimation'")
          }
          locd2hdv2 <- calcDhDv2(ZAL,w.resid=auglinmodblob$w.resid,
                                 auglinmodblob$wranefblob$w.ranef)
          aphls <- calc.p_v(mu=auglinmodblob$muetablob$mu,u_h=auglinmodblob$u_h,
                            dvdu=auglinmodblob$wranefblob$dvdu,
                            lambda_est=loc_lambda_est,phi_est=phi_est,
                            d2hdv2=locd2hdv2,cum_n_u_h=cum_n_u_h,
                            lcrandfamfam=lcrandfamfam,processed=processed,
                            family=family,prior.weights=prior.weights,returnLad=FALSE)
          if (ncol(X.Re)>0) {
            hessnondiag <- crossprod(ZAL,sweep(X.Re,MARGIN=1,auglinmodblob$w.resid,`*`))  
            Md2hdbv2 <- rbind(cbind(ZtWZ(X.Re,auglinmodblob$w.resid), t(hessnondiag)),
                              cbind(hessnondiag, - locd2hdv2)) 
            ladbv <- LogAbsDetWrap(Md2hdbv2/(2*pi))
          } else { ## fit ML: p_bv=p_v hence d2hdpbv reduces to d2hdv2
            Md2hdbv2 <- - d2hdv2 
            ladbv <- calcpv$lad
          }
          REMLcrit <- aphls$hlik-ladbv/2
          return(REMLcrit)
        } ## currently this refits the fixed effects together with the other params... probably not optimal
        ####  
        lowerb <- upperb <- matrix(NA,nrow=Xi_ncol,ncol=Xi_ncol)
        diag(lowerb) <- log(sqrt(1e-08))
        diag(upperb) <- log(sqrt(1e08))
        lowerb[2,1] <-   -(1-1e-08)
        upperb[2,1] <-   (1-1e-08)
        init <- attr(prev_LMatrices[[rt]],"par")
        if (is.null(init)) {
          init <- (upperb+lowerb)/2
          diag(init) <- 0
          init <- init[lowerbloc]        
        }
        upperb <- upperb[lowerbloc]
        lowerb <- lowerb[lowerbloc]
        parscale <- (upperb-lowerb)        
        ################# OPTIM
        optr <- optim(init,objfn,lower=lowerb,upper=upperb,method="L-BFGS-B",
                      control=list(parscale=parscale,fnscale=-1))
        ################# 
        ## reproduces representation in objfn
        COVcorr <- makeLcovLt(optr$par)
        blob <- selfAdjointSolverCpp(COVcorr) ## COVcorr= blob$u %*% diag(blob$d) %*% t(blob$u) ## 
        blib <- Matrix::expand(Matrix::lu(t(blob$u))) ## use pivoting in lu as a useful permutation...
        blob <- list(u=t(as.matrix(with(blib,L %*% U))),d=blob$d[blib$P@perm]) ## + jolies façon de permuter $d ?
        loc_lambda_est[u.range] <- rep(blob$d,rep(blocksize,Xi_ncol)) 
        loc_lambda_est[loc_lambda_est<1e-08] <- 1e-08 
        Lcompact <- blob$u  #
        next_LMatrix <- makelong(Lcompact) ## il faut updater pour estimer les ranef correctement...
        attr(next_LMatrix,"Lcompact") <- Lcompact ## kept for updating in next iteration and for output
        attr(next_LMatrix,"par") <- optr$par ## kept for updating in next iteration and for output
        attr(next_LMatrix,"ranefs") <- attr(ZAlist,"ranefs")[rt]
      } else next_LMatrix <- NULL
      next_LMatrices[[rt]] <- next_LMatrix
    } ## loop on rt = ranefs
    return(list(next_LMatrices=next_LMatrices,next_lambda_est=loc_lambda_est,
                latest.unique.cov=optr$par[2]))
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
      if (is.null(u_h)) {u_h <- u_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam)}
      lambda_est <- resize.lambda(lambda.Fix)
      eta <- eta + ZAL %*% etaFix$v_h ## updated at each iteration
    } ## FREQS
    ## conversion to mean of response variable (COUNTS for binomial)
    muetablob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen,priorWeights=prior.weights) 
    mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
    w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
    if (models[[1]]=="etaHGLM") { ## linear predictor for mean with ranef
      wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## no fit, likelihood computation
      dvdu <- wranefblob$dvdu
      w.ranef <- wranefblob$w.ranef
      d2hdv2 <- calcDhDv2(ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    }
    return(list(APHLs=calc.p_v(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,d2hdv2=d2hdv2,
                               cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,processed=processed,
                               family=family,prior.weights=prior.weights))) ### RETURN !! ## FR->FR but p_bv is not returned.
  } 
  
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  
  `provide.resglm` <- function() { ## family,y,pforpv,off,prior.weights
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
    resglm <- glm(locform,family=family,offset=off,weights=prior.weights) 
    if (pforpv>0) {
      ## Two potential problems (1) NA's pour param non estimables (cas normal); 
      ## (2) "glm.fit: fitted probabilities numerically 0 or 1 occurred" which implies separation or large offset
      if (max(abs(c(coefficients(resglm))),na.rm=TRUE)>1e10) { ## na.rm v1.2 
        message("(!) Apparent divergence of estimates in a *glm* analysis of the data.")
        message("    Check your data for separation or bad offset values.")
        stop("    I exit.") 
      } 
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
          fv <- fitted(resglm)
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
  
  
  ##############################################################################################
  ######### Initial estimates for mu by GLM ####################
  if ( ( pforpv>0 && is.null(init.HLfit$fixef)) || is.null(phi.Fix) || is.null(init.HLfit$v_h) || is.null(lambda.Fix) ) { 
    ## all cases where an initial resglm is needed (even when pforpv=0, may be needed to provide init phi or init lambda)
    resglm <- provide.resglm()   
  }
  beta_eta <- numeric(pforpv)
  if (pforpv>0) { 
    beta_eta <- init.HLfit$fixef
    if (is.null(beta_eta) ) {
      beta_eta<-c(coefficients(resglm)) ## this may include NA's. Testcase: HLfit(Strength ~ Material*Preheating+Method,data=weld)
      names(beta_eta) <- unlist(lapply(names(beta_eta),substring,first=5)) ## removes "X.pv" without guessing any order or length
    } 
    beta_eta[names(etaFix$beta)] <- etaFix$beta
  } 
  if (any(is.na(beta_eta))) {   
    ## FR->FR das preprocess en utilisant lm ? mais l'interet et de montrerles NA explicites dans la sortie comme par glm()
    XpvOri <- X.pv
    namesOri <- colnames(XpvOri)
    pforpvori <- pforpv
    beta_etaOri <- beta_eta
    names(beta_etaOri) <- namesOri 
    validbeta <- which(!is.na(beta_eta))
    beta_eta <- beta_eta[validbeta]
    X.pv <- X.pv[,validbeta,drop=FALSE]
    pforpv <- ncol(X.pv)
    if (ncol(X.Re)>0)   X.Re <- X.Re[,validbeta,drop=FALSE]    
  } else beta_etaOri <- NULL
  if (!is.null(control.HLfit$intervalInfo)) {
    parmcol <- attr(control.HLfit$intervalInfo$parm,"col")
    beta_eta[parmcol] <- control.HLfit$intervalInfo$init 
  }  
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
    psi_M <- unlist(lapply(seq(nrand), function(it) {
      lpsi <- switch(lcrandfamfam[it], 
                     gaussian = 0,
                     gamma = 1, 
                     beta = 1/2, 
                     "inverse.gamma" = 1
      )
      rep(lpsi,vec_n_u_h[it])
    })) ## rand.families is processed$ and thus it is unsafe to set psi_M as an attribute to it...
    v_h <- etaFix$v_h
    if (is.null(v_h) ) v_h <- init.HLfit$v_h
    if (is.null(v_h) ) {
      v_h <- unlist(lapply(seq(nrand), function(it){
        u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
        rand.families[[it]]$linkfun(psi_M[u.range]) ## v as link(mean(u)) 
      }))
    }
    u_h <- u_h_from_v_h(v_h,rand.families=rand.families,cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam) 
    checkv_h <- attr(u_h,"v_h") ## verif is v_h was corrected 
    if (!is.null(checkv_h)) v_h <- checkv_h
    init.lambda <- generateInitLambda()
    ## one could imagine fixing some of then but not others...  
    lambda_est <- resize.lambda(init.lambda)
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
  muetablob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen,priorWeights=prior.weights) 
  mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
  dmudeta<- muetablob$dmudeta ## if Bin/Pois, must be O(n)
  Vmu <- muetablob$Vmu ## if Bin/Pois, O(n)
  w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
  
  if (models[[1]]=="etaHGLM") {
    wranefblob <- updateW_ranefS(cum_n_u_h,rand.families,lambda_est,u_h,v_h) ## initilization !
    w.ranef <- wranefblob$w.ranef
    dlogWran_dv_h <- wranefblob$dlogWran_dv_h ## (1/Wran)(dWran/dv_h), the thing between P and K2 in the d3 coef. See LeeL12 appendix
    dvdu <- wranefblob$dvdu
  }
  #betaV <- c(beta_eta,v_h) 
  conv.phi <- FALSE; conv.lambda <- FALSE; conv.corr <- FALSE
  if (models[[1]]=="etaHGLM") {
    Sig <- ZWZt(ZAL,1/w.ranef) + diag(1/w.resid) ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid) 
    d2hdv2 <- calcDhDv2(ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    #Sig <- sweep(ZAL,MARGIN=2,1/w.ranef,`*`)  %*% t(ZAL) + diag(1/w.resid) ## ZAL %*% diag(1/w.ranef) %*% t(ZAL) + diag(1/w.resid) 
    #d2hdv2 <- - sweep(t(ZAL),MARGIN=2,w.resid,`*`) %*% ZAL - diag(w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    qr.d2hdv2 <- NULL
    OO1 <- matrix(0,cum_n_u_h[nrand+1L],pforpv)
    TT <- cbind(rbind(X.pv,OO1),attr(ZAL,"ZALI"))  ## aug design matrix
    if ( distinct.X.ReML ) {  
      OO1leve <- matrix(0,cum_n_u_h[nrand+1L],ncol(X.Re))
      TTleve <- cbind(rbind(X.Re,OO1leve),attr(ZAL,"ZALI"))
    }
    if (length(etaFix$beta)==ncol(X.pv) && !is.null(etaFix$v_h)) {
      maxit.mean <- 0 ## used in test near the end...
    } else if ( LMMbool ) {
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
    qr.XtX <- QRwrap(crossprodCpp(X.pv)) ## qr(t(X.pv)%*%X.pv) 
    SEMargs$beta_eta <- beta_eta
    betalambda <- do.call(SEMbetalambda,SEMargs)
    beta_eta <- betalambda$beta_eta
    lambda_est <- resize.lambda(betalambda$lambda)
    u_h <- v_h <- betalambda$vs
    logLapp <- betalambda$logLapp
    attr(logLapp,"seInt") <- betalambda$seInt ## may be NULL
    tXinvS <- NULL
  } else while ( TRUE ) { ##the main loop with steps: new linear predictor, new leverages, new phi, new w.resid, new lambda, new fn(lambda)
    if (models[[1]]=="etaHGLM") {
      ##############################
      auglinmodblob <- auglinmodfit(TT=TT,ZAL=ZAL,lambda_est=lambda_est,wranefblob=wranefblob,
                                    d2hdv2=d2hdv2,qr.d2hdv2=qr.d2hdv2,w.resid=w.resid,beta_eta=beta_eta,
                                    maxit.mean=maxit.mean,eta=eta,u_h=u_h,v_h=v_h,Sig=Sig,
                                    control.HLfit=control.HLfit,
                                    X.pv=X.pv,etaFix=etaFix,
                                    cum_n_u_h=cum_n_u_h,psi_M=psi_M,
                                    muetablob=muetablob,family=family,
                                    prior.weights=prior.weights,phi_est=phi_est,verbose=verbose,
                                    ranFix=ranFix,corrPars=corrPars,
                                    processed=processed
                                    ) ## HL(.,.) estim of beta, v for given lambda,phi
      ##############################
      beta_eta <- auglinmodblob$beta_eta
      v_h <- auglinmodblob$v_h
      u_h <- auglinmodblob$u_h
      eta <- auglinmodblob$eta
      wranefblob <- auglinmodblob$wranefblob
      w.ranef <- wranefblob$w.ranef ; dlogWran_dv_h <- wranefblob$dlogWran_dv_h ; dvdu <- wranefblob$dvdu
      muetablob <- auglinmodblob$muetablob
      mu <- muetablob$mu
      dmudeta <- muetablob$dmudeta
      Vmu <- muetablob$Vmu
      w.resid <- auglinmodblob$w.resid
      Sig <- auglinmodblob$Sig
      d2hdv2 <- auglinmodblob$d2hdv2
      qr.d2hdv2 <- auglinmodblob$`qr.d2hdv2`
      wAugX <- auglinmodblob$wAugX
      tXinvS <- auglinmodblob$tXinvS
      sqrt.ww <- auglinmodblob$sqrt.ww
      innerj <- auglinmodblob$innerj
      levQ <- auglinmodblob$levQ
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
          muetablob <- muetafn(family=family,eta=eta,BinomialDen=BinomialDen,priorWeights=prior.weights) 
          mu <- muetablob$mu ## if Bin/Pois, O(n): facteur BinomialDen dans la transfo mu -> eta ## nonstandard mu des COUNTS
          dmudeta <- muetablob$dmudeta
          Vmu <- muetablob$Vmu ## if Bin/Pois, O(n)
          ## update functions of v_h -> blob
          w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases
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
          wAugXleve <- sweepZ1W(TTleve,sqrt.ww)
          hatval <- leverages(wAugXleve)
        } else { ## basic REML, leverages from the same matrix used for estimation of betaV (even simply V)
          if (.spaMM.data$options$USEEIGEN) {
            if (is.null(levQ)) {
              if (FALSE)  {
                ## wAugX updated not only by change in lambda, phi but also GLM weights -> leverage comput difficult to optimize  
                ## the following could be useful if the GLM wights are unity, phiScal et lamScal...
                # mais il manque bouts pour pour produire <u> et <d> | unperturbed RpR = u d u'                
                #                hatval <- LevPerturbedQCpp(perturbedwAugX=wAugX,pforREML=ncol(X.Re),
                #                                           RpRu = <u>,RpRd=<d>,lamOverLam0=lambda/lambda0,phiOverPhi0=phi/phi0)
              } else hatval <- leverages(wAugX)
            } else hatval <- rowSums(levQ^2) ## if we have levQ, we use it   
          } else hatval <- rowSums(qr.qy(auglinmodblob$SQR, 
                                         diag(1, nrow = nrow(wAugX), ncol = ncol(wAugX)))^2) ## == but faster than rowSums(qr.Q(SQR)^2) !
        }
        if (any(abs(hatval) > 1 - 1e-8)) {
          hatval <- ifelse(abs(hatval) > 1 - 1e-8, 1 - 1e-8,hatval)
          warningList$leveLam1 <-TRUE
        }
        lev_phi <- hatval[1:nobs] ## for the error residuals (phi)
        lev_lambda <- hatval[(nobs+1L):(nobs+cum_n_u_h[nrand+1L])]  ## for the ranef residuals (lambda)
      }
    } else { ## GLM
      if ( distinct.X.ReML ) { 
        wAugXleve <- sweep(X.Re,MARGIN=1,sqrt(w.resid),`*`) # rWW%*%X.Re
        lev_phi <- leverages(wAugXleve)
      } else { ## basic REML, leverages from the same matrix used for estimation of beta
        wAugX <- sweep(X.pv,MARGIN=1,sqrt(w.resid),`*`) # rWW %*% X.pv 
        lev_phi <- leverages(wAugX)
      }
    }
    d2mudeta2 <- d2mudeta2fn(link=family$link,mu=mu,eta=eta,BinomialDen=BinomialDen)
    if (verbose["trace"]) {print(paste("beta=",paste(signif(beta_eta,4),collapse=", ")),quote=F)}
    if (HL[2]>0) { ## LeeN01 HL(.,1) ie the + in 'EQL+'
      ## (0): previous hat matrix -> p, notEQL -> tilde(p), (1): full correction -> q 
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
              if (.spaMM.data$options$USEEIGEN) {
                dvdloglamMat <- pseudoSolvediag(d2hdv2,as.vector(neg.d2f_dv_dloglam)) ## FR->FR dangereux car contient (Eigen:)solve(R)
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
            dleve <- ((hatval * dlW_deta_or_v) %*% attr(ZAL,"ZALI") ) %*% dvdloglamMat # (r+n) . (r+n)Xr . rXr = r (each element is a sum over r+n terms= a trace)
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
            dleve <- ((hatval * dlW_deta_or_v) %*% attr(ZAL,"ZALI")) %*% dvdlogphiMat # (r+n) . (r+n)Xr . rXn = n (each element is a sum over r+n terms= a trace)
            lev_phi <- lev_phi - as.vector(dleve)  
          } 
        }
      } 
    }
    if (HL[2]>1) {
      stop("Need a_i correction in Table 7 of NohL07 ie derivatives of second order correction wrt dips param.")
    }
    if (HL[3]!=0 ) {## HL(.,.,1) ie , p_bv(h), not EQL p_bv(q+), LeeNP p89; distinction does not arise for PQL <=> Gaussian ranefs...  
      # lambda
      if (models[[1]]=="etaHGLM" && is.null(lambda.Fix)) ## d h/ d !log! lambda coorection     
        lev_lambda <- lev_lambda + corr.notEQL.lambda(nrand,cum_n_u_h,lambda_est,lcrandfamfam) 
      # phi hence not poiss,binom:
      if (family$family=="Gamma" && is.null(phi.Fix) ) { ## d h/ d !log! phi correction (0 for gauss. resid. error). Not tied to REML
        phiscaled <- phi_est*prior.weights ## 08/2014 ## trick for still using deviances residuals in the Gamma GLM
        lev_phi <- lev_phi +  1+2*(log(phiscaled)+digamma(1/phiscaled))/phiscaled ## LNP p. 89 and as in HGLMMM IWLS_Gamma
      }    
    }
    ## updated residuals from updated mu must be used (LeeNP p.161) [not so in dhglmfit !!]
    deviance_residual <- family$dev.resids(y,mu,wt=1) 
    ######### Dispersion Estimates for phi #####################
    resid.family <- processed$resid.family
    phifam <- GammaForDispGammaGLM(resid.family$link)
    if (is.null(phi.Fix)) { ## if phi is estimated (phi.Fix set to 1 for Poisson, Binomial)
      ## leverages have been computed before the  inner loop, which did not change the design matrices 
      Offset_disp <- attr(resid.predictor,"offsetObj")$vector 
      lev_phi <- pmin(lev_phi, 1 - 1e-8)
      locmethod <- "glm" ## the method should not matter here
      if (models[[3]]=="phiScal") { ## 
        if (attr(resid.predictor,"offsetObj")$nonZeroInfo) {
          resglm_phi <- dispGammaGLM(dev.res=deviance_residual*prior.weights,lev=lev_phi,
                                     X=X_disp,offset=Offset_disp,family=phifam,method=locmethod)
          if (! is.null(locw <- resglm_phi$warnmess)) warningList$innerPhiGLM <- locw
          beta_phi <- resglm_phi$beta_disp
          next_phi_est <- resglm_phi$next_disp_est
        } else { ## one case where we can easily avoid an explicit call to a glm (but one will be used to compute SEs later) 
          next_phi_est <- sum(deviance_residual*prior.weights)/sum(1-lev_phi) ## NOT in linkscale
          beta_phi <- c("(Intercept)"=phifam$linkfun(next_phi_est)) ## linkscale value
          resglm_phi <- NULL
        }
        if (verbose["trace"]) {print(paste("phi_est=",signif(next_phi_est,4)),quote=F)}
      } else if (models[[3]]=="phiGLM") { ## there is a phi predictor to estimate but no ranef in this predictor
        resglm_phi <- dispGammaGLM(dev.res=deviance_residual*prior.weights,lev=lev_phi,
                                   X=X_disp,offset=Offset_disp,family=phifam,method=locmethod)
        if (! is.null(locw <- resglm_phi$warnmess)) warningList$innerPhiGLM <- locw
        # glm(I(deviance_residual/(1-lev_phi))~X_disp-1,weights=(1-lev_phi)/2,family=Gamma(log))
        beta_phi <- resglm_phi$beta_disp
        next_phi_est <- resglm_phi$next_disp_est
      } else if (models[[3]]=="phiHGLM") { ## random effect(s) in predictor for phi
        stop("random effects in predictor or residual variance (phi) not yet implemented")
        ## there is a template code with comments in version 260812 of HLfit
        reshglm_phi <- list()
      } 
      #lowphi <- which(next_phi_est < 1e-08); next_phi_est[lowphi] <- 1e-08 ## 2014/09/04 better correction in calc.p_v
      if (all(abs(next_phi_est-phi_est) < conv.threshold* (phi_est+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.phi <- TRUE ## 'weak convergence'... 
      } else conv.phi <- FALSE
    } else {conv.phi <- TRUE} ## there is a phi.Fix, already -> phi_est
    ######### Dispersion Estimates for lambda #####################
    if (models[[1]]=="etaHGLM" && is.null(lambda.Fix)) { ## lambda must be estimated     
      if (any(abs(lev_lambda) > 1 - 1e-8)) { ## abs... not commented when written...
        lev_lambda <- ifelse(abs(lev_lambda) > 1 - 1e-8, 1 - 1e-8,lev_lambda)
        warningList$leveLam1 <- TRUE
      }
      ## Build pseudo response for lambda GLM/HGLM
      resp_lambda <- matrix(0,cum_n_u_h[nrand+1L],1L)
      #########################
      for (it in 1:nrand) {
        u.range <- (cum_n_u_h[it]+1L):cum_n_u_h[it+1L]
        resp_lambda[u.range] <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=1) ## must give d1 in table p 989 de LeeN01
      }
      ## then analyse pseudoresponse
      if (all(models[[2]]=="lamScal")) { ## simplified estimate
        if (any(attr(X_lamres,"Xi_cols")>1)) {
          ## handling correlation in random slope models # slmt pr gaussian ranefs, verif dans preprocess
          LMatricesBlob <- makeCovEst(u_h,ZAlist=ZAlist,cum_n_u_h=cum_n_u_h
                                      ,X_lamres=X_lamres,prev_LMatrices=next_LMatrices,
                                      userLfixeds=attr(ZALlist,"userLfixeds"),
                                      hessUL=ZtWZ(X.Re,w.resid),
                                      hessFac=sweep(X.Re,MARGIN=1,w.resid,`*`),
                                      w.resid=w.resid,
                                      processed=processed                            
          )            
          next_LMatrices <- LMatricesBlob$next_LMatrices ## a list of matrices with NULL elements for non-random-slope terms
          next_lambda_est <- LMatricesBlob$next_lambda_est ## a full-length vector with values only in the appropriate u ranges 
          ## only for testing convergence: 
          cov12_est <- next_cov12_est
          next_cov12_est <- LMatricesBlob$latest.unique.cov
        } else next_lambda_est <- numeric(length(u_h)) ## next_LMatrices remains an empty list()
        ## then fill all missing values i.e. for terms without random slope 
        for (it in seq_len(nrand)) {
          if (attr(X_lamres,"Xi_cols")[it]==1) {
            u.range <- (cum_n_u_h[it]+1L):(cum_n_u_h[it+1L])
            unique.lambda <- sum(resp_lambda[u.range])/sum(1-lev_lambda[u.range]) ## NOT in linkscale 
            unique.lambda <- pmax(unique.lambda,1e-8) # FR->FR still corrected
            unique.lambda <- pmin(unique.lambda,.spaMM.data$options$maxLambda)  
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
        ## FR->FR code non operationnel
        resglm_lambda <- dispGammaGLM(dev.res=resp_lambda,lev=lev_lambda,X=X_lamres)
        if (! is.null(locw <- resglm_lambda$warnmess)) warningList$innerLamGLM <- locw
        next_lambda_est <- resglm_lambda$next_disp_est ## $fitted.value is NOT in linkscale, contrairement a $coefficients
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
    if (! is.null(corr_est)) { ## this code does not apply for the random slope model
      corrEst.args$formula <- Predictor(formula=corrEstBlob$corrEst.form,offset=off + X.pv %*% beta_eta) ## FR->FR ugly input for offset
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
      next_corr_est <- pff$corrPars[corrNames_in_init_HLfit] ## rho,nu,  pas trRho, trNu 
      #FR->FR maybe conv_threshold a bit strict here...
      if (all(abs(log(unlist(next_corr_est)/unlist(corr_est))) < conv.threshold) ) { ## 
        conv.corr <- TRUE ## this is the simplest, best case. ## but if slow geometric decrease to 0, this is never true 
      } else if (all(abs(unlist(next_corr_est)-unlist(corr_est)) < conv.threshold* (unlist(corr_est)+0.1)) ) { ## ie 1e-6 ~ 1e-5*(1e-6+0.1)
        conv.corr <- TRUE ## 'weak convergence'... 
      } else conv.corr <- FALSE
    } else {
      if (!is.null(next_cov12_est)) {
        if (iter>1 && abs(cov12_est-next_cov12_est) < conv.threshold ) { 
          conv.corr <- TRUE 
        } else conv.corr - FALSE       
      } else conv.corr <- TRUE
##      print(cov12_est-next_cov12_est)
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
        w.resid <- calc.w.resid(muetablob$GLMweights,phi_est) ## 'weinu', must be O(n) in all cases; blob was updated when eta was
      }
      if (length(next_LMatrices)>0 || ! is.null(corr_est) ) {
        if (length(next_LMatrices)>0) {
          ZALlist <- compute.ZALlist(LMatrix=next_LMatrices,ZAlist=ZAlist,Groupings=Groupings)
        } else if (! is.null(corr_est)) {
          corr_est <- next_corr_est 
          LMatrix <- attr(pff$predictor,"LMatrix")
          ZALlist <- compute.ZALlist(LMatrix=LMatrix,ZAlist=ZAlist,Groupings=Groupings)
        }
        ZAL <- post.process.ZALlist(ZALlist,predictor=predictor)
        TT <- cbind(rbind(X.pv,OO1),attr(ZAL,"ZALI"))  ## aug design matrix
        if (distinct.X.ReML) TTleve <- cbind(rbind(X.Re,OO1leve),attr(ZAL,"ZALI"))
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
     #      beta_se <- rep(Inf,pforpv) ## maybe...
      beta_cov <- matrix(Inf,ncol=pforpv,nrow=pforpv) ## maybe...
    } else {
      beta_cov <- try(solve(tXinvS%*%X.pv),silent=TRUE) ## solve(small matrix !)
      if (class(beta_cov)=="try-error") {
        #        beta_se <- rep(Inf,pforpv) ## maybe...
        beta_cov <- matrix(Inf,ncol=pforpv,nrow=pforpv) ## maybe...
      } else {
        #        beta_se <- diag(beta_cov)
        if (any(diag(beta_cov)<0)) { ## divergence of tXinvS%*%X.pv leads to negative variance estimates
         #          beta_se <- rep(Inf,pforpv) ## maybe... 
          beta_cov <- matrix(Inf,ncol=pforpv,nrow=pforpv) ## maybe...
        } #else beta_se <- sqrt(beta_se)
      }
    }
    #    if (any(is.infinite(beta_se))) {
    #      message("Suspected divergence of lambda estimates. Check model formula (wrong offset for example),")
    #      message(" otherwise try increasing control.HLfit$iter.mean.dispVar")
    #    }
  } else {#beta_se <- NULL
    beta_cov <- NULL
  } 
  ######################
  ######################
  ######################
  ##### LAMBDA
  if (HL[1]!="SEM" && models[[1]]=="etaHGLM" && is.null(lambda.Fix)) {
    if (all(models[[2]]=="lamScal")) { ## there is a single X_lamres for all lambda's
      ## to compute the se we need the GLM residuals etc. So if the GLM has not been previously used it's better to use it here
      resglm_lambda <- dispGammaGLM(dev.res=resp_lambda,lev=lev_lambda,X=X_lamres)
      linkscale.lambda <- resglm_lambda$beta_disp
      p_lambda <- ncolX_lam <- ncol(X_lamres) ## to be modifie
      lambda_se <- summary(resglm_lambda$resglm,dispersion=1)$coefficients[(ncolX_lam+1L):(2L*ncolX_lam)]       
      if (! is.null(next_LMatrices)) {
        ## lignes suiv supposent que L_matrix decrit random slope model
        p_corr <- sum(unlist(lapply(next_LMatrices,function(mat) {
          dimL <- nrow(attr(mat,"Lcompact"))
          (dimL-1)*dimL/2
        })))
        p_lambda <- p_lambda+p_corr
      }
    } else {
      stop("From HLfit: 'lamHGLM' and 'lamGLM' not fully implemented.")
      ## there is a template code with comments in version 260812 of HLfit
    }
  } else p_lambda <- 0       
  ##### PHI
  if ( is.null(phi.Fix)) {
    if (models[[3]]=="phiHGLM") {
      ## there is a template code with comments in version 260812 of HLfit
      stop("HGLM for phi not implemented")
    } else {
      if (is.null(resglm_phi)) {
        resglm_phi <- dispGammaGLM(dev.res=deviance_residual*prior.weights, 
                                   lev=lev_phi, X=X_disp, offset=Offset_disp, family=phifam)
      }
      phi_se <- summary(resglm_phi$resglm,dispersion=1)$coefficients[(p_phi+1L):(2L*p_phi)]       
      ## note dispersion set to 1 to match SmythHV's 'V_1' method, which for a log link has steps:
      #SmythHVsigd <- as.vector(sqrt(2)*phi_est);SmythHVG <- as.vector(phi_est); tmp <- SmythHVG / SmythHVsigd ## tmp is here sqrt(2) !
      #if (length(tmp)>1) {SmythHVZstar <- diag(tmp) %*% X_disp} else SmythHVZstar <- tmp * X_disp
      #SmythHVcovmat <- solve(ZtWZ(SmythHVZstar,(1-lev_phi))); phi_se <- sqrt(diag(SmythHVcovmat)) print(phi_se)
    }
  } 
  ########## LIKELIHOODS
  #    theta<-theta.mu.canonical(mu/BinomialDen,family$family)  
  if (HL[1]=="SEM") {
    APHLs <- list(logLapp=logLapp) ## keeps attributes
  } else {
    if (models[[1]]=="etaHGLM" && pforpv==0) { 
      d2hdv2 <- calcDhDv2(ZAL,w.resid,w.ranef) ##  - t(ZAL) %*% diag(w.resid) %*% ZAL - diag(w.ranef)
    }
    calcpv <- calc.p_v(mu=mu,u_h=u_h,dvdu=dvdu,lambda_est=lambda_est,phi_est=phi_est,d2hdv2=d2hdv2,
                       cum_n_u_h=cum_n_u_h,lcrandfamfam=lcrandfamfam,processed=processed,
                       family=family,prior.weights=prior.weights,returnLad=TRUE)
    if (models[[1]] != "etaHGLM" && models[3] != "phiHGLM") { ## ie GLM, not HGLM
      ## note that p_v=p_bv here, whether an REML estimation of phi is used or not... 
      ml <- calcpv$clik ## vanilla likelihood
      d2hdx2 <- - ZtWZ(X.Re,w.resid)  ## t(X.Re)%*%Wresid%*%X.Re ## X should be the one for leverages
      lad <- LogAbsDetWrap(d2hdx2/(2*pi))
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
          hessnondiag <- crossprod(ZAL,sweep(X.Re,MARGIN=1,w.resid,`*`))  
          Md2hdbv2 <- rbind(cbind(ZtWZ(X.Re,w.resid), t(hessnondiag)),
                          cbind(hessnondiag, - d2hdv2)) 
          ladbv <- LogAbsDetWrap(Md2hdbv2/(2*pi))
          if (AIC) { ## diff de d2hdbv2 slmt dans dernier bloc (FR->FR AIC on REML ????)
            Md2clikdbv2 <- rbind(cbind(ZtWZ(X.Re,w.resid), t(hessnondiag)),
                                 cbind(hessnondiag, ZtWZ(ZAL,w.resid)))            
          }
        } else { ## fit ML: p_bv=p_v hence d2hdpbv reduces to d2hdv2
          Md2hdbv2 <- - d2hdv2 
          ladbv <- calcpv$lad
          if (AIC) Md2clikdbv2 <-  ZtWZ(ZAL,w.resid) ## for AIC
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
      eigvals <- eigen(Md2hdbv2/(2*pi),only.values = T)$values
      eigvals <- pmax(eigvals,1e-12)
    }
    if (is.nan(ladbv)) ladbv <- sum(log(eigvals)) 
    p_bv <- hlik-(hv10+hv20+ladbv/2)  
    if ( ! is.null(calcpv$second.corr)) p_bv <- p_bv + calcpv$second.corr
    dAIC <- -2*p_bv + 2 * (p_lambda+p_phi) ## HaLM (10) focussed for dispersion params
    if ( AIC ) {
      # a debugging issue is that options(error=recover) acts before tryCatch gets the return value
      # from its first argument. So a tryCatch on solve is not a good idea.
      if (min(eigvals)>1e-11) {
        qr.Md2hdbv2 <- QRwrap(Md2hdbv2)
        pd <- sum(diag(solveWrap.matrix(qr.Md2hdbv2,Md2clikdbv2,stop.on.error=stop.on.error)))
        if (class(pd)=="try-error") {
          warning("Computation of cAIC/GoF df's failed because the 'd2hdbv2' matrix appears singular")
          pd <- NA
        }
      } else pd <- Inf
      GoFdf <- nobs - pd
      ## eqs 4,7 in HaLM07
      cAIC <- -2*clik + 2*pd ## 
      # there is also a "focussed" AIC in HaLM07 that would be   
      # - 2 p_bv + 2 * <number of dispersion parameters> (ie lambda,phi,nu,rho...)
      ## that would be used to select among different dispersion models 
      # discussion in section 7 of the paper suggests using an AIC based on p_v for selection among different fixed effect component models
      # but Yu and Yau then suggest caicc <- -2*clik + ... where ... involves d2h/db d[disp params] and d2h/d[disp params]2
    } else {cAIC <-NULL; GoFdf <- NULL}
    if (models[[1]] != "etaHGLM") {
      APHLs <- list(p_v=ml,p_bv=p_bv) ## FR->FR rename ?
    } else APHLs <- c(calcpv,list(p_bv=p_bv))
    APHLs$cAIC <- cAIC
    APHLs$dAIC <- dAIC
    APHLs$GoFdf <- GoFdf    
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
  res$X.pv <- X.pv
  res$ranFix <- ranFix ## currently as a uniform template consistent with projected changes ; excpt that lamFix, phiFix info is now in lambda.object, etc
  corrPars[corrNames_in_init_HLfit] <- corr_est[corrNames_in_init_HLfit]
  res$corrPars <- corrPars 
  ## FR->FR il serait logique ? de regrouper $ranFix et $corrPars dans la sortie ? Diffcile car corrPars inclut fixed and variable corr pars
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
    beta_etaOri[names(beta_eta)] <- beta_eta ## keeps the original NA's
    res$fixef <- beta_etaOri
    #    se <- beta_etaOri ## also put NA's
    #    se[names(beta_eta)] <- beta_se ## assumes that beta_eta and beta_se havesame order (the latter has no names)
    #    res$fixef_se <- se
    betaOri_cov <- matrix(0,ncol=ncol(XpvOri),nrow=ncol(XpvOri),dimnames=list(rownames=namesOri,colnames=namesOri))
    betaOri_cov[names(beta_eta),names(beta_eta)] <- beta_cov
    res$beta_cov <- betaOri_cov
  } else {
    names(beta_eta) <- colnames(X.pv)
    res$fixef <- beta_eta
  #    res$fixef_se <- beta_se
    res$beta_cov <- beta_cov
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
  if (is.null(phi.Fix)) {  
    res$resid.family <- resid.family  ## summary will use link and linkinv...
  }
  if ( distinct.X.ReML ) res$X.Re <- X.Re
  ###################
  ## ALL other LAMBDA returns
  ###################
  res$rand.families <- rand.families 
  ##
  res$ranef <- structure(u_h,cum_n_u_h=cum_n_u_h) ## FR->FR added cum_n_u_h attribute 11/2014: slightly duplicates info in lambda object
  res$v_h <- v_h
  res$w.resid <- w.resid ## useful to reconstruct Sig in predVar
  if (models[[1]]=="etaHGLM") res$w.ranef <- w.ranef ## useful to reconstruct Sig in predVar
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
      namesTerms <- attr(ZAlist,"namesTerms") ## for each random term, the names of the coefficients fitted
      if (HL[1]=="SEM") {
        lambda.object <- list(linkscale.lambda=log(print_lambda),lambda_se=NA,namesX_lamres="(Intercept)",namesTerms=namesTerms)
      } else {
        namesnames <- unlist(lapply(names(namesTerms),function(st) {
          if (nchar(st)>10) st <- paste(substr(st,0,9),".",sep="")
          st
        }))
        names(namesTerms) <- make.unique(namesnames,sep=".") ## makes group identifiers unique (names of coeffs are unchanged)
        lambda.object <- list(linkscale.lambda=linkscale.lambda,lambda_se=lambda_se,
                              namesX_lamres=colnames(X_lamres),namesTerms = namesTerms)
        attr(lambda.object,"warning") <- resglm_lambda$warnmess ## may be NULL
        if (! is.null(next_LMatrices)) {
          res$cov.mats <- lapply(next_LMatrices,function(mat) {
            ZWZt(attr(mat,"Lcompact"),exp(linkscale.lambda))
          })
        }
      }
    } else {
      lambda.object <- list(lambda.fix=lambda.Fix) ## la distinction lambda.fix doit être maintenue (summary, df des LRTs...)
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
  ## phi_est comes from the iterative algo, not from additional GLM for SE
  if (models[[3]]=="phiScal") {res$phi <- phi_est[1]} else res$phi <- phi_est
  if (is.null(phi.Fix)) {
    names(beta_phi) <- unlist(lapply(names(beta_phi),substring,first=2)) ## removes "X" without guessing any order or length
    phi.object <- list(beta_phi=beta_phi,phi_se=phi_se,namesX_disp=colnames(X_disp))
    attr(phi.object,"warning") <- resglm_phi$warnmess ## may be NULL
    res$phi.object <- phi.object
  } else {
    res$phi.object <- list(phi.Fix=phi.Fix)
  }
  ###################
  ## something to be checked
  ###################
  ## if (max.iter<1L) return(res) ## FR->FR  !! NOT class==HLfit !! FR->FR 06/2014: no longer clear what for: devel code ? 
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
  if (! is.null(locw <- warningList$innerPhiGLM)) {
    warningList$innerPhiGLM <- paste("'",locw,"' in some sub-final iteration(s) of phi estimation;", sep="")
  }
  if (! is.null(locw <- warningList$innerLamGLM)) {
    warningList$innerLamGLM <- paste("'",locw,"' in some sub-final iteration(s) of lambda estimation;", sep="")
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
  res$spaMM.version <- .spaMM.data$Constants$Version
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
      abyss <- sapply(length(seriousWarnings),function(i) {
        warning(paste("In HLfit :\n",seriousWarnings[[i]],sep=""),call.=FALSE)}) 
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

