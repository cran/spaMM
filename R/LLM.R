# cf .calc_d2mudeta2 with different arguments
.D2muDeta2 <- function(link) switch(link,
                                    "log" = function (eta) pmax(exp(eta), .Machine$double.eps),
                                    "identity" = function (eta) rep.int(0, length(eta)),
                                    "sqrt" = function (eta) rep.int(2, length(eta)),
                                    "logit" = function(eta) {
                                      expeta <- exp(eta)
                                      expeta*(1-expeta)/(1+expeta)^3
                                    }, 
                                    "probit" = function(eta) - eta * dnorm(eta, 0,1), 
                                    "cloglog" = function(eta) {
                                      eta <- pmin(eta,700) ## as in binomial(cloglog)$mu.eta
                                      expeta <- exp(eta)
                                      (1-expeta)*exp(eta-expeta)
                                    }, 
                                    "cauchit" = function(eta) { -2 *eta/(pi * (1+eta^2)^2)},
                                    "inverse" = function(eta) { 2/eta^3 }, # for gaussian(inverse) -> does not mean -1/...
                                    "loglambda" = function(eta) {stop("this function should not be called")},
                                    stop("link not yet handled in .D2muDeta2() [but easy to fix]")
)

.D3muDeta3 <- function(link) switch(link,
                                    "log" = function (eta) pmax(exp(eta), .Machine$double.eps),
                                    "identity" = function (eta) rep.int(0, length(eta)),
                                    "sqrt" = function (eta) rep.int(0, length(eta)),
                                    "logit" = function(eta) {
                                      eta <- pmin(eta,700) ## as in binomial(cloglog)$mu.eta
                                      expeta <- exp(eta)
                                      expeta*(1-4*expeta+expeta^2)/(1+expeta)^4
                                    }, 
                                    "probit" = function(eta) (eta^2-1) * dnorm(eta, 0,1), 
                                    "cloglog" = function(eta) {
                                      expeta <- exp(eta)
                                      (1-3*expeta + expeta^2)*exp(eta-expeta)
                                    }, 
                                    "cauchit" = function(eta) { (-2+6*eta^2)/(pi * (1+eta^2)^3)},
                                    "inverse" = function(eta) { -6/eta^4 }, # for gaussian(inverse) -> does not mean -1/...
                                    "loglambda" = function(eta) {stop("this function should not be called")},
                                    stop("link not yet handled in .D3muDeta3() [but easy to fix]")
)

# fixed-effect models only
.do_damped_WLS_llm <- function(grad, damping, negHess, X.pv, clik, family, old_beta_eta, phi_est, off, processed, verbose) {
  restarted <- FALSE
  dampingfactor <- 2
  while(TRUE) { 
    dampDpD <- damping * diag(negHess)
    nc <- ncol(negHess)
    diagPos <- seq.int(1L,nc^2,nc+1L)
    negHess[diagPos] <- negHess[diagPos] + dampDpD
    dbetaV <- as.vector(solve(negHess,grad)) # perturbation of solve(mH, grad) = - solve(D2logLDbeta2) . DlogLDbeta
    levMblob <- .eval_gain_clik_LevM(LevenbergMstep_result=list(dbetaV=dbetaV, dampDpD=dampDpD, rhs=grad),
                                       X.pv=X.pv, clikold=clik, family=family,
                                       coefold=old_beta_eta,
                                       phi_est=phi_est, offset=off,
                                       processed=processed)
    gainratio <- levMblob$gainratio
    conv_crit <- max(levMblob$conv_clik, abs(grad)/(1+clik))
    if (is.nan(gainratio)) {
      break
    } else if (gainratio>0) { ## success
      damping <- damping * max(1/3,1-(2*gainratio-1)^3)  
      # damping <- max(damping, .Machine$double.eps)
      dampingfactor <- 2
      break
    } else if (dampingfactor>4 ## iter not accessible for test
               && gainratio==0) { # apparently flat deviance
      if (conv_crit < 1e-8) { # we were at optimum
        damping <- 1e-7
        if (verbose["trace"]) cat("#")
        break ## apparently flat dev
      } else { ## measurable gradient, but we don't move => too high damping ? (if damping from previous LevM step was too high)
        if ( ! restarted) { # condition to avoid infinite loop
          damping <- 1e-7
          dampingfactor <- 2
          restarted <- TRUE
          if (verbose["trace"]) cat("-")
          # and continue # but this allows an  infinite loop
        } else {
          # hmm; well, break and diagnose...
          if (verbose["trace"]) cat("?")
          break
        }
      }
    } else { ## failure: increase damping and continue iterations
      damping <- dampingfactor*damping
      dampingfactor <- dampingfactor*2
    } 
    if (damping>1e10) break # stop("reached damping=1e10")
  } ## while TRUE
  RESU <- levMblob
  RESU$damping <- damping
  RESU$dbetaV <- dbetaV
  return(RESU)
}

.intervalStep_llm <- function(old_beta,X.pv, dlogcLdeta, d2logcLdeta2, currentlik,for_intervals,currentDy) {
  ## voir code avant 18/10/2014 pour une implem rustique de VenzonM pour debugage  
  ## somewhat more robust algo (FR->FR: still improvable ?), updates according to a quadratic form of lik near max
  ## then target.dX = (current.dX)*sqrt(target.dY/current.dY) where dX,dY are relative to the ML x and y 
  ## A nice thing of this conception is that if the target lik cannot be approached, 
  ##   the inferred x converges to the ML x => this x won't be recognized as a CI bound (important for locoptim) 
  # if (is.null(names(szAug))) stop("Programming error: 'szAug' must have names")
  if (is.null(colnames(X.pv))) stop("Programming error: 'X.pv' must have colnames")
  parmcol_X <- for_intervals$parmcol_X
  beta <- rep(NA,length(old_beta))
  if (currentDy <0) { 
    beta[parmcol_X] <- old_beta[parmcol_X]
  } else {
    currentDx <- (old_beta[parmcol_X]-for_intervals$MLparm)
    targetDy <- (for_intervals$fixeflik-for_intervals$targetlik)
    Dx <- currentDx*sqrt(targetDy/currentDy)
    ## pb is if Dx=0 , Dx'=0... and Dx=0 can occur while p_v is still far from the target, because other params have not converged.
    ## FR->FR patch:
    if (currentDy<targetDy) { ## we are close to the ML: we extrapolate a bit more confidently
      min_abs_Dx <- for_intervals$asympto_abs_Dparm/1000
    } else min_abs_Dx <- 1e-6 ## far from ML: more cautious move our of Dx=0
    Dx <- sign(currentDx)*max(abs(Dx),min_abs_Dx)
    beta[parmcol_X] <- for_intervals$MLparm + Dx 
  }
  #locsXaug <- sXaug[,-(parmcol_X),drop=FALSE]
  locX <- X.pv[,-(parmcol_X),drop=FALSE]
  
  #locszAug <- as.matrix(szAug-sXaug[,parmcol_X]*beta[parmcol_X])
  fullnegHess <- crossprod(X.pv, .Dvec_times_m_Matrix( - d2logcLdeta2, X.pv))
  grad <- X.pv %*% dlogcLdeta
  locgrad <- as.matrix(grad -fullnegHess[,parmcol_X]*beta[parmcol_X])
  
  #beta[-(parmcol_X)] <- get_from_MME(locsXaug,szAug=locszAug) 
  locnegHess <- fullnegHess[-(parmcol_X),-(parmcol_X)] # crossprod(locX, .Dvec_times_m_Matrix( - d2logcLdeta2[-(parmcol_X)], locX))
  loccholH <- chol(locnegHess)
  beta[-(parmcol_X)] <- backsolve(loccholH, backsolve(loccholH, locgrad, transpose = TRUE)) # solve(negHess
  return(list(beta=beta)) # levQ is presumably always dense
}


# fixed-effect mean response
# uses gradient and negHess, while .calc_etaGLMblob uses z1 and w_resid
# Maybe not quite different otherwise (__F I X M E__ merge ?? pbbly not a good idea)
# condition for either seems to be if(obsInfo) => $obs methods used (so requested, and non canonical)
.calc_etaLLMblob <- function(processed, muetablob, 
                             mu=muetablob$mu, eta=muetablob$sane_eta, 
                             old_beta_eta, ## scaled since X.pv is scaled; same for return beta_eta. An init.HLfit$fixef would be (.)/attr(spaMM:::.scale(zut$X.pv),"scaled:scale")
                             w.resid, # ignored!
                             phi_est, 
                             off=processed$off, 
                             maxit.mean, 
                             verbose, 
                             for_intervals=NULL,
                             Xtol_rel=processed$spaMM_tol$Xtol_rel) {
  BinomialDen <- processed$BinomialDen
  X.pv <- processed$AUGI0_ZX$X.pv
  y <- processed$y
  family <- processed$family
  LM_called <- FALSE
  damping <- 1e-7 ## as suggested by Madsen-Nielsen-Tingleff... # Smyth uses abs(mean(diag(XtWX)))/nvars
  newclik <- .calc_clik(mu=mu,phi_est=phi_est,processed=processed) ## handles the prior.weights from processed
  dlogL_blob <- .calc_dlogL_blob(eta, mu, y, weights=processed$prior.weights, family, phi=phi_est, muetaenv=muetablob,
                                 BinomialDen=BinomialDen)
  for (innerj in seq_len(maxit.mean)) {
    ## breaks when Xtol_rel is reached
    clik <- newclik
    # Historical oddity: the fit has worked with code which was OK for solving, but not for CI as the CI code suppresses 
    # a column of the design matrix, which is not sufficient on the premultiplied (scaled X) system.

    negHess <- crossprod(X.pv, .Dvec_times_m_Matrix( - dlogL_blob$d2logcLdeta2, X.pv))
    
    # names(szAug) <- colnames(X.pv) ## also important for intervalStep_glm
    ## simple QR solve with LevM fallback
    if ( ! is.null(for_intervals) || ! LM_called) {
      if ( ! is.null(for_intervals)) {
        currentDy <- (for_intervals$fixeflik-newclik)
        # The following check was previously performed at each iteration of etaxLM_fn: fixed-effect mean response but possibly mixed-effect residual response.
        # However, as long as the the random effects of the residual model have not converged (along with other estimates), likelihood calculations based of the predicted phi
        # are not appropriate for comparison. Hence one would have to condition the check on ! any(processed$models[["phi"]]=="phiHGLM")
        # But then the check is not very useful. It should be elsewhere, in cases where the v_h have converged but not other params.
        # = > REMOVE
        # if (currentDy < -1e-4 &&
        #     (is.null(bestlik <- processed$envir$confint_best$lik) || newclik > bestlik)) {
        #   if (is.null(bestlik)) {
        #     locmess <- paste("A higher",names(for_intervals$fixeflik),"was found than for the original fit.",
        #                      "\nThis suggests the original fit did not fully maximize",names(for_intervals$fixeflik),
        #                      "\nExpect more information at end of computation.")
        #     message(locmess)
        #   }
        #   processed$envir$confint_best$lik <- newclik
        #   processed$envir$confint_best$beta_eta <- .unscale(X.pv, old_beta_eta)
        # }
        intervalBlob <- .intervalStep_llm(old_beta=old_beta_eta, X.pv=X.pv, dlogcLdeta=dlogL_blob$dlogcLdeta, d2logcLdeta2=dlogL_blob$d2logcLdeta2, 
                                          currentlik=newclik, for_intervals=for_intervals,currentDy=currentDy) 
        beta_eta <- intervalBlob$beta
      } else if ( ! LM_called)  {
        beta_eta <- .llm_step_solver(dlogL_blob, eta, offset=off, X=X.pv, start=NULL, control=list(), damping=0L)
      }
      names(beta_eta) <- colnames(X.pv)
      # # PROBLEM is that NaN/Inf test does not catch all divergence cases so we need this :
      eta <- off + drop(X.pv %*% beta_eta) ## updated at each inner iteration
      muetablob <- .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed, phi_est=phi_est) 
      newclik <- .calc_clik(mu=muetablob$mu, phi_est=phi_est,processed=processed) 
    }  
    if ( is.null(for_intervals) &&
         (LM_called || # always use LevM when it has been called before: motivated by spaMM_glm example, cf email Alex, 01/06/2022, 14:13
          newclik < clik-1e-5 || anyNA(beta_eta) || any(is.infinite(beta_eta))) ) { 
      ## more robust LevM
      LM_called <- TRUE
      # a pure crossprod(*M*atrix) would return grad as a (1-col) Matrix, not well handled by .do_damped_WLS_llm
      grad <- .crossprod(X.pv, dlogL_blob$dlogcLdeta, eval_dens = FALSE, as_mat=TRUE) # H beta + grad
      
      damped_WLS_blob <- .do_damped_WLS_llm(grad, # 1-col *m*atrix
                                             damping, negHess, X.pv, clik, family, old_beta_eta, phi_est, off, processed, verbose)
      beta_eta <- damped_WLS_blob$beta 
      eta <- damped_WLS_blob$eta #off + drop(X.pv %*% beta_eta) ## updated at each inner iteration
      muetablob <- damped_WLS_blob$muetablob # .muetafn(eta=eta,BinomialDen=BinomialDen,processed=processed, phi_est=phi_est) 
      newclik <- damped_WLS_blob$clik
      damping <- damped_WLS_blob$damping
      dbetaV <- damped_WLS_blob$dbetaV
    } else dbetaV <- beta_eta - old_beta_eta
    dlogL_blob <- .calc_dlogL_blob(eta, mu=muetablob$mu, # muCOUNT 
                                   y, weights=processed$prior.weights, family, phi=phi_est, muetaenv=muetablob,
                                   BinomialDen=BinomialDen)
    if (verbose["trace"]) {
      print(paste0("Inner iteration ",innerj))
      print_err <- c(beta_eta=beta_eta)
      if (innerj>1L) print_err <- c(norm.dbetaV=sqrt(sum(dbetaV^2)),print_err)
      print(print_err)
      print("================================================")
    } 
    #if(innerj>600) browser()
    if (maxit.mean>1 && (damping>1e100 || mean(abs(dbetaV)) < Xtol_rel)) break
    #### done with one inner iteration
    old_beta_eta <- beta_eta
  } ## end for (innerj in 1:maxit.mean)
  names(beta_eta) <- colnames(X.pv)
  return(list(eta=muetablob$sane_eta, muetablob=muetablob, beta_eta=beta_eta, w.resid= - dlogL_blob$d2logcLdeta2, innerj=innerj,
              sXaug=structure(NA,info="no sXaug returned by .calc_etaLLMblob().") 
              ))
}
