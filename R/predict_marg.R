.integrate_marg_lowup <- function(cum_nactive, sdvars, famfams, quantil) {
  
  len <- length(cum_nactive)-1L
  bounds <- rep(NA_real_,len)
  for (it in seq_len(len)) {
    u.range <- (cum_nactive[it]+1L):cum_nactive[it+1L]
    sdvar <- sdvars[[it]]
    bounds[u.range] <- switch(famfams[it],
                              gaussian = qnorm(quantil,sd=sdvar),
                              gamma = qgamma(quantil,shape=1/sdvar,scale=sdvar),
                              beta = qbeta(quantil,1/(2*sdvar),1/(2*sdvar)),
                              "inverse.gamma" = 1/rgamma(quantil,shape=1+1/sdvar,scale=sdvar), ## yields inverse gamma (1+1/object$lambda,1/object$lambda)
                              stop("(!) density for given rand.family not yet implemented")
    )
  }
  return(bounds)
}

.integrate_marg <- function(integrand, newZAXlist, eta_fix, sdvars, nactives, newinold, famfams, object, ii) {
  for (oldrd in seq_len(length(newZAXlist))) {
    if (oldrd %in% newinold) {
      active_levels <- which(newZAXlist[[oldrd]][ii,]!=0)
      newZAXlist[[oldrd]] <- newZAXlist[[oldrd]][ii,active_levels] ## allows drop !
      nactives[oldrd] <- length(active_levels)
    } 
  }
  newZAXrow <- do.call(c,newZAXlist) 
  cum_nactive <- cumsum(c(0,nactives))
  loq <- 1e-7
  for (it in seq_len(length(nactives))) {
    u.range <- (cum_nactive[it]+1L):cum_nactive[it+1L]
    lower <- .integrate_marg_lowup(cum_nactive, sdvars, famfams, quantil=loq)
    upper <- .integrate_marg_lowup(cum_nactive, sdvars, famfams, quantil=1-loq)
  }
  if (length(lower)>1L) {
    warning("Experimental code, will give poor results. Check the integration error in attr(.,'integrate_info')")
    # Vectorize vectorizes a fn(scalar), but does not "matricize" a fn(vector)
    integrand_v <- function(V, eta_fix, newZAXrow, cum_nactive, object) {
      matrix(apply(V, 2L, integrand, eta_fix=eta_fix, newZAXrow=newZAXrow, cum_nactive=cum_nactive, object=object), ncol = ncol(V))
    }
    locarglist <- list(f=integrand_v, lower=lower, upper=upper, method="hcubature", nVec=2L,
                       eta_fix=eta_fix, newZAXrow=newZAXrow, cum_nactive=cum_nactive, object=object)
      resu <- .do_call_wrap("cubintegrate",locarglist, pack="cubature")
  } else {
    integrand <- Vectorize(integrand,vectorize.args = "u")
    resu <- integrate(f=integrand,lower=lower, upper=upper, eta_fix=eta_fix, newZAXrow=newZAXrow, cum_nactive=cum_nactive,
                      object=object #, sdvars=sdvars, famfams=famfams, zero_truncated=zero_truncated
    ) 
    resu$call <- NULL ## call is problem for next line
  }
  data.frame(resu[]) ## for easy binding; [] removes the "integrate" class for which there is no data.frame method
}

# gives nonsense in first serious test case. Not sure why, bt integration error diverges.
# cf private test-predict-marginal.R
.predict_marg <- function(object, newdata, re.form, control) {
  if (is.null(newdata)) newdata <- object$data # quick way of forcing computation of $etaFix (FIXME)
  delayedAssign("invCov_oldLv_oldLv_list", .get_invColdoldList(object, control=control))
  new_X_ZACblob <- .calc_new_X_ZAC(object=object, newdata=newdata, re.form = re.form,
                                   variances=list(residVar=FALSE, cov=FALSE), invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list) 
  newinold <- new_X_ZACblob$newinold ## says which ranef is kept by re.form
  newZAXlist <- new_X_ZACblob$newZACpplist
  lcrandfamfam <- attr(object$rand.families,"lcrandfamfam")
  nrand <- length(newZAXlist)
  nactives <- rep(0L, nrand)
  sdvars <- vector("list", nrand)
  famfams <- rep( NA_character_, nrand)
  for (oldrd in seq_len(nrand)) {
    if (oldrd %in% newinold) {
      lmatrix <- object$strucList[[oldrd]]
      if (! is.null(lmatrix)) { ## spatial or random-coef
        newZAXlist[[oldrd]] <- t(solve(lmatrix,t(newZAXlist[[oldrd]])))  
      }
      famfams[oldrd] <- lcrandfamfam[oldrd]
      if (famfams[oldrd]=="gaussian") {
        sdvars[[oldrd]] <- sqrt(object$lambda.object$lambda_list[[oldrd]]) 
      } else sdvars[[oldrd]] <- object$lambda.object$lambda_list[[oldrd]]
    } else newZAXlist[[oldrd]] <- list(NULL)
  }
  eta_fix <- new_X_ZACblob$eta_fix
  #
  integrand <- function(u, eta_fix, newZAXrow, cum_nactive, object
                        #, sdvars, famfams, zero_truncated ## defined n envir where integrand is defined 
  ) { 
    v <- u_densities <- NA*u ## using u densities since integration variable is u
    for (it in seq_len(length(cum_nactive)-1L)) {
      u.range <- (cum_nactive[it]+1L):cum_nactive[it+1L]
      v[u.range] <- object$rand.families[[it]]$linkfun(u[u.range])
      sdvar <- sdvars[[it]]
      u_densities[u.range] <- switch(famfams[it],
                                     gaussian = dnorm(u[u.range],mean=0,sd=sdvar),
                                     gamma = dgamma(u[u.range],shape=1/sdvar,scale=sdvar),
                                     beta = rbeta(u[u.range],1/(2*sdvar),1/(2*sdvar)),
                                     "inverse.gamma" = 1/rgamma(u[u.range],shape=1+1/sdvar,scale=sdvar), ## yields inverse gamma (1+1/object$lambda,1/object$lambda)
                                     stop("(!) density for given rand.family not yet implemented")
      )
    }
    eta <- eta_fix + sum(newZAXrow*v)
    fv <- .fv_linkinv(eta=eta, family=object$family, families=object$families)
    return(fv * prod(u_densities)) ## predicted value|u * density(u)
  }
  #
  resu <- vector("list",length(eta_fix))
  for (ii in seq(nrow(newdata))) {
    resu[[ii]] <- .integrate_marg(integrand, newZAXlist, eta_fix=eta_fix[ii], sdvars=sdvars, 
                                  object=object, nactives=nactives, newinold=newinold, famfams=famfams, ii=ii)
  }
  resu <- do.call(rbind,resu)
  resu <- structure(matrix(resu[,1L],ncol=1L,dimnames=list(rownames(newdata),NULL)),
                    integrate_info=resu[,-1L], class=c("predictions","matrix"))
  return(resu)
}
