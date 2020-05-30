# derived from .calc_dvdlogphiMat_new()
# possible input:
# non-NULL lhs, no explicit rhs => rhs=t(lhs)
# non-NULL lhs and explicit rhs
# NULL lhs and explicit rhs => no premultiplication
.calc_lhs_inv_MD2hdv2_rhs <- function(lhs, rhs=t(lhs),
                                   sXaug,d2hdv2_info=NULL ## either one
) {
  if (is.null(d2hdv2_info)) { 
    inv_D2hdv2_rhs <- get_from_MME(sXaug,"solve_d2hdv2",B=rhs) 
  } else if (inherits(d2hdv2_info,"qr") || inherits(d2hdv2_info,"sparseQR") ) { # call on dense-prec fits may reach here. 
    inv_D2hdv2_rhs <- solve(d2hdv2_info, as.matrix(rhs))        
  } else if (is.environment(d2hdv2_info)) { # but in fact this is not called in post-fit code
    rhs <- d2hdv2_info$chol_Q %*% rhs
    rhs <- solve(d2hdv2_info$G_CHMfactor,rhs)
    inv_D2hdv2_rhs <- - .crossprod(d2hdv2_info$chol_Q,rhs) # inv_D2hdv2 in - t(chol_Q).inv(G).chol_Q
  } else { ## call on an spprec fit may reach here as .calc_d2hdv2_info -> -.crossprodCpp(as.matrix(factor_inv_Md2hdv2), yy = NULL) provides inv_d2hdv2 as a matrix.
    if (inherits(d2hdv2_info,"dsCMatrix")) {
      d2hdv2_info <- as(d2hdv2_info,"dgeMatrix") ## more efficient if inv_d2hdv2 is math-dense
    }
    inv_D2hdv2_rhs <- d2hdv2_info %*% rhs # is indeed inv_D2hdv2 %*% rhs     
  }
  if ( ! is.null(lhs)) {
    return( - lhs %id*id% inv_D2hdv2_rhs) 
  } else return( - inv_D2hdv2_rhs) # don't forget '-' for MD2hdv2
}


# V, section 13.2.2 de mes notes;  template was .calc_Evar; 
.calc_Var_given_fixef <- function(object, new_X_ZACblob, covMatrix, fix_X_ZAC.object=NULL) {
  newinold <- new_X_ZACblob$newinold # should be trivial since get_cPredVar() has no re.from argument
  if ( ! (newnrand <- length(newinold))) { ## fixed-effect model: no d2hdv2 hence no variance due to it
    nr <- nrow(new_X_ZACblob$newX.pv)
    if (covMatrix) {
      return(matrix(0,nrow=nr,ncol=nr))
    } else return(rep(0,nr))
  }
  
  cov_newLv_oldv_list <- new_X_ZACblob$cov_newLv_oldv_list
  cov_newLv_newLv_list <- new_X_ZACblob$cov_newLv_newLv_list
  diag_cov_newLv_newLv_list <- new_X_ZACblob$diag_cov_newLv_newLv_list
  newZAlist <- new_X_ZACblob$newZAlist # the one of the newdata
  ZAL <- get_ZALMatrix(object) # the one of the fitted data
  cov_fixLv_oldv_list <- fix_X_ZAC.object$cov_newLv_oldv_list
  fixZAlist <- fix_X_ZAC.object$newZAlist
  d2hdv2_info <- .calc_d2hdv2_info(object, ZAL) 
  Vlist <- vector("list",newnrand)
  for (new_rd in seq_along(Vlist)) {
    #old_rd <- newinold[new_rd]
    # In IMRF case there should be no new L (only a new A) 
    lhs <- newZAlist[[new_rd]] %ZA*gI% cov_newLv_oldv_list[[new_rd]]
    if (covMatrix) {
      if ( ! is.null(fixZAlist[[new_rd]])) { # ie rhs != t(lhs)
        rhs <- t(fixZAlist[[new_rd]] %ZA*gI% cov_fixLv_oldv_list[[new_rd]])
        terme <- .calc_lhs_inv_MD2hdv2_rhs(lhs=lhs, d2hdv2_info=d2hdv2_info, rhs=rhs)
      } else terme <- .calc_lhs_inv_MD2hdv2_rhs(lhs=lhs, d2hdv2_info=d2hdv2_info)
      Vlist[[new_rd]] <- as.matrix(terme)
    } else {
      if ( ! is.null(fixZAlist[[new_rd]])) { # ie rhs != t(lhs)
        rhs <- t(fixZAlist[[new_rd]] %ZA*gI% cov_fixLv_oldv_list[[new_rd]])
        terme <- rowSums(lhs * t(.calc_lhs_inv_MD2hdv2_rhs(lhs=NULL, d2hdv2_info=d2hdv2_info, rhs=rhs)))
      } else {
        rhs <- t(lhs)
        terme <- colSums(rhs * .calc_lhs_inv_MD2hdv2_rhs(lhs=NULL, d2hdv2_info=d2hdv2_info, rhs=rhs))
      }
      Vlist[[new_rd]] <- terme 
    }
  }
  V <- Reduce("+",Vlist)
  return(V)
}

get_cPredVar <- function(pred_object, newdata=NULL, nsim, seed=NULL, type="residual", variances=NULL, 
                         nb_cores=NULL, fit_env=NULL, sim_object=pred_object) {
  if ( requireNamespace("foreach", quietly = TRUE)) { # foreach is in Suggests
    nsim <- as.integer(nsim) 
    if (nsim<1L) {
      warning("'nsim' must be at least 1.")
      return(list())
    }
    sim_ys <- simulate(sim_object, type=type, nsim=nsim, seed=seed) # the bootstrap reproduces the source of error we wish to measure (cf cAIC)
    fit_call <- getCall(pred_object)
    fittingFunction <- paste(fit_call[[1]]) 
    ctrl_opt <- .update_control(fit_call=fit_call, optim_boot=.spaMM.data$options$optim_boot, from_fn=fittingFunction) # need .safe_opt when newinits are at bound.
    newinits <- get_inits_from_fit(from=pred_object) 
    names_not_cP <- setdiff(names(newinits$init),"corrPars")
    if (fittingFunction=="fitme") {
      new_args <- list(init=newinits$init[names_not_cP], init.HLfit=newinits$init.HLfit, control=ctrl_opt)
      if (! is.null(newinits$init[["corrPars"]])) new_args$fixed <- .modify_list(eval(fit_call$fixed),newinits$init["corrPars"])
    } else if (fittingFunction=="corrHLfit") {
      new_args <- list(init.corrHLfit=newinits$init.corrHLfit[names_not_cP], init.HLfit=newinits$init.HLfit, control.corrHLfit=ctrl_opt)
      if (! is.null(newinits$init.corrHLfit[["corrPars"]])) new_args$ranFix <- .modify_list(eval(fit_call$ranFix),newinits$init.corrHLfit["corrPars"])
    } else new_args <- list(init.HLfit=newinits$init.HLfit) ## could include HLCor with inner estimation of corrPars 
    simuland <- function(y, ...) { # always a dot args in the fn for dopar()
      hpff <- do.call(update_resp, c(list(object=pred_object, newresp = y),new_args)) 
      return(attr(predict(hpff,newdata=newdata, variances=list(naive=TRUE,cov=variances$cov)),"naive")) 
    }
    #
    bootreps <- dopar(newresp = sim_ys,nb_cores = nb_cores,fn = simuland, fit_env = fit_env)
    naive <- attr(predict(pred_object,newdata=newdata,variances=list(naive=TRUE,cov=variances$cov)),"naive") 
    bias <- rowMeans(bootreps)-naive
    SEs <- apply(bootreps,1L,sd)/sqrt(nsim)
    resu <- get_predVar(pred_object,newdata=newdata)-bias # implicitly
    return(structure(resu,info=list2env(list(SEs=SEs, bias=bias, naive=naive))))
  } else {stop(paste("'foreach' required but not available.",sep=""))}
}

