.get_new_damping <- function(damping, which_LevMar_step) {
  if (which_LevMar_step=="strict_v|b") {
    damping <- "not used"
  } else if (damping>2e-7) {
    damping <- damping/(log2(damping)-log2(1e-7))
  } else damping <- 1e-7
  return(damping)
}

.calc_summand_gainratio <- function(processed, which_LevMar_step, LevMarblob, seq_n_u_h, ZAL_scaling, gainratio_grad) {
  if (processed$HL[1L]==1L) { ## ML fit 
    if (which_LevMar_step=="v_b") { 
      tempodvhbeta <- LevMarblob$dVscaled_beta
      tempodvhbeta[seq_n_u_h] <- tempodvhbeta[seq_n_u_h]*ZAL_scaling
      summand <- tempodvhbeta*(gainratio_grad+ LevMarblob$dampDpD * tempodvhbeta)
    } else if (which_LevMar_step=="b_&_v_in_b") {  # summand is computed after v_h ha been modified too
      # dbeta <- LevMarblob$dVscaled_beta[-seq_n_u_h]
      # summand <- dbeta*(gainratio_grad[-seq_n_u_h]+ LevMarblob$dampDpD * dbeta) 
      #### soame as for "v_b", *assuming* that we have put back the v_h changes in $dVscaled_beta
      tempodvhbeta <- LevMarblob$dVscaled_beta
      tempodvhbeta[seq_n_u_h] <- tempodvhbeta[seq_n_u_h]*ZAL_scaling
      summand <- tempodvhbeta*(gainratio_grad+ LevMarblob$dampDpD * tempodvhbeta)
    } else if (which_LevMar_step=="b") {
      summand <- LevMarblob$dbeta*(gainratio_grad[-seq_n_u_h]+ LevMarblob$dampDpD * LevMarblob$dbeta) 
    } else if (which_LevMar_step=="v") { ## v_h estimation given beta (FIXME can surely be made more exact)
      tempodvh <- LevMarblob$dVscaled*ZAL_scaling
      summand <- tempodvh*(gainratio_grad[seq_n_u_h]+ LevMarblob$dampDpD * tempodvh) 
    } else if (which_LevMar_step=="b_from_v_b") {
      dbeta <- LevMarblob$dVscaled_beta[-seq_n_u_h]
      summand <- dbeta*(gainratio_grad[-seq_n_u_h]+ LevMarblob$dampDpD * dbeta) 
    } else stop("Unhandled 'which_LevMar_step'")
  } else { ## joint hlik maximization
    summand <- LevMarblob$dVscaled_beta*(gainratio_grad+ LevMarblob$dampDpD * LevMarblob$dVscaled_beta) 
  }
}

.calc_summand_gainratio_spprec <- function(processed, which_LevMar_step, LevMarblob, seq_n_u_h, ZAL_scaling, gainratio_grad) {
  if (processed$HL[1L]==1L) { ## ML fit 
    if (which_LevMar_step=="v_b") { 
      tempodvhbeta <- unlist(LevMarblob[c("dVscaled","dbeta_eta")])
      tempodvhbeta[seq_n_u_h] <- tempodvhbeta[seq_n_u_h]*ZAL_scaling
      summand <- tempodvhbeta*(gainratio_grad+ LevMarblob$dampDpD * tempodvhbeta)
    } else if (which_LevMar_step=="b_&_v_in_b") {  # summand is computed after v_h ha been modified too
      # dbeta <- LevMarblob$dVscaled_beta[-seq_n_u_h]
      # summand <- dbeta*(gainratio_grad[-seq_n_u_h]+ LevMarblob$dampDpD * dbeta) 
      #### soame as for "v_b", *assuming* that we have put back the v_h changes in $dVscaled_beta
      tempodvhbeta <- unlist(LevMarblob[c("dVscaled","dbeta_eta")])
      tempodvhbeta[seq_n_u_h] <- tempodvhbeta[seq_n_u_h]*ZAL_scaling
      summand <- tempodvhbeta*(gainratio_grad+ LevMarblob$dampDpD * tempodvhbeta)
    } else if (which_LevMar_step=="b") {
      summand <- LevMarblob$dbeta_eta*(gainratio_grad[-seq_n_u_h]+ LevMarblob$dampDpD * LevMarblob$dbeta_eta) 
    } else if (which_LevMar_step=="v") { ## v_h estimation given beta (FIXME can surely be made more exact)
      tempodvh <- LevMarblob$dVscaled*ZAL_scaling
      summand <- tempodvh*(gainratio_grad[seq_n_u_h]+ LevMarblob$dampDpD * tempodvh) 
    } else if (which_LevMar_step=="b_from_v_b") {
      summand <- LevMarblob$dbeta_eta*(gainratio_grad[-seq_n_u_h]+ LevMarblob$dampDpD * LevMarblob$dbeta_eta) 
    } else stop("Unhandled 'which_LevMar_step'")
  } else { ## joint hlik maximization
    tempodvhbeta <- unlist(LevMarblob[c("dVscaled","dbeta_eta")])
    tempodvhbeta[seq_n_u_h] <- tempodvhbeta[seq_n_u_h]*ZAL_scaling
    summand <- tempodvhbeta*(gainratio_grad+ LevMarblob$dampDpD * tempodvhbeta)
  }
}

.wrap_do_damped_WLS_outer <- function(damped_WLS_fn, LevM_HL11_method, which_LevMar_step, old_relV_beta, constant_v_infer_args, 
                                      looseness, damping, rescue,
                                      ...) { ## I cannot always list(...) bc it contains promises that are not always defined
  loc_LevMar_step <- which_LevMar_step
  if (which_LevMar_step=="V_IN_B" || LevM_HL11_method[["b_step"]]=="v_in_b") { # i.e. we have previous inferred the need for the nested procedure
    ## each step of the damping loop  updates b and includes a v_h_IRLS
    loc_LevMar_step <- "b_&_v_in_b"
    v_infer_args <- constant_v_infer_args
    v_infer_args$looseness <- looseness  
  } else if (which_LevMar_step=="strict_v|b") { # i.e. we have previous inferred the need for the nested procedure
    ## damping <- Inf means that there is only one step of the damping loop which does not update b but includes a v_h_IRLS.
    v_infer_args <- constant_v_infer_args
    v_infer_args$looseness <- .spaMM.data$options$spaMM_tol$loose_resc  
    damping <- Inf
  } else v_infer_args <- NULL
  # run the damping loop in all cases
  damped_WLS_blob <- structure(
    damped_WLS_fn(v_infer_args=v_infer_args, which_LevMar_step=loc_LevMar_step, damping=damping, outer=TRUE,
                  stylefn=switch(loc_LevMar_step,
                                 v=.spaMM.data$options$stylefns$vloop,
                                 "strict_v|b"=.spaMM.data$options$stylefns$strictv,
                                 .spaMM.data$options$stylefns$betaloop ),
                  ...),
    step=which_LevMar_step
  )
  # optional *single* "v_in_b" IRLS for given beta REPLACES the previous one
  if (rescue && is.null(v_infer_args) && # loc_LevMar_step != "b_&_v_in_b" && ## not already "strict_v|b", V_IN_B or v_in_b
      ((breakcond <- damped_WLS_blob$breakcond)=="stuck_obj" || breakcond=="div_gain") ## suspect fit: if they are *both* low, the fit is presumably good
  ) { # coefficients do not move despite evidence that they should:
    # then we know that the heuristic procedure fails, and we will use the rigorous one;
    # But before that we perform a correct HL11 fit for the current beta (v_in_b truncated by damping=Inf argument):
    v_infer_args <- constant_v_infer_args
    v_infer_args$looseness <- .spaMM.data$options$spaMM_tol$loose_resc 
    ## this use the damped_WLS_fn function but in a way that does not call anny outer damping loop
    damped_WLS_blob <- 
      damped_WLS_fn(v_infer_args=v_infer_args, 
                    which_LevMar_step="strict_v|b", # for v_in_b ! # "v_b" neither handled nor really meaningful when damping=Inf
                    outer=TRUE,
                    damping=Inf, ## there are comparisons of damping to numerical values
                    stylefn=.spaMM.data$options$stylefns$rescue,
                    ...)
    attr(damped_WLS_blob,"step") <- "rescue"
  }
  return(damped_WLS_blob)
}

.calc_Xscal_newscaled <- function(newXscal, newZAL_scaling, ZAL, which_i_affected_rows, n_u_h, seq_n_u_h, processed) {
  if (TRUE) { ## alternative clause shows the meaning, but this version is distinctly faster. 
    scaledZAL <- .m_Matrix_times_Dvec(ZAL, newZAL_scaling)
    if (inherits(ZAL,"Matrix")) {
      # Next line assumes Xscal (IO_ZX) has no @diag component
      newXscal@x[which_i_affected_rows] <- scaledZAL@x ## should create an error if some elements are stored in @diag
    } else {
      newXscal[n_u_h+seq(nrow(scaledZAL)),seq_n_u_h] <- scaledZAL
    }
  } else newXscal <- .make_Xscal(ZAL, ZAL_scaling = newZAL_scaling, AUGI0_ZX=processed$AUGI0_ZX)
  return(newXscal) 
}

# checks and add attribute but does NOT re-evaluate the condition
.diagnose_coeff_not_moving <- function(coeff_not_moving, 
                                       relV_beta, # may actually be ofor only beta or only v_h
                                       damped_WLS_blob, innerj, damping, is_HL1_1, oldAPHLs, 
                                       Ftol, 
                                       trace, LevenbergM, stylefn) {
  info <- list()  
  if (is.na(coeff_not_moving)) { # diagnose the NA
    if (anyNA(relV_beta)) {
      if ( ! is.null(damped_WLS_blob)) {
        message(paste("innerj=",innerj,"damping=",damping,"lik=",damped_WLS_blob$lik))
        stop("Numerical problem despite Levenberg algorithm being used: complain.")
      } else stop("Numerical problem: try control.HLfit=list(LevenbergM=TRUE)")
    } else stop("Error in evaluating break condition")
  } 
  if (coeff_not_moving) { # DIAGNOSE this case
    # for ML, the test depends on values updated when ( hlik_stuck || ! need_v_step) and thus involves checks of the likelihood;
    # for PQL, we need a likelihood check:
    if ( ! is_HL1_1 ) { ## PQL in particular
      if ( ! is.null(damped_WLS_blob)) {
        info$hlik_stuck <- (damped_WLS_blob$APHLs$hlik < oldAPHLs$hlik + Ftol)
        if (info$hlik_stuck) { if (trace>1L) {cat(stylefn(" break bc 'coeff_not_moving' && 'hlik_stuck'\n"))} }
      } else {if (trace>1L) {
        info$hlik_stuck <- "Not assessed bc 'damped_WLS_blob' is NULL"
        cat(stylefn(" break bc 'coeff_not_moving' && 'is.null(damped_WLS_blob)'\n"))} 
      }
    } else {
      if ( ! is.null(damped_WLS_blob)) {
        # We end here if p_v is exactly not moving
        info$p_v_stuck <- (damped_WLS_blob$APHLs$p_v < oldAPHLs$p_v + Ftol) 
        if (info$p_v_stuck) {
          if (trace>1L) {cat(stylefn(" break bc 'coeff_not_moving' && 'p_v_stuck'\n"))} 
        }
      } else {
        info$p_v_stuck <- "Not assessed bc 'damped_WLS_blob' is NULL"
        if (trace>1L) {cat(stylefn("break bc 'coeff_not_moving'\n"))} ## sufficient condition in ML case 
      }
    }
  }
  return(structure(coeff_not_moving,info=info))
}

.wrap_v_h_IRLS <- function(v_h, beta_eta, seq_n_u_h, GLMMbool, wranefblob, 
                           constant_u_h_v_h_args, updateW_ranefS_constant_arglist, v_infer_args, Trace, IRLS_fn) {
  if (GLMMbool) {
    u_h <- v_h 
    newwranefblob <- wranefblob ## keep input wranefblob since GLMM and lambda_est not changed
  } else {
    u_h_v_h_from_v_h_args <- c(constant_u_h_v_h_args,list(v_h=v_h))
    u_h <- do.call(".u_h_v_h_from_v_h",u_h_v_h_from_v_h_args)
    if ( ! is.null(attr(u_h,"v_h"))) { ## second test = if u_h_info$upper.v_h or $lower.v_h non NULL
      v_h <- attr(u_h,"v_h")
    }
    ## update functions u_h,v_h
    newwranefblob <- do.call(".updateW_ranefS",c(updateW_ranefS_constant_arglist,list(u_h=u_h,v_h=v_h)))
  } 
  #
  v_infer_args$wranefblob <- newwranefblob
  v_infer_args$beta_eta <- beta_eta
  v_infer_args$v_h <- v_h
  v_infer_args$u_h <- u_h
  v_h_blob <- do.call(IRLS_fn, v_infer_args) ## blue + loop of damping loops -> loop of underline greens...
  return(v_h_blob)
}