.pot4improv <- function(which_LevMar_step, sXaug, gainratio_grad, seq_n_u_h) {
  switch(which_LevMar_step,
         "v_b"= get_from_MME(sXaug=sXaug, which="Mg_solve_g", B=gainratio_grad), # 0 at full solution (althoug B is not full gradient everywhere)
         "v"= get_from_MME(sXaug=sXaug, which="Mg_invH_g", B=gainratio_grad[seq_n_u_h]),  # 0 on the manifold
         "b"= { 
           B=gainratio_grad[-seq_n_u_h]
           if (length(B)) {
             get_from_MME(sXaug=sXaug, which="Mg_invXtWX_g", B=B) # hum not 0 zero at full solution (bc using only block of matrix) nor at hlik maximum (since B has full logdet gradient) 
           } else {0}
         }, 
         "b_&_v_in_b"= get_from_MME(sXaug=sXaug, which="Mg_solve_g", B=gainratio_grad), # again, 0 at full solution
         "b_from_v_b"= get_from_MME(sXaug=sXaug, which="Mg_solve_g", B=gainratio_grad) # again, 0 at full solution
  )
}

.loc_pot_tol <- function(which_LevMar_step, spaMM_tol) {
  switch(which_LevMar_step,
         "v_b"= spaMM_tol$b_pot_tol,
         "v"= spaMM_tol$v_pot_tol,
         "b"= spaMM_tol$b_pot_tol, 
         "b_&_v_in_b"= spaMM_tol$b_pot_tol,
         "b_from_v_b"= spaMM_tol$b_pot_tol
  )
}

.wrap_wrap_v_h_IRLS <- function(IRLS_fn, v_h, beta_eta, seq_n_u_h, GLMMbool, wranefblob, processed, updateW_ranefS_subarglist, v_infer_args, Trace) {
  #v_infer_args$maxit.mean <- ceiling(v_total_maxit_mean/5)
  v_h_blob <- .wrap_v_h_IRLS(v_h=v_h , 
                             beta_eta=beta_eta, seq_n_u_h, GLMMbool, wranefblob, 
                             processed$reserve$constant_u_h_v_h_args, updateW_ranefS_subarglist=updateW_ranefS_subarglist, 
                             v_infer_args, Trace, IRLS_fn=IRLS_fn)
  if (v_h_blob$break_info$IRLS_breakcond=="maxit") { # problematic failure: we need to do something
    #v_infer_args$maxit.mean <- ceiling(v_total_maxit_mean*4/5)
    psi_M <- rep(attr(processed$rand.families,"unique.psi_M"),diff(processed$cum_n_u_h))
    v_h_blob <- .wrap_v_h_IRLS(v_h=psi_M, 
                               beta_eta=beta_eta, seq_n_u_h, GLMMbool, wranefblob, 
                               processed$reserve$constant_u_h_v_h_args, updateW_ranefS_subarglist=updateW_ranefS_subarglist, 
                               v_infer_args, Trace, IRLS_fn=IRLS_fn)
  } 
  v_h_blob
}

.cat_break_info <- function(v_h_blob, stylefn_v, stylefn) {
  break_info <- v_h_blob$break_info
  cat(stylefn_v(paste0("v_h IRLS returns max(|grad|): v=",.prettify_num(break_info$maxs_grad[1L]), # grad when v_h IRLS exits
                       " beta=",.prettify_num(break_info$maxs_grad[2L]),
                       " after ",v_h_blob$innerj, # number of iterations of v_h IRLS
                       " iter"
  )))
  break_info$maxs_grad <- NULL
  for (st in names(break_info)) {
    if (is.numeric(stinfo <- break_info[[st]])) {
      cat(stylefn_v(paste0(", ",st,"=",.prettify_num(stinfo))))
    } else cat(stylefn_v(paste0(", ",st,"=",stinfo)))
  }
  cat(stylefn(";"))
}

.diagnose_conv_problem_LevM <- function(beta_cov_info, # __F I X M E__ redefine it to use tcrossfac_beta_v_cov rather than it tcrossprod? But RSpectra may be more efficient on symmetric matrices...
                                        w.resid, processed) {
  condnum <- NULL
  if ( ncol(beta_cov_info$tcrossfac_beta_v_cov)<2000L) { # for larger matrices the crossprod itself may be slow, perhaps the slowest step ?
    tc <- tcrossprod(beta_cov_info$tcrossfac_beta_v_cov) # hm. it's fairly dense
    if (inherits(tc,"sparseMatrix")) decomp <- .try_RSpectra(tc, symmetric=TRUE) # 1000 -> 0.28s
    if (is.null(decomp)) { # RSpectra was not available or it failed or matrix was not sparse
      if (ncol(tc)<1000L) { # 1000 -> 0.5s
        # kappa() computes by default (an estimate of) the 2-norm condition number of a matrix or of 
        # the R matrix of a QR decomposition, perhaps of a linear fit. The 2-norm condition number can 
        # be shown to be the ratio of the largest to the smallest *non-zero* singular value of the matrix.
        condnum <- kappa(tc)
      } # else condnum remains NULL
    } else condnum <- decomp$eigrange[2]/decomp$eigrange[1]
  } # else condnum remains NULL
  if ( ( ! is.null(condnum)) && condnum>1e08) {
    processed$envir$PQLdivinfo$high_kappa$ranFixes <- c(processed$envir$PQLdivinfo$high_kappa$ranFixes,
                                                        list(processed$envir$ranFix))
  } else { # problem too large, or no singularity
    # "unknown" : ideally this does not occur, if previous tests such a pot_tol checks provide explanations
    # But if this occurs, there are several other possible things to check:
    # check damping value ? provide info about hlik difference and about breakcond ?
    processed$envir$PQLdivinfo$unknown$ranFixes  <- c(processed$envir$PQLdivinfo$unknown$ranFixes,
                                                      list(processed$envir$ranFix))
  }
}

