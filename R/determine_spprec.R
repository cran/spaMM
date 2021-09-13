.calc_denseness <- function(X, relative=FALSE) {
  if ( ! ( inherits(X,"ddiMatrix") || inherits(X,"dtCMatrix") || 
           inherits(X,"dsCMatrix") || inherits(X,"dgCMatrix")) ) X <- drop0(X) # result always Csparse
  if (inherits(X,"ddiMatrix")) {
    absolute <- ncol(X)
  } else {
    absolute <- length(X@x) # this counts only once the symmetric elements
    if (methods::.hasSlot(X, "diag")) absolute <- absolute + ncol(X) # possible for dtC (but not dsc)
  }
  if (relative) {
    return(absolute/prod(dim(X)))
  } else return(absolute)
}

.determine_spprec <- local({
  ######## pb is that solve(BLOB$G_CHMfactor....) is not efficient,
  ##  (HLfit -> .get_hatvalues -> "hatval" is particularly inefficient and can take most of the time)
  ## Ideally we should identify those cases.
  ## Special case of LMM pure block effects : correl algo is always OK in that case but spprec can be marginally better (augZXy case!)
  ## 'useful_for_det_spprec' in test-for-scaled-spprec is an AR1 fit by augZXy. It is fast by spprec or not.
  ## Ideally we should identify those cases. but .preprocess_augZXy() is evaluated long after spprec (_F I X M E_?)
  warned1 <- warned2 <- FALSE
  function(ZAlist, 
           processed, ## single envir, not list of envirs
           init.HLfit=processed$init_HLfit,
           X.pv,
           #HLmethod, 
           fast=TRUE,
           nc = sum(unlist(lapply(ZAlist,ncol))),
           nr= nrow(ZAlist[[1]])
  ) {
    
    if (processed$HL[1L]=="SEM") return(FALSE)
    # fast shortcut : ## presumably efficient use of Matrix::qr by .sXaug_Matrix_QRP_CHM_scaled algo
    if (inherits(X.pv,"sparseMatrix") && ncol(X.pv)> sum(unlist(lapply(ZAlist,ncol)))) return(FALSE)
    
    Xi_cols <- attr(ZAlist, "Xi_cols")
    nc <- sum(unlist(lapply(ZAlist,ncol)))
    any_corrMatrix_is_spprec <- (any(sapply(processed$corr_info$corrMatrices,inherits,"precision"))) ## derives from inherits(corrMatrix,"precision")
    if (any_corrMatrix_is_spprec) {
      anyRandomSlope <- any(Xi_cols>1L) ## FIXME seems oK for later code but semantically sloppy, cf (X-1|id) terms)
      if (anyRandomSlope && processed$For!="fitme" && !warned1) {
        warned1 <<- TRUE
        message(paste("Sparse-precision method was implied by the 'covStruct' or 'corrMatrix' argument,\n", 
                      "but random-coefficient terms may not yet be fitted efficiently by ",processed$For,"() in this case.\n", 
                      "fitme() may be more efficient for this combination of corrMatrix and random-coefficient terms."))
      } else sparse_precision <- TRUE ## force sparse
    } else {
      sparse_precision <- .spaMM.data$options$sparse_precision ## global user control
    }
    ## best decision rule not obvious. For adjacency, trade off between repeated sparse Cholesky and a eigen() . 
    if (is.null(sparse_precision)) {
      exp_ranef_types <- attr(ZAlist,"exp_ranef_types")
      any_adj <- any(exp_ranef_types %in% c("SAR_WWt","adjacency"))
      if (any_adj) {
        inner_estim_adj_rho <- ( processed$For=="HLCor" || (.get_cP_stuff(init.HLfit,"rho",count=TRUE)))
        if (inner_estim_adj_rho) sparse_precision <- FALSE
      # } else if (all(exp_ranef_types=="(.|.)") # && 
      #           # HLmethod=="ML(0,0,1)" # processed value for PQL/L
      #            ) {
        ## see comments below
      #   sparse_precision <- nr<3*nc
      } 
      # WHEN I CHANGE THE CODE HERE, I MUST CHECK THE CONTENTS OF help("sparse_precision")
      if (is.null(sparse_precision)) {
        ## mrf <- fitme(migStatus ~ 1 + (1|pos) + multIMRF(1|longitude+latitude,margin=2,levels=2, coarse=4)... better in spprec
        any_IMRF <- any(exp_ranef_types== "IMRF")
        if (any_IMRF) {
          if (FALSE) { 
            G_diagnosis <- .provide_G_diagnosis(corr_info=processed$corr_info, ZAlist=ZAlist, fast=FALSE)
            sparse_precision <-  with(G_diagnosis, (dens_G_rel_ZL<1 && density_G*dens_G_rel_ZL<0.05))
            if (FALSE && ! sparse_precision ) {
              cat(crayon::red(unlist(G_diagnosis)))
              cat(crayon::red(sparse_precision))
            }
          } else sparse_precision <- TRUE  # always for IMRF
        } else {
          G_diagnosis <- .provide_G_diagnosis(corr_info=processed$corr_info, ZAlist=ZAlist)
          if (G_diagnosis$fast) { # as determined internally by .provide_G_diagnosis()
            # actually no G diagnosis ; instead compares ZL to a ZL_without_AR, 
            # which amounts to assume that the cost of spprec is that of ZL without AR 
            rel_ZAL_denseness <- G_diagnosis$denseness_via_ZL/G_diagnosis$denseness_noAR 
            crit <-  rel_ZAL_denseness*nr/nc # can reach high values (e.g. adjacency-long > 400)
            sparse_precision <- crit >.spaMM.data$options$spprec_threshold ## from numerical experiments on ohio
          } else {
            rel_ZAL_denseness <- G_diagnosis$denseness_via_ZL/G_diagnosis$denseness_noAR 
            crit <-  rel_ZAL_denseness*nr/nc # can reach high values (e.g. adjacency-long > 400)
            sparse_precision <- crit >.spaMM.data$options$spprec_threshold ## from numerical experiments on ohio
            # Gryphon has second criterion below 4e-5 and covfit in test-adjacency-corrMatrix has it >2e-3
            sparse_precision <-  with(G_diagnosis, (dens_G_rel_ZL<1 && density_G*dens_G_rel_ZL<2e-4)) 
            if (FALSE) {
              cat(crayon::yellow(unlist(G_diagnosis)))
              cat(crayon::yellow(sparse_precision))
            }
          }
          if ( FALSE ) {
            old_rel_ZAL_denseness <- (G_diagnosis$denseness_via_ZL-G_diagnosis$denseness_noAR)/(nc^2) # rel_ZAL_denseness=0 for pure block effects
            old_sparse_precision <- old_rel_ZAL_denseness*nc*nr>5000 ## from numerical experiments
            if (sparse_precision!=old_sparse_precision){ # changed spprec
              cat(crayon::green(c(old_rel_ZAL_denseness,rel_ZAL_denseness)),"\n")
              cat(crayon::green(c(old_rel_ZAL_denseness*nc*nr,crit)),"\n")
              cat(crayon::green(c(old_sparse_precision,sparse_precision)),"\n")
            } else if (old_sparse_precision) { # always spprec: useful to see criteria
              cat(crayon::blue(c(old_rel_ZAL_denseness,rel_ZAL_denseness)),"\n")
              cat(crayon::blue(c(old_rel_ZAL_denseness*nc*nr,crit)),"\n")
              cat(crayon::blue(c(old_sparse_precision,sparse_precision)),"\n")
            }
          } 
        }
      }
    } else if (sparse_precision && is.numeric(.getPar(init.HLfit,"rho"))) {
      stop("Conflicting option 'sparse_precision' and argument init.HLfit$rho")
    }
    return(sparse_precision)
  }
})
