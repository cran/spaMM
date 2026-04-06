.calc_denseness <- function(X, relative=FALSE) {
  #if (inherits(X,"dgeMatrix")) warning("Possibly inefficient use of 'dgeMatrix' caught by .calc_denseness.") 
  if ( ! inherits(X,c("ddiMatrix", "dtCMatrix", "dsCMatrix", "dgCMatrix")) ) 
    X <- drop0(X) # eg drop0(<dge>) result always Csparse 
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

.warn_once_not_fitme <- local({
  warned <- FALSE
  function(processed) {
    if (!warned) {
      warned <<- TRUE
      message(paste("Sparse-precision method was implied by the 'covStruct' or 'corrMatrix' argument,\n", 
                    "but random-coefficient terms may not yet be fitted efficiently by ",processed$For,"() in this case.\n", 
                    "fitme() may be more efficient for this combination of corrMatrix and random-coefficient terms."))
    }
  }
})

.determine_spprec_from_fast_G_diagnosis <- function(G_diagnosis, nc, nr, corr_info, ZAlist) {
  
  # The old crit uses cross_ZL to measure de **corr cost, presumably bc it was faster (?)
  dimZL <- G_diagnosis$dimZL
  denseness_cross_ZL <- G_diagnosis$rel_denseness_cross_ZL*dimZL[2]^2
  denseness_noAR <- G_diagnosis$rel_denseness_noAR*prod(dimZL)
  rel_ZAL_denseness <- denseness_cross_ZL/denseness_noAR 
  crit <-  rel_ZAL_denseness*(nr/nc)^(5/3)
  # tests for criterion: see test-determine_spprec.R
  sparse_precision <- crit >.spaMM.data$options$sp_alg_thresholds[["spprec"]] ## from numerical experiments on ohio
  sparse_precision <- sparse_precision || ( .spaMM.data$options$dec2spp &&
                                              # condition on nc a bit ad hoc but effective in # /covfit test
                                              # That case has rel_ZAL_denseness ~ 8 but nr/nc=0.5 
                                              nc>140L && 
                                              # .provide_QR_diagnosis() called within:
                                              .calc_default_QR_method(corr_info=corr_info, ZAlist=ZAlist, is_spprec=FALSE)=="dense"
                                            # note that .calc_default_QR_method() -> .provide_QR_diagnosis() -> .calc_QR_diagnosis()
                                            # is not called within .provide_G_diagnosis() -> .calc_G_diagnosis()
                                            # so it does not depend on a ZAlist locally .addrightcols_Z()-modified within .calc_G_diagnosis();
                                            # this seems OK: .addrightcols_Z() is applied on the final Z's only for spprec fits.
                                            # So we don't wan't is to be applied to the Z that would hypothetically be used for decorr/spcorr.
                                            # Further, .calc_default_QR_method() can use G_diagnosis...
  ) 
  sparse_precision
}

######## pb is that solve(BLOB$G_CHMfactor....) is not efficient,
##  (HLfit -> .get_hatvalues -> "hatval" is particularly inefficient and can take most of the time)
## Ideally we should identify those cases.
## Special case of LMM pure block effects : correl algo is always OK in that case but spprec can be marginally better (augZXy case!)
## 'useful_for_det_spprec' in test-for-scaled-spprec is an AR1 fit by augZXy. It is fast by spprec or not.
## Ideally we should identify those cases. but .preprocess_augZXy() is evaluated long after spprec (_F I X M E_?)
.determine_spprec <- function(ZAlist, 
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
  corr_info <- processed$corr_info
  any_corrMatrix_is_spprec <- (any(sapply(corr_info$corrMatrices,inherits,"precision")) ||
                                 any(.unlist(lapply(corr_info$corr_families,`[[`, "type"))=="precision")
  ) ## derives from inherits(corrMatrix,"precision")
  
  if (any_corrMatrix_is_spprec) {
    anyRandomSlope <- any(Xi_cols>1L) ## FIXME seems oK for later code but semantically sloppy, cf (X-1|id) terms)
    if (anyRandomSlope && processed$For!="fitme") {
      .warn_once_not_fitme(processed)
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
          G_diagnosis <- .provide_G_diagnosis(corr_info=corr_info, ZAlist=ZAlist)
          sparse_precision <-  with(G_diagnosis, (dens_G_rel_ZL<1 && density_G*dens_G_rel_ZL<0.05))
          if (FALSE && ! sparse_precision ) {
            cat(cli::col_red(unlist(G_diagnosis)))
            cat(cli::col_red(sparse_precision))
          }
        } else sparse_precision <- TRUE  # always for IMRF
      } else {
        G_diagnosis <- .provide_G_diagnosis(corr_info=corr_info, ZAlist=ZAlist)
        if (G_diagnosis$fast) {
          sparse_precision <- .determine_spprec_from_fast_G_diagnosis(
            G_diagnosis=G_diagnosis, nc=nc, nr=nr, corr_info=corr_info, ZAlist=ZAlist)
        } else {
          sparse_precision <- (names(which.min(G_diagnosis$costs))=="spprec")
        }
      }
    }
  } else if (sparse_precision && is.numeric(.getPar(init.HLfit,"rho"))) {
    stop("Conflicting option 'sparse_precision' and argument init.HLfit$rho")
  }
  return(sparse_precision)
}
