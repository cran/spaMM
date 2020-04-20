.calc_denseness <- function(X) {
  if ( ! ( inherits(X,"ddiMatrix") || inherits(X,"dtCMatrix") || 
           inherits(X,"dsCMatrix") || inherits(X,"dgCMatrix")) ) X <- drop0(X) # result always Csparse
  if (inherits(X,"ddiMatrix")) {
    return(ncol(X))
  } else {
    res <- length(X@x) # this counts only once the symmetric elements
    if (methods::.hasSlot(X, "diag")) res <- res + ncol(X) # possible for dtC (but not dsc)
    return(res)
  }
}

.determine_spprec <- local({
  warned1 <- warned2 <- FALSE
  function(ZAlist, 
           processed, ## single envir, not list of envirs
           init.HLfit,
           #HLmethod, 
           nc = sum(unlist(lapply(ZAlist,ncol))),
           nr= nrow(ZAlist[[1]])
  ) {
    Xi_cols <- attr(ZAlist, "Xi_cols")
    nc <- sum(unlist(lapply(ZAlist,ncol)))
    anyRandomSlope <- any(Xi_cols>1L) ## FIXME seems oK for later code but semantically sloppy, cf (X-1|id) terms)
    if (any(sapply(processed$corr_info$corrMatrices,inherits,"precision"))) { ## derives from inherits(corrMatrix,"precision")
      if (anyRandomSlope && processed$For!="fitme" && !warned1) {
        warned1 <<- TRUE
        message(paste("Sparse-precision method was implied by the 'covStruct' or 'corrMatrix' argument,\n", 
                      "but random-coefficient terms may not yet be fitted efficiently by ",processed$For,"() in this case.\n", 
                      "fitme() may be more efficient for this combination of corrMatrix and random-coefficient terms."))
      } else sparse_precision <- TRUE ## force sparse
      # } else if (any( ! sapply(processed$corr_info$corrMatrices,is.null))) { ## tests whether it's explicit (we already know it is not of class  "precision")
      #   if (identical(spaMM.getOption("sparse_precision"),TRUE)) {
      #     message("Sparse-precision method was requested, but this will be ignored as an explicit correlation matrix was provided")
      #   }
      #   sparse_precision <- FALSE
    } else {
      sparse_precision <- spaMM.getOption("sparse_precision") ## global user control
      # if (identical(sparse_precision,TRUE) && anyRandomSlope && !warned2) {
      #   warned2 <<- TRUE
      #   message(paste("Sparse-precision method was requested by spaMM.options(),\n", 
      #                 "but random-coefficient terms may not yet be fit efficiently by such methods.\n", 
      #                 "Set spaMM.options(sparse_precision=FALSE) ?"))
      # }
    }
    # if (processed$For %in% c("HLfit","corrHLfit")) {
    #   if (anyRandomSlope) {
    #     if (identical(sparse_precision,TRUE)) {
    #       # FIXME the ollowing restriction is slightly too strict: "HLfit","corrHLfit" work with fixed ranCoefs but there's something missing
    #       # for estimating the ranCoefs.
    #       message(paste0("Sparse-precision method was requested, but random-coefficient terms\n", 
    #                     "are not yet handled by ",processed$For,"() in this case. _Another method will be used._\n",
    #                     "Use fitme() if you really want to fit random-coefficient terms by sparse-precision algorithms."))
    #     } 
    #     sparse_precision <- FALSE
    #   }
    # }
    ## best decision rule not obvious. For adjacency, trade off between repeated sparse Cholesky and a eigen() . 
    if (is.null(sparse_precision)) {
      exp_ranef_types <- attr(ZAlist,"exp_ranef_types")
      #which_is_adj <- which(exp_ranef_types %in% c("SAR_WWt","adjacency"))
      #any_adj <- length(which_is_adj)>0L
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
      if (is.null(sparse_precision)) {
        # which_is_IMRF <- which(exp_ranef_types %in% c("IMRF"))
        # nc_adj <- sum(unlist(lapply(ZAlist[which_is_adj],ncol)))
        # if (nc_adj>.spaMM.data$options$spprec_threshold) { # scotlip better by dgC... 
        #   sparse_precision <- TRUE
        # } else ... 
        any_IMRF <- any(exp_ranef_types== "IMRF")
        if (any_IMRF) {
          sparse_precision <- TRUE
        } else {
          ## all scotlip, and 'ohio' LMMs (using augZXy), slower by spprec (v3.1.45)
          #
          ## pb with nonspatial (including "(.|.)") as for any other model is that G may be dense, and the solve(BLOB$G_CHMfactor) are not efficient,
          #  (HLfit -> .get_hatvalues -> "hatval" is particularly inefficient and can take most of the time)
          # while the Matrix_QRP_CHM code is efficient on sparse sXaug 
          # ie when (L of corrMatrix) is sparse. So sparse precision should be used when there is a corrMatrix and
          # the precision matrix is significantly _sparser_ than the corrMatrix
          # Good examples: test-Rasch; 
          #                tnb <- fitme(resp~1+(1|ID), data=lll,family=Tnegbin(2))
          #                nested AR1 in test-AR1 is also slow by spprec
          # largeLMM <- fitme(...) nearly as fast by spprec as dgC.
          ## there are some cases where spprec is faster (sub-data of twin study)
          ## but it's difficult to generalize.
          G_diagnosis <- .provide_G_diagnosis(corr_info=processed$corr_info, ZAlist=ZAlist)
          ## trying to assess the cost of fitting AR using correlation-based algos:
          rel_ZAL_denseness <- (G_diagnosis$len_crossZC-G_diagnosis$len_noAR)/(nc^2) # rel_ZAL_denseness=0 for pure block effects
          ## special case of pure block effects : dgeMatrix is always OK in that case:
          ## test cases are: 
          ##    Rasch GLMM (rel_G_denseness=0.00938...) faster by dgeMatrix
          ##    LMM version or Rasch is fast in spprec (=> use Rasch GLMM as test case for seeking improvements in spprec?)
          ##    big-ranefs GLMM (rel_G_denseness=4.292e-05) faster by dgeMatrix; was MUCH faster by spprec until I modified .damping_to_solve() in 2.7.26
          ##    test-large-LMM (rel_G_denseness=0.03715023) equally fast by both methods
          # condition1 <- (length(ZAlist)>1L && (! anyRandomSlope) && all(exp_ranef_types=="(.|.)"))
          # condition1 <- condition1 && ((processed$LMMbool && rel_G_denseness<0.05) || rel_G_denseness < 0.0002)
          sparse_precision <- rel_ZAL_denseness*nc*nr>.spaMM.data$options$spprec_threshold ## from numerical experiments 
          # __F I X M E__ this is still poor, not finding that spprec is good for 
          #mrf <- fitme(migStatus ~ 1 + (1|pos) + multIMRF(1|latitude+longitude,margin=2,levels=2, coarse=4) better in dgCMatrix
          #checkinput <-  corrHLfit(... better in dgCMatrix
          #adjfit <- (HGLM) better in dgCMatrix
          # adjacency(1|gridcode) better in dgCMatrix
        }
      }
    } else if (sparse_precision && is.numeric(.getPar(init.HLfit,"rho"))) {
      stop("Conflicting option 'sparse_precision' and argument init.HLfit$rho")
    }
    return(sparse_precision)
  }
})
