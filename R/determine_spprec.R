.determine_spprec <- function(ZAfix, ZAlist, 
                              processed, ## single envir, not list of envirs
                              init.HLfit) {
  Xi_cols <- attr(ZAlist, "Xi_cols")
  anyRandomSlope <- any(Xi_cols>1L) ## FIXME seems oK for later code but semantically sloppy, cf (X-1|id) terms)
  if (any(sapply(processed$corr_info$corrMatrices,inherits,"precision"))) { ## derives from inherits(corrMatrix,"precision")
    if (anyRandomSlope) {
      stop(paste("Sparse-precision method was implied by the 'covStruct' or 'corrMatrix' argument,\n", 
                 "but random-coefficient terms are not yet handled by ",processed$For,"() in this case.\n", 
                 "Use fitme() for this combination of corrMatrix and random-coefficient terms."))
    } else sparse_precision <- TRUE ## force sparse
  # } else if (any( ! sapply(processed$corr_info$corrMatrices,is.null))) { ## tests whether it's explicit (we already know it is not of class  "precision")
  #   if (identical(spaMM.getOption("sparse_precision"),TRUE)) {
  #     message("Sparse-precision method was requested, but this will be ignored as an explicit correlation matrix was provided")
  #   }
  #   sparse_precision <- FALSE
  } else {
    sparse_precision <- spaMM.getOption("sparse_precision") ## global user control
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
  ## best decision rule not obvious. For adjacency, trade off between repeated sparse Cholesky and a sym_eigen() . 
  if (is.null(sparse_precision)) {
    inner_estim_adj_rho <- (length(intersect(attr(ZAlist,"exp_ranef_types"),c("SAR_WWt","adjacency"))) && 
        ( processed$For=="HLCor" || (.get_cP_stuff(init.HLfit,"rho",count=TRUE))))
    if (inner_estim_adj_rho) {
      sparse_precision <- FALSE
    } else {
      ## pb with nonspatial (including "(.|.)") as for any other model is that G may be dense, and the solve(BLOB$G_CHMfactor) are not efficient,
      #  (HLfit -> .get_hatvalues -> "hatval" is particularly inefficient and can take most of the time)
      # while the Matrix_QRP_CHM code is efficient on sparse sXaug (cf test-Rasch as a good example)
      # ie when (L of corrMatrix) is sparse. So sparse precision should be used when there is a corrMatrix and
      # the precision matrix is significantly _sparser_ than the corrMatrix
      # (?fixme? ! not obvious from first attempt): rewrite part of sparse_precision code using solve(Q,... rather than solve(G,...)
      ## Examine density of Gmat
      if (TRUE) {
        G_denseness <- length(crossprod(ZAfix)@x)
        nc <- ncol(ZAfix)
        # I could increase G_denseness by small terms for each of "AR1", "SAR_WWt", "adjacency"
        # I could increase G_denseness by the denseness of the precision matrix for corrMatrix
        condition1 <- ( all(attr(ZAlist,"exp_ranef_types") %in% c("AR1", "SAR_WWt", "adjacency", "(.|.)")) && 
                          (G_denseness-nc)/(nc^2)< 0.0002
        ) ## FALSE for "Matern", "", etc.
      } else  {
        condition1 <- ( all(attr(ZAlist,"exp_ranef_types") %in% c("AR1", "SAR_WWt", "adjacency")) && 
                          nrow(ZAfix)/ncol(ZAfix)>3 ## tall ZAfix
        ) ## FALSE for "Matern", "", etc.
      }
      sparse_precision <- condition1 && (
        (
          processed$LMMbool &&
            1e-4*(ncol(ZAfix)^2 -160^2)+1e-7*(nrow(ZAfix)^3 -250^3)>0 ## from numerical experiments
        ) || (
          # the Poisson-gamma adjacency HGLM in 'old donttest examples' is slow in sparse...
          0.5e-3*(ncol(ZAfix)^2 -120^2)+4.5e-8*(nrow(ZAfix)^3 -200^3)>0 ## from numerical experiments
        ) ## slow test_all without such conditions...
      )
    }
  } else if (sparse_precision && is.numeric(.getPar(init.HLfit,"rho"))) {
    stop("Conflicting option 'sparse_precision' and argument init.HLfit$rho")
  }
  return(sparse_precision)
}

