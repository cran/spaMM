.determine_spprec <- local({
  warned1 <- warned2 <- FALSE
  function(ZAfix, ZAlist, 
                              processed, ## single envir, not list of envirs
                              init.HLfit) {
  Xi_cols <- attr(ZAlist, "Xi_cols")
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
      #
      # another stricking example is the test with tnb <- fitme(resp~1+(1|ID), data=lll,family=Tnegbin(2))
      # nested AR1 in test-AR1 is also slow by spprec
      #
      # ie when (L of corrMatrix) is sparse. So sparse precision should be used when there is a corrMatrix and
      # the precision matrix is significantly _sparser_ than the corrMatrix
      ## Examine density of Gmat
      if (TRUE) {
        G_denseness <- length(crossprod(ZAfix)@x)
        nc <- ncol(ZAfix)
        # I could increase G_denseness by small terms for each of "AR1", "SAR_WWt", "adjacency"
        # I could increase G_denseness by the denseness of the precision matrix for corrMatrix
        if (anyRandomSlope) {denseness_thr <- 0.0002} else  denseness_thr <- 0.0002
        exp_ranef_types <- attr(ZAlist,"exp_ranef_types")
        isMRFsl <- (exp_ranef_types %in% c("AR1", "SAR_WWt", "adjacency"
                                          # ,"(.|.)" # includes ranCoefs
                                          ))
        isMRFss <- (exp_ranef_types=="MRF")
        
        condition1a <- ( any(isMRFsl) && all(isMRFsl | exp_ranef_types=="(.|.)") && ## allows "(.|.)" only if also MRF present
                           (G_denseness-nc)/(nc^2)< denseness_thr
        ) ## FALSE for "Matern", "", etc.
        # alternative thrshold when there is at least one MRF, different threshold
        condition1b <- ( any(isMRFss) && all(isMRFsl | isMRFss | exp_ranef_types=="(.|.)") && ## allows "(.|.)" only if also MRF present
                           (G_denseness-nc)/(nc^2)< 0.02
        ) 
        condition1 <- (condition1a || condition1b)
      } else  {
        condition1 <- ( all(attr(ZAlist,"exp_ranef_types") %in% c("AR1", "SAR_WWt", "adjacency", "MFR")) && 
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
})
