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

.calc_G_denseness <- function(ZAlist) { # avoid cbinding the ZAlist
  cum_denseness <- 0L
  for (rd in seq_len(length(ZAlist))) {
    if (inherits(ZAlist[[rd]],"ddiMatrix")) { # no need to crossprod. Else, crossprod + .calc_denseness()
      cum_denseness <- cum_denseness + ncol(ZAlist[[rd]])
    } else cum_denseness <- cum_denseness + .calc_denseness(.crossprod(ZAlist[[rd]], allow_as_mat=FALSE))
    for (it in seq_len(rd-1L)) { # crossterms counting individuals that are jointed affected by ranef levels from two ranefs
      # .crossprod -> Matrix::crossprod() [or maybe some more complex approach]: ddi->ddi; dgC->dgC; dge->dpo (symmetric dense and storage doesn't use symmetry)
      cZA <- .crossprod(ZAlist[[rd]],ZAlist[[it]], allow_as_mat=FALSE) 
      cum_denseness <- cum_denseness + .calc_denseness(cZA) ## this counts only once the symmetric elements 
    }
  }
  return(cum_denseness)
}


.determine_spprec <- local({
  warned1 <- warned2 <- FALSE
  function(ZAlist, 
           processed, ## single envir, not list of envirs
           init.HLfit, 
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
        abs_G_denseness <- .calc_G_denseness(ZAlist) # length(crossprod(ZAfix)@x)
        # I could increase G_denseness by small terms for each of "AR1", "SAR_WWt", "adjacency"
        #   checkinput <- ... in test-AR1 is slower in sppprec and is run if G_thr <- 0.002
        # I could increase G_denseness by the denseness of the precision matrix for corrMatrix
        rel_G_denseness <- (abs_G_denseness-nc)/(nc^2) ## densness is relative do dim(G)=nc,nc
        exp_ranef_types <- attr(ZAlist,"exp_ranef_types")
        isMRFsl <- (exp_ranef_types %in% c("AR1", "SAR_WWt", "adjacency"
                                          # ,"(.|.)" # includes ranCoefs
                                          ))
        isIMRFss <- (exp_ranef_types=="IMRF")
        if (anyRandomSlope) {G_thr <- 0.0002} else  G_thr <- 0.0002
        condition1a <- ( any(isMRFsl) && all(isMRFsl | exp_ranef_types=="(.|.)") && ## allows "(.|.)" only if also MRFsl present
                           rel_G_denseness < G_thr
        ) ## FALSE for "Matern", "", etc.
        # alternative threshold when there is at least one IMRF, different threshold
        condition1b <- ( any(isIMRFss) && all(isMRFsl | isIMRFss | exp_ranef_types=="(.|.)") && ## allows "(.|.)" only if also IMRF present
                           rel_G_denseness < 0.05 # even with blackcap it's better to use spprec. 
                         # for IUCNI G_dens goes from 0.03+ to 0.02+ when  coarse from 4 to 10
        ) 
        condition1 <- (condition1a || condition1b)
      } else  {
        condition1 <- ( all(attr(ZAlist,"exp_ranef_types") %in% c("AR1", "SAR_WWt", "adjacency", "MFR")) && 
                          nr/nc>3 ## tall ZAfix
        ) ## FALSE for "Matern", "", etc.
      }
      sparse_precision <- condition1 && (
        (
          processed$LMMbool &&
            1e-4*(nc^2 -160^2)+1e-7*(nr^3 -250^3)>0 ## from numerical experiments
        ) || (
          # the Poisson-gamma adjacency HGLM in 'old donttest examples' is slow in sparse...
          0.5e-3*(nc^2 -120^2)+4.5e-8*(nr^3 -200^3)>0 ## from numerical experiments
        ) ## slow test_all without such conditions...
      )
    }
  } else if (sparse_precision && is.numeric(.getPar(init.HLfit,"rho"))) {
    stop("Conflicting option 'sparse_precision' and argument init.HLfit$rho")
  }
  return(sparse_precision)
}
})
