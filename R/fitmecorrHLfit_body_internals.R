.get_coordinates <- function(spatial.model, data) {
  if ( is.null(spatial.model)) {
    stop("An obsolete syntax for the adjacency model appears to be used.")
    ## coordinates <- c("x","y") ## backward compatibility... old syntax with (1|pos) and default values of the coordinates argument
  } else {
    if ( inherits(data,"list")) {
      dataForCheck <- data[[1]]
    } else dataForCheck <- data
    coordinates <- .extract_check_coords(spatial.model=spatial.model,datanames=names(dataForCheck))
  }
}

.get_cov_matrices__from_covStruct <- function(covStruct, corr.model, HLCor.args) {
  if ( ! is.null(covStruct)) covStruct <- .preprocess_covStruct(covStruct)
  resu <- list()
  if (corr.model  %in% c("corrMatrix")) {
    if (is.null(HLCor.args$corrMatrix)) resu$corrMatrix <- .get_corr_prec_from_covStruct(covStruct)
  } else if (corr.model  %in% c("SAR_WWt","adjacency")) { 
    ## adjMatrix should become weightMatrix
    if ( is.null(HLCor.args$adjMatrix) ) {
      resu$adjMatrix <- covStruct[["adjMatrix"]]
      if (is.null(resu$adjMatrix)) stop("missing 'adjMatrix' for adjacency model") ## or SAR_WWt model...
    }
  }
  if (length(resu)==0L) {
    return(NULL)
  } else return(resu)
}

.determine_sparse_precision <- function(processed, corr.model, init.HLfit) {
  if (.getProcessed(processed,"models",1L)[["eta"]]=="etaHGLM") {
    anyRandomSlope <- any(attr(.getProcessed(processed,"ZAlist",1L), "Xi_cols")>1L)
    sparse_precision <- spaMM.getOption("sparse_precision") ## user control
    if (is.null(sparse_precision)) {
      if (anyRandomSlope) sparse_precision <- FALSE
    } else if (sparse_precision && anyRandomSlope) {
      stop(paste("sparse precision was selected, but random-coefficient terms are not yet handled by sparse precision code."))
    }
    ## best decision rule not obvious. Trade off between repeated sparse Cholesky and a sym_eigen(). 
    if (is.null(sparse_precision)) {
      ZAfix <- .getProcessed(processed,"AUGI0_ZX$ZAfix",1L)
      if (corr.model %in% c("SAR_WWt","adjacency","AR1")) {
        sparse_precision <- (
          (
            ! is.numeric(init.HLfit$rho) ## init.HLfit$rho implies inner estim
          ) && (
            (
              .getProcessed(processed,"LMMbool",1L) 
              && 1e-4*(ncol(ZAfix)^2 -160^2)+1e-7*(nrow(ZAfix)^3 -250^3)>0 ## from numerical experiments
            ) || (
              # the Poisson-gamma adjacency HGLM in 'old donttest examples' i slow in sparse...
              0.5e-3*(ncol(ZAfix)^2 -120^2)+4.5e-8*(nrow(ZAfix)^3 -200^3)>0 ## from numerical experiments
            )
          )
        )  
      } else sparse_precision <- FALSE ## "Matern", "", etc.
      ## pb with nonspatial "" model is that G may be dense, and the solve(BLOB$G_CHMfactor) are not efficient,
      #  (HLfit -> .get_hatvalues -> "hatval" is particularly inefficient and can take most of the time)
      # while the Matrix_QRP_CHM code is efficient on sparse sXaug (cf test-Rasch as a good example)
      # ie when (L of corrMatrix) is sparse. So sparse precision should be used when 
      # the precision matrix is significantly _sparser_ than the corrMatrix
      # (?fixme? ! not obvious from first attempt): rewrite part of sparse_precision code using solve(Q,... rather than solve(G,...)
    } else if (sparse_precision && is.numeric(init.HLfit$rho)) {
      stop("Conflicting option 'sparse_precision' and argument init.HLfit$rho")
    }
  } else sparse_precision <- FALSE
  return(sparse_precision)
}

.provide_AR_factorization <- function(HLCor.args, sparse_precision, corr.model) {
  if (corr.model  %in% c("SAR_WWt")) {
    decomp <- eigen(HLCor.args$adjMatrix,symmetric=FALSE) ## FR->FR not RcppEigen-optimized
    return(list(u=decomp$vectors,d=decomp$values,u.=solve(decomp$vectors)))
  }
  # ELSE
  
  if (corr.model  %in% c("adjacency")) {
    ## eigenvalues needed in all cases for the bounds. Full decomp not always needed
    if (isSymmetric(HLCor.args$adjMatrix)) {
      if ( sparse_precision) { 
        decomp <- list(d=eigen(HLCor.args$adjMatrix,only.values = TRUE)$values) ## only eigenvalues
      } else {
        decomp <- sym_eigen(HLCor.args$adjMatrix)
      }
      return(decomp)
    } else stop("'adjMatrix' is not symmetric") ## => invalid cov mat for MVN
  }

}

.check_conflict_init_fixed <- function(fixed, init, errstring) {
  Fixrho <- .getPar(fixed,"rho")
  if ( (! is.null(.getPar(fixed,"nu"))) && (! is.null(init$nu)) ) {
    stop(paste("(!) 'nu'",errstring))    
  }
  if ( (! is.null(.getPar(fixed,"ARphi"))) && (! is.null(init$ARphi)) ) {
    stop(paste("(!) 'ARphi'",errstring))    
  }
  if ( (! is.null(.getPar(fixed,"Nugget"))) && (! is.null(init$Nugget)) ) {
    stop(paste("(!) 'Nugget'",errstring))    
  }
  if ( (! is.null(Fixrho)) && (! is.null(init$rho)) ) {
    stop(paste("(!) 'rho'",errstring))    
  } else return(max(length(Fixrho),length(init$rho)))
}


.provide_rho_mapping <- function(control.dist, coordinates, rho.size) {
  rho_mapping <- control.dist$rho.mapping ## may be NULL
  if (is.null(rho_mapping) ) { ## && length(coordinates)>1L ?
    if (length(coordinates)==rho.size) { ## rho.size comes from explicit rho from user
      rho_mapping <- seq_len(rho.size)           
      names(rho_mapping) <- coordinates
    } else if (length(rho.size)>1L) stop("'rho.mapping' missing with no obvious default from the other arguments.")
  } ## then (for given corr.model's) there is rho_mapping
  return(rho_mapping)
}

.calc_maxrange <- function(rho.size,distMatrix=NULL,uniqueGeo=NULL,rho_mapping,dist.method) {
  if (rho.size<2L) { ## can be 0 if no explicit rho in the input  
    if (inherits(distMatrix,"list")) {
      maxrange <- max(unlist(lapply(distMatrix, function(dd) max(c(-Inf,dd))))) ## les Inf to handle dist(0)...
      -min(unlist(lapply(distMatrix, function(dd) min(c( Inf,dd)))))
    } else maxrange <- max(distMatrix)-min(distMatrix)   
    if (maxrange<.Machine$double.eps) stop("Only one distance value: spatial model cannot be fitted.")
  } else { ## rho.size >1
    if (inherits(uniqueGeo,"list")) {
      maxrange <- lapply(unique(rho_mapping), function(idx) {
        ranges <- matrix(unlist(lapply(uniqueGeo, function(uu){
          if (nrow(uu)>1) {
            range(proxy::dist(uu[,rho_mapping==idx],method=dist.method))
          } else c(Inf,-Inf) ## encore des Inf to handle dist(0)...
        })),ncol=2)
        max(ranges[,2])-min(ranges[,1]) 
      })
    } else { ## single data set
      maxrange <- lapply(unique(rho_mapping), function(idx) {
        rng <- range(proxy::dist(uniqueGeo[,rho_mapping==idx],method=dist.method))
        rng[2]-rng[1] 
      })
    }  
    maxrange <- unlist(maxrange)
  }
  return(maxrange)
}

.do_TRACE <- function(level) {
  if (level) {
    suppressMessages(trace(HLfit_body,print=FALSE, 
                           tracer=quote(try({
                             rpblob <- .canonizeRanPars(ranFix,corr.model="",checkComplete = FALSE)
                             trueCorrpars <- rpblob$trueCorrpars
                             ntC <- names(trueCorrpars)
                             for (lit in seq_along(trueCorrpars)) cat(ntC[lit],"= ",paste(signif(trueCorrpars[[lit]],6),collapse=" ")," ")
                             ranPars <- rpblob$ranPars
                             ntC <- names(ranPars)
                             for (lit in seq_along(ranPars)) cat(ntC[lit],"= ",paste(signif(unlist(ranPars[[lit]]),6),collapse=" ")," ")
                           })),
                           exit=quote({
                             aphl <- unlist(res$APHLs[c("p_bv","p_v")])[1L] ## unlist drops a NULL p_bv
                             print(paste(names(aphl),"= ",signif(aphl,6),sep=""),quote=FALSE)
                           }))) # FIXME should print p_bv in REML fits ?
    suppressMessages (trace(.solve_IRLS_as_ZX,print=FALSE,tracer=quote(cat(">"))))
    suppressMessages (trace(.solve_IRLS_as_spprec,print=FALSE,tracer=quote(cat(">"))))
    suppressMessages (trace(spaMM.getOption("matrix_method"),print=FALSE,tracer=quote(cat("."))))
    suppressMessages(trace(spaMM.getOption("Matrix_method"),print=FALSE,tracer=quote(cat("."))))
    suppressMessages(trace("def_AUGI0_ZX_sparsePrecision",print=FALSE,tracer=quote(cat("."))))
    if (level>3L) {
      fn <- paste("get_from_MME",strsplit(spaMM.getOption("matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_EigenDense_QRP_Chol_scaled
      suppressMessages(trace(fn,print=FALSE,tracer=quote(cat(which))))
      fn <- paste("get_from_MME",strsplit(spaMM.getOption("Matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_Matrix_QRP_CHM_scaled
      suppressMessages(trace("fn",print=FALSE,tracer=quote(cat(which))))
      suppressMessages(trace("get_from_MME.AUGI0_ZX_sparsePrecision",print=FALSE,tracer=quote(cat(which))))
    }
  } else {
    suppressMessages(untrace(spaMM::HLfit_body))
    suppressMessages(try(untrace(.solve_IRLS_as_ZX),silent=TRUE)) ## untracing untraced internal functiosn fails
    suppressMessages(try(untrace(.solve_IRLS_as_spprec),silent=TRUE))
    suppressMessages(untrace(spaMM.getOption("matrix_method")))
    suppressMessages(untrace(spaMM.getOption("Matrix_method")))
    suppressMessages(untrace("def_AUGI0_ZX_sparsePrecision"))
    fn <- paste("get_from_MME",strsplit(spaMM.getOption("matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_EigenDense_QRP_Chol_scaled
    suppressMessages(untrace(fn))
    fn <- paste("get_from_MME",strsplit(spaMM.getOption("Matrix_method"),"def_")[[1L]][2],sep=".") ## get_from_MME.sXaug_Matrix_QRP_CHM_scaled
    suppressMessages(untrace(fn))
    suppressMessages(untrace("get_from_MME.AUGI0_ZX_sparsePrecision"))
  } 
} 

.init_precision_info <- function(processed) {
  envir <- processed$AUGI0_ZX$envir
  nranef <- length(envir$types)
  precisionFactorList <- vector("list",nranef) ## will contain diagonal matrices/info for non-trivial (diagonal) precision matrices 
  for (it in seq_len(nranef)) {
    if ( envir$types[it] %in% c("adjacency","corrMatrix","AR1") ) {
      ## terms of these types must have been dealt with by ad hoc code for each type
    } else if ( envir$types[it] %in% c("(.|.)") ) {
      nc <- ncol(processed$ZAlist[[it]])
      precisionFactorList[[it]] <- list(Qmat=Diagonal(n=nc), ## used independently of chol_Q_list, see precisionBlocks
                                        chol_Q=new("dtCMatrix",i= 0:(nc-1L), p=0:(nc), Dim=c(nc,nc),x=rep(1,nc)) )
      # All chol_Q's must be dtCMatrix so that bdiag() gives a dtCMatrix
    } else stop(paste("sparse precision was selected, but",envir$types[it],"terms are not yet handled by sparse precision code."))
    # Matern fails with Error in solve(Q_CHMfactor, system = "Lt") : object 'Q_CHMfactor' not found 
  }
  envir$precisionFactorList <- precisionFactorList
  ## environment modified, no return value
}

