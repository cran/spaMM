## alternative to lmerFactorList or so. Now see mkReTrms 

.calc_dataordered_levels <- function(txt,mf,type) {
  ## Standard (.|.) ranefs are NOT handled by this function but by code calling as.factor()
  if (type=="seq_len") { ## does not try to find redundant levels. Used by predict.HLfit() for spatial terms
    splt <- NULL
    dataordered_levels <- seq_len(nrow(mf))
  } else if (type==".ULI") { # spaMM's special ranef types that have no "nested nesting" : 
    # handles ( | ...+...) A N D importantly differs from the standard (.|.) code below,
    splt <- NULL
    aslocator <- parse(text=paste(".ULI(",gsub("\\+|:", "\\,", txt),")"))
    dataordered_levels <- eval(expr=aslocator,envir=mf) ## associates ordered levels 1...n to unique combin of rhs variables ordered as in the data.
  } else { ## handles "nested nesting" for AR1_sparse_Q || raneftype=="corrMatrix"
    splt <- strsplit(txt,c("%in%|:|\\+| "))[[1L]] ## things to be removed so that only variable names remain
    splt <- splt[splt!=""]
    if ( ! all(splt %in% names(mf)) ) stop(" ! all(splt %in% names(mf))")
    if (length(splt)==1L) {
      dataordered_levels <- mf[[splt[1L]]]  ## depending on the user, mf[[splt[1L]]] may be integer or factor...
    } else {
      dataordered_levels <- apply(mf[,splt],1,paste,collapse=":") ## paste gives a character vector, not a factor.
    } 
  } 
  return(list(factor=as.factor(dataordered_levels),splt=splt))
}

.calc_AR1_sparse_Q_ranges <- function(mf,dataordered_levels_blob) {
  res <- list()
  splt <- dataordered_levels_blob$splt
  if (length(splt)==1L) {
    levelrange <- range(as.integer(levels(dataordered_levels_blob$factor))) 
    AR1_block_n_u_h_s <- diff(levelrange)+1L
    uniqueGeo <- seq_levelrange <- levelrange[1L]:levelrange[2L]
    dim(uniqueGeo) <- c(length(uniqueGeo),1)
    #  this reorders levels differently to their order of appearance in the data, consistently with the code producing Q -> Lunique 
  } else {
    # We need integers instead of factors in e_uniqueGeo. We use level labels (levels()) 
    # rather than level indices (1...n) there so we should use the same here.
    # Without the next line of code, df[1,splt[-1L]]) would return the indices 1...n not the labels. 
    # Modifications of this code S H O U L D be tested on data where the nesting factors are coded as factors 
    # whose level names are  not 1...n (using subsetting from a large data frame is an appropriate way to 'mismatch' levels...)
    for (nam in splt[-1L]) {if (is.factor(fac <- mf[[nam]])) mf[[nam]] <- as.numeric(levels(fac))[fac]}
    # : see https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-an-integer-numeric-without-a-loss-of-information
    by_values <- by(mf,mf[,splt[-1L]],function(df) df[1,splt[-1L]]) ## (unique) combination of the nesting indices
    by_levels <- by(mf,mf[,splt[-1L]],getElement,name=splt[1L]) ## levels of nested (time) factor
    by_levelranges <- lapply(by_levels,range) 
    uniqueGeos <- vector("list",length(by_levels))
    AR1_block_n_u_h_s <- integer(length(by_levels))
    for (lit in seq_along(by_levels)) {
      seq_levelrange <- by_levelranges[[lit]][1L]:by_levelranges[[lit]][2L]
      uniqueGeos[[lit]] <- cbind(seq_levelrange,by_values[[lit]])
      AR1_block_n_u_h_s[lit] <- diff(by_levelranges[[lit]])+1L
      #by_levels[[lit]] <- paste(prefixes[[lit]],by_levels[[lit]],sep=":")
    }
    uniqueGeo <- do.call(rbind,uniqueGeos)
    seq_levelrange <- apply(uniqueGeo,1L,paste,sep="",collapse=":")
  }
  colnames(uniqueGeo) <- splt
  res$uniqueGeo <- uniqueGeo ## more values than in the data ## .get_dist_nested_or_not() expects a *numeric* uniqueGeo, including cols for nesting factor
  ## we need seq_levelrange to construct the factor
  res$seq_levelrange <- seq_levelrange
  ## we always need AR1_block_n_u_h_s later 
  res$AR1_block_n_u_h_s <- AR1_block_n_u_h_s
  return(res)
}

.spMMFactorList_locfn <- function(x,mf,rmInt,drop,sparse_precision,type=".ULI") {
  ## le bidule suivant evalue le bout de formule x[[3]] et en fait un facteur. 
  ## but fac may be any vector returned by the evaluation of x[[3]] in the envir mf
  rhs <- x[[3]]
  txt <- .DEPARSE(rhs) ## should be the rhs of (|) cleanly converted to a string by terms(formula,data) in HLframes
  ## converts '%in%' to ':' 
  if (length(grep("%in%",txt))>0) {
    splittxt <- strsplit(txt,"%in%")[[1]]
    rhs <- as.formula(paste("~",splittxt[1],":",splittxt[2]))[[2]]
    txt <- .DEPARSE(rhs)
  }
  #    if (length(grep("\\+",txt))>0) { ## coordinates is a vector with a single string; grep is 1 if  + was found in this single string and numeric(0) otherwise
  ## if sparse_precision is not yet determined
  #  this build the design matrix as if it was TRUE,
  #  but adds the info dataordered_levels that allows a later modif of the design matrix if sparse_precision is set to FALSE
  AR1_sparse_Q <- FALSE
  if (( ! is.null(raneftype <- attr(x,"type"))) && raneftype != "(.|.)"){ ## Any term with a 'spatial' keyword (incl. corrMatrix); cf comment in last case
    ## if sparse not yet determined for AR1, we generate the required info for sparse (and non-sparse) and thus assume AR1_sparse_Q: 
    if (nullsparse <- is.null(AR1_sparse_Q <- sparse_precision)) AR1_sparse_Q <- (raneftype=="AR1")  
    ## for AR1_sparse and corMatrix, we cannot use dummy levels as created by .ULI() of factor(). THe level names have special meaning
    #   matching a time concept, or user-provided names for the corrMatrix
    if (AR1_sparse_Q || raneftype=="corrMatrix") {
      dataordered_levels_blob <- .calc_dataordered_levels(txt=txt,mf=mf,type="mf")
    } else dataordered_levels_blob <- .calc_dataordered_levels(txt=txt,mf=mf,type=type)
    if (raneftype=="corrMatrix") {
      ff <- dataordered_levels_blob$factor
    } else if (AR1_sparse_Q) { 
      AR1_sparse_Q_ranges_blob <- .calc_AR1_sparse_Q_ranges(mf=mf,dataordered_levels_blob)
      ff <- factor(dataordered_levels_blob$factor,levels=AR1_sparse_Q_ranges_blob$seq_levelrange) ## rebuild a new factor with new levels
    } else { # other raneftype's: handles for ( | ...+...) A N D importantly differs from the standard (.|.) code below,
      # which creates a Zt matrix with rows (then ZA cols) reordered as the automatic levels of the factor
      # while the cov mats / LMatrix has the original order
      # In particular im <- as(ff... creates a non-diagonal matrix in the he standard (.|.) code to represent this reordering.
      ff <- dataordered_levels_blob$factor
    }
  } else if (length(grep("c\\(\\w*\\)",txt))>0) { ## c(...,...) was used (actually detects ...c(...)....) (but in which context ?)
    aslocator <-  parse(text=gsub("c\\(", ".ULI(", txt)) ## slow pcq ULI() est slow
    ff <- as.factor(eval(expr=aslocator,envir=mf))
  } else { ## standard ( | ) rhs 
    mfloc <- mf
    ## automatically converts grouping variables to factors as in lme4::mkBlist (10/2014)
    for (i in all.vars(rhs)) { if ( ! is.null(curf <- mfloc[[i]])) mfloc[[i]] <- as.factor(curf)}
    if (is.null(ff <- tryCatch(eval(substitute(as.factor(fac), list(fac = rhs)), mfloc),
                               error=function(e) NULL))) {
      message("couldn't evaluate grouping factor ", deparse(rhs)," within model frame:")
      if (length(grep("as.factor",rhs))>0) {
        stop("'as.factor' found in grouping factor term is not necessary and should be removed.",call.=FALSE)
      } else stop(" try adding grouping factor to data frame explicitly if possible",call.=FALSE)        
    }
    if (all(is.na(ff))) stop("Invalid grouping factor specification, ", deparse(rhs),call.=FALSE)
    if (drop) ff <- droplevels(ff)
    ## note additional code in lme4::mkBlist for handling lhs in particular
  }
  ## We have ff. 
  ## This is building Z(A) not Z(A)L hence reasonably sparse even in spatial models
  if (AR1_sparse_Q) {
    im <- sparseMatrix(i=as.integer(ff),j=seq(length(ff)),x=1L, # ~general code except that empty levels are not dropped
                       dimnames=list(levels(ff),NULL)) # names important for corrMatrix case at least
  } else im <- as(ff, "sparseMatrix") ##FR->FR slow step; but creates S4 objects with required slots
  if (!isTRUE(methods::validObject(im, test = TRUE))) {
    stop("invalid conditioning factor in random effect: ", format(rhs))
  }
  tempexp <- x[[2]] ## analyzing the LEFT hand side for non-trivial term (random-slope model)
  mm <- model.matrix(eval(substitute(~expr, list(expr = tempexp))), mf)
  if (rmInt) { ## remove intercept column
    if ( ! is.na(icol <- match("(Intercept)", colnames(mm)))) {
      if (ncol(mm) < 2) stop("lhs of a random-effects term cannot be an intercept only")
      mm <- mm[, -icol, drop = FALSE]
    }
  }
  ans <- list(f = ff, 
              A = do.call(rbind, rep(list(im),ncol(mm))),
              ## Zt is design obs <-> levels of ranef, either dgCMatrix (sparse) or dgeMatrix (dense)
              Zt = do.call(rbind, lapply(seq_len(ncol(mm)), .fillZtbloc, template=im, source=mm)),  ## mm stores (<this info>| ...) => numeric for random slope model 
              ST = matrix(0, ncol(mm), ncol(mm), dimnames = list(colnames(mm), 
                                                                 colnames(mm)))
              #lambda_X=mm ## design matrix for predictor of lambda
  )
  if (drop) {
    ans$A@x <- rep(0, length(ans$A@x))
    ans$Zt <- Matrix::drop0(ans$Zt)
  }
  if (identical(raneftype,"AR1")) {
    if (AR1_sparse_Q) { ## this is TRUE is sparse_precision has not yet been determined !
      ## Following is different from levels(dataordered_levels_blob$factor) which are reordered as character
      #  Effect in first fit in testèAR1, when spprec goes from NULL to FALSE
      ans$dataordered_unique_levels <- unique(as.character(dataordered_levels_blob$factor)) ## allow reformatting for ! sparse prec
      rownames(ans$Zt) <- AR1_sparse_Q_ranges_blob$seq_levelrange ## allow reformatting for ! sparse prec
      # ! ! ! caveat when changing the name of the following elements here, to change it elsewhere ! ! !
      ans$AR1_block_n_u_h_s <- AR1_sparse_Q_ranges_blob$AR1_block_n_u_h_s ## required for t_chol_Q computation
      ans$uniqueGeo <- AR1_sparse_Q_ranges_blob$uniqueGeo 
    } else {
      splt <- strsplit(txt,c("%in%|:|\\+| "))[[1L]]
      ans$uniqueGeo <- .calcUniqueGeo(data=mf[,splt,drop=FALSE])
    }
  }
  ans
}

.fillZtbloc <- function(col,source,template) {
  template@x <- source[,col]
  template
}


.spMMFactorList <- function (formula, mf, rmInt, drop, sparse_precision=spaMM.getOption("sparse_precision"),
                             type=".ULI") {
  ## drop=TRUE elimine des niveaux spurious (test: fit cbind(Mate,1-Mate)~1+(1|Female/Male) ....)
  ## avec des consequences ultimes sur tailles des objets dans dispGammaGLM
  ranef_terms <- .findbarsMM(formula[[length(formula)]])
  if (!length(ranef_terms)) return(structure(list(),anyRandomSlope=FALSE))
  exp_ranef_terms <- .spMMexpandSlash(ranef_terms) ## this is what expands a nested random effect
  ## ELSE
  x3 <- lapply(exp_ranef_terms, `[[`,i=3)
  names(exp_ranef_terms) <- unlist(lapply(x3, .DEPARSE)) ## names are RHS of (.|.)
  #######
  fl <- lapply(exp_ranef_terms, .spMMFactorList_locfn,mf=mf,rmInt=rmInt,drop=drop,sparse_precision=sparse_precision,type=type)
  ZAlist <- vector("list",length(fl))
  ## Subject <- list(0) ## keep this as comment; see below
  namesTerms <- vector("list",length(fl))
  GrpNames <- names(exp_ranef_terms)
  for (i in 1:length(fl)) {
    ###################
    # Subject[[i]] <- as.factor(fl[[i]]$f) # levels of grouping var for all obs ('ff' returned by locfn)
    ## : Subject was used only for random slope model, where ncol(Design) != nlevels(Subject). I tried to get rid of this.
    ## see commented use of Subject in preprocess()
    ###################
    Zt <- fl[[i]]$Zt
    ## ALL ZAlist[[i]] are (dg)*C*matrix ie a Compressed *C*olumn Storage (CCS) matrix 
    ##  (see http://netlib.org/linalg/html_templates/node92.html#SECTION00931200000000000000)
    ##  @x must contain the nonzero elements (except diagonal elements if @diag represents them)
    ##  @i contains the row indices of nonzero elements (except diagonal elements if @diag represents them)
    ##  @p[c] must contain the index _in x_ of the first nonzero element of column c, x[p[c]] in col c and row i[p[c]])  
    if (.is_identity(Zt)) {
      ZAlist[[i]] <- Diagonal(n=ncol(Zt)) ## diagonal matrix (ddiMatrix) with @diag = "U"
      colnames(ZAlist[[i]]) <- rownames(Zt) ## used e.g. in prediction (no more the case ?)
      rownames(ZAlist[[i]]) <- colnames(Zt) ## with transposition col/row names
    } else ZAlist[[i]] <- t(Zt) ## 
    attr(ZAlist[[i]],"dataordered_unique_levels") <- fl[[i]]$dataordered_unique_levels ## for conversion to nore standard, ! sparse-precision, matrices 
    attr(ZAlist[[i]],"AR1_block_n_u_h_s") <- fl[[i]]$AR1_block_n_u_h_s ## for building t_chol_Q
    attr(ZAlist[[i]],"uniqueGeo") <- fl[[i]]$uniqueGeo ## computed for AR1 only: predict, at last, will use it after copy in geo_info
    nt <- colnames(fl[[i]]$ST) ## length=npar
    namesTerms[[i]] <- nt ## possibly several variables, eg intercept or slope... 
    names(namesTerms)[i] <- GrpNames[i] ## the name of the list member namesTerms[i]
  }
  ## One should not check .is_identity -> isDiagonal when 10000 points to predict... (FR->FR: modif def one of these functions ?)
  return(structure(ZAlist,  
                   exp_ranef_terms=exp_ranef_terms, ## matches ZAlist elements
                   exp_ranef_types=unlist(lapply(exp_ranef_terms,attr,which="type")), ## matches ZAlist elements
                   # The following is not equivalent to ranef_strings produced by .parseBars(predictor),
                   #   as the latter returns strings  "<type>(.|.)" with () and optional type, instead of _terms_ .|. without () nor type
                   # ranef_terms=ranef_terms, ## matches formula terms
                   # ranef_types=unlist(lapply(ranef_terms,attr,which="type")), ## matches formula terms
                   namesTerms=namesTerms, ## contains info for identifying random-coef terms
                   Xi_cols= unlist(lapply(namesTerms,length)))
  )
}
