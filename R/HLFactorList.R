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
      ## without the next line, apply(mf[,splt] -> as.matrix(mf) produces artefacts such as space characters. 
      # the same artefacts should then be produced by seq_levelrange <- apply(uniqueGeo,1L,paste0,collapse=":")
      # mf[splt] <- lapply(mf[splt],factor)
      dataordered_levels <- apply(mf[splt],1,paste,collapse=":") ## paste gives a character vector, not a factor.
      
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
    # We need numeric values instead of factors in e_uniqueGeo. We use level labels (levels()) 
    # rather than level indices (1...n) there, so we should use the same here.
    # Without the next line of code, df[1,splt[-1L]]) would return the indices 1...n not the labels. 
    # Modifications of this code S H O U L D be tested on data where the nesting factors are coded as factors 
    # whose level names are  not 1...n (using subsetting from a large data frame is an appropriate way to 'mismatch' levels...)
    for (nam in splt[-1L]) {if (is.factor(fac <- mf[[nam]])) mf[[nam]] <- as.character(levels(fac))[fac]}
    # : see https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-an-integer-numeric-without-a-loss-of-information
    by_values <- by(mf,mf[,splt[-1L]],function(df) df[1,splt[-1L]]) ## (unique) combination of the nesting indices
    by_levels <- by(mf,mf[,splt[-1L]],getElement,name=splt[1L]) ## levels of nested (time) factor
    by_levelranges <- lapply(by_levels,range) 
    uniqueGeos <- vector("list",length(by_levels))
    AR1_block_n_u_h_s <- integer(length(by_levels))
    for (lit in seq_along(by_levels)) {
      seq_levelrange <- by_levelranges[[lit]][1L]:by_levelranges[[lit]][2L]
      uniqueGeos[[lit]] <- data.frame(seq_levelrange, by_values[[lit]])
      AR1_block_n_u_h_s[lit] <- diff(by_levelranges[[lit]])+1L
      #by_levels[[lit]] <- paste(prefixes[[lit]],by_levels[[lit]],sep=":")
    }
    uniqueGeo <- do.call(rbind,uniqueGeos)  ## data.frame (v2.3.9)
    seq_levelrange <- apply(uniqueGeo,1L,paste0,collapse=":")
  }
  colnames(uniqueGeo) <- splt
  res$uniqueGeo <- uniqueGeo ## more values than in the data ## .get_dist_nested_or_not() expects a *numeric* uniqueGeo, including cols for nesting factor
  ## we need seq_levelrange to construct the factor
  res$seq_levelrange <- seq_levelrange
  ## we always need AR1_block_n_u_h_s later 
  res$AR1_block_n_u_h_s <- AR1_block_n_u_h_s
  return(res)
}

.calc_Z_model_matrix <- function(leftOfBar_terms, leftOfBar_mf, raneftype) {
  modmat <- model.matrix(leftOfBar_terms, leftOfBar_mf) ## contrasts.arg not immed useful, maybe later.
  #if (raneftype == "(.|.)") stop("this does not occur") # does not seem to occur here
  if (! (is.null(raneftype))) {  ## NULL for ranCoefs! 
    if (ncol(modmat)>1L) { 
      # allowed: the variable was logical, or numeric (not factor, for which (fac|.) as well as (0+fac|.) generates cols for each level of the factor)
      # if numeric, should have used (0+z|.)
      classe <- attr(attr(leftOfBar_mf,"terms"),"dataClasses")
      if (classe=="logical") {
        ## TRUE/FALSE has created an intercept column...
        modmat <- modmat[,colnames(modmat) != "(Intercept)",drop=FALSE]
      } else if (classe=="factor") {
        ## TRUE/FALSE factor has created an intercept column... 0+. syntax not useful here => remove col and check again
        modmat <- modmat[,colnames(modmat) != "(Intercept)",drop=FALSE]
        if (ncol(modmat)>1L) stop(paste0("Unhandled expression in ", raneftype,"(<factor>|.):\n",
                                         " only TRUE/FALSE factor is allowed; '0 + <factor>' syntax is not allowed."))
      } else if (classe=="numeric") { ## true for integer variables  
        stop(paste0("Unhandled expression in ", raneftype,"(<numeric>|.): use explicit '0 + .' syntax to remove Intercept."))
      } else 
        stop(paste0("Unhandled expression in ", raneftype, "(<LHS>|.) for this type of random effect"))
    }
  }
  return(modmat)
}


.calc_Zmatrix <- local({
  trivial_incidMat <- sparseMatrix(i=1L,j=1L,x=1L, dimnames=list("1",NULL)) 
  function(x,mf,
           rmInt, ## remove Intercept
           drop,sparse_precision,type=".ULI",
           cov_mat_info, old_leftOfBar_mf=NULL) {
    ## le bidule suivant evalue le bout de formule x[[3]] et en fait un facteur. 
    ## but fac may be any vector returned by the evaluation of x[[3]] in the envir mf
    rhs <- x[[3]]
    txt <- .DEPARSE(rhs) ## should be the rhs of (|) cleanly converted to a string by terms(formula,data) in HLframes
    ## converts '%in%' to ':' 
    if (length(grep("%in%",txt))) {
      splittxt <- strsplit(txt,"%in%")[[1]]
      rhs <- as.formula(paste("~",splittxt[1],":",splittxt[2]))[[2]]
      txt <- .DEPARSE(rhs)
    }
    #    if (length(grep("\\+",txt))) { ## coordinates is a vector with a single string; grep is 1 if  + was found in this single string and numeric(0) otherwise
    ## if sparse_precision is not yet determined
    #  this build the design matrix as if it was TRUE,
    #  but adds the info dataordered_levels that allows a later modif of the design matrix if sparse_precision is set to FALSE
    AR1_sparse_Q <- info_mat_is_prec <- FALSE
    raneftype <- attr(x,"type")
    #if (identical(raneftype, "(.|.)")) stop("this does not occur") # does not occurs here, as explained in calling fn, .calc_Zlist()
    if (( ! is.null(raneftype))){ ## Any term with a 'spatial' keyword (incl. corrMatrix); cf comment in last case
      ## if sparse not yet determined for AR1, we generate the required info for sparse (and non-sparse) and thus assume AR1_sparse_Q: 
      if (is.null(AR1_sparse_Q <- sparse_precision)) AR1_sparse_Q <- (raneftype=="AR1")  
      info_mat_is_prec <- (raneftype=="corrMatrix" && inherits(cov_mat_info,"precision")) 
      ## for AR1_sparse and corMatrix, we cannot use dummy levels as created by .ULI() of factor(). THe level names have special meaning
      #   matching a time concept, or user-provided names for the corrMatrix
      if (raneftype %in% c("Matern","Cauchy", "MRF")) {
        # for MRF Z matches geo to uniqueGeo and A matches uniqueGeo to nodes
        dataordered_levels_blob <- .calc_dataordered_levels(txt=txt,mf=mf,type=type) ## even in sparse case
      } else if (AR1_sparse_Q || raneftype=="corrMatrix") {
        dataordered_levels_blob <- .calc_dataordered_levels(txt=txt,mf=mf,type="mf")
      } else dataordered_levels_blob <- .calc_dataordered_levels(txt=txt,mf=mf,type=type)
      #
      if (raneftype %in% c("Matern","Cauchy", "MRF")) {
        ff <- dataordered_levels_blob$factor ## so that Z cols will not be reordered.
      } else if (raneftype=="corrMatrix") {
        if (info_mat_is_prec) { 
          ff <- factor(dataordered_levels_blob$factor, levels=colnames(cov_mat_info$matrix))
        } else {
          ff <- dataordered_levels_blob$factor
        }
      } else if (AR1_sparse_Q) { 
        AR1_sparse_Q_ranges_blob <- .calc_AR1_sparse_Q_ranges(mf=mf,dataordered_levels_blob)
        ff <- factor(dataordered_levels_blob$factor,levels=AR1_sparse_Q_ranges_blob$seq_levelrange) ## rebuild a new factor with new levels
        if (anyNA(ff)) {
          stop(paste("Levels of the factor for an AR1 random effect should take integer values\n",
                     "(for convenient use of sparse-precision methods).")
          )
        }
      } else { # other raneftype's: handles for ( | ...+...) A N D importantly differs from the standard (.|.) code below,
        # which creates a Zt matrix with rows (then ZA cols) reordered as the automatic levels of the factor
        # while the cov mats / LMatrix has the original order
        # In particular im <- as(ff... creates a non-diagonal matrix in the he standard (.|.) code to represent this reordering.
        ff <- dataordered_levels_blob$factor
      }
    } else if (length(grep("c\\(\\w*\\)",txt))) { ## c(...,...) was used (actually detects ...c(...)....) (but in which context ?)
      aslocator <-  parse(text=gsub("c\\(", ".ULI(", txt)) ## slow pcq ULI() est slow
      ff <- as.factor(eval(expr=aslocator,envir=mf))
    } else { ## standard ( | ) rhs 
      mfloc <- mf
      ## automatically converts grouping variables to factors as in lme4::mkBlist (10/2014)
      for (i in all.vars(rhs)) { if ( ! is.null(curf <- mfloc[[i]])) mfloc[[i]] <- as.factor(curf)}
      if (is.null(ff <- tryCatch(eval(substitute(as.factor(fac), list(fac = rhs)), mfloc),
                                 error=function(e) NULL))) {
        message("couldn't evaluate grouping factor ", deparse(rhs)," within model frame:")
        if (length(grep("as.factor",rhs))) {
          stop("'as.factor' found in grouping factor term is not necessary and should be removed.",call.=FALSE)
        } else stop(" try adding grouping factor to data frame explicitly if possible",call.=FALSE)        
      }
      if (all(is.na(ff))) stop("Invalid grouping factor specification, ", deparse(rhs),call.=FALSE)
      ## note additional code in lme4::mkBlist for handling lhs in particular
    }
    ## Done with ff. Now the incidence matrix: 
    if (drop && ! (AR1_sparse_Q || info_mat_is_prec))  ff <- droplevels(ff)
    if (nrow(mf)==1L && levels(ff)=="1") {
      im <- trivial_incidMat ## massive time gain when optimizing spatial point predictions
    } else im <- sparseMatrix(i=as.integer(ff),j=seq(length(ff)),x=1L, # ~ as(ff, "sparseMatrix") except that empty levels are not dropped
                       dimnames=list(levels(ff),NULL)) # names important for corrMatrix case at least
    # : this is faster than   im <- Matrix::fac2sparse(ff,drop.unused.levels = (drop && ! (AR1_sparse_Q || info_mat_is_prec)))
    if (!isTRUE(methods::validObject(im, test = TRUE))) {
      stop("invalid conditioning factor in random effect: ", format(rhs))
    }
    ## model matrix for LHS in [...](LHS|rhs) (Intercept if ...(1|.)) 
    tempexp <- x[[2]] ## LHS
    leftOfBar_form <- eval(substitute(~expr, list(expr = tempexp)))
    leftOfBar_terms <- terms(leftOfBar_form) ## Implicitly assumes Intercept is included
    # old_leftOfBar_mf for prediction:
    if ( length(old_leftOfBar_mf)) { ## excludes NULL, or 0-col data.frames as in Matern(1|.) in *OLD*_leftOfBar_mf 
      # new predvars set on 'mf' by new_mf_ranef <- .calc_newFrames_ranef(.)$mf
      ori_levels <- .getXlevels(attr(old_leftOfBar_mf,"terms"), old_leftOfBar_mf) 
      leftOfBar_mf <- model.frame(leftOfBar_terms, mf, xlev = ori_levels) 
    } else leftOfBar_mf <- model.frame(leftOfBar_terms, mf) ## Matern(1|.) => [0 col; nrow=nrow(mf)]
    # note the test of contrasts on predict with ranCoefs with factors, in test-ranCoefs.R
    modmat <- .calc_Z_model_matrix(leftOfBar_terms, leftOfBar_mf, raneftype) ## handles non-trivial LHS in e.g. Matern(LHS|rhs)
    if (rmInt) { ## remove intercept column
      if ( ! is.na(icol <- match("(Intercept)", colnames(modmat)))) {
        if (ncol(modmat) < 2) stop("lhs of a random-effects term cannot be an intercept only")
        modmat <- modmat[, -icol, drop = FALSE]
      }
    }
    ## This is building Z(A) not Z(A)L hence reasonably sparse even in spatial models
    ZA <- .calc_raw_ZA(incidMat=im, modmat)
    attr(ZA,"leftOfBar_mf") <- leftOfBar_mf
    attr(ZA, "namesTerm") <- colnames(modmat) ## length=npar
    if (identical(raneftype,"AR1")) {
      if (AR1_sparse_Q) { ## this is TRUE is sparse_precision has not yet been determined !
        ## Following is different from levels(dataordered_levels_blob$factor) which are reordered as character
        #  Effect in first fit in test-AR1, when spprec goes from NULL to FALSE
        attr(ZA,"dataordered_unique_levels") <- unique(as.character(dataordered_levels_blob$factor)) ## allow reformatting for ! sparse prec
        colnames(ZA) <- AR1_sparse_Q_ranges_blob$seq_levelrange ## allow reformatting for ! sparse prec
        # ! ! ! caveat when changing the name of the following elements here, to change it elsewhere ! ! !
        attr(ZA,"AR1_block_n_u_h_s") <- AR1_sparse_Q_ranges_blob$AR1_block_n_u_h_s ## required for t_chol_Q computation
        attr(ZA,"uniqueGeo") <- AR1_sparse_Q_ranges_blob$uniqueGeo 
      } else {
        splt <- strsplit(txt,c("%in%|:|\\+| "))[[1L]]
        attr(ZA,"uniqueGeo") <- .calcUniqueGeo(data=mf[,splt,drop=FALSE])
      }
    } 
    ZA
  }
})

# incidMat is a incidence matrix with one 1 per COL (it is transposed relative to Z) (~perm but not necess square) hence its @x has nrow(im) elements
# its # of col is the number of levels of the grouping variable.
## modmat stores (<this info>| ...) => numeric for random slope model  
.calc_raw_ZA <- function(incidMat, modmat) {
  if (ncol(modmat)==1L && (length(umm <- unique(modmat[,1L]))==1L) && umm==1) { # classic (1|.) case
    if (nrow(incidMat)>1 && .is_identity(incidMat)) { ## if nrow=1 we may be optimizing point predictions and Diagonal is not worth the cost. 
      ZA <- Diagonal(n=nrow(incidMat))
      colnames(ZA) <- rownames(incidMat) 
    } else ZA <- t(incidMat)
    # incidMat has no colnames and modmat does not provide names in the alternative general code
  } else { ## first conceived for ranCoefs, with ncol(modmat)>1L. But also handles e.g. Matern(not-1|.) 
    ZA <- vector("list",ncol(modmat))
    for (col in seq_len(ncol(modmat))) {
      ZA_col <- incidMat     # don't try to fill with template <- t(incidMat) as @x would no longer have the correct order
      ZA_col@x <- modmat[,col]
      ZA_col <-  drop0(ZA_col)
      ZA[[col]] <- t(ZA_col)
    }
    ZA <- do.call(cbind, ZA) ## colnames are repeated if modmat has several cols...
    if (.spaMM.data$options$Zcolsbyrows) ZA <- ZA[,as.integer(matrix(seq(ncol(ZA)),byrow=TRUE,ncol=nrow(incidMat)))]
  }
  return(ZA)
}



.calc_Zlist <- function (formula, mf, rmInt, drop, sparse_precision=spaMM.getOption("sparse_precision"),
                         type=".ULI", corr_info, 
                         old_ZAlist=NULL,newinold=NULL, ## for prediction
                         barlist ## missing and NULL ust be distinguished 
) {
  ## drop=TRUE elimine des niveaux spurious (test: fit cbind(Mate,1-Mate)~1+(1|Female/Male) ....)
  ## avec des consequences ultimes sur tailles des objets dans dispGammaGLM
  exp_ranef_terms <- .process_bars(formula[[length(formula)]], 
                                   barlist=barlist, # remains missing in the simulate() call
                                   expand=TRUE, which. = "exp_ranef_terms") # which does not provide the "(.|.)" type
  if (!length(exp_ranef_terms)) return(structure(list(),anyRandomSlope=FALSE))
  x3 <- lapply(exp_ranef_terms, `[[`,i=3)
  names(exp_ranef_terms) <- unlist(lapply(x3, .DEPARSE)) ## names are RHS of (.|.)
  #######
  ZAlist <- vector("list",length(exp_ranef_terms))
  for (lit in seq_along(exp_ranef_terms)) {
    ZAlist[[lit]] <- .calc_Zmatrix(exp_ranef_terms[[lit]], mf=mf, rmInt=rmInt,
                                   drop=drop, sparse_precision=sparse_precision, type=type, 
                                   cov_mat_info=corr_info$corrMatrices[[lit]],
                                   old_leftOfBar_mf = attr(old_ZAlist[[newinold[lit]]],"leftOfBar_mf"))
    ## ALL ZAlist[[i]] are either 
    #     diagonal matrix (ddiMatrix) with @diag = "U"
    #  or (dg)*C*matrix ie a Compressed *C*olumn Storage (CCS) matrix 
    ##  (see http://netlib.org/linalg/html_templates/node92.html#SECTION00931200000000000000)
    ##  @x must contain the nonzero elements (except diagonal elements if @diag represents them)
    ##  @i contains the row indices of nonzero elements (except diagonal elements if @diag represents them)
    ##  @p[c] must contain the index _in x_ of the first nonzero element of column c, x[p[c]] in col c and row i[p[c]])  
  }
  ## Subject <- list(0) ## keep this as comment; see below
  namesTerms <- vector("list",length(ZAlist))
  GrpNames <- names(exp_ranef_terms)
  for (i in seq_len(length(ZAlist))) {
    ###################
    # Subject[[i]] <- as.factor(fl[[i]]$f) # levels of grouping var for all obs ('ff' returned by locfn)
    ## : Subject was used only for random slope model, where ncol(Design) != nlevels(Subject). I tried to get rid of this.
    ## see commented use of Subject in preprocess()
    ###################
    namesTerms[[i]] <- attr(ZAlist[[i]],"namesTerm") ## possibly several variables, eg intercept or slope... 
    names(namesTerms)[i] <- GrpNames[i] ## the name of the list member namesTerms[i]
  }
  ## One should not check .is_identity -> isDiagonal when 10000 points to predict... (FR->FR: modif def one of these functions ?)
  return(structure(ZAlist,  
                   exp_ranef_terms=exp_ranef_terms, ## matches ZAlist elements
                   exp_ranef_types=unlist(lapply(exp_ranef_terms,attr,which="type")), ## matches ZAlist elements
                   namesTerms=namesTerms, ## contains info for identifying random-coef terms
                   Xi_cols= unlist(lapply(namesTerms,length)))
  )
}
