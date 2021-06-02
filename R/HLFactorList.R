## alternative to lmerFactorList or so. Now see mkReTrms 

.as_factor <- function(txt,mf,type) {
  ## Standard (.|.) ranefs are NOT handled by this function but by code calling as.factor()
  if (type=="seq_len") { ## does not try to find redundant levels. Used by predict.HLfit() for spatial terms
    splt <- NULL
    raw_levels <- seq_len(nrow(mf))
  } else if (type==".ULI") { # has been the default type of .calc_Zlist(). Now .ULI() no longer called through this case. 
    # handles ( | ...+...) but not nested groups, while the general alternative below may handle both (at least separately)
    # Used to make sure that results are data_ordered, but type="data_order" may now handle that. 
    splt <- NULL
    aslocator <- parse(text=paste(".ULI(",gsub("\\+|:", "\\,", txt),")"))
    raw_levels <- eval(expr=aslocator,envir=mf) ## associates ordered levels 1...n to unique combin of rhs variables ordered as in the data.
  } else { ## all other types such as "mf" or "data_order"; handles "nested nesting" for AR1 spprec || raneftype=="corrMatrix"
    splt <- strsplit(txt,c("%in%|:|\\+| "))[[1L]] ## things to be removed so that only variable names remain
    splt <- splt[splt!=""]
    if ( ! all(splt %in% names(mf)) ) stop(" ! all(splt %in% names(mf))")
    if (length(splt)==1L) {
      raw_levels <- mf[[splt[1L]]]  ## depending on the user, mf[[splt[1L]]] may be integer or factor...
    } else {
      #raw_levels <- apply(mf[splt],1,paste,collapse=":") ## paste gives a character vector, not a factor.
      # far-fetched code to the same effect (from plotrix::pasteCols)
      x <- t(mf[splt])
      pastestring <- paste("list(",paste("x","[",seq_len(nrow(x)),",]",sep="",collapse=","),")",sep="")
      raw_levels <- do.call(paste,c(eval(parse(text = pastestring)),sep=":"))
    } 
  } 
  if (type=="data_order") {
    return(list(factor=factor(raw_levels,levels=unique(raw_levels)),splt=splt))
  } else return(list(factor=as.factor(raw_levels),splt=splt)) # ".ULI" already data_ordered
}

.calc_AR1_sparse_Q_ranges <- function(mf,levels_blob) {
  res <- list()
  splt <- levels_blob$splt
  if (length(splt)==1L) {
    levelrange <- range(as.integer(levels(levels_blob$factor))) 
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

.calc_Z_model_matrix <- function(leftOfBar_terms, leftOfBar_mf, raneftype,lcrandfamfam) {
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
  } else if (ncol(leftOfBar_mf)==1L) { ## ncol=0L is for (1|.) ## single variable, but modmat may have an intercept col
    if (lcrandfamfam != "gaussian" && 
        (attr(attr(leftOfBar_mf,"terms"),"dataClasses"))=="numeric"
    ) { ## Gamma(wei|.)
      if (ncol(modmat)>1L) stop(paste0("Unhandled expression in ", raneftype,"(<numeric>|.): use explicit '0 + .' syntax to remove Intercept."))
      prior_lam_fac <- modmat[,1]^2
      modmat[] <- 1
      attr(modmat,"prior_lam_fac") <- prior_lam_fac
    }
  } 
  return(modmat)
}

.add_levels <- function(ff, adj_or_prec, old_ZA) {
  if (is.null(old_ZA)) {
    RHS_levels <- colnames(adj_or_prec) # may well be NULL (=> later error only if automatic matching of ZA and 'L' is not possible)
  } else RHS_levels <- colnames(old_ZA)
  if ( ! is.null(RHS_levels)) {  # factor() does not satisfactorily handle levels=NULL
    ff <- factor(ff, levels=RHS_levels)
  }
  ff
}

.calc_Zmatrix <- local({
  trivial_incidMat <- sparseMatrix(i=1L,j=1L,x=1L, dimnames=list("1",NULL)) 
  function(x, data, 
           rmInt, ## remove Intercept
           drop, sparse_precision, 
           levels_type, # note that "data_order", "seq_len" and ".ULI" are all data-ordered, and others ('"mf"') are not
           # "data_order" is faster than ".ULI" and seems to work. attr(fitobject,"info.uniqueGeo") is unaffected.
           # extra info avout levels, for preprocessing vs for post-fit (we cannot simply pass corr_info bc we would need a ranef index 'lit' too):
           corrMat_info, # is corr_info$corrMatrices[[lit]] 
           adjMatrix=NULL, # is corr_info$adjMatrices[[lit]] 
           old_ZA=NULL, # post_fit
           lcrandfamfam) {
    ## le bidule suivant evalue le bout de formule x[[3]] et en fait un facteur. 
    ## but fac may be any vector returned by the evaluation of x[[3]] in the envir 
    rhs <- x[[3]]
    txt <- .DEPARSE(rhs) ## should be the rhs of (|) cleanly converted to a string by terms(formula,data) in .get_terms_info()
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
    assuming_spprec <- FALSE
    raneftype <- attr(x,"type")
    #if (identical(raneftype, "(.|.)")) stop("this does not occur") # does not occurs here, as explained in calling fn, .calc_Zlist()
    if ( ! is.null(raneftype)) { ## Any term with a 'spatial' keyword (incl. corrMatrix); cf comment in last case
      ## if sparse not yet determined for AR1, we generate the required info for sparse (and non-sparse) and thus assume spprec: 
      if (is.null(assuming_spprec <- sparse_precision)) assuming_spprec <- (raneftype=="AR1")  
      ## for AR1_sparse and corrMatrix, we cannot use dummy levels as created by .ULI() of factor(). The level names have special meaning
      #   matching a time concept, or user-provided names for the corrMatrix.
      ## Further, we can drop rows/cols of a correlation matrix, but not of a precision matrix
      if (raneftype %in% c("Matern","Cauchy")) { ## even in sparse case
        # uses default levels_type
      } else if (raneftype =="IMRF") {
        # for IMRF Z matches geo to uniqueGeo and A matches uniqueGeo to nodes
        levels_type <- .spaMM.data$options$uGeo_levels_type # $uGeo_levels_type used to make sure 
        #                                               that same type is used in .calc_AMatrix_IMRF() -> .as_factor()
      } else if (assuming_spprec || raneftype=="corrMatrix") {
        levels_type <- "data_order" # at some point I may have changed the type (from ".ULI" to "mf"?)
        #                                       and not seen the effect on a corrMatrix example, now included in the tests
      } else { # e.g. ranefType="adjacency", *!*assuming_spprec (immediate in the tests)
        # uses this function's default levels_type
      }
      levels_blob <- .as_factor(txt=txt,mf=data,type=levels_type) # levelstype not further needed below
      #
      if (raneftype %in% c("Matern","Cauchy", "IMRF")) {
        ff <- levels_blob$factor ## so that Z cols will not be reordered.
      } else if (raneftype=="corrMatrix") {
        if ( ! is.null(corrMat_info) && inherits(corrMat_info,"precision")) { # PRE-FIT ONLY:
          # cf example with  covStruct=list(precision=as_precision(MLdistMat)):
          # The user provided a precision matrix for a ranef of type corrMatrix
          ## we have to keep all levels of the precision matrix even those absent from the data
          # this info is in corr_info$corrMatrices, not in $adjMatrices (a prec mat is not an adj mat)
          ff <- .add_levels(ff=levels_blob$factor, adj_or_prec=corrMat_info$matrix, old_ZA=NULL)
          drop <- FALSE
        } else {
          ff <- levels_blob$factor
        }
      } else if (raneftype=="adjacency") { ## pre-fit AND post-fit: we have to keep all levels even those absent from the data
        ff <- .add_levels(ff=levels_blob$factor, adj_or_prec=adjMatrix, old_ZA=old_ZA)
        drop <- FALSE
      } else if (assuming_spprec && raneftype=="AR1") { 
        AR1_sparse_Q_ranges_blob <- .calc_AR1_sparse_Q_ranges(mf=data,levels_blob) # we need all 'time steps' for AR1 by spprec
        ff <- factor(levels_blob$factor,levels=AR1_sparse_Q_ranges_blob$seq_levelrange) ## rebuild a new factor with new levels
        drop <- FALSE
        if (anyNA(ff)) {
          stop(
            paste("Something wrong in levels of the factor for AR1 effects.") # levels_blob$factor does not match AR1_sparse_Q_ranges_blob$seq_levelrange
          )
        }
      } else { # other non-NULL raneftype's (hence not (.|.) ranefs) 
        ff <- levels_blob$factor
      }
    } else if (length(grep("c\\(\\w*\\)",txt))) { ## c(...,...) was used (actually detects ...c(...)....) (but in which context ?)
      #levels_type <- ".ULI"
      aslocator <-  parse(text=gsub("c\\(", ".ULI(", txt)) ## slow pcq ULI() est slow
      ff <- as.factor(eval(expr=aslocator,envir=data))
    } else { ## standard ( | ) rhs: automatically converts grouping variables to factors as in lme4::mkBlist (10/2014)
      # This creates a Zt matrix with rows (then ZA cols) reordered as the automatic levels of the factor
      # Hence Z is not 'ordered' (eventually diagonal) if levels are not 'ordered' in the data.
      #levels_type <- "as.factor"
      mfloc <- data
      ## 
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
    ## If info_mat was corr then it must have the levels that a precision matrix would need
    ## If info_mat_is_prec we drop nothing
    ## if assuming_spprec (i.e. if spprec already determined, or AR1) we drop nothing.
    if (drop)  ff <- droplevels(ff)
    ## Done with ff. Now the incidence matrix: 
    if (nrow(data)==1L && levels(ff)=="1") {
      im <- trivial_incidMat ## massive time gain when optimizing spatial point predictions
    } else im <- sparseMatrix(i=as.integer(ff),j=seq(length(ff)),x=1L, # ~ as(ff, "sparseMatrix") except that empty levels are not dropped
                       dimnames=list(levels(ff),NULL)) # names important for corrMatrix case at least
    # : this is faster than   im <- Matrix::fac2sparse(ff,drop.unused.levels = (drop && ! (AR1_sparse_Q || info_mat_is_prec)))
    if (!isTRUE(methods::validObject(im, test = TRUE))) {
      stop("invalid conditioning factor in random effect: ", format(rhs))
    }
    ## model matrix for LHS in [...](LHS|rhs) (Intercept if ...(1|.)) 
    tempexp <- x[[2]] ## LHS
    if (grepl("mv(", .DEPARSE(tempexp), fixed=TRUE)) {
      # We construct a template Z matrix from a fake modmat. 
      # .calc_Z_model_matrix() will 'expand' the template and .correct_ZA_mv_ranCoefs() will fill the template.
      # model_ids <- eval(tempexp) # using def of mv() function... 
      model_ids <- sub("(mv)(\\([^|]+)","c\\2", .DEPARSE(tempexp))
      model_ids <- eval(parse(text=model_ids))
      submv <- sub("mv\\([^|]+",".mv", .DEPARSE(tempexp))
      leftOfBar_form <- as.formula(paste("~", submv))
      leftOfBar_terms <- terms(leftOfBar_form)
      leftOfBar_mf <- model.frame(leftOfBar_terms, data.frame(.mv=factor(model_ids)), xlev = attr(old_ZA,"LHS_levels"))
      dummymodmat <- .calc_Z_model_matrix(leftOfBar_terms, leftOfBar_mf, raneftype = NULL, lcrandfamfam = "gaussian")
      modmat <- matrix(1, nrow=nrow(data), ncol=ncol(dummymodmat)) 
      colnames(modmat) <- colnames(dummymodmat) # the model_ids 
      #levels_type <- "model_ids"
    } else {
      leftOfBar_form <- eval(substitute(~expr, list(expr = tempexp)))
      leftOfBar_terms <- terms(leftOfBar_form) ## Implicitly assumes Intercept is included
      # LHS_levels for prediction:
      # for poly() terms, there is a crucial difference between the data (with raw variables) 
      #   and the [value of model.frame(), with monomials]. Then one cannot call recursively model.frame() on a mf.  
      leftOfBar_mf <- model.frame(leftOfBar_terms, data, xlev = attr(old_ZA,"LHS_levels")) 
      # note the test of contrasts on predict with ranCoefs with factors, in test-ranCoefs.R
      modmat <- .calc_Z_model_matrix(leftOfBar_terms, leftOfBar_mf, raneftype, lcrandfamfam) ## handles non-trivial LHS in e.g. Matern(LHS|rhs)
    }
    if (rmInt) { ## remove intercept column
      if ( ! is.na(icol <- match("(Intercept)", colnames(modmat)))) {
        if (ncol(modmat) < 2L) stop("lhs of a random-effects term cannot be an intercept only")
        modmat <- modmat[, -icol, drop = FALSE]
      }
    }
    ## This is building Z not Z(A)L hence reasonably sparse even in spatial models
    Z_ <- .calc_raw_ZA(incidMat=im, modmat) ## modmat allows simple forms of heteroscedasticity of lambda.
    if ( length(leftOfBar_mf)) { ## excludes NULL, or 0-col data.frames as in Matern(1|.)
      attr(Z_,"LHS_levels") <- .getXlevels(attr(leftOfBar_mf,"terms"), leftOfBar_mf) 
      # "                      1", "2" even if the model matrix is of the form (Intercept, contrast)
    } # else attr(Z_,"LHS_levels") <- c() # leaving it NULL is risky 
    attr(Z_, "namesTerm") <- colnames(modmat) ## e.g., "(Intercept)" ".mv2" ;  length=npar # and, consistently, the number of models in an 'mv' term
    if (identical(raneftype,"AR1")) {
      if (assuming_spprec) { ## this is TRUE is sparse_precision has not yet been determined !
        ## Following is different from levels(levels_blob$factor) which are reordered as character
        #  Effect in first fit in test-AR1, when spprec goes from NULL to FALSE
        attr(Z_,"dataordered_unique_levels") <- unique(as.character(levels_blob$factor)) ## allow reformatting for ! sparse prec
        colnames(Z_) <- AR1_sparse_Q_ranges_blob$seq_levelrange ## allow reformatting for ! sparse prec
        # ! ! ! caveat when changing the name of the following elements here, to change it elsewhere ! ! !
        attr(Z_,"AR1_block_n_u_h_s") <- AR1_sparse_Q_ranges_blob$AR1_block_n_u_h_s ## required for t_chol_Q computation
        attr(Z_,"uniqueGeo") <- AR1_sparse_Q_ranges_blob$uniqueGeo 
      } else {
        splt <- strsplit(txt,c("%in%|:|\\+| "))[[1L]]
        attr(Z_,"uniqueGeo") <- .calcUniqueGeo(data=data[,splt,drop=FALSE])
      }
    } 
    attr(Z_,"prior_lam_fac") <- attr(modmat,"prior_lam_fac") 
    return(Z_)
  }
})

# incidMat is a incidence matrix with one 1 per COL (it is transposed relative to Z) (~perm but not necess square) hence its @x has nrow(im) elements
# its # of col is the number of levels of the grouping variable.
## modmat stores (<this info>| ...) => numeric for random slope model  
.calc_raw_ZA <- function(incidMat, modmat) {
  if (ncol(modmat)==1L && (length(umm <- unique(modmat[,1L]))==1L) && umm==1) { # classic (1|.) case
    if (nrow(incidMat)>1L && .is_identity(incidMat)) { ## nrow=1 may occur when optimizing an objective function (eg bboptim) and Diagonal is not worth the cost. 
      ZA <- Diagonal(n=nrow(incidMat))
      colnames(ZA) <- rownames(incidMat) 
    } else ZA <- t(incidMat)
    # incidMat has no colnames and modmat does not provide names in the alternative general code
    attr(ZA,"is_incid") <- TRUE # we use for predVar that it then has a diagonal tcrossprod
  } else { ## first conceived for ranCoefs, with ncol(modmat)>1L. But also handles e.g. Matern(not-1|.) .... Matern(female|.)
    ZA <- vector("list",ncol(modmat))
    for (col in seq_len(ncol(modmat))) {
      ZA_col <- incidMat     # don't try to fill with template <- t(incidMat) as @x would no longer have the correct order
      ZA_col@x <- modmat[,col]
      ZA_col <-  drop0(ZA_col)
      ZA[[col]] <- t(ZA_col)
    }
    ZA <- do.call(cbind, ZA) ## colnames are repeated if modmat has several cols...
    attr(ZA,"is_incid") <- FALSE ## ___FIXME___ hmmm could be set to TRUE when  modmat is an incid matrix ? 0+ case ?
  }
  return(ZA)
}

.calc_Zlist <- function (exp_ranef_terms, data, rmInt, drop, sparse_precision=spaMM.getOption("sparse_precision"),
                         levels_type="data_order", # cf comment for Matern case in .calc_ZMatrix(); there is one non-default use in post-fit code
                         # two alternative ways to provide info about levels in .calc_ZMatrix():
                         corr_info=NULL, # pre-fit
                         ZAlist_info=NULL, # post-fit: (subset by newinold of) object$ZAlist
                         lcrandfamfam # for heteroscedastic non-gaussian random effects
                         
) {
  ## drop=TRUE elimine des niveaux spurious (test: fit cbind(Mate,1-Mate)~1+(1|Female/Male) ....)
  ## avec des consequences ultimes sur tailles des objets dans dispGammaGLM
  if (!length(exp_ranef_terms)) return(structure(list(),anyRandomSlope=FALSE))
  x3 <- lapply(exp_ranef_terms, `[[`,i=3)
  names(exp_ranef_terms) <- unlist(lapply(x3, .DEPARSE)) ## names are RHS of (.|.)
  #######
  nrand <- length(exp_ranef_terms)
  Zlist <- structure(vector("list", nrand), names=seq_len(nrand))
  for (lit in seq_along(exp_ranef_terms)) {
    # if ( ! is.null(ZAlist_info)) {
    #   if ( attr(ZAlist_info,"exp_ranef_types")[lit] %in% c("Matern","Cauchy")
    #        && is.null(attr(ZAlist_info,"AMatrices")[[as.character(lit)]])) {
    #     # cases where we can avoid some level-matching checks; but _F I X M E_ can we generalize this? 
    #     levels_type <- "seq_len" # bc new_old corr mat will have names as seq_len, 
    #     # even though "data order" was used for the original Z matrix
    #   } else levels_type <- attr(ZAlist_info, "levels_types")[lit]
    # }
    Zlist[[lit]] <- .calc_Zmatrix(exp_ranef_terms[[lit]], data=data, rmInt=rmInt,
                                   drop=drop, sparse_precision=sparse_precision, levels_type=levels_type, 
                                   corrMat_info=corr_info$corrMatrices[[lit]],
                                  adjMatrix=corr_info$adjMatrices[[lit]],
                                  old_ZA=ZAlist_info[[lit]],
                                  lcrandfamfam=lcrandfamfam[[lit]])
    ## ALL Zlist[[i]] are either 
    #     diagonal matrix (ddiMatrix) with @diag = "U"
    #  or (dg)*C*matrix ie a Compressed *C*olumn Storage (CCS) matrix 
    ##  (see http://netlib.org/linalg/html_templates/node92.html#SECTION00931200000000000000)
    ##  @x must contain the nonzero elements (except diagonal elements if @diag represents them)
    ##  @i contains the row indices of nonzero elements (except diagonal elements if @diag represents them)
    ##  @p[c] must contain the index _in x_ of the first nonzero element of column c, x[p[c]] in col c and row i[p[c]])  
  }
  if ( ! is.null(ZAlist_info)) names(Zlist) <- names(ZAlist_info) ## bc .calc_ZAlist() matches AMatrices by names (itself useful for re.form) 
  ## Subject <- list(0) ## keep this as comment; see below
  namesTerms <- vector("list", nrand)
  #levels_types <- character(nrand)
  GrpNames <- names(exp_ranef_terms)
  for (rd in seq_len(nrand)) {
    ###################
    # Subject[[i]] <- as.factor(fl[[i]]$f) # levels of grouping var for all obs ('ff' returned by locfn)
    ## : Subject was used only for random slope model, where ncol(Design) != nlevels(Subject). I tried to get rid of this.
    ## see commented use of Subject in preprocess()
    ###################
    namesTerms[[rd]] <- attr(Zlist[[rd]],"namesTerm") ## possibly several variables, eg intercept or slope... 
    names(namesTerms)[rd] <- GrpNames[rd] ## the name of the list member namesTerms[i]
    #levels_types[rd] <- attr(Zlist[[rd]],"levels_type")
  }
  return(structure(Zlist,  
                   exp_ranef_terms=exp_ranef_terms, ## matches ZAlist elements
                   exp_ranef_types=attr(exp_ranef_terms,"type"), ## matches ZAlist elements
                   #levels_types=levels_types, ## matches ZAlist elements
                   namesTerms=namesTerms, ## contains info for identifying random-coef terms
                   Xi_cols= unlist(lapply(namesTerms,length)))
  )
}
