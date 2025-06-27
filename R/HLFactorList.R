.pasteCols <- function(x, collapse=":") { # fast version of apply(mf[splt],1,paste,collapse=":") ; cf from plotrix::pasteCols
  pastestring <- paste("list(",paste("x","[",seq_len(nrow(x)),",]",sep="",collapse=","),")",sep="")
  do.call(paste,c(eval(parse(text = pastestring)),sep=collapse)) 
}
# spaMM:::.pasteCols(t(blackcap)) same as apply(blackcap,1L,paste0, collapse=":") = apply(blackcap,1L,paste, collapse=":") but for names
# spaMM:::.pasteCols(t(blackcap), collapse=" ") same as apply(blackcap,1L,paste, collapse=" ") but for names

## alternative to lmerFactorList or so. Now see mkReTrms 
.as_factor <- function(txt=.DEPARSE(term), term, mf, type, has_.in.= inForm(txt, quote(`%in%`))) {
  ## Standard (.|.) ranefs are NOT handled by this function but by .rhs2factor -> base::as.factor()
  if (type=="seq_len") { ## does not try to find redundant levels. Used by predict.HLfit() for spatial terms
    splt <- NULL
    raw_levels <- seq_len(nrow(mf))
    # __F I X M E___ (was 'TAG') look for potential problem for future devel of <correlated>(mv(.)|.)
    return(list(factor=as.factor(raw_levels),splt=splt, has_fn_in_RHS=FALSE))
  } else if (type=="data_order" || ## all <keyword>() ranefs, including those with "nested nesting" for AR1 spprec || raneftype=="corrMatrix"
             type=="time_series") {
    ### 
    has_fn_in_RHS <- length(grep("^[\\w\\.]+\\(([\\w,\\(\\)\\s]*)\\)",txt, perl=TRUE))
    # setdiff(all.names(x), all.vars(x)) [https://stackoverflow.com/questions/14276728/finding-the-names-of-all-functions-in-an-r-expression]
    # would catch '+' as a function. 
    if (has_fn_in_RHS) {
      splt <- sub("^[\\w\\.]+\\(([\\w,\\(\\)\\s]*)\\)","\\1",txt, perl=TRUE) # takes the arguments 
      #   of whatever looks like a function call: "pair(a, b)" => "a, b"
      # The inner \\(\\) allows "pair(a(c), b(d))" => "a(c), b(d)". Not sure its useful.
      # The starting "^" is then necessary so that (non-first) outer match is substituted, not 
      #   the "first" (innermost, probably) match.
      # \\s for space and \\w are understood by perl=TRUE only. Otherwise use [:space:] and ...?
      splt <- strsplit(splt,c("%in%|:|\\+|, "))[[1L]] ## things to be removed so that only variable names remain
    } else splt <- strsplit(txt,c("%in%|:|\\+| "))[[1L]] ## things to be removed so that only variable names remain
    # splt <- strsplit(txt,c("%in%|:|\\+|-| "))[[1L]]  would allow '-' in an RHS but this is confusing (think about '-' in model formulas)
    splt <- splt[splt!=""]
    if ( ! all(splt %in% names(mf)) ) stop(" ! all(splt %in% names(mf))")
    if (has_.in.) {
      raw_levels <- .pasteCols(x=t(mf[splt]))
      # 
      umf_RHS <- t(unique(mf[splt])) 
      nest_row <- nrow(umf_RHS)
      raw_nested_levels <- .pasteCols(x=umf_RHS[-nest_row,,drop=FALSE])
      # I will need to go from active ZA cols to indices of unique Geo pos in data order... argh.
      # uGeo_ind_of_raw_nested_levels <- match(raw_nested_levels,unique(raw_nested_levels))
      umf_nestvar <- umf_RHS[nest_row,]
      nesting_levels <- unique(umf_nestvar)
      nested_LEVELS <- nested_Zcols <- full_LEVELS <- # uGeo_indices <- 
                                                      vector("list",length(nesting_levels))
      for (it in seq_along(nesting_levels)) {
        nested_Zcols[[it]] <- which(umf_nestvar==nesting_levels[it])
        nested_LEVELS[[it]] <- raw_nested_levels[nested_Zcols[[it]]] # e.g., geo coord for Matern; HAS TO match dimnames of a correlation matrix and lhs of colnames(ZA)
        #uGeo_indices[[it]] <- uGeo_ind_of_raw_nested_levels[nested_Zcols[[it]]]
        full_LEVELS[[it]] <- paste0(nested_LEVELS[[it]],":",nesting_levels[it])
      }
      names(nested_Zcols) <- nesting_levels
      ZA_order <- sort.list(.unlist(nested_Zcols)) # permutation to the order of colnames(ZA)
      return(list(factor=factor(raw_levels,levels=unique(raw_levels)), splt=splt, 
                  nested_Zcols=nested_Zcols,
                  nested_LEVELS=.unlist(nested_LEVELS)[ZA_order],
                  #uGeo_indices=.unlist(uGeo_indices)[ZA_order],
                  full_LEVELS=.unlist(full_LEVELS)[ZA_order],
                  has_fn_in_RHS=has_fn_in_RHS) 
      )
    } else { # raw_levels are the contents of x not some of its names...
      if (has_fn_in_RHS) {
        raw_levels <- sapply(splt, function(var) mf[[var]])
        raw_levels <- as.vector(raw_levels) # 'data order' read by cols
        # BUT IN PRINCIPLE the levels are re-created by the given RHS fn                                  
      } else if (length(splt)==1L) {
        raw_levels <- mf[[splt[1L]]]  ## depending on the user, mf[[splt[1L]]] may be integer or factor...
      } else {
        x <- t(mf[splt])
        # pastestring <- paste("list(",paste("x","[",seq_len(nrow(x)),",]",sep="",collapse=","),")",sep="")
        raw_levels <- .pasteCols(x=x) # do.call(paste,c(eval(parse(text = pastestring)),sep=":"))
      }
      return(list(factor=factor(raw_levels,levels=unique(raw_levels)),splt=splt,
                  has_fn_in_RHS=has_fn_in_RHS))
    } 
  } # I previously had type "mf" which produced the raw_levels as "data_order but returned 
    # list(factor=as.factor(raw_levels),splt=splt). By contrast, levels=. forces the order produced by unique().
}

.seq_levelrange <- function(RHS_info) {
  AR1_block_u_h_ranges <- RHS_info$AR1_block_u_h_ranges
  if ( is.null(by_values <- RHS_info$by_values) ) {
    levelrange <- AR1_block_u_h_ranges[[1L]]
    seq_levelrange <- levelrange[1L]:levelrange[2L]
    #  this reorders levels differently to their order of appearance in the data, consistently with the code producing Q -> Lunique 
  } else {
    uniqueGeos <- vector("list",length(AR1_block_u_h_ranges))
    for (lit in seq_along(AR1_block_u_h_ranges)) {
      seq_levelrange <- AR1_block_u_h_ranges[[lit]][1L]:AR1_block_u_h_ranges[[lit]][2L] # say 31 32 33 34 35 36 37 38 39 40 for by_values[[lit]]= '2'
      uniqueGeos[[lit]] <- data.frame(seq_levelrange, by_values[[lit]])
    }
    uniqueGeo <- do.call(rbind,uniqueGeos)  ## data.frame (v2.3.9)
    seq_levelrange <- .pasteCols(t(uniqueGeo)) # apply(uniqueGeo,1L,paste0,collapse=":")
    # "21:2" "22:2" "23:2" "24:2" "25:2" "26:2" "27:2" "28:2" "29:2" "30:2" "31:3" "32:3" "33:3" "34:3" "35:3" "36:3" "37:3" "38:3" "39:3" "40:3"
  }
  seq_levelrange # integer vec, or character vec from apply(.,paste0,collapse=":")
}

.calc_AR1_sparse_Q_ranges <- function(mf,RHS_info) {
  splt <- RHS_info$splt
  if (length(splt)==1L) {
    levelrange <- range(as.integer(levels(RHS_info$factor)))
    AR1_block_u_h_ranges <- list("NA"=levelrange)
    return(list(AR1_block_u_h_ranges=AR1_block_u_h_ranges))
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
    AR1_block_u_h_ranges <- lapply(by_levels,range) 
    names(AR1_block_u_h_ranges) <- by_values
    list(AR1_block_u_h_ranges=AR1_block_u_h_ranges, by_values=by_values) # by_values as names of ranges is not enough bc its type is not conserved.   
  }
}

.calc_Z_LHS_model_matrix <- function(leftOfBar_terms, leftOfBar_mf, raneftype, lcrandfamfam) {
  modmat <- model.matrix(leftOfBar_terms, leftOfBar_mf) ## contrasts.arg not immed useful, maybe later.
  
  # special handling of two cases: 
  #   TRUE|FALSE LHS (either logical, or factor with TRUE|FALSE levels),
  #   and non-gaussian ranef with numeric LHS. 
  # In the first case we removes cols of modmat (Note the tests to not remove cols in the case of factor with >2 levels, 
  #   for which (fac|.) as well as (0+fac|.) generates cols for each level of the factor).
  # In the second we check that LHS is 0+numeric.
  if ( ! (is.null(raneftype))) {
    # => quickly excludes (1|.) and simple ranCoefs. It would be nice to quickly exclude composite ranCoefs here too,
    # but passing the is_composite info (pre- and post-fit) seems a bit laborious. Instead we apply the tests on 'classe' 
    # which will de facto exclude composite ranCoefs. However some of the tests assume 
    # that 'classe' is a single value, as it should be for the target TRUE|FALSE cases, 
    # but which may be false for ranCoefs => additional test on length(classe).
    if (ncol(modmat)>1L) { # Don't put such condition(s) in parent test bc 
      # otherwise the 'parent alternative' for non-gaussian ranefs would be evaluated for such cases.
      classe <- attr(attr(leftOfBar_mf,"terms"),"dataClasses")
      if (length(classe)==1L) { # as explained above. Alternative would be to be able to test is_composite in this fn,
        # but .process_ranCoefs() depends on info provided later (by the ZAlist construction) 
        ##
        ## Cf test case     (fit1 <- fitme(migStatus ~ bool + corrMatrix(bool|loc), data=toy, corrMatrix=rcov))
        # ...(bool|.) => classe is "logical", there are no levels, colnames are "(Intercept)" "boolTRUE"
        ## With default contrasts:
        # ...(as.factor(bool)|.) => classe is "factor", levels are "FALSE" and "TRUE", colnames are "(Intercept)" and "as.factor(bool)TRUE" 
        # ...(as.character(bool)|.) => classe is "character", there are no levels, colnames are "(Intercept)" and "as.character(bool)TRUE" 
        ## With nondefault contrasts: the name of the second column changes, but there is still an "(Intercept)" col.
        ## So there seems to always be an "(Intercept)" column to remove   ***BUT***
        ## if <factor>-1 or 0+<factor> is used, there is not Intercept col, only indicators column for each level 
        ## Likewise if <logical>-1 is used, there is not Intercept col, only indicators column for each level 
        if (classe=="logical") { ## TRUE/FALSE (not factor) 
          TRUEcol <- grep(".+TRUE",colnames(modmat))
          modmat <- modmat[,TRUEcol,drop=FALSE]
        } else if (classe=="factor") {
          if (setequal(levels(leftOfBar_mf[,1L]),c("TRUE","FALSE"))) {
            if ("(Intercept)" %in% colnames(modmat)) { # <factor>-1 or 0+<factor> were NOT used
              modmat <- modmat[,colnames(modmat) != "(Intercept)",drop=FALSE]
            } else if (length(TRUEcol <- grep(".+TRUE",colnames(modmat)))) { # these should still be such a column (as implied by levels(leftOfBar_mf[,1L]))
              modmat <- modmat[,TRUEcol,drop=FALSE]
            } else warning("Unexpected case in .calc_Z_LHS_model_matrix(). Anything can happen...", immediate. = TRUE) # Reaylly unexpected, but better catch it...
          } # else the user used a non-logical factor, no reason to remove the intercept; OR the user messed up a F/T factor 
          #     by changing the levels of the logical factor (say levels(<factor>) <- c("non", "oui") with a strong french accent) 
          #     and it's not clear how to handle such mess. 
        } 
      }
    }
  } else { # raneftype is NULL... we have something like (0+wei|.) 
    if (ncol(leftOfBar_mf)==1L) { ## single variable; ncol=0L would be for (1|.) 
      if (lcrandfamfam != "gaussian" && 
          (attr(attr(leftOfBar_mf,"terms"),"dataClasses"))=="numeric" ) { 
        if (ncol(modmat)>1L) # we must guard against <Gamma ranef>(wei|), where modmat also has an intercept col
          stop(paste0("Unhandled expression in (<numeric>|.): use explicit '0 + .' syntax to remove Intercept."))
        prior_lam_fac <- modmat[,1]
        prior_lam_fac <- prior_lam_fac*prior_lam_fac
        modmat[] <- 1
        attr(modmat,"prior_lam_fac") <- prior_lam_fac
      }
    }
  } 
  return(modmat)
}

# .add_levels <- function(ff, RHS_levels) {
#   if (is.null(RHS_levels)) RHS_levels <- colnames(old_ZA) 
#   # !! for composite ranefs the RHS_info is important bc old_ZA will have repeated names
#   # may (both) well be NULL (=> later error only if automatic matching of ZA and 'L' is not possible)
#   if ( ! is.null(RHS_levels)) {  # factor() does not satisfactorily handle levels=NULL.
#     ff <- factor(ff, levels=RHS_levels) # this add the RHS_levels to the ff
#     # the opposite is not possible: 'levels' specifies the possible values. 
#     # new levels in 'ff' generate NA's.
#   }
#   ff
# }

## Simplified version of reformulas's (unexported) `%i%`.
# Cannot be named `%in%` which is a base fn... 
# ':' in formula term implies that all combinations of all levels of each 
# involved factor are generated; %i% removes unused ones.
`%i%` <- function (f1, f2, fix.order = TRUE) { # f1 and f2 should be factors
  f12 <- paste(f1, f2, sep = ":")
  u <- which(!duplicated(f12))
  levs_rank <- length(levels(f2)) * as.numeric(f1[u]) + as.numeric(f2[u])
  return(factor(f12, levels = (f12[u])[order(levs_rank)]))
}

# 'rename' reformulas::replaceForm for clarity
# reformulas v0.4.0 has an unexported replaceTerm which appears unused and identical to replaceForm...
.replaceTerm <- replaceForm

.rhs2factor <- function(data, rhs) { 
  ## standard ( | ) rhs: automatically converts grouping variables to factors as in lme4::mkBlist (10/2014)
  # This creates a Zt matrix with rows (then ZA cols) reordered as the automatic levels of the factor
  # Hence Z is not 'ordered' (eventually diagonal) if levels are not 'ordered' in the data.
  mfloc <- data
  ## 
  for (i in all.vars(rhs)) { if ( ! is.null(curf <- mfloc[[i]])) mfloc[[i]] <- as.factor(curf)}
  rhs <- .replaceTerm(rhs, quote(`:`), quote(`%i%`)) # see comment in def of %i%
  if (is.null(ff <- tryCatch(eval(substitute(factor(fac), list(fac = rhs)), mfloc),
                             error=function(e) NULL))) {
    message("couldn't evaluate grouping factor ", deparse(rhs)," within model frame:")
    if (length(grep("as.factor",rhs))) {
      stop("'as.factor' found in grouping factor term is not necessary and should be removed.",call.=FALSE)
    } else stop(" try adding grouping factor to data frame explicitly if possible",call.=FALSE)        
  }
  if (all(is.na(ff))) stop("Invalid grouping factor specification, ", deparse(rhs),call.=FALSE)
  ff
}

.trivial_incidMat <- sparseMatrix(i=1L,j=1L,x=1L, dimnames=list("1",NULL)) 

.get_levels_type <- function(corr_info, it, corr_family=corr_info$corr_families[[it]],
                              corr_type=corr_info$corr_types[[it]], default) {
  if ( identical(corr_type,"corrFamily")) {
    lty_cF <- corr_family$levels_type
    if (is.null(lty_cF)) { # Should not occur  
      # even when 'corrFamily(.|.)' is called in fitmv (cf next alternative)
      warning("No 'levels_type' in corrFamily, or earlier error?")
      ""
    } else if (lty_cF=="stub") { # can be reached by preprocessing of submodels of mv fits,
      default
    } else lty_cF    
  } else ""
}

.incidence_mat <- function(ff, drop, data, rhs) {
  ## If info_mat was corr then it must have the levels that a precision matrix would need
  ## If info_mat_is_prec we drop nothing
  ## if assuming_spprec (i.e. if spprec already determined, or AR1) we drop nothing.
  if (drop)  ff <- droplevels(ff)
  ## Done with ff. Now the incidence matrix: 
  if (nrow(data)==1L && 
      nlevels(ff)==1L && # in the adjacency case levels(ff) are all the original levels => need to avoid comparing with "1"
      levels(ff)=="1") { 
    im <- .trivial_incidMat ## precomputed => massive time gain when optimizing spatial point predictions
  } else {
    im <- sparseMatrix(i=as.integer(ff),j=seq(length(ff)),x=1L, # ~ as(ff, "sparseMatrix") except that empty levels are not dropped
                       dims=c(nlevels(ff),length(ff)), dimnames=list(levels(ff),NULL)) # names important for corrMatrix case at least
    # : this is faster than   im <- Matrix::fac2sparse(ff,drop.unused.levels = (drop && ! (AR1_sparse_Q || info_mat_is_prec)))
    # Matrix::sparse.model.matrix(< rhs formula>, data) might have been used in some cases? _F I X M E__
    if (!isTRUE(methods::validObject(im, test = TRUE))) 
      stop("invalid conditioning factor in random effect: ", format(rhs)) #_F I X M E__ find a more economical check ?
  }
  im
}

# incidMat is a incidence matrix with one 1 per COL (it is transposed relative to Z) (~perm but not necess square) hence its @x has nrow(im) elements
# its # of col is the number of levels of the grouping variable.
## modmat stores (<this LHS info>| ...) => numeric for random slope model  
.calc_raw_ZA <- function(incidMat, modmat, umm=unique(modmat[,1L])) {
  has1col <- ncol(modmat)==1L
  if (has1col && (length(umm)==1L) && umm==1) { # classic (1|.) case but still possibly in eg Matern(1|.) context
    if (nrow(incidMat)>1L && .is_identity(incidMat)) { ## nrow=1 may occur when optimizing an objective function (eg bboptim) and Diagonal is not worth the cost. 
      ZA <- Diagonal(n=nrow(incidMat))
      colnames(ZA) <- rownames(incidMat) 
    } else ZA <- t(incidMat)
    # incidMat has no colnames and modmat does not provide names in the alternative general code
    is_incid <- TRUE
    attr(is_incid,"is01col") <- FALSE # used for detection of possible sparsity so MUST indicates there are zero's
  } else { ## first conceived for ranCoefs, with ncol(modmat)>1L. But also handles e.g. Matern(0+verif|.) .... Matern(female|.)
    ZA <- vector("list",ncol(modmat))
    for (col in seq_len(ncol(modmat))) {
      ZA_col <- incidMat     # don't try to fill with template <- t(incidMat) as @x would no longer have the correct order
      ZA_col@x <- modmat[,col]
      ZA_col <-  drop0(ZA_col)
      ZA[[col]] <- t(ZA_col)
    }
    ZA <- do.call(cbind, ZA) ## colnames are repeated if modmat has several cols...
    is_incid <- has1col && (length(umm)==2L) && all(sort(umm)==c(0,1))
    # For ranCoefs, it is not incidMat;  neither for Matern(0+verif|.) as the elements are not 0 or 1.
    #  Matern(female|.) would be an exception but modmat has no is_incid attribute... B/C ratio not clear...  
    attr(is_incid,"is01col") <- is_incid # used for detection of possible sparsity so MUST indicates there are zero's
  }
  attr(ZA,"is_incid") <- is_incid # We use it in predVar for the property that it has a diagonal tcrossprod. 
  # The fact that the Elements of this diagonal are 1 might not be necessary there.
  # But there are some other uses of is_incid...
  attr(ZA, "namesTerm") <- colnames(modmat) ## e.g., "(Intercept)" ".mv2" ;  length=npar # and, consistently, the number of models in an 'mv' term
  return(ZA)
}

.expand_Z_mm_old <- function(LHS_modmat, LHS_factors, 
                         modmat1, modmat2, # levels for each partner
                         factor.) {
  asgn <- attr(LHS_modmat,"assign") 
  nblocks = ncol(LHS_modmat)
  colblocks = vector("list", nblocks)
  for (colit in seq_len(nblocks)) {
    term_col <- asgn[colit] 
    if (term_col==0L || ! LHS_factors[factor.,term_col]) { # fake Intercept || LHS col not affected by the factor,
      ## i.e. fake intercept, or col affected by another variable (say 'sex' in PAIR*sex term)
      ## term_col=0L generally represents the Intercept, and this is the way the LHS_modmat was produced.
      ## But here that col does NOT represent Intercept, but rather the first member of the pair.
      colblocks[[colit]] = .Dvec_times_Matrix(LHS_modmat[,colit], modmat1)
    } else colblocks[[colit]] = .Dvec_times_Matrix(LHS_modmat[,colit], modmat2)
  }
  do.call(cbind, colblocks)
}

.expand_Z_mm <- function(LHS_modmat, LHS_factors, 
                         modmats, # levels for each partner
                         factor.) {
  asgn <- attr(LHS_modmat,"assign") 
  nblocks <- ncol(LHS_modmat) # affected by interaction in LHS...
  colblocks <- vector("list", nblocks)
  lvlcnts <- 1L + integer(max(asgn))
  for (colit in seq_len(nblocks)) {
    term_col <- asgn[colit] 
    if (term_col==0L || ! LHS_factors[factor.,term_col]) { # fake Intercept || LHS col not affected by the factor,
      ## i.e. fake intercept, or col affected by another variable (say 'sex' in PAIR*sex term)
      ## term_col=0L generally represents the Intercept, and this is the way the LHS_modmat was produced.
      ## But here that col does NOT represent Intercept, but rather the first member of the pair.
      colblocks[[colit]] <- .Dvec_times_Matrix(LHS_modmat[,colit], modmats[[1L]])
    } else {
      lvlcnts[term_col] <- lvlcnts[term_col]+1L
      colblocks[[colit]] <- .Dvec_times_Matrix(LHS_modmat[,colit], modmats[[lvlcnts[term_col]]]) 
    }
  }
  do.call(cbind, colblocks)
}

.calc_Z_rhs_mm <- function(mmvars, data) {
  ID1 <- eval(parse(text=mmvars[1]), envir=data)
  ID2 <- eval(parse(text=mmvars[2]), envir=data) 
  # To enforce same levels for the two factors
  levs <- unique(c(as.character(ID1),as.character(ID2))) # as.character bc bug if mixing factor and integer here 
  ID1 <- factor(ID1,levels=levs) 
  ID2 <- factor(ID2,levels=levs)
  
  locdf <- data.frame(ID1=ID1, ID2=ID2)
  modmat1 <- (model.matrix(~0+ID1,locdf)) # for the levels of ID1
  modmat1 <- as(modmat1,"CsparseMatrix")
  colnames(modmat1) <- sub("^ID1","",colnames(modmat1))
  modmat2 <- (model.matrix(~0+ID2,locdf)) # for the levels of ID2
  modmat2 <- as(modmat2,"CsparseMatrix")
  colnames(modmat2) <- sub("^ID2","",colnames(modmat2)) # must be the same names as for modmat1
  # so that all blocks of columns of Z will be ordered identically. 
  list(modmat1 = modmat1, modmat2 = modmat2,
       attrs=list(spec_levs=levs, is_incid=FALSE) )
}

.mrdots <- function(..., lhs=NULL) {
  dotlist <- list(...)
  eval(str2lang(paste0("c(",paste(dotlist, collapse=","),")")))
}

.mrnamed <- function(..., lhs=NULL) {
  nmlist <- match.call()
  nmlist <- nmlist[c(1L,which(names(nmlist)!=""))]
  nmlist[[1L]] <- quote(list) 
  eval(nmlist) # evaluated list of names arguments, excluding the dots which are assumed not named
}


.calc_LHS_blob <- function(LHSexpr, oldZA, data, raneftype, lcrandfamfam,
                           rmInt) {
  LHSstr <- .DEPARSE(LHSexpr)
  if (grepl("mv(", LHSstr, fixed=TRUE)) {
    # We construct a template Z matrix from a fake modmat. 
    # .calc_Z_LHS_model_matrix() will 'expand' the template and .correct_ZA_mv_ranCoefs() will fill the template.
    #
    ## Telling apart submodel indices from names arguments.
    ## Which named args are handled in build in the formals of .mrdots() and .mrnamed()
    model_ids <- sub("(mv)(\\([^|]+)",".mrdots\\2", LHSstr) # converts to an expr such as 0+c(1,2)
    model_ids <- eval(str2lang(model_ids)) # evaluates c(1,2)...
    mrterm  <- attr(terms(as.formula(paste("~", LHSstr))),"term.labels")
    namedargs <- sub("(mv)(\\([^|]+)",".mrnamed\\2", mrterm) # converts the named args to an expression for a list
    namedargs <- eval(str2lang(namedargs)) # namedargs is list of named arguments 
    #
    lhs <- namedargs$lhs
    if (is.null(lhs)) {
      submv <- sub("mv\\([^|]+",".mv", LHSstr) # "0+.mv(1,2)", say
    } else submv <- lhs # "1", say
    leftOfBar_form <- as.formula(paste("~", submv))
    leftOfBar_terms <- terms(leftOfBar_form)
    leftOfBar_mf <- model.frame(leftOfBar_terms, data.frame(.mv=factor(model_ids)), xlev = attr(oldZA,"LHS_levels"))
    dummymodmat <- .calc_Z_LHS_model_matrix(leftOfBar_terms, leftOfBar_mf, raneftype = NULL, lcrandfamfam = "gaussian")
    LHS_modmat <- matrix(1, nrow=nrow(data), ncol=ncol(dummymodmat)) # assuming later .correct_ZA_mv_ranCoefs() call.
    colnames(LHS_modmat) <- colnames(dummymodmat) # matching the model_ids 
  } else if (grepl("mm(", LHSstr, fixed=TRUE)) {
    # We construct a template Z matrix from a fake modmat. 
    # .calc_Zmatrix() -> .expand_Z_mm() will expand according to levels 
    # of the factors for each member of the pair, using the modmat1 and modmat2 objects
    # (stored in the PAIRfn() environment when the PAIRfn() is used).
    
    mmvars <- sub("(.*)(mm)(\\([^|]+\\))(.*)","\\3", LHSstr)
    mmvars <- strsplit(sub("\\(([^|]+)\\)","\\1", mmvars),",")[[1]]
    mmvars <- sapply(mmvars, gsub, pattern=" ", replacement="")
    
    submm <- sub("mm\\([^|]+\\)",".mm", LHSstr) # replaces mm(...) by .mm 
    leftOfBar_form <- as.formula(paste("~ ", submm))
    ## If I used paste("~ 0 + ", submm) the effect would be (src/library/stats/src/model.c):
    ## /* If there is no intercept we look through the factor pattern */
    ## /* matrix and adjust the code for the *****first factor***** found so that */
    ##   /* it will be coded by dummy variables rather than contrasts. */
    ## Thus only the first two colnames of    
    #   model.matrix( as.formula(paste("~ 0 + ", submm)) , data=locdata)
    ## would be the correct ones.
    leftOfBar_terms <- terms(leftOfBar_form) 
    locdata <- vector("list", length(mmvars))
    for (mmbr in seq_along(locdata)) {
      locdata[[mmbr]] <- cbind(data,.mm=mmvars[mmbr])
    }
    locdata <- do.call(rbind, locdata)
    locdata$.mm  <- factor(locdata$.mm, levels=mmvars)
    leftOfBar_mf <- model.frame(leftOfBar_terms, locdata)
    LHS_modmat <- .calc_Z_LHS_modmat_mm(mmvars, locdata, leftOfBar_terms, 
                                        leftOfBar_mf, nr=nrow(data))
  } else {
    leftOfBar_form <- eval(substitute(~expr, list(expr = LHSexpr)))
    leftOfBar_terms <- terms(leftOfBar_form) ## Implicitly assumes Intercept is included;
    ## in particular when PAIR & PAIRfn() are used, PAIR is here only a numeric '1' constant and 
    ## a 0+PAIR generates a 0-row model frame. 
    # LHS_levels for prediction:
    # for poly() terms, there is a crucial difference between the data (with raw variables) 
    #   and the [value of model.frame(), with monomials]. Then one cannot call recursively model.frame() on a mf.  
    leftOfBar_mf <- model.frame(leftOfBar_terms, data, drop.unused.levels=TRUE, xlev = attr(oldZA,"LHS_levels")) 
    # note the test of contrasts on predict with ranCoefs with factors, in test-ranCoefs.R
    LHS_modmat <- .calc_Z_LHS_model_matrix(leftOfBar_terms, leftOfBar_mf, raneftype, 
                                           lcrandfamfam) ## handles non-trivial LHS in e.g. Matern(LHS|rhs)
  }
  #
  if (rmInt) { ## remove intercept column
    if ( ! is.na(icol <- match("(Intercept)", colnames(LHS_modmat)))) {
      if (ncol(LHS_modmat) < 2L) stop("lhs of a random-effects term cannot be an intercept only")
      LHS_modmat <- LHS_modmat[, -icol, drop = FALSE]
    }
  }
  
  LHS_blob <- list(LHS_modmat=LHS_modmat, leftOfBar_terms=leftOfBar_terms)
  if ( length(leftOfBar_mf)) { ## excludes NULL, or 0-col data.frames as in Matern(1|.)
    LHS_blob$LHS_levels <- .getXlevels(attr(leftOfBar_mf,"terms"), leftOfBar_mf)
    # : a list as defined for xlev argument of model.matrix().
  }  
  LHS_blob
}


.calc_Zmatrix <- function(x, # a term (element of exp_ranef_terms)
                          data, 
                          rmInt, ## remove Intercept
                          drop=TRUE, 
                          sparse_precision, 
                          levels_type, # note that "data_order" and "seq_len" are all data-ordered, 
                          # and other input values are not handled, but other values may be taken from cF_levels_type
                          # but now there is also cF_levels_type
                          
                          corr_info, lit,
                          oldZA=NULL, # post_fit
                          lcrandfamfam,
                          cF_levels_type=.get_levels_type(corr_info=corr_info, 
                                                          it=lit, default=levels_type)) {
  ## le bidule suivant evalue le bout de formule x[[3]] et en fait un facteur. 
  ## but fac may be any vector returned by the evaluation of x[[3]] in the envir 
  rhs <- x[[3]]
  has_.in. <- inForm(rhs, quote(`%in%`))
  ## if sparse_precision is not yet determined
  #  this build the design matrix as if it was TRUE,
  #  but adds the info dataordered_levels that allows a later modif of the design matrix if sparse_precision is set to FALSE
  assuming_spprec <- FALSE
  #                                          attr(x,"type") vs corr_info$corr_types[lit]:
  # For corrFamilies, e.g., "MaternIMRFa"     "MaternIMRFa" vs "corrFamily"
  # and for (.|.),                                     NULL vs NA
  # All code long based on attr(x,"type")... switching to corr_type would require further changes as e.g. .calc_Z_LHS_model_matrix also tests ( ! is.null(raneftype)) "## exclude (1|.) and ranCoefs! "
  # so currently we stick to attr(x,"type") and test raneftype %in% .spaMM.data$keywords$all_cF   (_F I X M E__?)
  raneftype <- attr(x,"type") 
  #if (identical(raneftype, "(.|.)")) stop("this does not occur") # does not occurs here, as explained in calling fn, .calc_Zlist()
  if ( ! is.null(raneftype)) { ## Any ranef term with a keyword (incl. corrMatrix, corrFamilies); cf comment in last case
    is_time_series_s.l. <- (raneftype=="AR1" || cF_levels_type=="time_series" )
    ## In the mv case, .calc_Zlist() is called 
    ## for each submodel while .preprocess_corrFamily has been called
    ## on a LOCAL covStruct contructed from the submodel formula.
    ## AR1(), ARp() or corrFamily(), will all have been processed 
    ## but corrFamily() will not have used the fitmv call's covStruct argument, 
    ## so appropriate levels_type info may be missing in that case,
    ## and the appropriate RHS_info may not be produced.
    ## Next, the levels_type from fitmv call's covStruct argument will be tested, 
    ## and then RHS_info will be needed. => To force its generation, I implemented the private hack 
    ## corrFamily(., levels_type="time_series")
    ## [not, doc'ed, but the doc correctly say that corrFamily() may not work in fitmv()]
    #
    old_levels_type <- attr(oldZA,"Z_levels_type")
    if (is.null(old_levels_type)) { # pre-fit; or post-fit fitted with old version of spaMM before 'old_levels_type' info has been introduced
      ## if sparse not yet determined for AR1, we generate the required info for sparse (and non-sparse) and thus assume spprec: 
      ## this block however does not correctly sets post-fit 'levels_type' to "data_order" for some composite nested ranefs,
      ## which is why old_levels_type is now used.
      if (is.null(assuming_spprec <- sparse_precision)) assuming_spprec <- is_time_series_s.l.
      ## for AR1_sparse and corrMatrix, we cannot use dummy levels as created by .ULI() of factor(). The level names have special meaning
      #   matching a time concept, or user-provided names for the corrMatrix.
      ## Further, we can drop rows/cols of a correlation matrix, but not of a precision matrix
      if (raneftype %in% c("Matern","Cauchy")) { ## even in sparse case, so this must be checked here, rather than be default final case
        # Do notchange the current .calc_Zlist()'s levels_type, typically "data_order" at this point.
      } else if (raneftype =="IMRF") {
        # for IMRF Z matches geo to uniqueGeo and A matches uniqueGeo to nodes
        levels_type <- .spaMM.data$options$uGeo_levels_type # $uGeo_levels_type used to make sure 
        #                                               that same type is used in .calc_AMatrix_IMRF() -> .as_factor()
      } else if (is_time_series_s.l. || raneftype %in% c("corrMatrix","adjacency")) {
        levels_type <- "data_order" # otherwise in prediction, any set of levels=location indices is reduced to 1 2 3... 
      } else if ( identical(corr_info$corr_types[[lit]],"corrFamily")) {
        levels_type <- cF_levels_type # notably, "time_series"
        # But not the previous test on assuming_spprec, which will typically catch the ARp case
        # in which case levels_type ("data_order") remains distinct from cF_levels_type ("time_series")   
      } else { # e.g. ranefType="adjacency", NOT assuming_spprec (immediate in the tests)
        # uses .calc_Zlist()'s default levels_type: "data_order"; or "seq_len" in post-fit calls (permuted newdata tests important here)
      }
    } else { # post-fit:
      if (raneftype %in% c("Matern","Cauchy") && ! has_.in.) {
        ## old comment below. But with  newdata, simulate.HLfit() -> 
        #  .calc_ZAlist_newdata() -> .calc_Zlist() with implicit default levels_type="data_order" 
        # so with newdata, this is "data_order" here. Test code: simulate(HLC, newdata=Loaloa[c(1:3,1:3)+0.1,])
        #
        ### OLD comment
        # keep the default post-fit levels_type, i.e. seq_len, as fast shortcut sufficient in that case
        # (permuted newdata tests important here)
      } else levels_type <- old_levels_type
    }
    
    RHS_info <- .as_factor(term=rhs,mf=data,type=levels_type,has_.in.=has_.in.) # levelstype not further needed below
    #
    if (raneftype %in% c("Matern","Cauchy", "IMRF")) {
      ff <- RHS_info$factor ## so that Z cols will not be reordered.
    } else if (raneftype=="corrMatrix") {
      if ( ! is.null(corrMat_info <- corr_info$corrMatrices[[lit]]) && inherits(corrMat_info,"precision")) { # PRE-FIT ONLY:
        # cf example with  covStruct=list(precision=as_precision(MLdistMat)):
        # The user provided a precision matrix for a ranef of type corrMatrix
        ## we have to keep all levels of the precision matrix even those absent from the data
        # this info is in corr_info$corrMatrices, not in $adjMatrices (a prec mat is not an adj mat)
        if ( ! is.null(RHS_levels <- colnames(corrMat_info$matrix))) {  
          RHS_info$factor <- factor(RHS_info$factor, levels=RHS_levels) # this add the RHS_levels to the ff
        }
        drop <- FALSE
      } # otherwise order in unchanged but cov_info_mat <- .subset_corrMatrix(...) will be used to reorder the corrMatrix to match ZA cols.
      ff <- RHS_info$factor
    } else if (raneftype=="adjacency") { ## pre-fit AND post-fit: we have to keep all levels even those absent from the data
      
      RHS_levels <- colnames(corr_info$adjMatrices[[lit]]) # in-fit
      # Post-fit we should presumably reconstruct the Z as in the fit, then the ZA,
      # then apply the ZA_update permutation?
      # The ZA colnames cannot be used anyway (they are typically permuted, and
      #   for composite ranefs, old_ZA has repeated colnames, whose order is
      #   determined by the permutation in .ZA_update(), where the cols for the ranCoefs blocks
      #   are completely scrambled).
      ##  presumptive____F I X M E____: Next line presumably fails if a non-trivial (permutation) A matrix was declared by the user... 
      ## But I have no test for this odd case.
      if (is.null(RHS_levels)) RHS_levels <- colnames(corr_info$kron_Y_Qmats[[lit]]) # composite post-fit
      # Next line should be clean for remaining case
      if (is.null(RHS_levels)) RHS_levels <- levels(attr(oldZA,"RHS_info")$factor) # of the original Z, as should be kept in ZA
      
      if ( ! is.null(RHS_levels)) {  # factor() does not satisfactorily handle levels=NULL.
        # this add the RHS_levels to the ff
        # the opposite is not possible: 'levels' specifies the possible values. 
        # new levels in 'ff' generate NA's.
        RHS_info$factor <- factor(RHS_info$factor, levels=RHS_levels) 
      }
      ff <- RHS_info$factor 
      drop <- FALSE
    } else if (assuming_spprec && is_time_series_s.l.) { 
      RHS_info <- c(RHS_info, .calc_AR1_sparse_Q_ranges(mf=data,RHS_info)) # info for 'all time steps' for t_chol_Q computation for AR1 by spprec
      seq_levelrange <- .seq_levelrange(RHS_info)
      ff <- factor(RHS_info$factor,levels=seq_levelrange) ## rebuild a new factor with new levels; seq_levelrange is  # integer vec, or character vec from apply(.,paste0,collapse=":")
      drop <- FALSE
      if (anyNA(ff)) {
        stop(
          paste("Something wrong in levels of the factor for AR1 effects.") # RHS_info$factor does not match  .seq_levelrange() result
        )
      }
      # the names of the raw Z_ should be <levels of nested var>:<levels of block index>  
    } else { # other non-NULL raneftype's (hence not (.|.) ranefs) 
      ff <- RHS_info$factor
    }
    has_fn_in_RHS <- RHS_info$has_fn_in_RHS
  } else {
    # .rhs2factor(data, rhs=block %in% year) would interpret rhs as a test (evaluated to FALSE), so:
    if (has_.in.) rhs <- .replaceTerm(rhs, quote(`%in%`), quote(`%i%`)) 
    ff <- .rhs2factor(data, rhs) # standard ( | ), including ranCoefs case
    is_time_series_s.l. <- has_fn_in_RHS <- FALSE
  }
  
  ## model matrix for LHS in [...](LHS|rhs) (Intercept if ...(1|.)) 
  LHS_blob <- .calc_LHS_blob(LHSexpr=x[[2]],
                             oldZA, data, raneftype, lcrandfamfam[[lit]], rmInt)
  LHS_modmat <- LHS_blob$LHS_modmat
  
  if (has_fn_in_RHS) {
    rhs_call <- rhs
    rhs_call[["only.vars"]] <- FALSE
    Z_blob <- eval(rhs_call, envir=data)
    if (grepl("mm(", .DEPARSE(x[[2]]), fixed=TRUE)) {
      if (TRUE) {
        Z_rhs_blob <- Z_blob$calc_Z_rhs()
        Z_ <- Z_blob$expand_Z(LHS_modmat=LHS_modmat,
                              LHS_factors=attr(LHS_blob$leftOfBar_terms,"factors"),
                              factor.=".mm")
      } else {
        Z_rhs_blob <- .calc_Z_rhs_mm(mmvars=attr(LHS_modmat,"mmvars"), data)
        Z_ <- .expand_Z_mm(LHS_modmat=LHS_modmat,
                           LHS_factors=attr(LHS_blob$leftOfBar_terms,"factors"),
                           modmats = Z_rhs_blob$modmats,
                           factor.=".mm")
      }
      for (st in names(Z_rhs_blob$attrs)) attr(Z_, st) <- Z_rhs_blob$attrs[[st]]
      attr(Z_,"mmvars") <- attr(LHS_modmat,"mmvars")
    } else { # 'PAIR' syntax
      # colnames(LHS_modmat) <- colnames(
      #   model.matrix( as.formula(substitute( ~0 + expr, list(expr = x[[2]]))) , 
      #                 data=data)      )
      asgn <- attr(LHS_modmat,"assign") 
      nblocks <- ncol(LHS_modmat)
      #
      # Rename columns. See comments in code for mm() syntax
      colnams <- colnames(LHS_modmat)
      levelsmm <- RHS_info$splt
      factor. <- attr(LHS_blob$leftOfBar_terms,"factors")[Z_blob$VAR,]
      for (colit in seq_len(nblocks)) {
        term_col <- asgn[colit]
        if (term_col==0L || ! factor.[term_col]) { 
          #
        } else {
          colnams[colit-1L] <- sub(paste0("^PAIR"),
                                   paste0("",levelsmm[1]),
                                   colnams[colit])
          colnams[colit] <- sub(paste0("^PAIR"),
                                paste0("",levelsmm[2]),
                                colnams[colit])
        }
      }
      colnames(LHS_modmat) <- colnams
      #
      data[[Z_blob$VAR]] <- factor(NA, levels=RHS_info$splt) # all $PAIR values become NA, which they should not be before.
      Z_ <- Z_blob$expand_Z(LHS_modmat=LHS_modmat, 
                            LHS_factors=attr(LHS_blob$leftOfBar_terms,"factors"),
                            factor.=Z_blob$VAR)
      for (st in names(Z_blob$attrs)) attr(Z_, st) <- Z_blob$attrs[[st]]
    }
    attr(Z_, "namesTerm") <- colnames(LHS_modmat) ## LHS terms
  } else {
    incid_RHS <- .incidence_mat(ff, drop, data, rhs) # allows simple forms of heteroscedasticity of lambda.
    Z_ <- .calc_raw_ZA(incidMat=incid_RHS, modmat=LHS_modmat) 
  }
  ## This is building Z not Z(A)L hence reasonably sparse even in spatial models
  
  attr(Z_,"LHS_levels") <- LHS_blob$LHS_levels
  if (is_time_series_s.l.) {
    ## Following is different from levels(RHS_info$factor) which are reordered as character
    #  Effect in first fit in test-AR1, when spprec goes from NULL to FALSE
    # Same as colnames but ordered as in the data, while colnames Z are ordered (blocking var ordered, and within-var ordered within block)
    RHS_info$dataordered_unique_levels <- unique(as.character(RHS_info$factor))
    attr(Z_,"RHS_info") <- RHS_info ## allow reformatting for ! sparse prec
  # } else if (identical(raneftype,"adjacency") && is.null(oldZA)) {
  # RHS_info$all_grid_levels <- levels(ff) # all levels even those not in the data
  } else if (identical(raneftype,"adjacency")) attr(Z_,"RHS_info") <- RHS_info
  
  attr(Z_,"Z_levels_type") <- levels_type # 'Z_...' as it describes cols of Z, rather than levels of a component factor 
  attr(Z_,"prior_lam_fac") <- attr(LHS_modmat,"prior_lam_fac") 
  if (has_.in. &&  ! is.null(raneftype)) {
    splits <- strsplit(unique(colnames(Z_)),":") # unique() necessary for the LHS case => Z_ has repeated colnames 
    for (it in seq_along(splits)) splits[[it]] <- tail(splits[[it]],1L)
    attr(Z_,"RHS_nesting_info") <- RHS_info[c("nested_Zcols", # indices of Z cols
                                              "nested_LEVELS", # names matching lhs of Z colnames
                                              # "uGeo_indices", 
                                              "full_LEVELS" # names matching Z colnames
                                              )] # dropping the "factor" element from the RHS_info, and a former blocksizes=table(.unlist(splits))
  }
  return(Z_)
}


.calc_Zlist <- function (exp_ranef_terms, data, rmInt, 
                         sparse_precision=.spaMM.data$options$sparse_precision,
                         levels_type="data_order", # cf comment for Matern case in .calc_ZMatrix(); there is one non-default use in post-fit code
                         # two alternative ways to provide info about levels in .calc_ZMatrix():
                         corr_info=NULL, # Only pre-fit until corrFamily was introduced
                         rd_in_mv=NULL, # for .calc_ZAlist_newdata_mv()
                         sub_oldZAlist=NULL, # post-fit: (subset by newinold of) object$ZAlist
                         lcrandfamfam # for heteroscedastic non-gaussian random effects
                         
) {
  ## drop=TRUE elimine des niveaux spurious (test: fit cbind(Mate,1-Mate)~1+(1|Female/Male) ....)
  ## avec des consequences ultimes sur tailles des objets dans dispGammaGLM
  if (!length(exp_ranef_terms)) return(structure(list(),anyRandomSlope=FALSE))
  x3 <- lapply(exp_ranef_terms, `[[`,i=3)
  names(exp_ranef_terms) <- unlist(lapply(x3, .DEPARSE)) ## names are RHS of (.|.)
  #######
  ids_in_mv <- seq_along(exp_ranef_terms)
  # in mv fits for the fit, this is called for each submodel, Zlist has submodel length.
  # Post-fit we must reconstruct a similar submodel list:
  # map_rd_mv[[2]]=c("1"=2) => rd_in_mv=c("1"=2) for submodel 2 reads as for this submodel the first Z is the second of the full-model Zlist.  
  if ( ! is.null(rd_in_mv)) {
    ids_in_mv <- rd_in_mv 
    type <- attr(exp_ranef_terms,"type")
    exp_ranef_terms <- exp_ranef_terms[ids_in_mv]
    attr(exp_ranef_terms,"type") <- type[ids_in_mv]
    # in that case the sub_oldZAlist is currently the object$ZAlist
    sub_oldZAlist <- sub_oldZAlist[ids_in_mv]  
  }
  nrand <- length(exp_ranef_terms)
  Zlist <- structure(vector("list", nrand), names=seq_len(nrand))
  for (rd in seq_along(ids_in_mv)) {
    Zlist[[rd]] <- .calc_Zmatrix(exp_ranef_terms[[rd]], data=data, rmInt=rmInt,
                                  sparse_precision=sparse_precision, 
                                  levels_type=levels_type, 
                                  corr_info=corr_info, lit=ids_in_mv[[rd]],
                                  oldZA=sub_oldZAlist[[rd]], 
                                  lcrandfamfam=lcrandfamfam)
    ## ALL Zlist[[i]] are either 
    #     diagonal matrix (ddiMatrix) with @diag = "U"
    #  or (dg)*C*matrix ie a Compressed *C*olumn Storage (CCS) matrix 
    ##  (see http://netlib.org/linalg/html_templates/node92.html#SECTION00931200000000000000)
    ##  @x must contain the nonzero elements (except diagonal elements if @diag represents them)
    ##  @i contains the row indices of nonzero elements (except diagonal elements if @diag represents them)
    ##  @p[c] must contain the index _in x_ of the first nonzero element of column c, x[p[c]] in col c and row i[p[c]])  
  }
  if ( ! is.null(sub_oldZAlist)) names(Zlist) <- names(sub_oldZAlist) ## bc .calc_ZAlist() matches AMatrices by names (itself useful for re.form) 
  ## Subject <- list(0) ## keep this as comment; see below
  namesTerms <- vector("list", nrand)
  GrpNames <- names(exp_ranef_terms)
  for (rd in seq_along(ids_in_mv)) {
    ###################
    # Subject[[i]] <- as.factor(fl[[i]]$f) # levels of grouping var for all obs ('ff' returned by locfn)
    ## : Subject was used only for random slope model, where ncol(Design) != nlevels(Subject). I tried to get rid of this.
    ## see commented use of Subject in preprocess()
    ###################
    namesTerms[[rd]] <- attr(Zlist[[rd]],"namesTerm") ## possibly several variables, eg intercept or slope... 
    names(namesTerms)[rd] <- GrpNames[rd] ## the name of the list member namesTerms[i]
  }
  return(structure(Zlist,  
                   exp_ranef_terms=exp_ranef_terms, ## matches ZAlist elements
                   exp_ranef_types=attr(exp_ranef_terms,"type"), ## matches ZAlist elements
                   namesTerms=namesTerms, ## contains info for identifying random-coef terms
                   Xi_cols= unlist(lapply(namesTerms,length)))
  )
}
