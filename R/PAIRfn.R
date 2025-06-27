# seek 'has_fn_in_RHS' for usage...

mmfn <- function(..., VAR="PAIR", only.vars=TRUE) {
  mmmc <- mc <- match.call()
  mmmc$only.vars <- NULL
  mmmc$VAR <- NULL
  nmm  <- ...length()
  members <- vector("list", nmm)
  varnames <- character(nmm) 
  for (it in seq_len(nmm)) {
    varnames[[it]] <- member <- deparse(mmmc[[it+1L]])
    members[[it]] <- as.character(eval(parse(text=member), envir=parent.frame()))# as.character bc bug if mixing factor and integer here 
  }
  levs <- unique(.unlist(members)) 
  for (it in seq_len(nmm)) members[[it]] <- factor(members[[it]],levels=levs) 
  names(members) <- varnames
  if (only.vars) { # trouble is that I need some info for length of value to be returned.
    # if (missing(ID1)) return(NULL) # but then VAR should be in the 'data'
    # corrMatrix(PAIRfn(ID) | PAIRfn(ID,mother)) # currently fails in substitute( ~0 + expr, list(expr = x[[2]]))
    # as the expression PAIRfn(ID) does not match a variable in the data
    LHS_mf <- do.call("cbind", members) # a nmm-col matrix, possibly including NA's
    return(LHS_mf)
  } else { # seek  rhs_call[["only.vars"]] to see how this is used.
    locdf <- as.data.frame(members)
    for (member in varnames) {
      modmat <- (model.matrix(as.formula(paste("~0+",member)),locdf)) # for the levels of ID1
      modmat <- as(modmat,"CsparseMatrix")
      colnames(modmat) <- sub(paste0("^",member),"",colnames(modmat)) # must be the same names for all member's
      # so that all blocks of columns of Z will be ordered identically. 
      members[[member]] <- modmat
    }
    names(members) <- paste0("modmat", seq_len(nmm))
    attrs <- list(call=mc, spec_levs=levs, is_incid=FALSE, VAR=VAR)
    
    calc_Z_rhs <- function() {
      c(members, list(attrs=attrs ))
    }
    
    expand_Z <- function(LHS_modmat, LHS_factors, factor.) {
      .expand_Z_mm(LHS_modmat=LHS_modmat, LHS_factors=LHS_factors, 
                   modmats=members, # levels for each partner
                   factor.=factor.)
    }
    resu <- list(attrs=attrs, calc_Z_rhs=calc_Z_rhs, expand_Z=expand_Z, VAR=VAR)
    return(resu)
  }
}

.PAIRfn_old <- function(ID1, ID2, VAR="PAIR", only.vars=TRUE) {
  mc <- match.call()
  var1 <- deparse(mc$ID1)
  var2 <- deparse(mc$ID2)
  ID1 <- eval(parse(text=var1), envir=parent.frame())
  ID2 <- eval(parse(text=var2), envir=parent.frame()) 
  # To enforce same levels for the two factors
  levs <- unique(c(as.character(ID1),as.character(ID2))) # as.character bc bug if mixing factor and integer here 
  ID1 <- factor(ID1,levels=levs) 
  ID2 <- factor(ID2,levels=levs)
  if (only.vars) { # trouble is that I need some info for length of value to be returned.
    # if (missing(ID1)) return(NULL) # but then VAR should be in the 'data'
    # corrMatrix(PAIRfn(ID) | PAIRfn(ID,mother)) # currently fails in substitute( ~0 + expr, list(expr = x[[2]]))
    # as the expression PAIRfn(ID) does not match a variable in the data
    LHS_mf <- cbind(ID1, ID2) # a 2-col matrix, possibly including NA's
    return(LHS_mf)
  } else { # seek  rhs_call[["only.vars"]] to see how this is used.
    locdf <- data.frame(ID1=ID1, ID2=ID2)
    modmat1 <- (model.matrix(~0+ID1,locdf)) # for the levels of ID1
    modmat1 <- as(modmat1,"CsparseMatrix")
    colnames(modmat1) <- sub("^ID1","",colnames(modmat1))
    modmat2 <- (model.matrix(~0+ID2,locdf)) # for the levels of ID2
    modmat2 <- as(modmat2,"CsparseMatrix")
    colnames(modmat2) <- sub("^ID2","",colnames(modmat2)) # must be the same names as for modmat1
    # so that all blocks of columns of Z will be ordered identically. 
    attrs <- list(call=mc, spec_levs=levs, is_incid=FALSE)
    
    calc_Z_rhs <- function() {
      list(modmat1 = modmat1, modmat2 = modmat2, attrs=attrs )
    }
    
    expand_Z <- function(LHS_modmat, LHS_factors, factor.) {
      .expand_Z_mm_old(LHS_modmat=LHS_modmat, LHS_factors=LHS_factors, 
                       modmat1=modmat1, modmat2=modmat2, # levels for each partner
                       factor.=factor.)
    }
    resu <- list(attrs=attrs, calc_Z_rhs=calc_Z_rhs, expand_Z=expand_Z, VAR=VAR)
    return(resu)
  }
}

PAIRfn <- mmfn 

.calc_Z_LHS_modmat_mm <- function(mmvars, locdata, leftOfBar_terms, leftOfBar_mf, nr) {

  LHS_modmat <- model.matrix(leftOfBar_terms, leftOfBar_mf) # source in https://svn.r-project.org/R/trunk/src/library/stats/src/model.c
  attrs <- attributes(LHS_modmat)
  asgn <- attrs$assign # "for each column in the [LHS_modmat] 
  # ... the term in the [interpreted] formula which gave rise to the column" 
  factor. <- attr(leftOfBar_terms,"factors")[".mm",]
  seq_nr <- seq_len(nr) 
  levelsmm <- levels(locdata$.mm)
  colnams <- colnames(LHS_modmat)
  lvlcnts <- integer(length(factor.))
  for (colit in seq_len(ncol(LHS_modmat))) {
    term_col <- asgn[colit]
    # factor.[asgn[colit]] indicates whether "the alternative level(s)" of .mm affects column 'colit'
    # Hence
    if (term_col==0L || ! factor.[term_col]) { 
      blockrows  <- seq_nr 
    } else {
      altlevel_it <- lvlcnts[term_col] <- lvlcnts[term_col]+1L
      blockrows  <- nr*(altlevel_it) + seq_nr 
      if (altlevel_it==1L) colnams[colit-1L] <- sub(paste0("^\\.mm",levelsmm[2]),
                               paste0("",levelsmm[1]),
                               colnams[colit])
      colnams[colit] <- sub(paste0("^\\.mm",levelsmm[altlevel_it+1L]),
                            paste0("",levelsmm[altlevel_it+1L]),
                            colnams[colit])
    }
    LHS_modmat[seq_nr,colit] <- LHS_modmat[blockrows,colit]
  }
  colnames(LHS_modmat) <- colnams
  LHS_modmat <- LHS_modmat[seq_nr,]
  ## It's dense but a sparse matrix format is expected in later use
  LHS_modmat <- as(as(LHS_modmat,"generalMatrix"),"CsparseMatrix")
  for (st in setdiff(names(attrs),c("dim","dimnames" ))) attr(LHS_modmat, st) <- attrs[[st]]
  attr(LHS_modmat,"mmvars") <- mmvars
  ## I used a two-level factor for something that is not really such a factor:
  ## indeed the first two columns of the LHS_modmat are identical, and 
  ## additional pairs of columns for interaction terms will also be identical; 
  ## but of course it's not the final Z_.
  #
  LHS_modmat
}

