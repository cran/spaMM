X.GCA <- function(term, contr="contr.treatment", ...) {
  
  mc <- match.call()
  term <- deparse(mc$term)
  varnames <- strsplit(term,":")[[1]]
  ID1 <- eval(parse(text=varnames[1L]), envir=parent.frame())
  ID2 <- eval(parse(text=varnames[2L]), envir=parent.frame()) # make sure these are factors, otherwise
  if (is.factor(ID1) || is.factor(ID2)) {
    # ID1>ID2 is possible if the  two factors are ordered, and have the same levels 
    Afactor <- factor(unique(c(ID1,ID2)))
    Alevels <- levels(Afactor)
    ID1 <- factor(ID1, levels=Alevels, ordered = TRUE)
    ID2 <- factor(ID2, levels=Alevels, ordered = TRUE)
    # ID1 <- factor(ID1, ordered=TRUE)
    # ID2 <- factor(ID2, ordered=TRUE)
  } # else ID1>ID2 is assumed to be feasible 
  Z2A <- data.frame(ID1=ID1, ID2=ID2)
  which_reord <- which(ID1>ID2) 
  Z2A[which_reord,c(1L,2L)] <- Z2A[which_reord,c(2L,1L)]
  if (is.null(xlev <- mc$spec_levs)) {
    Afactor <- factor(unique(unlist(Z2A)))
  } else Afactor <- factor(unique(unlist(Z2A)), levels=xlev)
  Alevels <- levels(Afactor)
  ID1 <- factor(Z2A[,1], levels=Alevels, ordered = TRUE)
  ID2 <- factor(Z2A[,2], levels=Alevels, ordered = TRUE) # make sure these are factors, otherwise
  # (addfit <- fitme(y ~x +X.GCA(id1:id2), data=dyaddf)) fails: .GetValidData_info() -> model.frame() -> tries to assign contrasts to id1 not a factor
  contrasts(ID1) <- contr
  contrasts(ID2) <- contr
  Z2A <- data.frame(ID1=ID1, ID2=ID2)

  # if (length(Alevels)==1L) {
  #   Amatrix <- .trivial_incidMat * 0   # This case not really tested yet
  # } else {
  Amatrix <- model.matrix(~ID1,Z2A) + model.matrix(~ID2,Z2A) ## use Matrix::sparse.model.matrix here ? 
  # }
  Amatrix <- Amatrix[, -1L, drop=FALSE]
  if (contr=="contr.sum") {
    mlevels <- Alevels[-length(Alevels)] # cf lmDiallel::GCA
  } else mlevels <- Alevels[-1]
  colnames(Amatrix) <- mlevels
  resu <- as.matrix(Amatrix)
  mc$contr <- contr # keep evaluated version for later reuse 
  attr(resu,"call") <- mc
  attr(resu,"spec_levs") <- Alevels
  class(resu) <- c("factorS",class(resu))
  return(resu)
}

X.antisym <- function(term, contr="contr.treatment", ...) {
  mc <- match.call()
  term <- deparse(mc$term)
  varnames <- strsplit(term,":")[[1]]
  ID1 <- eval(parse(text=varnames[1L]), envir=parent.frame())
  ID2 <- eval(parse(text=varnames[2L]), envir=parent.frame())   
  if (is.null(xlev <- mc$spec_levs)) {
    Afactor <- factor(unique(c(ID1,ID2))) # primary fit
  } else Afactor <- factor(unique(c(ID1,ID2)), levels=xlev) # post-fit call with levels info from primary fit
  #
  Alevels <- levels(Afactor) # ALL levels of the factors used to construct the design matrices:
  ID1 <- factor(ID1, levels=Alevels) # same levels in same order for both ID1 and ID2... 
  ID2 <- factor(ID2, levels=Alevels) # ...so that the columns of the two design matrices will match.
  contrasts(ID1) <- contr
  contrasts(ID2) <- contr
  Z2A <- data.frame(ID1=ID1, ID2=ID2)
  
  # if (length(Alevels)==1L) {
  #   Amatrix <- .trivial_incidMat * 0   # This case not really tested yet
  # } else {
    Amatrix <- model.matrix(~ID1,Z2A) - model.matrix(~ID2,Z2A) ## use Matrix::sparse.model.matrix here ? 
  # }
  Amatrix <- Amatrix[, -1L, drop=FALSE]
  if (contr=="contr.sum") {
    mlevels <- Alevels[-length(Alevels)] # cf lmDiallel::GCA
  } else mlevels <- Alevels[-1]
  colnames(Amatrix) <- mlevels
  
  resu <- as.matrix(Amatrix)
  mc$contr <- contr # keep evaluated version for later reuse of the call
  attr(resu,"spec_levs") <- Alevels # keep factor information
  attr(resu,"call") <- mc # processed call for easy later re-evaluation post-fit (with factor-levels information added then) 
  class(resu) <- c("factorS",class(resu)) # "factorS" class so that the "spec_levs" attribute will be used by prediction procedure.
  return(resu)
}

# For post-fit modification of the model frame created for new data:
.hack.model.frame <- function(mf, data, spec_levs) {
  # the displayed result of model.frame() looks like a matrix but its structure shows a list...
  for (it in seq_along(mf)) {
    if(inherits(mf[[it]],"factorS")) {
      mc <- attr(mf[[it]], "call")
      mc$spec_levs <- spec_levs[[names(mf[it])]] # names(mf[it]) is the term label
      mf[[it]] <- eval(mc, envir=data)      
      attr(attr(mf,"terms"),"dataClasses")[it] <- stats::.MFclass(mf[[it]])
    }
  }
  mf
}



