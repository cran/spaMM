compare.model.structures <-
function(object,object2) {
  if (inherits(object,"HLfitlist") || inherits(object2,"HLfitlist")) {
    stop("compare.model.structures does not yet work on HLfitlist objects")
  }
  X1 <- colnames(object$X)
  X2 <- colnames(object2$X)
  if (length(X1)==0) {
    REML1 <- NULL ## compatible with both ML or REML tests
  } else REML1 <- (object$APHLs$p_v != object$APHLs$p_bv)
  if (length(X2)==0) {
    REML2 <- NULL ## idem
  } else REML2 <- (object2$APHLs$p_v != object2$APHLs$p_bv)
  REML <- unique(REML1,REML2)
  meth1 <- object$HL
  meth2 <- object2$HL
  if (! identical(meth1,meth2) || length(REML)>1 ) {
    stop("object fitted by different methods cannot be compared")
  }
  dX12 <- setdiff(X1,X2)
  dX21 <- setdiff(X2,X1)
  if (length(dX12)>0 && length(dX21)>0) {
    stop("Non-nested fixed-effect models")
  } else if (length(dX12)>0) {
    Xnest <- "2in1"
  } else if (length(dX21)>0) {
    Xnest <- "1in2"
  } else Xnest <- NULL
  ranterms1 <- attr(object$ZAlist,"ranefs")
  ranterms2 <- attr(object2$ZAlist,"ranefs")
  randist1 <- lapply(object$rand.families,function(v) paste(paste(v)[1:2],collapse="")) ## makes a string from each $family and $link 
  randist2 <- lapply(object2$rand.families,function(v) paste(paste(v)[1:2],collapse="")) ## makes a string from each $family and $link 
  ranterms1 <- paste(ranterms1,randist1) ## joins each term and its distrib
  ranterms2 <- paste(ranterms2,randist2) ## joins each term and its distrib
  dR12 <- setdiff(ranterms1,ranterms2)
  dR21 <- setdiff(ranterms2,ranterms1)
  if (length(dR12)>0 && length(dR21)>0) { 
    stop("Non-nested random-effect models")
  } else if (length(dR12)>0) {
    Rnest <- "2in1"
  } else if (length(dR21)>0) {
    Rnest <- "1in2"
  } else Rnest <- NULL
  nest <- c(Xnest,Rnest)
  unest <- unique(nest)
  if (length(nest)==0) { ## NULL,NULL
    stop("models to not appear different from each other") 
  } else if (length(unest)==2) {
    stop("Models not nested (opposite nestings for fixed and random terms). ")
  } else {
    df1 <- length(X1)
    df2 <- length(X2)
    if (!is.null(Rnest)) {
      lambda.object <- object$lambda.object
      if (!is.null(lambda.object)) df1 <- df1+nrow(lambda.object$linkscale.lambda)
      cov.mats <- object$cov.mats
      if ( ! is.null(cov.mats)) {
        nrows <- unlist(lapply(cov.mats,nrow))
        df1 <- df1+sum(nrows*(nrows-1)/2)
      }
      lambda.object <- object2$lambda.object
      if (!is.null(lambda.object)) df2 <- df2+nrow(lambda.object$linkscale.lambda)
      cov.mats <- object2$cov.mats
      if ( ! is.null(cov.mats)) {
        nrows <- unlist(lapply(cov.mats,nrow))
        df2 <- df2+sum(nrows*(nrows-1)/2)
      }
    }
    if (unest=="1in2") {
      fullm <- object2
      nullm <- object
      df <- df2-df1
    } else {
      fullm <- object
      nullm <- object2
      df <- df1-df2
    }
    if (length(nest)==2) {
      message("Nested models differing both by in their fixed and in their random terms. ")
      message("Tentatively using marginal likelihood to compare them... ")
      testlik <- "p_v" 
    } else {
      if (is.null(Rnest)) { ## fixed effect test 
        if (REML) {
          ## checking the comparability of REML fits
          if ( ! is.null(fullm$X.Re) ) {
            df.f.Re <-ncol(fullm$X.Re)
          } else df.f.Re <-ncol(fullm$X)
          if ( ! is.null(nullm$X.Re) ) {
            df.n.Re <-ncol(nullm$X.Re)
          } else df.n.Re <-ncol(nullm$X)
          if ( df.f.Re !=  df.n.Re ) {
            warning("LRT comparing REML fits with different designs is highly suspect")
          }
        }
        testlik <- "p_v"
      } else { ## random effect test
        if ( ! REML) warning("ML fits used to compare different random-effects models...")
        testlik <- "p_bv" ## used in both case, identical to p_v in the non-standard case
        stop("The two models have identical fixed-effect formulas\n and cannot yet be compared properly by this function.")
        ## need to take into account correlations in random slope models for example
      }
    } 
  }
  return(list(fullm=fullm,nullm=nullm,testlik=testlik,df=df))
}
