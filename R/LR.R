.compare_model_structures <- function(object,object2) {
  if (inherits(object,"HLfitlist") || inherits(object2,"HLfitlist")) {
    stop("This does not yet work on HLfitlist objects")
  }
  X1 <- colnames(object$`X.pv`)
  X2 <- colnames(object2$`X.pv`)
  if (length(X1)==0L) {
    REML1 <- NULL ## compatible with both ML or REML tests
  } else REML1 <- (object$APHLs$p_v != object$APHLs$p_bv)
  if (length(X2)==0L) {
    REML2 <- NULL ## idem
  } else REML2 <- (object2$APHLs$p_v != object2$APHLs$p_bv)
  REML <- unique(c(REML1,REML2))
  meth1 <- object$HL
  meth2 <- object2$HL
  if (! identical(object$family[c("family","link")],object2$family[c("family","link")] ) ) {
    stop("Models are not nested (distinct families)")
  }
  if (! identical(meth1,meth2) || length(REML)>1 ) {
    stop("object fitted by different methods cannot be compared")
  }
  if ( ! is.null(X1)) X1 <- sapply(strsplit(X1,':'),function(x) paste(sort(x),collapse=':')) ## JBF 2015/02/23: sort variables in interaction terms before comparison
  if ( ! is.null(X2)) X2 <- sapply(strsplit(X2,':'),function(x) paste(sort(x),collapse=':'))
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
  if (length(nest)==0L) { ## NULL,NULL
    stop("Fixed-effect specifications do not appear different from each other.") 
  } else if (length(unest)==2) {
    stop("Models not nested (opposite nestings for fixed and random terms). ")
  } else {
    df1 <- length(X1)
    df2 <- length(X2)
    if (!is.null(Rnest)) {
      lambda.object <- object$lambda.object
      if (!is.null(lambda.object)) df1 <- df1+length(unlist(lambda.object$coefficients_lambdaS))
      cov.mats <- object$cov.mats
      if ( ! is.null(cov.mats)) {
        nrows <- unlist(lapply(cov.mats,NROW))
        df1 <- df1+sum(nrows*(nrows-1)/2)
      }
      lambda.object <- object2$lambda.object
      if (!is.null(lambda.object)) df2 <- df2+length(unlist(lambda.object$coefficients_lambdaS))
      cov.mats <- object2$cov.mats
      if ( ! is.null(cov.mats)) {
        nrows <- unlist(lapply(cov.mats,NROW))
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
          if ( ! is.null(fullm$distinctX.Re) ) {
            df.f.Re <- ncol(fullm$distinctX.Re)
          } else df.f.Re <- ncol(fullm$`X.pv`)
          if ( ! is.null(nullm$distinctX.Re) ) {
            df.n.Re <- ncol(nullm$distinctX.Re)
          } else df.n.Re <- ncol(nullm$`X.pv`)
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
  return(list(fullfit=fullm,nullfit=nullm,test_obj=testlik,df=df))
}

LRT <- function(object,object2,boot.repl=0,nb_cores=NULL) { ## compare two HM objects
  if (nrow(object$data)!=nrow(object2$data)) {
    stop("models were not both fitted to the same size of dataset.")
  }
  info <- .compare_model_structures(object,object2)
  nullfit <- info$nullfit
  fullfit <- info$fullfit
  test_obj <- info$test_obj
  df <- info$df
  LRTori <- 2*(logLik(fullfit,which=test_obj)-logLik(nullfit,which=test_obj))
  pvalue <- 1-pchisq(LRTori,df=df) ## but not valid for testing null components of variance
  resu <- list(nullfit=nullfit,fullfit=fullfit,basicLRT = data.frame(chi2_LR=LRTori,df=df,p_value=pvalue)) ## format appropriate for more tests  
  if (boot.repl>0) {
    if (boot.repl<100) message("It is recommended to set boot.repl>=100 for Bartlett correction")
    aslistfull <- as.list(getCall(fullfit)) 
    ## problem is for corrHLfit etc this is the call of the final HLfit call with $processed and a lot of missing original arguments  
    aslistfull$processed <- NULL ## may capture bugs 
    aslistnull <- as.list(getCall(nullfit))
    aslistnull$processed <- NULL ## may capture bugs
    simbData <- nullfit$data
    if (tolower(nullfit$family$family)=="binomial") {
      cbf <- .check_binomial_formula(nullfit=nullfit, data=fullfit$data, fullfit=fullfit)
      cbindTest <- cbf$cbindTest
      if (cbindTest) {
        nnegname <- cbf$nnegname
        nposname <- cbf$nposname
        aslistfull$formula <- cbf$full_formula
        aslistnull$formula <- cbf$null_formula
      }
    } else cbindTest <- FALSE
    eval_replicate <- function(newy,only_vector=TRUE) { ## only_vector controls how to handle errors
      if (cbindTest) {
        simbData[[nposname]] <- newy
        simbData[[nnegname]] <- .get_BinomialDen(nullfit)  - newy
      } else {simbData[[as.character(nullfit$predictor[[2L]])]] <- newy} ## allows y~x syntax for binary response
      ## analyze under both models
      aslistfull$data <- simbData
      aslistnull$data <- simbData
      fullfit <- (eval(as.call(aslistfull)))
      if (inherits(fullfit,"try-error")) {
        if (only_vector) {
          return(c(NA,NA))
        } else return(fullfit)
      }
      nullfit <- try(eval(as.call(aslistnull)))
      if (inherits(nullfit,"try-error")) {
        if (only_vector) {
          return(c(NA,NA))
        } else return(nullfit)
      }
      ## return pair of likelihoods
      return(c(logLik(fullfit,which=test_obj),logLik(nullfit,which=test_obj)))
    }
    bootLs <- .eval_boot_replicates(eval_replicate=eval_replicate,boot.repl=boot.repl,nullfit=nullfit,nb_cores=nb_cores,
                                    aslistfull=aslistfull, aslistnull=aslistnull,simbData=simbData)
    colnames(bootLs) <- paste(c("full.","null."),test_obj,sep="")
    bootdL <- bootLs[,1]-bootLs[,2]
    meanbootLRT <- 2*mean(bootdL)
    resu <- c(resu,list(rawBootLRT = data.frame(chi2_LR=LRTori,df=df,p_value=(1+sum(bootdL>=LRTori/2))/(boot.repl+1)))) ## format appropriate for more tests  
    LRTcorr <- LRTori*df/meanbootLRT
    resu <- c(resu,list(BartBootLRT = data.frame(chi2_LR=LRTcorr,df=df,p_value=1-pchisq(LRTcorr,df=df)))) ## format appropriate for more tests  
    bootInfo <- list(meanbootLRT = meanbootLRT,bootreps = bootLs)
    resu <- c(resu,list(bootInfo=bootInfo)) ## keeps the sublist structure, which is not compatible with hglmjob.R...  
  }
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}

## anova treated as alias for LRT
anova.HLfit <- function(object,object2,...) {
  LRT(object,object2,...)
}

