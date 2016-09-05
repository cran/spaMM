## derived from lmerFactorList or so. Now see mkReTrms 

spMMFactorList <- function (formula, mf, rmInt, drop) {
  ## drop=TRUE elimine des niveaux spurious (test: fit cbind(Mate,1-Mate)~1+(1|Female/Male) ....)
  ## avec des consequences ultimes sur tailles des objets dans dispGammaGLM
  bars <- spMMexpandSlash(findbarsMM(formula[[length(formula)]])) ## this is what expands a nested random effect
  if (!length(bars)) stop("No random effects terms specified in formula")
  names(bars) <- unlist(lapply(bars, function(x) DEPARSE(x[[3]])))
  #######
  locfn <- function(x) {
    ## le bidule suivant evalue le bout de formule x[[3]] et en fait un facteur. 
    ## but fac may be any vector returned by the evaluation of x[[3]] in the envir mf
    rhs <- x[[3]]
    txt <- DEPARSE(rhs) ## should be the rhs of (|) cleanly converted to a string by terms(formula,data) in HLframes
    ## converts '%in%' to ':' 
    if (length(grep("%in%",txt))>0) {
      splittxt <- strsplit(txt,"%in%")[[1]]
      rhs <- as.formula(paste("~",splittxt[1],":",splittxt[2]))[[2]]
      txt <- DEPARSE(rhs)
    }
    #    if (length(grep("\\+",txt))>0) { ## coordinates is a vector with a single string; grep is 1 if  + was found in this single string and numeric(0) otherwise
    if (! is.null(raneftype <- attr(x,"type"))){ ## Any term with a 'spatial' keyword (incl. corrMatrix); cf comment in last case
      if (raneftype=="ar1") {
        obslevels <- eval(expr=parse(text=paste(txt)),envir=mf) ## keeps original levels, contrary to ULI
        levelrange <- range(obslevels)
        ff <- factor(obslevels,levels=levelrange[1L]:levelrange[2L])
      } else if (raneftype=="corrMatrix") {
        splt <- strsplit(txt,c("%in%|\\+| "))[[1L]]
        splt <- splt[splt!=""]
        if (length(splt)==1L && splt %in% names(mf)) { 
          ff <- as.factor(mf[[splt]]) 
          attr(ff,"generator") <- "as.factor"
        } else {
          aslocator <- parse(text=paste("ULI(",gsub("\\+", "\\,", txt),")")) ## handles for (| ...+...) 
          ff <- as.factor(eval(expr=aslocator,envir=mf))
          attr(ff,"generator") <- "ULI"
        }   
      } else {
        ## note that later ZALlist has an attr(attr(.,"ranefs"),"type") where non-spatial types are NOT null
        ## evaluating expression(ULI(...)) serves as a way to identify unique combinations of its arguments
        aslocator <- parse(text=paste("ULI(",gsub("\\+", "\\,", txt),")")) ## handles for (| ...+...) 
        ff <- as.factor(eval(expr=aslocator,envir=mf))
        ## FR->FR the general code has thesame final effect ??  
      }
    } else if (length(grep("c\\(\\w*\\)",txt))>0) { ## c(...,...) was used (actually detects ...c(...)....) 
      aslocator <-  parse(text=gsub("c\\(", "ULI(", txt)) ## slow pcq ULI() est slow
      ff <- as.factor(eval(expr=aslocator,envir=mf))
    } else { ## standard ( | ) rhs 
      ## Je pensais que "(if a single variable, it does not matter whether it is spatial or not )"
      # mais ce n'est pas vrai: s'il y a une seule variable spatiale,
      # le fait d'en faire un facteur va mener a creer une matrice Zt (puis ZA) avec des cols réordonnées comme les niveaux du factor
      # alors que cov mats / LMatrix ne sont pas réordonnées...
      # le probleme devient visible avec im <- as(ff... qui cree une matrice non diagsi on est passé par ce bloc
      mfloc <- mf
      ## automatically converts grouping variables to factors as in lme4::mkBlist (10/2014)
      makeFac <- function(x) if (!is.factor(x)) factor(x) else x
      for (i in all.vars(rhs)) {
        if (!is.null(curf <- mfloc[[i]]))
          mfloc[[i]] <- makeFac(curf)
      }
      if (is.null(ff <- tryCatch(eval(substitute(makeFac(fac),
                                                 list(fac = rhs)), mfloc),
                                 error=function(e) NULL))) {
        message("couldn't evaluate grouping factor ",
             deparse(rhs)," within model frame:")
        if (length(grep("as.factor",rhs))>0) {
          stop("'as.factor' found in grouping factor term is not necessary and should be removed.",call.=FALSE)
        } else stop(" try adding grouping factor to data ",
             "frame explicitly if possible",call.=FALSE)        
      }
      if (all(is.na(ff)))
        stop("Invalid grouping factor specification, ",
             deparse(rhs),call.=FALSE)
      if (drop) ff <- droplevels(ff)
      ## note additional code in lme4::mkBlist for handling lhs in particular
    }
    ## *******back to code common to all cases********. 
    ## This is building Z(A) not Z(A)L hence reasonably sparse even in spatial models
    if (identical(raneftype,"ar1")) {
      im <- sparseMatrix(i=as.integer(ff),j=seq(length(ff)),x=1L) # ~general code except that empty levels are not dropped
    } else im <- as(ff, "sparseMatrix") ##FR->FR slow step; but creates S4 objects with required slots
    ##(lme4) Could well be that we should rather check earlier .. :
    if (!isTRUE(validObject(im, test = TRUE))) {
      stop("invalid conditioning factor in random effect: ", format(rhs))
    }
    tempexp <- x[[2]] ## analyzing the LEFT hand side for non-trivial term (random-slope model)
    mm <- model.matrix(eval(substitute(~expr, list(expr = tempexp))), mf)
    if (rmInt) {
      if (is.na(icol <- match("(Intercept)", colnames(mm)))) {
        ## break  ##FR "break may be used in wrong context: no loop is visible"
        ## but the break was in lmer code ! ... introduced the 'else'  
      } else {
        if (ncol(mm) < 2) 
          stop("lhs of a random-effects term cannot be an intercept only")
        mm <- mm[, -icol, drop = FALSE]
      }
    }
    ##FR at this point the code diverges from lmerFactorList  
    ans <- list(f = ff, 
                A = do.call(rbind, lapply(seq_len(ncol(mm)),function(j) im)),
                ## Zt is design obs <-> levels of ranef, either dgCMatrix (sparse) or dgeMatrix (dense)
                Zt = do.call(rbind, lapply(seq_len(ncol(mm)), 
                                           function(j) {
                                             im@x <- mm[, j] ## mm stores (<this info>| ...) => numeric for random slope model 
                                             im
                                           })), 
                ST = matrix(0, ncol(mm), ncol(mm), dimnames = list(colnames(mm), 
                                                                   colnames(mm))),
                lambda_X=mm ## design matrix for predictor of lambda
    )
    if (drop) {
      ans$A@x <- rep(0, length(ans$A@x))
      ans$Zt <- drop0(ans$Zt)
    }
    ans
  }
  #######
  fl <- lapply(bars, locfn)
  Design <- list(0)
  ## Subject <- list(0) ## keep this as comment; see below
  namesTerms <- list(0)
  GrpNames <- names(bars)
  for (i in 1:length(fl)) {
    ###################
    # Subject[[i]] <- as.factor(fl[[i]]$f) # levels of grouping var for all obs ('ff' returned by locfn)
    ## : Subject was used only for random slope model, where ncol(Design) != nlevels(Subject). I tried to get rid of this.
    ## see commented use of Subject in preprocess()
    ###################
    ## is as.factor necess here ?? check that fl[[i]]$f is already factor in all cases. 
    Zt <- fl[[i]]$Zt
    ## ALL Design[[i]] are (dg)*C*matrix ie a Compressed *C*olumn Storage (CCS) matrix 
    ##  (see http://netlib.org/linalg/html_templates/node92.html#SECTION00931200000000000000)
    ##  @x must contain the nonzero elements
    ##  @i must contain the row indices of nonzero elements
    ##  @p[c] must contain the index _in x_ of the first nonzero element of column c, x[p[c]] in col c and row i[p[c]])  
    if (is.identity(Zt)) {
      Design[[i]] <- Diagonal(n=ncol(Zt)) ## diagonal matrix (ddiMatrix) with @diag = "U"
      colnames(Design[[i]]) <- rownames(Zt) ## used e.g. in prediction (no more the case ?)
      rownames(Design[[i]]) <- colnames(Zt) ## with transposition col/row names
    } else Design[[i]] <- t(Zt) ## 
    nt <- colnames(fl[[i]]$ST) ## length=npar
    namesTerms[[i]] <- nt ## possibly several variables, eg intercept or slope... 
    names(namesTerms)[i] <- GrpNames[i] ## the name of the list member namesTerms[i]
  }
  ## One should not check is.identity -> isDiagonal when 10000 points to predict... (FR->FR: modif def one of these functions ?)
  for (iMat in seq_len(length(Design))) {
    attr(Design[[iMat]],"nlevels") <- ncol(Design[[iMat]])
    attr(Design[[iMat]],"colnames") <- colnames(Design[[iMat]])
    attr(Design[[iMat]],"generator") <- attr(fl[[iMat]]$f,"generator") ## ( NULL except for corrMatrix )
  }
  list(Design = Design, #Subject = Subject, 
       namesTerms=namesTerms)
}

getgroups <- function(formlist,data) { ## reduction extreme de nlme:::getGroups.data.frame
  vlist <- lapply(formlist, function(x) { 
    val <- eval(x[[length(x)]], data)
    if (length(val) == 1) {
      return(as.factor(rep(val, nrow(data))))
    }
    else {
      return(as.factor(val)[drop = TRUE])
    }
  })
  if (length(vlist) == 1) 
    return(vlist[[1]])
  value <- do.call("data.frame", vlist)
  return(value) 
}
