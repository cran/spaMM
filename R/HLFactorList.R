##' Create the list of model matrices from the random-effects terms in
##' the formula and the model frame.
##' 
##' @param formula model formula
##' @param mf model frame
##' @param rmInt logical scalar - should the `(Intercept)` column
##'        be removed before creating Zt
##' @param drop logical scalar indicating if elements with numeric
##'        value 0 should be dropped from the sparse model matrices 
##'

## FR->FR this is no longer such a function in lme4...but see mkReTrms next time 


spMMFactorList <- function (formula, mf, rmInt, drop) {
  ##(lme4) record dimensions and algorithm settings
  ##(lme4) create factor list for the random effects
  ## drop=TRUE elimine des niveaux spurious (test: fit cbind(Mate,1-Mate)~1+(1|Female/Male) ....)
  ## avec des consequences ultimes sur tailles des objets dans dispGammaGLM
  bars <- spMMexpandSlash(findbarsMM(formula[[length(formula)]])) ## lme4::: refs removed
  if (!length(bars)) stop("No random effects terms specified in formula")
  names(bars) <- unlist(lapply(bars, function(x) DEPARSE(x[[3]]))) 
  #######
  locfn <- function(x) {
    ## le bidule suivant evalue le bout de formule x[[3]] et en fait un facteur. Useless for spatial effects like longitude + latitude
    ## but fac may be any vector returned by the evaluation of x[[3]] in the envir mf
    rhs <- x[[3]]
    txt <- DEPARSE(rhs) ## should be the rhs of (|) cleanly converted to a string by terms(formula,data) in HLframes
    ## converts '%in%' to ':' 
    if (length(grep("%in%",txt))>0) {
      splittxt <- strsplit(txt,"%in%")[[1]]
      rhs <- as.formula(paste("~",splittxt[1],":",splittxt[2]))[[2]]
      txt <- DEPARSE(rhs)
    }
    if (length(grep("\\+",txt))>0) { ## coordinates is a vector with a single string; grep is 1 if  + was found in this single string and numeric(0) otherwise
      aslocator <- parse(text=paste("ULI(",gsub("\\+", "\\,", txt),")"))
      ## evaluating expression(ULI(...)) serves as a way to identify unique combinations of its arguments
      ff <- as.factor(eval(expr=aslocator,envir=mf))
      ## old trick using ULI to interpret %in%  
      #} else if (length(grep("%in%",txt))>0) { ## 
      #  aslocator <- parse(text=paste("ULI(",gsub("%in%", "\\,", txt),")")) ## spaMM hack using ULI 
      #  ff <- as.factor(eval(expr=aslocator,envir=mf))
    } else if (length(grep("c\\(\\w*\\)",txt))>0) { ## c(...,...) was used (actually detects ...c(...)....) 
      aslocator <-  parse(text=gsub("c\\(", "ULI(", txt)) ## slow pcq ULI() est slow
      ff <- as.factor(eval(expr=aslocator,envir=mf))
    } else { ## standard ( | ) rhs (if a single variables, it does not matter whether it is spatial or not )
      mfloc <- mf
      ## automatically converts grouping variables to factors as in lme4::mkBlist (10/2014)
      makeFac <- function(x) if (!is.factor(x)) factor(x) else x
      for (i in all.vars(rhs)) {
        if (!is.null(curf <- mfloc[[i]]))
          mfloc[[i]] <- makeFac(curf)
      }
      ## ff <- eval( expr = substitute(as.factor(fac)[, drop = drop], list(fac = rhs)), envir = mfloc)  ## older code (from lme4::MkReTrm ? )
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
    im <- as(ff, "sparseMatrix") ##FR->FR slow step; but creates S4 objects with slots as assumed in following code
    ##(lme4) Could well be that we should rather check earlier .. :
    if (!isTRUE(validObject(im, test = TRUE))) {
      stop("invalid conditioning factor in random effect: ", format(rhs))
    }
    tempexp <- x[[2]] ## analyzing the LEFT hand side for non-trivial term (random-slope model)
    mm <- model.matrix(eval(substitute(~expr, list(expr = tempexp))), mf)
    if (rmInt) {
      if (is.na(icol <- match("(Intercept)", colnames(mm)))) {
        ## break  ##FR "break may be used in wrong context: no loop is visible"
        ## but the break was in lmer code ! 
        ##FR: introduced the 'else'  
      } else {
        if (ncol(mm) < 2) 
          stop("lhs of a random-effects term cannot be an intercept only")
        mm <- mm[, -icol, drop = FALSE]
      }
    }
    ##FR at this point the code diverges from lmerFactorList  
    ans <- list(f = ff, 
                A = do.call(rBind, lapply(seq_len(ncol(mm)),function(j) im)),
                Zt = do.call(rBind, lapply(seq_len(ncol(mm)), 
                                           function(j) {
                                             im@x <- mm[, j]
                                             im
                                           })), ## obs <-> levels of ranef
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
  termsModels <- c()
  Design <- list(0)
  Subject <- list(0)
  namesTerms <- list(0)
  lambda_Xs <- list(0)
  GrpNames <- names(bars)
  for (i in 1:length(fl)) {
    termsModels[i] <- "lamScal" ## FR->FR tempo fix because it's not here that this should be determined
    Subject[[i]] <- as.factor(fl[[i]]$f) # levels of grouping var for all obs
    Design[[i]] <- as.matrix(Matrix::t(fl[[i]]$Zt)) ## nobs *(nr*npar) matrix => used to compute *ZAL* not a model matrix from a formula for lambda 
    ## colnames cannot match the names of uniqueGeo <- unique(data[,coordinates,drop=F]) because input argument do not contain info about the rownames of the data
    nt <- colnames(fl[[i]]$ST) ## length=npar
    namesTerms[[i]] <- nt ## eg intercept of slope... possibly several variables
    names(namesTerms)[i] <- GrpNames[i] ## eg intercept of slope... possibly several variables
  }
  for (iMat in seq_len(length(Design))) {
    if (is.identity(Design[[iMat]])) class(Design[[iMat]]) <- c("identityMatrix",class(Design[[iMat]]))
    attr(Design[[iMat]],"nlevels") <- ncol(Design[[iMat]])
    attr(Design[[iMat]],"colnames") <- colnames(Design[[iMat]])
  }
  ## to each (.|.) corresponds a Grouping (not necess distinct) and an element of namesTerms each identifying one or more terms
  list(Design = Design, Subject = Subject, Groupings = GrpNames,namesTerms=namesTerms,termsModels=termsModels
       #,lambda_Xs=lambda_Xs
  )
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
