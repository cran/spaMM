

ULI <- function(...) {  ## Unique Location Index; '...' are simply names of variables in the data
  ## the behaviour of unique(<several columns>) is to compare character representations
  ## ... but not the default behaviour is a single column
  ## => we need to standardize this 
  ## If this is modified then the computation ofuniqueGeo in HLCor needs to be modified 
  redondGeo <- cbind(...) ## always a matrix
  redondFac <- apply(redondGeo,1,paste,collapse=" ") ## always characters whatever the number of columns 
  redondFac <- as.integer(as.factor(redondFac)) ## as.factor effectively distinguishes unique character strings 
  uniqueFac <- unique(redondFac) ## seems to preserve order ## unique(<integer>) has unambiguous behaviour
  sapply(redondFac, function(v) { which(v==uniqueFac) })
}
# old usage:     corrHLfit(migStatus ~ 1+ (1|ULI(latitude,longitude)),data=blackcap,objective="p_v")

### Utilities for parsing the mixed model formula
## Functions more of less distantly derived from lme4 version of findbars

##' Determine random-effects expressions from a formula

# different ! :
# findbars(~ 1 + (1|batch/cask))
# findbarsMM(~ 1 + (1|batch/cask))


findbarsMM <-  function (term) {
    if (is.name(term) || !is.language(term)) 
      return(NULL)
    if (term[[1]] == as.name("(")) 
      return(findbarsMM(term[[2]]))
    if (!is.call(term)) 
      stop("term must be of class call")
    if (term[[1]] == as.name("|")) 
      return(term)
    if (length(term) == 2L) { ## 
      term1 <- as.character(term[[1]])
      if (term1 %in% c("adjacency","Matern","AR1","corrMatrix","ar1")) {
        res <- findbarsMM(term[[2]])
        attr(res,"type") <- term1
        return(res) 
      } else return(findbarsMM(term[[2]])) 
    }
    c(findbarsMM(term[[2]]), findbarsMM(term[[3]]))
  }

## Returns a formula without the non-spatial random-effects terms.
#currently not used
noNonSpatialbarsMM <- function (term) {
  if (!("|" %in% all.names(term))) 
    return(term)
  if (is.call(term) && term[[1]] == as.name("|")) 
    return(NULL) ## removes (|) but not Matern(|)
  if (length(term) == 2L) {
    term1 <- as.character(term[[1]])
    if (term1 %in% c("adjacency","Matern","AR1","corrMatrix","ar1")) {
      return(term) 
    } else return(noNonSpatialbarsMM(term[[2]])) 
  }
  nb2 <- noNonSpatialbarsMM(term[[2]])
  nb3 <- noNonSpatialbarsMM(term[[3]])
  if (is.null(nb2)) 
    return(nb3)
  if (is.null(nb3)) 
    return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

## Remove the random-effects terms from a mixed-effects formula
nobarsMM <- function (term) { ## different from lme4::nobars
  nb <- nobarsMM_(term)
  # if (is(term, "formula") && length(term) == 3 && is.symbol(nb)) { # wong bc is.symbol(cbind(...)) is FALSE
    if (is(term, "formula") && length(term) == 3 && ! inherits(nb,"formula")) {
      #nb <- reformulate("0", response = deparse(nb)) ## wrong for nb=cbind(...)
      nb <- as.formula(paste(deparse(nb),"~ 0")) 
  }
  # HLframes -> if (inherits(fixef.form, "formula")) { ... model.frame(...)) does not handle the formula  '~ 0'. Maybe the test should be modified
  # if (is.null(nb)) {
  #   nb <- if (is(term, "formula")) 
  #     ~0
  #   # else 0 # would convert Matern(...) en Matern(0)
  # }
  nb
}
  
nobarsMM_ <- function (term) { ## compare to lme4:::nobars_
    if (!("|" %in% all.names(term))) 
      return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
      return(NULL)
    if (length(term) == 2) {
      nb <- nobarsMM_(term[[2]])
      if (is.null(nb)) 
        return(NULL)
      term[[2]] <- nb
      return(term)
    }
    nb2 <- nobarsMM_(term[[2]])
    nb3 <- nobarsMM_(term[[3]])
    if (is.null(nb2)) 
      return(nb3)
    if (is.null(nb3)) 
      return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
  }

noOffset <- function (term) { 
  nb <- noOffset_(term)
  ## cf comment in nobarsMM
  if (is(term, "formula") && length(term) == 3  && ! inherits(nb,"formula")) {
    #nb <- reformulate("0", response = deparse(nb))
    nb <- as.formula(paste(deparse(nb),"~ 0")) 
  }
  ## cf comment in nobarsMM
  # if (is.null(nb)) {
  #   nb <- if (is(term, "formula")) 
  #     ~0
  #   # else 0 ## would convert y~x+offset(...) into y~x+0 
  # }
  nb
}


noOffset_ <- function (term)   {
  if (!("offset" %in% all.names(term))) 
    return(term)
  if (length(term) == 2) { ## this is the case if this term is  offset(...)
    term1 <- as.character(term[[1]])
    if (term1=="offset") {
      return(NULL) 
    }  
    nb <- noOffset_(term[[2]])
    if (is.null(nb)) 
      return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- noOffset_(term[[2]])
  nb3 <- noOffset_(term[[3]])
  if (is.null(nb2)) {
    #if (inherits(term,"formula")) {   ## code for autonomous noOffset fn up to 07/2016
    #  return(as.formula(paste(term[[1]],nb3))) 
    #} else 
    return(nb3)
  } 
  if (is.null(nb3)) {
    #if (inherits(term,"formula")) {   ## idem
    #  return(as.formula(paste(nb2,"~ 0"))) ## return(nb2) would only return the LHS, not a formula 
    #} else 
    return(nb2)
  } 
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

nobarsNooffset <- function (term) {
  nb <- nobarsMM(term)
  noOffset(nb)
}

## spaces should be as in parseBars because terms can be compared as strings in later code
findSpatial <- function (term, which=c("adjacency","Matern","AR1","corrMatrix","ar1")) { ## derived from findbars
  if (inherits(term,"formula") && length(term) == 2L) 
    return(c(NULL,findSpatial(term[[2]],which=which)))  
  ##       c(NULL,...) ensures that the result is always a list of language objects [critical case: term = ~ Matern() , without LHS nor other RHS terms]
  if (is.name(term) || !is.language(term)) 
    return(NULL)
  if (term[[1]] == as.name("(")) ## i.e. (blob) but not ( | )  [update formula can add parenthese aroundspatial terms...]
    return(findSpatial(term[[2]],which=which)) 
  if (!is.call(term)) 
    stop("term must be of class call")
  if (term[[1]] == as.name("|")) ## i.e. ( | ) expression
    return(NULL)
  if (length(term) == 2L) { 
    term1 <- as.character(term[[1]])
    if (term1 %in% which) {
      return(term) 
    } else return(NULL) 
  }
  c(findSpatial(term[[2]],which=which), findSpatial(term[[3]],which=which))
}


## findSpatialOrNot replaced by parseBars 11/2014; old code in a R.txt file
## spaces should be as in findSpatial because terms can be compared as strings in later code
parseBars <- function (term) { ## derived from findbars, ... return strings
  if (inherits(term,"formula") && length(term) == 2L) ## added 2015/03/23
    return(parseBars(term[[2]]))  ## process RHS
  if (is.name(term) || !is.language(term)) 
    return(NULL)
  if (term[[1]] == as.name("(")) ## i.e. (blob) but not ( | ) 
  {return(parseBars(term[[2]]))} 
  if (!is.call(term)) 
    stop("term must be of class call")
  if (term[[1]] == as.name("|")) { ## i.e.  .|. expression
    term <- spMMexpandSlash(term) ## 06/2015
    bT <- paste("(",c(term),")")
    attr(bT,"type") <- rep("(.|.)",length(bT)) ## Random-slope models are not best identified here [(X2-1|block) is not randomslope].
    return(bT)
  } ## le c() donne .|. et non |..
  if (length(term) == 2) {
    term1 <- as.character(term[[1]])
    if (term1 %in% c("adjacency","Matern","AR1","corrMatrix","ar1")) {
      res <- paste(c(term))
      attr(res,"type") <- term1
      return(res) 
    } else return(NULL) 
  }
  bT2 <- parseBars(term[[2]])
  bT3 <- parseBars(term[[3]])
  res <-c(bT2, bT3)
  attr(res,"type") <- c(attr(bT2,"type"),attr(bT3,"type"))
  return(res)
}


findOffset <- function (term) { ## derived from findbars
  if (is.name(term) || !is.language(term)) 
    return(NULL)
  if (term[[1]] == as.name("(")) ## i.e. (blob) but not ( | )
    return(findOffset(term[[2]])) 
  if (!is.call(term)) 
    stop("term must be of class call")
  if (term[[1]] == as.name("|")) ## i.e. ( | ) expression
    return(NULL)
  if (length(term) == 2) {
    term1 <- as.character(term[[1]])
    if (term1=="offset") {
      return(term) 
    } else if (term1=="~") { 
      return(findOffset(term[[2]]))
    } else return(NULL) 
  }
  c(findOffset(term[[2]]), findOffset(term[[3]]))
}


asStandardFormula <- function(formula) {
  aschar <- DEPARSE(formula)
  aschar <- gsub("adjacency(","(",aschar,fixed=T)
  aschar <- gsub("Matern(","(",aschar,fixed=T)
  aschar <- gsub("AR1(","(",aschar,fixed=T)
  aschar <- gsub("ar1(","(",aschar,fixed=T)
  aschar <- gsub("corrMatrix(","(",aschar,fixed=T)
  as.formula(aschar)
}

## function that handles prior.weights too:
getValidData <- function(formula,resid.formula=NULL,data,
                         callargs=list() ## expects a list from callargs, ie match.call of the parent frame, 
                                         ## rather than an element, to avoid premature eval of refs to data variables 
                         ) {
  envform <- environment(formula)
  formula <- asStandardFormula(formula) ## removes spatial tags
  frame.form <- subbarsMM(formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
  if (!is.null(resid.formula)) {
    resid.formula <- asStandardFormula(resid.formula) ## removes spatial tags
    frame.resid.form <- subbarsMM(resid.formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
    frame.form <- paste(DEPARSE(frame.form),"+",DEPARSE(frame.resid.form[[2]]))
  }
  frame.form <- as.formula(frame.form)
  environment(frame.form) <- envform
  mf <- match.call()
  m <- match(c("data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$weights <- callargs$prior.weights
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf$formula <- frame.form
  mf <- eval(mf)
  return(mf) ## data.frame with many attributes
}

## cf model.frame.default from package stats mais ne traite pas les effets alea !
HLframes <- function (formula, data,fitobject=NULL) {
  ## m gives either the position of the matched term in the matched call 'mc', or 0
  formula <- asStandardFormula(formula) ## strips out the spatial information, retaining the variables
  if (is.character(formula[[2]])) { ## implies that the var designated by a string (phi, for example) should not be in the data frame 
    respname <- formula[[2]]
    if (is.null(data[[respname]])) {
      validname <- respname
    } else {
      datanames <- names(data)
      allphis <- datanames[which(substr(datanames,0,nchar(respname))==respname)] ## "phi..."
      allphis <- substring(allphis,nchar(respname)+1) ## "..."
      allphis <- as.numeric(allphis[which( ! is.na(as.numeric(allphis)))  ]) ## as.numeric("...")
      if (length(allphis) == 0) {
        num <- 0
      } else num <- max(allphis)+1
      validname <-paste ( respname , num , sep="") 
    }
    data[validname] <- 1 ## adds a column $phi of 1 
    formula[[2]] <- as.name(validname) ## now the formula is standard
  }
  ####### first construct a mf for all variables required by the formula (NA would be removed here if they had not been by a previous call to validData)
  mf <- call("model.frame",data=data) ## it adds the formula argument below....
  frame.form <- subbarsMM(formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
  environment(frame.form) <- environment(formula)
  mf$formula <- frame.form
  mf$drop.unused.levels <- TRUE
  fe <- mf ## copy before further modif of mf
  mf <- eval(mf) ## will contain vars required for fixef and for ranef
  Y <- model.response(mf, "any")
  Y <- as.vector(Y) ## problem: Y is array1d here if input data frame contains array1d
  ####### Then constructs the design X by evaluating the model frame (fe) with fe$formula <- fixef.form
  fixef.form <- nobarsMM_(formula) 
  fixef.form <- noOffset_(fixef.form) 
  fixef_levels <- NULL
  fixef_terms <- NULL
  if (inherits(fixef.form, "formula")) {
    if ( ! is.null(fitobject)) { ## call for prediction
      fixef_terms <- fitobject$fixef_terms
      if (is.null(fixef_terms)) { 
        X <- matrix(nrow=nrow(mf),ncol=0L) ## model without fixed effects, not even an Intercept 
      } else if (!is.empty.model(fixef_terms)) {
        Terms <- delete.response(fixef_terms)
        fixef_mf <- model.frame(Terms, data, xlev = fitobject$fixef_levels) ## xlev gives info about the original levels
        X <- model.matrix(Terms, fixef_mf, contrasts.arg=attr(fitobject$X.pv,"contrasts")) ## contrasts.arg not immed useful, maybe later.
      } else {
        mess <- pastefrom("no variables identified. Check model formula.",prefix="(!) From ")
        stop(mess)
      }
    } else {
      fe$formula <- fixef.form
      fe <- eval(fe)
      fixef_terms <- attr(fe, "terms")
      if (fixef_terms[[length(fixef_terms)]]==0) { 
        X <- matrix(nrow=nrow(mf),ncol=0L) ## model without fixed effects, not even an Intercept 
      } else if (!is.empty.model(fixef_terms)) {
        X <- model.matrix(fixef_terms, mf, contrasts) ## contrasts is a function from the stats package
      } else {
        mess <- pastefrom("no variables identified. Check model formula.",prefix="(!) From ")
        stop(mess)
      }
      fixef_levels <- .getXlevels(fixef_terms, fe) ## added 2015/12/09 useful for predict
    }
  } else {
    X <- matrix(nrow=nrow(data), ncol=0L) ## NROW(Y) =0 if formula has no LHS, yielding inappropriate X
  }
  storage.mode(X) <- "double" ## otherwise X may be logi[] rather than num[] in particular when ncol=0
  list(Y = Y, 
       X = X, 
       wts = NULL, 
       off = NULL, 
       mf = mf, 
       fixef_levels = fixef_levels, ## added 2015/12/09 useful for predict
       fixef_terms = fixef_terms ## added 2015/12/09 useful for predict
       #,  fixef = fixef ## removed 2015/12/09 no clear use
  )
}
