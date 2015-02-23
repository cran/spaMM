ULI<- function(...) {  ## Unique Location Index; '...' are simply names of variables in the data
 redondGeo<-cbind(...) ## always a matrix
 uniqueGeo<-unique(redondGeo) ## seems to preserve order
 ##FR->FR and I should in principle use a comparison of character (see ?unique)
 #### a 1 * nrow(redondGeo) matrix which nth element gives the index of the unique location in uniqueGeo :
 # designRU <- apply(redondGeo,1,function(v) {which(apply(v==t(uniqueGeo),2,all))}) ## this is slow; alternative using proxy:
 bla <- proxy::dist(uniqueGeo,redondGeo,method="Manhattan")
 designRU <- apply(bla,2,function(v){which(v==0L)})
 as.vector(designRU) ## has no names  
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
      if (term1 %in% c("adjacency","Matern","corMatern","AR1","corrMatrix")) {
        res <- findbarsMM(term[[2]])
        attr(res,"type") <- term1
        return(res) 
      } else return(findbarsMM(term[[2]])) 
    }
    c(findbarsMM(term[[2]]), findbarsMM(term[[3]]))
  }




##' Remove the random-effects terms from a mixed-effects formula,
nobarsMM <-
  function (term) 
  {
    if (!("|" %in% all.names(term))) 
      return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
      return(NULL)
    if (length(term) == 2) {
      nb <- nobarsMM(term[[2]])
      if (is.null(nb)) 
        return(NULL)
      term[[2]] <- nb
      return(term)
    }
    nb2 <- nobarsMM(term[[2]])
    nb3 <- nobarsMM(term[[3]])
    if (is.null(nb2)) 
      return(nb3)
    if (is.null(nb3)) 
      return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
  }

noOffset <-
  function (term) 
  {
    if (!("offset" %in% all.names(term))) 
      return(term)
    if (length(term) == 2) { ## this is the case if this term is  offset(...)
      term1 <- as.character(term[[1]])
      if (term1=="offset") {
        return(NULL) 
      }  
      nb <- noOffset(term[[2]])
      if (is.null(nb)) 
        return(NULL)
      term[[2]] <- nb
      return(term)
    }
    nb2 <- noOffset(term[[2]])
    nb3 <- noOffset(term[[3]])
    if (is.null(nb2)) {
      if (inherits(term,"formula")) {
        return(as.formula(paste(term[[1]],nb3))) 
      } else return(nb3)
    } 
    if (is.null(nb3)) {
      if (inherits(term,"formula")) {
        return(as.formula(paste(nb2,"~ 0"))) ## return(nb2) would only return the LHS, not a formula 
      } else return(nb2)
    } 
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
  }


nobarsNooffset <-
  function (term) 
  {
    if (!("|" %in% all.names(term) || "offset" %in% all.names(term))) 
      return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
      return(NULL)
    if (length(term) == 2) {
      term1 <- as.character(term[[1]])
      if (term1=="offset") {
        return(NULL) 
      }  
      nb <- nobarsNooffset(term[[2]])
      if (is.null(nb)) 
        return(NULL)
      term[[2]] <- nb
      return(term)
    }
    nb2 <- nobarsNooffset(term[[2]])
    nb3 <- nobarsNooffset(term[[3]])
    if (is.null(nb2)) 
      return(nb3)
    if (is.null(nb3)) 
      return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
  }

## spaces should be as in parseBars because terms can be compared as strings in later code
findSpatial <- function (term) { ## derived from findbars
  if (inherits(term,"formula") && length(term) == 2L) 
    return(findSpatial(term[[2]]))  ## process RHS
  if (is.name(term) || !is.language(term)) 
    return(NULL)
  if (term[[1]] == as.name("(")) ## i.e. (blob) but not ( | )  [update formula can add parenthese aroundspatial terms...]
    return(findSpatial(term[[2]])) 
  if (!is.call(term)) 
    stop("term must be of class call")
  if (term[[1]] == as.name("|")) ## i.e. ( | ) expression
    return(NULL)
  if (length(term) == 2L) { 
    term1 <- as.character(term[[1]])
    if (term1 %in% c("adjacency","Matern","corMatern","AR1","corrMatrix")) {
      return(term) 
    } else return(NULL) 
  }
  c(findSpatial(term[[2]]), findSpatial(term[[3]]))
}

## folwoing is confusing; type attribute of findbarsMM should be sufficient to identify spatial terms
# findNonSpatial <- function(formula) {
#   aschar <- DEPARSE(formula)
#   aschar <- gsub("adjacency([^\\|]+\\|[^)]+)","(0",aschar)
#   aschar <- gsub("Matern([^\\|]+\\|[^)]+)","(0",aschar)
#   aschar <- gsub("corMatern([^\\|]+\\|[^)]+)","(0",aschar)
#   aschar <- gsub("AR1([^\\|]+\\|[^)]+)","(0",aschar)
#   aschar <- gsub("corrMatrix([^\\|]+\\|[^)]+)","(0",aschar)
#   formula <- as.formula(aschar)
#   findbarsMM(formula)
# }

## findSpatialOrNot replaced by parseBars 11/2014; old code in a R.txt file
## spaces should be as in findSpatial because terms can be compared as strings in later code
parseBars <- function (term) { ## derived from findbars, ... return strings
  if (is.name(term) || !is.language(term)) 
    return(NULL)
  if (term[[1]] == as.name("(")) ## i.e. (blob) but not ( | ) 
  {return(parseBars(term[[2]]))} 
  if (!is.call(term)) 
    stop("term must be of class call")
  if (term[[1]] == as.name("|")) { ## i.e.  .|. expression
    bT <- paste("(",c(term),")")
    attr(bT,"type") <- "(.|.)" ## FR->FR but I need to be able to detect random-slope models
    return(bT)
  } ## le c() donne .|. et non |..
  if (length(term) == 2) {
    term1 <- as.character(term[[1]])
    if (term1 %in% c("adjacency","Matern","corMatern","AR1","corrMatrix")) {
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
    } else return(NULL) 
  }
  c(findOffset(term[[2]]), findOffset(term[[3]]))
}


asStandardFormula <- function(formula) {
  aschar <- DEPARSE(formula)
  aschar <- gsub("adjacency(","(",aschar,fixed=T)
  aschar <- gsub("Matern(","(",aschar,fixed=T)
  aschar <- gsub("corMatern(","(",aschar,fixed=T)
  aschar <- gsub("AR1(","(",aschar,fixed=T)
  aschar <- gsub("corrMatrix(","(",aschar,fixed=T)
  as.formula(aschar)
}

validData <- function(formula,resid.formula=NULL,data) {
  formula <- asStandardFormula(formula) ## removes spatial tags
  frame.form <- subbarsMM(formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
  if (!is.null(resid.formula)) {
    resid.formula <- asStandardFormula(resid.formula) ## removes spatial tags
    frame.resid.form <- subbarsMM(resid.formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
    frame.form <- paste(DEPARSE(frame.form),"+",DEPARSE(frame.resid.form[[2]]))
  }
  frame.form <- as.formula(frame.form)
  environment(frame.form) <- environment(formula)
  mf <- call("model.frame",data=data) ## it adds the formula argument below....
  mf$formula <- frame.form
  mf$drop.unused.levels <- TRUE
  mf <- eval(mf) ## data.frame with many attributes
  return(mf)
}

## cf model.frame.default from package stats mais ne traite pas les effets alea !
HLframes <- function (formula, data) {
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
    mf <- eval(mf)
    Y <- model.response(mf, "any")
    if (! is.null(Y)) {
      ## if binomial Y (may be) a numeric vector and length(dim(Y)) = length(NULL) = 0 
      ## if poisson Y (may be) an integer(!) vector and length(dim(Y)) = length(NULL) = 0
      Y <- as.matrix(Y) ## also useful in binomial case because preprocess tests ncol(Y) later. FR->FR Revise...
    }
    ####### Then constructs the design X by evaluating the model frame (fe) with fe$formula <- fixef.form
    fixef.form <- nobarsNooffset(formula) ## nobars removes the (...|...) terms...
    if (inherits(fixef.form, "formula")) {
       fe$formula <- fixef.form
       fe <- eval(fe)
       mt <- attr(fe, "terms")
       if (mt[[length(mt)]]==0) { 
    	   X <- matrix(nrow=nrow(mf),ncol=0L) ## model without fixed effects, not even an Intercept 
       } else if (!is.empty.model(mt)) {
         X <- model.matrix(mt, mf, contrasts) ## contrasts is a function from the stats package
       } else {
         mess <- pastefrom("no variables identified. Check model formula.",prefix="(!) From ")
         stop(mess)
       }
    } else X <- matrix(, NROW(Y), 0)
    storage.mode(X) <- "double" ## not clear what for...
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL
    ## creates problem with non-gaussian ranef... because the format of <rand.family>$<member fn>(<matrix>) may be vector or matrix depending on <rand.family>
    ## on the other hand X.pv %*% <vector> is a matrix, so either we try to convert everything to vector (hum) or we just care where it matters...
    list(Y = Y, X = X, wts = NULL, off = NULL, 
        mf = mf, fixef = fixef)
}
