.ULI <- function(...) {
  redondGeo <- cbind(...) ## NOT always a matrix: if ... is a model frame, there may be $ poly(x, 2, raw = TRUE): 'poly' num  inside
  if (ncol(redondGeo)==0L) return(rep(1L,nrow(redondGeo))) ## .generateInitPhi with constant predictor led here
  if (nrow(redondGeo)==1L) return(1L) ##  trivial case where the forllowingcode fails
  ## redondGeo is NOT always a matrix: if ... is a model frame, there may be <$ poly(x, 2, raw = TRUE): 'poly' num>  inside
  ## poly may be a two-col matrix and ncol(redoncGeo) is *1*
  redondFac <- apply(redondGeo,2L,factor,labels="") # not cute use of labels... 
  redondFac <- apply(redondFac,1L,paste,collapse=":") ## paste factors
  #redondFac <- as.character(factor(redondFac))
  uniqueFac <- unique(redondFac) ## seems to preserve order ## unique(<integer>) has unambiguous behaviour
  uniqueIdx <- seq(length(uniqueFac))
  names(uniqueIdx) <- uniqueFac
  return(uniqueIdx[redondFac])
}

.ULI_failure <- function(...) { # doc for the code of .ULI()...
  redondGeo <- cbind(...) ## NOT always a matrix: if ... is a model frame, there may be $ poly(x, 2, raw = TRUE): 'poly' num  inside
  if (ncol(redondGeo)==0L) return(rep(1L,nrow(redondGeo))) 
  if (nrow(redondGeo)==1L) return(1L) 
  for (colit in seq_len(ncol(redondGeo))) redondGeo[,colit] <- factor(redondGeo[,colit],labels="") # fails in poly() case...
  redondFac <- character(nrow(redondGeo))
  for (rowit in seq_len(nrow(redondGeo))) redondFac[rowit] <- paste(redondGeo[rowit,],collapse=":") 
  uniqueFac <- unique(redondFac) 
  uniqueIdx <- seq(length(uniqueFac))
  names(uniqueIdx) <- uniqueFac
  return(uniqueIdx[redondFac])
}


### Utilities for parsing the mixed model formula

## Remove the random-effects terms from a mixed-effects formula
.stripRanefs <- function (term) { ## different from lme4::nobars
  nb <- .stripRanefs_(term)
  if (is(term, "formula") && length(term) == 3 && ! inherits(nb,"formula")) {
    nb <- as.formula(paste(deparse(nb),"~ 0")) ## apparently needed bc model.frame(...)) does not handle the formula  '~ 0' (?)
  }
  nb
}
  
.stripRanefs_ <- function (term) { ## compare to lme4:::nobars_ ; 'term is formula or any term, recursively
  if (!("|" %in% all.names(term))) 
    return(term)
  if (is.call(term) && term[[1]] == as.name("|")) 
    return(NULL)
  if (term[[1]] == as.name("IMRF") || term[[1]] == as.name("multIMRF")) {
    return(NULL)
  }
  if (length(term) == 2) {
    nb <- .stripRanefs_(term[[2]])
    if (is.null(nb)) 
      return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- .stripRanefs_(term[[2]])
  nb3 <- .stripRanefs_(term[[3]])
  if (is.null(nb2)) 
    return(nb3)
  if (is.null(nb3)) 
    return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

## lme4::subbars "Substitute the '+' function for the '|' function in a mixed-model formula, recursively (hence the argument name term). This provides a formula suitable for the current model.frame function."
## Original version shows handling of '||'
.subbarsMM <- function (term) {   
  if (is.name(term) || !is.language(term)) 
    return(term)
  if (term[[1]] == as.name("IMRF") || term[[1]] == as.name("multIMRF")) {
    return(.subbarsMM(term[[2]]))
  }
  if (length(term) == 2) {
    term[[2]] <- .subbarsMM(term[[2]])
    return(term)
  }
  # stopifnot(length(term) >= 3) ## inefficient
  if (is.call(term) && term[[1]] == as.name("|")) 
    term[[1]] <- as.name("+")
  for (j in 2:length(term)) term[[j]] <- .subbarsMM(term[[j]])
  term
}



.stripOffset <- function (term) { 
  nb <- .stripOffset_(term)
  ## cf comment in nobarsMM
  if (is(term, "formula") && length(term) == 3  && ! inherits(nb,"formula")) {
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

.explicit_0_term <- function(x) { ## not yet used
  check_x <- .stripRanefs_(x)
  check_x <- .stripOffset_(check_x) 
  if ( ! inherits(check_x,"formula")) {
    x_len <- length(x)
    rhs <- paste0("0 + ",x[x_len])
    if (x_len > 2L) {
      res <- formula(paste("lhs ~", rhs))
      res[[2L]] <- x[[2L]]
      return(res)
    } else return(formula(paste("~", rhs)))
  } else return(x)
}

.stripOffset_ <- function (term)   {
  if (!("offset" %in% all.names(term))) 
    return(term)
  if (length(term) == 2) { ## this is the case if this term is  offset(...)
    term1 <- as.character(term[[1]])
    if (term1=="offset") {
      return(NULL) 
    }  
    nb <- .stripOffset_(term[[2]])
    if (is.null(nb)) 
      return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- .stripOffset_(term[[2]])
  nb3 <- .stripOffset_(term[[3]])
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

# To match the lme4 'syntactic sugar' (.||.):
.expandDoubleVert <- function(term) {
  frml <- formula(substitute(~x, list(x = term[[2]])))
  newtrms <- paste0("0+", attr(terms(frml), "term.labels"))
  if (attr(terms(frml), "intercept") != 0) 
    newtrms <- c("1", newtrms)
  nterms <- length(newtrms)
  collapsand <- character(nterms)
  for (it in seq_len(nterms)) collapsand[it] <- paste0(newtrms[it], "|", deparse(term[[3]]))
  as.formula(paste("~(", paste(collapsand, collapse = ")+("), ")"))[[2]]
}

.expandDoubleVerts <- function (term) {
  if (!is.name(term) && is.language(term)) {
    if (term[[1]] == as.name("(")) {
      term[[2]] <- .expandDoubleVerts(term[[2]])
    }
    stopifnot(is.call(term))
    if (term[[1]] == as.name("||")) 
      return(.expandDoubleVert(term))
    term[[2]] <- .expandDoubleVerts(term[[2]])
    if (length(term) == 3) 
      term[[3]] <- .expandDoubleVerts(term[[3]])
  }
  term
}

.parseBars <- function (term) { ## derived from findbars, ... 
  if (inherits(term,"formula") && length(term) == 2L) ## added 2015/03/23
    return(.parseBars(term[[2]]))  ## process RHS
  if (is.name(term) || !is.language(term)) 
    return(NULL)
  if (term[[1]] == as.name("(")) ## i.e. (blob) but not ( | ) 
  {return(.parseBars(term[[2]]))} 
  if (!is.call(term)) 
    stop("term must be of class call")
  if (term[[1]] == as.name("|")) { ## i.e.  .|. expression
    attr(term,"type") <- rep("(.|.)",length(paste(c(term)))) ## paste is the simplest way to count terms (rather than elements of a single term)
    return(term)
  } ## le c() donne .|. et non |..
  if (term[[1]] == as.name("IMRF")) {
    pars <- paste0(paste0(names(term),"=",term)[-c(1,2)], collapse=",") 
    pars <- eval(parse(text = paste("list(",pars,")")))
    if ( ! is.null(pars$spde)) {
      warning("'spde' argument is obsolete. Use 'model' instead.")
      pars$model <- pars$spde
      pars$spde <- NULL
    }
    if ( ! is.null(pars$model)) {
      pars$no <- FALSE
      # it's hard to guess the alpha parameter from the object!
      if (pars$model$param.inla[["B0"]][1,3]!=0) {
        stop("Apparently fractional SPDE model, not implemented in spaMM.")
      } else {
        SPDE_alpha <- pars$model$param.inla[["B1"]][1,3]
        if ( ! (SPDE_alpha  %in% c(1,2))) stop("Unrecognized model from inla.spde2. Contact the spaMM maintainer to extend spaMM.")
        pars$SPDE_alpha <- SPDE_alpha
      }
    } else if (is.null(pars$no)) pars$no <- TRUE
    if (is.null(pars$ce)) pars$ce <- TRUE
    attr(term,"type") <- structure("IMRF", pars=pars) # pars as attr to type avoid problems in building the return value.
    return(term) # (lhs|rhs) (language, with the ()); or character string.
  }
  if (term[[1]] == as.name("multIMRF")) { # useful for .findSpatial() before preprocessing, as in filled.mapMM()
    attr(term,"type") <- as.character(term[[1]])
    return(term) 
  }
  if (length(term) == 2) { # no extra args ater the grouping rhs !
    term1 <- as.character(term[[1]])
    if (term1 %in% c("adjacency","Matern","Cauchy","AR1","corrMatrix")) {
      attr(term,"type") <- term1
      return(term) 
    } else return(NULL) 
  }
  # This point may not be reached by the outer .parseBars() call, so the following is not consistently terminal code.
  bT2 <- .parseBars(term[[2]])
  bT3 <- .parseBars(term[[3]])
  res <- c(bT2, bT3) ## if the bT's are language objects, c() creates a *list* of them
  attr(res,"type") <- c(attr(bT2,"type"),attr(bT3,"type"))
  # res is a list or a vector, depending on as_character
  # If a list, each element has a type (conversely, vector elements cannot have attributes except a name)
  # The list, or the vector, *has* a type attribute. It's a full vector, hence it has "(.|.)" for elements that have no type
  return(res)
}

.lhs_rhs_bars <- function(barlist) { ## replicates the old .findbarsMM result, but starting from a .parseBars() result
  ## this means that a .|. term loses its "(.|.)" attribute, 
  ##        and that a Matern(.|.) becomes .|. without the Matern(), but keeping its "Matern" attribute.
  res <- barlist
  for (it in seq_along(res)) {
    if ( ! is.null(type <- attr(res[[it]],"type"))) { 
      if (type == "(.|.)") {
        attr(res[[it]],"type") <- NULL
      } else res[[it]] <- structure(.parseBars(res[[it]][[2]]), type=type)} # 1 | grp is unchanged, and 
  }
  return(structure(res,type=attr(barlist,"type"))) # list attribute modified, 2019/3/6
}

.as_char_bars <- function(barlist) { # character vector from language list
  if (is.null(barlist)) return(NULL)
  # ELSE
  res <- character(length(barlist))
  for (it in seq_along(res)) {
    if ( ! is.null(type <- attr(barlist[[it]],"type"))) { 
      if (type == "(.|.)") {
        res[it] <- paste("(",c(barlist[[it]]),")")
      } else res[it] <- paste(c(barlist[[it]]))
    }
  }
  return(structure(res,type=attr(barlist,"type")))
}

.process_bars <- function(formula, barlist, expand= (which.==""),
                          as_character= (which.==""), # return strings vs. full <keyword>(.|.) terms
                          which.="") {
  if (missing(barlist)) barlist <- .parseBars(formula)
  # NULL barlist is a result from a previous.parseBars / .process_bars call !   
  if (expand) {
    barlist <- .spMMexpandSlash(barlist)
  }
  attr(barlist,"type") <- unlist(lapply(barlist,attr,which="type")) # before the element attributes are hacked by .lhs_rhs_bars()
  if (which.=="exp_ranef_terms") {
    return(.lhs_rhs_bars(barlist)) # keeps the "type" attribute of the list
  } else if (as_character) {
    return(.as_char_bars(barlist)) # keeps the "type" attribute of the list
  } else return(barlist)
} # with a "type" attribute to the list 

## spaces should be as in parseBars because terms can be compared as strings in later code
.findSpatial <- function (term, barlist, expand=FALSE, as_character=FALSE, nomatch=NULL) { 
  res <- .process_bars(term, 
                       barlist=barlist, ## may be missing
                       expand=expand, as_character=as_character)
  if (is.null(nomatch)) {
    res <- res[attr(res,"type")!="(.|.)"] ## nomatch=NULL removes terms of type "(.|.)" 
  } else res[attr(res,"type")=="(.|.)"] <- nomatch ## typically sets to nomatch=NA the terms of type "(.|.)" 
  return(res)
}

.findOffset <- function (term) { ## derived from findbars
  if (is.name(term) || !is.language(term)) 
    return(NULL)
  if (term[[1]] == as.name("(")) ## i.e. (blob) but not ( | )
    return(.findOffset(term[[2]])) 
  if (!is.call(term)) 
    stop("term must be of class call")
  if (term[[1]] == as.name("|")) ## i.e. ( | ) expression
    return(NULL)
  if (length(term) == 2) {
    term1 <- as.character(term[[1]])
    if (term1=="offset") {
      return(term) 
    } else if (term1=="~") { 
      return(.findOffset(term[[2]]))
    } else return(NULL) 
  }
  c(.findOffset(term[[2]]), .findOffset(term[[3]]))
}


.asNoCorrFormula <- function(formula) {
  aschar <- .DEPARSE(formula)
  aschar <- gsub("adjacency(","(",aschar,fixed=TRUE)
  aschar <- gsub("Matern(","(",aschar,fixed=TRUE)
  aschar <- gsub("Cauchy(","(",aschar,fixed=TRUE)
  aschar <- gsub("AR1(","(",aschar,fixed=TRUE)
  aschar <- gsub("corrMatrix(","(",aschar,fixed=TRUE)
  #                     (  1st group )  (2nd group)        \1: only first group is retained
  aschar <- gsub("multIMRF\\((.+?\\|[^,]+?), ([^)]+?)\\)", "\\(\\1\\)", aschar, fixed = FALSE) # removes IMRF and the parameters
  aschar <- gsub("IMRF\\((.+?\\|[^,]+?), ([^)]+?)\\)", "\\(\\1\\)", aschar, fixed = FALSE) # removes IMRF and the parameters
  ## thank you https://regex101.com/ ... ? (:lazy matching) is essential
  as.formula(aschar)
}

## function that handles prior.weights too:
.getValidData <- function(formula,resid.formula=NULL,data,
                         callargs=list() ## expects a list from callargs, ie match.call of the parent frame, 
                                         ## rather than an element, to avoid premature eval of refs to data variables 
                         ) {
  ## 
  envform <- environment(formula)  
  ## The aim is to reduce environment(formula) to almost nothing as all vars should be in the data. 
  #  The revised Infusion code shows how to use the 'data' for prior.weights in programming, using eval(parse(text=paste(<...>,priorwName,<...>)))
  # As long as prior.weights was still using the environment, fitting functions still had
  #     mc <- oricall
  #     oricall$formula <- .preprocess_formula(formula) 
  #  (stripping of mc$formula being done only later by .as_predictor() )
  #  where we otherwise have the opposite order
  #     oricall$formula <- .preprocess_formula(formula) 
  #     mc <- oricall
  formula <- .asNoCorrFormula(formula) ## removes spatial tags
  frame.form <- .subbarsMM(formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
  if (!is.null(resid.formula)) {
    resid.formula <- .asNoCorrFormula(resid.formula) ## removes spatial tags
    frame.resid.form <- .subbarsMM(resid.formula) ## this comes from lme4 and converts (...|...) terms to some "+" form
    frame.form <- paste(.DEPARSE(frame.form),"+",.DEPARSE(frame.resid.form[[2]]))
  }
  check <- grep('$',frame.form,fixed=TRUE)
  if (length(check)) {
    message("'$' detected in formula: suspect and best avoided. 
            In particular, one should never need to specify the 'data' in the 'formula'. 
            See help('good-practice') for explanations.")
  }
  frame.form <- as.formula(frame.form)
  environment(frame.form) <- envform
  mf <- match.call()
  m <- match(c("data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  verif <- try(eval(callargs$prior.weights,data),silent=TRUE)
  if (inherits(verif,"try-error")) stop("All variables should be in the 'data', including those for prior weights.")
  mf$weights <- callargs$prior.weights
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- get("model.frame", asNamespace("stats"))
  mf$formula <- frame.form
  mf <- eval(mf)
  return(mf) ## data.frame with many attributes
}

.find_validname <- function(formula, data) {
  ## implies that the var designated by a string (phi, for example) should not be in the data frame 
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
    validname <-paste0( respname , num) 
  }
  return(validname)
}

## cf model.frame.default from package stats mais ne traite pas les effets alea !
.HLframes <- function (formula, data) {
  ## m gives either the position of the matched term in the matched call 'mc', or 0
  formula <- .asNoCorrFormula(formula) ## strips out the spatial information, retaining the variables
  if (is.character(formula[[2]])) {
    validname <- .find_validname(formula, data)
    data[validname] <- 1 ## adds a column $phi of 1 
    formula[[2]] <- as.name(validname) ## now the formula is standard
  }
  ####### first construct a mf for all variables required by the formula (NA would be removed here if they had not been by a previous call to validData)
  frame.form <- .subbarsMM(formula) ## this comes from lme4 and converts (...|...) terms to some "+" form 
  environment(frame.form) <- environment(formula)
  mf_call <- call("model.frame", data=data, formula=frame.form,  drop.unused.levels=TRUE) # language object reused later
  full_frame <- eval(mf_call) ## data.frame with all vars required for fixef and for ranef, and a "terms" attribute
  res <- list(mf=full_frame)
  #
  Y <- model.response(full_frame, "any")
  if ( ! is.null(Y)) { ## exclude this cases which occurs when no LHS formula is handled
    if (is.factor(Y)) { 
      Y <- (Y != levels(Y)[1L]) ## cf ?glm: 'success' is interpreted as the factor not having the first level
    } ## indeed as in  binomial()$initialize
    if (NCOL(Y)==1L) { ## test on NCOL excludes the case cbind-LHS -> matrix Y (with any array1d stuff already erased)
      Y <- as.matrix(Y) ## to cope with array1d, and see $y <- ... code in .preprocess
    }
  }
  res$Y <- Y
  #
  # compare to .calc_newFrames_fixed() to keep changes in code consistent
  fixef_off_form <- .stripRanefs_(formula) 
  if (inherits(fixef_off_form, "formula")) {
    # to keep info about offset in $fixef_off_terms
    mf_call$formula <- fixef_off_form
    feo_frame <- eval(mf_call)
    res$fixef_off_terms <- attr(feo_frame, "terms") ## added 2015/12/09 useful for .calc_newFrames()
    #
    fixef_form <- .stripOffset_(fixef_off_form) # formula if something remains (e.g, ". ~ 0") after the offset has been removed
    if (inherits(fixef_form, "formula")) { # something remains after the offset has been removed (e.g., explicit '0')
      ## Produce X design matrix:
      fe_call <- mf_call
      fe_call$formula <- fixef_form
      fe_frame <- eval(fe_call)
      fixef_terms <- attr(fe_frame, "terms")
      if (is.null(fixef_terms)) {
        message("Note: formula without explicit fixed effects is interpreted
              as formula without fixed effects [i.e., as .~ 0 + (random effect)].
              Use ~ 1 + (random effect) rather than ~ (random effect) to include an Intercept.")
      }
      if (fixef_terms[[length(fixef_terms)]]==0L) { ## check that the fixef are only an explicit '0' (odd that it compares to 0L, but it does)
        res$X <- matrix(nrow=nrow(full_frame),ncol=0L) ## model without fixed effects, not even an Intercept 
      } else {
        res$X <- model.matrix(fixef_terms, full_frame, contrasts.arg = NULL) ## always valid, but slower
      } 
      res$fixef_levels <- .getXlevels(fixef_terms, fe_frame) ## added 2015/12/09 useful for .calc_newFrames()
    } else { ## only an offset in formula, not even an explicit 0: .stripOffset_(fixef_off_form) produced a 'name'
      message("Note: formula without explicit fixed effects is interpreted
              as formula without fixed effects [i.e.,  as .~ 0 + offset + (random effect)].
              Use ~ 1 + offset + (random effect) rather than ~ offset + (random effect) to include an Intercept.")
      res$X <- matrix(nrow=nrow(full_frame),ncol=0L) ## model without fixed effects, not even an Intercept 
    }
  } else {
    res$X <- matrix(nrow=nrow(data), ncol=0L) ## Not using NROW(Y) which is 0 if formula has no LHS errr... RHS ?
  }
  ####### Then constructs the design X by evaluating the model frame (fe) with fe$formula <- fixef.form
  storage.mode(res$X) <- "double" ## otherwise X may be logi[] rather than num[] in particular when ncol=0
  return(res) # full_frame, fixef_off_terms, fixef_levels, X, Y
}

.calc_newpredvars <- function(oldterms, formula) {
  predvars <- attr(oldterms, "predvars") ## coming from attr(HLframes$mf, "terms")
  if ( ! is.null(predvars)) { 
    newterms <- terms(formula)
    newtermnames <- rownames(attr(newterms,"factors"))
    oldtermnames <- rownames(attr(oldterms,"factors"))
    textnewform <- .DEPARSE(formula)
    newpredvars_factors <- character(length(newtermnames))
    for (it in seq_len(length(newtermnames))) {
      oldpos <- which(oldtermnames==newtermnames[it])
      newpredvars_factors[it] <- .DEPARSE(predvars[-1L][[oldpos]])
    }
    #
    newpredvars_offset <- character(length(attr(newterms,"offset"))) ## template with empty string(s) ""
    offsetpos <- attr(newterms,"offset")
    for (it in seq_len(length(offsetpos))) newpredvars_offset[it] <- .DEPARSE(attr(newterms,"variables")[[offsetpos[it]+1L]])
    #
    newpredvars <- unique(c(newpredvars_factors,newpredvars_offset))
    newpredvars <- paste0("list(",paste(newpredvars,collapse=","),")")
    return(parse(text=newpredvars)) ## contains poly(., coefs) information,
  } else return(NULL)
}

.calc_newFrames_fixed <- function (formula, data, fitobject, need_allFrames=TRUE) {
  ## X may or may not contain offset info, which should not be used (see .newEtaFix()) 
  #  but fixef_mf should contain such info bc .newEtaFix calls off <- model.offset( newMeanFrames$mf)
  fixef_off_terms <- .get_fixef_off_terms(fitobject) 
  if (is.null(formula)) {
    X <- matrix(nrow=nrow(data),ncol=0L) ## model without fixed effects, not even an Intercept 
    fixef_mf <- NULL 
  } else { 
    fixef_off_form <- .stripRanefs_(formula) 
    if (inherits(fixef_off_form, "formula")) {
      if (is.character(formula[[2L]])) fixef_off_form <- fixef_off_form[-2L] ## something like ".phi" ....
      Terms <- terms(fixef_off_form)
      Terms <- stats::delete.response(Terms)
      attr(Terms,"predvars") <- .calc_newpredvars(fixef_off_terms, fixef_off_form) ## for poly()
      fixef_form <- .stripOffset_(fixef_off_form) # formula if something remains after the offset has been removed
      if ( ! inherits(fixef_form, "formula")) { ## only an offset in formula, not even an explicit 0: .stripOffset_(fixef_off_form) produced a 'name'
        attr(Terms,"intercept") <- 0L # removes offset that terms() assumes if there is no explicit '0'.
      }
      # handles offset:  (without the small shortcut used in .HLframes())
      fixef_mf <- model.frame(Terms, data, xlev = fitobject$HLframes$fixef_levels) ## xlev gives info about the original levels
      # :here for a poly(age,.) Terms and age=Inf in the 'data', fixef_mf had zero rows and linkinv will fail on numeric(0) 'eta'
      if (nrow(fixef_mf)!=nrow(data)) {
        if (any(pb <- which( ! sapply(lapply(data,is.finite), all)))) {
          stop(paste0("NA/NaN/Inf in 'data' for fixed-effects prediction: check variable(s) '", paste(names(pb), collapse="', '"),"'."))
        } else stop("nrow(fixef_mf)!=nrow(data) for undetermined reason") 
      }
      X <- model.matrix(Terms, fixef_mf, contrasts.arg=attr(fitobject$X.pv,"contrasts")) 
      ## : original contrasts definition is used to define X cols that match those of original X, whatever was the contrast definition when the model was fitted
    } else {
      fixef_mf <- NULL
      X <- matrix(nrow=nrow(data), ncol=0L) ## Not using NROW(Y) which is 0 if formula has no LHS
    }
  }
  storage.mode(X) <- "double" ## otherwise X may be logi[] rather than num[] in particular when ncol=0
  return(list(X = X, mf = fixef_mf)) 
}

.calc_newFrames_ranef <- function (formula, data, fitobject) {
  formula <- .asNoCorrFormula(formula) ## strips out the correlation information, retaining the ranefs as (.|.)
  if (is.character(formula[[2L]])) formula <- formula[-2L] ## something like ".phi" ....
  plusForm <- .subbarsMM(formula) ## this comes from lme4 and converts (.|.) terms to (.+.) form 
  environment(plusForm) <- environment(formula)
  Terms <- terms(plusForm) ## assumes an Intercept implicitly
  Terms <- stats::delete.response(Terms)
  attr(Terms,"predvars") <- .calc_newpredvars(fitobject$HLframes$all_terms, Terms) ## for poly in ranefs ? which never worked
  mf <- model.frame(Terms, data, drop.unused.levels=TRUE) 
  return(list(mf = mf))
}
