.ULI <- function(...) { # Try to rmeove its last uses ? _F I X M E_
  redondGeo <- cbind(...) ## NOT always a matrix: if ... is a model frame, there may be $ poly(x, 2, raw = TRUE): 'poly' num  inside
  if (ncol(redondGeo)==0L) return(rep(1L,nrow(redondGeo))) ## .generateInitPhi [removed] with constant predictor led here
  if (nrow(redondGeo)==1L) return(1L) ##  trivial case where the following code fails
  ## redondGeo is NOT always a matrix: if ... is a model frame, there may be <$ poly(x, 2, raw = TRUE): 'poly' num>  inside
  ## poly may be a two-col matrix and ncol(redondGeo) is *1*. Then
  # for (colit in seq_len(ncol(redondGeo))) redondGeo[,colit] <- factor(redondGeo[,colit],labels="") 
  ## fails in poly() case...
  redondFac <- apply(redondGeo,2L,factor,labels="") # not cute use of labels... 
  redondFac <- .pasteCols(t(redondFac)) # apply(redondFac,1L,paste,collapse=":") ## paste factors
  #redondFac <- as.character(factor(redondFac))
  uniqueFac <- unique(redondFac) ## seems to preserve order ## unique(<integer>) has unambiguous behaviour
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

# cf .remove_all_fixef() for opposite effect  
.stripRanefs_ <- function (term) { ## compare to lme4:::nobars_ ; 'term is formula or any term, recursively
  if ( ! ("|" %in% all.names(term))) return(term) ## no bar => not a ranef term
  termname <- term[[1L]]
  if (is.call(term) && termname == as.vector("|", "symbol")) return(NULL) # (.|.) ranef
  if (as.vector(termname, "character") %in% .spaMM.data$keywords$all_ranefs) return(NULL) # rather explicit
  if (termname == as.vector("multIMRF")) return(NULL) # "multIMRF" is not (formally declared as) a ranef keyword so special handling
  if (length(term) == 2L) { # 
    # this would strip Matern(1 | longitude + latitude),  that has 2 elements  Matern  and  1 | longitude + latitude,
    # if this term was not already caught by check against .spaMM.data$keywords$all_ranefs.
    # This operated before .spaMM.data$keywords$all_ranefs was introduced to deal with user-defined ranefs with additional params (length>2) 
    nb <- .stripRanefs_(term[[2L]])
    if (is.null(nb)) 
      return(NULL)
    term[[2L]] <- nb
    return(term)
  }
  nb2 <- .stripRanefs_(term[[2L]])
  nb3 <- .stripRanefs_(term[[3L]])
  if (is.null(nb2)) 
    return(nb3)
  if (is.null(nb3)) 
    return(nb2)
  term[[2L]] <- nb2
  term[[3L]] <- nb3
  term
}

# Correct, used until v3.6.15, but operations on 'symbols' in .subbarsMM() are now sufficient => obsolete
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
  # Handling the 'virtual' factor represented by an mv() expression; necess for .preprocess -> .GetValidData_info() as the virtual levels are not in the data.
  aschar <- gsub("mv\\([^|]+", "1 ", aschar, fixed=FALSE) # finds mv(.... up to the closing ')' included and replaces it by '1' 
  as.formula(aschar, env=environment(formula))
}

## lme4::subbars "Substitute the '+' function for the '|' function in a mixed-model formula, recursively (hence the argument name term). This provides a formula suitable for the current model.frame function."
## Original version shows handling of '||'
.subbarsMM <- function (term) {   
  if (is.name(term) || !is.language(term)) return(term) #  single variable (such as the LHS of a formula)
  if (as.vector(term[[1L]], "character") %in% .spaMM.data$keywords$all_keywords)  # handles "multIMRF" and "mv", not simply ranefs
             return(.subbarsMM(term[[2L]]))
  if (length(term) == 2L) { 
    # an unregistered corrFamily reaches here ARp(1 | time) would become ARp(1 + time)
    # but also other syntaxes understood by R such as  I(prop.ag/10). So it's not easy to distinguish valid and invalid terms.
    # The unregistered cF generates an error in .GetValidData_info -> model.frame... which catches the list returned by execution of the constructor
    term[[2L]] <- .subbarsMM(term[[2L]])
    return(term)
  }
  # stopifnot(length(term) >= 3) ## inefficient
  if (is.call(term) && term[[1L]] == as.vector("|", "symbol"))  term[[1L]] <- as.vector("+", "symbol")
  # if term is   1 + <other formula terms>, term[[1]] is  +  and we reach here:
  for (j in 2L:length(term)) term[[j]] <- .subbarsMM(term[[j]])
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

.param.inla_dgT2dgC <- function(param.inla) {
  if (inherits(M0 <- param.inla$M0,"dgTMatrix")) {
    param.inla$M0 <- as(M0,"CsparseMatrix")
  }
  if (inherits(M1 <- param.inla$M1,"dgTMatrix")) {
    param.inla$M1 <- as(M1,"CsparseMatrix")
  }
  if (inherits(M2 <- param.inla$M2,"dgTMatrix")) {
    param.inla$M2 <- as(M2,"CsparseMatrix")
  }
  param.inla
}

.process_IMRF_bar <- function(term, env) {
  pars <- paste0(paste0(names(term),"=",term)[-c(1,2)], collapse=",") 
  pars <- try(eval(parse(text = paste("list(",pars,")")), envir=env))
  if (inherits(pars,"try-error")) {
    stop("Hmmm... it looks like you should put some object(s) in control.HLfit$formula_env.")
  }
  if ( ! is.null(pars$spde)) {
    warning("'spde' argument is obsolete. Use 'model' instead.")
    pars$model <- pars$spde
    pars$spde <- NULL
  }
  if ( ! is.null(pars$model)) {
    pars$no <- FALSE
    # it's hard to guess the alpha parameter from the object!
    # diagnose the object
    if (pars$model$mesh$manifold=="R1") {dim_ <- 1} else dim_ <- 2
    # test for # inla.spde.matern with default parameters B.tau and B.kappa
    is_default_spde2_matern <- (diff(range((pars$model$param.inla$BLC[c("tau.1","kappa.1"),] - matrix(c(0,0,1,0,0,1),ncol=3))))==0)
    if (is_default_spde2_matern) {
      # 'pars' looks like the result of inla.spde2.matern() with default B.tau
      if (pars$model$param.inla[["B0"]][1,3]!=0) {
        pars$SPDE_alpha <- 2+pars$model$param.inla$B0[1,3] # 2+ seem OK in 1D
      } else {
        # ' only the case d=2 and alpha= 1 or 2' is implemented in spaMM => only in that case the following code is true
        pars$SPDE_alpha <- pars$model$param.inla[["B1"]][1,3]
        if ( ! (pars$SPDE_alpha  %in% c(1,2))) stop("Unrecognized model from inla.spde2. Contact the spaMM maintainer to extend spaMM.")
      }
    } else {
      #  'pars' looks like the result of inla.spde2.pcmatern() (which uses non-default B.tau)
      pars$SPDE_alpha <- dim_/2 + pars$model$param.inla[["BLC"]][["tau.1",2]] 
    }
    # print((c(dim_,pars$SPDE_alpha)))
  } else if (is.null(pars$no)) pars$no <- TRUE
  if (is.null(pars$ce)) pars$ce <- TRUE
  pars$model$param.inla <- .param.inla_dgT2dgC(pars$model$param.inla)
  attr(term,"type") <- structure("IMRF", pars=pars) # pars as attr to type avoid problems in building the return value.
  return(term) # (lhs|rhs) (language, with the ()); or character string.
}

## derived from findbars, ... 
# Must add a type attr to each formula term: 
.parseBars <- function (term, # term or formula
                        env=environment(term) # parent call typically has a formula as argument, so the default env is that of the formula,
                        #                       which is then passed recursively by explicit env=env 
                        ) { 
  if (inherits(term,"formula") && length(term) == 2L) {## corrected 2021/03/07
    resu <- .parseBars(term[[2L]], env=env) ## process RHS
    if (is.call(resu)) { # top call spaMM:::.parseBars(~(1|a)) leads here (no LHS, single ranef). resu is not (yet) a list here;
      #  Then .process_bars => names(barlist) <- seq_along(barlist) on a call leads to nonsense and bug for get_predVar(, re.form=~(1|a))
      # So we need to convert to list.  (testing if resu is a list in wrong: resu may be NULL)
      resu <- list(resu)
      attr(resu,"type") <- attr(resu[[1L]],"type") 
    }
    return(resu)
  }
  if (is.name(term) || !is.language(term)) 
    return(NULL)
  termname <- term[[1L]]
  if (termname == as.vector("(","symbol")) return(.parseBars(term[[2L]], env=env)) ## i.e. (blob) but not ( | ) 
  if (!is.call(term)) stop("term must be of class call")
  if (termname == as.vector("|","symbol")) { ## i.e.  .|. expression
    attr(term,"type") <- rep("(.|.)",length(paste(c(term)))) ## paste is the simplest way to count terms (rather than elements of a single term)
    return(term)
  } ## le c() donne .|. et non |..
  if (as.vector(termname,"character") %in% .spaMM.data$keywords$all_ranefs) {
    attr(term,"type") <- as.character(term[[1L]])
    if (termname == as.vector("IMRF","symbol")) {
      return(.process_IMRF_bar(term, env=env)) # 
    } else return(term) 
  } else if (length(term) == 2L) return(NULL) # this occurs. Maybe the extra arguments in IMRF? 
  # # terms are of class call, so expectedly ',' separates different elements in a term
  # # E.g., if not handled specifically, <somekeyword>(1,2|.) would be parsed as term[[2]]=1, term[[3]]=2|batch
  # # That means that the test length(term) == 2 is not appropriate for terms with some ','
  # # <somekeyword>(1 2|.) does not work natively: R doesn't know what to do about the 2 (unexpected numeric constant).
  # if (term[[1L]] == as.name("multIMRF") ) { # useful for .findSpatial() before preprocessing, as in filled.mapMM()
  #   return(term) 
  # }

  # This point may not be reached by the outer .parseBars() call, so the following is not consistently terminal code.
  bT2 <- .parseBars(term[[2L]], env=env)
  bT3 <- .parseBars(term[[3L]], env=env)
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
  if ( ! is.null(barlist)) {
    names(barlist) <- seq_along(barlist)
    # NULL barlist is a result from a previous.parseBars / .process_bars call !   
    if (expand) {
      barlist <- .spMMexpandSlash(barlist)
    }
    attr(barlist,"type") <- unlist(lapply(barlist,attr,which="type")) # before the element attributes are hacked by .lhs_rhs_bars()
  }
  if (which.=="exp_ranef_terms") {
    return(.lhs_rhs_bars(barlist)) # keeps the "type" attribute of the list # but removes the type keyword from the term
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

#mv <- function(...) return(unlist(list(...)))

if (FALSE) { # seems correct, but ultimately not needed
  .asNoPoly <- function(formula) {
    aschar <- .DEPARSE(formula)
    while(TRUE) {
      newaschar <- gsub("poly\\(cbind\\((.+?),(.+?)\\)(.+?)\\)", "poly(cbind(\\1+\\2)\\3)", aschar, fixed = FALSE)
      if (newaschar==aschar) break
      aschar <- newaschar 
    } # replaces ',' separating variables in cbind() by '+'
    aschar <- gsub("poly\\(cbind\\((.+?)\\),(.+?)\\)", "poly\\(\\1\\)", aschar, fixed = FALSE) # keeps only the argument within cbind()
    aschar <- gsub("poly\\(([^,]+?),(.+?)\\)", "\\1", aschar, fixed = FALSE) # or keeps only the variable
    as.formula(aschar, env=environment(formula))
  }
}

.GetValidData_info <- function(formula, resid.formula=NULL, data, prior.weights,
                                   ... # to ignore a lot of other arguments of the call from which the .GetValidData_info() call is constructed.
) {
  ## As in .preprocess_pw(): 
  # if user entered prior.weights=pw, model.frame will handle substitute(prior.weights), a 'symbol' (or 'name'), we substitute
  # if user entered prior.weights=1/varx, model.frame will handle substitute(prior.weights) also.
  # if user entered quote(1/varx), model.frame will handle eval(substitute(prior.weights)), a 'call'.
  prior.weights <- substitute(prior.weights)
  if ( (inherits(prior.weights,"call") && prior.weights[[1L]] == "quote") ) prior.weights <- eval(prior.weights)
  #
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
  #formula <- .asNoCorrFormula(formula) ## removes spatial tags
  frame.form <- .subbarsMM(formula) ## this converts (...|...) terms to some "+" form
  if (!is.null(resid.formula)) { 
    frame.resid.form <- .subbarsMM(resid.formula) 
    frame.form <- paste(.DEPARSE(frame.form),"+",.DEPARSE(frame.resid.form[[2]])) ## only good to select rows 
    ## (eg, may have offset terms from the two formulas, summed by model.offset() ! )
  }
  check <- grep('$',frame.form,fixed=TRUE)
  if (length(check)) {
    message("'$' detected in formula: suspect and best avoided. 
            In particular, one should never need to specify the 'data' in the 'formula'. 
            See help('good-practice') for explanations.")
  }
  check <- grep('$',frame.form,fixed=TRUE)
  if (length(check)) {
    message("'$' detected in formula: suspect and best avoided. 
            In particular, one should never need to specify the 'data' in the 'formula'. 
            See help('good-practice') for explanations.")
  }
  check <- grep('c(',frame.form,fixed=TRUE)
  if (length(check)) {
    if ((inherits(frame.form,"formula") && check[1L]==length(frame.form)-1L) ||
        is.character(frame.form) && check[1L]==1L
    )
      warning("'c(' detected in formula: did you mean cbind() for binomial response or for poly()?", immediate.=TRUE)
    # warning useful as eval(mf) generates an obscure error.
  }
  frame.form <- as.formula(frame.form)
  environment(frame.form) <- envform
  # model.frame calls terms(formula, data) [if the 'formula' is not already terms]. 
  # Here we call terms(formula) without data to detect use of '.'
  chkterms <- try(terms(frame.form), silent=TRUE)
  if (inherits(chkterms,"try-error") && length(grep("'.'",attr(chkterms,"condition")$message))) {
    warning("It looks like there is a '.' in the RHS of the formula.\n Fitting may be successful, but post-fit functions such as predict() will fail.", immediate.=TRUE)
    call2mf <- call("model.frame", data = data, formula= frame.form, drop.unused.levels = TRUE, weights=prior.weights)
  } else call2mf <- call("model.frame", data = data, formula= chkterms, drop.unused.levels = TRUE, 
                         weights=prior.weights # we don't want to evaluate it since we want its variables to be included in the data.
                         )
  resu <- try(eval(call2mf), silent=TRUE) ## # "directly calling model.frame failed to handle the prior weights (tried again 03/2022)."
                                          ## => A bit opaque. Problem handling promises in a programming context where prior weights
                                          ## are not specified by the user but by spaMM procedures ?
  if (inherits(resu,"try-error")) { # => try to diagnose a case that is otherwise difficult for even the developer to guess.
    verif <- try(eval(prior.weights,data),silent=TRUE)
    if (inherits(verif,"try-error")) {
      stop("All variables should be in the 'data', including those for prior weights.")
    } else return(verif) 
  }
  return(list(rownames=rownames(resu), weights=model.weights(resu))) # so that valid rows and weights always have the same length.
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

.sanitize_Y <- local({
  #int_warned <- FALSE
  function(y, famfam) {
    if ( famfam %in% c("binomial","poisson","COMPoisson","negbin1","negbin2", "betabin")) { # COUNTS
      ## the response variable should always be Counts
      safe_y <- as.integer(y+0.5) # non-negative values only # non-array from array, hence:
      if (NCOL(y)) dim(safe_y) <- dim(y)
      if (max(abs(y-safe_y))>1e-05) {
        anynegy <- any(y<0L)
        if (anynegy) { # at this point, there are 'large' negative values
          stop(paste0("negative values not allowed for the '",famfam,"' family"))
        } else stop("response variable should be integral values.")
      } else {
        # if ( ! int_warned) {
        #   int_warned <<- TRUE
        #   message("Response converted to integer for integral-response families")
        # }
        y <- safe_y # silent sanitizing # tiny negative values would stop() later
      }
    } else if (famfam=="Gamma") {
      Gamma_min_y <- .spaMM.data$options$Gamma_min_y
      is_low_y <- (y < Gamma_min_y)
      if (any(is_low_y)) {
        #y[which(is_low_y)] <- Gamma_min_y
        warning(paste0("Found Gamma response < (Gamma_min_y=",Gamma_min_y,") . Troubles may happen."), immediate. = TRUE)
      }
      is_high_y <- (y > 1/Gamma_min_y)
      if (any(is_high_y)) {
        #y[which(is_low_y)] <- Gamma_min_y
        warning(paste0("Found Gamma response > (1/Gamma_min_y=",1/Gamma_min_y,") . Troubles may happen."), immediate. = TRUE)
      }
    }else if (famfam=="beta_resp") {
      if (any(y < 0 | y > 1)) {
        stop("Found Beta responses outside valid (0,1) range.")
      }
      beta_min_y <- .spaMM.data$options$beta_min_y
      if (any(y < beta_min_y | y > 1 - beta_min_y)) {
        #y[which(is_low_y)] <- Gamma_min_y
        warning(paste0("Found Beta responses close to 0 or 1 by less than ",beta_min_y,") . Troubles may happen."), immediate. = TRUE)
      }
    }
    return(y)
  }
})

.get_Y <- function(full_frame, famfam) {
  Y <- model.response(full_frame, "any")
  if ( ! is.null(Y)) { ## exclude this cases which occurs when no LHS formula is handled
    if (is.factor(Y)) { 
      Y <- (Y != levels(Y)[1L]) ## cf ?glm: 'success' is interpreted as the factor not having the first level
    } ## indeed as in  binomial()$initialize
    if (NCOL(Y)==1L) { ## test on NCOL excludes the case cbind-LHS -> matrix Y (with any array1d stuff already erased)
      Y <- as.matrix(Y) ## to cope with array1d, and see $y <- ... code in .preprocess
      respname <- colnames(full_frame)[1]
    } else respname <- colnames(full_frame[[1]])[1] ## binomial, full_frame[[1]] is a matrix with cols for npos and nneg
    Y <- .sanitize_Y(Y, famfam) # sanitize Y, rather than processed$y which is not used by .get_inits_by_xLM()
    Y <- structure(Y, respname=respname) # respname drops automatically when processed$y is evaluated
  }
  Y
}

## cf model.frame.default from package stats mais ne traite pas les effets alea !
.get_terms_info <- function (formula, data, famfam, weights=NULL) {
  if (is.character(formula[[2]])) { # Interestingly, not needed for .getValidData(), as the latter is not called on the resid.formula...
    validname <- .find_validname(formula, data)
    data[validname] <- 1 ## adds a column $phi of 1 
    formula[[2]] <- as.name(validname) ## now the formula is standard
  }
  ####### first construct a mf for all variables required by the formula (NA would be removed here if they had not been by a previous call to validData)
  frame.form <- .subbarsMM(formula) ## converts (...|...) terms to some "+" form, from a similar concept in lme4 
  environment(frame.form) <- environment(formula)
  mf_call <- call("model.frame", data=data, formula=frame.form,  drop.unused.levels=TRUE, weights=weights) # language object reused later
  full_frame <- eval(mf_call) ## data.frame with all vars required for fixef and for ranef, and a "terms" attribute
  # 
  res <- list(mf=full_frame)   ## MODEL FRAME 
  #
  res$Y <- .get_Y(full_frame, famfam)
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
      fe_call <- mf_call # uses the data, not a model frame !
      fe_call$formula <- fixef_form
      fe_frame <- eval(fe_call)
      res$fixef_terms <- fixef_terms <- attr(fe_frame, "terms") ## return fixef_terms for processing data with updated response.
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
      #
      specials_levels <- list()
      labs <- names(fe_frame)
      for (it in seq_along(fe_frame)) specials_levels[[labs[it]]] <- attr(fe_frame[[labs[it]]],"spec_levs")
      res$specials_levels <- specials_levels
      #
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
  return(res) # $mf=full_frame, $fixef_off_terms, $fixef_levels, $X, $Y
}

.calc_newpredvars <- function(fitobject, formula) {
  predvars <- .get_from_data_attrs(fitobject, which="fixefpredvars") 
  if ( ! is.null(predvars)) { 
    newterms <- terms(formula)
    newtermnames <- rownames(attr(newterms,"factors")) # or paste() or .DEPARSE() on attr(,"variables") ?
    oldtermnames <- .get_from_data_attrs(fitobject, which="fixefvarnames")
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
