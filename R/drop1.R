# remove specific fixefs. 
# cf buildmer::remove.terms for another approach, not suitable here 
# (here: no concern with marginality not problem with _apparent_ marginality).
#
# stats::terms.formula() first drops parentheses from simple (.|.) ranefs (pfff) but not prefixed ones.
#   It then calls a local fixFormulaObject(terms) functions that puts the parentheses back...
# We need to put back them but not on the prefixed terms.
#
# There are distinct issues using terms.formula(): 
# (1) it does not handle correctly random effect
# cf terms( y ~ foo(1|x) + (1|id) + a + b) => "foo(1 | x)" "1 | id"     "a"          "b"  
# I defined a term.HLfit() function that does not use the formula(object)
# Here I use terms.formula() but I must call .stripRanefs() first.
# (2) it returns an object inheriting from class "formula" but `[``
# does not work as form a basic formula.
#
# One should check that any change is backward-compatible with 'mandrills' tests.
#
.remove_fixef <- function (formula, # ideally formula(<fit object>) incl. possible offsets
                           remove, # regressors to be removed (>1 is possible, as vectors of strings) OR formula
                                   # argument better obtained from another formula, as attr(terms(.),"term.labels")
                           nofixef= NULL, keep_intercept=TRUE, keep_offset=TRUE) {
  if (is.null(nofixef)) nofixef <- .stripFixefs(formula, keep_offset)
  # Case where 'remove' is a formula: convert it to regressor names.
  if (inherits(remove,"formula")) {
    if (length(remove)==3L) remove <- remove[-2L]
    remove <- .stripRanefs_(remove)
    if (is.null(remove)) {
      warning("'remove' argument does not contain fixed effects.", immediate. = TRUE)
      return(formula)
    }
    remove <- attr(terms(remove),"term.labels")
  }

  if (inherits(formula,"terms")) {
    # Then formula[-2] does not give the desired result...
    stop("'formula' argument is already a 'terms' object.")
  } else {
    if (length(formula)==3L) {
      termobj <- .stripRanefs_(formula[-2])
    } else termobj <- .stripRanefs_(formula)
  }
  if (is.null(termobj)) {
    warning("'formula' argument does not contain fixed effects.", immediate. = TRUE)
    return(formula)
  }
  termobj <- terms(termobj)
  termlabs <- attr(termobj,"term.labels")
  ind <- grep("|", termlabs, fixed = TRUE)
  if (length(ind)) termlabs[ind] <- .process_bars(formula) # .process_bars() -> .as_char_bars() to handle parentheses.
  
  remove_list <- strsplit(remove,":")
  for (it in seq_along(remove_list))  { 
    nterms <- length(termlabs)
    remove_i <- remove_list[[it]]
    if (length(remove_i)>1L) {
      remove_i <- .gen_all_perms(remove_i)
      remove_i <- sapply(remove_i, paste, collapse = ":")
    } 
    termlabs <-  setdiff(termlabs, remove_i)
    if (length(termlabs) == nterms) 
      warning(paste("term(s)",remove[it],
                    "aleady removed or not in terms implied by 'formula'."),
              immediate. = TRUE)
  }
  
  if (length(termlabs)) {
    newform <- paste(termlabs,collapse=" + ")
    if (keep_intercept && attr(termobj,"intercept")) newform <- paste("1", newform, sep="+")
  } else if (keep_intercept && attr(termobj,"intercept")) {
    newform <- "1"
  } else newform <- NULL
  
  # => This newform is only a charstring for the RHS, for fixed effects. Add random effects:
  if ( ! is.null(nofixef) &&
      ( ! (rhs_nofixef <- deparse(nofixef[[length(nofixef)]]))=="0") # nofixef != ~ 0
     ) newform <- paste(c(newform, rhs_nofixef), collapse="+")
  if (is.null(newform)) newform <- "0"
  # Convert to length-2 or length-3 formula:
  newform <- paste("~", newform)
  if (length(formula)==3L) newform <- paste(deparse(formula[[2L]]), newform)
  as.formula(newform)
}

# Create API version while keeping internal version used in possibly public code. 
remove_fixef <- .remove_fixef 

# Minimal marginality check from fixed effects. See buildmer::remove.terms for more general stuff
.is_marginal <- function(remove, have) {
  forbidden <- if (!all(have == "1")) {"1"} else NULL
  for (x in have) {
    x.star <- gsub(":", "*", x)
    partterms <- attr(terms(stats::as.formula(paste0("~", x.star))), "term.labels")
    forbidden <- c(forbidden, partterms[partterms != x])
  }
  !remove %in% forbidden
}

.preprocess_scope <- function(scope=NULL, object, check, tl= attr(terms(object), "term.labels")) {
  if (is.null(scope)) {
    scope <- drop.scope(object)
    if (is.null(check)) check <- FALSE
  } else {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(object, scope)), 
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L)) 
      stop("scope is not a subset of term labels")
    if (is.null(check)) check <- TRUE
  }
  if (check) {
    marg_ok <- .is_marginal(scope, tl)
    if (any( ! marg_ok)) {
      message("Some test(s) removed as not satisfying marginality condition. See drop1.HLfit()'s 'check_marg' argument.")
      scope <- scope[marg_ok]
    }
  }
  scope
}



.drop1.lm <- function(object, scope=NULL, check=NULL, scale = 0, all.cols = TRUE, test = c("none", 
                                                              "Chisq", "F"), k = 2, ...) {
  test <- match.arg(test)
  if (test=="F") {
    w <- weights(object, type="prior")
    ssr <- sum(if (is.null(w)) residuals.HLfit(object,type="response")^2 else w * residuals.HLfit(object,type="response")^2)
    mss <- sum(if (is.null(w)) fitted(object)^2 else w * fitted(object)^2)
    if (ssr < 1e-10 * mss) 
      warning("F-tests on an essentially perfect fit are unreliable")
  }
  x <- model.matrix(object)
  offset <- model.offset(model.frame(object))
  iswt <- !is.null(wt <- weights(object, type="prior"))
  n <- nrow(x)
  tl <- attr( terms(object), "term.labels")
  scope <- .preprocess_scope(scope, object=object, check=check, tl=tl)
  ndrop <- match(scope, tl)
  ns <- length(scope)
  rdf <- df.residual(object)
  chisq <- deviance(object)
  dfs <- RSS <- numeric(ns)
  y <- object$y
  na.coef <- seq_along(object$coefficients)[!is.na(object$coefficients)]
  asgn <- attr(x, "assign")
  for (i in seq_len(ns)) {
    ii <- seq_along(asgn)[asgn == ndrop[i]]
    jj <- setdiff(if (all.cols) { seq(ncol(x)) } else na.coef, ii)
    z <- if (iswt) {
      lm.wfit(x[, jj, drop = FALSE], y, wt, offset = offset)
    } else stats::lm.fit(x[, jj, drop = FALSE], y, offset = offset)
    dfs[i] <- z$rank
    oldClass(z) <- "lm"
    RSS[i] <- deviance(z)
  }
  #
  scope <- c("<none>", scope)
  dfs <- c(object$dfs$pforpv, dfs)
  RSS <- c(chisq, RSS)
  if (scale > 0) {
    aic <- RSS/scale - n + k * dfs
  } else aic <- n * log(RSS/n) + k * dfs
  dfs <- dfs[1L] - dfs
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, `Sum of Sq` = c(NA, RSS[-1L] - 
                                                RSS[1L]), RSS = RSS, AIC = aic, row.names = scope, check.names = FALSE)
  if (scale > 0) 
    names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
  #
  if (test == "Chisq") {
    dev <- aod$"Sum of Sq"
    if (scale == 0) {
      dev <- n * log(RSS/n)
      dev <- dev - dev[1L]
      dev[1L] <- NA
    } else dev <- dev/scale
    df <- aod$Df
    nas <- !is.na(df)
    dev[nas] <- pchisq(dev[nas], df[nas], lower.tail = FALSE)
    aod[, "Pr(>Chi)"] <- dev
  }
  else if (test == "F") {
    dev <- aod$"Sum of Sq"
    dfs <- aod$Df
    rdf <- df.residual(object)
    rms <- aod$RSS[1L]/rdf
    Fs <- (dev/dfs)/rms
    Fs[dfs < 1e-04] <- NA
    P <- Fs
    nas <- !is.na(Fs)
    P[nas] <- pf(Fs[nas], dfs[nas], rdf, lower.tail = FALSE)
    aod[, c("F value", "Pr(>F)")] <- list(Fs, P)
  }
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)), 
            if (scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}


.drop1.glm <- function (object, scope=NULL, check=NULL, scale = 0, test = c("none", "Rao", "LRT", 
                                             "Chisq", "F"), k = 2, ...) {
  test <- match.arg(test)
  doscore <- !is.null(test) && test == "Rao"
  
  if (test == "Chisq") test <- "LRT"
  x <- model.matrix(object)
  n <- nrow(x)
  termsobj <- terms(object)
#  tl <- attr(termsobj, "term.labels")
  scope <- .preprocess_scope(scope, object=object, check=check,  tl=attr(terms(termsobj), "term.labels"))
  ns <- length(scope)
  rdf <- df.residual(object)
  chisq <- logLik(object)
  dfs <- aics <- dev <- numeric(ns)
  score <- numeric(ns)
  y <- object$y
  if (is.null(y)) {
    y <- model.response(model.frame(object))
    if (!is.factor(y)) 
      storage.mode(y) <- "double"
  }
  wt <- object$prior.weights 
  form <- formula(object)
  nofixef <- .stripFixefs(form, keep_offset=TRUE)
  for (i in seq_along(scope)) {
    # newform <- drop.terms(termsv, drop_ids[i]) # fails when there is a single term
    newform <- .remove_fixef(form, scope[i], nofixef=nofixef, keep_intercept = TRUE)
    refit <- update(object, formula.= newform) 
    dfs[i] <- refit$dfs$pforpv
    aics[i] <- AIC(refit,verbose=FALSE)[[1]] # from p_v, even for REML fits
    dev[i] <- deviance(refit)
    if (doscore) {
      r <- residuals(refit, type="working")
      w <- weights(refit, type="working")
      zz <- glm.fit(x, r, w)
      score[i] <- zz$null.deviance - zz$deviance
    }
  }
  scope <- c("<none>", scope)
  dfs <- c(object$dfs$pforpv, dfs)
  dev <- c(deviance(object),dev)
  if (doscore) {
    score <- c(NA, score)
  }
  fam <- object$family$family
  dispersion <- if (is.null(scale) || scale == 0) 
    dispersion <- residVar(object,"fit")
  else scale
  fam <- object$family$family
  loglik <- if (fam == "gaussian") {
    if (scale > 0) 
      dev/scale - n
    else n * log(dev/n)
  }
  else dev/dispersion
  dfs <- dfs[1L] - dfs
  dfs[1] <- NA
  aics <- c(AIC(object,verbose=FALSE)[[1]],aics) # from p_v, even for REML fits
  aod <- data.frame(#logL = logliks, 
    Df = dfs, Deviance=dev, AIC = aics, row.names = scope, 
                    check.names = FALSE)
  if (all(is.na(aics))) aod <- aod[, -3]
  if (test == "LRT") {
    dev <- pmax(0, loglik - loglik[1L]) # don't try to use the logLik() from a REML fit...
    dev[1L] <- NA
    nas <- !is.na(dev)
    LRT <- if (dispersion == 1) 
      "LRT"
    else "scaled dev."
    aod[, LRT] <- dev
    dev[nas] <- pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
    aod[, "Pr(>Chi)"] <- dev
  }
  else if (doscore) {
    dev <- pmax(0, score)
    nas <- !is.na(dev)
    SC <- if (dispersion == 1) 
      "Rao score"
    else "scaled Rao sc."
    dev <- dev/dispersion
    aod[, SC] <- dev
    dev[nas] <- pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
    aod[, "Pr(>Chi)"] <- dev
  }
  else if (test == "F") {
    fam <- object$family$family
    if (fam == "binomial" || fam == "poisson") 
      warning(gettextf("F test assumes 'quasi%s' family", 
                       fam), domain = NA)
    rms <- dev[1L]/rdf
    dev <- pmax(0, dev - dev[1L])
    dfs <- aod$Df
    rdf <- df.residual(object)
    Fs <- (dev/dfs)/rms
    Fs[dfs < 1e-04] <- NA
    P <- Fs
    nas <- !is.na(Fs)
    P[nas] <- pf(Fs[nas], dfs[nas], rdf, lower.tail = FALSE)
    aod[, c("F value", "Pr(>F)")] <- list(Fs, P)
  }
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)), 
            if (!is.null(scale) && scale > 0) paste("\nscale: ", 
                                                    format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}


.drop1_fallback <- function(object, scope=NULL, check=NULL, check_time, ...) {
  REML <- (object$APHLs$p_v != object$APHLs$p_bv)
  if (REML) warning("LRTs comparing REML fits with different fixed-effect conditions are highly suspect", 
                    immediate.=TRUE)
  form <- formula(object)
  termsobj <- terms(.stripRanefs(form))
  scope <- .preprocess_scope(scope, object=object, check=check, 
                             tl=attr(termsobj, "term.labels"))
  if (is_long <- (fit_time <- how(object, verbose=FALSE)$fit_time)*length(scope)>check_time ) {
    message(paste0("Fitting the original model took ",fit_time,"s and drop1() may take a few times longer."))
  }
  progbar <- (is_long && length(scope>2L))
  basicLRTs <- vector("list", length(scope)) # __F I X M E___ other variants of LR test? bootstrap, etc
  names(basicLRTs) <- scope
  nofixef <- .stripFixefs(form, keep_offset=TRUE)
  if (progbar) cat("\nProgress: ")
  for (i in seq_along(scope)) {
    if (progbar) cat(".")
    # newform <- buildmer::remove.terms(object$predictor, scope[i],check=check) # quite different approach:
    newform <- .remove_fixef(form, scope[i], nofixef=nofixef, keep_intercept = TRUE)
    refit <- update(object, formula.= newform) 
    lrt <- LRT(object, refit, ...) 
    basicLRTs[[i]] <- lrt$basicLRT
  }
  if (progbar) cat("\n")
  basicLRTs <- do.call(rbind,basicLRTs)
  head <- c("Likelihood-ratio tests for single-term deletions", "\nModel:", deparse(formula(object)))
  class(basicLRTs) <- c("anova", "data.frame") # but the names of the data.frame do not those for which stats:::print.anova has specific actions.
  attr(basicLRTs, "heading") <- head
  basicLRTs
}



drop1.HLfit <- function(object, scope=NULL, method="", check_marg = NULL, check_time=60, ...) { # there may also be a 'check_deriv' argument to be passed to as_LMLT
  if (method != "LRT") {
    models <- object$models[c("eta","phi")]
    if (length(models$phi)==1L && models$phi %in% c("phiScal","")) {
      if (models$eta=="etaGLM") { 
        if (object$family$family=="gaussian" && object$family$link=="identity") {
          return(.drop1.lm(object, scope, check=check_marg, ...))
        } else return(.drop1.glm(object, scope, check=check_marg, ...))
      } else if (object$family$family=="gaussian" && object$family$link=="identity") { # LMM
        if (requireNamespace("lmerTest",quietly=TRUE)) {
          scope <- .preprocess_scope(scope, object=object, check=check_marg)
          if ((fit_time <- how(object, verbose=FALSE)$fit_time)*length(scope)>check_time ) {
            message(paste0("Fitting the original model took ",fit_time,"s and drop1() may take a few times longer."))
          }
          lmlt <-  as_LMLT(object, ...)
          return(drop1(lmlt, scope=scope, ...)) 
        } else if ( ! identical(spaMM.getOption("lmerTest_warned"),TRUE)) {
          message("If the lmerTest package were installed, a drop1 single-deletions table could be computed.")
          .spaMM.data$options$lmerTest_warned <- TRUE
        } 
      } 
    }
  }
  # Fallback if no earlier return:
  return(.drop1_fallback(object=object, scope=scope, check=check_marg, check_time=check_time, ...))
}

drop1.LMLT <- function(object, scope, ...) { 
  if (is.null(getClassDef("LMLT", where = .spaMM.data$class_cache, inherits = FALSE))) {
    # is as_LMLT has not been called since restarting R session-> lmerTest presumably not loaded
    if (requireNamespace("lmerTest",quietly=TRUE)) {
      # Hack to define object not of class LMLT (-> infinite recursion) but with similar "contains", using only default coerce methods
      setClass(Class="LMLT", contains = c("LMLTslots","lmerModLmerTest"), where=.spaMM.data$class_cache) 
      setClass(Class="LMLTinternal", contains = c("LMLTslots","lmerModLmerTest"), where=.spaMM.data$class_cache) 
    } else message("If the lmerTest package were installed, a drop1 single-deletions table could be computed.")
  } 
  object <- as(as(object,"LMLTslots"),"LMLTinternal")
  drop1(object, scope, ...) # calling the lmerTest:: method.
}
