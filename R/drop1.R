
# cf .stripRanefs() for opposite effect 
.remove_all_fixef <- function (term, ## 'term is formula or any term, recursively. Offset only handled as formula term
                                     ## so <HLfit obj>$predictor is a suitable object in contrast to term objects. 
                               keep_offset) { 
  
  if (length(term) == 1L) { # '~', or single fixef term
    if (term == as.vector("~", "symbol")) return(term) 
    return(NULL) 
  }
  termname <- term[[1L]]
  #
  if (keep_offset && termname == as.vector("offset")) return(term) 
  #
  # if (is.call(term) && termname == as.vector("|", "symbol")) return(term) # .|. part of a (.|.) ranef 
  if (as.vector(termname, "character") %in% .spaMM.data$keywords$all_ranefs) return(term) # rather explicit
  if (termname == as.vector("multIMRF")) return(term) # "multIMRF" is not (formally declared as) a ranef keyword so special handling
  #
  if (length(term) == 2L) { # 
    # this would strip Matern(1 | longitude + latitude),  that has 2 elements  Matern  and  1 | longitude + latitude,
    # if this term was not already caught by check against .spaMM.data$keywords$all_ranefs.
    if (term[[1L]] == as.vector("(", "symbol") && 
        term[[2L]][[1L]] == as.vector("|", "symbol") ) { #(.|.) 
      return(term)
    } else if (term[[1L]] == as.vector("~", "symbol")) {
      nb <- .remove_all_fixef(term[[2L]], keep_offset=keep_offset)
      if (is.null(nb)) {
        return( ~ 0 )
      } else return(NULL)
    } else return(NULL) # I() or log() ...fixef
  }
  #
  # length(term) == 3L
  nb3 <- .remove_all_fixef(term[[3L]], keep_offset=keep_offset)
  if (term[[1L]] == as.vector("~", "symbol")) {
    attributes(term) <- NULL # remove attributes of the input term object
    if ( is.null(nb3)) nb3 <- 0
    nb2 <- term[[2]] # nothing to remove on a formula LHS
  } else { 
    nb2 <- .remove_all_fixef(term[[2L]], keep_offset=keep_offset)
    if ( is.null(nb3)) return(nb2)
    if ( is.null(nb2)) return(nb3) # say original term was x + (1|g) => without this return, + (NULL, (1|g)) yields (1|g) + (1|g)
  }
  term[[2L]] <- nb2
  term[[3L]] <- nb3
  term
}

# remove specific ranefs; matching is by term labels, 
# A new formula is build from term labels, and from the random effects extracted by .remove_all_fixef()
.remove_fixef <- function (formula, # ideally <HLfit obj>$predictor so that offset is retained
                           remove, # argument better obtained from another formula, as attr(terms(.),"term.labels")
                           nofixef= NULL, keep_intercept=TRUE, keep_offset=TRUE) {
  if (is.null(nofixef)) nofixef <- .remove_all_fixef(formula, keep_offset)
  termobj <- terms(formula)
  termlabs <- attr(termobj,"term.labels")
  termlabs <- setdiff(termlabs, remove)
  if (length(termlabs)) {
    newform <- paste(termlabs,collapse=" + ")
    if (keep_intercept && attr(termobj,"intercept")) newform <- paste("1", newform, sep="+")
  } else if (keep_intercept && attr(termobj,"intercept")) {
    newform <- "1"
  } else newform <- "0"
  if (!is.null(nofixef) &&
      ( ! (rhs_nofixef <- deparse(nofixef[[length(nofixef)]]))=="0")
     ) newform <- paste(newform, "+", rhs_nofixef)
  newform <- paste("~", newform)
  if (length(formula)==3L) newform <- paste(deparse(formula[[2L]]), newform)
  as.formula(newform)
}

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
      message("Some test(s) removed as not satisfying marginalit condition. See 'check' argument.")
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
  nofixef <- .remove_all_fixef(object$predictor, keep_offset=TRUE)
  for (i in seq_along(scope)) {
    # newform <- drop.terms(termsv, drop_ids[i]) # fails when there is a single term
    newform <- .remove_fixef(termsobj, scope[i], nofixef=nofixef, keep_intercept = TRUE)
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
  termsobj <- terms(object)
  scope <- .preprocess_scope(scope, object=object, check=check, tl=attr(terms(termsobj), "term.labels"))
  if (is_long <- (fit_time <- how(object, verbose=FALSE)$fit_time)*length(scope)>check_time ) {
    message(paste0("Fitting the original model took ",fit_time,"s and drop1() may take a few times longer."))
  }
  progbar <- (is_long && length(scope>2L))
  basicLRTs <- vector("list", length(scope)) # ___F I X M E___ other variants of LR test? bootstrap, etc
  names(basicLRTs) <- scope
  nofixef <- .remove_all_fixef(object$predictor, keep_offset=TRUE)
  if (progbar) cat("\nProgress: ")
  for (i in seq_along(scope)) {
    if (progbar) cat(".")
    # newform <- buildmer::remove.terms(object$predictor, scope[i],check=check) # quite different approach:
    newform <- .remove_fixef(termsobj, scope[i], nofixef=nofixef, keep_intercept = TRUE)
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
          return(.drop1.lm(object, scope, ...))
        } else return(.drop1.glm(object, scope, ...))
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
