.calc_CARdispGammaGLM <- function(data,family=spaMM_Gamma("inverse"),lev=NULL,
                                lambda_rd, control) {
  if ( ! is.null(lev)) {
    resp <- data$resp/(1-lev) 
    resp[resp==0L] <- 1e-100
    resp[resp>1e150] <- 1e150 ## v1.2 fix for glm -> glm.fit -> .Call(C_Cdqrls, x[good, , drop = FALSE]*w problem
    data$resp <- resp
    weights <- (1-lev)/2 
    loclist <- list(data=data, family=family, weights=weights)
  } else loclist <- list(data=data, family=family)
  if (is.na(lambda_rd)) {
    loclist$formula <- resp~1+adjd
    #loclist$start <- c(1/mean(data$resp),0)
  } else {
    loclist$formula <- resp ~ adjd -1
    loclist$offset <- rep(1/lambda_rd,nrow(data))
    loclist$start <- c(0)
  }
  loclist$control <- control
  suppressWarnings(do.call("spaMM_glm",loclist))
}

# a specialized version of glm, which further constructs the response value internally; the formula has no lhs 
.calc_dispGammaGLM <- function (formula, 
            dev.res, ## info not in the 'data' argument;
            # In the case where there was a prior_lam_fac is the 'design' for non-ranCoef (wei-1|.), the unique_lambda was obtained by
            # computing the dev.res as fn of such weights, rather than through a non-unit design X. The present code will replicate this, as it will use 
            # a using X, and the same dev.res passed through resp_lambda to the present 'dev.res'. Cf the block with
            # safe_dev.res_it <- rand.families[[it]]$dev.resids(u_h[u.range],psi_M[u.range],wt=wt) ## must give d1 in table p 989 de LeeN01
            # ...      
            # resp_lambda[u.range] <- safe_dev.res_it                                                       (all correct but quite contrived!)
            lev,
            data, ## provide only the predictor variables ! alternatively could proceed as calc_CARdispGammaGLM
            family = spaMM_Gamma(link=log), #inverse link is used in .calc_CARdispGammaGLM()
            offset=rep.int(0, NROW(dev.res)),
            na.action, start = NULL, etastart, mustart, 
            control, 
            #try = FALSE, 
            ## more args from glm() def:
            subset, ##  not used but rethink.
            model = TRUE, ## whether to return the mf, cf end of code
            x = FALSE, ## whether to return the design matrix, cf end of code
            y = TRUE, ## whether to include y... idem 
            contrasts = NULL,
            method="glm.fit",
            ...) {
  Y <- as.numeric(dev.res/((1-lev))) 
  ## potential problem of any corrections is that v_h and lambda estimates may be inconsistent 
  Y[Y==0L] <- 1e-100 ## v_h ~0 occurs eg Poisson, y ~ 1+ ranef(lambda->very small)  ... and COMPoisson
  Y[Y>1e150] <- 1e150 ## v1.2 fix for glm -> glm.fit -> .Call(C_Cdqrls, x[good, , drop = FALSE]*w problem
  Y <- exp(.sanitize_eta_log_link(log(Y), max=30,y=Y)) # previous two lines important for the y=Y argument, which is itself important.
  #
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset","X"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame) 
  if (missing(data)) data <- environment(formula)
  mf <- eval(mf, parent.frame()) ## formula argument is required to eval mf => for mt too
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf, contrasts) ## from glm() code
  #
  offset <- as.vector(model.offset(mf))
  intercept <- FALSE ## previously tried attr(mt, "intercept") > 0L but this is wrong
  #
  weights <- (1-lev)/2 # glm() code : weights <- as.vector(stats::model.weights(mf))
  fit <- NULL
  if ( ! is.null(start)) {
    tryfit <- .tryCatch_W_E(eval(call(method, 
                                      x = X, y = Y, weights = weights, offset = offset, family = family, 
                                      start=start, etastart=etastart, # two mandatory arguments for method= .glm_reformat(); 
                                      control = control, intercept = intercept))) 
    fit <- tryfit$value
    warnmess <- fit$warning$message
    #if (inherits(fit,"error") || is.na(fit$null.deviance)) { ## is.na(...): dgamma problem but returns a fit
    if (inherits(fit,"error")) fit <- NULL # spaMM_glm.fit could be directly called as it should work in all cases (and 'method'=.glm_reformat) has long been called 
    # with incomplete arguments so it failed systematically, so that spaMM_glm.fit) was called systematically).
  }
  if (is.null(fit)) { 
    # Gamma() with small response values... 
    ## we need a valid GLM: we fit again with a more controlled method
    control$maxit <- 1000L ## convergence is very slow when fitting y ~ 1 where most y values are <1e-12 and one is ~1e-10. 
    # Further, if all y are very small, a strict control may give spurious warnings:
    if (all(Y<1e-12)) control$epsilon <- 1e-7 # rather than 1e-8
    fit <- spaMM_glm.fit(x = X, y = Y, weights = weights, offset = offset, family = family, 
                     control = control, intercept = intercept)
    warnmess <- NULL
  }
  if (model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x) fit$x <- X
  if (!y) fit$y <- NULL
  fit <- c(fit, list(call = match.call(), formula = formula, terms = mt, 
                     data = data, offset = offset, control = control, method= "glm.fit",
                     contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf),
                     warnmess=warnmess
                     ))
  class(fit) <- c(fit$class, c("glm", "lm"))
  fit
}