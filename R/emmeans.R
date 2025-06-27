# don't need to be formally exported as methods 

recover_data.HLfit <- function (object, frame = model.frame(object), ...) {
  fcall <- getCall(object)
  pwts <- object$prior.weights
  if (attr(pwts,"is_unit")) pwts <- NULL
  recover_data_call <- get("recover_data.call", asNamespace("emmeans"), inherits=FALSE) 
  # should fail in standard way if not installed, so no wrapping by .get_wrap("recover_data.call", pack="emmeans")
  recover_data_call(fcall, delete.response(terms(object)), object$na.action,
                    frame = frame, pwts = pwts, ...)
}

# Conceived to replicate a subset of what glmmTMB:::emm_basis.glmmTMB() does
# So the code of the latter function may serve as 'doc' 
# for the design decisions of the present one.
emm_basis.HLfit <- local({
  warned_lm_compar <- FALSE
  warned_experimental  <- FALSE
  function(object, trms, xlev, grid, ...) {
    if ( ! warned_experimental) {
      warning("emmeans() method for HLfit objects is experimental.")
      warned_experimental  <<- TRUE
    }
    m <- model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    contrasts.arg <- attr(model.matrix(object),"contrasts")
    X <- model.matrix(trms, m, contrasts.arg = contrasts.arg) 
    bhat <- fixef(object) 
    Xmat <- model.matrix(trms, data=object$data)     
    if ( ! warned_lm_compar && 
         ! .REMLmess(object,return_message=FALSE) &&
         length(object$rand.families)==0L && # not mixed
         object$family$family=="gaussian" && ## on a mv-GLM, evaluates to FALSE
         .DEPARSE(.get_phiform(object))=="~1" # not heteroscedastic
    ) {
      warning(paste0("The model should have been fitted by REML if you want results consistent\n", 
                     "with standard F-test theory (and with results derived from 'lm' fit objects)."))
      warned_lm_compar <<- TRUE
    }
    ## Use of .get_new_X_ZAC_blob was suggested by my interpretation of an emmeans vignette.
    ## (cf ref to parallels with predict code)
    ## But glmmTMB's approach is much simpler...
    # new_X_ZACblob <- .get_new_X_ZAC_blob(object, newdata=NULL, re.form=NULL, variances=list(), 
    #                                      #invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list, control=control,
    #                                      na.action=na.pass)
    # 
    # new_X_ZACblob <- .get_new_X_ZAC_blob(object, newdata=NULL, re.form=NULL, variances=list(), 
    #                                      #invCov_oldLv_oldLv_list=invCov_oldLv_oldLv_list, control=control,
    #                                      na.action=na.pass)
    # 
    V <-  vcov(object)
    class(V) <- c("matrix", "array") # it already inherits from these, but emmeans does not handle the additional class
    nbasis <-  matrix(NA) 
    #
    if ( 
      (inherits(object,"fitmv") && all(sapply(object$families, `[[`, x="family")=="gaussian")) ||
      object$family$family == "gaussian"
    ) {
      dfargs <- list(df = df.residual(object))
      dffun <- function(k, dfargs) dfargs$df
    } else {
      dffun <- function(k, dfargs) Inf
      dfargs <- list()
    }
    #
    list(X = X, bhat = bhat, nbasis = nbasis, V = V,                
         dffun = dffun, dfargs = dfargs)
  }
})
