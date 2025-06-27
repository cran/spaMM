# quick base R reformatting
reshape2long  <- function(mndata, types) {
  template <- rep(0L,length(types))
  names(template) <- types
  vvv <- lapply(seq_len(nrow(mndata)),
                function(i) {
                  NA_template <- template
                  v  <- mndata[i,]
                  NA_template[is.na(v)] <- NA_integer_ 
                  vv <- lapply(types, function(typ) {
                    if (is.na(nr <- v[,typ])) {
                      NULL
                    } else {
                      NA_template[typ]  <- 1L
                      v[,types] <- NA_template[types]
                      v[rep(1,nr),]
                    }
                  })
                  do.call(rbind,vv)
                }
  )
  do.call(rbind,vvv)
}

.get_validrownames <- function(mc, data) {
  mc[["data"]] <- data
  mc[["verbose"]]["getCall"] <- TRUE
  mc[[1L]] <-  get("fitmv", asNamespace("spaMM"), inherits=FALSE)  
  mvp <- eval(mc,parent.frame()) # first fit# using input 'init' if any
  attr(mvp$processed$data,"validrownames") 
}

# experimental... and not used in effective code so far.
.get_processed_call <- function(mc, ...) {
  dotlist <- list(...)
  for (st in names(dotlist)) 
  mc[[st]] <- dotlist[[st]]
  mc[["verbose"]]["getCall"] <- TRUE # This avoids the fit (since we want $processed)
  mc <- eval(mc,parent.frame()) 
  #
  # if mc[[1]] was fitme/fitmv then ..._body() returned a call with $processed, 
  # and the mc[[1]] function attached a $call to it.
  # Then if we eval() a new call built from this return value
  # the new call may have an argument $call that may be evaluated when the dotlist is analyzed
  # This should be avoided! (argh). So
  mc[["call"]] <- NULL
  #
  mc[["verbose"]]["getCall"] <- FALSE
  mc[["processed"]]$"verbose"[["getCall"]] <- FALSE
  mc
}

.get_muP_template <- function(data, types, validrownames) {
  dim_data <- dim(data[,types])
  muP_template <- rep(NA_real_, prod(dim_data))
  dim(muP_template) <- dim_data
  rownames(muP_template) <- rownames(data)
  for (it in seq_along(validrownames)) muP_template[validrownames[[it]],it] <- 0
  muP_template
}

pois4mlogit <- function(submodels, data, to.long=FALSE,
                        init=list(), control=list(), # names(init) possibly used in all iterations...
                        ..., next_inits=c("ranPars","v_h","fixef"), 
                        types, n_iter=1000L, tol=1e-5, fac=1, progress=FALSE) {
  update_fitmv_body <- control$update_fitmv_body
  if (is.null(update_fitmv_body)) update_fitmv_body <- FALSE
  fac_is_1  <- abs(fac-1) < 100*.Machine$double.eps
  time1 <- Sys.time()
  mc <- match.call(expand.dots = TRUE) 
  mc["to.long"] <- NULL
  mc["types"] <- NULL
  mc["n_iter"] <- NULL
  mc["tol"] <- NULL
  mc["progress"] <- NULL
  mc["next_inits"] <- NULL
  mc["fac"] <- NULL
  # if (missing(types)) { # formula is not equiv to submodel 
  #   types <- sapply(submodels, function(form) deparse(form[[1]][[2]]))
  # }
  # Need a first .dynoffset before preprocessing:  
  null_user_dynoffset <- is.null(data$.dynoffset)
  if (null_user_dynoffset) data$.dynoffset <- 0 # only for .get_validrownames() call
  # Build template matrix with NAs for missing response or predictors per submodel.
  validrownames <- .get_validrownames(mc, data) # preprocessing to diagnose the data. (nested fitmv call with full preprocessing)
  muP_template <- .get_muP_template(data, types, validrownames)
  data <- data[rowSums( ! is.na(muP_template)) > 1L,] # more than one class must have information
  if (to.long) {
    # if length(setdiff(types, colnames(data))) stop("responses are not variables in the data.frame: provide a 'types' argument")
    data <- reshape2long(mndata=data, types = types)
    log_mnsizes <- 0
    mnsizes <- 1L
  } else {
    # mnsizes <- eval(expr=str2lang(paste(types, collapse="+")), data) # not requiring the same test as in to.long case
    mnsizes <- rowSums(data[,types],na.rm = TRUE)
    if (any(mnsizes==0L)) {
      data <- data[mnsizes!=0L,]
      mnsizes <- rowSums(data[,types],na.rm = TRUE)
    }  
    log_mnsizes <- log(mnsizes)
  }
  if (null_user_dynoffset) data$.dynoffset <- log_mnsizes
  
  mc[["data"]] <- data
  mc[[1L]] <- get("fitmv", asNamespace("spaMM"), inherits=FALSE)  
  if (update_fitmv_body) {
    mc <- .get_processed_call(mc=mc, data=data) 
    processed <- mc[["processed"]] 
    data  <- processed$data # get the processed data
    mc[[1L]] <- get("fitmv_body", asNamespace("spaMM"), inherits=FALSE)  
    names(mc)[names(mc)=="fixed"] <- "fixedS" 
  } 
  mvp <- eval(mc,parent.frame()) # first fit# using input 'init' if any
  
  # Build new template matrix with NAs for missing response values (relevant for dynoffset computation)
  validrownames <- attr(mvp$data,"validrownames") 
  muP_template <- .get_muP_template(data, types, validrownames)
  #
  muP_template[ ! is.na(muP_template)] <- mvp$eta
  muP_template <- muP_template - data$.dynoffset
  muP <- exp(muP_template) # 
  output.dynoffset <- log_mnsizes-log(rowSums(muP,na.rm = TRUE))
  names.init <- names(init)
  prevmsglength <- 0L
  oldScrit <- oldOcrit <- Inf
  oldlogL <- -Inf
  logL <- logLik(mvp)
  d_off_fac <- 1
  notwarned <- progress>=0L
  successful.input.dynoffset  <- data$".dynoffset"
  for (it in seq_len(n_iter)) {
    data$".dynoffset" <- d_off_fac*output.dynoffset +(1-d_off_fac)*successful.input.dynoffset# \sum^J
    if ("ranPars" %in% next_inits) {
      newinits <- get_inits_from_fit(mvp, to_fn="fitmv_body")[["init"]]
    } else  newinits <- get_inits_from_fit(mvp, to_fn="fitmv_body")[["init"]][names.init]
    if (update_fitmv_body) { # direct update on fitmv_body call
      processed$data  <- data
      # The processed offset is stored in processed$off, which has to be recomputed
      # (and its mv-length is distinct from that of data$".dynoffset" )
      processed$off  <- model.offset.HLfit(mvp, data=data)
      mc[["processed"]]  <- processed
      # Initialize next outer optim: 
      mc["init"] <- list(newinits) 
      # The current port_env (with scaled values) will be used to initialize the next "inner fit",
      # provided we remove any original init.HLfit ...
      mc["init.HLfit"] <- list(NULL) 
      # ... but we also have to hack the $port_env$objective, otherwise the (irrelevant) final logL 
      # of the previous fit with different .dynoffset would control further updating:
      processed$port_env$objective  <- -Inf
      mvp <- eval(mc,parent.frame())
    } else {
      init.HLfit <- list()
      # $processed is recreated in each iter so we use init.HLfit for the first "inner fit" of the call
      if ("v_h" %in% next_inits) init.HLfit$v_h <- unlist(ranef(mvp))
      if ("fixef" %in% next_inits) init.HLfit$fixef <- fixef(mvp)
      mvp <- update(mvp, data=data, init.HLfit=init.HLfit, init=newinits) 
    }
       
    #
    # Scrit computed before correction of muP_template by .dynoffset
    muP_template[ ! is.na(muP_template)] <- mvp$eta
    Scrit <- rowSums(exp(muP_template), na.rm=TRUE)
    Scrit <- mean(abs(Scrit-mnsizes))
    #
    # Ocrit use new offsets computed from (old-offset)-included muP's
    muP_template <- muP_template - data$.dynoffset
    muP <- exp(muP_template) # 
    output.dynoffset <- log_mnsizes-log(rowSums(muP,na.rm = TRUE))
    Ocrit <- mean(abs(output.dynoffset-data$.dynoffset))
    #
    
    logL <- logLik(mvp)
    dlogL <- logL-oldlogL
    if (dlogL > 0 || fac_is_1) {
      if ( ! fac_is_1) d_off_fac <- d_off_fac*fac
      successful.input.dynoffset <- data$".dynoffset"
      oldlogL <- logL
    } else { # no progress
      if (d_off_fac > 1) { # problem due to too large 'd_off_fac'
        logL <- oldlogL
      } else { # possible problem, but it's still best to accept step
        successful.input.dynoffset <- data$".dynoffset"
        oldlogL <- logL
      }
      d_off_fac <- max(1,d_off_fac*1/3) # do not allow d_off_fac to vanish 
    }
    
    if (progress>1L) prevmsglength <- overcat(paste0(it,": logL:",signif(logL,5),
                                                     " Ocrit: ",signif(Ocrit,3L),
                                                     " Scrit: ",signif(Scrit,3L),"         "), 
                                              prevmsglength)
    if (notwarned && abs(dlogL) < tol) {
      # if a>0 in (Scrit=a+b lambda_b^t) Scrit_lam_b is O(lambda_b^t), t Scrit_lam_b should still vanish
      # Ideally, a=0, Scrit_lam_b is O(lambda_b), t Scrit_lam_b will diverge 
      Scrit_lam_b  <- 1- Scrit/oldScrit # ideally large. So we test if Scrit_lam_b is small, as this is suspect.
      Ocrit_lam_b  <- 1- Ocrit/oldOcrit # same.
      # Rkably one of the identifiable 'pollen' fits shows these two lam_b's 
      #   staying large and ~constant for some time. Only the logL improves.
      if ( abs(Ocrit_lam_b)<tol && abs(Scrit_lam_b)<tol) { 
        warning("Something suspect... maybe unidentifiable model,\n    e.g. with an intercept shared among submodels?",
                immediate. = TRUE)
        notwarned <- FALSE
      }
    }
    oldScrit <- Scrit
    oldOcrit <- Ocrit
    # oldlogL <- logL
    cond <- Ocrit<tol && Scrit<tol
    if (cond) break
  } 
  if (progress>=0L) {
    if ( ! cond) {
      warning(
        paste("pois4mlogit() fit did not converge in",n_iter,
              "iterations (Ocrit: ",signif(Ocrit,3L),
              ", Scrit: ",signif(Scrit,3L),")"),
        immediate. = TRUE)
    } else if (progress>1L) { # case with overcat's
      cat("\n")
    } else print(paste("Fit converged in", it,"iterations."), quote=FALSE)
  }
  # mvp$APHLs$mnConst <- lgamma(nrow(data)+1L)
  if ( ! inherits(mvp,"HLfitlist") && ! is.call(mvp) ) {
    mvp$call[[1]] <- quote(spaMM::pois4mlogit)
    mvp$mnsizes <- mnsizes
    mvp$info <- list(Ocrit=Ocrit, Scrit=Scrit, it=it)
    mvp$how$fnname <- "pois4mlogit"
    fit_time <- .timerraw(time1) 
    mvp$how$fit_time <- structure(fit_time,
                              message="Please use how(<fit object>)[['fit_time']] to extract this information cleanly.")
  }
  class(mvp) <- c("pois4mlogit",class(mvp))
  mvp 
}

predict.pois4mlogit <- function(object, newdata=NULL, ...) {
  .predict.pois4mlogit(object=object, newdata=newdata, ...) # this hides the private 'check' argument
}

.predict.pois4mlogit.check <- function(object, mnsizes, ...) {
  ndata <- object$data
  ndata$".dynoffset" <- ndata$".dynoffset" - log(mnsizes) # quick correction without predict.HLfit call.
  # This code was first conceived to avoid the double predict call in .predict.pois4mlogit(), 
  # but may be confusing because the predictions then 'exactly' sum to 1 only if convergence was 'exact'.
  predict.HLfit(object, newdata=ndata, ...)
}

.predict.pois4mlogit <- function(object, newdata=NULL, ..., check=FALSE) {
  if (check) {
    # This code may be useful to provide some measure of convergence inaccuracy on prediction.
    if ( ! is.null(newdata)) stop("check=TRUE does not handle 'newdata'." ) 
    if ( is.null(mnsizes <- object$mnsizes)) stop("Old object, $mnsizes missing")
    .predict.pois4mlogit.check(object, mnsizes=mnsizes, ...)
  } else { # default case: double predict call, the first without the '...'
    # This double call makes sure that frequencies sum to 1, even if convergence was not 'exact'. 
    if (is.null(newdata)) newdata <- object$data
    newdata$".dynoffset" <- 0
    pred <- predict.HLfit(object, newdata=newdata) 
    pred <- matrix(pred, ncol=length(formula(object)))
    rowsums <- rowSums(pred)
    newdata$".dynoffset" <- - log(rowsums)
    predict.HLfit(object, newdata=newdata, ...)
  }
}

