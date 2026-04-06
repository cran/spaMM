transffit <- function(object, lower=-2, upper=2, verbose=FALSE,
                      extracts= quote({y <- object$y}),
                      updates= quote({
                        if (lambda==0) {
                          newy <- log(y)
                        } else newy <- (y^lambda-1)/lambda
                        refit <- update_resp(object, newresp=newy)
                      }),
                      logDetJac= quote({(lambda-1)*sum(log(y)) })
) {
  if ( ! missing(updates) && missing(logDetJac)) 
    message("'updates' argument specified: maybe 'logDetJac' needed too?")
  eval(extracts)
  refit <- NULL
  objfn <- function(lambda, return_fit=FALSE) {
    eval(updates)
    if (return_fit) {
      return(refit)
    } else {
      logdetjac <- eval(logDetJac)
      resu <- logLik(refit)+logdetjac
      if (verbose) print(resu)
      return( - resu)
    }
  }
  if (length(lower)==1L) {
    optr <- optimize(f=objfn, interval=c(lower,upper), return_fit=FALSE)
    solution <- optr$minimum
  } else {
    optr <- .safe_opt(init=(lower+upper)/2, objfn=objfn,lower=lower, upper=upper,
                      LowUp=list(lower=lower, upper=upper), verbose=verbose, return_fit=FALSE)
    solution <- optr$solution
  }
  return(list(optr=optr,
              fit=objfn(lambda = solution, return_fit = TRUE)))
}
