.get_bobyqa_controls <- function(init, upper, lower, maxeval_corr=.spaMM.data$options$maxeval_corr) {
  bobyqa_controls <- .spaMM.data$options$bobyqa
  if (is.null(bobyqa_controls$npt)) bobyqa_controls$npt <- 2*length(init)+1L
  if (is.null(bobyqa_controls$rhobeg)) bobyqa_controls$rhobeg <- .spaMM.data$options$bobyqa_rhofn(lower,upper)
  if (is.null(bobyqa_controls$rhoend)) bobyqa_controls$rhoend <- bobyqa_controls$rhobeg*1e-8 # bobyqa's default is rhobeg*1e-6 which is unsafe in the test_rC_transf sph case
  if (is.null(bobyqa_controls$maxfun)) {
    bobyqa_controls$maxfun <- max(2*( eval(.spaMM.data$options$maxeval,envir=list(initvec=init)))*maxeval_corr,
                                  1+10*length(init)^2) # bobyqa will complain if not > second value
  }
  bobyqa_controls
}

.get_nloptr_controls <- function(init, LowUp, maxeval_corr=.spaMM.data$options$maxeval_corr) {
  nloptr_controls <- .spaMM.data$options$nloptr
  if (is.null(nloptr_controls$maxeval)) nloptr_controls$maxeval <-  eval(.spaMM.data$options$maxeval,list(initvec=init))*maxeval_corr
  if (is.null(nloptr_controls$xtol_abs)) nloptr_controls$xtol_abs <- eval(.spaMM.data$options$xtol_abs, 
                                                                          list(LowUp=LowUp, rC_transf=.spaMM.data$options$rC_transf))
  if (is.null(nloptr_controls$xtol_abs)) nloptr_controls$xtol_abs <- 1e-12
  #nloptr_controls$print_level <- 3L # can be controlled by spaMM.options()!
  nloptr_controls
}

.safe_opt <- function(init, objfn, lower, upper, verbose, maxeval_corr=.spaMM.data$options$maxeval_corr, 
                      recheck_at_bound=.spaMM.data$options$recheck_at_bound, 
                      adjust_init=list(), 
                      LowUp,  
                      ...) { # minimization
  names_init <- names(init) # may be lost in later operations
  prevmin <- Inf
  delayedAssign("bobyqa_controls", .get_bobyqa_controls(init, upper, lower, maxeval_corr))
  delayedAssign("nloptr_controls", .get_nloptr_controls(init, LowUp, maxeval_corr))
  use_bobyqa <- (min(c(init-lower,upper-init))<1e-4)
  while (TRUE) {
    if (use_bobyqa) {
      if (verbose) cat("bobyqa: ")
      # This only bc bobyqa is more sensitive to the 14th decimal than nloptr at the boundaries  
      bobyqa_margin <- .spaMM.data$options$bobyqa_margin
      margin <- (upper-lower)*bobyqa_margin # test_rC_transf (sph) was a test of (this together with rhoend) but in that case adjust_init is a better fix
      margin <- pmin(bobyqa_margin,margin) # handles infinite ranges (but not only)
      init <- pmax(lower+margin,pmin(upper-margin,init))
      # And this is a more substantial adjustment at the margin
      if ( ! is.null(adjust_init$lower)) init <- pmax(adjust_init$lower, init)
      if ( ! is.null(adjust_init$upper)) init <- pmin(adjust_init$upper, init)
      #
      optr <- minqa::bobyqa(par=init, fn=objfn, lower=lower, upper=upper, control=bobyqa_controls, ...) ## does not use gradients
      if (optr$fval > prevmin-1e-8) { # i.e. progress is at most 1e-8
        if (optr$fval > prevmin) optr <- prev_optr # bobyqa may return much worse than initial value! if the objfn diverges at the bounds
        if (inherits(optr,"bobyqa")) { # may be FALSE if prev_optr was brought back.
          optr$solution <- optr$par ## let's stick to the nloptr() return elements (solution, objective)
          optr$objective <- optr$fval ## let's stick to the nloptr() return elements (solution, objective)
        }
        break ## no measurable minimization progress in <= maxeval iterations (indeed, stop if small regression)
      }
      if (optr$ierr==1L) { #maximum number of function evaluations exceeded
        init <- optr$par
        prev_optr <- optr
        prevmin <- prev_optr$fval
        use_bobyqa <- (min(c(init-lower,upper-init))<1e-4)
        #if (verbose) message("maxeval reached in bobyqa(); optimizer called again until apparent convergence of objective.") 
      } else { # conversion of bobyqa output to return format
        optr$solution <- optr$par ## let's stick to the nloptr() return elements (solution, objective)
        optr$objective <- optr$fval ## let's stick to the optim() return elements (par, value)
        break 
      }
    } else {
      if (verbose) cat("nloptr: ")
      optr <- nloptr::nloptr(x0=init,eval_f=objfn,lb=lower,ub=upper,
                             opts=nloptr_controls, ...)
      loc_ftol <- max(1e-8, optr$options$ftol_abs)
      if (optr$objective > prevmin-loc_ftol) {
        break ## no measurable minimization progress in <= maxeval iterations (indeed, stop if small regression)
      }
      init <- optr$solution
      at_bound <- (min(c(init-lower,upper-init))<1e-4)
      if (( recheck_at_bound && at_bound) || optr$status==5L) { ## 5 => termination bc maxeval has been reached
        # loop again
        prev_optr <- optr
        prevmin <- prev_optr$objective
        use_bobyqa <- at_bound
        #if (at_bound && verbose) message("maxeval reached in nloptr(); optimizer called again until apparent convergence of objective.") 
        #if (verbose && optr$status==5L) message("nloptr() went close to bound; calling bobyqa...") 
      } else break 
    }
    if (verbose>1L) print(c(objective=prevmin,next_init=init))
  }
  names(optr$solution) <- names_init
  return(optr) # use nloptr format (solution, objective) for return
}