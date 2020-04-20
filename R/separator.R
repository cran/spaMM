.do_call_wrap <- local(
  {
    warned_dcw <- list()
    function(chr_fnname,arglist, pack="e1071", info_mess) {
      if (length(grep(pack,packageDescription("spaMM")$Imports))) {
        ## then the necessary functions are imported-from in the NAMESPACE  
        do.call(chr_fnname,arglist) 
      } else if (length(grep(pack,packageDescription("spaMM")$Suggests))) {
        ## then the necessary functions are not imported-from in the NAMESPACE  (and the package must be written in an appropriate way)
        if ( requireNamespace(pack, quietly = TRUE)) {
          myfun <- get(chr_fnname, asNamespace(pack)) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
          do.call(myfun,arglist) 
        } else {
          if ( ! identical(warned_dcw[[pack]],TRUE)) {
            if (pack=="e1071") message("If the 'e1071' package were installed, spaMM could check separation in binary regression problem.")
            if (pack=="cubature") message("If the 'cubature' package were installed, spaMM could compute a requested marginal prediction.")
            warned_dcw[[pack]] <<- TRUE
          }
          return(NULL)
        }
      } else { ## package not declared in DESCRIPTION
        if (suppressWarnings(do.call("require",list(package=pack, quietly = TRUE)))) { ## 'quietly' only inhibits the startup message
          do.call(chr_fnname,arglist) 
        } else {
          if ( ! identical(warned_dcw[[pack]],TRUE)) {
            if (pack=="e1071") message("If the 'e1071' package were installed, spaMM could check separation in binary regression problem.")
            if (pack=="cubature") message("If the 'cubature' package were installed, spaMM could compute a requested marginal prediction.")
            if (pack=="pracma") message(info_mess)
            warned_dcw[[pack]] <<- TRUE
          }
          return(NULL)
        }
      }
    }
  }
)

## See also Zorn, A Solution to Separation in Binary Response Models, Political Analysis (2005) 13:157-170
## see also http://www.ats.ucla.edu/stat/mult_pkg/faq/general/complete_separation_logit_models.htm
## there is a Firth method for dealing with this (and a package, brglm... brglm(locform,offset=Offset) ).
.test_and_find_sep <- function(x, y, solver, tol = 1e-3, verbose=TRUE, ...) {
  n <- nrow(x)
  p <- ncol(x)
  
  y.bar <- -sign(y - 0.5)
  x.bar <- y.bar * x
  
  op <- OP(colSums(x.bar))
  constraints(op) <- L_constraint(L = x.bar, dir = leq(n), rhs = double(n))
  bounds(op) <- V_bound(lb = rep.int(-Inf, p), ub = rep.int(Inf, p))
  
  time1 <- Sys.time()
  s <- ROI_solve(op, solver = solver, ...)
  sep_time <- .timerraw(time1)
  if (sep_time>1) message(paste0("Re-checking separation for binomial-response model took ",sep_time," s."))
  separation <- as.logical(solution(s, "status_code"))
  if (separation) { # some problem
    if (solution(s, "msg")$status==6) { # unbounded: separation
      bounds(op) <- V_bound(lb = rep.int(-1, p), ub = rep.int(1, p))
      
      s <- ROI_solve(op, solver = solver, ...)
      beta <- solution(s, force = TRUE)
      separating.terms <- dimnames(x)[[2]][abs(beta) > tol]
      if(length(separating.terms)) {
        mess <- paste("The following terms are causing separation among the sample points:",
                      paste(separating.terms, collapse = ", "))
        if (verbose) message(paste(mess,
                                   "\n\tsome estimates of fixed-effect coefficients could be practically infinite,",
                                   "\n\tcausing numerical issues in various functions."))
        warning(mess)
      } 
    } else { # unidentified problem
      warning(paste("ROI solver returned",s$status$msg$message,"during check for separation."),immediate. = TRUE)
      separation <- FALSE
    }
  }
  return(separation)
}

is_separated <- local(
  {
    warned_is <- FALSE
    function(x,y, verbose=TRUE, solver=spaMM.getOption("sep_solver")) {
      pb_size <- (1e-5*prod(dim(x)))
      if (pb_size> spaMM.getOption("separation_max")) {
        message(paste("Increase spaMM.options(separation_max=<.>) to at least", ceiling(pb_size),
                      "if you want to check separation (see 'help(separation)')."))
        return(FALSE)
      } else separation <- .test_and_find_sep(x, y, verbose=verbose, solver=solver)
    }
  }
)

is_separated.formula <- function(formula, ..., separation_max=spaMM.getOption("separation_max"),
                                 solver=spaMM.getOption("sep_solver")) { 
  # is_separated.formula -> .preprocess(,For="is_separated") -> is_separated(X.pv, as.numeric(y))
  mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(mat_sqrt)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(fitme)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in is_separated() call."))
  # 
  FHF <- formals(HLfit) ## makes sure about default values 
  names_FHF <- names(FHF)
  names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
  FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
  preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(.preprocess)))] 
  preprocess.formal.args$For <- "is_separated"
  preprocess.formal.args$family <- binomial() 
  preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
  preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
  preprocess.formal.args$ranFix <- mc$fixed ## because preprocess expects ranFix # note that etaFix is handled in the FHF's
  preprocess.formal.args$adjMatrix <- mc$adjMatrix ## because adjMatrix not in formals(HLfit)
  preprocess.formal.args$corrMatrix <- mc$corrMatrix ## because corrMatrix not in formals(HLfit)    #
  preprocess.formal.args$covStruct <- mc$covStruct ## because covStruct not in formals(HLfit)    #
  preprocess.formal.args$method <- mc$method ## forces evaluation
  preprocess.formal.args$init <- mc$init ## because init not in formals(HLfit)    #
  #
  oldopt <- spaMM.options(separation_max=separation_max, sep_solver=solver)
  isSeparated <- do.call(.preprocess,preprocess.formal.args,envir=parent.frame(1L))
  spaMM.options(oldopt)
  return(isSeparated)
}




