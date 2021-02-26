.do_call_wrap <- local(
  {
    warned_dcw <- list()
    inla_already <- FALSE
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
        if (pack=="INLA") {
          success <- suppressMessages(do.call("require",list(package=pack))) # messge might suggest that the fit is by INLE...
          if ( ! inla_already) {
            message("INLA::inla.spde.make.A() will be used to construct the IMRF models fitted by spaMM.")
            inla_already <<- TRUE
          }
        } else success <- suppressWarnings(do.call("require",list(package=pack, quietly = TRUE)))
        if (success) { ## 'quietly' only inhibits 'Le chargement a necessite le package :'... 
          do.call(chr_fnname,arglist) 
        } else {
          if ( ! identical(warned_dcw[[pack]],TRUE)) {
            if (pack=="INLA") message("If the 'INLA' package were installed, spaMM could use INLA:::inla.spde.make.A().")
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
  if (solver=="svm") {
    varcols <- logical(ncol(x))
    for (it in seq_len(ncol(x))) varcols[it] <- (diff(range(x[,it])) > .Machine$double.eps ^ 0.5)
    if (any(varcols)) {
      time1 <- Sys.time()
      if (inherits(x,"sparseMatrix")) x <- as.matrix(x) ## bc next line ->  model.frame.default...
      arglist <- list(formula= y~x[,varcols,drop=FALSE],type='C-classification', kernel='linear')
      svmfit <- .do_call_wrap("svm",arglist=arglist)
      sep_time <- .timerraw(time1)
      if (sep_time>1) message(paste0("Checking separation for binomial-response model took ",sep_time," s."))
      if ( ! is.null(svmfit)) {
        zut <- cbind(y=y,fv=svmfit$fitted)
        #return( ! any(svmfit$fitted!=y)) ## wrong ,this tested for perfect prediction of all response levels
        ## Better but (still) fails to detect quasi-separation:
        separation <- (any(by(zut,zut[,1],function(v) length(unique(v[,2])))==1L)) ## this tests for perfect prediction of some response level
      } else separation <- FALSE
    } else separation <- FALSE
  } else { # ROI-based methods
    n <- nrow(x)
    p <- ncol(x)
    
    y.bar <- -sign(y - 0.5)
    x.bar <- y.bar * x
    
    op <- OP(colSums(x.bar))
    if (inherits(x.bar,"sparseMatrix")) x.bar <- as.matrix(x.bar) ## bc L_constraint does not handle sparse...
    constraints(op) <- L_constraint(L = x.bar, dir = leq(n), rhs = double(n))
    bounds(op) <- V_bound(lb = rep.int(-Inf, p), ub = rep.int(Inf, p))
    
    time1 <- Sys.time()
    s <- ROI_solve(op, solver = solver, ...)
    sep_time <- .timerraw(time1)
    if (sep_time>1) message(paste0("checking separation for binomial-response model took ",sep_time," s."))
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
  }
  return(separation)
}

is_separated <- local({
    warned_is <- FALSE
    function(x,y, verbose=TRUE, solver=spaMM.getOption("sep_solver")) {
      if (solver!="svm") {
        if (requireNamespace("ROI.plugin.glpk",quietly=TRUE)) {
          pb_size <- (1e-5*prod(dim(x)))
        } else {
          if ( ! warned_is) {
            message(paste0("If the 'ROI.plugin.glpk' package were installed,\n",
                           "spaMM could properly check (quasi-)separation in binary regression problem.\n",
                           "See help('external-libraries') if you have troubles installing 'ROI.plugin.glpk'."))
            warned_is <<- TRUE
          }
          pb_size <- length(y)/100
          solver <- "svm"
        }
      }
      if (pb_size> spaMM.getOption("separation_max")) {
        message(paste("Increase spaMM.options(separation_max=<.>) to at least", ceiling(pb_size),
                      "if you want to check separation (see 'help(separation)')."))
        return(FALSE)
      } else separation <- .test_and_find_sep(x, y, verbose=verbose, solver=solver)
    }
})

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
  preprocess_args <- .get_inits_preprocess_args(For="is_separated") 
  # the doc says about '...' "possibly other arguments of a fitme call" but many will be ignored and the following code cannot be as in fitme()
  names_nondefault  <- intersect(names(mc),names(preprocess_args)) ## mc including dotlist
  preprocess_args[names_nondefault] <- mc[names_nondefault] # Notably handling etaFix, one of the few relevant arguments for is_separated()  
  preprocess_args$predictor <- mc$formula ## because preprocess stll expects $predictor 
  #
  oldopt <- spaMM.options(separation_max=separation_max, sep_solver=solver)
  isSeparated <- do.call(.preprocess, preprocess_args, envir=parent.frame(1L))
  spaMM.options(oldopt)
  return(isSeparated)
}




