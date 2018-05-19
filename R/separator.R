.do_call_wrap <- local(
  {
    warned_dcw <- FALSE
    function(chr_fnname,arglist, pack="e1071") {
      if (length(grep(pack,packageDescription("spaMM")$Imports))) {
        ## then the necessary functions must be imported-from in the NAMESPACE  
        do.call(chr_fnname,arglist) 
      } else if (length(grep(pack,packageDescription("spaMM")$Suggests))) {
        ## then the necessary functions cannot be imported-from in the NAMESPACE  (and the package must be written in an appropriate way)
        if ( requireNamespace(pack, quietly = TRUE)) {
          myfun <- get(chr_fnname, asNamespace(pack)) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
          do.call(myfun,arglist) 
        } else {
          if ( ! warned_dcw) {
            message("If the'",pack,"'package were installed, spaMM could check separation in binary regression problem.")
            warned_dcw <<- TRUE
          }
          return(NULL)
        }
      } else { ## package not declared in DESCRIPTION
        if (do.call("require",list(package=pack, quietly = TRUE))) {
          do.call(chr_fnname,arglist) 
        } else {
          if ( ! warned_dcw) {
            message("If the'",pack,"'package were installed, spaMM could check separation in binary regression problem.")
            warned_dcw <<- TRUE
          }
          return(NULL)
        }
      }
    }
  }
)



.is_separated <- local(
  {
    warned_is <- FALSE
    function(x,y) {
      if (requireNamespace("lpSolveAPI",quietly=TRUE)) {
        ## test for and/or find the direction of separation
        ## x a design matrix and y a 0-1 binary response vector
        separation <- .separator(x, as.numeric(y), purpose = "test")$separation
        if(separation) {
          message("Separation exists among the sample points.\n\tThis model cannot be fit by maximum likelihood.")
          message("The following terms are causing separation among the sample points:")
          beta <- .separator(x, as.numeric(y), purpose = "find")$beta
          separating.terms <- dimnames(x)[[2]][abs(beta) > 1e-09]
          if(length(separating.terms)) {
            separating.terms <- paste(separating.terms, collapse = ", ")
            message(paste(separating.terms,"\n"))
            warning(paste("separation due to", separating.terms))
          }
        }
        return(separation)
      } else {
        if ( ! warned_is) {
          message("If the 'lpSolveAPI' package were installed, spaMM could properly check (quasi-)separation in binary regression problem.")
          warned_is <<- TRUE
        }
        varcols <- logical(ncol(x))
        for (it in seq_len(ncol(x))) varcols[it] <- (diff(range(x[,it])) > .Machine$double.eps ^ 0.5)
        if (any(varcols)) {
          if (length(y)> spaMM.getOption("separation_max")) {
            message(paste("Increase spaMM.options(separation_max=<.>) to at least",length(y),
                          "if you want to check separation (see 'help(separation)')."))
            return(FALSE)
          } else {
            if (inherits(x,"sparseMatrix")) x <- as.matrix(x) ## bc next line ->  model.frame.default...
            arglist <- list(formula= y~x[,varcols,drop=FALSE],type='C-classification', kernel='linear')
            svmfit <- .do_call_wrap("svm",arglist=arglist)
            if ( ! is.null(svmfit)) {
              zut <- cbind(y=y,fv=svmfit$fitted)
              #return( ! any(svmfit$fitted!=y)) ## wrong ,this tested for perfect prediction of all response levels
              ## Better but (still) fails to detect quasi-separation:
              return(any(by(zut,zut[,1],function(v) length(unique(v[,2])))==1L)) ## this tests for perfect prediction of some response level
            } else return(FALSE)
          }
        } else return(FALSE)
      }
    }
  }
)


## assuming  Imports lpSolveAPI (>= 5.5.0.14)
## code derived from the glm() function in the safeBinaryRegression package
## does not import all package to prevent interference between glm fns and different error handling 

##FR See also Zorn, A Solution to Separation in Binary Response Models, Political Analysis (2005) 13:157-170
##FR see also http://www.ats.ucla.edu/stat/mult_pkg/faq/general/complete_separation_logit_models.htm
##FR there is a Firth method for dealing with this (and a package, brglm... brglm(locform,offset=Offset) ).
.separator <- function(x, y, method = c("primal", "dual"), purpose = c("test", "find"),
                      tol = 1e-3)
{
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  dimnames(x) <- NULL
  
  y.bar <- -sign(y - 0.5)
  x.bar <- y.bar * x
  
  ans <- list()
  
  if(method == "primal" && purpose == "test") {
    lp <- lpSolveAPI::make.lp(n, p)
    for(j in 1:p) status <- lpSolveAPI::set.column(lp, j, x.bar[, j])
    status <- lpSolveAPI::set.rhs(lp, rep(0.0, n))
    status <- lpSolveAPI::set.constr.type(lp, rep(1, n))
    status <- lpSolveAPI::set.objfn(lp, -colSums(x.bar))
    status <- lpSolveAPI::set.bounds(lp, lower = rep(-Inf, p), upper = rep(Inf, p))
    control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "max",
                          simplextype = c("primal", "primal"))
    status <- lpSolveAPI::solve.lpExtPtr(lp)
    
    if(status == 0)
      ans$separation <- FALSE
    else if(status == 3)
      ans$separation <- TRUE
    else {
      stop("unexpected result for primal test.")
    }
  }
  
  if(method == "primal" && purpose == "find") {
    lp <- lpSolveAPI::make.lp(n, p)
    for(j in 1:p) status <- lpSolveAPI::set.column(lp, j, x.bar[, j])
    status <- lpSolveAPI::set.rhs(lp, rep(0.0, n))
    status <- lpSolveAPI::set.constr.type(lp, rep(1, n))
    status <- lpSolveAPI::set.objfn(lp, -colSums(x.bar))
    status <- lpSolveAPI::set.bounds(lp, lower = rep(-1, p), upper = rep(1, p))
    control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "max",
                          simplextype = c("primal", "primal"))
    status <- lpSolveAPI::solve.lpExtPtr(lp)
    
    if(status != 0) {
      stop("unexpected result for primal test.")
    }
    beta <- lpSolveAPI::get.variables(lp)
    
    if(sum(abs(beta)) > tol)
      ans$separation <- TRUE
    else
      ans$separation <- FALSE
    
    ans$beta <- beta
  }
  
  if(method == "dual" && purpose == "test") {
    lp <- lpSolveAPI::make.lp(p, n)
    for(j in 1:n)  status <- lpSolveAPI::set.column(lp, j, x.bar[j, ])
    status <- lpSolveAPI::set.rhs(lp, -colSums(x.bar))
    status <- lpSolveAPI::set.constr.type(lp, rep(3, p))
    status <- lpSolveAPI::set.objfn(lp, rep(0.0, n))
    status <- lpSolveAPI::set.bounds(lp, lower = rep(0.0, n), upper = rep(Inf, n))
    control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "min",
                          simplextype = c("primal", "primal"))
    status <- lpSolveAPI::solve.lpExtPtr(lp)
    
    if(status == 0)
      ans$separation <- FALSE
    else if(status == 2)
      ans$separation <- TRUE
    else {
      stop("unexpected result for dual test.")
    }
  }
  
  if(method == "dual" && purpose == "find") {
    lp <- lpSolveAPI::make.lp(p, n + 2*p)
    for(j in 1:n)
      status <- lpSolveAPI::set.column(lp, j, x.bar[j, ])
    for(j in 1:p)
      status <- lpSolveAPI::set.column(lp, n+j, -1.0, j)
    for(j in 1:n)
      status <- lpSolveAPI::set.column(lp, n+p+j, 1.0, j)
    b <- -colSums(x.bar)
    status <- lpSolveAPI::set.rhs(lp, b)
    status <- lpSolveAPI::set.constr.type(lp, rep(3, p))
    status <- lpSolveAPI::set.objfn(lp, rep(c(0.0, 1.0), c(n, 2*p)))
    status <- lpSolveAPI::set.bounds(lp, lower = rep(0.0, n + 2*p), upper = rep(Inf, n + 2*p))
    control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "min",
                          simplextype = c("primal", "primal"))
    basis <- 1:p
    basis[b >= 0.0] <- basis[b >= 0.0] + p
    status <- lpSolveAPI::set.basis(lp, -(n + p + basis))
    status <- lpSolveAPI::solve.lpExtPtr(lp)
    
    beta <- lpSolveAPI::get.dual.solution(lp)[2:(p+1)]
    
    if(sum(abs(beta)) > tol)
      ans$separation <- TRUE
    else
      ans$separation <- FALSE
    
    ans$beta <- beta
  }
  
  ans
}

