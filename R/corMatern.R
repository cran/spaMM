corMatern <- function(value = c(1, 0.5), form = ~1, nugget=FALSE, nuScaled=FALSE, 
                      metric = c("euclidean", "maximum", "manhattan"), fixed = FALSE) {
    attr(value, "formula") <- form
    attr(value, "nugget") <- nugget
    attr(value, "nuScaled") <- nuScaled
    attr(value, "metric") <- match.arg(metric)
    attr(value, "fixed") <- fixed
    class(value) <- c("corMatern", "corStruct") ## corSpatial removed
    value
}


getCovariate.corMatern <- function(object, form = formula(object), data) {
  if (is.null(covar <- attr(object, "covariate"))) { # need to calculate it
    if (missing(data)) {
      stop("need data to calculate covariate")
    }
    covForm <- nlme::getCovariateFormula(form)
    if (length(all.vars(covForm)) > 0) { # covariate present
      if (attr(terms(covForm), "intercept") == 1) {
        covForm <- eval(parse(text = paste("~", .DEPARSE(covForm[[2]]),"-1",sep="")))
      }
      covar <- as.data.frame(unclass(model.matrix(covForm,
                                                  model.frame(covForm, data, drop.unused.levels = TRUE))))
    } else {
      covar <- NULL
    }
    if (!is.null(nlme::getGroupsFormula(form))) { # by groups
      grps <- getGroups(object, data = data)
      if (is.null(covar)) {
        distmats <- lapply(split(grps, grps), function(x) (proxy::dist(1:length(x))))
      } else {
        locfn <- function(el, metric) {
          el <- as.matrix(el)
          if (nrow(el) > 1) {
            (proxy::dist(el, metric))
          } else {
            numeric(0)
          }
        }
        distmats <- lapply(split(covar, grps), locfn, metric = attr(object, "metric"))
      }
      covar <- lapply(distmats,as.vector)
      covar <- covar[sapply(covar, length) > 0]  # no 1-obs groups ## takes only list elements of nonzero length
    } else {				# no groups
      if (is.null(covar)) {
        distmats <- (proxy::dist(1:nrow(data)))
      } else {
        distmats <- (proxy::dist(as.matrix(covar),
                                 method = attr(object, "metric")))
      }
      covar <- lapply(distmats,as.vector)
    }
    if (any(unlist(covar) == 0)) {
      stop("cannot have zero distances in \"corMatern\"")
    }
    attr(covar,"distmats") <- distmats
  }
  covar
}

Initialize.corMatern <- function(object, data, ...) {
  if (!is.null(attr(object, "minD"))) { #already initialized
    return(object)
  }
  
  form <- formula(object)
  ## obtaining the groups information, if any
  if (!is.null(nlme::getGroupsFormula(form))) {
    attr(object, "groups") <- getGroups(object, form, data = data)
    attr(object, "Dim") <- nlme::Dim(object, attr(object, "groups"))
  } else {
    # no groups
    attr(object, "Dim") <- nlme::Dim(object, as.factor(rep(1, nrow(data))))
  }
  ## obtaining the covariate(s)
  ## is this where the distance matrix is actually computed ?
  attr(object, "covariate") <- getCovariate.corMatern(object, data = data)
  
  nug  <- attr(object, "nugget")
  #nusc <- attr(object, "nuScaled")
  
  val <- as.vector(object)            # how many parameters?
  if (length(val) > 0) {		# is initialized
    if (val[1] <= 0) {                  # test for values of range
      stop("'range' must be > 0 in \"corMatern\" initial value")
    }
    if (val[2] <= 0) {                  # test for values of nu
      stop("'nu' must be > 0 in \"corMatern\" initial value")
    }    
    if (nug) {				# with nugget effect
      if (length(val) == 2) {		# assuming nugget effect not given
        val <- c(val, 0.1)		# setting it to 0.1
      } else {
        if (length(val) != 3) {
          stop("initial value for \"corMatern\" parameters of wrong dimension")
        }
      }
      if ((val[3] <= 0) || (val[3] >= 1)) {
        stop("initial value of nugget ratio must be between 0 and 1")
      }
    } else {				# only range and nu
      if (length(val) != 2) {
        stop("initial value for \"corMatern\" parameters of wrong dimension")
      }
    }
  } else { ## val of length 0
    val <- min(unlist(attr(object, "covariate"))) / 0.9
    val <- c(val,0.5) ## FR 18/03/13
    if (nug) val <- c(val, 0.1)
  }
  val[1] <- log(val[1]) ## log(range=1/rho)
  val[2] <- log(val[2]) ## log(nu)    
  if (nug) val[3] <- log(val[3]/(1 - val[3]))
  oldAttr <- attributes(object)
  object <- val
  attributes(object) <- oldAttr
  attr(object, "minD") <- min(unlist(attr(object, "covariate")))
  attr(object, "factor") <- corFactor.corMatern(object)
  attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")
  object
} ## Initialize.corMatern


logDet.corMatern <- function(object, covariate = getCovariate(object), ...) {
  if (!is.null(aux <- attr(object, "logDet"))) {
    return(aux)
  }
  if (is.null(aux <- attr(object, "factor"))) {
    ## getting the transpose sqrt factor
    aux <- corMatrix(object, covariate = covariate, corr = FALSE)
  }
  if (is.null(aux1 <- attr(aux, "logDet"))) {
    ## checking for logDet attribute; if not present, get corr matrix
    aux <- corMatrix(object, covariate)
    if (data.class(aux) == "list") {    # by group
      sum(log(abs(unlist(lapply(aux, function(el) svd(el)$d)))))/2
    } else {
      sum(log(abs(svd(aux)$d)))/2
    }
  } else {
    -aux1
  }
}


corMatrix.corMatern <-  function(object, covariate = getCovariate(object), 
                                corr = TRUE, ## returns matrix, else returna cholesky factor + the logDet  
                                ...) {
  if (data.class(covariate) == "list") { ## groups are defined
    if (is.null(names(covariate))) {
      names(covariate) <- 1:length(covariate)
    }
    corD <- nlme::Dim(object, rep(names(covariate),
                                  unlist(lapply(covariate,
                                                function(el) round((1 + sqrt(1 + 8 * length(el)))/2)))))
  } else { ## no groups
    corD <- nlme::Dim(object, rep(1, round((1 + sqrt(1 + 8* length(covariate)))/2)))
  }
  
  par <- as.vector(object)
  rho <- exp(-par[1]);
  nu <- exp(par[2]);
  if (attr(object, "nugget")) {
    aux <- exp(par[3]);
    Nugget <- 1 - 1 / (1.0 + aux) 
  } else Nugget <- 0.
  lD <- NULL
  val <- vector("list",length(covariate))
  distmats <- attr(covariate,"distmats")
  val <- vector("list",length(distmats))
  if ( ! corr) lD <- numeric(length(covariate))
  for (it in seq_along(distmats)) {
    tmp <- MaternCorr(distmats[[it]], nu=nu, rho=rho,Nugget=Nugget)
    tmp <- as.matrix(tmp)
    diag(tmp) <- 1
    if (corr) { ## returns a list of correlation matrices
      val[[it]] <- tmp
    } else {
      val[[it]] <- solve(designL.from.Corr(tmp)) # t(solve(chol(tmp)))
      lD[it] <- determinant(val[[it]])$modulus
      lD <- sum(lD)
    }
  }
  if (corD[["M"]] > 1) {
    val <- lapply(val, function(el) {
      nel <- round(sqrt(length(el)))
      array(el, c(nel, nel))
    })
    names(val) <- names(corD[["len"]])
    val <- as.list(val)
  } else {
    val <- array(val, c(corD[["N"]], corD[["N"]]))
  }
  attr(val, "logDet") <- lD
  val
} ## corMatrix.corMatern

coef.corMatern <- function(object, unconstrained = TRUE, ...) {
  ##cat("tata",object[1],object[2],"\n")
  if (attr(object, "fixed") && unconstrained) {
    return(numeric(0))
  }
  val <- as.vector(object)
  if (length(val) == 0) {               # uninitialized
    return(val)
  }
  if (!unconstrained) {
    val <- exp(val)
    if (attr(object, "nugget")) val[3] <- val[3]/(1+val[3])
  }
  if (attr(object, "nugget")) names(val) <- c("range", "nu", "nugget")
  else names(val) <- c("range","nu")
  val
}


"coef<-.corMatern" <- function(object, ..., value) {
  if (length(value) != length(object)) {
    stop("cannot change the length of the parameter after initialization")
  }
  #  cat("tutu",object[1],"->",value[1],object[2],"->",value[2],"\n")
  object[] <- value
  aux <- corFactor.corMatern(object, ...)
  attr(object, "factor") <- aux
  attr(object, "logDet") <- - attr(aux,"logDet")
  object
} ## "coef<-.corMatern"

corFactor.corMatern <- function(object, ...) {
  corD <- nlme::Dim(object)
  ## parameter assumed in unconstrained form = log transf for unconstrained range
  par <- as.vector(object)
  rho <- exp(-par[1]);
  nu <- exp(par[2]);
  if (attr(object, "nugget")) {
    aux <- exp(par[3]);
    Nugget <- 1 - 1 / (1.0 + aux) 
  } else Nugget <- 0.
  distmats <- attr(attr(object,"covariate"),"distmats")
  lD <- 0
  val <- vector("list",length(distmats))
  for (it in seq_along(distmats)) {
    tmp <- MaternCorr(distmats[[it]], nu=nu, rho=rho,Nugget=Nugget)
    tmp <- as.matrix(tmp)
    diag(tmp) <- 1
    tmp <- solve(designL.from.Corr(tmp)) #  t(solve(chol(tmp)))
    lD <- lD + determinant(tmp)$modulus[[1L]]
    val[[it]] <- as.vector(tmp) 
  }
  val <- unlist(val)
  attr(val, "logDet") <- lD
  return(val)
} ## corFactor.corMatern


.Dim.corMatern <- function(object, groups, ...) {
  if (missing(groups)) return(attr(object, "Dim"))
  ugrp <- unique(groups)
  groups <- factor(groups, levels = ugrp)
  len <- table(groups)
  val <- list(N = length(groups),
              M = length(len),
              maxLen = max(len),
              sumLenSq = sum(len^2),
              len = len,
              start = NA,
              spClass=0)
  val[["start"]] <- c(0, 4, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
  val
}

recalc.corMatern <- function(object, conLin, ...) {
  ## parameter assumed in unconstrained form = log transf for unconstrained range
  rho <- exp(-object[1]);
  nu <- exp(object[2]);
  if (attr(object, "nugget")) {
    aux <- exp(object[3]);
    Nugget <- 1 - 1 / (1.0 + aux) # /* 1 - nugget */
  } else Nugget <- 0
  covariate <- getCovariate.corMatern(object)
  distmats <- attr(covariate,"distmats")
  val <- vector("list",length(distmats))
  cumranges <- cumsum(c(0,unlist(lapply(distmats,nrow))))
  for (it in seq_along(distmats)) {
    tmp <- MaternCorr(distmats[[it]], nu=nu, rho=rho,Nugget=Nugget)
    tmp <- as.matrix(tmp)
    diag(tmp) <- 1
    tmp <- solve(designL.from.Corr(tmp)) #  t(solve(chol(tmp)))
    yrange <- (cumranges[it]+1L):(cumranges[it+1L])
    conLin[["Xy"]][yrange,] <- tmp %*% conLin[["Xy"]][yrange,,drop=FALSE]
  }
  conLin[["logLik"]] <- conLin[["logLik"]] +attr(attr(object,"factor"),"logDet") 
  conLin
} ## recalc.corMatern



Variogram.corMatern <- function(object, distance = NULL, sig2 = 1, length.out = 50, FUN, ...) {
  if (is.null(distance)) {
    rangeDist <- range(unlist(getCovariate(object)))
    distance <- seq(rangeDist[1], rangeDist[2], length = length.out)
  }
  params <- coef(object, unconstrained = FALSE)
  if (length(params) == 2) {            # no nugget effect
    range  <- params[1]
    nu   <- params[2]
    nugg <- 0
  } else {                              # nugget effect
    range  <- params[1]
    nu   <- params[2]
    nugg <- params[3]
  }
  par <- as.vector(object)
  rho <- exp(-par[1]);
  nu <- exp(par[2]);
  zz <- 1-MaternCorr(distance, nu=nu, rho=rho,Nugget=0)
  val <- data.frame(variog = sig2 * (nugg + (1 - nugg) * zz),
                    dist = distance)
  class(val) <- c("Variogram", "data.frame")
  val
} ## Variogram.corMatern

if (FALSE) {
  library(nlme)
  set.seed(1)
  d <- data.frame(x = rnorm(50), y = rnorm(50))
  gls(distance ~ Sex, data=Orthodont, correlation = corExp(1,form = ~1 | Subject)) #This works
  
  # use corMatern(form = ~age | Subject) to have distances according to age ! 
  str(csage <- corMatern(form = ~age | Subject)) 
  str(csage <- Initialize(csage, data = Orthodont[1:8,])) 
  corFactor(csage)
  
  ## corExp
  str(ce1 <- corExp(1,form = ~1 | Subject)) 
  locdata <- Orthodont[0+(1:8),]
  str(ce1o <- Initialize(ce1, data = locdata)) 
  corMatrix(ce1o)
  set.seed(123); 
  Xy <- matrix(runif(nrow(locdata)))
  recalc(ce1o,list(Xy=Xy,logLik=0))
  coef(ce1o) <- -log(2) # - log(rho)
  Variogram(ce1o)

  ## ~ cor exp via corMatern
  str(cs1 <- corMatern(c(1,0.5),form = ~1 | Subject)) 
  str(cs1n <- Initialize.corMatern(cs1, data =locdata)) 
  corMatrix.corMatern(cs1n)
  recalc.corMatern(cs1n,list(Xy=Xy,logLik=0))
  Variogram.corMatern(cs1n)
  coef(cs1n) <- c(log(2),log(0.5)) # log(rho), log(nu)

}



