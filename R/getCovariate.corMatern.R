getCovariate.corMatern <-
function(object, form = formula(object), data)
{
  if (is.null(covar <- attr(object, "covariate"))) { # need to calculate it
    if (missing(data)) {
      stop("need data to calculate covariate")
    }
    covForm <- getCovariateFormula(form)
    if (length(all.vars(covForm)) > 0) { # covariate present
      if (attr(terms(covForm), "intercept") == 1) {
	covForm <-
          eval(parse(text = paste("~", deparse(covForm[[2]]),"-1",sep="")))
      }
      covar <-
          as.data.frame(unclass(model.matrix(covForm,
                                             model.frame(covForm, data,
                                                         drop.unused.levels = TRUE))))
    } else {
      covar <- NULL
    }

    if (!is.null(getGroupsFormula(form))) { # by groups
      grps <- getGroups(object, data = data)
      if (is.null(covar)) {
	covar <- lapply(split(grps, grps),
                        function(x) as.vector(dist(1:length(x))))
      } else {
	covar <- lapply(split(covar, grps),
			function(el, metric) {
                          el <- as.matrix(el)
                          if (nrow(el) > 1) {
                            as.vector(dist(el, metric))
                          } else {
                            numeric(0)
                          }
			}, metric = attr(object, "metric"))
      }
      covar <- covar[sapply(covar, length) > 0]  # no 1-obs groups
    } else {				# no groups
      if (is.null(covar)) {
	covar <- as.vector(dist(1:nrow(data)))
      } else {
	covar <- as.vector(dist(as.matrix(covar),
                                method = attr(object, "metric")))
      }
    }
    if (any(unlist(covar) == 0)) {
      stop("cannot have zero distances in \"corMatern\"")
    }
  }
  covar
}
