\name{good-practice}
\alias{good-practice}
\alias{formula_env}

\title{Clear and trustworthy formulas}

\description{
Base fitting functions in R will seek variables in the environment where the \code{formula} was defined was made (i.e., typically in the global environment), if they are not in the \code{data}. This increases the memory size of fit objects (as the formula and attached environment are part of such objects). This also easily leads to errors (see example in the discussion of \code{\link{update.HLfit}}). Indeed Chambers (2008, p.221), after describing how the environment is defined, comments that \dQuote{Where clear and trustworthy software is a priority, I would personally avoid such tricks. Ideally, all the variables in the model frame should come from an explicit, verifiable data source...}. Fitting functions in spaMM try to adhere to such a principle, as they assume by default that all variables from the \code{formula} should be in the \code{data} argument (and then, \bold{one never needs to specify \dQuote{\code{data$}} in the \code{formula}.}. The variables defining the \code{prior.weights} should also be in the \code{data}. However, variables used in other arguments such as \code{ranFix} are looked up neither in the data nor in the formula environment, but in the calling environment as usual.

spaMM implements this by default by stripping the formula environment from any variable. It is also possible to assign a given environment to the formula, through the control \code{control.HLfit$formula_env}: see Examples. However, the search mechanism of R is such that variables present in the \code{formula} but not in the data nor in the formula environment will still be sought in the global environment, so bugs are not entirely preventable.
}
\references{
Chambers J.M. (2008) Software for data analysis: Programming with R. Springer-Verlag New York
}
\examples{
set.seed(123)
d2 <- data.frame(y = seq(10)/2+rnorm(5)[gl(5,2)], x1 = sample(10), grp=gl(5,2), seq10=seq(10))
# Using only variables in the data:
HLfit(y ~ x1 + seq10+(1|grp), data = d2, 
      control.HLfit = list(formula_env=list2env(list(data=d2))))
# 
HLfit(y ~ x1 + I(seq_len(nrow(data)))+(1|grp), data = d2, 
      control.HLfit = list(formula_env=list2env(list(data=d2))))
#
# Internal imprementation exploits partial matching of argument names
#  so that this can also be in 'control' if 'control.HLfit' is absent:      
fitme(y ~ x1 + I(seq_len(nrow(data)))+(1|grp), data = d2, 
      control = list(formula_env=list2env(list(data=d2))))
}