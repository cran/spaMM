\name{good-practice}
\alias{good-practice}
\alias{formula_env}

\title{Clear and trustworthy formulas and prior weights}

\description{
Base fitting functions in R will seek variables in the environment where the \code{formula} was defined (i.e., typically in the global environment), if they are not in the \code{data}. This increases the memory size of fit objects (as the formula and attached environment are part of such objects). This also easily leads to errors (see example in the discussion of \code{\link{update.HLfit}}). Indeed Chambers (2008, p.221), after describing how the environment is defined, comments that \dQuote{Where clear and trustworthy software is a priority, I would personally avoid such tricks. Ideally, all the variables in the model frame should come from an explicit, verifiable data source...}. Fitting functions in \pkg{spaMM} try to adhere to such a principle, as they assume by default that all variables from the \code{formula} should be in the \code{data} argument (and then, \bold{one never needs to specify \dQuote{\code{data$}} in the \code{formula}.}. \pkg{spaMM} implements this by default by stripping the formula environment from any variable. It is also possible to assign a given environment to the formula, through the control \code{control.HLfit$formula_env}: see Examples. 
% However, the search mechanism of R is such that variables present in the \code{formula} but not in the data nor in the formula environment will still be sought in the global environment, so bugs are not entirely preventable. => .GetValidData_info() ?

The variables defining the \code{prior.weights} should also be in the \code{data}. However, the implementation of the \code{prior.weights} argument has limitations that can be overcome by using the more recently introduced \code{weights.formula} argument of \pkg{spaMM} fitting functions (see Examples, where this is also compared  with \code{stats::lm}'s handling of its \code{weights} argument).  
% \code{data$} should be OK in prior.weights and can make programming easier when a fit is performed in a function from a package (else R CMD check would fail to ascertain the variable)

However, variables used in other arguments such as \code{ranFix} are looked up neither in the data nor in the formula environment, but in the calling environment as usual.

}
\references{
Chambers J.M. (2008) Software for data analysis: Programming with R. Springer-Verlag New York
}
\examples{
#######  Controlling the formula environment

set.seed(123)
d2 <- data.frame(y = seq(10)/2+rnorm(5)[gl(5,2)], x1 = sample(10), grp=gl(5,2), seq10=seq(10))
# Using only variables in the data: basic usage
# HLfit(y ~ x1 + seq10+(1|grp), data = d2)
# is practically equivalent to
HLfit(y ~ x1 + seq10+(1|grp), data = d2, 
      control.HLfit = list(formula_env=list2env(list(data=d2))))
# 
# The 'formula_env' avoids the need for the 'seq10' variable:
HLfit(y ~ x1 + I(seq_len(nrow(data)))+(1|grp), data = d2, 
      control.HLfit = list(formula_env=list2env(list(data=d2))))
#
# Internal imprementation exploits partial matching of argument names
#  so that this can also be in 'control' if 'control.HLfit' is absent:      
fitme(y ~ x1 + I(seq_len(nrow(data)))+(1|grp), data = d2, 
      control = list(formula_env=list2env(list(data=d2))))
      
####### Prior-weights misery

data("hills", package="MASS")

(fit <- lm(time ~ dist + climb, data = hills, weights=1/dist^2)) 
# same as
(fit <- fitme(time ~ dist + climb, data = hills, prior.weights=1/dist^2, method="REML")) 

# possible calls:
(fit <- fitme(time ~ dist + climb, data = hills, prior.weights=quote(1/dist^2))) 
(fit <- fitme(time ~ dist + climb, data = hills, prior.weights= 1/hills$dist^2)) 
(fit <- fitme(time ~ dist + climb, data = hills, weights.form= ~ 1/dist^2)) 
(fit <- fitme(time ~ dist + climb, data = hills, weights.form= ~ I(1/dist^2))) 

# Also syntactically correct since 'dist' is found in the data:
(fit <- fitme(time ~ dist + climb, data = hills, weights.form= ~ rep(2,length(dist)))) 

#### Programming with prior weights:

## Different ways of passing prior weights to fitme() from another function:

wrap_as_form <- function(weights.form) {
  fitme(time ~ dist + climb, data = hills, weights.form=weights.form) 
}

wrap_as_pw <- function(prior.weights) {
  fitme(time ~ dist + climb, data = hills, prior.weights=prior.weights) 
}

wrap_as_dots <- function(...) {
  fitme(time ~ dist + climb, data = hills,...) 
}


## Similarly for lm:

wrap_lm_as_dots <- function(...) {
  lm(time ~ dist + climb, data = hills, ...) 
}

wrap_lm_as_arg <- function(weights) {
  lm(time ~ dist + climb, data = hills,weights=weights) 
}

## Programming errors with stats::lm():

pw <- rep(1e-6,35) # or even NULL

(fit <- wrap_lm_as_arg(weights=pw)) # catches weights from global envir!
(fit <- lm(time ~ dist + climb, data = hills, weights=pw)) # idem!

(fit <- lm(time ~ dist + climb, data = hills, 
           weights=hills$pw)) # fails silently - no $pw in 'hills'
(fit <- wrap_lm_as_dots(weights=hills$pw)) # idem!
(fit <- wrap_lm_as_arg(weights=hills$pw)) # idem!

## Safer spaMM results:

try(fit <- wrap_as_pw(prior.weights= pw)) # correctly catches problem
try(fit <- wrap_as_dots(prior.weights=hills$pw)) # correctly catches problem
(fit <- wrap_as_dots(prior.weights=1/dist^2)) # correct
(fit <- wrap_as_dots(prior.weights=quote(1/dist^2))) # correct

## But 'prior.weights' limitations: 

try(fit <- wrap_as_pw(prior.weights= 1/hills$dist^2)) # fails (stop)
try(fit <- wrap_as_pw(prior.weights= 1/dist^2)) # fails (stop)
try(fit <- wrap_as_pw(prior.weights= quote(1/dist^2))) # fails (stop)

## Problems all solved by using 'weights.form':

try(fit <- wrap_as_form(weights.form= ~ pw)) # correctly catches problem
(fit <- wrap_as_form(weights.form= ~1/dist^2)) # correct
(fit <- wrap_as_form(weights.form= ~1/hills$dist^2)) # correct
(fit <- wrap_as_dots(weights.form= ~ 1/dist^2)) # correct

rm("pw")

      
      
}