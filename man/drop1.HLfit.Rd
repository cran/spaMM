\name{drop1.HLfit}
\alias{drop1}
\alias{drop1.HLfit}
\alias{drop1.LMLT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drop all possible single fixed-effect terms from a model
}
\description{
Drop single terms from the model. The \code{drop1} method for \pkg{spaMM} fit objects is conceived to replicate the functionality, output format and details of pre-existing methods for similar models. Results for LMs and GLMs should replicate base R \code{drop1} results, with some exceptions:\cr
* somewhat stricter default check of non-default \code{scope} argument;\cr
* Because the dispersion estimates for Gamma GLMs differ between \code{stats::glm} and \pkg{spaMM} fits (see Details in \code{\link{method}}), some tests may differ too; results from \pkg{spaMM} REML fits being closer than ML fits to those from \code{glm()} fits;\cr 
* AIC values reported in tables are always the marginal AIC as computed by \code{AIC.HLfit}, while \code{drop1.glm} may report confusing (to me, at least) values (see \code{\link{AIC.HLfit}}) for reasons that seem to go beyond differences in dispersion estimates.

\bold{For LMMs}, ANOVA tables are provided by interfacing \code{lmerTest::anova} (with non-default \code{type}).

For other classes of models, a table of likelihood ratio tests is returned, each test resulting from a call to \code{\link{LRT}}. 
}
\usage{
\method{drop1}{HLfit}(object, scope, method="", check_marg=NULL, check_time=60, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{Fit object returned by a \pkg{spaMM} fitting function.}
  \item{scope}{
  Default \dQuote{scope} (terms to be tested, specified as a formula, see Examples) is determined by applying \code{stats::drop.scope} on fixed-effect terms.
  Non-default scope can be specified a formula giving the terms to be considered for dropping. It is also possible to specify them as a character vector, but then one has to make sure that the elements are consistent with term labels produced by \code{terms}, as inconsistent elements will be ignored.
}
\item{method}{
  Only non-default value is \code{"LRT"} which forces evaluation of a likelihood ratio tests by \code{\link{LRT}}, instead of specific methods for specific classes of models.
} 
\item{check_marg}{
  NULL or boolean: whether effects should be checked for marginality. By default, this check is performed when a non-default scope is specified, and then no test is reported for terms that do not satisfy the marginality condition. If \code{check_marg=FALSE}, marginality is not checked and tests are always performed.
} 
\item{check_time}{
  numeric: whether to output some information when the execution time of \code{drop1} may be of the order of the time specified by \code{check_time}, or more. This is checked only when random effect are present. Such output can thus be suppressed by \code{check_time=Inf}.
} 
\item{\dots}{Further arguments passed from or to methods, or to \code{LRT}.}
}
\details{
As for the ANOVA-table functionality, it has been included here mainly to provide access to F tests (including, for LMMs, the \dQuote{Satterthwaite method}, using pre-existing procedures as template or backend for expediency and familiarity. The procedures for specific classes of models have various limitations, e.g., none of them handle models with variable dispersion parameter. For classes of models not well handled by these procedures (by design or due to the experimental nature of the recent implementation), \code{method="LRT"} can still be applied (and will be applied by default for GLMMs). 
}
\value{
The return format is that of the function called (\code{lmerTest::drop1} for LMMs), or emulated (base \code{drop1} methods for LMs or GLMs), or is a data frame whose rows are each the result of calling \code{LRT}.
}
\seealso{
\code{\link{as_LMLT}} for the interface to \code{lmerTest::drop1}.
}
\examples{
data("wafers")
#### GLM

wfit <- fitme(y ~ X1+X2+X1*X3+X2*X3+I(X2^2),  family=Gamma(log), data=wafers)
drop1(wfit, test = "F")              
drop1(wfit, test = "F", scope= ~ X1 +  X1 * X3 )  # note the message!

#### LMM
if(requireNamespace("lmerTest", quietly=TRUE)) {
  lmmfit <- fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),data=wafers)
  drop1(lmmfit) # => Satterthwaite method here giving p-values quite close to 
                #    traditional t-tests given by:
  summary(lmmfit, details=list(p_value=TRUE)) 
}

#### GLMM 
wfit <- fitme(y ~ X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), family=Gamma(log),
              rand.family=inverse.Gamma(log), resid.model = ~ X3+I(X3^2) , data=wafers)
drop1(wfit)              
drop1(wfit, scope= ~ X1 +  X1 * X3 )  # note the message!
}
