\name{as_LMLT}
\alias{as_LMLT}
\alias{LMLTslots}
\alias{class:LMLTslots}
\alias{LMLTslots-class}
\alias{model.matrix.LMLTslots} 
%
\title{
Conversion to input for procedures from lmerTest package
}
\description{
The \code{lmerTest::contest} function, \code{drop1} and \code{anova} methods implement a number of tests for linear mixed models, e.g. using effective degrees of freeedom based on (a generalization of) Satterthwaite's method. These tests can be performed using \pkg{spaMM} fits through the conversion of the fit object, by the \code{as_LMLT} function, to an ad-hoc format acceptable as input to \pkg{lmerTest}'s internal procedures. The separately documented \code{\link{drop1.HLfit}} and (optionally) \code{\link{anova.HLfit}} methods, when called on a single LMM fit object, perform the conversion by \code{as_LMLT} and call \code{drop1} or \code{anova} methods defined by \pkg{lmerTest}.

Only the tests using \pkg{lmerTest}'s default method \code{ddf="Satterthwaite"} are formally supported, as the converted object do not have the required format for the other methods. Only LMMs are handled by \pkg{lmerTest}, and residual-dispersion models are not yet handled by the conversion. However, the conversion extends \pkg{lmerTest}'s functionality by handling all random-effect parameters handled by \code{\link{numInfo}}, therefore including (e.g.) spatial-correlation parameters not handled by \pkg{lme4}.  

}
\usage{
as_LMLT(fitobject, nuisance=NULL, verbose=TRUE, transf=TRUE, check_deriv=NULL, ...)
}
\arguments{
\item{fitobject}{
  Object of class \code{HLfit} resulting from the fit of a linear mixed model (LMM).
}
\item{nuisance}{
  A list of fitted values of parameters that affect the distribution of the test of fixed effects, in the format of the \code{fixed} argument of the \code{fitme} function. If NULL (default), then the list is constructed from the fitted values of the random-effect parameters and of \code{phi} (residual dispersion parameter). 
  The \code{nuisance} argument is better ingored unless the extractor he construct the default value fails in some way.
}
\item{verbose}{
  boolean: controls printing of the message that shows the \code{unlist}ed value of the \code{nuisance} list.
}
\item{transf}{
  boolean: whether to evaluate numerical derivatives on a transformed parameter scale, or not (may affect numerical precision).
}
\item{check_deriv}{See same-named argument of \code{\link{numInfo}}} 
\item{\dots}{Other arguments that may be needed by some method (currently ignored).}
}
\value{
The value is returned invisibly. It is an S4 object of class \code{"LMLT"} with slots matching those required in objects of S4 class \code{"lmerModLmerTest"} when used by package \pkg{lmerTest} with \code{ddf="Satterthwaite"} (many additional slots of a formal \code{"lmerModLmerTest"} object are missing). The additional \code{nuisance} slot contains the \code{nuisance} list.

}
\references{
Alexandra Kuznetsova, Per B. Brockhoff and Rune H. B. Christensen (2017) lmerTest Package: Tests in Linear Mixed Effects Models. Journal of Statistical Software, 82(13), 1–26. doi:10.18637/jss.v082.i13
}
\examples{
## Reproducing an example from the doc of lmerTest::contest.lmerModLmerTest,
#    using a spaMM fit as input.
\dontrun{
  data("sleepstudy", package="lme4")
  
  ## The fit:
  spfit <- fitme(Reaction ~ Days + I(Days^2) + (1|Subject) + (0+Days|Subject),
             sleepstudy, method="REML")
  
  ## Conversion:           
  spfit_lmlt <-  as_LMLT(spfit) 

  ## Functions from package lmerTest can then be called on this object:
  lmerTest::contest(spfit_lmlt, L=diag(3)[2:3, ])   # Test of 'Days + I(Days^2)'.
  #
  anova(spfit_lmlt, type="1")         # : using lmerTest:::anova.lmerModLmerTest()
  drop1(spfit_lmlt)                   # : using lmerTest:::drop1.lmerModLmerTest()

}
}
