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
The \code{lmerTest::contest} function and \code{::anova} method implement a number of tests for linear mixed models, e.g. using effective degrees of freeedom based on (a generalization of) Satterthwaite's method. These tests can be performed using \pkg{spaMM} fits through the conversion of the fit object, by the \code{as_LMLT} function. Only the tests using \pkg{lmerTest}'s default method \code{ddf="Satterthwaite"} are formally supported, as the converted object do not have therequired format for the other methods. Only LMMs are handled by \pkg{lmerTest}, and residual-dispersion models are not yet handled by the conversion.  
}
\usage{
as_LMLT(fitobject, nuisance=NULL, verbose=TRUE, transf=TRUE, ...)
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

}
}
