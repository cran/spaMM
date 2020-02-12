\name{eval_replicate}
\alias{eval_replicate}
\alias{.eval_replicate2}
\title{
Evaluating bootstrap replicates
}
\description{
\code{eval_replicate} is the default \code{simuland} function applied to simulated bootstrap samples by likelihood-ratio testing functions (\code{fixedLRT, LRT, anove.HLfit}). This documentation presents the requirements and possible features of this function and of possible user-defined alternatives. 

An alternative function \code{spaMM:::.eval_replicate2} is also provided. It is slower, as it refits the models compared with different initial values for random-effect parameters, which is useful in some difficult cases where initial values matter. The \code{eval_replicate} function may also refit the \dQuote{full} models with different initial values when the logLik of the refitted full model is subtstantially lower than that of the refitted null model. \dQuote{Substantially} means that a tolerance of \code{1e-04} is applied to account for inaccuracies of numerical maximization. 
}
\usage{
eval_replicate(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a response vector on which a previously fitted model may be refitted.   
}
}
\value{ 
A vector of the form \code{c(full=logLik(}<refitted full model>\code{),null=logLik(}<refitted null model>\code{)}; or possibly in debugging contexts, a list with the same elements each with some additional information provided as attribute.
}
\details{
likelihood-ratio testing functions have a \code{debug.} argument whose effect depends on the \code{simuland} function. The default behaviour is thus defined by \code{eval_replicate}, as: if \code{debug.=TRUE}, upon error in the fitting procedures, \code{dump.frames} will be called, in which case \bold{a dump file will be written on disk}; and a \bold{list} with debugging information will be returned (so that, say, \code{pbapply} will not return a matrix). This behaviour may change in later versions, so non-default \code{debug.} values should not be used in reproducible code. In serial computation, \code{debug.=2} may induce a \code{stop}; this should not happen in parallel computation because the calling functions check \code{debug.} values to prevent it.

Essential information such as the originally fitted models is passed to the function not as arguments but through its environment, which is controlled by the calling functions (see the \code{eval_replicate} source code to know which are these arguments). Users should thus not assume that they can control their own \code{simuland} function's environment as this environment will be altered.

Advanced users can define their own \code{simuland} function. The \code{eval_replicate} source code provides a template showing how to use the function's environment. The Example below illustrates another approach augmenting \code{eval_replicate}. A further example is provided in the file \code{tests/testthat/test-LRT-boot.R}, using \ldots to pass additional arguments beyond response values. 
}
\seealso{
Calling functions \code{\link{fixedLRT}, \link{LRT}}.
}
\examples{
\dontrun{
# Simple wrapper enhancing the default 'simuland'
#  with a call to some obscure option, and dealing with 
#  the need to pass the environment assigned to 'simuland'
eval_with_opt <- function(y) { 
  spaMM.options(some_obscure_option="some_obscure_value")
  eval_rep <- spaMM:::.eval_replicate
  environment(eval_rep) <- parent.env(environment()) # passing the environment
  eval_rep(y)
}
}
}
