\name{anova}
\alias{anova}
\alias{anova.HLfit}
\alias{anova.LMLT}
\title{
  ANOVA tables (and likelihood ratio tests).
}
\description{
The \code{anova} method for fit objects from \pkg{spaMM} has two uses: if a single fit object is provided, ANOVA tables may be returned, with specific procedures for univariate-response LMs, GLMs and LMMs (see Details). Alternatively, if a second fit object is provided (\code{object2} argument), \code{anova} performs as an alias for \code{\link[spaMM]{LRT}}.  
}
\usage{
\method{anova}{HLfit}(object, object2, type = "2", method="", ...)
}
\arguments{
\item{object}{Fit object returned by a \pkg{spaMM} fitting function.}
\item{object2}{Optional second model fit to be compared to the first (their order does not matter, except in non-standard cases where the second model is taken as the null one and a message is issued).}
\item{type}{ANOVA type for LMMs, as interpreted by \pkg{lmerTest}. Note that the default (single-term deletion ANOVA) differs from that of \pkg{lmerTest}.} 
%
\item{method}{
  Only non-default value is \code{"t.Chisq"} which forces evaluation of a table of chi-squared tests for each fixed-effect term, using the classical \dQuote{Wald} test (see Details). Further methods are available through the \dots for specific classes of models.
} 
%
\item{\dots}{
  Further arguments, passed to \code{\link{spaMM_boot}} (e.g., for parallelization) in the case of LRTs (i.e., when \code{object2} is provided). For ANOVA tables, arguments of functions \code{\link[stats]{anova.lm}} \code{\link[stats]{anova.glm}}, 
  and \code{\link{as_LMLT}}, % including, notably, its \code{check} argument
  respectively for LMs, GLMs and LMMs, may be handled (e.g. the \code{test} argument for \code{anova.glm}).
}
}
\details{
The ANOVA-table functionality has been included in \pkg{spaMM} mainly to provide access to F tests (including, for LMMs, the \dQuote{Satterthwaite method} as developed by Fai and Cornelius, 1996), using pre-existing procedures as template or backend for expediency and familiarity:\cr 
  \enumerate{
 \item ANOVA tables \bold{for LMs and GLMs} have been conceived to replicate the functionality, output format and details of base R \code{anova}, and therefore replicate some of their limitations, e.g., they only perform sequential analysis (\dQuote{type 1}) in the same way as \code{anova.lm} and \code{anova.glm}. However, a difference occurs for Gamma GLMs, because the dispersion estimates for Gamma GLMs differ between \code{stats::glm} and \pkg{spaMM} fits (see Details in \code{\link{method}}). Therefore, F tests and Mallows' Cp differ too; results from \pkg{spaMM} REML fits being closer than ML fits to those from \code{glm()} fits;\cr \item \bold{For LMMs}, ANOVA tables are provided by interfacing \code{lmerTest::anova} (with non-default \code{type}). This procedure should handle all types of LMMs that can be fitted by \pkg{spaMM}; yet, the displayed information should be inspected to check that some fitted random-effect parameters are not ignored when computing information for the Satterthwaite method. 
\item  For fitted models that do not lay within previous categories, such as \bold{GLMMs}, models with a \bold{residual-dispersion} submodel, and \bold{multivariate-response} models, a table of tests for single-term deletions using the classical \dQuote{Wald} chi-squared test based on coefficient values and their conditional standard error estimates will be returned. LRTs (moreover, with bootstrap correction) are more reliable than such tests and, as calling them requires a second model to be explicitly specified, they may also help users thinking about the hypothesis they are testing. 
  }

}
\value{
The return format is that of the function called (\code{lmerTest::anova} for LMMs) or emulated (for LMs or GLMs). For GLMs, by default no test is reported, as has been the default for \code{anova.glm} before R 4.4.0.
} %- value
\references{
Fai AH, Cornelius PL (1996). Approximate F-tests of multiple degree of freedom hypotheses in generalised least squares analyses of unbalanced split-plot experiments.
Journal of Statistical Computation and Simulation, 54(4), 363-378. \doi{10.1080/00949659608811740}
}
\seealso{
\code{\link[spaMM]{LRT}},\cr
\code{\link{as_LMLT}} for the interface to \code{lmerTest::anova},\cr
and \code{\link{summary.HLfit}}\code{(.,details=list(<true|"Wald">))} for reporting the p-value for each t-statistic in the summary table for fixed effects, either by Student's t distribution, or by the approximation of t^2 distribution by the Chi-squared distribution (\dQuote{Wald's test}).
}
\examples{
data("wafers")
## Gamma GLMM with log link
m1 <- HLfit(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),family=Gamma(log),
          resid.model = ~ X3+I(X3^2) ,data=wafers,method="ML")
m2 <- update(m1,formula.= ~ . -I(X2^2))
#
anova(m1,m2) # LRT

## 'anova' (Wald chi-squared tests...) for GLMM or model with a 'resid.model'
anova(m1) 

## ANOVA table for GLM
# Gamma example, from McCullagh & Nelder (1989, pp. 300-2), as in 'glm' doc:
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100),
    lot1 = c(118,58,42,35,27,25,21,19,18),
    lot2 = c(69,35,26,21,18,16,13,12,12))
spglm <- fitme(lot1 ~ log(u), data = clotting, family = Gamma, method="REML")
anova(spglm, test = "F") 
anova(spglm, test = "Cp") 
anova(spglm, test = "Chisq")
anova(spglm, test = "Rao") 

## ANOVA table for LMM
if(requireNamespace("lmerTest", quietly=TRUE)) {
  lmmfit <- fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch),data=wafers)
  print(anova(lmmfit)) # => Satterthwaite method, here giving p-values 
                       #   quite close to traditional t-tests given by:
  summary(lmmfit, details=list(p_value=TRUE))
}
}
\keyword{ htest }
