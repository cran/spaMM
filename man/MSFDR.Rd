\name{MSFDR}
\alias{MSFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple-Stage False Discovery Rate procedure
}
\description{
This implements the procedure described by Benjamini and Gavrilov (2009) for model-selection \bold{of fixed-effect terms} based on False Discovery Rate (FDR) concepts. It uses forward selection based on penalized likelihoods. The penalization for the number of parameters is distinct from that in Akaike's Information Criterion, and variable accross iterations of the algorithm (but functions from the \code{stats} package for AIC-based model-selection are still called, so that some screen messages refer to AIC). 
}
\usage{
MSFDR(nullfit, fullfit, q = 0.05, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nullfit}{
An ML fit to the minimal model to start the forward selection from; an object of class \code{HLfit}.
}
  \item{fullfit}{
An ML fit to the maximal model; an object of class \code{HLfit}.
}
  \item{q}{
Nominal error rate of the underlying FDR procedure (expected proportion of
incorrectly rejected null out of the rejected). Benjamini and Gavrilov (2009) recommend \code{q=0.05} on the basis of minimizing mean-squared prediction error in various simulation conditions considering only linear models. 
}
  \item{verbose}{
Whether to print information about the progress of the procedure.
}
}
\value{
The fit of the final selected model; an object of class \code{HLfit}.
}
\references{
A simple forward selection procedure based on false discovery rate control. Ann. Appl. Stat, 3, 179-198 (2009).
}
\examples{
if (spaMM.getOption("example_maxtime")>1.4) {
data("wafers")
nullfit <- fitme(y~1+(1|batch), data=wafers,family=Gamma(log))
fullfit <- fitme(y ~X1+X2+X1*X3+X2*X3+I(X2^2)+(1|batch), data=wafers, family=Gamma(log))
MSFDR(nullfit=nullfit,fullfit=fullfit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

