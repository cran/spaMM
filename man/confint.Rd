\name{confint.HLfit}
\alias{confint.HLfit}
\alias{confint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals for fixed-effect parameters 
}
\description{
This computes confidence intervals for a given fixed effect parameter, based on the \code{p_v}-based approximation 
of the profile likelihood ratio for this parameter. The profiling is other all other fitted parameters: other fixed effects, as well as variances of random effects and spatial correlations if these were fitted.  
}
\usage{
\method{confint}{HLfit}(object, parm, level=0.95, verbose=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{HLfit}, such as return object of HLfit, HLCor or corrHLfit calls;}
  \item{parm}{The name of a parameter to be fitted, or its position in the the object's \code{$fixef} vector. 
  Valid names are those of the object's \code{$fixef}; }
  \item{level}{The coverage of the interval;}
  \item{verbose}{whether to print the interval or not. As the function returns its more extensive results invisibly, 
  this printing is the only visible output;}
  \item{\dots}{Additional arguments (maybe not used, but conforming to the generic definition of \code{confint}).}
}
\value{ A list including the confidence \code{interval} for the target parameter, and the 
fits \code{lowerfit} and \code{upperfit} giving the profile fits at the confidence bounds.
This is returned invisibly.
}
\examples{
% checked in test-confint.R
\dontrun{
 data("wafers")
 wfit <- HLfit(y ~X1+(1|batch),family=Gamma(log),data=wafers,HLmethod="ML")
 confint(wfit,"X1")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
