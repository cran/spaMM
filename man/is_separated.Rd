\name{is_separated}
\alias{is_separated}
\alias{separation}
\title{
Checking for (quasi-)separation in binomial-response model.
}
\description{
Separation occurs in binomial response models when a combination of the predictor variables perfectly predict a level of the response. In such a case the estimates of the coefficients for these variables diverge to (+/-)infinity, and the numerical algorithms typically fail. To anticipate such a problem, the fitting functions in \code{spaMM} try to check for separation by default, using the \code{lpSolveAPI} package which can also detect some borderline cases (\dQuote{quasi-separation}). If this package is not available, \code{spaMM} tries to use the \code{e1071} package (in a way which will not detect quasi-separation). The check may take much time, and is skipped if the \dQuote{problem size} exceeds a threshold defined by \code{spaMM.options(separation_max=<.>)}, in which case a message will tell users by how much they should increase \code{separation_max} to force the check (the definition of the \dQuote{problem size} differs betwwen the two methods and may be complicated). 

\code{is_separated} is a convenient interface to procedures from the \code{lpSolveAPI} package, which can be called explicitly by the user to check bootstrap samples (see Example in \code{\link{anova}}). 
}
\usage{
is_separated(x, y, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Design matrix for fixed effects.
}
  \item{y}{
  Numeric response vector
}
  \item{verbose}{
  Whether to print some messages or not.
}
}
\value{
Returns a boolean; \code{TRUE} means there is (quasi-)separation.
}
\references{
Konis, K. 2007. Linear Programming Algorithms for Detecting Separated Data in Binary Logistic Regression Models. DPhil Thesis, Univ. Oxford.
}
\seealso{
See also the 'safeBinaryRegression' package.
}
\keyword{regression}% use one of  RShowDoc("KEYWORDS")
