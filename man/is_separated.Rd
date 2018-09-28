\name{is_separated}
\alias{is_separated}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checking for (quasi-)separation in binomial-response model.
}
\description{
This is a convenient interface to procedures from the lpSolveAPI package, if this package is installed (otherwise a cruder approach and possibly flawed will be used), to test for (quasi-)separation. This is used by default by the fitting funtions, but can also be called explicitly by the user to check bootstrap samples (see Example in \code{\link{anova}}). 
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
