\name{is_separated}
\alias{is_separated}
\alias{is_separated.formula}
\alias{separation}
\title{
Checking for (quasi-)separation in binomial-response model.
}
\description{
Separation occurs in binomial response models when a combination of the predictor variables perfectly predict a level of the response. In such a case the estimates of the coefficients for these variables diverge to (+/-)infinity, and the numerical algorithms typically fail. To anticipate such a problem, the fitting functions in \code{spaMM} try to check for separation by default. The check may take much time, and is skipped if the \dQuote{problem size} exceeds a threshold defined by \code{spaMM.options(separation_max=<.>)}, in which case a message will tell users by how much they should increase \code{separation_max} to force the check (its exact meaning and default value are subject to changes without notice but the default value aims to correspond to a separation check time of the order of 1s on the author's computer). 

\code{is_separated} is a convenient interface to procedures from the \code{ROI} package, which can be called explicitly by the user to check bootstrap samples (see Example in \code{\link{anova}}). 
\code{is_separated.formula} is a variant (not yet a formal S3 method) that performs the same check, but using arguments similar to those of \code{fitme(., family=binomial())}.
}
\usage{
is_separated(x, y, verbose = TRUE, solver=spaMM.getOption("sep_solver"))
is_separated.formula(formula, ..., separation_max=spaMM.getOption("separation_max"),
                     solver=spaMM.getOption("sep_solver"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Design matrix for fixed effects.
}
  \item{y}{
  Numeric response vector
}
\item{formula}{
  A model formula
}
\item{...}{
  \code{data} and possibly other arguments of a \code{fitme} call. \code{family} is ignored if present.  
}
\item{separation_max}{
  numeric: non-default value allow for easier local control of this spaMM option.  
}
\item{solver}{
  character: name of linear programming solver used to assess separation; passed to \code{\link[ROI]{ROI_solve}}'s \code{solver} argument. One can select other solvers if the corresponding ROI plugin is installed. 
}
\item{verbose}{
  Whether to print some messages or not.
}
}
\value{
Returns a boolean; \code{TRUE} means there is (quasi-)separation.
}
\references{
The method accessible by \code{solver="glpk"} implements algorithms described by

Konis, K. 2007. Linear Programming Algorithms for Detecting Separated Data in Binary Logistic Regression Models. DPhil Thesis, Univ. Oxford. \url{https://ora.ox.ac.uk/objects/uuid:8f9ee0d0-d78e-4101-9ab4-f9cbceed2a2a}.
}
\seealso{
See also the 'safeBinaryRegression' and 'detectseparation' package.
}
\examples{
set.seed(123)
d <- data.frame(success = rbinom(10, size = 1, prob = 0.9), x = 1:10)
is_separated.formula(formula= success~x, data=d) # FALSE
is_separated.formula(formula= success~I(success^2), data=d) # TRUE
}
\keyword{regression}% use one of  RShowDoc("KEYWORDS")
