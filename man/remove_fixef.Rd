\name{remove_fixef}
\alias{remove_fixef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove specified fixed-effect terms from a model formula 
}
\description{
As the Title says.
}
\usage{
remove_fixef(formula, remove, nofixef = NULL, keep_intercept = TRUE, keep_offset = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
The model formula from which terms should be removed.
}
  \item{remove}{
The fixed-effect terms to be removed, specified either as a formula from which fixed-effect terms will be deduced (any random effet in such a formula will be ignored), or as a vector of character strings. The variable names should match those in the  \code{formula} argument (not to be confused with the labels of fitted coefficients of fixed effect, when these differ, as for factor terms). Permutation of variables in interaction terms are handled (see Examples).
}
  \item{nofixef}{
Optional \code{\link{expression}} or formula allowing control of terms included in the returned formula, beyond fixed effects not first removed according to the \code{remove} argument (see Details).
By default, these additional terms are all the offsets and random effects from the original formula.
See Details for non-default usage.
}
  \item{keep_intercept}{
Boolean: whether to keep an original intercept in the returned formula (does not add one if not originally present).
}
  \item{keep_offset}{
Boolean: whether to keep the offsets in the returned formula.
}
}
\details{
The following toy example explains the results with a non-default \code{nofixef}:\cr
\preformatted{
  form <- ~ offset() + a*b + (1|grp)
  remove_fixef(form, "a:b", nofixef=form)
}
which gives \code{~1 + a + b + offset(o) + a * b + (1 | grp)}.

First, fixed effect terms are extracted from the \code{formula}. Next, the terms implied by \code{remove} are removed, yielding an expression containing only other fixed-effect terms (here, \code{1 + a + b}). Then, this expression is joined to the terms included in \code{nofixef}. The \code{nofixef} argument by default contains only the offset and random effect terms from the original formula (hence its name), but this example shows that anything in \code{nofixef}, including fixed effects, is included in the result. 
}
\value{
A model formula. The response variable is retained in it if it was present in the input formula.
}
\examples{
remove_fixef( y ~ a*b, ~ a*b) #  y ~ 1, different from:
# remove_fixef(y ~ a*b, "a*b") # Warning and nothing removed 
# (because a*b not in terms(formula), which instead contains a, b and a:b)

remove_fixef( y ~ a*b, ~ a*b) # => y ~ 1, same result as:
remove_fixef( y ~ a*b, ~ b*a)   

remove_fixef( y ~ a*b, "a:b") # Same result as:
remove_fixef( y ~ a*b, "b:a") # a:b removed

# Remove all fixed effects term by setting 'remove' equal to 'formula':
form <- ~ offset() + a*b + (1|grp)
remove_fixef(form, form)
remove_fixef(form, form, keep_intercept = FALSE, keep_offset = FALSE)
}
