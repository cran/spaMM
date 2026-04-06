\name{X2X}
\alias{genX2X}
\alias{X2X}
\title{
Fitting fixed-effects coefficients shared between submodels
}
\description{
\code{X2X} is an argument of the \code{\link{fitmv}} function, whose use is illustrated in the Examples below. By providing an \code{X2X} matrix \bold{M}, the default fixed-effect design matrix \bold{X} of the model is replaced by \bold{XM}. The fixed-effect term \bold{X}\eqn{\beta} of the linear predictor is then modified to \bold{XM}\eqn{\beta^*} for a new vector \eqn{\beta^*} of fixed-effect coefficients, distinct from the vector \eqn{\beta} of coefficients of the default fit.

In the examples, the vector \eqn{\beta} of fixed-effect coefficients of the default fit, without \code{X2X} argument, has four distinct elements, including one Intercept for each of the two sub-models. The \bold{M} matrix is used to declare that the Intercept coefficient is identical in the two sub-models, so \bold{M} has three columns and four rows. 

\bold{M} must have column names, labeling the \eqn{\beta^*} coefficients.

Because specifying the \bold{M} matrix can be laborious and error-prone, the helper function \code{genX2X} is provided. The typical way of using it, illustrated below, is \code{X2X=genX2X(matches=\var{list})}, with the second argument of the function missing. This argument will be provided internally, so the user does not have to find its correct form. Without this argument, the function returns its own call, which is not of any immediate interest for the user.   
}
\usage{
genX2X(matches, names_ori)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matches}{
A named list, whose elements are vectors of names (character strings). The list names are those of some of the coefficients to be fitted(notably, the names to be given to the new, shared coefficients); and the vector elements are some names of coefficients of the original model (the one without an \code{X2X} argument) that should be constrained to be identical to the corresponding shaed coefficient.
}
  \item{names_ori}{
Either missing, or the names of all coefficients of the default model (the one without \code{X2X} argument). 
}
}
\value{
\code{genX2X} returns a matrix if \code{names_ori} is provided, and returns the call to the \code{genX2X} function is this argument is missing.
}
\examples{
### Data preparation
data(clinics)
climv <- clinics
(fitClinics <- HLfit(cbind(npos,nneg)~treatment+(1|clinic),
                     family=binomial(),data=clinics))
set.seed(123)
climv$np2 <- simulate(fitClinics, type="residual")

### fits

## Default fit without 'X2X' argument
(mvfit <- fitmv(
   submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                  mod2=list(formula=np2~treatment+(1|clinic),
                            family=poisson(), fixed=list(lambda=c("1"=1)))), 
   data=climv))
   
## Fits with 'X2X' argument
   
# Suppose we want to fit the same intercept for the two submodels 
# (there may be cases where this is meaningful, even if not here).
# The original fit has four coefficients corresponding to four columns 
# of fixed-effect design matrix:

head(design_X <- model.matrix(mvfit))
#      (Intercept)_1 treatment_1 (Intercept)_2 treatment_2
# [1,]             1           1             0           0
#      ...
# where '_1' or '_2' identifies the submodel to which each coefficient belongs.

# The three coefficients of the intended model are (say) 
# "(Intercept)" "treatment_1" "treatment_2"
# We build a matrix that relates the original 4 coefficients to these 3 ones:

X_4to3 <- 
  matrix(c(1,0,0,
           0,1,0,
           1,0,0,
           0,0,1), nrow=4, ncol=3, byrow=TRUE,
         dimnames=list(NULL, c("(Intercept)","treatment_1","treatment_2")))
                   
# defined such that design_X \%*\% X_4to3 will be the design matrix 
# for the intended model, and the single "(Intercept)" coefficient 
# of the three-parameter model will operate as a shared estimate 
# of the "(Intercept)_1" and "(Intercept)_2" coefficients
# of the original 4-coefficients model, as intended. 

# To define such matrices, it is *strongly advised* to either fit 
# the unconstrained model first, and to examine the structure 
# of its model matrix (as shown above), or to use genX2X(). 

# The 'X2X' argument provides the matrix:

(mvfit3m <- fitmv(
    submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                   mod2=list(formula=np2~treatment+(1|clinic),
                             family=poisson(), fixed=list(lambda=c("1"=1)))), 
    X2X = X_4to3,
    data=climv))

# => the column names of 'X_4to3' are the fixed-effect names in all output.

# Alternatively, the 'X2X' argument provides a genX2X() call:

(mvfit3g <- fitmv(
    submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                   mod2=list(formula=np2~treatment+(1|clinic),
                             family=poisson(), fixed=list(lambda=c("1"=1)))), 
    X2X = genX2X(list("(Intercept)"=c("(Intercept)_1","(Intercept)_2"))),
    data=climv))

# => the last two fits are equivalent (although the order of the coefficients may differ).
# The internally produced 'X2X' matrix is that provided by 
genX2X(list("(Intercept)"=c("(Intercept)_1","(Intercept)_2")),
       names_ori=colnames(design_X) )

}
