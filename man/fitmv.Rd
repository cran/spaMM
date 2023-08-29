\name{fitmv}
\alias{fitmv}
\alias{X2X}
\title{
Fitting multivariate responses
}
\description{
This function extends the \code{fitme} function to fit a joint model for different responses (following possibly different response families) sharing some random-effects, including a new type of random effect defined to exhibit correlations across different responses (see \code{\link{mv}}).
It is also possible to declare shared fixed-effect coefficients among different submodels, using the \code{X2X} argument.  
Only a few features available for analysis of univariate response may not yet work (see Details). 
}
\usage{
fitmv(submodels, data, fixed=NULL, init=list(), lower=list(), upper=list(), 
      control=list(), control.dist = list(), method="ML", init.HLfit=list(), 
      X2X=NULL, ...)
}
\arguments{
  \item{submodels}{
A list of sublists each specifying a model for each univariate response. The names given to each submodel in the main list are currently ignored. The names and syntax of elements within each sublist are those of a \code{fitme} call. In most cases, each sublist should not contain arguments whose names are those of formal arguments of \code{fitmv} itself (with the possible exception for \code{fixed}). 

\code{prior.weights} (or better, \code{weights.form}), if any, should be specified as part of a submodel. 
}
  \item{data}{
A data frame containing the variables in the response and the model formulas.
}
  \item{fixed}{A list of fixed values of the parameters controlling random effects. The syntax is that of the same argument in \code{fitme} (the optional \code{fixed} argument in each sublist of \code{submodels} may also be used but this feature may be confusing). Fixed \code{phi} values must be specified as a list, e.g., \code{fixed=list(phi=list("2"=0.1))} to set the value for the second submodel.   
}
  \item{init, lower, upper}{
Lists of initial values or bounds. The syntax is that of the same arguments in \code{fitme}. In these lists, random effects should be indexed according to their order of appearance in the total model (see Details). Any \code{init}, \code{lower}, or \code{upper} in a sublist of \code{submodels} will be ignored.   
}
  \item{control}{A list of control parameters, with possible elements as described for \code{\link{fitme}}} 
  \item{control.dist}{ See \code{control.dist} in \code{\link{HLCor}} }
  \item{method}{Character: the fitting method to be used, such as \code{"ML"}, \code{"REML"} or \code{"PQL/L"}. \code{"ML"} is the default, as for \code{fitme} and in contrast to \code{"REML"} for the other fitting functions. Other possible values of \code{HLfit}'s \code{method} argument are handled.
}
\item{init.HLfit}{ See identically named \code{\link{HLfit}} argument. }
\item{X2X}{ NULL, or a matrix \bold{M} by which one can specify, as \eqn{\beta=}\bold{M}\eqn{\beta^*}, fixed effects \eqn{\beta} with some coefficients shared between submodels, e.g. as shown in the \dQuote{Shared fixed effect} Example, where \eqn{\beta^*} has three distinct elements, and \eqn{\beta} has four elements including identical Intercept coefficients among the two submodels. The fixed-effect term \bold{X}\eqn{\beta} of the linear predictor thus takes the form \bold{XM}\eqn{\beta^*}, meaning that the default design matrix of the model \bold{X} is replaced by \bold{XM}. \bold{M} must have column names, labeling the \eqn{\beta^*} coefficients.}
    \item{\dots}{
Optional arguments passed to (or operating as if passed to) \code{\link{HLCor}}, \code{\link{HLfit}} or  \code{\link{mat_sqrt}}, for example  \code{control.HLfit} or the \code{covStruct}, \code{distMatrix}, \code{corrMatrix} or \code{adjMatrix} arguments of \code{HLCor}. }
}
\details{
\bold{Matching random effects across submodels, and referring to them};\cr
Random effects are recognized as identical across submodels by matching the formula terms. As shown in the Examples, if the two models formulas share the \code{(1|clinic)} term, this term is recognized as a single random effect shared between the two responses. But the \code{(1|clinic)} and \code{(+1|clinic)} terms are recognized as distinct random effects. In that case, the \code{init} argument \code{init=list(lambda=c('1'=1,'2'=0.5))} is shown to refer to these by names \code{1,2}... where the order is defined as the order of first appearance of the terms across the model formulas in the order of the \code{submodels} list. 
Alternatively, the syntax \code{fixed=list(lambda=c('clinic.1'=0.5,'clinic'=1))} works: this syntax makes order of input irrelevant but assumes that the user guesses names correctly (these are typically the names that appear in the summary of lambda values from the fit object or, more programmatically, \cr
\code{names(<fit object>$lambda.object$print_namesTerms)}). Finally, fixed values of parameters can \bold{also} be specified through each sub-model, with indices referring to the order of random effects with each model.

The matching of random-effect terms occurs after expansion of \code{\link{multIMRF}} terms, if any. This may have subtle consequences if two multIMRF terms differ only by their number of levels, as some of the expanded IMRF terms are then shared. 

\bold{Capacities and limitations}:\cr
Practically all features of models that can be fitted by \code{fitme} should be available: this includes all combinations of GLM response families, residual dispersion models, and all types of random-effect terms, whether autocorrelated or not. Among the arguments handled through the \ldots, \code{covStruct}, \code{distMatrix}, \code{corrMatrix} should be effective; \code{control.HLfit$LevenbergM} and \code{verbose=c(TRACE=TRUE)} will work but some other controls available in \code{fitme} may not. Usage of the \code{\link{REMLformula}} argument is restricted as it cannot be used to specify a non-standard REML correction (but the more useful \code{\link{keepInREML}} attribute for fixed fixed-effect coefficients is handled).

The \code{\link{multi}} family-like syntax for multinomial models should not be used, but \code{fitmv} could provide other means to model multinomial responses.

Most post-fit functions work, at least with default arguments. This includes point \code{predict}ion and prediction variances calculations sensu lato, including with \code{newdata}; but also \code{simulate}, \code{spaMM_boot}, \code{confint}, \code{anova}, \code{update_resp}, and \code{update}. The \code{re.form} argument now works for \code{predict} and \code{simulate}. Bootstrap computation may require special care for models where simulation of one response variable may depend on drawns of another one (see Hurdle model example in the \dQuote{Gentle introdution} to \pkg{spaMM},
\url{https://gitlab.mbb.univ-montp2.fr/francois/spamm-ref/-/blob/master/vignettePlus/spaMMintro.pdf}).

Prediction functions try to handle most forms of missing information in \code{newdata} (including information missing for a residual-dispersion model when predcitiosn fro mit are needed: see Examples). As information may be missing for some submodels but not others, different numbers of predictions are then returned for different submodels. As for univariate-response models, \code{predict} will return point predictions as a single 1-column matrix, here concatenating the prediction results of the different submodels. The \code{nobs} attribute specifies how may values pertain to each submodel. 
  

Some plotting functions may fail. \code{update.formula} fails (see \code{\link{update_formulas}} for details). \code{terms} returns a list, which is not usable by other base R functions. \code{stats::step} is a good example of resulting limitations, as it is currently unable to perform any sensible operation on \code{fitmv} output. \code{spaMM::MSFDR} which calls \code{stats::step} likewise fails. 
\code{multcomp::glht} fails. % one mcp2matrix2() wherein factor_contrasts(model) surely does not return correct info
% get_RLRsim_args also look into family...

A perhaps not entirely satisfying feature is that \code{simulate} by default stacks the results of simulating each submodel in a single vector. Some non-trivial reformatting may then be required to include such simulation results in a suitable \code{newdata} data frame with (say) sufficient information for prediction of all responses. The syntax \cr
\code{update_resp(<fit>, newresp = simulate(<fit>, ...), evaluate = FALSE)$data}\cr 
may be particularly useful to reformat simulation results in this perspective. 

\bold{Which arguments belong to \code{submodels}?}:\cr 
Overall, arguments specifying individuals submodels should go into \code{submodels}, while other arguments of \code{fitmv} should be those potentially affecting several submodels (notably, random-effect structures, \code{lower}, and \code{upper}) and fitting controls (such as \code{init} and \code{init.HLfit}). One rarely-used exception is \code{\link{REMLformula}} which controls the fitting method but should be specified through the \code{submodels}.    

The function proceeds by first preprocessing all submodels independently, before merging the resulting information by matching random effects across submodels. The merging operation includes some checks of consistency across submodels, implying that redundant arguments may be needed across submodels (e.g. specifying twice a non-default \code{rand.family} for a random effect shared by two submodels). 

%Some fits may by default be slower than is possible by playing with control arguments \code{init} and \code{init.HLfit}.
}
\value{
A (single) list of class \code{HLfit}, as returned by other fitting functions in \pkg{spaMM}. The main difference is that it contains a \code{families} element describing the response families, instead of the \code{family} elements of fitted objects for univariate response.
}
\seealso{
See further examples in \code{\link{mv}} (modelling correlated random effects over the different submodels),
and \code{\link{residVar}}.
}
\examples{
### Data preparation
data(clinics)
climv <- clinics
(fitClinics <- HLfit(cbind(npos,nneg)~treatment+(1|clinic),
                     family=binomial(),data=clinics))
set.seed(123)
climv$np2 <- simulate(fitClinics, type="residual")
#
### fits

#### Shared random effect
(mvfit <- fitmv(
   submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                  mod2=list(formula=np2~treatment+(1|clinic),
                            family=poisson(), fixed=list(lambda=c("1"=1)))), 
   data=climv))

# Two univariate-response independent fits because random effect terms are distinct
# (note how two lambda values are set; same syntax for 'init' values):   
(mvfitind <- fitmv(
   submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                  mod2=list(formula=np2~treatment+(+1|clinic),family=poisson())), 
   data=climv, fixed=list(lambda=c('1'=1,'2'=0.5)))) # '1': (1|clinic); '2': (+1|clinic)

#### Specifying fixed (but not init) values in submodels is also possible (maybe not a good idea)
# (mvfitfix <- fitmv(
#   submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),
#                            family=binomial(),fixed=list(lambda=c('1'=1))), # '1': (1|clinic) 
#                  mod2=list(formula=np2~treatment+(+1|clinic),family=poisson(),
#                            fixed=list(lambda=c('1'=0.5)))),              # '2': (+1|clinic) 
#   data=climv))

#### Shared fixed effect
# Suppose we want to fit the same intercept for the two submodels 
# (there may be cases where this is meaningful, even if not here).
# The original fit has four coefficients corresponding to four columns 
# of fixed-effect design matrix:

head(design_X <- model.matrix(mvfit))
#      (Intercept)_1 treatment_1 (Intercept)_2 treatment_2
# [1,]             1           1             0           0
#      ...

# The three coefficients of the intended model are (say) 
# "(Intercept)" "treatment_1" "treatment_2"
# We build a matrix that relates the original 4 coefficients to these 3 ones:

X_4to3 <- 
  matrix(c(1,0,0,
           0,1,0,
           1,0,0,
           0,0,1), nrow=4, ncol=3, byrow=TRUE,
         dimnames=list(NULL, c("(Intercept)","treatment_1","treatment_2")))
                   
# defined such that design_X \%*\% X_4to3 will be the design matrix for the intended model, 
# and the single "(Intercept)" coefficient of the three-parameter model will operate as 
# a shared estimate of the "(Intercept)_1" and "(Intercept)_2" coefficients
# of the original 4-coefficients model, as intended. 

# The new fit is obtained by providing the matrix as the 'X2X' argument:

(mvfit3 <- fitmv(
    submodels=list(mod1=list(formula=cbind(npos,nneg)~treatment+(1|clinic),family=binomial()),
                   mod2=list(formula=np2~treatment+(1|clinic),
                             family=poisson(), fixed=list(lambda=c("1"=1)))), 
    X2X = X_4to3,
    data=climv))

# => the column names of 'X_4to3' are the fixed-effect names in all output.

#### Prediction with a residual-dispersion model
set.seed(123)
beta_dat <- data.frame(y=runif(100),grp=sample(2,100,replace = TRUE), x_het=runif(100),
                       y2=runif(100))
(mvfit <- fitmv(list(list(y ~1+(1|grp), family=beta_resp(), resid.model = ~x_het),
                   list(y2 ~1+(1|grp), family=beta_resp())), 
              data= beta_dat))
              
misspred <- beta_dat[1:3,]
misspred$x_het[1] <- NA # missing info for residual variance of first submodel

## => prediction missing when this info is needed:
#
length(predict(mvfit, newdata=misspred)) # 6 values: missing info not needed for point predictions
length(get_residVar(mvfit, newdata=misspred)) # 5 values  
length(get_respVar(mvfit, newdata=misspred)) # 5 values
#  Missing info not needed for predVar (**as opposed to respVar**)
length(get_predVar(mvfit, newdata=misspred)) # 6 values
#
# Same logic for interval computations:
#
dim(attr(predict(mvfit, newdata=misspred, intervals="respVar"),"intervals")) # 5,2  
dim(attr(predict(mvfit, newdata=misspred, intervals="predVar"),"intervals")) # 6,2  
#
# Same logic for simulate():
#
length(simulate(mvfit, newdata=misspred)) # 5 as simulation requires residVar


}


\keyword{ model }
