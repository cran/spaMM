\name{spaMM}
% ----- \Rdversion{1.1}
\alias{spaMM}
\alias{spaMM-package}

\docType{package}

\title{Inference in mixed models, in particular spatial GLMMs}

\description{Fits a range of mixed-effect models, including those with spatially correlated random effects. The random effects are either Gaussian (which defines GLMMs), or other distributions (which defines the wider class of hierarchical GLMs), or simply absent (which makes a LM or GLM).  Multivariate-response models can be fitted by the \code{\link{fitmv}} function. Other models can be fitted by \code{\link{fitme}} (the most general function). Also available are previously conceived fitting functions \code{\link{HLfit}} (sometimes faster, for non-spatial models), \code{\link{HLCor}} (sometimes faster, for conditional-autoregressive models and fixed-correlation models), and \code{\link{corrHLfit}} (now of lesser interest). Additional functions are available such as \code{\link{fixedLRT}} for likelihood-ratio testing, \code{\link{simulate}} and \code{\link{predict}}.

Both maximum likelihood (ML) and restricted likelihood (REML) can be used for linear mixed models, and extensions of these methods using Laplace approximations are used for non-Gaussian random response. Several variants of these methods discussed in the literature are included (see Details in \code{\link{HLfit}}), the most notable of which may be \dQuote{PQL/L} for binary-response GLMMs (see Example for \code{\link{arabidopsis}} data). PQL methods implemented in spaMM are closer to (RE)ML methods than those implemented in \code{MASS::glmmPQL}. 
}

\details{
  The standard response families \code{gaussian}, \code{binomial}, \code{poisson}, and \code{Gamma} are handled, as well as negative binomial (see \code{\link{negbin}}), zero-truncated poisson and negative binomial, and Conway-Maxwell-Poisson response (see \code{\link{Tpoisson}}, \code{\link{Tnegbin}} and \code{\link{COMPoisson}}). A \code{multi} family look-alike is also available for \code{\link{multinomial}} response, with some constraints. 
  
  The variance parameter of residual error is denoted \eqn{\phi} (\code{phi}): this is the residual variance for gaussian response, but for Gamma-distributed response, the residual variance is \eqn{\phi}\eqn{\mu^2} where \eqn{\mu} is expected response. A fixed-effects linear predictor for \eqn{\phi}, modeling heteroscedasticity, can be considered (see Examples). 
  
  The package fits models including several nested or crossed random effects, including autocorrelated ones with the following correlation models: \code{\link{Matern}}, \code{\link{Cauchy}}, interpolated Markov Random Fields (\code{\link{IMRF}}), first-order autoregressive (\code{\link{AR1}}), conditional autoregressive as specified by an \code{\link{adjacency}} matrix, or any fixed correlation matrix (\code{\link{corrMatrix}}). GLMMs and HGLMs are fit via Laplace approximations for (1) the marginal likelihood with respect to random effects and (2) the restricted likelihood (as in REML), i.e. the likelihood of random effect parameters given the fixed effect estimates. 
  
     All handled models can be formulated in terms of a linear predictor of the traditional form \code{offset}+ \bold{X}\eqn{\beta} + \bold{Z b}, where \bold{X} is the design matrix of fixed effects, \eqn{\beta} (\code{beta}) is a vector of fixed-effect coefficients, \bold{Z} is a \dQuote{design matrix} for the random effects (which is instead denoted \bold{M}=\bold{ZAL} elsewhere in the package documentation), and \bold{b} a vector of random effect values. The general structure of \bold{Mb} is described in \code{\link{random-effects}}.
     
  Gaussian and non-gaussian random effects can be fitted. Different \strong{gaussian} random-effect terms are handled, with the following effects:\cr 
\preformatted{
* (1|<RHS>), for non-autocorrelated random effects as in lme4;
* (<LHS>|<RHS>), for random-coefficient terms as in lme4, *and 
   additional terms depending on the <LHS> type* (further detailed below);
* (<LHS> || <RHS>) is interpreted as in lme4: any such term is immediately 
   converted to ( (1|<RHS>) + (0+<LHS>|<RHS>) ). It should be counted as two 
   random effects for all purposes (e.g., for fixing the variances of the 
   random effects). However, this syntax is useless when the LHS includes a 
   factor (see help('lme4::expandDoubleVerts')).
* <prefix>(1|<RHS>), to specify autocorrelated random effects, 
   e.g. Matern(1|long+lat). 
* <prefix>(<LHS>|<RHS>), where the <LHS> can be used to alter the 
   autocorrelated random effect as detailed below. 
}   

Different LHS types of \strong{gaussian} \code{(<LHS>|<RHS>)} random-effect terms are handled, with the following effects:
\preformatted{
* <logical> (TRUE/FALSE): affects only responses for which <LHS> is TRUE; % equiv O + as.numeric(<logical>)
* <factor built from a logical>: same a <logical> case;
* <factor not built from a logical>: random-coefficient term as in lme4;
* 0 + <factor not built from a logical>: same but contrasts are not used;
* factors specified by the mv(...) expression, generate random-coefficient 
  terms specific to multivariate-response models fitted by fitmv() (see 
  help("mv")). 0 + mv(...) has the expected effect of not using contrasts; 
* <numeric> (but not '0+<numeric>'): random-coefficient term as in lme4, 
  with 2*2 covariance matrix of effects on Intercept and slope;
* 0 + <numeric>: no Intercept so no covariance matrix (random-slope-only 
   term);
* the "experimental" dummy() specifier described in the lme4 documentation
  has been found to work with spaMM too. 
}
The '0 + <numeric>' effect is achieved by direct control of the elements of the incidence matrix \bold{Z} through the \code{<LHS>} term: for numeric \code{z}, such elements are multiplied by \code{z} values, and thus provide a variance of order O(\code{z} \bold{squared}).

If one wishes to fit uncorrelated group-specific random-effects with distinct variances for different groups or for different response variables, three syntaxes are thus possible. The most general, suitable for fitting several variances (see {help("GxE")} for an example), is to fit a (0 + <factor>| <RHS>) random-coefficient term with correlation(s) fixed to 0. Alternatively, one can define \bold{numeric} (0|1) variables for each group (as \code{as.numeric(<boolean for given group membership>)}), and use each of them in a \verb{0 + <numeric>} LHS (so that the variance of each such random effect is zero for response not belonging to the given group). Using \code{dummy} factors for each level of the original factor can achieve the same effect without coding new variables in the data. These syntaxes are illustrated in the \dQuote{Elementary multivariate-response} example below for two groups (males and females). 

The following \strong{gaussian} \code{<prefix>(<LHS>|<RHS>)} random-effect terms are handled: \cr 
\verb{   } * Since version 3.9.0 it is possible to fit \code{corrMatrix(<LHS>|<RHS>)} terms that combines features of any of the above LHS types and of a random effect with given \code{corrMatrix}. See \code{\link{composite-ranef}} for details.\cr 
\verb{   } * It is possible to fit some random-effect terms involving other (\code{<prefix>}) correlation models and a non-\sQuote{\code{1}} LHS. Only the LHS types that do not imply estimating a random-coefficient covariance matrix are currently handled. For example, independent Matérn effects can be fitted for males and females by using the syntax \code{Matern(male|.) + Matern(female|.)}, where \code{male} and \code{female} are TRUE/FALSE factors; and \code{Matern(0+<numeric>|.)} represents an autocorrelated random-slope (only) term (or, equivalently, a direct specification of heteroscedasticity of the Matérn random effect). By contrast, \code{Matern(z|.)} (with effect similar to that described in \code{\link{composite-ranef}}) is not yet implemented.  
\verb{   } * In these autocorrelated random effects the \code{<RHS>} is generally restricted to the simplest form, except for AR1 and geostatistical models where the \code{<.> \%in\% <.>} form of nested random effect is allowed (see example in \code{\link{MaternCorr}}). 
% cf fns .calc_Z_model_matrix() and .calc_raw_ZA().
% tests: test-spaMM.R and in tests_private/test-syntax-Matern-RHS.R

The syntax \code{(z-1|.)}, for \bold{numeric} \code{z} only, can also be used to fit \strong{some heteroscedastic non-Gaussian} random effects. For example, a Gamma random-effect term \code{(wei-1|block)} specifies an heteroscedastic Gamma random effect \eqn{u} with constant mean 1 and variance \code{wei^2} \eqn{\lambda}, where \eqn{\lambda} is still the estimated variance parameter. See Details of \code{\link{negbin}} for a possible application. Here, this effect is not implemented through direct control of \bold{Z} (multiplying the elements of an incidence matrix \bold{Z} by \code{wei}), as this would have a different effect on the distribution of the random effect term. \code{(z|.)} is not defined for \emph{non-Gaussian} random effects. It could mean that a correlation structure between random intercepts and random slopes for (say) Gamma-distributed random effects is considered, but such correlation structures are not well-specified by their correlation matrix. 

}
\author{
\code{spaMM} was initially published by \enc{François}{Francois} Rousset and Jean-Baptiste Ferdy, and is continually developed by F. Rousset and tested by Alexandre Courtiol.  
}
\references{
Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalized linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.

Rousset F., Ferdy, J.-B. (2014) Testing environmental and genetic effects in the presence of spatial autocorrelation. Ecography, 37: 781-790.
    \doi{10.1111/ecog.00566}

}
\seealso{
The \code{test} directory of the package provides many additional examples of spaMM usage beyond those from the formal documentation. 
}

\examples{
data("wafers")
data("scotlip") ## loads 'scotlip' data frame, but also 'Nmatrix'

##     Linear model
fitme(y ~ X1, data=wafers)

##     GLM
fitme(y ~ X1, family=Gamma(log), data=wafers)
fitme(cases ~ I(log(population)), data=scotlip, family=poisson)

##     Non-spatial GLMMs
fitme(y ~ 1+(1|batch), family=Gamma(log), data=wafers)
fitme(cases ~ 1+(1|gridcode), data=scotlip, family=poisson)
#
# Random-slope model (mind the output!)        
fitme(y~X1+(X2|batch),data=wafers, method="REML")

## Spatial, conditional-autoregressive GLMM
if (spaMM.getOption("example_maxtime")>2) {   
  fitme(cases ~ I(log(population))+adjacency(1|gridcode), data=scotlip, family=poisson, 
        adjMatrix=Nmatrix) # with adjacency matrix provided by data("scotlip")
} 
# see ?adjacency for more details on these models 

## Spatial, geostatistical GLMM: 
# see e.g. examples in ?fitme, ?corrHLfit, ?Loaloa, or ?arabidopsis;
# see examples in ?Matern for group-specific spatial effects.

##     Hierachical GLMs with non-gaussian random effects
 data("salamander")
if (spaMM.getOption("example_maxtime")>1) {   
 # both gaussian and non-gaussian random effects
 fitme(cbind(Mate,1-Mate)~1+(1|Female)+(1|Male),family=binomial(),
        rand.family=list(gaussian(),Beta(logit)),data=salamander)
 
 # Random effect of Male nested in that of Female:
 fitme(cbind(Mate,1-Mate)~1+(1|Female/Male),
       family=binomial(),rand.family=Beta(logit),data=salamander)
 # [ also allowed is cbind(Mate,1-Mate)~1+(1|Female)+(1|Male \%in\% Female) ]
}

##    Modelling residual variance ( = structured-dispersion models)    
# GLM response, fixed effects for residual variance 
fitme( y ~ 1,family=Gamma(log),
      resid.model = ~ X3+I(X3^2) ,data=wafers)
#
# GLMM response, and mixed effects for residual variance
if (spaMM.getOption("example_maxtime")>1.5) {   
  fitme(y ~ 1+(1|batch),family=Gamma(log),
        resid.model = ~ 1+(1|batch) ,data=wafers)
}

## Elementary multivariate-response model (see fitmv() for a more general procedure) 
# Data preparation
fam <- rep(c(1,2),rep(6,2)) # define two biological 'families'
ID <- gl(6,2) # define 6 'individuals'
resp <- as.factor(rep(c("x","y"),6)) # distinguishes two responses per individual
set.seed(123)
toymv <- data.frame(
  fam = factor(fam), ID = ID, resp = resp, 
  y = 1 + (resp=="x") + rnorm(4)[2*(resp=="x")+fam] + rnorm(12)[6*(resp=="x")+as.integer(ID)],
  respX = as.numeric(resp=="x"),
  respY = as.numeric(resp=="y")
)
#
# fit response-specific variances of random effect and residuals: 
(fitme(y ~ resp+ (0+respX|fam)+ (0+respY|fam), resid.model = ~ 0+resp ,data=toymv))
#
# Same result by different syntaxes:
#
(fitme(y ~ resp+ (0+resp|fam), resid.model = ~ 0+resp ,data=toymv, 
       fixed=list(ranCoefs=list("1"=c(NA,0,NA)))))

# or by the dummy() specifier from lme4:
# (fitme(y ~ resp+ (0+dummy(resp,"x")|fam)+ (0+dummy(resp,"y")|fam), 
#        resid.model = ~ 0+resp ,data=toymv))
}
\keyword{ package }
\keyword{spatial}
% should add \keyword{smooth} later... and htest ? models ? see RShowDoc("KEYWORDS")
