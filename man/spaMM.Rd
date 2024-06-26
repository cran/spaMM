\name{spaMM}
% ----- \Rdversion{1.1}
\alias{spaMM}
\alias{spaMM-package}

\docType{package}

\title{Inference in mixed models, in particular spatial GLMMs}

\description{Fits a range of mixed-effect models, including those with spatially correlated random effects. The random effects are either Gaussian (which defines GLMMs), or other distributions (which defines the wider class of hierarchical GLMs), or simply absent (which makes a LM or GLM).  Multivariate-response models can be fitted by the \code{\link{fitmv}} function. Other models can be fitted by \code{\link{fitme}}. Also available are previously conceived fitting functions \code{\link{HLfit}} (sometimes faster, for non-spatial models), \code{\link{HLCor}} (sometimes faster, for conditional-autoregressive models and fixed-correlation models), and \code{\link{corrHLfit}} (now of lesser interest). A variety of post-fit procedures are available for prediction, simulation and testing (see, e.g., \code{\link{fixedLRT}}, \code{\link{simulate}} and \code{\link{predict}}).

A variety of special syntaxes for fixed effects, such as \code{\link{poly}}, \code{splines::}\code{\link[splines]{ns}} or \code{\link[splines]{bs}}, or \code{lmDiallel::GCA}, may be handled natively although some might not be fully handled by post-fit procedures such as \code{predict}. \code{poly} is fully handled. \code{lmDiallel::GCA} cannot be due to its inherent limitations, but see \code{\link{X.GCA}} for a more functional alternative for diallel/multi-membership fixed-effect terms. Note that packages implementing these terms must be attached to the search list as \code{::} will not be properly understood in a \code{formula}. % ___F I X M E____ provide informative messages ?   

Both maximum likelihood (ML) and restricted likelihood (REML) can be used for linear mixed models, and extensions of these methods using Laplace approximations are used for non-Gaussian random response. Several variants of these methods discussed in the literature are included (see Details in \code{\link{HLfit}}), the most notable of which may be \dQuote{PQL/L} for binary-response GLMMs (see Example for \code{\link{arabidopsis}} data). PQL methods implemented in spaMM are closer to (RE)ML methods than those implemented in \code{MASS::glmmPQL}. 
}

\details{
  The standard response families \code{gaussian}, \code{binomial}, \code{poisson}, and \code{Gamma} are handled, as well as negative binomial (see \code{\link{negbin1}} and \code{\link{negbin2}}), beta (\code{\link{beta_resp}}), beta-binomial (\code{\link{betabin}}), zero-truncated poisson and negative binomial and Conway-Maxwell-Poisson response (see \code{\link{Tpoisson}}, \code{\link{Tnegbin}} and \code{\link{COMPoisson}}). A \code{multi} family look-alike is also available for \code{\link{multinomial}} response, with some constraints. 
  
  The variance parameter of residual error is denoted \eqn{\phi} (\code{phi}): this is the residual variance for gaussian response, but for Gamma-distributed response, the residual variance is \eqn{\phi}\eqn{\mu^2} where \eqn{\mu} is expected response. A (possibly mixed-effects) linear predictor for \eqn{\phi}, modeling heteroscedasticity, can be considered (see Examples). 
  
  The package fits models including several nested or crossed random effects, including autocorrelated ones. An interface is being developed allowing users to implement their own parametric correlation models (see \code{\link{corrFamily}}), beyond the following ones which are built in \pkg{spaMM}:\cr 
\verb{ } * geostatistical (\code{\link{Matern}}, \code{\link{Cauchy}}),\cr 
\verb{ } * interpolated Markov Random Fields (\code{\link{IMRF}}, \code{\link{MaternIMRFa}}),\cr 
\verb{ } * autoregressive time-series (\code{\link{AR1}}, \code{\link{ARp}}, \code{\link{ARMA}}),\cr 
\verb{ } * conditional autoregressive as specified by an \code{\link{adjacency}} matrix,\cr 
\verb{ } * pairwise interactions with individual-level random effects, such as diallel experiments (\code{\link{diallel}}),\cr 
\verb{ } * or any fixed correlation matrix (\code{\link{corrMatrix}}). 
  
  GLMMs and HGLMs are fit via Laplace approximations for (1) the marginal likelihood with respect to random effects and (2) the restricted likelihood (as in REML), i.e. the likelihood of random effect parameters given the fixed effect estimates. All handled models can be formulated in terms of a linear predictor of the traditional form \code{offset}+ \bold{X}\eqn{\beta} + \bold{Z b}, where \bold{X} is the design matrix of fixed effects, \eqn{\beta} (\code{beta}) is a vector of fixed-effect coefficients, \bold{Z} is a \dQuote{design matrix} for the random effects (which is instead denoted \bold{M}=\bold{ZAL} elsewhere in the package documentation), and \bold{b} a vector of random effect values. The general structure of \bold{Mb} is described in \code{\link{random-effects}}.
     
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
* <logical> (TRUE/FALSE): affects only responses for which <LHS> is TRUE. 
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
}
The '0 + <numeric>' effect is achieved by direct control of the elements of the incidence matrix \bold{Z} through the \code{<LHS>} term: for numeric \code{z}, such elements are multiplied by \code{z} values, and thus provide a variance of order O(\code{z} \bold{squared}).

If one wishes to fit uncorrelated group-specific random-effects with distinct variances for different groups or for different response variables, three syntaxes are thus possible. The most general, suitable for fitting several variances (see \code{\link{GxE}} for an example), is to fit a (0 + <factor>| <RHS>) random-coefficient term with correlation(s) fixed to 0. Alternatively, one can define \bold{numeric} (0|1) variables for each group (as \code{as.numeric(<boolean for given group membership>)}), and use each of them in a \verb{0 + <numeric>} LHS (so that the variance of each such random effect is zero for response not belonging to the given group). See \code{\link{lev2bool}} for various ways of specifying such indicator variables for several levels.

\strong{Gaussian} \code{<prefix>(<LHS not 1>|<RHS>)} random-effect terms may be handled, with two main cases depending on the LHS type, motivated by the following example: independent Matérn effects can be fitted for males and females by using the syntax \code{Matern(male|.) + Matern(female|.)}, where \code{male} and \code{female} are TRUE/FALSE (or a factor with TRUE/FALSE levels). In contrast to a \code{(male|.)} term, no random-coefficient correlation matrix is fitted. However, for some other types of RHS, one can fit \emph{composite random effects} combining a random-coefficient correlation matrix and the correlation model defined by the \dQuote{prefix}. This combination is defined in \code{\link{composite-ranef}}. This leads to the following distinction:\cr 
\verb{   } * The terms are *not* composite random effects when the non-\sQuote{\code{1}} LHS type is boolean or factor-from-boolean, a just illustrated, but also \code{0+<numeric>}: for example, \code{Matern(0+<numeric>|.)} represents an autocorrelated random-slope (only) term or, equivalently, a direct specification of heteroscedasticity of the Matérn random effect.\cr
% cf fns .calc_Z_LHS_model_matrix() and .calc_raw_ZA().
% tests: test-spaMM.R and in tests_private/test-syntax-Matern-RHS.R; MaternCorr.Rd
\verb{   } * By contrast, \code{Matern(<numeric>|.)} implies estimating a random-coefficient covariance matrix and thus defines a composite random effects, as does an LHS that is a factor constructed from numeric or character levels.\cr 
Composite random effects can be fitted in principle for all \dQuote{prefixes}, including for \code{<\link{corrFamily}>} terms. In practice, this functionality has been checked for \code{Matern}, \code{corrMatrix}, \code{AR1} and the \code{ARp}-corrFamily term. In these terms, the \code{<.> \%in\% <.>} form of nested random effect is allowed. 
% test-composite, and much more in test-composite-extra.R incl. a lot of mv models.
% tested up to Matern(grp|latitude %in% grp) and corrMatrix([idem]) in test-composite-extra, composite-ranef.Rd, test-syntax-Matern-RHS.R.
% For the non-composite random effects, the \code{<.> \%in\% <.>} is tested by AR1(0+age|age %in% Subject) in test-AR1.R,
% and by Matern(sex|latitude %in% sex) in test-syntax-Matern-RHS.R.

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
See the \code{test} directory of the package for many additional examples of \pkg{spaMM} usage beyond those from the formal documentation. 

See \code{\link{fitme}} for multivariate-response models. 
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

}
\keyword{ package }
\keyword{spatial}
% should add \keyword{smooth} later... and htest ? models ? see RShowDoc("KEYWORDS")
