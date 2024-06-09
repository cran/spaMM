\name{composite-ranef}
\alias{composite-ranef}
\title{
Composite random effects
}
\description{
An example of a composite random effect is \code{corrMatrix(sex|pair)}. It combines features of a random-coefficient model \code{(sex|pair)} and of an autocorrelated random effect \code{corrMatrix(.|.)}. The random-coefficient model is characterized by a \eqn{2*2} covariance matrix \bold{C} for the random effects \eqn{u_{1,pair}} and \eqn{u_{2,pair}} both affecting each of the two sexes for each \code{pair}, and the \code{corrMatrix} random effect assumes that elements of each of the two vectors \eqn{u_i=(u_{i,pair})} for \code{pair}=1,...,\eqn{P} are correlated according to a given \eqn{P*P} correlation matrix \bold{A}. Then the composite random effect is defined as the one with 2\eqn{P*}2\eqn{P} covariance matrix \code{\link{kronecker}}(\bold{C},\bold{A}).

The definition of composite random effects through the \code{kronecker} product may be motivated and understood in light of a quantitative-genetics application (see help("Gryphon") for an example). In this context the two response variables are typically two individual traits. Each trait is affected by two sets of genes, the effect of each set being represented by a gaussian random effect (\code{u_1} or \code{u_2}). The effect of genetic relatedness on the correlation of random effects \code{u_i,ID} among individuals \code{ID} within each set \eqn{i} of genes is described by the corrMatrix \bold{A}. The effects on the two traits for each individual are interpreted as different linear combinations of these two random effects (the coefficients of these linear combinations determining the \bold{C} matrix). Under these assumptions the correlation matrix of the responses (in order (trait, individual)=(1,1)...(1,ID)... (2,1)...(2,ID)...) is indeed \code{\link{kronecker}}(\bold{C},\bold{A}).

Composite random effects are not restricted to \code{corrMatrix} terms, and can also be fitted for multivariate-response models. For example, \code{Matern(mv(1,2)|longitude+latitude)} terms can be fitted, in which case the correlation model is still defined through the Kronecker product, where \bold{A} will be a (fitted) Matérn correlation matrix, and \bold{C} will be the correlation matrix of the fitted random-coefficient model for the \code{\link{mv}} virtual factor for multivariate response. 

Implementation and tests of composite random effects is work in progress, with the following ones having been tested: \code{corrMatrix}, \code{Matern}, \code{Cauchy}, \code{adjacency}, \code{IMRF}, \code{AR1}, and to a lesser extent \code{MaternIMRFa}. Fits of other composite terms may fail. Even if they succeed, not all post-fit procedures may be operational: in particular, prediction (and then, simulation) with \code{newdata} may not yet work. 
% * Tests of composite non-mv
%   In test-composite-extra.R:
%     tests of point predictions for AR1, ARp, Matern, corrMatrix, and less strict for IMRF (but see below) and MaternIMRFa
%     tests of get_predVar for AR1, ARp, at least
%     Warning+stop for missing code for corrFamily ? 
% * Tests of composite mv
%   In test-composite-extra.R:
%     test of corrMatrix(0+mv(...)|...)
%   In test-mv-extra.R:
%     test of adjacency(0+mv(...)|...)
%     test of IMRF(0+mv(...)|...)
%
%  predVar inaccurate as stated: 
Further, as for random-coefficient terms in univariate-response models, some components of the computed prediction variance depend on a poorly characterized approximation, yielding different results for different fitting algorithms (see Details in \code{\link{predVar}}).

The summary of the model provides fitted parameters for \bold{A} if this matrix derives from a parametric correlation model (e.g., \code{Matern}), and a description of the \bold{C} matrix where it is viewed as a covariance matrix, parametrized by its variances and its correlation coefficient(s). In a standard random-coefficient model these variances would be those of the correlated  random effects (see \code{\link{summary.HLfit}}). In the composite random-effect model this is not necessarily so as the variance of the correlated random effects also depend on the variances implied by the \bold{A} matrix, which are not necessarily 1 if \bold{A} is a covariance matrix rather than simply a correlation matrix.    

A \code{<prefix>(<LHS>|<RHS>)} term is *not* a composite random effect when the LHS in boolean, a factor from boolean, or \dQuote{0+numeric}. See the Matérn examples, and the corrMatrix \dQuote{\code{<LHS> is 0+numeric}} example, below.
}
\examples{
if (spaMM.getOption("example_maxtime")>1.8) {

########################
#### corrMatrix examples
########################

## Toy data preparation

data("blackcap")
toy <- blackcap
toy$ID <- gl(7,2)
grp <- rep(1:2,7)
toy$migStatus <- toy$migStatus +(grp==2)
toy$loc <- rownames(toy) # to use as levels matching the corrMatrix dimnames

toy$grp <- factor(grp)
toy$bool <- toy$grp==1L
toy$boolfac <- factor(toy$bool)
toy$num <- seq(from=1, to=2, length.out=14)

## Build a toy corrMatrix as perturbation of identity matrix:
n_rhs <- 14L
eps <- 0.1
set.seed(123)
rcov <- ((1-eps)*diag(n_rhs)+eps*rWishart(1,n_rhs,diag(n_rhs)/n_rhs)[,,1])
# eigen(rcov)$values
colnames(rcov) <- rownames(rcov) <- toy$loc # DON'T FORGET NAMES


##### Illustrating the different LHS types

### <LHS> is logical (TRUE/FALSE) => No induced random-coefficient C matrix; 
#   corrMatrix affects only responses for which <LHS> is TRUE:
#
(fit1 <- fitme(migStatus ~ bool + corrMatrix(bool|loc), data=toy, corrMatrix=rcov))
#
# Matrix::image(get_ZALMatrix(fit1))


### <RHS> is a factor built from a logical => same a 'logical' case above:
#
(fit2 <- fitme(migStatus ~ boolfac + corrMatrix(boolfac|loc), data=toy, corrMatrix=rcov))
#
# Matrix::image(get_ZALMatrix(fit2))


### <RHS> is a factor not built from a logical: 
# (grp|.) and (0+grp|.) lead to equivalent fits of the same composite model, 
#   but contrasts are not used in the second case and the C matrices differ,
#   as for standard random-coefficient models.
#
(fit1 <- fitme(migStatus ~ grp +  corrMatrix(grp|loc), data=toy, corrMatrix=rcov))
(fit2 <- fitme(migStatus ~ grp +  corrMatrix(0+grp|loc), data=toy, corrMatrix=rcov))
# 
# => same fits, but different internal structures:
Matrix::image(fit1$ZAlist[[1]]) # (contrasts used) 
Matrix::image(fit2$ZAlist[[1]]) # (contrasts not used)
# Also compare ranef(fit1) versus ranef(fit2) 
#
#
## One can fix the C matrix, as for standard random-coefficient terms 
#
(fit1 <- fitme(migStatus ~ grp +  corrMatrix(0+grp|loc),data=toy, corrMatrix=rcov, 
               fixed=list(ranCoefs=list("1"=c(1,0.5,1)))))
#       
# same result without contrasts hence different 'ranCoefs':             
#
(fit2 <- fitme(migStatus ~ grp +  corrMatrix(grp|loc), data=toy, corrMatrix=rcov, 
               fixed=list(ranCoefs=list("1"=c(1,-0.5,1)))))


### <LHS> is numeric (but not '0+numeric'):
# composite model with C being 2*2 for Intercept and numeric variable
#
(fitme(migStatus ~ num +  corrMatrix(num|loc), data=toy, corrMatrix=rcov))

### <LHS> is 0+numeric: no random-coefficient C matrix 
#  as the Intercept is removed, but the correlated random effects 
#  arising from the corrMatrix are multiplied by sqrt(<numeric variable>)
#
(fitme(migStatus ~ num +  corrMatrix(0+num|loc), data=toy, corrMatrix=rcov))


### <LHS> for multivariate response (see help("Gryphon") for more typical example)
## More toy data preparation for multivariate response
ch <- chol(rcov)
set.seed(123)
v1 <- tcrossprod(ch,t(rnorm(14,sd=1)))
v2 <- tcrossprod(ch,t(rnorm(14,sd=1)))
toy$status <- 2*v1+v2
toy$status2 <- 2*v1-v2

## Fit:
fitmv(submodels=list(mod1=list(status ~ 1+ corrMatrix(0+mv(1,2)|loc)),
                     mod2=list(status2 ~ 1+ corrMatrix(0+mv(1,2)|loc))), 
      data=toy, corrMatrix=rcov)
      
##################################################
#### Matern examples: sex-dependent random effects
##################################################

if (spaMM.getOption("example_maxtime")>2) {
data(Leuca)
subLeuca <- Leuca[c(1:10,79:88),] # subset for faster examples

# The random effects in the following examples are composite because 'sex' is not 
# boolean nor factor from boolean. If 'Matern(factor(female)|x+y)' were used, the effect 
# would be the same 'Matern(female|x)', fitting 

  fitme(fec_div ~ 1 + Matern(sex|x+y),data=subLeuca)  # => builds a random effect
# correlated across sexes, from 2 independent realizations u_1 and u_2 of 20 values 
# (for the 20 locations in the data). In a (sex|x) term the 20 values would be 
# independent from each other within each u_i. In the Matern(sex|x+y) such 20 values 
# are autocorrelated within each u_i.

# For pedagogic purposes, one can also fit
  fitme(fec_div ~ 1 + Matern(sex|x + y \%in\% sex),data=subLeuca)
# which again builds a random effect from 2 independent realizations 
# u_1 and u_2, but each u_i now contains two realizations u_i1 and u_i2 of 10 values,
# autocorrelated within each u_ij following the Matern model, 
# but independent between u_i1 and u_i2. As a result, the overall random effect in each sex, 
# v_m or v_f, is a weighted sum of two sex-specific Matern random effect, 
# so that v_m and v_f are independent from each other. 
}

}
}
