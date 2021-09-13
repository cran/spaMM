\name{composite-ranef}
\alias{composite-ranef}
\title{
Composite random effects
}
\description{
An example of a composite random effect is \code{corrMatrix(sex|pair)}. It combines features of a random-coefficient model \code{(sex|pair)} and of a random effect \code{corrMatrix(1|pair)}. The random-coefficient model is  characterized by a \eqn{2*2} covariance matrix \bold{C} for the random effects \eqn{u_{1,pair}} and \eqn{u_{2,pair}} both affecting each of the two sexes for each \code{pair}, and the \code{corrMatrix} random effect assumes that elements of each of the two vectors \eqn{u_i=(u_{i,pair})} for \code{pair}=1,...,\eqn{P} are correlated according to a given \eqn{P*P} correlation matrix \bold{A}. Then the composite random effect is defined as the one with 2\eqn{P*}2\eqn{P} covariance matrix \code{\link{kronecker}}(\bold{C},\bold{A}).

Composite random effects can also be fitted for multivariate-response models, e.g. \code{corrMatrix(mv(1,2)|ID)} for two responses variables for each individual \code{ID}. 

The definition of composite random effects through the \code{kronecker} product may be motivated and understood in light of a quantitative-genetics application (see help("Gryphon") for an example). In this context the two response variables are typically two individual traits. Each trait is affected by two sets of genes, the effect of each set being represented by a gaussian random effect (\code{u_1} or \code{u_2}). The effect of genetic relatedness on the correlation of random effects \code{u_i,ID} among individuals \code{ID} within each set \eqn{i} of genes is described by the corrMatrix \bold{A}. The effects on the two traits for each individual are interpreted as different linear combinations of these two random effects (the coefficients of these linear combinations determining the \bold{C} matrix). Under these assumptions the correlation matrix of the responses (in order (trait, individual)=(1,1)...(1,ID)... (2,1)...(2,ID)...) is indeed \code{\link{kronecker}}(\bold{C},\bold{A}).

The summary of the model provides a description of the \bold{C} matrix in terms of its variances and its correlation coefficient(s)) when \bold{C} is viewed as a covariance matrix. In a standard random-coefficient model these variances are those of the correlated  random effects (see \code{\link{summary.HLfit}}). In the composite random-effect model this is not necessarily so as the variance of the correlated random effects also depend on the variances implied by the \bold{A} matrix, which are not necessarily 1 if \bold{A} is a covariance matrix rather than simply a correlation matrix.    

In a \code{corrMatrix(<LHS>|<RHS>)} term the type (logical, factor...) of the <LHS> has an effect identical to its effect in a non-composite \code{(<LHS>|<RHS>)} term, as described in \code{\link{spaMM}}. In particular, in some cases no random-coefficient correlation matrix \bold{C} is implied: see the Examples below.
}
\examples{
if (spaMM.getOption("example_maxtime")>1.8) {
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


}
}
