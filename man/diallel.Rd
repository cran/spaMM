\name{diallel}
\alias{ranGCA}
\alias{diallel}
\title{
Random-effect structures for diallel experiments and other symmetric dyadic interactions
}
\description{
In diallel experiments, the phenotypes of offspring from crosses where the same parent can play the role of mother as well as that of father is analyzed. Then one must find a way to fit a model accounting for the fact that the same individual (with same genetic effect) appears in two variables (mother genotype and father genotype). 

This concern arises more generally in the case where interactions occur between pairs of individuals, where one analyzes an outcome for each interaction (here one response value for each pairwise interaction, rather than two individual response values). If individual-level random effects of the form (1|ID1)+ (1|ID2) were included in the model formula, this would result in different variances being fitted for each random effect (breaking the assumption of symmetry), and the value of the random effect would differ for an individual whether it appears as a level of the first random effect or of the second (which is also inconsistent with the idea that the random effect represents a property of the individual). Instead one should include a \code{ranGCA(1|ID1+ID2)} term in the model formula.

Although this random-effect structure may be used in many different contexts (including social network analysis, or \dQuote{round robin} experiments in psychology), its present name refers to the semantics established for diallel experiments (e.g., Lynch & Walsh, 1998, p. 611), because it is not easy to find a more general yet intuitive semantics. In this context the additive genetic effects of each parent's genotypes are described as \dQuote{general combining abilities} (GCA). In case of non-additive effect, the half-sib covariance is not half the full-sib covariance and this is represented by another effect described as \dQuote{specific combining abilities} (SCA) . More generally, think about GCAs as random effects associated with each member of a dyad and having additive effects on the linear predictor, and about GCAs + SCAs as a way of representing non-additive effects of each member.    

\verb{ * }The \code{ranGCA} \code{\link{corrFamily}} constructor defines a random effect which is the sum of the random effects contributed by each parent.\cr
\verb{ * }The \code{diallel} corrFamily constructor defines a random effect with a one-parameter covariance structure that can represent GCAs and SCAs, and more generally non-additive effects of each individual. This extended model is formulated as an autocorrelated random effect, which is related to the additive case as follows. A level of the random effect is assigned for each unordered pair of individuals, and a correlation \eqn{\rho} is assigned between levels for pairs sharing one individual. In the context of a diallel experiment this amounts to define a synthetic \dQuote{mid-parent} random effect \dQuote{received} by the offspring, with distinct levels for each unordered parental pair, and this assigns a correlation \eqn{\rho} between effects received by half-sibs (one shared parent). \eqn{\rho} corresponds to var(GCA)/[2*var(GCA)+var(SCA)] and is necessarily \eqn{\le 0.5}. 
The \code{diallel} model should actually be fitted for \eqn{\rho < 0.5} only; \code{ranGCA} fits the exact case \eqn{\rho= 0.5}.
%\eqn{\rho=0.5} actually generates a singular correlation matrix and the fit may then fail, or succeed with warnings, depending on the exact algorithm used internally by \pkg{spaMM}, 
% ... depending on .calc_denseness(cov_info_mat,relative=TRUE) < 0.15  => CHMfactor computed or not.  
%... Indeed, with \eqn{n} individuals, the distinct values of the random effect for up to \eqn{n(n-1)/2} pairs are all linearly related to the \eqn{n} individual values when \eqn{\rho=0.5}, so the rank of the correlation matrix is at most \eqn{n} even if it as more levels (more pairs).  

\code{diallel} fits can be slow for large data if the correlation matrix is large, as it can have a fair proportion of nonzero elements. 
There may also be identifiability issues for variance parameters: in a LMM as shown in the examples, there will be three parameters for the random variation (\code{phi}, \code{lambda} and \code{rho}) but only two can be estimated if only one observation is made for each dyad.  

GCAs and SCAs can also be fitted as fixed effects. \pkg{spaMM} has no specific syntax for that purpose, but it seems that the fixed-effect terms defined in the \pkg{lmDiallel} package (Onofri & Terzaroli, 2021) work in a formula for a \pkg{spaMM} fit. 

If the names \code{ranGCA} and \code{diallel} sound inappropriate for your context of application, you can register these terms under a new name using \code{\link{register_cF}}.
}
\usage{
## formula terms:
# ranGCA(1| <.> + <.>) 
# diallel(1|  <.> + <.>, tpar, fixed = NULL, public = NULL)
## where the <.> are two factor identifiers, ** whose levels
## must be identical when representing the same individual **

## corrFamily constructors:
ranGCA() # no argument
diallel(tpar, fixed = NULL, public = NULL)
}
\arguments{
  \item{tpar}{
  Numeric: template value of the correlation coefficient for pairs sharing one individual. 
}
  \item{fixed}{
NULL or fixed value of the correlation coefficient.
}
\item{public}{NULL, or an environment. When an empty environment is provided, a template \code{CorNA} for the correlation matrix (with NA's in place of \eqn{rho}) will be copied therein, for inspection at user level.}
}
\value{
Both functions return corrFamily descriptors whose general format is described in \code{\link{corrFamily}}. The one produced by \code{ranGCA} is atypical in that only its \code{Af} element is non-trivial.
}
\references{
Lynch, M., Walsh, B. (1998) Genetics and analysis of quantitative traits. Sinauer, Sunderland, Mass.

Onofri A., Terzaroli N. (2021) lmDiallel:  Linear Fixed Effects Models for Diallel Crosses. Version 0.9.4. \url{https://cran.r-project.org/package=lmDiallel}.
}
\examples{

#### Simulate dyadic data

set.seed(123)
nind <- 10
x <- runif(nind^2) 
id12 <- expand.grid(id1=seq(nind),id2=seq(nind))
id1 <- id12$id1
id2 <- id12$id2
u <-  rnorm(50,mean = 0, sd=0.5)

## additive individual effects:
y <-  0.1 + 1*x + u[id1] +  u[id2] + rnorm(nind^2,sd=0.2)

## Same with non-additive individual effects:
dist.u <- abs(u[id1] -  u[id2])
z <- 0.1 + 1*x + dist.u + rnorm(nind^2,sd=0.2)

dyaddf <- data.frame(x=x, y=y, z=z, id1=id1,id2=id2)
# : note that this contains two rows per dyad, which avoids identifiability issues.

# Enforce that interactions are between distinct individuals (not essential for the fit):
dyaddf <- dyaddf[- seq.int(1L,nind^2,nind+1L),] 


# Fits:
#
(addfit <- fitme(y ~x +ranGCA(1|id1+id2), data=dyaddf))
#
# practically equivalent to:
#
(fitme(y ~x +diallel(1|id1+id2, fixed=0.49999), data=dyaddf))

(distfit <- fitme(z ~x +diallel(1|id1+id2), data=dyaddf)) 
     
}
