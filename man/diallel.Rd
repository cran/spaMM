\name{diallel}
\alias{ranGCA}
\alias{antisym}
\alias{diallel}
\title{
Random-effect structures for diallel experiments and other dyadic interactions
}
\description{
\code{ranGCA} and \code{diallel} are random-effect structures designed to represent the effet of symmetric interactions between pairs of individuals (order of individuals in the pair does not matter), while \code{antisym} represents anti-symmetric interactions (the effect of reciprocal ordered pairs on the outcome are opposite, as in the so-called Bradley-Terry models). These random-effect structures all account for multiple membership, i.e., the fact that the same individual may act as the first or the second individual among different pairs, or even within one pair if this makes sense).   

More formally, the outcome of an interaction between a pair \eqn{i,j} of agents is subject to a symmetric overall random effect \eqn{v_{ij}} when the effect  \dQuote{on} individual \eqn{i} (or viewed from the perspective of individual \eqn{i}) equals the effect on \eqn{j}: \eqn{v_{ij}=v_{ji}}. This may result from the additive effect of individual random effects \eqn{v_{i}} and \eqn{v_{j}}:  \eqn{v_{ij}=v_i+v_j}, but also from non-additive effects \eqn{v_{ij}=v_i+v_j+a_{ij}} if the interaction term  \eqn{a_{ij}} is itself symmetric (\eqn{a_{ij}=a_{ji}}). \code{ranGCA} and \code{diallel} effects represent such symmetric effects, additive or non-additive respectively, in a model formula (see Details for the semantic origin of these names and how they can be changed). Conversely, antisymmetry is characterized by \eqn{v_{ij}=v_i-v_j=-v_{ji}} and is represented by the \code{antisym} formula term.  

If individual-level random effects of the form (1|ID1)+ (1|ID2) were included in the model formula instead of \code{ranGCA(1|ID1+ID2)} for symmetric additive interactions, this would result in different variances being fitted for each random effect (breaking the assumption of symmetry), and the value of the random effect would differ for an individual whether it appears as a level of the first random effect or of the second (which is also inconsistent with the idea that the random effect represents a property of the individual). 

When \code{ranGCA} or \code{antisym} random effects are fitted, the individual effects are inferred. By contrast, when a \code{diallel} random effect is fitted, an autocorrelated random effect \eqn{v_{ij}} is inferred for each \bold{unordered} pair (no individual effect is inferred), with correlation \eqn{\rho} between levels for pairs sharing one individual. This correlation parameter is fitted and is constrained by \eqn{\rho < 0.5} (see Details). \code{ranGCA} is equivalent to the case \eqn{\rho= 0.5}. \code{diallel} fits can be slow for large data if the correlation matrix is large, as this matrix can have a fair proportion of nonzero elements. 
There may also be identifiability issues for variance parameters: in a LMM as shown in the examples, there will be three parameters for the random variation (\code{phi}, \code{lambda} and \code{rho}) but only two can be estimated if only one observation is made for each dyad.  


%\eqn{\rho=0.5} actually generates a singular correlation matrix and the fit may then fail, or succeed with warnings, depending on the exact algorithm used internally by \pkg{spaMM}, 
% ... depending on .calc_denseness(cov_info_mat,relative=TRUE) < 0.15  => CHMfactor computed or not.  
%... Indeed, with \eqn{n} individuals, the distinct values of the random effect for up to \eqn{n(n-1)/2} pairs are all linearly related to the \eqn{n} individual values when \eqn{\rho=0.5}, so the rank of the correlation matrix is at most \eqn{n} even if it as more levels (more pairs).  
}
\details{
Although the symmetric random-effect structures may be used in many different contexts (including social network analysis, or \dQuote{round robin} experiments in psychology; another posibly relevant literature keyword here \dQuote{multi membership}), their present names refer to the semantics established for diallel experiments (e.g., Lynch & Walsh, 1998, p. 611), because it is not easy to find a more general yet intuitive semantics. If the names \code{ranGCA} and \code{diallel} sound inappropriate for your context of application, you can declare and use an alternative name for them, taking advantage of the fact that they are random-effect structures defined through \code{\link{corrFamily}} constructors, which are functions named as the formula term. For example, \code{symAdd(1|ID1+ID2)} can be used in a model formula after the following two steps:\cr
\preformatted{
# Define the 'symAdd' corrFamily constructor (a function) by copy:
symAdd <- ranGCA     
# Associate the 'symAdd' function to 'symAdd' formula terms:
register_cF("symAdd")  
}

In diallel experiments, one analyzes the phenotypes of offspring from multiple crosses among which the mother in a cross can be the father in another, so this is an example of multiple-membership. The additive genetic effects of each parent's genotypes are described as \dQuote{general combining abilities} (GCA). In case of non-additive effect, the half-sib covariance is not half the full-sib covariance and this is represented by the interaction \eqn{a_{ij}} described as \dQuote{specific combining abilities} (SCA). The sum of GCA and SCA defines a synthetic random effect \dQuote{received} by the offspring, with distinct levels for each unordered parental pair, and with correlation \eqn{\rho} between effects received by half-sibs (one shared parent). \eqn{\rho} corresponds to var(GCA)/[2*var(GCA)+var(SCA)] and is necessarily \eqn{\le 0.5}. 

See the \code{\link{X.GCA}} documentation for similar constructs for fixed effects.
}
\usage{
## formula terms:

# ranGCA(1| <.> + <.>) 
# antisym(1| <.> + <.>) 
# diallel(1|  <.> + <.>, tpar, fixed = NULL, public = NULL)

## where the <.> are two factor identifiers, ** whose levels
## must be identical when representing the same individual **

## corrFamily constructors:
ranGCA() # no argument
antisym() # no argument
diallel(tpar, fixed = NULL, public = NULL)
}
\arguments{
  \item{tpar}{
  Numeric: template value of the correlation coefficient for pairs sharing one individual. 
}
  \item{fixed}{
NULL or fixed value of the correlation coefficient.
}
\item{public}{NULL, or an environment. When an empty environment is provided, a template \code{CorNA} for the correlation matrix (with NA's in place of \eqn{\rho}) will be copied therein, for inspection at user level.}
}
\value{
The functions return corrFamily descriptors whose general format is described in \code{\link{corrFamily}}. The ones produced by \code{ranGCA} and \code{antisym} are atypical in that only their \code{Af} element is non-trivial.
}
\references{
Lynch, M., Walsh, B. (1998) Genetics and analysis of quantitative traits. Sinauer, Sunderland, Mass.
}
\examples{

#### Simulate dyadic data

set.seed(123)
nind <- 10       # Beware data grow as O(nind^2)
x <- runif(nind^2) 
id12 <- expand.grid(id1=seq(nind),id2=seq(nind))
id1 <- id12$id1
id2 <- id12$id2
u <-  rnorm(nind,mean = 0, sd=0.5)

## additive individual effects:
y <-  0.1 + 1*x + u[id1] +  u[id2] + rnorm(nind^2,sd=0.2)

## Same with non-additive individual effects:
dist.u <- abs(u[id1] -  u[id2])
z <- 0.1 + 1*x + dist.u + rnorm(nind^2,sd=0.2)

## anti-smmetric individual effects:
t <-  0.1 + 1*x + u[id1] - u[id2] + rnorm(nind^2,sd=0.2)

dyaddf <- data.frame(x=x, y=y, z=z, t=t, id1=id1,id2=id2)
# : note that this contains two rows per dyad, which avoids identifiability issues.

# Enforce that interactions are between distinct individuals (not essential for the fit):
dyaddf <- dyaddf[- seq.int(1L,nind^2,nind+1L),] 


# Fits:

(addfit <- fitme(y ~x +ranGCA(1|id1+id2), data=dyaddf))
#
# practically equivalent to:
#
(fitme(y ~x +diallel(1|id1+id2, fixed=0.49999), data=dyaddf))

(antifit <- fitme(t ~x +antisym(1|id1+id2), data=dyaddf))

(distfit <- fitme(z ~x +diallel(1|id1+id2), data=dyaddf)) 
     
}
