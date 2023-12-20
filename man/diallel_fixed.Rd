\name{X.GCA}
\alias{X.GCA}
\alias{X.antisym}
\title{
Fixed-effect terms for dyadic interactions
}
\description{
\code{X.GCA} and \code{X.antisym} are functions which, when called in a model formula, stand for terms designed to represent the effet of symmetric interactions between pairs of individuals (order of individuals in the pair does not matter). \code{antisym} likewise represents anti-symmetric interactions (the effect of reciprocal ordered pairs on the outcome are opposite, as in the so-called Bradley-Terry models). These constructs all account for multiple membership, i.e., the fact that the same individual may act as the first or the second individual among different pairs, or even within one pair if this makes sense).   

The outcome of an interaction between a pair \eqn{i,j} of agents is subject to a symmetric overall effect \eqn{a_{ij}} when the effect  \dQuote{on} individual \eqn{i} (or viewed from the perspective of individual \eqn{i}) equals the effect on \eqn{j}: \eqn{a_{ij}=a_{ji}}. This may result from the additive effect of individual effects \eqn{a_{i}} and \eqn{a_{j}}:  \eqn{a_{ij}=a_i+a_j}. A \code{X.GCA} call represents such symmetric additive effects. Conversely, antisymmetry is characterized by \eqn{a_{ij}=a_i-a_j=-a_{ji}} and is represented by a \code{X.antisym} call. See the \code{\link{diallel}} documentation for similar constructs for random effects, for additional comments on semantics (e.g. about \dQuote{GCA}), and for further references. 

If individual-level factors ID1 + ID2 were included in a formula for dyadic interactions, this would result in different coefficients being fitted for the same level in each factor. By contrast, the present constucts ensure that a single coefficient is fitted for the same-named levels of factors ID1 and ID2. 
}
\usage{
X.GCA(term, contr="contr.treatment", ...) 
X.antisym(term, contr="contr.treatment", ...) 
}
\arguments{
  \item{term}{
an expression of the form <.>:<.> where each <.> represents a factor (or a variable that will automaticall be converted to a factor) present in the \code{data} of the fitting function call. 
}
  \item{contr}{
The \code{\link{contrasts}} used. Only the default and \code{"contr.sum"} are implemented.
}
  \item{\dots}{
  For programming purposes, not documented.
}
}
\value{
The functions return design matrices with additional class \code{"factorS"} and attributes \code{"call"} and \code{"spec_levs"}.
}
\details{
The fixed-effect terms (\code{GCA(Par1,Par2)}, etc) from the \pkg{lmDiallel} package (Onofri & Terzaroli, 2021), defined by functions returning design matrices for usage with \code{stats:lm}, work in a formula for a \pkg{spaMM} fit, and users can define use such functions as templates for additional functions that will work similarly. However, not all post-fit functions will handle terms defined in this way well: checking robustness of \code{predict} on small and permuted \code{newdata}, as shown in the Examples, is a good test. Such problems happen because the formula-handling machinery of R handles terms represented by either a matrix or a factor, while both the model matrix and the factor information used to construct dyadic-interaction terms are needed to correctly predict, with new data, from models including such terms.  

The presently designed functions are defined to solve this issue. By using such functions as template, users can define additional functions with the same return format (as further explained in the documented source code of \code{X.antisym}), which will allow them to perform correct predictions from fitted models.    
}
\references{
Onofri A., Terzaroli N. (2021) lmDiallel:  Linear fixed effects models for diallel crosses. Version 0.9.4. \url{https://cran.r-project.org/package=lmDiallel}.
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

## anti-smmetric individual effects:
t <-  0.1 + 1*x + u[id1] - u[id2] + rnorm(nind^2,sd=0.2)

dyaddf <- data.frame(x=x, y=y, t=t, id1=id1,id2=id2)
# : note that this contains two rows per dyad, which avoids identifiability issues.

# Enforce that interactions are between distinct individuals (not essential for the fit):
dyaddf <- dyaddf[- seq.int(1L,nind^2,nind+1L),] 


# Fits:

(addfit <- fitme(y ~x +X.GCA(id1:id2), data=dyaddf))

(antifit <- fitme(t ~x +X.antisym(id1:id2), data=dyaddf))

if (FALSE) { #### check of correct handling by predict():

  # First scramble the data so that input factors are in no particular order
  set.seed(123)
  dyaddf <- dyaddf[sample(nrow(dyaddf)),]
  
  addfiti <- fitme(y ~x +X.GCA(id1:id2), data=dyaddf)
  foo <- rev(2:4)
  p1 <- predict(addfiti)[foo]
  p2 <- predict(addfiti, newdata=dyaddf[foo,])
  diff(range(p1-p2))<1e-10    # must be TRUE
}
}
