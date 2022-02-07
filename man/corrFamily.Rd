\name{corrFamily}
\alias{corrFamily}
\title{
Defining parametric correlation families
}
\description{
The \code{corrFamily} keyword can be used to declare and fit correlated random effects belonging to a user-defined correlation model (i.e., a parametric family of correlation matrices), or a covariance model. This feature is experimental in the sense that it has been tested in a limited number of applications and a limited number of matrix formats (see recommendations below). In the long term it should facilitate the addition of new \dQuote{canned} correlation models in \pkg{spaMM}.

In the model formula, the random effect is declared as a term of the form\cr
\code{corrFamily(1|<grouping factor>)}; and the \code{covStruct} argument is used to provide 
\itemize{
\item a function \code{f} that returns a correlation (or covariance) matrix from an input parameter vector (see below further information on class and attributes of the returned matrix). 

\item a template parameter vector \code{tpar}, also with some constraints described below. 
}

\code{lower} and \code{upper} value must be specified for the correlation parameters to be estimated.
%\preformatted{
%list(
%  corrFamily=list( # defines possibly non-zero elements of the correlation matrix:
%    map=Matrix::drop0(toeplitz(seq(0,6))),    # defines its variable elements
%    fixed=Matrix::.symDiagonal(7)    # defines its constant non-zero elements
%))   
%}
%In this example, \code{fixed} is an identity matrix: it serves to assign a constant value 1 to the main diagonal of the correlation matrix. The \code{map} has constant integer values along the other diagonals, 1 along the first super/sub diagonals, 2 along the second diagonals, etc. This defines the first parameter as the element of the first diagonals, the second parameter as the elements of the second diagonals, etc. Note that \code{toeplitz} is filled with \code{0} on the main diagonal, then these \code{0} are \code{Matrix::drop0}()ed so that the resulting sparse matrix has elements only in variable positions of the correlation matrix to be constructed.   

Using \code{corrFamily} terms will be more efficient than the equally general method, of maximizing an objective function  of the correlation parameters that calls (say) \code{fitme()} on a model including a \code{corrMatrix} itself function of the correlation parameters. But this may still be inefficient if a few issues are ignored. The Details should be read before conceiving a \code{corrFamily} model. 
}
\usage{
## 'covStruct' element
# corrFamily=list(f=<.>, tpar=<.>)

## 'lower', 'upper', 'init' and 'fixed' values in fitting-function call
## are numeric vectors in the 'corrPars', provided 
## following the general syntax of these arguments; e.g.,
# lower=list(corrPars=list("1"= <numeric vector> )) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{f}{(required): function returning the correlation (or covariance) matrix, given its first argument, a parameter vector. The function value should
(1) be of constant class for all parameter values. For families of mathematically sparse matrices, the \code{CsparseMatrix} class is recommended (and more specifically the \code{dsCMatrix} class since the matrix is symmetric); (2) have row names that match the levels of the grouping factor  (the nested random effect Example shows the code needed when this nested effect is defined from two variables).}
\item{tpar}{(required): a feasible argument of \code{f}. If it is a named vector, its names will be used for the parameters, otherwise default names \code{"p1"}, \code{"p2"}... are provided. In order to favor the automatic selection of suitable algorithms, \code{tpar} should be chosen so that \code{f(tpar)} is \bold{least} sparse in the correlation family, in terms of its sparsity and of the sparsity of its inverse. A \code{tpar} yielding an identity matrix is often a \bold{*bad*} template as least sparse correlation matrices and their inverses are denser for most families except diagonal ones.}
%  \item{map}{(required): symmetric sparse matrix (from a \pkg{Matrix} class) without any explicit zero.}
%  \item{fixed}{While this argument is formally optional, 
%  it is required for correlation matrices whose diagonal must be fixed to 1. This case can be specified as \code{fixed="unit diag"}. More generally,
%  a sparse matrix (from a \pkg{Matrix} class) without any explicit zero. }
%Both \code{map} and \code{fixed} will be to sparse symmetric matrices (\code{dsCMatrix} class) internally by \code{Matrix::forceSymmetric}, so asymmetries in them will be ignored.
\item{lower,upper}{(required): complete vector of lower and upper bounds. Parameter names (as determined by the \code{tpar} argument) are not required but will be used if present.} 
\item{init}{(optional): a vector of initial values. If it is incomplete (not providing values for all parameters), names of the specified initial parameter values must be provided.} 
%(given the distinctive use of \code{tpar}, \code{spaMM} never tries to use it as initial value).
\item{fixed}{(optional): vector of fixed values. Currently only a complete vector is formally handled and there is no check that the vector is complete.}
}
\details{
\itemize{
\item \pkg{spaMM} will regularize invalid or nearly-singular correlation or covariance matrices internally if the correlation function has not done so already, but it it better to control this in the correlation function. The \code{\link{regularize}} convenience function is available for that purpose, but parametrizations that avoid the need for regularization are even better, since fitting models with nearly-singular correlation matrices is prone to various difficulties (The Toeplitz example below is poor from that perspective, although it is good to illustrate potential problems).

%If a candidate correlation matrix is found not to have positive eigenvalues, a warning is issued and constant correction is added to the diagonal so that the corrected matrix is positive definite (though nearly singular). 

\item Users should make sure that any regularized matrix still belongs to the intended parametric family of matrices, and they should keep in mind that the spaMM output will show the input parameters of the unregularized matrix, not the parameters of the regularized one (e.g., in the Toeplitz example below, the fitted matrix is a regularized Toepliz matrix with slightly different coefficients than the input parameters).  

And for efficiency,

 \item Let us repeat that the correlation function should return matrices of constant class, and in sparse format when the matrices are indeed mathematically sparse. For mathematically dense matrices (as in the Toeplitz example below), the \code{dsyMatrix} class may be suitable.
 
 \item Let us repeat that in order to favor the automatic selection of suitable algorithms, \code{tpar} should be chosen so that \code{f(tpar)} is \bold{least} sparse in the correlation family. For matrices of \code{CsparseMatrix}, a check is implemented to catch wrong choices of \code{tpar}.
 
 \item For challenging problems (large data, many parameters...) it may pay to optimize a bit the correlation function. The Example of nested effects with heterogenous variance below illustrates a possible trick. In the same cases,  It may also pay to try the alternative \code{\link{algebra}}ic methods, by first comparing speed of the different methods (\code{control.HLfit=list(algebra= <"spprec"|"spcorr"|"decorr">)}) for given correlation parameter values, rather than to assume that \pkg{spaMM} will find the best method (even if it often does so).
}
}
\examples{
if (spaMM.getOption("example_maxtime")>2 &&
      requireNamespace("agridat", quietly = TRUE)) {

data("onofri.winterwheat", package="agridat")

##### Fitting a Toeplitz correlation model for temporal correlations #####

# A Toeplitz correlation matrix of dimension d*d has d-1 parameters 
# (by symmetry, and with 1s on the main diagonal). These d-1 parameters 
# can be fitted as follows:

Toepfn <- function(v) {
  toepmat <- Matrix::forceSymmetric(toeplitz(c(1,v))) # dsyMatrix
  # Many of the matrices in this family are not valid correlation matrices;
  #   the regularize() function is handy here:
  toepmat <- regularize(toepmat, maxcondnum=1e12)
  # And don't forget the rownames!
  rownames(toepmat) <- unique(onofri.winterwheat$year)
  toepmat
}

(Toepfit <- spaMM::fitme(
  yield ~ gen + corrFamily(1|year), data=onofri.winterwheat, method="REML",
  covStruct=list(corrFamily=list(f=Toepfn, tpar=rep(1e-4,6))), 
         # (Note the gentle warning if one instead uses tpar=rep(0,6) here)  
  lower=list(corrPars=list("1"=rep(-0.999,6))), 
  upper=list(corrPars=list("1"=rep(0.999,6))))) 

# The fitted matrix is (nearly) singular, and was regularized:

eigen(Corr(Toepfit)[[1]])$values

# which means that the returned likelihood may be inaccurate, 
# and also that the actual matrix elements differ from input parameters:

Corr(Toepfit)[[1]][1,-1]  

### The usual rules for specifying covStruct, 'lower', 'upper' and 'init' apply
# here when the corrFamily term is the second random-effect:

(Toep2 <- spaMM::fitme(
        yield ~ 1 + (1|gen) + corrFamily(1|year), data=onofri.winterwheat, method=method,
        covStruct=list("1"=NULL, corrFamily=list(f=Toepfn, tpar=rep(1e-4,6))), 
        , init=list(corrPars=list("2"=rep(0.1,6))),
        lower=list(corrPars=list("2"=rep(-0.999,6))), 
        upper=list(corrPars=list("2"=rep(0.999,6)))))

##### Fitting one variance among years per each of 8 genotypes. #####

# First, note that this can be *more efficiently* fitted by another syntax:

### Fit as a constrained random-coefficient model: 
    
# Diagonal matrix of NA's, represented as vector for its lower triangle:
ranCoefs_for_diag <- function(nlevels) { 
  vec <- rep(0,nlevels*(nlevels+1L)/2L)
  vec[cumsum(c(1L,rev(seq(nlevels-1L)+1L)))] <- NA
  vec
} 
    
(by_rC <- spaMM::fitme(yield ~ 1 + (0+gen|year), data=onofri.winterwheat, method="REML",
                       fixed=list(ranCoefs=list("1"=ranCoefs_for_diag(8)))))
                         
### Fit as a corrFamily model:   

gy_levels <- paste0(gl(8,1,length =56,labels=levels(onofri.winterwheat$gen)),":",
                        gl(7,8,labels=unique(onofri.winterwheat$year)))
                        
# A log scale is often suitable for variances, hence is used below;

# a correct but crude implementation of the model is
diagf <- function(logvar) {
  corr_map <- kronecker(Matrix::.symDiagonal(n=7),diag(x=exp(logvar)))
  rownames(corr_map) <- gy_levels
  corr_map
}

# but we can minimize matrix operations as follows:

corr_map <- Matrix::forceSymmetric(kronecker(Matrix::.symDiagonal(n=7),diag(x=seq(8))))
rownames(corr_map) <- gy_levels

diagf <- function(logvar) {
  corr_map@x <- exp(logvar)[corr_map@x]
  corr_map
}                 # (and this returns a dsCMatrix)

(by_cF <- spaMM::fitme(
  yield ~ 1 + corrFamily(1|gen \%in\% year), data=onofri.winterwheat, method="REML",
  covStruct=list(corrFamily = list(f=diagf, tpar=rep(1,8))), 
  fixed=list(lambda=1),            # Don't forget to fix this redundant parameter!
  # init=list(corrPars=list("1"=rep(log(O.1),8))), # 'init' optional 
  lower=list(corrPars=list("1"=rep(log(1e-6),8))), # 'lower' and 'upper' required
  upper=list(corrPars=list("1"=rep(log(1e6),8)))))

exp(get_ranPars(by_cF)$corrPars[[1]]) # fitted variances 
  
}
}
