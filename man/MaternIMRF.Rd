\name{MaternIMRFa}
\alias{MaternIMRFa}
\title{
corrFamily constructor for Interpolated Markov Random Field (IMRF) covariance structure approximating a 2D Matern correlation model.
}
\description{
Reimplements the \code{\link{IMRF}} correlation model approximating a Matern correlation function, through a \code{\link{corrFamily}} constructor. This allows the efficient joint estimation of the \code{alpha} parameter of the approximating Markov random field (in principle related to the smoothness parameter of the \code{\link{Matern}} correlation function) together with its \code{kappa} parameter. By contrast, random effects terms specified as \code{IMRF(1| . , model = <INLA::inla.spde2.matern result>)}  assume a fixed \code{alpha}.

Using this feature requires that the not-on-CRAN package \pkg{INLA} (\url{https://www.r-inla.org}) is installed so that \code{INLA::inla.spde2.matern} can be called for each \code{alpha} value. 
}
\usage{
# corrFamily constructor:
MaternIMRFa(mesh, tpar = c(alpha = 1.25, kappa = 0.1), fixed = NULL)
}
\arguments{
  \item{mesh}{
An \code{inla.mesh} object as produced by \code{INLA::inla.mesh.2d}.
}
and consistently with the general format of \code{\link{corrFamily}} constructors:
  \item{tpar}{
  Named numeric vector: template values of the parameters of the model. Better not modified unless you know what you are doing. 
}
  \item{fixed}{NULL or numeric vector, to fix the parameters of this model.}
}
\value{
A list suitable as input in the \code{covStruct} argument, with the following elements: 
\item{f }{function returning a precision matrix for the random effect in mesh vertices;} 
\item{tpar }{template parameter vector (see general requirements of a \code{\link{corrFamily}} descriptor);} 
\item{Af}{function returning a matrix that implements the prediction of random effect values in data locations by interpolation of values in \code{mesh} locations (similarly to \code{INLA::inla.spde.make.A});} 
\item{type}{specifies that the matrix returned by \code{Cf} is a precision matrix rather than a correlation matrix;}  
and possibly other elements which should not be considered as stable features of the return value.

}
\references{
Lindgren F., Rue H., Lindström J. (2011) An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73: 423-498. \doi{10.1111/j.1467-9868.2011.00777.x}
}
\examples{
\dontrun{

if(requireNamespace("INLA", quietly = TRUE)) {

data("Loaloa")

mesh <- INLA::inla.mesh.2d(loc = Loaloa[, c("longitude", "latitude")],
                     max.edge = c(3, 20)) 

### Fit with fixed alpha                     

(fit_MaternIMRF <- fitme(
  cbind(npos,ntot-npos) ~ elev1 + elev2 + elev3 + elev4 + maxNDVI1 + 
        seNDVI +  MaternIMRFa(1|longitude+latitude, mesh, fixed=c(alpha=1.05)),
  family=binomial(), 
  data=Loaloa, verbose=c(TRACE=interactive())) )
   
# For data sets with a small number of locations (as here), fitting
# the Matern model as follows is faster than fitting its MaternIMRFa approximation. 
# Here this appears more than twofold faster

fit_Matern <- fitme(
  cbind(npos,ntot-npos) ~ elev1 + elev2 + elev3 + elev4 + maxNDVI1 + 
        seNDVI + Matern(1|longitude+latitude), 
  fixed=list(nu=0.05),                          # in principle similar to alpha=0.05
  data=Loaloa,family=binomial())   
  
  
  
### Same with variable alpha                     

(fit_MaternIMRF <- fitme(
  cbind(npos,ntot-npos) ~ elev1 + elev2 + elev3 + elev4 + maxNDVI1 + 
        seNDVI + MaternIMRFa(1|longitude+latitude, mesh),         
  family=binomial(), 
  data=Loaloa, verbose=c(TRACE=interactive())) )

# Comparable Matern fit:
fit_Matern <- fitme(
  cbind(npos,ntot-npos) ~ elev1 + elev2 + elev3 + elev4 + maxNDVI1 + 
        seNDVI + Matern(1|longitude+latitude), 
  init=list(nu=0.25), lower=list(nu=0), upper=list(nu=1),
  data=Loaloa,family=binomial())   
  
# Note that the fitted nu and alpha parameters do not quite match each other,
# and that the IMRF likelihood does not really approximate the Matern likelihood.
# Mesh design would also be seen to matter.   

} else print("INLA must be installed to run this example.")

}
}