\name{multIMRF}
\alias{IMRF}
\alias{multIMRF}
\alias{str.inla.spde2}
\alias{small_spde}
\title{
Interpolated Markov Random Field models
}
\description{

\code{IMRF} is a syntax to specify random-effect terms of the forms considered by Lindgren et al. (2011) or Nychka et al. (2015, 2019).
For example, using \code{IMRF} with its \code{model} argument provides good approximations of random effects with \code{Matern} correlation structure win smoothness=1. 

The random effects considered here all involve a multivariate Gaussian random effect over a lattice, from which the random-effect value in any spatial position is determined by interpolation of values on the lattice. \bold{IMRF} stands for \bold{I}nterpolated \bold{M}arkov \bold{R}andom \bold{F}ield because the specific process considered on the lattice is currently known as a Gaussian Markov Random Field (see the Details for further information). Lindgren et al. considered irregular lattices that can be specified by the \code{model} argument, while Nychka et al. considered regular grids that can be specified by the other arguments. 

The \code{multIMRF} syntax implements the multiresolution model of Nychka et al. Any \code{multIMRF} term in a formula is immediately converted to \code{IMRF} terms, which should be counted as distinct random effects for all purposes (e.g., for fixing the variances of other random effects).  However, the arguments that control \code{multIMRF} terms are lists with names referring to successive \code{multIMRF} terms in the un-expanded formula, not to successive random-effect terms in the expanded formula.  
}
\details{
Gaussian Markov Random Field (MRF) and conditional autoregressive models are essentially the same thing, apart from details of specification. % Rue Held p. 28
\code{\link{adjacency}} and \code{\link{AR1}} random effects can be seen as specific MRFs.  
The common idea is the Markov-like property that the distribution of each element \eqn{b_i} of the random-effect \bold{b}, given values of a few specific elements (the \dQuote{neighbours} of \eqn{i}), is independent of other elements (i.e., of non-neighbours). The non-zero non-diagonal elements of a precision matrix characterize the neighbours. 

Given the inferred vector \bold{b} of values of the MRF on the lattice, the interpolation of the MRF in any focal point is of the form \bold{Ab} where each row of \bold{A} weights the values of \bold{b} according to the position of the focal point relative to the vertices of the lattice. Following the original publications, for regular grids (NULL \code{model}), the weights are computed as <Wendland function>(<scaled Euclidean distances beween focal point and vertices>); and for grids given by \code{model=<inla.spde2 object>}, the non-zero weights are the barycentric coordinates of the focal point in the enclosing triangle from the mesh triangulation (points from outside the mesh would have zero weights, so the predicted effect \bold{Ab=0}). 

The \code{IMRF} model defines both a lattice in space, the precision matrix for a Gaussian MRF over this lattice, and the \bold{A} weights. The full specification of the MRF on \bold{irregular lattices} is complex. The \eqn{\kappa} parameter considered by \code{spaMM} is the \eqn{\kappa} considered by Lindgren et al.; only the case \emph{d=2} and \eqn{alpha=} 1 or 2, approximating a Matérn correlation model with \eqn{nu=0} or \eqn{nu=1} (\eqn{alpha=nu + d/2}), is currently implemented in \code{spaMM}.  

For the MRFs on \bold{regular grids} implemented here, the precision matrix is defined (up to a variance parameter) as \bold{M'M} where the diagonal elements \eqn{m_{ii}} of \bold{M} are 4+\eqn{\kappa^2} and the \eqn{m_{ij}} for the four nearest neighbours are -1 (note that \bold{M'M} involves more than these four neighbours). 

The precision matrix defined in this way is the inverse of an heteroscedastic covariance matrix \bold{C}, but by default a normalization is applied so that the random effect is homoscedastic. As for other random effects, the variance is further controlled by a multiplicative factor \eqn{\lambda}.  The \bold{ormalization} is as follows: the design matrix of the random effect term is viewed as \bold{WAL} where \bold{W} is a diagonal normalization matrix, \bold{A} is the above-described weight matrix, and \bold{L} is a \dQuote{square root} of \bold{C}. If no normalization is applied, the covariance matrix of the random effect is of the form \eqn{\lambda}\bold{ALL'A'}, which is heteroscedastic; \eqn{\lambda} may then be quite different from the marginal variance of the random effect, and is difficult to describe in a simple way. Hence, by default,  \bold{W} is defined such that \bold{WALL'A'W'} has unit diagonal; then, \eqn{\lambda} is the marginal variance of the random effect. 

By default, the IMRF lattice is rectangular (currently the only option) and is made of a core lattice, to which margins of \code{margin} steps are added on each side. The core lattice is defined as follows: in each of the two spatial dimensions, the range of axial coordinates is determined. The largest range is divided in \code{nd-1} steps, determining \code{nd} values and step length \eqn{L}. The other range is divided in steps of the same length \eqn{L}. If it extends over (say) \eqn{2.5 L}, a grid of 2 steps and 3 values is defined, and by default centered on the range (the extreme points therefore typically extend slightly beyond the grid, within the first of the additional steps defined by the \code{margin}; if not centered, the grid start from the lower coordinate of the range).

\code{multIMRF} implements multilevel IMRFs. It defines a sequence of IMRFs, with progressively finer lattices, a common \eqn{\kappa} value \code{hy_kap} for these IMRFs, and a single variance parameter \code{hy_lam} that determines \eqn{\lambda} values decreasing by a factor of 4 for successive IMRF terms. By default, each component \code{IMRF} is normalized independently as described above (as in Nychka et al. 2019), and \code{hy_lam} is the sum of the variances of these terms (e.g., if there are three levels and \code{hy_lam=1}, the successive variances are (1,1/4,1/16)/(21/16) ). The \code{nd} of the first IMRF is set to the \code{coarse} value, and its lattice is defined accordingly. If \code{coarse=4}  and \code{margin=5}, a grid of 14 coordinates is therefore defined over the largest range. In the second IMRF, the grid spacing is halved, so that new steps are defined halfway between the previous ones (yielding a grid of 27 stepin the widest range). The third IMRF proceeds from the second in the same way, and so on.

To control initial or fixed values of \code{multIMRF} \eqn{\kappa} and variance parameters, the \code{hyper} syntax shown in the Examples should be used. \code{hyper} possible elements are named \code{"1"}, \code{"2"}, ... referring to successive \code{multIMRF} terms in the formula.
}
\usage{
# IMRF( 1 | <coordinates>, model, nd, m, no, ce, ...) 
# multIMRF( 1 | <coordinates>, levels, margin, coarse=10L, 
#            norm=TRUE, centered=TRUE ) 
}
% IMRF( 1 | <coordinates>, model, nd, m, no, ce, ...) defaults implemented in .parseBars()
% multIMRF( 1 | . ) defaults in .expand_multIMRF() default args
\arguments{
  \item{model}{An \code{inla.spde2} object as produced by \code{INLA::inla.spde2.matern} (see Examples below, and \url{http://www.r-inla.org} for further information).
}
  \item{levels}{integer; 
Number of levels in the hierarchy, i.e. number of component IMRFs.
}
  \item{margin, m}{integer; 
  width of the margin, as a number of additional grid points on each side (applies to all levels of the hierarchy).
}
  \item{coarse}{integer; 
  number of grid points (excluding the margins) per dimension for the coarsest IMRF. The number of grids steps nearly doubles with each level of the hierarchy (see Details). 
}
  \item{nd}{integer; 
  number of grid steps (excluding the margins) per dimension for the given IMRF. 
}
  \item{norm, no}{Boolean; whether to apply normalization (see Details), or not.
}
  \item{centered, ce}{Boolean; whether to center the grid in all dimensions, or not.
}
\item{\ldots}{Not documented, for programming purposes} % l used in .assign_geoinfo_and_LMatrices_but_ranCoefs()
}
\references{
D. Nychka, S. Bandyopadhyay, D. Hammerling, F. Lindgren, S. Sain (2015)
A multiresolution gaussian process model for the analysis of large spatial datasets.
Journal of Computational and Graphical Statistics 24 (2), 579-599. \doi{10.1080/10618600.2014.914946}

D. Nychka, D. Hammerling, Mitchel. Krock, A. Wiens (2018) Modeling and emulation of nonstationary Gaussian fields. Spat. Stat. 28: 21-38. \doi{10.1016/j.spasta.2018.08.006}

%Krainski et al. (2018) Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. \url{https://becarioprecario.bitbucket.io/spde-gitbook}

Lindgren F., Rue H., Lindström J. (2011) An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73: 423-498. \doi{10.1111/j.1467-9868.2011.00777.x}
}
\examples{
data("blackcap") ## toy examples; but IMRF may be useful only for much larger datasets

########################## Irregular lattice specified by 'model':

if (spaMM.getOption("example_maxtime")>1.6) {

data("small_spde") ## load object of class 'inla.spde2', created and saved by :
  # spd <- sp::SpatialPointsDataFrame(coords = blackcap[, c("longitude", "latitude")],
  #                            data = blackcap)
  # small_mesh <- INLA::inla.mesh.2d(loc = INLA::inla.mesh.map(sp::coordinates(spd)), 
  #                           max.n=100, # only for demonstration purposes
  #                           max.edge = c(3, 20)) 
  # small_spde <- INLA::inla.spde2.matern(small_mesh)
  # save(small_spde, file="small_spde.RData", version=2)
  
fit_SPDE <- fitme(migStatus ~ means + IMRF(1|longitude+latitude, model=small_spde), 
                  data=blackcap)
}

########################## Regular lattices:   
  # Using 'hyper' to control fixed hyper-parameters
  (mrf <- fitme(migStatus ~ 1 + (1|pos) + 
                            multIMRF(1|latitude+longitude,margin=5,levels=2), 
                data=blackcap, method="ML", fixed =list(phi=1,lambda=c("1"=0.5),
                hyper=list("1"=list(hy_kap=0.1,hy_lam=1)))) )
  if (spaMM.getOption("example_maxtime")>5) {
    # Using 'hyper' to control initial hyper-parameters
    (mrf <- fitme(migStatus ~ 1 + multIMRF(1|latitude+longitude,margin=5,levels=2),
                  data=blackcap, method="ML",fixed =list(phi=1),
                  init=list(hyper=list("1"=list(hy_kap=0.1,hy_lam=1)))) )
    # *Independent* IMRF terms (often giving dubious results)
    (mrf <- HLCor(migStatus ~ 1 + IMRF(1|latitude+longitude,margin=5, nd=4L)
                                + IMRF(1|latitude+longitude,margin=5, nd=7L),
            data=blackcap, HLmethod="ML",
            ranPars=list(phi=1,lambda=c(1/4,1/16),
                         corrPars=list("1"=list(kappa=0.1),"2"=list(kappa=0.1)))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }% use one of  RShowDoc("KEYWORDS")
\keyword{ spatial }% __ONLY ONE__ keyword per line
