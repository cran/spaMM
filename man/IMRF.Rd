\name{multIMRF}
\alias{IMRF}
\alias{multIMRF}
\alias{hyper}
\alias{str.inla.spde2}
\alias{inla.spde2.matern}
\alias{inla.spde2.pcmatern}
\alias{small_spde}
\title{
Interpolated Markov Random Field models
}
\description{

\pkg{spaMM} can fit random-effect terms of the forms considered by Lindgren et al. (2011) or Nychka et al. (2015, 2018).
The random effects considered here all involve a multivariate Gaussian random effect over a lattice, from which the random-effect value in any spatial position is determined by interpolation of values on the lattice. \bold{IMRF} stands for \bold{I}nterpolated \bold{M}arkov \bold{R}andom \bold{F}ield because the specific process considered on the lattice is currently known as a Gaussian Markov Random Field (see the Details for further information). Lindgren et al. considered irregular lattices designed to approximate of the Matern correlation model with fixed smoothness <= 2, while Nychka et al. considered regular grids. 

The correlation model of Lindgren et al. (2011) can be fitted by \pkg{spaMM} by declaring an \code{IMRF} random effect term in the model formula, with a \code{model} argument in the right-hand side whose value is the result of \code{INLA::inla.spde2.matern} (or \code{INLA::inla.spde2.pcmatern}) for given smoothness. The \pkg{spaMM} functions for such a fit do not call \pkg{INLA} functions. Alternatively, the same model with variable smoothness can be fitted by declaring a \code{corrFamily} term whose structure is described through the \code{\link{MaternIMRFa}} function, whose respective documentations should be considered for more details. In the latter case \code{INLA::inla.spde2.matern} is called internally by \code{spaMM}. The correlation models thus defined are fitted by the same methods as other models in \pkg{spaMM}.     

Regular lattices can also be declared as an \code{IMRF} term (with arguments distinct from \code{model}). The \code{multIMRF} syntax implements the multiresolution model of Nychka et al. Any \code{multIMRF} term in a formula is immediately converted to \code{IMRF} terms for regular grids wih different step sizes. This has distinct implications for controlling the parameters of these or other random effects in the model by \code{init} or \code{fixed} values: see Details if you need such control. 
}
\details{
\bold{Formulation of the covariance models:}

Gaussian Markov Random Field (MRF) and conditional autoregressive models are essentially the same thing, apart from details of specification. % Rue Held p. 28
\code{\link{adjacency}} and \code{\link{AR1}} random effects can be seen as specific MRFs.  
The common idea is the Markov-like property that the distribution of each element \eqn{b_i} of the random-effect \bold{b}, given values of a few specific elements (the \dQuote{neighbours} of \eqn{i}), is independent of other elements (i.e., of non-neighbours). The non-zero non-diagonal elements of a precision matrix characterize the neighbours. 

Given the inferred vector \bold{b} of values of the MRF on the lattice, the interpolation of the MRF in any focal point is of the form \bold{Ab} where each row of \bold{A} weights the values of \bold{b} according to the position of the focal point relative to the vertices of the lattice. Following the original publications,\cr 
{ * } for grids given by \code{model=<inla.spde2 object>}, the non-zero weights are the barycentric coordinates of the focal point in the enclosing triangle from the mesh triangulation (points from outside the mesh would have zero weights, so the predicted effect \bold{Ab=0});\cr { * }for regular grids (NULL \code{model}), the weights are computed as <Wendland function>(<scaled Euclidean distances between focal point and vertices>). 

The \code{IMRF} model defines both a lattice in space, the precision matrix for a Gaussian MRF over this lattice, and the \bold{A} matrix of weights. The full specification of the MRF on \bold{irregular lattices} is complex. The \eqn{\kappa} (\code{kappa}) parameter considered by \code{spaMM} is the \eqn{\kappa} scale parameter considered by Lindgren et al and comparable to the \eqn{\rho} scale factor of the Matérn model. The \eqn{\alpha} argument of the \code{INLA::inla.spde2.matern} controls the smoothness of the approximated Matern model, as \eqn{\alpha=\nu + d/2}) where \eqn{d} is the dimension of the space.   
Correlation models created by \code{INLA::inla.spde2.pcmatern} are handled so as to give the same correlation values as when \code{INLA::inla.spde2.matern} is used with the same \code{mesh} and \code{alpha} argument (thus, the extra functionalities of \dQuote{\code{pc}}\code{matern} are ignored). 

Not all options of the INLA functions may be compatible or meaningful when used with spaMM (only the effects of \code{alpha} and \code{cutoff} have been checked). 
% cf bit of code in .calc_IMRF_Qmat()

\bold{Normalization}:

For the MRFs on default \bold{regular grids} (missing \code{model} argument), the precision matrix is defined (up to a variance parameter) as \bold{M'M} where the diagonal elements \eqn{m_{ii}} of \bold{M} are 4+\eqn{\kappa^2} and the \eqn{m_{ij}} for the four nearest neighbours are -1 (note that \bold{M'M} involves more than these four neighbours). The precision matrix defined in this way is the inverse of an heteroscedastic covariance matrix \bold{C}, but (following Nychka et al.) by default a normalization is applied so that the random effect in each data position is homoscedastic (the precision matrix for the latent effect in grid positions is not modified, but it is the \bold{A} matrix of weights which is is modified). As for other random effects, the variance is further controlled by a multiplicative factor \eqn{\lambda}. 

Without normalization, the covariance matrix of the random effect in data locations is \eqn{\lambda}\bold{ALL'A'} (\bold{A} being the above-described weight matrix, and \bold{L} is a \dQuote{square root} of \bold{C}), and  \bold{AL} is the original \dQuote{design matrix} of the random effect. \eqn{\lambda} may then be quite different from the marginal variance of the random effect, and is difficult to describe in a simple way. 
For normalization, \bold{A} is modified as \bold{WA} where \bold{W} is a diagonal matrix such that \bold{WAL} is a correlation matrix (\bold{WALL'A'W'} has unit diagonal); then, \eqn{\lambda} is the marginal variance of the random effect. 


For irregular grids specified using the \code{model} argument, the precision matrix described by this object is also the inverse of an heteroscedastic covariance matrix, but here (again following original publicatiosn such as Lindgren at al. 2011) the normalization is not applied by default (and was not even an option before version 4.3.23). But for ease of presentation and interpretation, if for no other reason, the normalized model may be preferable. 

\bold{Details for rectangular grids:} 

By default (meaning in particular that \code{model} is not used to specify a lattice defined by the INLA procedures), the IMRF lattice is rectangular (currently the only option) and is made of a core lattice, to which margins of \code{margin} steps are added on each side. The core lattice is defined as follows: in each of the two spatial dimensions, the range of axial coordinates is determined. The largest range is divided in \code{nd-1} steps, determining \code{nd} values and step length \eqn{L}. The other range is divided in steps of the same length \eqn{L}. If it extends over (say) \eqn{2.5 L}, a grid of 2 steps and 3 values is defined, and by default centered on the range (the extreme points therefore typically extend slightly beyond the grid, within the first of the additional steps defined by the \code{margin}; if not centered, the grid start from the lower coordinate of the range).

\code{multIMRF} implements multilevel IMRFs. It defines a sequence of IMRFs, with progressively finer lattices, a common \eqn{\kappa} value \code{hy_kap} for these IMRFs, and a single variance parameter \code{hy_lam} that determines \eqn{\lambda} values decreasing by a factor of 4 for successive IMRF terms. By default, each component \code{IMRF} is normalized independently as described above (as in Nychka et al. 2019), and \code{hy_lam} is the sum of the variances of these terms (e.g., if there are three levels and \code{hy_lam=1}, the successive variances are (1,1/4,1/16)/(21/16) ). The \code{nd} of the first IMRF is set to the \code{coarse} value, and its lattice is defined accordingly. If \code{coarse=4}  and \code{margin=5}, a grid of 14 coordinates is therefore defined over the largest range. In the second IMRF, the grid spacing is halved, so that new steps are defined halfway between the previous ones (yielding a grid of 27 step in the widest range). The third IMRF proceeds from the second in the same way, and so on.

To control initial or fixed values of \code{multIMRF} \eqn{\kappa} and variance parameters, which are hyper-parameter controlling several \code{IMRF} terms, the \code{hyper} syntax shown in the Examples should be used. \code{hyper} is a nested list whose possible elements are named \code{"1"}, \code{"2"}, ... referring to successive \code{multIMRF} terms in the input formula, not to successive random effect in the expanded formula with distinct IMRF terms (see Examples). But the different IMRF terms should be counted as distinct random effects when controlling other parameters (e.g., for fixing the variances of other random effects).
}
\usage{
# IMRF( 1 | <coordinates>, model, nd, m, no, ce, ...) 
# multIMRF( 1 | <coordinates>, levels, margin, coarse=10L, 
#            norm=TRUE, centered=TRUE ) 
}
% IMRF( 1 | <coordinates>, model, nd, m, no, ce, ...) defaults implemented in .parseBars()
% multIMRF( 1 | . ) defaults in .expand_multIMRF() default args
\arguments{
  \item{model}{An \code{inla.spde2} object as produced by \code{INLA::inla.spde2.matern} or\cr 
  \code{INLA::inla.spde2.pcmatern} (see Examples below, and \url{https://www.r-inla.org} for further information).
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

if (spaMM.getOption("example_maxtime")>6) {

data("blackcap") ## toy examples; but IMRF may be useful only for much larger datasets
## and when using the 'cutoff' parameter of INLA::inla.mesh.2d()

########################## Irregular lattice specified by 'model':
#
data("small_spde") ## load object of class 'inla.spde2', created and saved by :
  # spd <- sp::SpatialPointsDataFrame(coords = blackcap[, c("longitude", "latitude")],
  #                            data = blackcap)
  # small_mesh <- INLA::inla.mesh.2d(loc = INLA::inla.mesh.map(sp::coordinates(spd)), 
  #                           max.n=100, # only for demonstration purposes
  #                           max.edge = c(3, 20)) 
  # small_spde <- INLA::inla.spde2.matern(small_mesh)
  # save(small_spde, file="small_spde.RData", version=2)
#  
fit_SPDE <- fitme(migStatus ~ means + IMRF(1|longitude+latitude, model=small_spde), 
                  data=blackcap)
                  
########################## Regular lattices:   
# 
#Using 'hyper' to control fixed hyper-parameters
#
(mrf <- fitme(migStatus ~ 1 + (1|pos) + 
                          multIMRF(1|longitude+latitude,margin=5,levels=2), 
              data=blackcap, fixed=list(phi=1,lambda=c("1"=0.5),
              hyper=list("1"=list(hy_kap=0.1,hy_lam=1)))) )
              
# Using 'hyper' to control initial hyper-parameters
#
(mrf <- fitme(migStatus ~ 1 + multIMRF(1|longitude+latitude,margin=5,levels=2),
                data=blackcap, method="ML", fixed =list(phi=1),
                init=list(hyper=list("1"=list(hy_kap=0.1,hy_lam=1)))) )
                
# *Independent* IMRF terms with default rectangular lattice (often giving dubious results)
#
(mrf <- fitme(migStatus ~ 1 + IMRF(1|longitude+latitude,margin=5, nd=4L)
                              + IMRF(1|longitude+latitude,margin=5, nd=7L),
          data=blackcap,  
          fixed=list(phi=1,lambda=c(1/4,1/16),
                       corrPars=list("1"=list(kappa=0.1),"2"=list(kappa=0.1)))))
                    
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }% use one of  RShowDoc("KEYWORDS")
\keyword{ spatial }% __ONLY ONE__ keyword per line
