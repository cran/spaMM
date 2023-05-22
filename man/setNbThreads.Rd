\name{setNbThreads}
\alias{setNbThreads}
\title{
Parallel computations in fits
}
\description{
A few steps of fitting can be parallelized. Currently it is possible to control the use of multiple threads by OpenMP by the Eigen library. By default only one thread will be used, but this may be modified by using \code{control.HLfit$NbThreads} in a fitting function's arguments, as in\cr
\preformatted{
avail_thr <- parallel::detectCores() - 1L 
fitme(., control.HLfit=list(NbThreads=max(avail_thr, 1L)))
}
This control is distinct from that of post-fit steps such as bootstraps where some parallel computations are controlled e.g. the \code{nb_cores} argument of \code{spaMM_boot}. In cases where post-fits computation imply refits of models (as is typical of parametric bootstraps), the two parallelizations should not be combined, and the \pkg{spaMM} code for post-fit operations will in principle automatically take care of this.  

According to \url{https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#OpenMP-support}, using openMP may decrease the precision of some computations, and may be less efficient under Windows; and according to \url{https://eigen.tuxfamily.org/dox/TopicMultiThreading.html} only a few Eigen computations will benefit from such parallelisation, mainly the dense matrix products. \pkg{spaMM} will \emph{suggest} using parallelisation when random effects have many levels and dense-correlation methods are selected (see \code{\link{algebra}}), that is mainly for geostatiscal models with many locations. Speed gains appear moderate, as the slowest steps are not parallelized.      
}
