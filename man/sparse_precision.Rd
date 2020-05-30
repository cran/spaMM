\name{sparse_precision}
\alias{sparse_precision}

\title{Sparse_precision algorithm}
\description{
spaMM includes fitting algorithms formulated in terms of the correlation matrix of random effects, and algorithms formulated in terms of the precision matrix (i.e. inverse covariance matrix) for random effects.  Formulation of random effects in terms of their precision or of their correlation matrices are meaningful for Gaussian random effects, but beyond this both classes of algorithms work when the model include random effects with non-Gaussian distribution and no intrinsic correlation structure, and therefore for the full class of HGLMs. 

The algorithms based on precision matrices may be more efficient when the precision matrix is sparse but the correlation matrix is dense. However, spaMM does not yet select the fastest algorithm by default, and the default choice has changed over versions without being properly documented here. A non-default choice of fitting algorithm can be selected by using \code{spaMM.options(sparse_precision= <TRUE|FALSE>)}. Currently it is selected by default in two cases (with exceptions indicated by specific messages): (1) for \code{\link{IMRF}} random effects, but not for other conditional autoregressive models (with a random effect of the form \code{adjacency(1|<grouping factor>)}); and (2) when the \code{\link{covStruct}} syntax is used to provide a fixed precision matrix (see \code{\link{pedigree}} for an example). 
}

