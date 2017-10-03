\name{sparse_precision}
\alias{sparse_precision}

\title{Sparse_precision algorithm}
\description{
A fitting algorithm efficient for random effects with sparse precision matrix (i.e. inverse covariance matrix) is implemented. It is used by default only in two cases: for conditional autoregressive models (with a random effect of the form \code{adjacency(1|<grouping factor>)}), and when the \code{\link{covStruct}} syntax is used to provide a fixed precision matrix (see \code{\link{pedigree}} for an example). A non-default choice of fitting algorithm can be selected in this and other models by using \code{spaMM.options(sparse_precision= <TRUE|FALSE>)} with often poor results.\cr
A precision matrix is meaningful for a Gaussian random effect, but beyond this the algorithm works for HGLMs, i.e. the model may include another random effect with non-Gaussian distribution.
}

