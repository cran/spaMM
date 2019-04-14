\name{ZAXlist}
\alias{ZAXlist}
\alias{class:ZAXlist}
\alias{ZAXlist-class}
\alias{missingOrNULL}
\alias{class:missingOrNULL}
\alias{missingOrNULL-class}
\docType{methods}
\alias{\%*\%-methods}
\alias{crossprod-methods}
\alias{tcrossprod-methods}
\alias{\%*\%,ZAXlist,Matrix-method}
\alias{\%*\%,ZAXlist,matrix-method}
\alias{\%*\%,ZAXlist,numeric-method}
\alias{\%*\%,numeric,ZAXlist-method}
\alias{crossprod,ZAXlist,Matrix-method}
\alias{crossprod,ZAXlist,matrix-method}
\alias{crossprod,ZAXlist,numeric-method}
\alias{tcrossprod,ZAXlist,missingOrNULL-method}
\title{ZAXlist class and (cross) products for ZAL matrix}
\description{
A \code{ZAXlist} object is a representation of the \dQuote{\code{ZAL}} matrix as a list of descriptors of each \code{ZAL} block for each random effect.

This documentation is for development purposes and may be incomplete. The objects and methods are not part of the programming interface and are subject to modification without notice.
}
\usage{
# new("ZAXlist", LIST=.) 
}
% perhaps define a constructor ?
\section{Slots}{
  \describe{
    \item{\code{LIST}:}{A list whose each block is either a \code{(M|m)atrix}, 
or a list with two elements (and addtional class \code{ZA_QCHM}): \code{ZA}, and the \code{\link[Matrix]{Cholesky}} factor \code{Q_CHMfactor} of the precision matrix (\code{L=solve(Q_CHMfactor,system="Lt")}).}
  }
}