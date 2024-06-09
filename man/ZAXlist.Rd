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
\alias{t.ZAXlist}
\alias{crossprod,ZAXlist,Matrix-method}
\alias{crossprod,ZAXlist,matrix-method}
\alias{crossprod,ZAXlist,numeric-method}
\alias{tcrossprod,ZAXlist,missingOrNULL-method}
\alias{Kronfacto}
\alias{class:Kronfacto}
\alias{Kronfacto-class}
\alias{dim.Kronfacto}
%\alias{\%*\%,Kronfacto,Matrix-method}
%\alias{\%*\%,Kronfacto,matrix-method}
\alias{\%*\%,Kronfacto,numeric-method}
%\alias{\%*\%,numeric,ZAXlist-method}
\alias{crossprod,Kronfacto,Matrix-method}
\alias{crossprod,Kronfacto,matrix-method}
\alias{crossprod,Kronfacto,numeric-method}

\title{S4 classes for structured matrices}
\description{
A \code{ZAXlist} object is a representation of the \dQuote{\code{ZAL}} matrix as an S4 class holding a list of descriptors of each \code{ZAL} block for each random effect.

A \code{Kronfacto} object is a representation of a Kronecker product as an S4 class holding its factors. Methods defined for this class may avoid the computation of the Kronecker product as an actual matrix of large dimensions.

This documentation is for development purposes and may be incomplete. The objects and methods are not part of the programming interface and are subject to modification without notice.
}
\usage{
# new("ZAXlist", LIST=.) 
# new("Kronfacto", BLOB=.) 
}
% perhaps define a constructor ?
\section{Slots}{
  \describe{
    \item{\code{LIST}:}{A list whose each element is a ZAL block represented as either a \code{(M|m)atrix}, 
or a list with two elements (and additional class \code{ZA_QCHM}): \code{ZA}, and the \code{\link[Matrix]{Cholesky}} factor \code{Q_CHMfactor} of the precision matrix (\code{L=solve(Q_CHMfactor,system="Lt")}).}
    \item{\code{BLOB}:}{An environment holding \code{lhs} and \code{rhs}, the factors of the Kronecker product, and other objects initialized as promises. See the source code of the non-exported \code{.def_Kranfacto} constructor for further information.}
  }
}
%\seealso{
%See also source of spaMM:::.def_Kronfacto(lhs,rhs)
%}
