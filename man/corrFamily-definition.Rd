\name{corrFamily-definition}
\alias{corrFamily-definition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
corrFamily definition
}
\description{
Tentative formal rules for definition of corrFamily descriptors (work in progress). This is likely to repeat and extend information partially given in \code{\link{corrFamily}} and \code{\link{corrFamily-design}} documentations.
\describe{

User-level rules (not relevant fo  corrFamily descriptors internally modified during a fit):
  \item{tpar}{
  Should always be present. For trivial parameterless cases (e.g. \code{ranGCA}), it should be \code{numeric(0L)}, not \code{NULL}.
}
  \item{Cf}{
  function; should always be present. For trivial uncorrelated random effects (e.g. \code{ranGCA}, 
  where only the \code{Af} function carries the information for the model), 
  it should return an identity matrix, not \code{NULL}, % can this be avoided ? That woudl be perhaps awkward (cf unused level removal - rownames generally expected)
  with row names to be matched to the column names of the \bold{Z} matrix for the random effect.
}
\item{calc_moreargs}{
optional function. If present, it should have formal arguments including at least \code{corrfamily} and \ldots. 
}
  \item{Af}{
  function; optional. If present, it should have row names to be matched to the column names of the \bold{Z} matrix for the random effect, and also needs column names if it is to be matched with the row names of a correlation matrix (or its inverse).
}
  \item{initialize}{
Optional function.  If present, should have formal arguments including at least \code{Zmatrix} and \ldots.  

In predefined corrFamily constructors, variables created by \code{initialize} for use by \code{Cf} or \code{Af} should be declared (typically as \code{NULL}) in the body of the constructor, so that R CMD check does not complain.
}
\item{public}{
An environment. \code{initialize} may write into it. It might also read into it, for example read the result of a long previous computation by \code{initialize} during a previous fit, though this opens the door to various errors. 
}

}

}
