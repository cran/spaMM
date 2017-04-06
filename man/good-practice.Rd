\name{good-practice}
\alias{good-practice}

\title{Clear and trustworthy formulas}

\description{
Base fitting functions in R will seek variables in the global environment (or more generally in the environment where a call to \code{`~`} was made, defining the model formula) if they are not in the \code{data}. This easily leads to errors (see example in the discussion of \code{\link{update.HLfit}}). Indeed Chambers (2008, p.221), after describing how the environment is defined, comments that \dQuote{Where clear and trustworthy software is a priority, I would personally avoid such tricks. Ideally, all the variables in the model frame should come from an explicit, verifiable data source...}. Hence, the main fitting functions in \code{spaMM} depart from the sloppy practice. They strip the formula environment from any variable, and seek all variables from the \code{formula} in the data frame given by their \code{data} argument. \bold{One never needs to specify the \code{data} in the \code{formula}.} The variables defining the \code{prior.weights} should also be in the \code{data}. Variables used in other arguments such as \code{ranFix} are looked up neither in the data nor in the formula environment, but in the calling environment as usual.
}
\references{
Chambers J.M. (2008) Software for data analysis: Programming with R. Springer-Verlag New York
}