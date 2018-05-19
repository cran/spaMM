\name{separation}
\alias{separation}

\title{Checking separation in binomial-response models}

\description{Separation occurs in binomial response models when a combination of the predictor variables perfectly predict a level of the response. In such a case the estimates of the coefficients for these variables diverge to (+/-)infinity, and the numerical algorithms typically fail. To anticipate such a problem, the fitting functions in \code{spaMM} try to check for separation by default, using the \code{lpSolveAPI} package which can also detect some borderline cases (\dQuote{quasi-separation}). It this package is not available, \code{spaMM} tries to use the \code{e1071} package (in a way which will not detect quasi-separation), except for large data sets where this may take time, in which case a message notifies this to the user. The threshold size of the data is set by \code{spaMM.options(separation_max=<.>)}}

