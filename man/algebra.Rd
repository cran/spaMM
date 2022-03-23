\name{algebra}
\alias{algebra}
\alias{sparse_precision}

\title{Control of matrix-algebraic methods}
\description{
Autocorrelated gaussian random effects can be specified in terms of their covariance matrix, or in terms of the precision matrix (i.e. inverse covariance matrix). In a pre-processing step, spaMM may assess whether such precision matrices are sparse but the correlation matrix is dense, and if so, it may use \dQuote{sparse-precision} algorithms efficient for this case. If the precision matrix does not appear sufficiently sparser than the correlation matrix, correlation matrices are used, and they can themselves be sparse or dense, with distinct algebraic methods used in each case. 

For example, when the model includes a corrMatrix term specified by a covariance matrix, the precision matrix may be computed to assess its sparseness. The Example below illustrates a case where detecting sparsity of the precision matrix allows a faster fit. However, such a comparison of correlation and precision matrices takes time and is not performed for all types of random-effect structures. Instead, some fast heuristics may be used (see Details). The default selection of methods may not always be optimal, and may be overcome by using the \code{control.HLfit} argument of the fitting function (or by \code{spaMM.options()}, see Details). In particular one can use either \code{control.HLfit=list(sparse_precision= <TRUE|FALSE>)} or 
\code{control.HLfit=list(algebra= <"spprec"|"spcorr"|"decorr">)} with the obvious expected effects.

Such control may be useful when you already know that the precision matrix is sparse (as spaMM may even kindly remind you of, see Example below). In that case, it is also efficient to specify the precision matrix directly (see Example in \code{\link{Gryphon}}), as spaMM then assumes that sparse-precision methods are better without checking the correlation matrix.

Such control may also be useful when the correlation matrix is nearly singular so that computation of its inverse fails. This may occur if the model is poorly specified, but also occurs sometimes for valid correlation models because inversion of large matrices though Cholesky methods is not numerically accurate enough. In the latter case, you may be directed to this documentation by an error message, and specifying \code{sparse_precision= FALSE} may be useful.

}
\details{
Currently the sparse-precision methods are selected by default in two cases (with possible exceptions indicated by specific messages): (1) for models including \code{\link{IMRF}} random effects; and (2) when the \code{\link{corrMatrix}} (or \code{\link{covStruct}}) syntax is used to provide a fixed precision matrix. Further, for models including autoregressive terms other than IMRF (i.e., adjacency, AR1), sparse-precision methods may or may not be selected on a simple heuristic based on the likely structure of the correlation matrix. 

Algebraic methods can be controlled globally over all further fits by using\cr 
\code{spaMM.options(sparse_precision= <TRUE|FALSE>)}\cr 
and, among the correlation-based methods,\cr 
\code{spaMM.options(QRmethod= <"sparse"|"dense">)}\cr 
to select \code{"spcorr"} vs. \code{"decorr"} methods. Fit-specific controls (by \code{control.HLfit}) override these global ones.
}
\seealso{\code{\link{pedigree}}}
\examples{
if (spaMM.getOption("example_maxtime")>6) {
 data("Gryphon")
 
 gry_df <- fitme(BWT ~ 1 + corrMatrix(1|ID), corrMatrix = Gryphon_A, 
                 data = Gryphon_df, method = "REML")
 how(gry_df)
               
 # => Note the message about 'Choosing matrix methods...'. 
 # Using control.HLfit=list(algebra="spprec") would indeed 
 # save the time used to select this method.
 
 # Conversely, using a correlation-based method would be a waste of time:
 
 gry_dn <- fitme(BWT ~ 1 + corrMatrix(1|ID), corrMatrix = Gryphon_A, 
               data = Gryphon_df, method = "REML",
               control.HLfit=list(sparse_precision=FALSE))
 how(gry_dn) # forced dense-correlation methods, which is slower here.
}

}
