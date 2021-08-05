\name{sparse_precision}
\alias{sparse_precision}

\title{Control of usage of sparse_precision methods}
\description{
Autocorrelated gaussian random effects can be specified in terms of their covariance matrix, or in terms of the precision matrix (i.e. inverse covariance matrix). In a pre-processing step, spaMM may assess whether such precision matrices are sparse but the correlation matrix is dense, and if so, it may use \dQuote{sparse-precision} algorithms efficient for this case. However this assessment takes times and is not systematically performed. Instead, some fast heuristics are used (see Details).

The default selection (or not) of sparse-precision methods may not always be optimal, and may be overcome by using \code{control.HLfit=list(sparse_precision= <TRUE|FALSE>)} (\code{control.HLfit} argument of fitting function) or \code{spaMM.options(sparse_precision= <TRUE|FALSE>)}. This control may also be useful in other cases, notably when a large \code{corrMatrix} is specified and either (1) you know whether the precision matrix will be sparse (in that case, it might even be better to specify the precision matrix directly), or dense; 
or (2) the correlation matrix is nearly singular so that computation of its inverse fails. This may occur if the model is poorly specified, but also occurs sometimes for valid correlation models because inversion of large matrices though Cholesky methods is not numerically accurate enough. In the latter case, you may be directed to this documentation by an error message, and specifying \code{sparse_precision= FALSE} may be useful.
}
\details{Currently the sparse-precision methods are selected by default in two cases (with possible exceptions indicated by specific messages): (1) for models including \code{\link{IMRF}} random effects; and (2) when the \code{\link{covStruct}} syntax is used to provide a fixed precision matrix (see \code{\link{pedigree}} for an example). Further, for models including autoregressive terms other than IMRF (i.e., adjacency, AR1), sparse-precision methods may or may not be selected on a simple heuristic based on the likely structure of the correlation matrix. When the model includes a corrMatrix term specified by a covariance matrix, the (costly) computation of the precision matrix is performed to assess its sparseness. If it is sparse enough, the corrMatrix term is handled simlarly to an adjacency one. 
The Example illustrates a case where detecting sparsity of the precision matrix deduced from a corrMatrix allows a faster fit.
}
\seealso{\code{\link{pedigree}}}
\examples{
if (spaMM.getOption("example_maxtime")>6) {
 data("Gryphon")
 gry_df <- fitme(BWT ~ 1 + corrMatrix(1|ID), corrMatrix = Gryphon_A, 
               data = Gryphon_df, method = "REML")
 how(gry_df) # Shows default fitting method
 gry_dn <- fitme(BWT ~ 1 + corrMatrix(1|ID), corrMatrix = Gryphon_A, 
               data = Gryphon_df, method = "REML",
               control.HLfit=list(sparse_precision=FALSE))
 how(gry_dn) # forced dense-correlation methods, which is slower here.
}

}
