\name{PLS-internals}
\alias{PLS-internals}

%\alias{Rcpp_chol_R}

%\alias{lmwithQR}
%\alias{lmwithQRP}
%\alias{lmwith_sparse_QRp}
%\alias{lmwith_sparse_LDLp}
%\alias{lmwith_sparse_LLp}
%\alias{sparse_LDLp_from_XtX}
%\alias{sparse_LLp_from_XtX}
%\alias{LevMar_cpp}
%\alias{fit_as_ZX}
%\alias{fit_as_sparsePrecision}
%\alias{calc_APHLs_from_ZX}
%\alias{calc_sXaug_Re}
%\alias{create_sXaug_EigenDense_QR_scaled}
%\alias{create_sXaug_baseDense_QRP_scaled}
%\alias{create_sXaug_EigenSparse_QR_scaled}
%\alias{create_sXaug_EigenSparse_QRP_scaled}
%\alias{create_sXaug_EigenSparse_LDLP_scaled}
%\alias{create_sXaug_EigenSparse_LDL_scaled}
%\alias{create_sXaug_EigenSparse_LLP_scaled}
\alias{def_AUGI0_ZX_sparsePrecision}
%\alias{def_sXaug_EigenDense_QR_scaled}
%\alias{def_sXaug_EigenDense_QRP_scaled}
\alias{def_sXaug_EigenDense_QRP_Chol_scaled}
%\alias{def_sXaug_baseDense_QRP_scaled}
%\alias{def_sXaug_EigenSparse_QR_scaled}
%\alias{def_sXaug_EigenSparse_QRP_scaled}
%\alias{def_sXaug_EigenSparse_LDLP_scaled}
%\alias{def_sXaug_EigenSparse_LDL_scaled}
%\alias{def_sXaug_EigenSparse_LLP_scaled}
%\alias{def_sXaug_Matrix_QRP_scaled}
%\alias{def_sXaug_Matrix_cholP_scaled}
\alias{def_sXaug_Matrix_QRP_CHM_scaled}
\alias{get_from_MME}
\alias{get_from_MME.AUGI0_ZX_sparsePrecision}
%\alias{get_from_MME.sXaug_EigenDense_QR_scaled}
%\alias{get_from_MME.sXaug_EigenDense_QRP_scaled}
\alias{get_from_MME.sXaug_EigenDense_QRP_Chol_scaled}
%\alias{get_from_MME.sXaug_baseDense_QRP_scaled}
%\alias{get_from_MME.sXaug_EigenSparse_QR_scaled}
%\alias{get_from_MME.sXaug_EigenSparse_QRP_scaled}
%\alias{get_from_MME.sXaug_EigenSparse_LDLP_scaled}
%\alias{get_from_MME.sXaug_EigenSparse_LDL_scaled}
%\alias{get_from_MME.sXaug_EigenSparse_LLP_scaled}
%\alias{get_from_MME.sXaug_Matrix_QRP_scaled}
%\alias{get_from_MME.sXaug_Matrix_cholP_scaled}
\alias{get_from_MME.sXaug_Matrix_QRP_CHM_scaled}
\alias{get_from_MME.default}
\alias{get_from_MME.sparseMatrix}
\alias{get_from_MME_default}
\alias{get_from_MME_default.matrix}
\alias{get_from_MME_default.Matrix}
%\alias{make_Xscal}

\title{Internal, Experimental Functions for procedure using the ((I,0),(Z,X)) block-order}
\description{
  Internal spaMM functions
}
\details{
These functions use the block order introduced by Bates and DebRoy (2004) and also used in \code{lme4}.
  These are not to be called by the user, or are waiting for documentation to be written.  
}
\references{
Bates D.M., DebRoy S. Linear mixed models and penalized least squares. Journal of Multivariate Analysis 91: 1–17.
}
\keyword{ internal }
