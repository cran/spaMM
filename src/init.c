#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bessel_lnKnu_e(void *, void *, void *, void *, void *, void *);
extern void matern_cor(void *, void *, void *, void *, void *);
extern void matern_factList(void *, void *, void *, void *, void *, void *, void *, void *);
extern void matern_matList(void *, void *, void *, void *, void *, void *, void *);
extern void matern_recalc(void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP spaMM_crossprodCpp(SEXP, SEXP);
extern SEXP spaMM_e_LevenbergMsolveCpp(SEXP, SEXP, SEXP);
extern SEXP spaMM_LevenbergMsolveCpp(SEXP, SEXP, SEXP);
extern SEXP spaMM_leverages(SEXP);
extern SEXP spaMM_LevMar_cpp(SEXP, SEXP);
extern SEXP spaMM_LevPerturbedQCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP spaMM_lmwith_sparse_LDLp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP spaMM_lmwith_sparse_LLp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP spaMM_lmwith_sparse_QRp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP spaMM_lmwithQR(SEXP, SEXP, SEXP, SEXP);
extern SEXP spaMM_lmwithQRP(SEXP, SEXP, SEXP, SEXP);
extern SEXP spaMM_LogAbsDetCpp(SEXP);
//extern SEXP spaMM_pseudoSolvediag(SEXP, SEXP);
extern SEXP spaMM_Rcpp_chol_R(SEXP);
extern SEXP spaMM_Rcpp_COMP_Z(SEXP, SEXP, SEXP, SEXP);
extern SEXP spaMM_Rcpp_d2hdv2(SEXP, SEXP, SEXP);
extern SEXP spaMM_Rcpp_Sig(SEXP, SEXP, SEXP);
extern SEXP spaMM_RcppChol(SEXP);
extern SEXP spaMM_selfAdjointSolverCpp(SEXP);
extern SEXP spaMM_sparse_LDLp_from_XtX(SEXP, SEXP);
extern SEXP spaMM_sparse_LLp_from_XtX(SEXP, SEXP);
extern SEXP spaMM_sweepZ1W(SEXP, SEXP);
extern SEXP spaMM_tcrossprodCpp(SEXP, SEXP);
extern SEXP spaMM_ZtWZ(SEXP, SEXP);
extern SEXP spaMM_ZWZt(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"bessel_lnKnu_e",  (DL_FUNC) &bessel_lnKnu_e,  6},
  {"matern_cor",      (DL_FUNC) &matern_cor,      5},
  {"matern_factList", (DL_FUNC) &matern_factList, 8},
  {"matern_matList",  (DL_FUNC) &matern_matList,  7},
  {"matern_recalc",   (DL_FUNC) &matern_recalc,   9},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"spaMM_crossprodCpp",         (DL_FUNC) &spaMM_crossprodCpp,         2},
  {"spaMM_e_LevenbergMsolveCpp", (DL_FUNC) &spaMM_e_LevenbergMsolveCpp, 3},
  {"spaMM_LevenbergMsolveCpp",   (DL_FUNC) &spaMM_LevenbergMsolveCpp,   3},
  {"spaMM_leverages",            (DL_FUNC) &spaMM_leverages,            1},
  {"spaMM_LevMar_cpp",           (DL_FUNC) &spaMM_LevMar_cpp,           2},
  {"spaMM_LevPerturbedQCpp",     (DL_FUNC) &spaMM_LevPerturbedQCpp,     6},
  {"spaMM_lmwith_sparse_LDLp",   (DL_FUNC) &spaMM_lmwith_sparse_LDLp,   5},
  {"spaMM_lmwith_sparse_LLp",    (DL_FUNC) &spaMM_lmwith_sparse_LLp,    5},
  {"spaMM_lmwith_sparse_QRp",    (DL_FUNC) &spaMM_lmwith_sparse_QRp,    5},
  {"spaMM_lmwithQR",             (DL_FUNC) &spaMM_lmwithQR,             4},
  {"spaMM_lmwithQRP",            (DL_FUNC) &spaMM_lmwithQRP,            4},
  {"spaMM_LogAbsDetCpp",         (DL_FUNC) &spaMM_LogAbsDetCpp,         1},
//  {"spaMM_pseudoSolvediag",      (DL_FUNC) &spaMM_pseudoSolvediag,      2},
  {"spaMM_Rcpp_chol_R",          (DL_FUNC) &spaMM_Rcpp_chol_R,          1},
  {"spaMM_Rcpp_COMP_Z",          (DL_FUNC) &spaMM_Rcpp_COMP_Z,          4},
  {"spaMM_Rcpp_d2hdv2",          (DL_FUNC) &spaMM_Rcpp_d2hdv2,          3},
  {"spaMM_Rcpp_Sig",             (DL_FUNC) &spaMM_Rcpp_Sig,             3},
  {"spaMM_RcppChol",             (DL_FUNC) &spaMM_RcppChol,             1},
  {"spaMM_selfAdjointSolverCpp", (DL_FUNC) &spaMM_selfAdjointSolverCpp, 1},
  {"spaMM_sparse_LDLp_from_XtX", (DL_FUNC) &spaMM_sparse_LDLp_from_XtX, 2},
  {"spaMM_sparse_LLp_from_XtX",  (DL_FUNC) &spaMM_sparse_LLp_from_XtX,  2},
  {"spaMM_sweepZ1W",             (DL_FUNC) &spaMM_sweepZ1W,             2},
  {"spaMM_tcrossprodCpp",        (DL_FUNC) &spaMM_tcrossprodCpp,        2},
  {"spaMM_ZtWZ",                 (DL_FUNC) &spaMM_ZtWZ,                 2},
  {"spaMM_ZWZt",                 (DL_FUNC) &spaMM_ZWZt,                 2},
  {NULL, NULL, 0}
};

void R_init_spaMM(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
