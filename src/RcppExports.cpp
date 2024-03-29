// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rcpp_COMP_Z
SEXP Rcpp_COMP_Z(int moment, double nu, double lambda, int maxn);
RcppExport SEXP _spaMM_Rcpp_COMP_Z(SEXP momentSEXP, SEXP nuSEXP, SEXP lambdaSEXP, SEXP maxnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type moment(momentSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type maxn(maxnSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_COMP_Z(moment, nu, lambda, maxn));
    return rcpp_result_gen;
END_RCPP
}
// COMP_Z_integrand
SEXP COMP_Z_integrand(Rcpp::NumericVector z, Nullable<NumericVector> eta, Nullable<NumericVector> lambda, double nu, int moment, double logScaleFac);
RcppExport SEXP _spaMM_COMP_Z_integrand(SEXP zSEXP, SEXP etaSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP momentSEXP, SEXP logScaleFacSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type moment(momentSEXP);
    Rcpp::traits::input_parameter< double >::type logScaleFac(logScaleFacSEXP);
    rcpp_result_gen = Rcpp::wrap(COMP_Z_integrand(z, eta, lambda, nu, moment, logScaleFac));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_COMP_Z_asympto
SEXP Rcpp_COMP_Z_asympto(double nu, double pow_lam_nu);
RcppExport SEXP _spaMM_Rcpp_COMP_Z_asympto(SEXP nuSEXP, SEXP pow_lam_nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type pow_lam_nu(pow_lam_nuSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_COMP_Z_asympto(nu, pow_lam_nu));
    return rcpp_result_gen;
END_RCPP
}
// lmwith_sparse_LDLp
SEXP lmwith_sparse_LDLp(SEXP XX, SEXP yy, bool returntQ, bool returnR, bool pivot);
RcppExport SEXP _spaMM_lmwith_sparse_LDLp(SEXP XXSEXP, SEXP yySEXP, SEXP returntQSEXP, SEXP returnRSEXP, SEXP pivotSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yy(yySEXP);
    Rcpp::traits::input_parameter< bool >::type returntQ(returntQSEXP);
    Rcpp::traits::input_parameter< bool >::type returnR(returnRSEXP);
    Rcpp::traits::input_parameter< bool >::type pivot(pivotSEXP);
    rcpp_result_gen = Rcpp::wrap(lmwith_sparse_LDLp(XX, yy, returntQ, returnR, pivot));
    return rcpp_result_gen;
END_RCPP
}
// lmwith_sparse_LLp
SEXP lmwith_sparse_LLp(SEXP XX, SEXP yy, bool returntQ, bool returnR, bool pivot);
RcppExport SEXP _spaMM_lmwith_sparse_LLp(SEXP XXSEXP, SEXP yySEXP, SEXP returntQSEXP, SEXP returnRSEXP, SEXP pivotSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yy(yySEXP);
    Rcpp::traits::input_parameter< bool >::type returntQ(returntQSEXP);
    Rcpp::traits::input_parameter< bool >::type returnR(returnRSEXP);
    Rcpp::traits::input_parameter< bool >::type pivot(pivotSEXP);
    rcpp_result_gen = Rcpp::wrap(lmwith_sparse_LLp(XX, yy, returntQ, returnR, pivot));
    return rcpp_result_gen;
END_RCPP
}
// update_R_in_place
SEXP update_R_in_place(SEXP RD);
RcppExport SEXP _spaMM_update_R_in_place(SEXP RDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type RD(RDSEXP);
    rcpp_result_gen = Rcpp::wrap(update_R_in_place(RD));
    return rcpp_result_gen;
END_RCPP
}
// lmwithQR
SEXP lmwithQR(SEXP XX, SEXP yy, bool returntQ, bool returnR);
RcppExport SEXP _spaMM_lmwithQR(SEXP XXSEXP, SEXP yySEXP, SEXP returntQSEXP, SEXP returnRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yy(yySEXP);
    Rcpp::traits::input_parameter< bool >::type returntQ(returntQSEXP);
    Rcpp::traits::input_parameter< bool >::type returnR(returnRSEXP);
    rcpp_result_gen = Rcpp::wrap(lmwithQR(XX, yy, returntQ, returnR));
    return rcpp_result_gen;
END_RCPP
}
// lmwithQRP
SEXP lmwithQRP(SEXP XX, SEXP yy, bool returntQ, bool returnR);
RcppExport SEXP _spaMM_lmwithQRP(SEXP XXSEXP, SEXP yySEXP, SEXP returntQSEXP, SEXP returnRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yy(yySEXP);
    Rcpp::traits::input_parameter< bool >::type returntQ(returntQSEXP);
    Rcpp::traits::input_parameter< bool >::type returnR(returnRSEXP);
    rcpp_result_gen = Rcpp::wrap(lmwithQRP(XX, yy, returntQ, returnR));
    return rcpp_result_gen;
END_RCPP
}
// lmwith_sparse_QRp
SEXP lmwith_sparse_QRp(SEXP XX, SEXP yy, bool returntQ, bool returnR, bool COLAMDO);
RcppExport SEXP _spaMM_lmwith_sparse_QRp(SEXP XXSEXP, SEXP yySEXP, SEXP returntQSEXP, SEXP returnRSEXP, SEXP COLAMDOSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yy(yySEXP);
    Rcpp::traits::input_parameter< bool >::type returntQ(returntQSEXP);
    Rcpp::traits::input_parameter< bool >::type returnR(returnRSEXP);
    Rcpp::traits::input_parameter< bool >::type COLAMDO(COLAMDOSEXP);
    rcpp_result_gen = Rcpp::wrap(lmwith_sparse_QRp(XX, yy, returntQ, returnR, COLAMDO));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_chol_R
SEXP Rcpp_chol_R(SEXP AA);
RcppExport SEXP _spaMM_Rcpp_chol_R(SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_chol_R(AA));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_dense_cbind_mat_mat
NumericMatrix Rcpp_dense_cbind_mat_mat(NumericMatrix a, NumericMatrix b);
RcppExport SEXP _spaMM_Rcpp_dense_cbind_mat_mat(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_dense_cbind_mat_mat(a, b));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_dense_cbind_mat_vec
NumericMatrix Rcpp_dense_cbind_mat_vec(NumericMatrix a, NumericVector b);
RcppExport SEXP _spaMM_Rcpp_dense_cbind_mat_vec(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_dense_cbind_mat_vec(a, b));
    return rcpp_result_gen;
END_RCPP
}
// RcppMatrixCb2
Eigen::SparseMatrix<double> RcppMatrixCb2(Eigen::MappedSparseMatrix<double>& matrix1, Eigen::MappedSparseMatrix<double>& matrix2);
RcppExport SEXP _spaMM_RcppMatrixCb2(SEXP matrix1SEXP, SEXP matrix2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MappedSparseMatrix<double>& >::type matrix1(matrix1SEXP);
    Rcpp::traits::input_parameter< Eigen::MappedSparseMatrix<double>& >::type matrix2(matrix2SEXP);
    rcpp_result_gen = Rcpp::wrap(RcppMatrixCb2(matrix1, matrix2));
    return rcpp_result_gen;
END_RCPP
}
// lowertri
SEXP lowertri(NumericMatrix A, bool diag, NumericVector value);
RcppExport SEXP _spaMM_lowertri(SEXP ASEXP, SEXP diagSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< bool >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(lowertri(A, diag, value));
    return rcpp_result_gen;
END_RCPP
}
// uppertri
SEXP uppertri(NumericMatrix A, bool diag, NumericVector value);
RcppExport SEXP _spaMM_uppertri(SEXP ASEXP, SEXP diagSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< bool >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(uppertri(A, diag, value));
    return rcpp_result_gen;
END_RCPP
}
// rC_inv_chol_cpp
SEXP rC_inv_chol_cpp(SEXP trCoef);
RcppExport SEXP _spaMM_rC_inv_chol_cpp(SEXP trCoefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type trCoef(trCoefSEXP);
    rcpp_result_gen = Rcpp::wrap(rC_inv_chol_cpp(trCoef));
    return rcpp_result_gen;
END_RCPP
}
// C_calc_cov_from_ranCoef
SEXP C_calc_cov_from_ranCoef(SEXP ranCoef, Nullable<IntegerVector> Xi_ncol);
RcppExport SEXP _spaMM_C_calc_cov_from_ranCoef(SEXP ranCoefSEXP, SEXP Xi_ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ranCoef(ranCoefSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type Xi_ncol(Xi_ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(C_calc_cov_from_ranCoef(ranCoef, Xi_ncol));
    return rcpp_result_gen;
END_RCPP
}
// makelong
SEXP makelong(NumericMatrix Lcompact, int Lsize);
RcppExport SEXP _spaMM_makelong(SEXP LcompactSEXP, SEXP LsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Lcompact(LcompactSEXP);
    Rcpp::traits::input_parameter< int >::type Lsize(LsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(makelong(Lcompact, Lsize));
    return rcpp_result_gen;
END_RCPP
}
// makelong2
SEXP makelong2(NumericMatrix Lcompact, int Lsize);
RcppExport SEXP _spaMM_makelong2(SEXP LcompactSEXP, SEXP LsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Lcompact(LcompactSEXP);
    Rcpp::traits::input_parameter< int >::type Lsize(LsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(makelong2(Lcompact, Lsize));
    return rcpp_result_gen;
END_RCPP
}
// C_dispInv
SEXP C_dispInv(NumericVector x);
RcppExport SEXP _spaMM_C_dispInv(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_dispInv(x));
    return rcpp_result_gen;
END_RCPP
}
// logit
SEXP logit(NumericVector mu);
RcppExport SEXP _spaMM_logit(SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(logit(mu));
    return rcpp_result_gen;
END_RCPP
}
// is_evaluated
bool is_evaluated(Symbol name, Environment env);
RcppExport SEXP _spaMM_is_evaluated(SEXP nameSEXP, SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Symbol >::type name(nameSEXP);
    Rcpp::traits::input_parameter< Environment >::type env(envSEXP);
    rcpp_result_gen = Rcpp::wrap(is_evaluated(name, env));
    return rcpp_result_gen;
END_RCPP
}
// is_promise2
bool is_promise2(Symbol name, Environment env);
RcppExport SEXP _spaMM_is_promise2(SEXP nameSEXP, SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Symbol >::type name(nameSEXP);
    Rcpp::traits::input_parameter< Environment >::type env(envSEXP);
    rcpp_result_gen = Rcpp::wrap(is_promise2(name, env));
    return rcpp_result_gen;
END_RCPP
}
// nuln_plus_bessel_lnKnu
NumericVector nuln_plus_bessel_lnKnu(Rcpp::NumericVector x, double nu);
RcppExport SEXP _spaMM_nuln_plus_bessel_lnKnu(SEXP xSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(nuln_plus_bessel_lnKnu(x, nu));
    return rcpp_result_gen;
END_RCPP
}
// set_thread_nbr
int set_thread_nbr(int thr);
RcppExport SEXP _spaMM_set_thread_nbr(SEXP thrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type thr(thrSEXP);
    rcpp_result_gen = Rcpp::wrap(set_thread_nbr(thr));
    return rcpp_result_gen;
END_RCPP
}
// rankinfo
SEXP rankinfo(SEXP x, SEXP tol);
RcppExport SEXP _spaMM_rankinfo(SEXP xSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(rankinfo(x, tol));
    return rcpp_result_gen;
END_RCPP
}
// leverages
SEXP leverages(SEXP XX);
RcppExport SEXP _spaMM_leverages(SEXP XXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type XX(XXSEXP);
    rcpp_result_gen = Rcpp::wrap(leverages(XX));
    return rcpp_result_gen;
END_RCPP
}
// sparse_cwiseprod
SEXP sparse_cwiseprod(SEXP AA, SEXP BB);
RcppExport SEXP _spaMM_sparse_cwiseprod(SEXP AASEXP, SEXP BBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    Rcpp::traits::input_parameter< SEXP >::type BB(BBSEXP);
    rcpp_result_gen = Rcpp::wrap(sparse_cwiseprod(AA, BB));
    return rcpp_result_gen;
END_RCPP
}
// sweepZ1W
SEXP sweepZ1W(SEXP ZZ, SEXP WW);
RcppExport SEXP _spaMM_sweepZ1W(SEXP ZZSEXP, SEXP WWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ZZ(ZZSEXP);
    Rcpp::traits::input_parameter< SEXP >::type WW(WWSEXP);
    rcpp_result_gen = Rcpp::wrap(sweepZ1W(ZZ, WW));
    return rcpp_result_gen;
END_RCPP
}
// ZWZt
SEXP ZWZt(SEXP ZZ, SEXP WW);
RcppExport SEXP _spaMM_ZWZt(SEXP ZZSEXP, SEXP WWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ZZ(ZZSEXP);
    Rcpp::traits::input_parameter< SEXP >::type WW(WWSEXP);
    rcpp_result_gen = Rcpp::wrap(ZWZt(ZZ, WW));
    return rcpp_result_gen;
END_RCPP
}
// ZtWZ
SEXP ZtWZ(SEXP ZZ, SEXP WW);
RcppExport SEXP _spaMM_ZtWZ(SEXP ZZSEXP, SEXP WWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ZZ(ZZSEXP);
    Rcpp::traits::input_parameter< SEXP >::type WW(WWSEXP);
    rcpp_result_gen = Rcpp::wrap(ZtWZ(ZZ, WW));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_d2hdv2
SEXP Rcpp_d2hdv2(SEXP ZZ, SEXP WA, SEXP WB);
RcppExport SEXP _spaMM_Rcpp_d2hdv2(SEXP ZZSEXP, SEXP WASEXP, SEXP WBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ZZ(ZZSEXP);
    Rcpp::traits::input_parameter< SEXP >::type WA(WASEXP);
    Rcpp::traits::input_parameter< SEXP >::type WB(WBSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_d2hdv2(ZZ, WA, WB));
    return rcpp_result_gen;
END_RCPP
}
// RcppChol
SEXP RcppChol(SEXP AA);
RcppExport SEXP _spaMM_RcppChol(SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(RcppChol(AA));
    return rcpp_result_gen;
END_RCPP
}
// crossprodCpp_d
SEXP crossprodCpp_d(SEXP Mat, SEXP yy);
RcppExport SEXP _spaMM_crossprodCpp_d(SEXP MatSEXP, SEXP yySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Mat(MatSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yy(yySEXP);
    rcpp_result_gen = Rcpp::wrap(crossprodCpp_d(Mat, yy));
    return rcpp_result_gen;
END_RCPP
}
// tcrossprodCpp
SEXP tcrossprodCpp(SEXP Mat, SEXP yy);
RcppExport SEXP _spaMM_tcrossprodCpp(SEXP MatSEXP, SEXP yySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Mat(MatSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yy(yySEXP);
    rcpp_result_gen = Rcpp::wrap(tcrossprodCpp(Mat, yy));
    return rcpp_result_gen;
END_RCPP
}
// LevenbergMsolveCpp
SEXP LevenbergMsolveCpp(SEXP AA, SEXP rrhhss, SEXP dd);
RcppExport SEXP _spaMM_LevenbergMsolveCpp(SEXP AASEXP, SEXP rrhhssSEXP, SEXP ddSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    Rcpp::traits::input_parameter< SEXP >::type rrhhss(rrhhssSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dd(ddSEXP);
    rcpp_result_gen = Rcpp::wrap(LevenbergMsolveCpp(AA, rrhhss, dd));
    return rcpp_result_gen;
END_RCPP
}
// LogAbsDetCpp
SEXP LogAbsDetCpp(SEXP AA);
RcppExport SEXP _spaMM_LogAbsDetCpp(SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(LogAbsDetCpp(AA));
    return rcpp_result_gen;
END_RCPP
}
// selfAdjointSolverCpp
SEXP selfAdjointSolverCpp(SEXP AA);
RcppExport SEXP _spaMM_selfAdjointSolverCpp(SEXP AASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    rcpp_result_gen = Rcpp::wrap(selfAdjointSolverCpp(AA));
    return rcpp_result_gen;
END_RCPP
}
// dgCprod
SEXP dgCprod(SEXP AA, SEXP BB);
RcppExport SEXP _spaMM_dgCprod(SEXP AASEXP, SEXP BBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    Rcpp::traits::input_parameter< SEXP >::type BB(BBSEXP);
    rcpp_result_gen = Rcpp::wrap(dgCprod(AA, BB));
    return rcpp_result_gen;
END_RCPP
}
// dgCcrossprod
SEXP dgCcrossprod(SEXP AA, SEXP BB, bool as_mat);
RcppExport SEXP _spaMM_dgCcrossprod(SEXP AASEXP, SEXP BBSEXP, SEXP as_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    Rcpp::traits::input_parameter< SEXP >::type BB(BBSEXP);
    Rcpp::traits::input_parameter< bool >::type as_mat(as_matSEXP);
    rcpp_result_gen = Rcpp::wrap(dgCcrossprod(AA, BB, as_mat));
    return rcpp_result_gen;
END_RCPP
}
// dgCtcrossprod
SEXP dgCtcrossprod(SEXP AA, SEXP BB);
RcppExport SEXP _spaMM_dgCtcrossprod(SEXP AASEXP, SEXP BBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    Rcpp::traits::input_parameter< SEXP >::type BB(BBSEXP);
    rcpp_result_gen = Rcpp::wrap(dgCtcrossprod(AA, BB));
    return rcpp_result_gen;
END_RCPP
}
// crossprod_not_dge
SEXP crossprod_not_dge(SEXP AA, SEXP BB, bool eval_dens, bool as_mat, bool keep_names);
RcppExport SEXP _spaMM_crossprod_not_dge(SEXP AASEXP, SEXP BBSEXP, SEXP eval_densSEXP, SEXP as_matSEXP, SEXP keep_namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    Rcpp::traits::input_parameter< SEXP >::type BB(BBSEXP);
    Rcpp::traits::input_parameter< bool >::type eval_dens(eval_densSEXP);
    Rcpp::traits::input_parameter< bool >::type as_mat(as_matSEXP);
    Rcpp::traits::input_parameter< bool >::type keep_names(keep_namesSEXP);
    rcpp_result_gen = Rcpp::wrap(crossprod_not_dge(AA, BB, eval_dens, as_mat, keep_names));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_crossprod
SEXP Rcpp_crossprod(SEXP AA, SEXP BB, bool eval_dens, bool as_mat, bool keep_names);
RcppExport SEXP _spaMM_Rcpp_crossprod(SEXP AASEXP, SEXP BBSEXP, SEXP eval_densSEXP, SEXP as_matSEXP, SEXP keep_namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    Rcpp::traits::input_parameter< SEXP >::type BB(BBSEXP);
    Rcpp::traits::input_parameter< bool >::type eval_dens(eval_densSEXP);
    Rcpp::traits::input_parameter< bool >::type as_mat(as_matSEXP);
    Rcpp::traits::input_parameter< bool >::type keep_names(keep_namesSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_crossprod(AA, BB, eval_dens, as_mat, keep_names));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_Csum
SEXP Rcpp_Csum(SEXP AA, SEXP BB);
RcppExport SEXP _spaMM_Rcpp_Csum(SEXP AASEXP, SEXP BBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AA(AASEXP);
    Rcpp::traits::input_parameter< SEXP >::type BB(BBSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_Csum(AA, BB));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_backsolve
SEXP Rcpp_backsolve(SEXP r, SEXP x, bool upper_tri, bool transpose);
RcppExport SEXP _spaMM_Rcpp_backsolve(SEXP rSEXP, SEXP xSEXP, SEXP upper_triSEXP, SEXP transposeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type r(rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type upper_tri(upper_triSEXP);
    Rcpp::traits::input_parameter< bool >::type transpose(transposeSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_backsolve(r, x, upper_tri, transpose));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_backsolve_M_M
SEXP Rcpp_backsolve_M_M(SEXP r, SEXP x, bool upper_tri, bool transpose);
RcppExport SEXP _spaMM_Rcpp_backsolve_M_M(SEXP rSEXP, SEXP xSEXP, SEXP upper_triSEXP, SEXP transposeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type r(rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type upper_tri(upper_triSEXP);
    Rcpp::traits::input_parameter< bool >::type transpose(transposeSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_backsolve_M_M(r, x, upper_tri, transpose));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_chol2solve
SEXP Rcpp_chol2solve(SEXP r, SEXP x);
RcppExport SEXP _spaMM_Rcpp_chol2solve(SEXP rSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type r(rSEXP);
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_chol2solve(r, x));
    return rcpp_result_gen;
END_RCPP
}
// adhoc_shermanMstep_sp
SEXP adhoc_shermanMstep_sp(SEXP AAinv, SEXP uu);
RcppExport SEXP _spaMM_adhoc_shermanMstep_sp(SEXP AAinvSEXP, SEXP uuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type AAinv(AAinvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type uu(uuSEXP);
    rcpp_result_gen = Rcpp::wrap(adhoc_shermanMstep_sp(AAinv, uu));
    return rcpp_result_gen;
END_RCPP
}
// adhoc_shermanM_sp
SEXP adhoc_shermanM_sp(SEXP QQt, SEXP iindic);
RcppExport SEXP _spaMM_adhoc_shermanM_sp(SEXP QQtSEXP, SEXP iindicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type QQt(QQtSEXP);
    Rcpp::traits::input_parameter< SEXP >::type iindic(iindicSEXP);
    rcpp_result_gen = Rcpp::wrap(adhoc_shermanM_sp(QQt, iindic));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spaMM_Rcpp_COMP_Z", (DL_FUNC) &_spaMM_Rcpp_COMP_Z, 4},
    {"_spaMM_COMP_Z_integrand", (DL_FUNC) &_spaMM_COMP_Z_integrand, 6},
    {"_spaMM_Rcpp_COMP_Z_asympto", (DL_FUNC) &_spaMM_Rcpp_COMP_Z_asympto, 2},
    {"_spaMM_lmwith_sparse_LDLp", (DL_FUNC) &_spaMM_lmwith_sparse_LDLp, 5},
    {"_spaMM_lmwith_sparse_LLp", (DL_FUNC) &_spaMM_lmwith_sparse_LLp, 5},
    {"_spaMM_update_R_in_place", (DL_FUNC) &_spaMM_update_R_in_place, 1},
    {"_spaMM_lmwithQR", (DL_FUNC) &_spaMM_lmwithQR, 4},
    {"_spaMM_lmwithQRP", (DL_FUNC) &_spaMM_lmwithQRP, 4},
    {"_spaMM_lmwith_sparse_QRp", (DL_FUNC) &_spaMM_lmwith_sparse_QRp, 5},
    {"_spaMM_Rcpp_chol_R", (DL_FUNC) &_spaMM_Rcpp_chol_R, 1},
    {"_spaMM_Rcpp_dense_cbind_mat_mat", (DL_FUNC) &_spaMM_Rcpp_dense_cbind_mat_mat, 2},
    {"_spaMM_Rcpp_dense_cbind_mat_vec", (DL_FUNC) &_spaMM_Rcpp_dense_cbind_mat_vec, 2},
    {"_spaMM_RcppMatrixCb2", (DL_FUNC) &_spaMM_RcppMatrixCb2, 2},
    {"_spaMM_lowertri", (DL_FUNC) &_spaMM_lowertri, 3},
    {"_spaMM_uppertri", (DL_FUNC) &_spaMM_uppertri, 3},
    {"_spaMM_rC_inv_chol_cpp", (DL_FUNC) &_spaMM_rC_inv_chol_cpp, 1},
    {"_spaMM_C_calc_cov_from_ranCoef", (DL_FUNC) &_spaMM_C_calc_cov_from_ranCoef, 2},
    {"_spaMM_makelong", (DL_FUNC) &_spaMM_makelong, 2},
    {"_spaMM_makelong2", (DL_FUNC) &_spaMM_makelong2, 2},
    {"_spaMM_C_dispInv", (DL_FUNC) &_spaMM_C_dispInv, 1},
    {"_spaMM_logit", (DL_FUNC) &_spaMM_logit, 1},
    {"_spaMM_is_evaluated", (DL_FUNC) &_spaMM_is_evaluated, 2},
    {"_spaMM_is_promise2", (DL_FUNC) &_spaMM_is_promise2, 2},
    {"_spaMM_nuln_plus_bessel_lnKnu", (DL_FUNC) &_spaMM_nuln_plus_bessel_lnKnu, 2},
    {"_spaMM_set_thread_nbr", (DL_FUNC) &_spaMM_set_thread_nbr, 1},
    {"_spaMM_rankinfo", (DL_FUNC) &_spaMM_rankinfo, 2},
    {"_spaMM_leverages", (DL_FUNC) &_spaMM_leverages, 1},
    {"_spaMM_sparse_cwiseprod", (DL_FUNC) &_spaMM_sparse_cwiseprod, 2},
    {"_spaMM_sweepZ1W", (DL_FUNC) &_spaMM_sweepZ1W, 2},
    {"_spaMM_ZWZt", (DL_FUNC) &_spaMM_ZWZt, 2},
    {"_spaMM_ZtWZ", (DL_FUNC) &_spaMM_ZtWZ, 2},
    {"_spaMM_Rcpp_d2hdv2", (DL_FUNC) &_spaMM_Rcpp_d2hdv2, 3},
    {"_spaMM_RcppChol", (DL_FUNC) &_spaMM_RcppChol, 1},
    {"_spaMM_crossprodCpp_d", (DL_FUNC) &_spaMM_crossprodCpp_d, 2},
    {"_spaMM_tcrossprodCpp", (DL_FUNC) &_spaMM_tcrossprodCpp, 2},
    {"_spaMM_LevenbergMsolveCpp", (DL_FUNC) &_spaMM_LevenbergMsolveCpp, 3},
    {"_spaMM_LogAbsDetCpp", (DL_FUNC) &_spaMM_LogAbsDetCpp, 1},
    {"_spaMM_selfAdjointSolverCpp", (DL_FUNC) &_spaMM_selfAdjointSolverCpp, 1},
    {"_spaMM_dgCprod", (DL_FUNC) &_spaMM_dgCprod, 2},
    {"_spaMM_dgCcrossprod", (DL_FUNC) &_spaMM_dgCcrossprod, 3},
    {"_spaMM_dgCtcrossprod", (DL_FUNC) &_spaMM_dgCtcrossprod, 2},
    {"_spaMM_crossprod_not_dge", (DL_FUNC) &_spaMM_crossprod_not_dge, 5},
    {"_spaMM_Rcpp_crossprod", (DL_FUNC) &_spaMM_Rcpp_crossprod, 5},
    {"_spaMM_Rcpp_Csum", (DL_FUNC) &_spaMM_Rcpp_Csum, 2},
    {"_spaMM_Rcpp_backsolve", (DL_FUNC) &_spaMM_Rcpp_backsolve, 4},
    {"_spaMM_Rcpp_backsolve_M_M", (DL_FUNC) &_spaMM_Rcpp_backsolve_M_M, 4},
    {"_spaMM_Rcpp_chol2solve", (DL_FUNC) &_spaMM_Rcpp_chol2solve, 2},
    {"_spaMM_adhoc_shermanMstep_sp", (DL_FUNC) &_spaMM_adhoc_shermanMstep_sp, 2},
    {"_spaMM_adhoc_shermanM_sp", (DL_FUNC) &_spaMM_adhoc_shermanM_sp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_spaMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
