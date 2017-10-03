#include "spaMM_linear.h"
#include "PLS.h"

using namespace Rcpp;
using namespace std;

//using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Upper;
using Eigen::VectorXd;

bool print_sparse_QR=false;

// [[Rcpp::export(.sparse_LDLp_from_XtX)]]
SEXP sparse_LDLp_from_XtX( SEXP XX, bool pivot ){
  if (pivot) {
    return(sparse_LDL_from_XtX_oT<Eigen::COLAMDOrdering<int> >( XX, pivot ));
  } else return(sparse_LDL_from_XtX_oT<Eigen::NaturalOrdering<int> >( XX, pivot )); 
}

// [[Rcpp::export(.lmwith_sparse_LDLp)]]
SEXP lmwith_sparse_LDLp( SEXP XX, SEXP yy, bool returntQ, bool returnR, bool pivot ){
  if (pivot) {
    return(lmwith_sparse_LDL_oT<Eigen::COLAMDOrdering<int> >( XX,  yy, returntQ,
                                                              returnR, pivot ));
  } else return(lmwith_sparse_LDL_oT<Eigen::NaturalOrdering<int> >( XX,  yy, returntQ,
                                                                    returnR, pivot )); 
}

// [[Rcpp::export(.sparse_LLp_from_XtX)]]
SEXP sparse_LLp_from_XtX( SEXP XX, bool pivot ){
  if (pivot) {
    return(sparse_LL_from_XtX_oT<Eigen::COLAMDOrdering<int> >( XX, pivot ));
  } else return(sparse_LL_from_XtX_oT<Eigen::NaturalOrdering<int> >( XX, pivot )); 
}

// [[Rcpp::export(.lmwith_sparse_LLp)]]
SEXP lmwith_sparse_LLp( SEXP XX, SEXP yy, bool returntQ, bool returnR, bool pivot ){
  if (pivot) {
    return(lmwith_sparse_LL_oT<Eigen::COLAMDOrdering<int> >( XX,  yy, returntQ,
                                                              returnR, pivot ));
  } else return(lmwith_sparse_LL_oT<Eigen::NaturalOrdering<int> >( XX,  yy, returntQ,
                                                                    returnR, pivot )); 
}


// [[Rcpp::export(.lmwithQR)]]
SEXP lmwithQR( SEXP XX, SEXP yy, bool returntQ, bool returnR ){
  if (printDebug)   Rcout <<"debut lmwithQR()"<<std::endl;
  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const Eigen::HouseholderQR<MatrixXd> QR(X);
  List resu=List::create();
  if (! Rf_isNull(yy)) {
    const Map<VectorXd> y(as<Map<VectorXd> >(yy));
    resu["coef"] = VectorXd(QR.solve(y));
  }
  if (returntQ) {
    MatrixXd Q(QR.householderQ());
    resu["t_Q_scaled"] = MatrixXd(Q.leftCols(X.cols()).transpose()); 
  }
  if (returnR) {
    const int r(X.cols()); // r(QRP.rank());
    resu["R_scaled"] = MatrixXd(QR.matrixQR().topLeftCorner(r,r).triangularView<Upper>()); 
  }
  if (printDebug)   Rcout <<"fin lmwithQR()"<<std::endl;
  return resu;
}

// [[Rcpp::export(.lmwithQRP)]]
SEXP lmwithQRP( SEXP XX, SEXP yy, bool returntQ, bool returnR ){
  if (printDebug)   Rcout <<"debut lmwithQRP()"<<std::endl;
  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const Eigen::ColPivHouseholderQR<MatrixXd> QRP(X);
  List resu=List::create();
  if (! Rf_isNull(yy)) {
    const Map<VectorXd> y(as<Map<VectorXd> >(yy));
    resu["coef"] = VectorXd(QRP.solve(y));
  }
  if (returntQ) {
    MatrixXd Q(QRP.householderQ());
    resu["t_Q_scaled"] = Q.leftCols(X.cols()).transpose(); 
  }
  if (returnR) {
    resu["perm"] = QRP.colsPermutation().indices();
    const int r(X.cols()); // r(QRP.rank());
    // oddly, topLeftCorner(r,r) seems required for QRP but not for QP 
    resu["R_scaled"] = MatrixXd(QRP.matrixQR().topLeftCorner(r,r).triangularView<Upper>()); 
  }
  if (printDebug)   Rcout <<"fin lmwithQRP()"<<std::endl;
  return resu;
} 

/* 
 * SEXP lmwith_sparse_QRp( SEXP XX, SEXP yy, bool returntQ, bool returnR )
 * could be defined similarly to lmwith_sparse_QRp, with another ordering than Eigen::NaturalOrdering<int>,
 * e.g. COLAMDOrdering. However, see Notes: leverages(X[,B]) are not deducible from Q[,B] 
 * => problem for Pdiag calculation in particular.
*/

// [[Rcpp::export(.lmwith_sparse_QRp)]]
SEXP lmwith_sparse_QRp( SEXP XX, SEXP yy, 
                       bool returntQ, // I N H I B I T E D
                       bool returnR, bool pivot ){
  if (pivot) {
    return(lmwith_sparse_QR_oT<Eigen::COLAMDOrdering<int> >( XX,  yy, 
                               returntQ, // I N H I B I T E D
                              returnR, pivot ));
  } else return(lmwith_sparse_QR_oT<Eigen::NaturalOrdering<int> >( XX,  yy, 
                                                                 returntQ, // I N H I B I T E D
                                                                 returnR, pivot )); 
}

// [[Rcpp::export(.Rcpp_chol_R)]]
SEXP Rcpp_chol_R( SEXP AA ){ // upper tri as in R chol()
  if (printDebug)   Rcout <<"debut Rcpp_chol_R()"<<std::endl;
  const Eigen::LLT<MatrixXd> llt(as<Map<MatrixXd> >(AA));
  int indic=0;
  switch(llt.info()) {
  case Eigen::Success   : indic=1 ; break;
  case Eigen::NumericalIssue : indic=2 ; break;
  case Eigen::NoConvergence : indic=3 ; break;
  case Eigen::InvalidInput : indic=4 ; break;
  }
  List out = List::create(Named("R") = MatrixXd(llt.matrixU()),
                          Named("Status")= indic); // Status==1 means Success. Note that higher values are not success but as booleans, are interpreted as TRUE 
  out.attr("class") = "Rcpp_chol_R";
  if (printDebug)   Rcout <<"fin Rcpp_chol_R()"<<std::endl;
  return out;
} // such that A = R'.R as in R's chol()

// only in if (FALSE) {...} blocks of code: 
// [[Rcpp::export(.LevMar_cpp)]] 
SEXP LevMar_cpp( SEXP AA, SEXP LMrhs ){ // 
  if (printDebug)   Rcout <<"begin levMar_cpp"<<std::endl;
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  const Map<VectorXd> LM_rhs(as<Map<VectorXd> >(LMrhs));
  // usng generic QR
  // more efficient algo using givens rotations: https://eigen.tuxfamily.org/dox/unsupported/LMqrsolv_8h_source.html
  // cf aussi ancienne approche util. LDL sur AtAdDpD
  const Eigen::HouseholderQR<MatrixXd> QR_A(A);
  const int r(A.cols()); // r(QRP.rank());
  MatrixXd R = MatrixXd(QR_A.matrixQR().topLeftCorner(r,r));
  List resu=List::create();
  // first compute t(t(rhs).inv(R)) = inv(Rt).rhs
  VectorXd dVscaled_beta = R.triangularView<Upper>().solve<Eigen::OnTheRight>(LM_rhs.transpose()).transpose();
  // then inv(R) . inv(Rt).rhs  = inv(RtR).rhs
  resu["dVscaled_beta"] = R.triangularView<Upper>().solve( dVscaled_beta );
  resu["rhs"] = LMrhs;
  if (printDebug)   Rcout <<"end levMar_cpp"<<std::endl;
  return(resu);
}

// DEVEL code
// http://stackoverflow.com/questions/25147577/get-information-about-a-promise-without-evaluating-the-promise
/*
  // [[Rcpp::export]]
bool is_promise2(Symbol name, Environment env) {
  SEXP object = Rf_findVar(name, env);
  return (TYPEOF (object) == PROMSXP);
}

// [[Rcpp::export]]
SEXP promise_code(Symbol name, Environment env) {
  SEXP object = Rf_findVar(name, env);
  return PRCODE(object);
}
// [[Rcpp::export]]
SEXP promise_value(Symbol name, Environment env) {
  SEXP object = Rf_findVar(name, env);
  return PRVALUE(object);
}
*/
