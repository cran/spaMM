#include "spaMM_linear.h"
#include "PLS.h"
//#include <Eigen/Jacobi> 

using namespace Rcpp;
using namespace std;

//using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Upper;
using Eigen::VectorXd;

bool print_sparse_QR=false;

// [[Rcpp::export(.lmwith_sparse_LDLp)]]
SEXP lmwith_sparse_LDLp( SEXP XX, SEXP yy, bool returntQ, bool returnR, bool pivot ){
  if (pivot) {
    return(lmwith_sparse_LDL_oT<Eigen::COLAMDOrdering<int> >( XX,  yy, returntQ,
                                                              returnR, pivot ));
  } else return(lmwith_sparse_LDL_oT<Eigen::NaturalOrdering<int> >( XX,  yy, returntQ,
                                                                    returnR, pivot )); 
}

// [[Rcpp::export(.lmwith_sparse_LLp)]]
SEXP lmwith_sparse_LLp( SEXP XX, SEXP yy, bool returntQ, bool returnR, bool pivot ){
  if (pivot) {
    return(lmwith_sparse_LL_oT<Eigen::COLAMDOrdering<int> >( XX,  yy, returntQ,
                                                              returnR, pivot ));
  } else return(lmwith_sparse_LL_oT<Eigen::NaturalOrdering<int> >( XX,  yy, returntQ,
                                                                    returnR, pivot )); 
}

int givens(std::vector<double>& resu, double& a, double& b) { // GolubL algo 5.1.3
  double tau, f;
  if (b==0.0) { 
    resu[0]=1.0;
    resu[1]=0.0;
  } else if (fabs(b)> fabs(a)) {
    tau=-a/b;
    f=1.0/sqrt(1.0+tau*tau);
    resu[0]=f*tau;
    resu[1]=f;
  } else {
    tau=-b/a;
    f=1.0/sqrt(1.0+tau*tau);
    resu[0]=f;
    resu[1]=f*tau;
  }
  return(0);
}

// [[Rcpp::export(.update_R_in_place)]]
SEXP update_R_in_place(SEXP RD) {
  Map<MatrixXd> A(as<Map<MatrixXd> >(RD));
  std::vector<double> ROT;
  ROT.resize(2);
  double x1,x2;
  Eigen::Index k, n_col=A.cols();
  for (Eigen::Index dt=0; dt<n_col; dt++) { 
/* perform the n(n-1)/2 Givens rotations suggested by Mor\'e 77 and slighlty mor detailed in NocedalW p.265
 #iteration dt=0 sets zeros on the diagonal of the diagonal block, but introduces non-zero 
 # on the upper-right triangle of what was the diagonal block
 #iteration dt=1 then fill the first super diagonal with zeros, still altering the triangle above it
 #iteration dt=2 then fill the next super diagonal with zeros, still altering the triangle above it
 #... iterate over super diagonals until the upper triangle is done. */
    for (Eigen::Index it=n_col-1; it>=dt; it--) {
      k=n_col-dt+it;
      givens(ROT, A(it,it), A(k,it));
      for (Eigen::Index jj=it; jj<n_col; jj++) {
        x1=A(it,jj);
        x2=A(k,jj);
        A(it,jj)=ROT[0]*x1-ROT[1]*x2;
        A(k,jj)=ROT[1]*x1+ROT[0]*x2;
      }
    }
  }
  return(wrap(A.topRows(n_col)));
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
    MatrixXd thinQ(X.rows(),X.cols());
    thinQ.setIdentity();
    thinQ= QR.householderQ() * thinQ; // v3.6.8 markedly faster than older code, see https://forum.kde.org/viewtopic.php?f=74&t=106635 ... but still comparatively slow.
    resu["t_Q_scaled"] = MatrixXd(thinQ.transpose()); 
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
    MatrixXd thinQ(X.rows(),X.cols());
    thinQ.setIdentity();
    thinQ= QRP.householderQ() * thinQ; 
    resu["t_Q_scaled"] = thinQ.transpose(); 
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
 * e.g. COLAMDOrdering. However, see Notes: leverages(X[,B]) are not deducible from Q[,B] 
 * => problem for Pdiag calculation in particular.
*/

// [[Rcpp::export(.lmwith_sparse_QRp)]]
SEXP lmwith_sparse_QRp( SEXP XX, SEXP yy, 
                       bool returntQ, // I N H I B I T E D
                       bool returnR, bool COLAMDO=true ){
  // R is pivoted in both cases !
  if (COLAMDO) {
    return(lmwith_sparse_QR_oT<Eigen::COLAMDOrdering<int> >( XX,  yy, 
                               returntQ, // I N H I B I T E D
                              returnR ));
  } else return(lmwith_sparse_QR_oT<Eigen::NaturalOrdering<int> >( XX,  yy, 
                                                                 returntQ, // I N H I B I T E D
                                                                 returnR )); 
}

// there are also .wrap_Ltri_t_chol() calling .RcppChol() which return the transpose.
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

// DEVEL code
// http://stackoverflow.com/questions/25147577/get-information-about-a-promise-without-evaluating-the-promise
// ## we can get info about the promise by using pryr::uneval(resid.model)
// ## class(uneval(resid.model)) may be "bytecode" if the corrHLfit argument was already an evaluated promise
// ## otherwise class(uneval(resid.model)) may be a "call".
// ## if needed, uneval is minimal R code  + the C++ code in promise.cpp: promise_code(resid_model) might suffice.
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

//https://stackoverflow.com/questions/31913437/r-fast-cbind-matrix-using-rcpp
// [[Rcpp::export(.Rcpp_dense_cbind_mat_mat)]]
NumericMatrix Rcpp_dense_cbind_mat_mat(NumericMatrix a, NumericMatrix b) {
  int acoln = a.ncol();
  int bcoln = b.ncol();
  NumericMatrix out = no_init_matrix(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      out(_, j) = b(_, j - acoln);
    }
  }
  return out;
}

// [[Rcpp::export(.Rcpp_dense_cbind_mat_vec)]]
NumericMatrix Rcpp_dense_cbind_mat_vec(NumericMatrix a, NumericVector b) {
  int acoln = a.ncol();
  NumericMatrix out = no_init_matrix(a.nrow(), acoln + 1);
  for (int j = 0; j < acoln ; j++)  out(_, j) = a(_, j);
  out(_, acoln) = b;
  return out;
}
