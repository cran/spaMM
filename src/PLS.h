#ifndef _PLS_H_
#define _PLS_H_

#include "spaMM_linear.h"
using namespace Rcpp;
using namespace std;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Upper;
using Eigen::VectorXd;

extern bool print_sparse_QR;

template<typename OrderingType>
SEXP sparse_LDL_from_XtX_oT( SEXP XX, bool pivot ){
  if (printDebug || print_sparse_QR)   Rcout <<"begin sparse_LDLp_from_XtX()"<<std::endl;
  const Eigen::MappedSparseMatrix<double> XtX(as<Eigen::MappedSparseMatrix<double> >(XX));
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double, Eigen::ColMajor> , Eigen::Lower, OrderingType > LDLp(XtX);
  List resu=List::create();
  resu["D_scaled"] = LDLp.vectorD();
  resu["XtX"] = XtX; // gloups
  if (pivot) {
    resu["sortPerm"] = LDLp.permutationP().indices();
    resu["perm"] = LDLp.permutationPinv().indices();
  } // more direct way to test the OrderingType ? 
  Eigen::SparseMatrix<double> U = LDLp.matrixU();  
  resu["U_scaled"] = U;
  if (printDebug || print_sparse_QR)   Rcout <<"end sparse_LDLp_from_XtX()"<<std::endl;
  return resu;
}

template<typename OrderingType>
SEXP sparse_LL_from_XtX_oT( SEXP XX, bool pivot ){
  if (printDebug || print_sparse_QR)   Rcout <<"begin sparse_LLp_from_XtX()"<<std::endl;
  const Eigen::MappedSparseMatrix<double> XtX(as<Eigen::MappedSparseMatrix<double> >(XX));
  Eigen::SimplicialLLT< Eigen::SparseMatrix<double, Eigen::ColMajor> , Eigen::Lower, OrderingType > LLp(XtX);
  List resu=List::create();
  resu["XtX"] = XtX; // gloups
  if (pivot) {
    resu["sortPerm"] = LLp.permutationP().indices();
    resu["perm"] = LLp.permutationPinv().indices();
  } // more direct way to test the OrderingType ? 
  Eigen::SparseMatrix<double> R = LLp.matrixU();  
  resu["R_scaled"] = R;
  if (printDebug || print_sparse_QR)   Rcout <<"end sparse_LLp_from_XtX()"<<std::endl;
  return resu;
}


template<typename OrderingType>
SEXP lmwith_sparse_LDL_oT( SEXP XX, SEXP yy, 
                           bool returntQ, // I G N O R E D but for consistent interface (cf get_from_default.Matrix)
                           bool returnR, bool pivot ){
  if (printDebug || print_sparse_QR)   Rcout <<"begin lmwith_sparse_LDL_oT()"<<std::endl;
  const Eigen::MappedSparseMatrix<double> X(as<Eigen::MappedSparseMatrix<double> >(XX));
  int nc=X.cols();
  Eigen::SparseMatrix<double> XtX(nc,nc); // resize necessary before rankUpdate
  XtX= X.transpose() * X;
  //https://forum.kde.org/viewtopic.php?f=74&t=95552&sid=a52d201181d55e50c3105c1213caddaf
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double, Eigen::ColMajor> , Eigen::Lower, OrderingType > LDLp(XtX);
  List resu=List::create();
  if (! Rf_isNull(yy)) {
    const Map<VectorXd> y(as<Map<VectorXd> >(yy));
    //VectorXd coef = LDLp.solve( (y.transpose() * X).transpose() );
    resu["coef"] = LDLp.solve( X.transpose() * y );
  }
  if (returnR) {
    resu["D_scaled"] = LDLp.vectorD();
    resu["XtX"] = XtX; // gloups
    if (pivot) {
      resu["sortPerm"] = LDLp.permutationP().indices();
      resu["perm"] = LDLp.permutationPinv().indices();
    } // more direct way to test the OrderingType ? 
    Eigen::SparseMatrix<double> U = LDLp.matrixU();  
    resu["U_scaled"] = U;
  }
  if (printDebug || print_sparse_QR)   Rcout <<"end lmwith_sparse_LDL_oT()"<<std::endl;
  return resu;
}

template<typename OrderingType>
SEXP lmwith_sparse_LL_oT( SEXP XX, SEXP yy, 
                           bool returntQ, // I G N O R E D but for consistent interface (cf get_from_default.Matrix)
                           bool returnR, bool pivot ){
  if (printDebug || print_sparse_QR)   Rcout <<"begin lmwith_sparse_LL_oT()"<<std::endl;
  const Eigen::MappedSparseMatrix<double> X(as<Eigen::MappedSparseMatrix<double> >(XX));
  int nc=X.cols();
  Eigen::SparseMatrix<double> XtX(nc,nc); // resize necessary before rankUpdate
  XtX= X.transpose() * X;
  //https://forum.kde.org/viewtopic.php?f=74&t=95552&sid=a52d201181d55e50c3105c1213caddaf
  Eigen::SimplicialLLT< Eigen::SparseMatrix<double, Eigen::ColMajor> , Eigen::Lower, OrderingType > LLp(XtX);
  List resu=List::create();
  if (! Rf_isNull(yy)) {
    const Map<VectorXd> y(as<Map<VectorXd> >(yy));
    //VectorXd coef = LDLp.solve( (y.transpose() * X).transpose() );
    resu["coef"] = LLp.solve( X.transpose() * y );
  }
  if (returnR) {
    resu["XtX"] = XtX; // gloups
    if (pivot) {
      resu["sortPerm"] = LLp.permutationP().indices();
      resu["perm"] = LLp.permutationPinv().indices();
    } // more direct way to test the OrderingType ? 
    Eigen::SparseMatrix<double> R = LLp.matrixU();  
    resu["R_scaled"] = R;
  }
  if (printDebug || print_sparse_QR)   Rcout <<"end lmwith_sparse_LL_oT()"<<std::endl;
  return resu;
}


template<typename OrderingType>
SEXP lmwith_sparse_QR_oT( SEXP XX, SEXP yy, 
                        bool returntQ, // I N H I B I T E D
                        bool returnR, bool pivot ){
  if (printDebug || print_sparse_QR)   Rcout <<"debut lmwith_sparse_QRp()"<<std::endl;
  const Eigen::MappedSparseMatrix<double> X(as<Eigen::MappedSparseMatrix<double> >(XX));
  Eigen::SparseQR< Eigen::SparseMatrix<double, Eigen::ColMajor> ,  OrderingType > QRp(X);
  List resu=List::create();
  if (! Rf_isNull(yy)) {
    const Map<VectorXd> y(as<Map<VectorXd> >(yy));
    resu["coef"] = VectorXd(QRp.solve(y));
  }
  if (false && returntQ) { // F A L S E: inhibit Q computation
    if(false) {
      /*
      * cf https://forum.kde.org/viewtopic.php?f=74&t=117500&p=293605#p292860:
      *    http://eigen.tuxfamily.org/bz/show_bug.cgi?id=596
      * This compounds bug(s) that makes returning Q difficult, 
      * and a very slow computation of Q!
      */
      if (print_sparse_QR)   Rcout <<"slow...";
      Eigen::SparseMatrix<double> Q;
      Q = QRp.matrixQ(); // clumsy (other way?)
      if (print_sparse_QR)   Rcout <<" ...done ";
      // SparseMatrix( transpose()) in that order to return a dg*C*Matrix
      resu["tQ"] = Eigen::SparseMatrix<double> (Q.leftCols(X.cols()).transpose());
    } else {
      // naive, but much faster! But still quite slow relative to a dense solver!
      Eigen::SparseMatrix<double> I(X.rows(),X.rows());
      I.setIdentity();
      resu["tQ"] = QRp.matrixR().topLeftCorner(X.cols(),X.cols()) * QRp.solve(I);
    }
  }
  if (returnR) {
    if (pivot) resu["P"] = QRp.colsPermutation().indices(); // more direct way to test the OrderingType ? 
    const int r(X.cols()); // r(QRP.rank());
    Eigen::SparseMatrix<double> R(r,r);
    R = QRp.matrixR().topLeftCorner(r,r);  
    resu["R"] = R;
  }
  if (printDebug || print_sparse_QR)   Rcout <<"fin lmwith_sparse_QRp()"<<std::endl;
  return resu;
}

#endif
