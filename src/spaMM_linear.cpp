// cf spaMM/devel/RcppDevel.R to generate this
#include "spaMM_linear.h"

using namespace Rcpp;
using namespace std;

using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

bool printDebug=false;

// [[Rcpp::export(.rankinfo)]] // called only if spaMM.getOption("rankMethod") is NULL (non-default)
SEXP rankinfo( SEXP x, SEXP tol){
  const Map<MatrixXd> X(as<Map<MatrixXd> >(x));
  const double threshold(as<double>(tol));
  Eigen::ColPivHouseholderQR<MatrixXd> QR(X);
  QR.setThreshold(threshold);
  return List::create(Named("pivot")=QR.colsPermutation().indices(), 
                      // following https://stackoverflow.com/questions/26385561/how-to-get-matrix-pivot-index-in-c-eigen,
                      // but this does not seem to be equivalent to qr()$pivot
                      Named("rank")=int(QR.rank()),
                      Named("method")=".rankinfo");
}

// [[Rcpp::export(.leverages)]]
SEXP leverages( SEXP XX ){
  if (printDebug)   Rcout <<"debut leverages()"<<std::endl;
  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const int c(X.cols()); 
  const Eigen::HouseholderQR<MatrixXd> QR(X);
  MatrixXd Qthin(MatrixXd(QR.householderQ()).leftCols(c)); /* householderQ SLOW and memory-inefficient (full Q)
  The main weakness of Eigen's QR is that Qthin computation requires Qfull computation*/
  // transpose -> 1-row matrix OK for return as VectorXd -> R vector 
  if (printDebug)   Rcout <<"fin leverages()"<<std::endl;
  return wrap(VectorXd(Qthin.cwiseProduct(Qthin).rowwise().sum().transpose())); //returns vector of leverages rather than thin Q matrix
}

// [[Rcpp::export(.Rcpp_sweepZ1W)]]
SEXP sweepZ1W( SEXP ZZ, SEXP WW ){
  if (printDebug)   Rcout <<"debut sweepZ1W()"<<std::endl;
  const Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
  const Map<VectorXd> W(as<Map<VectorXd> >(WW));
  if (printDebug)   Rcout <<"fin sweepZ1W()"<<std::endl; // almost...
  return wrap(W.asDiagonal() *Z);
}

/* attempt to add the colnames here by the following code, which was perhaps correct 
 but failed (as does the simpler version) on my first test because the matrix was integer... 
 Poor error message, waste of time!
 
 
 if (printDebug)   Rcout <<"debut sweepZ1W()"<<std::endl;
 const Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
 const Map<VectorXd> W(as<Map<VectorXd> >(WW));
 // Consistent with R's diag(<vector>) %*% X that keeps only the X colnames:
 SEXP dimnames = Rf_getAttrib(ZZ, R_DimNamesSymbol);
 SEXP ZZ_colnames;
 bool hascolnames =  (! Rf_isNull(dimnames));
 if (hascolnames) {
 ZZ_colnames = VECTOR_ELT(dimnames, 1);
 hascolnames = ( ! Rf_isNull(dimnames));
 }
 if (hascolnames) {
 NumericMatrix Xout(wrap(W.asDiagonal() *Z)); 
 colnames(Xout) = ZZ_colnames;
 if (printDebug)   Rcout <<"fin sweepZ1W()"<<std::endl;
 return(Xout);
 } else {
 if (printDebug)   Rcout <<"fin sweepZ1W()"<<std::endl; // hmmm
 return(wrap(W.asDiagonal() *Z));
 }
*/

/*
variations on 
ZWZt <- function(Z,W) {base::tcrossprod(Z %*% diag(W),Z)}  ## naive
ZtWZ <- function(Z,W) {base::crossprod(Z*W,Z)}  ## naive
or
ZWZt <- function(Z,W) {base::tcrossprod(sweep(Z,MARGIN=2,W,`*`),Z)}  ## still SLOW
ZtWZ <- function(Z,W) {base::crossprod(Z,sweep(Z,MARGIN=1,W,`*`))}  ## still SLOW
It's interesting to replicate the benchmark at http://stackoverflow.com/questions/6886416/automatic-multiplication-between-vector-and-matrix
The results may be quite different from those shown there.
*/

// [[Rcpp::export(.ZWZt)]]
SEXP ZWZt( SEXP ZZ, SEXP WW ){
  if (printDebug)   Rcout <<"debut ZWZt()"<<std::endl;
  const Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
  const Map<VectorXd> W(as<Map<VectorXd> >(WW));
  //return wrap(Z * W.asDiagonal() *Z.adjoint());
  VectorXd sqrtW = W.array().sqrt();
  MatrixXd swZ = Z * sqrtW.asDiagonal();
  const int r(Z.rows());
  swZ = MatrixXd(r, r).setZero().selfadjointView<Lower>().rankUpdate(swZ);// as in tcrossprod
  if (printDebug)   Rcout <<"fin ZWZt()"<<std::endl;
  return wrap(swZ); 
}

// [[Rcpp::export(.ZtWZ)]]
SEXP ZtWZ( SEXP ZZ, SEXP WW ){
  if (printDebug)   Rcout <<"debut ZtWZ()"<<std::endl;
  const Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
  const int c(Z.cols());
  if (c==0) return(wrap(MatrixXd(0,0))); // which the alternative code would do, but with an UBSAN issue.
  const Map<VectorXd> W(as<Map<VectorXd> >(WW));
  VectorXd sqrtW = W.array().sqrt();
  MatrixXd swZ = sqrtW.asDiagonal() *Z;
  // *this.rankUpdate(u) performs this= this+ (alpha=1) uu*
  swZ = MatrixXd(c, c).setZero().selfadjointView<Lower>().rankUpdate(swZ.transpose());// as in crossprod
  if (printDebug)   Rcout <<"fin ZtWZ()"<<std::endl;
  return wrap(swZ); 
}

// [[Rcpp::export(.Rcpp_d2hdv2)]]
SEXP Rcpp_d2hdv2( SEXP ZZ, SEXP WA, SEXP WB ){
  if (printDebug)   Rcout <<"debut Rcpp_d2hdv2()"<<std::endl;
  const Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
  if (Z.cols()==0) return(wrap(MatrixXd(0,0))); // which the alternative code would do, but with an UBSAN issue.
  const Map<VectorXd> wa(as<Map<VectorXd> >(WA));
  const Map<VectorXd> wb(as<Map<VectorXd> >(WB));
  VectorXd sqrtW = wa.array().sqrt();
  MatrixXd swZ = sqrtW.asDiagonal() *Z;
  const int c(swZ.cols());
  // *this.rankUpdate(u) performs this= this+ (alpha=1) uu*
  swZ = MatrixXd(c, c).setZero().selfadjointView<Lower>().rankUpdate(swZ.transpose());// as in crossprod
  swZ.diagonal() += wb;
  if (printDebug)   Rcout <<"fin Rcpp_d2hdv2()"<<std::endl;
  return wrap( - swZ); 
}


// [[Rcpp::export(.RcppChol)]]
SEXP RcppChol( SEXP AA ){ // returns t(base::chol)
  if (printDebug)   Rcout <<"debut RcppChol()"<<std::endl;
  const Eigen::LLT<MatrixXd> llt(as<Map<MatrixXd> >(AA));
  int indic=0;
  switch(llt.info()) {
  case Eigen::Success   : indic=1 ; break;
  case Eigen::NumericalIssue : indic=2 ; break;
  case Eigen::NoConvergence : indic=3 ; break;
  case Eigen::InvalidInput : indic=4 ; break;
  }
  List out = List::create(Named("L") = MatrixXd(llt.matrixL()),
                          Named("Status")= indic); // Status==1 means Success. Note that higher values are not success but as booleans, are interpreted as TRUE 
  out.attr("class") = "RcppChol";
  if (printDebug)   Rcout <<"fin RcppChol()"<<std::endl;
  return out;
} // such that A = L.Lt ie standard, NOT R's chol()

// https://gist.github.com/bobthecat/6509321 // wrapped bu .crossprod. Remain clearly faster for large matrices 09/2018 but not 3*3 matrix...
// [[Rcpp::export(.crossprodCpp)]]
SEXP crossprodCpp(SEXP Mat, SEXP yy){
  if (printDebug)   Rcout <<"debut crossprodCpp()"<<std::endl;
  const Map<MatrixXd> A(as<Map<MatrixXd> >(Mat));
  MatrixXd tAA;
  if (Rf_isNull(yy)) {
    const int c(A.cols());
    tAA= MatrixXd(c, c).setZero().selfadjointView<Lower>().rankUpdate(A.transpose());  
  } else {
    const Map<MatrixXd> y(as<Map<MatrixXd> >(yy));
    tAA = A.transpose() * y;
  }  
  if (printDebug)   Rcout <<"fin crossprodCpp()"<<std::endl;
  return wrap(tAA);
} // correspond bien a crossprod()

// [[Rcpp::export(.tcrossprodCpp)]]
SEXP tcrossprodCpp(SEXP Mat, SEXP yy){
  if (printDebug)   Rcout <<"debut tcrossprodCpp()"<<std::endl;
  const Map<MatrixXd> A(as<Map<MatrixXd> >(Mat));
  MatrixXd AAt;
  if (Rf_isNull(yy)) {
    const int r(A.rows());
    AAt = MatrixXd(r, r).setZero().selfadjointView<Lower>().rankUpdate(A);
  } else {
    const Map<MatrixXd> y(as<Map<MatrixXd> >(yy));
    AAt = A * y.transpose();
  }  
  if (printDebug)   Rcout <<"fin tcrossprodCpp()"<<std::endl;
  return wrap(AAt);
}

// [[Rcpp::export(.LevenbergMsolveCpp)]]
SEXP LevenbergMsolveCpp( SEXP AA, SEXP rrhhss, SEXP dd ){ // 
  if (printDebug)   Rcout <<"debut LevenbergMsolveCpp()"<<std::endl;
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  const Map<VectorXd> rhs(as<Map<VectorXd> >(rrhhss));
  const double damping(as<double>(dd));
  const int nc(A.cols());
  MatrixXd AtAdDpD(MatrixXd(nc, nc).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint()));
  const VectorXd dampDpD(damping * VectorXd(AtAdDpD.diagonal()));
  AtAdDpD += dampDpD.asDiagonal();
  const Eigen::LDLT<MatrixXd> PQR(AtAdDpD); // (LDLt is a version of cholesky) faster than QR
  VectorXd dbetaV = PQR.solve(rhs);
  if (printDebug)   Rcout <<"fin LevenbergMsolveCpp()"<<std::endl;
  return List::create(Named("dbetaV") = dbetaV,Named("dampDpD")=dampDpD);
}


// [[Rcpp::export(.LogAbsDetCpp)]]
SEXP LogAbsDetCpp( SEXP AA ) {
  if (printDebug)   Rcout <<"debut LogAbsDetCpp()"<<std::endl;
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  Eigen::PartialPivLU<MatrixXd> LU(A);
  if (printDebug)   Rcout <<"fin LogAbsDetCpp()"<<std::endl;
  return wrap(LU.matrixLU().diagonal().array().abs().log().sum()); // valid because L is *unit*-lower-triangular ...
}

// contrary to svd(), the following function ensures that for symm matrices A=u.d.t(u) [with svd, v  can be -u and d can be the opposite of the eigenvalues]
// [[Rcpp::export(.selfAdjointSolverCpp)]]
SEXP selfAdjointSolverCpp( SEXP AA ){ 
  if (printDebug)   Rcout <<"debut selfAdjointSolverCpp()"<<std::endl;
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(A); 
  if (printDebug)   Rcout <<"fin selfAdjointSolverCpp()"<<std::endl;
  return List::create(Named("vectors") = es.eigenvectors(),Named("values")=es.eigenvalues());
}
