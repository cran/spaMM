// cf spaMM/devel/RcppDevel.R to generate this
#include "spaMM_linear.h"

using namespace Rcpp ;

using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Upper;
using Eigen::VectorXd;

// [[Rcpp::export]]
SEXP lmwithQ( SEXP XX, SEXP yy ){
  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const Map<VectorXd> y(as<Map<VectorXd> >(yy));
//int n = X.rows();
//if (y.size() != n) throw invalid_argument("size mismatch");
  const Eigen::HouseholderQR<MatrixXd> QR(X);
  return List::create(Named("coef") = VectorXd(QR.solve(y)),
                    Named("Q") = MatrixXd(QR.householderQ()));
}

// [[Rcpp::export]]
SEXP leverages( SEXP XX ){
  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const int c(X.cols()); 
  const Eigen::HouseholderQR<MatrixXd> QR(X);
  MatrixXd Qthin(MatrixXd(QR.householderQ()).leftCols(c)); // householderQ SLOW: main bottleneck for GLMMs
  return wrap(VectorXd(Qthin.cwiseProduct(Qthin).rowwise().sum())); //returns vector of leverages rather than thin Q matrix
}

// [[Rcpp::export]]
SEXP Rcpp_qr_Q( SEXP XX){
  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const Eigen::HouseholderQR<MatrixXd> QR(X);
  return wrap(MatrixXd(QR.householderQ())); // returns thin Q matrix
}

// [[Rcpp::export]]
SEXP Rcpp_QR( SEXP XX){
  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const Eigen::HouseholderQR<MatrixXd> QR(X);
  List out = List::create(Named("Q") = MatrixXd(QR.householderQ()), // thin Q matrix
                      Named("R") = MatrixXd(QR.matrixQR().triangularView<Upper>())); 
  out.attr("class") = "Rcpp-QR";
  return(out);
}

// [[Rcpp::export]]
SEXP sweepZ1W( SEXP ZZ, SEXP WW ){
Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
const Map<VectorXd> W(as<Map<VectorXd> >(WW));
return wrap(W.asDiagonal() *Z);
}

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

// [[Rcpp::export]]
SEXP ZWZt( SEXP ZZ, SEXP WW ){
Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
const Map<VectorXd> W(as<Map<VectorXd> >(WW));
//MatrixXd resu=Z * W.asDiagonal() *Z.adjoint(); 
return wrap(Z * W.asDiagonal() *Z.adjoint());
}

// [[Rcpp::export]]
SEXP ZtWZ( SEXP ZZ, SEXP WW ){
Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
const Map<VectorXd> W(as<Map<VectorXd> >(WW));
//MatrixXd resu=Z.adjoint() * W.asDiagonal() *Z; 
return wrap(Z.adjoint() * W.asDiagonal() *Z);
}


// [[Rcpp::export]]
SEXP RcppChol( SEXP AA ){
  const Eigen::LLT<MatrixXd> llt(as<Map<MatrixXd> >(AA));
  int indic=0;
  switch(llt.info()) {
     case Eigen::Success   : indic=1 ; break;
     case Eigen::NumericalIssue : indic=2 ; break;
     case Eigen::NoConvergence : indic=3 ; break;
     case Eigen::InvalidInput : indic=4 ; break;
  }
  return List::create(Named("L") = MatrixXd(llt.matrixL()),
                      Named("Status")= indic); // Status==1 (==TRUE) means Success
} // such that A = L.Lt ie standard, NOT R's chol()

// https://gist.github.com/bobthecat/6509321
// [[Rcpp::export]]
SEXP crossprodCpp(SEXP Mat){
const Map<MatrixXd> A(as<Map<MatrixXd> >(Mat));
const int n(A.cols());
MatrixXd AtA(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.transpose()));
return wrap(AtA);
} // correspond bien a crossprod()

// [[Rcpp::export]]
SEXP tcrossprodCpp(SEXP Mat){
const Map<MatrixXd> A(as<Map<MatrixXd> >(Mat));
const int n(A.rows());
MatrixXd tAA(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A));
return wrap(tAA);
}

// [[Rcpp::export]] 
SEXP pseudoSolvediag( SEXP XX, SEXP bb ){ // pseudoSolvediag(X,b) should == solve.qr(qr(X),diag(b))
  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const Map<VectorXd> b(as<Map<VectorXd> >(bb));
  const int c(X.cols());
  const Eigen::HouseholderQR<MatrixXd> PQR(X);
  const MatrixXd pseudoInvX(PQR.matrixQR().topRows(c).triangularView<Upper>().solve(MatrixXd(PQR.householderQ()).leftCols(c).adjoint()));
return wrap(pseudoInvX * b.asDiagonal());
}

// [[Rcpp::export]]
SEXP LevenbergMsolveCpp( SEXP AA, SEXP rrhhss, SEXP dd ){ // 
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  const Map<VectorXd> rhs(as<Map<VectorXd> >(rrhhss));
  const double damping(as<double>(dd));
  const int n(A.cols());
  MatrixXd AtAdDpD(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint()));
  const VectorXd dampDpD(damping * VectorXd(AtAdDpD.diagonal()));
  AtAdDpD += dampDpD.asDiagonal();
  const Eigen::HouseholderQR<MatrixXd> PQR(AtAdDpD);
  VectorXd dbetaV(PQR.solve(rhs));
  return List::create(Named("dbetaV") = dbetaV,Named("dampDpD")=dampDpD);
}


// [[Rcpp::export]]
SEXP LogAbsDetCpp( SEXP AA ) {
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  Eigen::PartialPivLU<MatrixXd> LU(A);
  return wrap(LU.matrixLU().diagonal().array().abs().log().sum()); // valid because L is *unit*-lower-triangular ...
}

// contrary to svd(), the following function ensures that for symm matrices A=u.d.t(u) [with svd, v  can be -u and d can be the opposite of the eigenvalues]
// [[Rcpp::export]]
SEXP selfAdjointSolverCpp( SEXP AA ){ 
 const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
 Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
 es.compute(A); 
 return List::create(Named("u") = es.eigenvectors(),Named("d")=es.eigenvalues());
}

//FR->FR unfortunately useless:
// [[Rcpp::export]]
SEXP LevPerturbedQCpp( SEXP perturbedwAugX, SEXP pforREML, SEXP RpRu, SEXP RpRd, SEXP lamOverLam0, SEXP phiOverPhi0 ){
 const Map<MatrixXd> A(as<Map<MatrixXd> >(perturbedwAugX));
 const unsigned int p(as<unsigned int>(pforREML));
 const  Map<MatrixXd> RpRU(as< Map<MatrixXd> >(RpRu));
 const  Map<VectorXd> RpRdiag(as< Map<VectorXd> >(RpRd));
 const double a(sqrt(as<double>(lamOverLam0)));
 const double b(sqrt(as<double>(phiOverPhi0)));
 double psi = b/a; psi *=psi; psi -=1; // (b/a)^2 -1
 VectorXd W(RpRdiag.size());
 for (int it=0; it<W.size();it++) {W(it) = b*b/(RpRdiag(it)+psi);}
 MatrixXd BB = RpRU * W.asDiagonal() * RpRU.transpose();
 Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
 es.compute(BB.topLeftCorner(p,p));
 MatrixXd gauche = BB.leftCols(p) * es.eigenvectors();
 VectorXd diag(es.eigenvalues().size());
 for (int it=0; it<diag.size();it++) {diag(it) = 1/((es.eigenvalues())(it)-b*b/psi);}
 MatrixXd perturbedRRp = (BB - gauche * diag.asDiagonal() * gauche.adjoint());
 VectorXd leverages(A.rows());
 for (int it=0; it<leverages.size();it++) {leverages(it) = A.row(it) * perturbedRRp * A.row(it).transpose();}
 return(wrap(leverages));
}
