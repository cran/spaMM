// cf spaMM/devel/RcppDevel.R to generate this
#include "spaMM_linear.h"

using namespace Rcpp;
using namespace std;

using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Upper;
using Eigen::VectorXd;
// for sparse QR:
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::COLAMDOrdering;

bool printDebug=false;

// [[Rcpp::export(.Rcpp_as_dgCMatrix)]]
SEXP Rcpp_as_dgCMatrix_( SEXP XX_ ){ 
  // corrected from http://gallery.rcpp.org/articles/sparse-matrix-coercion/ where a NULL "dimnames" attribute is not correctly handled.
  typedef Eigen::SparseMatrix<double> SpMat;   
  typedef Eigen::Map<Eigen::MatrixXd> MapMatd; // Input: must be double
  MapMatd X(Rcpp::as<MapMatd>(XX_));
  SpMat Xsparse = X.sparseView();              // Output: sparse matrix
  S4 Xout(wrap(Xsparse));                      // Output: as S4 object
  SEXP dimnames = Rf_getAttrib(XX_, R_DimNamesSymbol);
  if (Rf_isNull(dimnames)) { 
    Xout.slot("Dimnames") = List::create(R_NilValue,R_NilValue);
  } else {
    NumericMatrix Xin(XX_);                      // Copy dimnames
    Xout.slot("Dimnames") = clone(List(Xin.attr("dimnames")));
  }
  return(Xout);
}

// [[Rcpp::export(.rankinfo)]]
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
  MatrixXd Qthin(MatrixXd(QR.householderQ()).leftCols(c)); // householderQ SLOW: main bottleneck for GLMMs
  // transpose -> 1-row matrix OK for return as VectorXd -> R vector 
  if (printDebug)   Rcout <<"fin leverages()"<<std::endl;
  return wrap(VectorXd(Qthin.cwiseProduct(Qthin).rowwise().sum().transpose())); //returns vector of leverages rather than thin Q matrix
}

// [[Rcpp::export(.Rcpp_sweepZ1W)]]
SEXP sweepZ1W( SEXP ZZ, SEXP WW ){
  if (printDebug)   Rcout <<"debut sweepZ1W()"<<std::endl;
  const Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
  const Map<VectorXd> W(as<Map<VectorXd> >(WW));
  if (printDebug)   Rcout <<"fin sweepZ1W()"<<std::endl;
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

// (fixme?) This function does not seem currently used
// [[Rcpp::export(.Rcpp_Sig)]]
SEXP Rcpp_Sig( SEXP ZZ, SEXP WA, SEXP WB ){
  if (printDebug)   Rcout <<"debut Rcpp_Sig()"<<std::endl;
  const Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
  const Map<VectorXd> wa(as<Map<VectorXd> >(WA));
  const Map<VectorXd> wb(as<Map<VectorXd> >(WB));
  VectorXd sqrtW = wa.array().sqrt();
  MatrixXd swZ = Z * sqrtW.asDiagonal();
  const int r(swZ.rows());
  swZ = MatrixXd(r, r).setZero().selfadjointView<Lower>().rankUpdate(swZ);
  swZ.diagonal() += wb;
  if (printDebug)   Rcout <<"fin Rcpp_Sig()"<<std::endl;
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
SEXP RcppChol( SEXP AA ){ // returns t(R::chol)
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

// https://gist.github.com/bobthecat/6509321
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

// // [[Rcpp::export]] 
// SEXP pseudoSolvediag( SEXP XX, SEXP bb ){ // pseudoSolvediag(X,b) should == solve.qr(qr(X),diag(b))
//   if (printDebug)   Rcout <<"debut pseudoSolvediag()"<<std::endl;
//   const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
//   const Map<VectorXd> b(as<Map<VectorXd> >(bb));
//   const int c(X.cols());
//   const Eigen::HouseholderQR<MatrixXd> PQR(X);
//   const MatrixXd pseudoInvX(PQR.matrixQR().topRows(c).triangularView<Upper>().solve(MatrixXd(PQR.householderQ()).leftCols(c).adjoint()));
//   if (printDebug)   Rcout <<"fin pseudoSolvediag()"<<std::endl;
//   return wrap(pseudoInvX * b.asDiagonal());
// }

// [[Rcpp::export(.e_LevenbergMsolveCpp)]]
SEXP e_LevenbergMsolveCpp( SEXP AA, SEXP wwAugz, SEXP dd ){ // 
  if (printDebug)   Rcout <<"debut e_LevenbergMsolveCpp()"<<std::endl;
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  const Map<VectorXd> LM_wAugz(as<Map<VectorXd> >(wwAugz));
  const double damping(as<double>(dd));
  const int nc(A.cols());
  const int nr(A.rows());
  VectorXd dampDpD = damping * A.cwiseProduct(A).colwise().sum(); //damping * colSums(wAugX*wAugX)
  VectorXd corrD = dampDpD.cwiseSqrt(); // sqrt(dampDpD);
  //
  VectorXd z(nr+nc);
  z.fill(0);
  const Eigen::HouseholderQR<MatrixXd> QR_A(A); //QR_A <- Rcpp_QR(wAugX)
  z.head(nr) = MatrixXd(QR_A.householderQ()).transpose() * LM_wAugz; // t(QR_A$Q) %*% LM_wAugz    
  //
  MatrixXd RD(nr+nc, nc);
  RD << MatrixXd(QR_A.matrixQR().triangularView<Upper>()),
        MatrixXd(corrD.asDiagonal()); // Matrix() essential! //  (rbind(QR_A$R,diag(corrD,nrow=length(corrD))))
  const Eigen::HouseholderQR<MatrixXd> QR_RD(RD);
  z = MatrixXd(QR_RD.householderQ()).leftCols(nc).transpose() * z;  //t(QR_QR_R$Q[,1:nc]) %*% z
  VectorXd dbetaV = QR_RD.matrixQR().topRows(nc).triangularView<Upper>().solve(z); //backsolve(QR_QR_R$R[1:nc,],z)
  if (printDebug)   Rcout <<"fin e_LevenbergMsolveCpp()"<<std::endl;
  return List::create(Named("dbetaV") = dbetaV,Named("dampDpD")=dampDpD);
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
  return List::create(Named("u") = es.eigenvectors(),Named("d")=es.eigenvalues());
}

// unfortunately useless: // but some improper declarations, would need debugging)
// [[Rcpp::export(.LevPerturbedQCpp)]]
SEXP LevPerturbedQCpp( SEXP perturbedwAugX, SEXP pforREML, SEXP RpRu, SEXP RpRd, SEXP lamOverLam0, SEXP phiOverPhi0 ){
  const Map<MatrixXd> A(as<Map<MatrixXd> >(perturbedwAugX));
  const unsigned int p(as<unsigned int>(pforREML));
  const  Map<MatrixXd> RpRU(as< Map<MatrixXd> >(RpRu));
  const  Map<VectorXd> RpRdiag(as< Map<VectorXd> >(RpRd));
  const double a(sqrt(as<double>(lamOverLam0)));
  const double b(sqrt(as<double>(phiOverPhi0)));
  double psi = b/a; psi *=psi; psi -=1; // (b/a)^2 -1
  int n(RpRdiag.size());
  VectorXd W(n);
  for (int it=0; it<W.size();it++) {W(it) = b*b/(RpRdiag(it)+psi);}
  MatrixXd BB(n,n); 
  BB = RpRU * W.asDiagonal() * RpRU.transpose();  // FR->FR can be improved ? 
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
