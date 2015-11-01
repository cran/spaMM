// cf spaMM/devel/RcppDevel.R to generate this
#include "spaMM_linear.h"

using namespace Rcpp ;

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
SEXP lmwithSparseQ( SEXP XX, SEXP yy ){
  const MappedSparseMatrix<double> X(as<MappedSparseMatrix<double> >(XX));
  const Map<VectorXd> y(as<Map<VectorXd> >(yy));
  SparseQR<SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > spQR(X);
  // cf http://eigen.tuxfamily.org/bz/show_bug.cgi?id=836
  //SparseMatrix<double> Q_ap; 
  SparseMatrix<double> Q_ap; 
  Q_ap = spQR.matrixQ(); // AP = Q_ap R_ap here, ie A = Q_ap.tP . P.R_ap.tP 
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(spQR.colsPermutation());
  return List::create(Named("coef") = VectorXd(spQR.solve(y)),
                      Named("P") = perm.indices(), // A= $Q_ap $R_ap[,sort.list($P)] as in ?qr.X = $Q_ap[,sort.list($P)] $R_ap[sort.list($P),sort.list($P)]
                      Named("Q_ap") =MatrixXd(Q_ap)); // full Q // returning as MatrixXd seems faster ??
//                      Named("Q_ap") =Q_ap); // full Q 
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
  return wrap(MatrixXd(QR.householderQ())); // returns full Q matrix
}

// [[Rcpp::export]]
SEXP Rcpp_QR( SEXP XX){
  const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
  const Eigen::HouseholderQR<MatrixXd> QR(X);
  List out = List::create(Named("Q") = MatrixXd(QR.householderQ()), // full Q matrix (cf leverages compute for thin Q)
                      Named("R") = MatrixXd(QR.matrixQR().triangularView<Upper>())); 
  out.attr("class") = "Rcpp_QR";
  return(out);
}

// [[Rcpp::export]]
SEXP Rcpp_sparseQR( SEXP XX){
  // well, hell http://eigen.tuxfamily.org/bz/show_bug.cgi?id=706 - Allowing SparseQR and SparseCholesky to work with MappedSparseMatrix
  // mais en remplaçant   const MappedSparseMatrix<double> X(as<MappedSparseMatrix<double> >(XX));   on ne résout pas le pb: 
  // (Error in validObject(x) : invalid class “dgCMatrix” object: row indices are not sorted within columns) on $R (depend de la taille de XX)
  const MappedSparseMatrix<double> X(as<MappedSparseMatrix<double> >(XX));

  SparseQR<SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > spQR(X);
  // cf http://eigen.tuxfamily.org/bz/show_bug.cgi?id=836
  SparseMatrix<double> Q_ap; 
  Q_ap = spQR.matrixQ(); // AP = Q_ap R_ap here, ie A = Q_ap.tP . P.R_ap.tP 
  // = Q_a R_a for Q_a= Q_ap.tP and R_a= P.R_ap.tP          (1) "no pivoting" method
  // or = Q_ap Z_a for Q_a = Q_ap and Z_a= R_ap.tP           (2)       Z_a not unpper triangular
  // since we will use Q_a t(Q_a)= Q_ap t(Q_ap) we dont need to compute Q_a;
  // with R_a (but not Z_a) we could use backsolve (matrix, not Matrix!) which would be useful if R is not sparse...
  // we would further need to have Q.tP as Q * tP and <<*>> is defined in principle for any matrix type, but
  // <<SparseQRMatrixQReturnType should inherit EigenBase and provide an evalTo member to be converted to a plain SparseMatrix.>>
  // http://eigen.tuxfamily.org/bz/show_bug.cgi?id=596 ; which means that it doesn t inherits
  // de meme je n arrive pas a premul par perm une matrice declaree comme SparseMatrix<double>...  mais je m en passe
  const int r(spQR.rank());
  SparseMatrix<double> R_ap = (spQR.matrixR().topLeftCorner(r,r));  
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(spQR.colsPermutation());
  // pour comparaison onpourrait exporter
  //SparseMatrix<double> AP;
  //AP = (X * spQR.colsPermutation()).sparseView();
  List out = List::create(Named("Q_ap") = Q_ap, // full Q matrix !
                          Named("P") = perm.indices(), // A= $Q_ap $R_ap[,sort.list($P)] as in ?qr.X = $Q_ap[,sort.list($P)] $R_ap[sort.list($P),sort.list($P)]
                          Named("R_ap") = R_ap); 
  out.attr("class") = "Rcpp_sparseQR";
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
//return wrap(Z * W.asDiagonal() *Z.adjoint());
VectorXd sqrtW=W.array().sqrt();
MatrixXd swZ= Z * sqrtW.asDiagonal();
const int n(swZ.rows());
swZ = MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(swZ);// as in tcrossprod
return wrap(swZ); 
}

// [[Rcpp::export]]
SEXP ZtWZ( SEXP ZZ, SEXP WW ){
Map<MatrixXd> Z(as<Map<MatrixXd> >(ZZ));
const Map<VectorXd> W(as<Map<VectorXd> >(WW));
//return wrap(Z.adjoint() * W.asDiagonal() *Z);
VectorXd sqrtW=W.array().sqrt();
MatrixXd swZ= sqrtW.asDiagonal() *Z;
const int n(swZ.cols());
swZ = MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(swZ.transpose());// as in crossprod
return wrap(swZ); 
}


// [[Rcpp::export]]
SEXP RcppChol( SEXP AA ){ // returns t(R::chol)
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
  return out;
} // such that A = L.Lt ie standard, NOT R's chol()

// https://gist.github.com/bobthecat/6509321
// [[Rcpp::export]]
SEXP crossprodCpp(SEXP Mat){
const Map<MatrixXd> A(as<Map<MatrixXd> >(Mat));
const int n(A.cols());
MatrixXd tAA(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A.transpose()));
return wrap(tAA);
} // correspond bien a crossprod()

// [[Rcpp::export]]
SEXP tcrossprodCpp(SEXP Mat){
const Map<MatrixXd> A(as<Map<MatrixXd> >(Mat));
const int n(A.rows());
MatrixXd AAt(MatrixXd(n, n).setZero().selfadjointView<Lower>().rankUpdate(A));
return wrap(AAt);
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
  //const Eigen::HouseholderQR<MatrixXd> PQR(AtAdDpD);
  const Eigen::LDLT<MatrixXd> PQR(AtAdDpD); // (LDLt is a version of cholesky) faster than QR
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

// unfortunately useless:
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
 MatrixXd BB = RpRU * W.asDiagonal() * RpRU.transpose();  // FR->FR can be improved ? 
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
