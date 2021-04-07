#include "internals.h"
using namespace Rcpp;

// [[Rcpp::export(`.lower.tri<-`)]]
SEXP lowertri(NumericMatrix A, bool diag, NumericVector value) { // for use as .lower.tri(., diag=..) <- ...
  int i, k=0, nc=A.ncol();
  for (int j=0; j<nc; j++) {
    if (diag) {
      A(j,j)=value[k];
      k++;
    }
    for (i=j+1; i<nc; i++) {
      A(i,j)=value[k];
      k++;
    }
  }
  return wrap(A);
}

// [[Rcpp::export(`.upper.tri<-`)]]
SEXP uppertri(NumericMatrix A, bool diag, NumericVector value) { // for use as .upper.tri(., diag=..) <- ...
  int i, k=0, nc=A.ncol();
  for (int j=0; j<nc; j++) {
    for (i=0; i<j; i++) {
      A(i,j)=value[k];
      k++;
    }
    if (diag) {
      A(j,j)=value[k];
      k++;
    }
  }
  return wrap(A);
}

// [[Rcpp::export(.rC_inv_chol_cpp)]]
SEXP rC_inv_chol_cpp(SEXP trCoef) { // C++ version of ..ranCoefsInv(., rC_transf="chol")
  // .calc_cov_from_trRancoef() part:
  VectorXd trRancoef(as<VectorXd>(trCoef));
  RObject ncol_attr = Rf_getAttrib(trCoef, wrap("Xi_ncol"));
  int Xi_ncol;
  if (is<IntegerVector>(ncol_attr)) { 
    IntegerVector zut=as<IntegerVector>(ncol_attr);
    Xi_ncol=zut[0];
  } else Xi_ncol=std::floor(std::sqrt(trRancoef.size()*2));
  
  MatrixXd covmat=MatrixXd(Xi_ncol, Xi_ncol).setZero();
  int k=0, i;
  for (int j=0; j<Xi_ncol; j++) {
    for (i=0; i<j; i++) {
      covmat(i,j)=trRancoef[k];
      k++;
    }
    covmat(j,j)=std::max(trRancoef[k],1e-32);
    k++;
  }
  covmat = MatrixXd(Xi_ncol, Xi_ncol).setZero().selfadjointView<Lower>().rankUpdate(covmat.transpose());
  
  //// .smooth_regul()
  //  .eigen_sym():
  Function f("eigen");   
  List es=f(wrap(covmat), Named("symmetric")=true);
  VectorXd values(as<VectorXd>(es["values"]));
  Map<MatrixXd> vectors(as<Map<MatrixXd> >(es["vectors"]));
  double v;
  const double epsi=1e-08;
  
  for(int it=0; it < values.size(); it++ ) {
    v=values(it);
    if (v<2*epsi) values(it)=4*epsi*v*v+epsi; // inverse is w[w<2*epsi] <- sqrt( (4*epsi)*(w[w<2*epsi]-epsi) )
  }
  VectorXd sqrtW = values.array().sqrt();
  MatrixXd rancoefs = vectors * sqrtW.asDiagonal();
  const int r(vectors.rows());
  rancoefs = MatrixXd(r, r).setZero().selfadjointView<Lower>().rankUpdate(rancoefs);// as in tcrossprod
  //return List::create(Named("ZWZt") = swZ,Named("vectors") = vectors,Named("values")=es["values"], Named("d_regul")=values);
  
  // lambdas <- rancoefs[diagPos];     torancoefs <- sqrt(1/lambdas)
  VectorXd lambdas=rancoefs.diagonal();
  VectorXd torancoefs=lambdas.array().sqrt().inverse();
  
  //    rancoefs <- torancoefs * rancoefs * rep(torancoefs, each = Xi_ncol) # cf cov2cor()
  rancoefs = torancoefs.asDiagonal() * rancoefs * torancoefs.asDiagonal();
  rancoefs.diagonal()=lambdas;
  
  k=0;
  for (int j=0; j<Xi_ncol; j++) {
    trRancoef[k]=rancoefs(j,j);
    k++;
    for (i=j+1; i<Xi_ncol; i++) {
      trRancoef[k]=rancoefs(i,j);
      k++;
    }
  }
  
  SEXP sexp = Rcpp::wrap(trRancoef); // path from Eigen object to Rcpp object
  Rcpp::NumericVector out(sexp); // Rcpp object => attributes can be added
  out.attr("transf")= trCoef;
  out.attr("Xi_ncol")=Xi_ncol;
  return wrap(out);
}

// [[Rcpp::export(.C_calc_cov_from_ranCoef)]]
SEXP C_calc_cov_from_ranCoef(SEXP ranCoef, Nullable<IntegerVector> Xi_ncol = R_NilValue) {
  Rcpp::IntegerVector zut;
  if (Xi_ncol.isNull()) {
    RObject ncol_attr = Rf_getAttrib(ranCoef, wrap("Xi_ncol"));
    zut=as<Rcpp::IntegerVector>(ncol_attr);
  } else zut=IntegerVector(Xi_ncol);
  int nc=zut[0];
  Map<VectorXd> rC(as<Map<VectorXd> >(ranCoef)); // fails if input is integer vector
  int i, k=0;
  MatrixXd compactcovmat=MatrixXd(nc,nc).setIdentity();
  VectorXd sigmas=VectorXd(nc);
  for (int j=0; j<nc; j++) {
    sigmas[j]=std::min(1e6,sqrt(rC[k]));
    k++;
    for (i=j+1; i<nc; i++) {
      compactcovmat(i,j)=rC[k];
      compactcovmat(j,i)=rC[k];
      k++;
    }
  }
  compactcovmat = sigmas.asDiagonal() * compactcovmat * sigmas.asDiagonal();  
  
  return wrap(compactcovmat);
} 

// [[Rcpp::export(.C_makelong)]]
SEXP makelong(NumericMatrix Lcompact, int Lsize) { // first filling of matrix, without template
  //NumericMatrix Lcompact(lcompact);
  //int Lsize = Rcpp::as<int>(longsize);
  
  Eigen::SparseMatrix<double> longLv(Lsize,Lsize); // resize necessary before rankUpdate
  int nc = Lcompact.ncol();
  int n_levels=Lsize/nc, ubeg1, ubeg2, jt, lt;
  //longLv.setZero();
  longLv.reserve(Eigen::VectorXi::Constant(Lsize,nc));
  double Cij;
  for (int it=0; it<nc; it++) {
    ubeg1 = it*n_levels;
    for (jt=0; jt<nc; jt++) {
      Cij=Lcompact(it,jt);
      ubeg2 = jt*n_levels;
      for (lt=0; lt<n_levels; lt++) {
        longLv.insert(ubeg1+lt,ubeg2+lt)=Cij;
      }
    }
  }
  longLv.makeCompressed(); // essential for correct conversion to dgC
  return wrap(longLv);
} 


// [[Rcpp::export(.C_makelong2)]]
SEXP makelong2(NumericMatrix Lcompact, int Lsize) { // Pedagogical: same result as makelong() but using more general algo. 
  //NumericMatrix Lcompact(lcompact);
  //int Lsize = Rcpp::as<int>(longsize);
  
  typedef Eigen::Triplet<double> T;
  
  Eigen::SparseMatrix<double> longLv(Lsize,Lsize); 
  int nc = Lcompact.ncol();
  int n_levels=Lsize/nc, ubeg1, ubeg2, jt, lt;
  double Cij;
  
  std::vector<T> tripletList;
  tripletList.reserve(Lsize*nc);
  for (int it=0; it<nc; it++) {
    ubeg1 = it*n_levels;
    for (jt=0; jt<nc; jt++) {
      Cij=Lcompact(it,jt);
      ubeg2 = jt*n_levels;
      for (lt=0; lt<n_levels; lt++) {
        tripletList.push_back(T(ubeg1+lt,ubeg2+lt,Cij));
      }
    }
  }
  
  longLv.setFromTriplets(tripletList.begin(), tripletList.end()); // directly in compressed mode
  return wrap(longLv);
}

// [[Rcpp::export(.C_dispInv)]]
SEXP C_dispInv(NumericVector x) {
  const double xref=5e-5,xreff=1e-1,xm=1e-3, thr=std::log(xref+xm), fac=(xm+xref)/(xm+xreff);
  double z;
  for (int it=0; it<x.size(); it++) {
    z=x[it];
    if (z<thr) {
      x[it]=std::exp(z)-xref;
    } else {
      z=std::exp((z-thr)*fac);
      x[it]=xm*z+xreff*(z-1);
    }
  }
  return(wrap(x));
}

///// logit 
// cf r-source/src/library/stats/src/family.c 
static inline double x_d_omx(double x) {return x/(1 - x);}

// [[Rcpp::export(.logit)]]
SEXP logit(NumericVector mu) {
  for (int it=0; it<mu.size(); it++) mu[it]=log(x_d_omx(mu[it]));
  return(wrap(mu));
}

// [[Rcpp::export(.is_evaluated)]]
bool is_evaluated(Symbol name, Environment env) {
  SEXP object = Rf_findVar(name, env);
  return PRVALUE(object) != R_UnboundValue;
}

// [[Rcpp::export(.is_promise)]]
bool is_promise2(Symbol name, Environment env) {
  SEXP object = Rf_findVar(name, env);
  return (TYPEOF (object) == PROMSXP);
}
