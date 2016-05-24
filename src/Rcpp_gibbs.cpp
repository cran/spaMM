#include <Rcpp.h>
#include <RcppEigen.h>
#include "Ziggurat.h"  

using namespace Rcpp;
using namespace std;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::export]]
SEXP Rcpp_gibbs_iter(SEXP fix, Rcpp::IntegerVector y, 
                     double lambda, int n_u_h,  
                     SEXP randbGivenObs, 
                     bool ZAisI, double eps,
                     Rcpp::List CondNorm, 
                     SEXP ZA, Rcpp::List forV, SEXP LHSCorrblob, SEXP condL) {
  const Map<VectorXd> fix_(as<Map<VectorXd> >(fix));
  // The maps provide light access without copy, but the accessed copy may be an arglist member for the call.
  // hence its better to always explicitly return values than to play dangerously with in place modifications.
  Map<VectorXd> randbGivenObs_(as<Map<VectorXd> >(randbGivenObs));
  Map<MatrixXd> ZA_(NULL,0,0);
  Map<MatrixXd> u_(NULL,0,0);
  Map<VectorXd> d_(NULL,0);  
  Map<MatrixXd> condLvReg_(NULL,0,0);
  Map<MatrixXd> sqrtCondCovLv_(NULL,0,0);
  Map<MatrixXd> LHSCorrblob_(NULL,0,0);
  Map<MatrixXd> condL_(NULL,0,0);
  if (ZAisI) {
    new (&condLvReg_) Map<MatrixXd>(as<Map<MatrixXd> >(CondNorm["condLvReg"]));
    new (&sqrtCondCovLv_) Map<MatrixXd>(as<Map<MatrixXd> >(CondNorm["sqrtCondCovLv"]));
  } else {
    new (&ZA_) Map<MatrixXd>(as<Map<MatrixXd> >(ZA));
    new (&u_) Map<MatrixXd>(as<Map<MatrixXd> >(forV["u"]));
    new (&d_) Map<VectorXd>(as<Map<VectorXd> >(forV["d"]));
    new (&LHSCorrblob_) Map<MatrixXd>(as<Map<MatrixXd> >(LHSCorrblob));
    new (&condL_) Map<MatrixXd>(as<Map<MatrixXd> >(condL));
  }
  VectorXd randcond_bMean(n_u_h);
  VectorXd randcond_brand(n_u_h);
  VectorXd rnorms(n_u_h);
  int nobs=fix_.size();
  VectorXd locv(nobs);
  VectorXd rand_augY(nobs);
  NumericVector rand_Esp_augY(nobs);
  double mu=0,pn=0;
  GetRNGstate();
  if (ZAisI) {
    rand_Esp_augY = fix_ + randbGivenObs_  ;
  } else rand_Esp_augY = fix_ + ZA_ * randbGivenObs_;  
  for (int i=0; i<rand_Esp_augY.size(); i++) { 
    mu = rand_Esp_augY[i];
    if (y[i]>0) mu=-mu;
    pn = R::unif_rand() * R::pnorm(0.0,mu,1.0,1,0);
    if (pn<=0) pn = eps;
    if (pn>=1) pn = 1.0-eps;
    rand_augY[i] = R::qnorm(pn,mu,1.0,1,0);
    if (y[i]>0) rand_augY[i]=-rand_augY[i];
  }
  if(ZAisI) {
    randcond_bMean = condLvReg_ * (rand_augY-fix_);
    for (int i =0; i<rnorms.size(); i++) rnorms[i]= zigg.norm();
    randcond_brand = sqrtCondCovLv_ * rnorms;
  } else { 
    locv = rand_augY - fix_;
    locv = (locv * u_);
    for (int i =0; i<locv.size(); i++) {locv[i] /= (1/lambda + d_[i]);}
    randcond_bMean = LHSCorrblob_ * locv;
    for (int i =0; i<rnorms.size(); i++) rnorms[i]= zigg.norm();
    randcond_brand = condL_ * rnorms;
  }
  randbGivenObs_ = randcond_bMean + randcond_brand;
  PutRNGstate();  
  return(List::create(_["randbGivenObs"]=randbGivenObs_,
                      _["rand_augY"]=rand_augY,
                      _["randcond_bMean"]=randcond_bMean)); 
}

// [[Rcpp::export]]
SEXP Rcpp_gibbs_debug(SEXP fix, Rcpp::IntegerVector y, 
                      Rcpp::List decomp, bool estimCARrho, 
                      double lambda, int n_u_h,  
                      SEXP randbGivenObs, 
                      bool ZAisI, double eps,
                      Rcpp::List CondNorm, 
                      SEXP ZA, Rcpp::List forV, SEXP LHSCorrblob, SEXP condL) {
  const Map<VectorXd> fix_(as<Map<VectorXd> >(fix));
  const Map<MatrixXd> decomp_u_(as<Map<MatrixXd> >(decomp["u"]));
  const Map<VectorXd> decomp_d_(as<Map<VectorXd> >(decomp["d"]));  
  // The maps provide light access without copy, but the accessed copy may be an arglist member for the call.
  // hence its better to always explicitly return values than to play dangerously with in place modifications.
  Map<VectorXd> randbGivenObs_(as<Map<VectorXd> >(randbGivenObs));
  Map<MatrixXd> ZA_(NULL,0,0);
  Map<MatrixXd> u_(NULL,0,0);
  Map<VectorXd> d_(NULL,0);  
  Map<MatrixXd> condLvReg_(NULL,0,0);
  Map<MatrixXd> sqrtCondCovLv_(NULL,0,0);
  Map<MatrixXd> LHSCorrblob_(NULL,0,0);
  Map<MatrixXd> condL_(NULL,0,0);
  if (ZAisI) {
    new (&condLvReg_) Map<MatrixXd>(as<Map<MatrixXd> >(CondNorm["condLvReg"]));
    new (&sqrtCondCovLv_) Map<MatrixXd>(as<Map<MatrixXd> >(CondNorm["sqrtCondCovLv"]));
  } else {
    new (&ZA_) Map<MatrixXd>(as<Map<MatrixXd> >(ZA));
    new (&u_) Map<MatrixXd>(as<Map<MatrixXd> >(forV["u"]));
    new (&d_) Map<VectorXd>(as<Map<VectorXd> >(forV["d"]));
    new (&LHSCorrblob_) Map<MatrixXd>(as<Map<MatrixXd> >(LHSCorrblob));
    new (&condL_) Map<MatrixXd>(as<Map<MatrixXd> >(condL));
  }
  VectorXd randcond_bMean(n_u_h);
  VectorXd randcond_brand(n_u_h);
  VectorXd rnorms(n_u_h);
  int nobs=fix_.size();
  VectorXd locv(nobs);
  VectorXd rand_augY(nobs);
  NumericVector rand_Esp_augY(nobs);
  double mu=0,pn=0;
  GetRNGstate();
  if (ZAisI) {
    rand_Esp_augY = fix_ + randbGivenObs_  ;
  } else rand_Esp_augY = fix_ + ZA_ * randbGivenObs_;  
  for (int i=0; i<rand_Esp_augY.size(); i++) { 
    mu = rand_Esp_augY[i];
    if (y[i]>0) mu=-mu;
    pn = R::unif_rand() * R::pnorm(0.0,mu,1.0,1,0);
    if (pn<=0) pn = eps;
    if (pn>=1) pn = 1.0-eps;
    rand_augY[i] = R::qnorm(pn,mu,1.0,1,0);
    if (y[i]>0) rand_augY[i]=-rand_augY[i];
  }
  if(ZAisI) {
    randcond_bMean = condLvReg_ * (rand_augY-fix_);
    for (int i =0; i<rnorms.size(); i++) rnorms[i]= zigg.norm();
    randcond_brand = sqrtCondCovLv_ * rnorms;
  } else { 
    locv = rand_augY - fix_;
    locv = (locv * u_);
    for (int i =0; i<locv.size(); i++) {locv[i] /= (1/lambda + d_[i]);}
    randcond_bMean = LHSCorrblob_ * locv;
    for (int i =0; i<rnorms.size(); i++) rnorms[i]= zigg.norm();
    randcond_brand = condL_ * rnorms;
  }
  randbGivenObs_ = randcond_bMean + randcond_brand;
  // D E B U G but valid
  VectorXd tmp(n_u_h); //dummy vector
  VectorXd condv2_or_w2(n_u_h);
  tmp = decomp_u_.transpose() * randbGivenObs_;
  if ( ! estimCARrho) {
    for (int i =0; i<tmp.size(); i++) tmp[i] /= sqrt(decomp_d_[i]);
    tmp = decomp_u_ * tmp;
  } 
  for (int i =0; i<tmp.size(); i++) condv2_or_w2[i] = tmp[i]*tmp[i]; 
  //
  PutRNGstate();  
  return(List::create(_["randbGivenObs"]=randbGivenObs_,
                      _["condv2_or_w2"]=condv2_or_w2,
                      _["rand_augY"]=rand_augY,
                      _["randcond_bMean"]=randcond_bMean)); 
}




// [[Rcpp::export]]
SEXP Rcpp_gibbs(int ngibbs, Rcpp::IntegerVector gibbsSample, bool estimCARrho, Rcpp::List decomp,
                SEXP fix, Rcpp::IntegerVector y, 
                double lambda, int n_u_h,  
                bool ZAisI, double eps,
                Rcpp::List CondNorm, 
                SEXP ZA, Rcpp::List forV, SEXP LHSCorrblob, SEXP condL) {
  const Map<VectorXd> fix_(as<Map<VectorXd> >(fix));
  // The maps provide light access without copy, but the accessed copy may be an arglist member for the call.
  // hence its better to always explicitly return values than to play dangerously with in place modifications.
  const Map<MatrixXd> decomp_u_(as<Map<MatrixXd> >(decomp["u"]));
  const Map<VectorXd> decomp_d_(as<Map<VectorXd> >(decomp["d"]));  
  // Conditionally used variables:
  Map<MatrixXd> ZA_(NULL,0,0);
  Map<MatrixXd> forV_u_(NULL,0,0);
  Map<VectorXd> forV_d_(NULL,0);  
  Map<MatrixXd> condLvReg_(NULL,0,0);
  Map<MatrixXd> sqrtCondCovLv_(NULL,0,0);
  Map<MatrixXd> LHSCorrblob_(NULL,0,0);
  Map<MatrixXd> condL_(NULL,0,0);
  std::sort(gibbsSample.begin(), gibbsSample.end()); // prerequisite for using binary_search later
  //
  MatrixXd tmp(0,0); //dummy matrix
  MatrixXd locv(0,0); //dummy matrix
  VectorXd rnorms(n_u_h);
  VectorXd sumcond_bMeans=VectorXd::Zero(n_u_h);
  VectorXd sum_condv2_or_w2=VectorXd::Zero(n_u_h);
  VectorXd randcond_brand(n_u_h);  
  VectorXd randbGivenObs=VectorXd::Zero(n_u_h);
  VectorXd randcond_bMean(n_u_h);
  int nobs=fix_.size();
  VectorXd rand_augY(nobs);
  VectorXd sum_rand_augY_Means=VectorXd::Zero(nobs);
  VectorXd rand_Esp_augY(nobs);
  double mu=0,pn=0;
  //
  if (ZAisI) {
    new (&condLvReg_) Map<MatrixXd>(as<Map<MatrixXd> >(CondNorm["condLvReg"]));
    new (&sqrtCondCovLv_) Map<MatrixXd>(as<Map<MatrixXd> >(CondNorm["sqrtCondCovLv"]));
    if (condLvReg_.cols()!=n_u_h) { Rcout<<"condLvReg_.cols()!=n_u_h"<<endl; return(List::create(_["error"]=1)); }
    if (condLvReg_.rows()!=n_u_h) { Rcout<<"condLvReg_.rows()!=n_u_h"<<endl; return(List::create(_["error"]=2)); }
    if (sqrtCondCovLv_.cols()!=n_u_h) { Rcout<<"sqrtCondCovLv_.cols()!=n_u_h"<<endl; return(List::create(_["error"]=3)); }
    if (sqrtCondCovLv_.rows()!=n_u_h) { Rcout<<"sqrtCondCovLv_.rows()!=n_u_h"<<endl; return(List::create(_["error"]=4)); }
  } else {
    new (&forV_u_) Map<MatrixXd>(as<Map<MatrixXd> >(forV["u"]));
    new (&forV_d_) Map<VectorXd>(as<Map<VectorXd> >(forV["d"]));
    new (&ZA_) Map<MatrixXd>(as<Map<MatrixXd> >(ZA));
    new (&LHSCorrblob_) Map<MatrixXd>(as<Map<MatrixXd> >(LHSCorrblob));
    new (&condL_) Map<MatrixXd>(as<Map<MatrixXd> >(condL));
    if (forV_u_.cols()!=nobs) { Rcout<<"forV_u_.cols()!=nobs"<<endl; return(List::create(_["error"]=5)); }
    if (forV_u_.rows()!=nobs) { Rcout<<"forV_u_.rows()!=nobs"<<endl; return(List::create(_["error"]=6)); }
    if (forV_d_.size()!=nobs) { Rcout<<"forV_d_.size()!=nobs"<<endl; return(List::create(_["error"]=7)); }
    if (ZA_.cols()!=n_u_h) { Rcout<<"ZA_.cols()!=n_u_h"<<endl; return(List::create(_["error"]=8)); }
    if (ZA_.rows()!=nobs) { Rcout<<"ZA_.rows()!=nobs"<<endl; return(List::create(_["error"]=9)); }
    if (LHSCorrblob_.cols()!=nobs) { Rcout<<"LHSCorrblob_.cols()!=nobs"<<endl; return(List::create(_["error"]=10)); }
    if (LHSCorrblob_.rows()!=n_u_h) { Rcout<<"LHSCorrblob_.rows()!=n_u_h"<<endl; return(List::create(_["error"]=11)); }
    if (condL_.cols()!=n_u_h) { Rcout<<"condL_.cols()!=n_u_h"<<endl; return(List::create(_["error"]=12)); }
    if (condL_.rows()!=n_u_h) { Rcout<<"condL_.rows()!=n_u_h"<<endl; return(List::create(_["error"]=13)); }
  }
  if (y.size()!=nobs) { Rcout<<"y.size()!=nobs"<<endl; return(List::create(_["error"]=14)); }
  if (decomp_u_.cols()!=n_u_h) { Rcout<<"decomp_u_.cols()!=n_u_h"<<endl; return(List::create(_["error"]=15)); }
  if (decomp_u_.rows()!=n_u_h) { Rcout<<"decomp_u_.rows()!=n_u_h"<<endl; return(List::create(_["error"]=16)); }
  if (decomp_d_.size()!=n_u_h) { Rcout<<"decomp_d_.size()!=n_u_h"<<endl; return(List::create(_["error"]=17)); }
  //
  GetRNGstate();
  for (int g=0;g<ngibbs;g++) {
    if (ZAisI) {
      rand_Esp_augY = fix_ + randbGivenObs  ;
    } else rand_Esp_augY = fix_ + ZA_ * randbGivenObs;  
    for (int i=0; i<rand_Esp_augY.size(); i++) { 
      mu = rand_Esp_augY[i];
      if (y[i]>0) mu=-mu;
      pn = R::unif_rand() * R::pnorm(0.0,mu,1.0,1,0);
      if (pn<=0) pn = eps;
      if (pn>=1) pn = 1.0-eps;
      rand_augY[i] = R::qnorm(pn,mu,1.0,1,0);
      if (y[i]>0) rand_augY[i]=-rand_augY[i];
    }
    if(ZAisI) {
      randcond_bMean = condLvReg_ * (rand_augY-fix_);
      for (int i =0; i<rnorms.size(); i++) rnorms[i]= zigg.norm();
      randcond_brand = sqrtCondCovLv_ * rnorms;
    } else { 
      locv = rand_augY - fix_; // nobs x 1
      locv = locv.transpose() * forV_u_; // all (1 x nobs).(nobs x nobs) = (1 x nobs) 
      for (int i =0; i<locv.cols(); i++) {locv(0,i) /= (1/lambda + forV_d_[i]);} // forV_d_ also nobs
      randcond_bMean = LHSCorrblob_ * locv.transpose(); // (n_u_h x nobs).(nobs x 1) = n_u_h x 1
      for (int i =0; i<rnorms.size(); i++) rnorms[i]= zigg.norm(); // n_u_h x 1
      randcond_brand = condL_ * rnorms;  // n_u_h x 1
    }
 //   Rcout << randcond_brand[0] << " " << randcond_brand[55]<<endl;
    randbGivenObs = randcond_bMean + randcond_brand;  //used to compute condv2 or condw2, ./.
    // ./. initiate next iter of this loop and retuned to initiate next iter of SEM main loop 
    if (std::binary_search(gibbsSample.begin(), gibbsSample.end(), g)) {
      sumcond_bMeans += randcond_bMean;
      sum_rand_augY_Means += rand_augY;
      tmp = randbGivenObs.transpose() * decomp_u_; // without the transpose, tmp becomes wrong-dimensioned
      tmp.transposeInPlace(); // back to 1-col vector t(u).v // dont use tmp = tmp.transpose() !
      //tmp = decomp_u_.transpose() * randbGivenObs;
      if ( ! estimCARrho) {
        for (int i =0; i<tmp.rows(); i++) tmp(i,0) /=  sqrt(decomp_d_[i]);
        tmp = decomp_u_ * tmp;
      } 
      for (int i =0; i<tmp.rows(); i++) sum_condv2_or_w2[i] += tmp(i,0)*tmp(i,0); 
      // v2 iid; w2 i_but not_id . for estimation of lambda and rho by GLM
    }
  }
  PutRNGstate();  
  return(List::create(_["sumcond_bMeans"]=sumcond_bMeans,
                      _["sum_rand_augY_Means"]=sum_rand_augY_Means,
                      _["sum_condv2_or_w2"]=sum_condv2_or_w2,
                      _["rand_augY"]=rand_augY // debug
           )); 
}

/*
library(spaMM)
data(scotlip)
 essai <- 1
 set.seed(124)
ldl <- selfAdjointSolverCpp(Nmatrix)
Lmat <- ldl$u %*% diag(sqrt(1/(1-0.17*ldl$d)))
lp <- 0.1 + 3* Lmat %*% rnorm(ncol(Lmat)) ## single intercept beta =0.1; lambda=3
resp <- rbinom(ncol(Lmat),1,1/(1+exp(-lp)))
donn <- data.frame(npos=resp,nneg=1-resp,gridcode=scotlip$gridcode)
 set.seed(essai)
 essai <- essai+1
 HLCor(cbind(npos,nneg)~1 +adjacency(1|gridcode),
      adjMatrix=Nmatrix,family=binomial(probit),data=donn,HLmethod="SEM")
 

*/
