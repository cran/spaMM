#include <Rcpp.h>
#include <cmath>
#include <numeric>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export(.Rcpp_COMP_Z)]]
SEXP Rcpp_COMP_Z(int moment,double nu, double lambda, int maxn) {
  //COMP_Z <- function(eta,nu,lambda=exp(eta),maxn=.COMP_maxn(lambda,nu)){
  double logScaleFac=0;
  double scaled=1;
  if (nu<1e-08) {
    switch(moment) {
    case 0: scaled=1/(1-lambda);break;
    case 1: scaled=lambda/pow(1.-lambda,2);break;
    case 2: scaled=lambda*(1+lambda)/pow(1.-lambda,3);break;
    case 3: scaled=lambda*(1+4*lambda+lambda*lambda)/pow(1.-lambda,4);break;
    }
  } else {
    int floorn = floor(maxn);
    double epsn = maxn - floorn;
    double jd; 
    vector<double> facs;
    if (moment==0) {
      facs.resize(floorn+2);
      facs[0] = 1;
      for(unsigned int i=1;i < facs.size();i++) facs[i]=lambda/pow(i,nu);
    } else {
      facs.resize(floorn+1);// implicit first term is 0 does not contribute to sums...
      facs[0] = lambda;
      for(unsigned int i=1;i < facs.size();i++) {
        jd=i+1.; // cast i before computing ratio...
        facs[i]=pow(jd/(jd-1),moment)*(lambda/pow(jd,nu)); 
      }
    }
    vector<double> cumprodfacs(facs.size());
    std::partial_sum(facs.begin(), facs.end(), cumprodfacs.begin(), multiplies<double>()); 
    (cumprodfacs.back()) *= epsn;
    scaled = std::accumulate(cumprodfacs.begin(), cumprodfacs.end(), double(0));
    if ( ! R_FINITE(scaled)) { // Writing R xtensions 6.4 & 1.6.4
      for(vector<double>::iterator it=facs.begin();it !=facs.end();it++) (*it)=log(*it);
      vector<double> cumsumlogfacs(facs.size());
      std::partial_sum(facs.begin(), facs.end(), cumsumlogfacs.begin());
      logScaleFac = *(std::max_element(cumsumlogfacs.begin(),cumsumlogfacs.end()));
      for(vector<double>::iterator it=cumsumlogfacs.begin();it !=cumsumlogfacs.end();it++) (*it) -= logScaleFac;
      for(vector<double>::iterator it=cumsumlogfacs.begin();it !=cumsumlogfacs.end();it++) (*it)=exp(*it);
      scaled = std::accumulate(cumsumlogfacs.begin(), cumsumlogfacs.end(), double(0));
    }
  }
  NumericVector comp_z= NumericVector::create(_["logScaleFac"]=logScaleFac, _["scaled"]=scaled);
  return(comp_z);
}

// [[Rcpp::export(.COMP_Z_integrand)]]
SEXP COMP_Z_integrand(Rcpp::NumericVector z, 
                      Nullable<NumericVector> eta = R_NilValue, Nullable<NumericVector> lambda = R_NilValue , // allows _missing_ values 
                      double nu=0.5, //needs default value since after arg with default value
                      int moment=0, double logScaleFac=0) {
  double eta_;
  double lambda_;
  if (eta.isNull()) {
    lambda_ = (NumericVector(lambda))[0];
    eta_ = log(lambda_);
  } else eta_ = (NumericVector(eta))[0];
  Rcpp::NumericVector logint = moment*log(z)+ z*eta_ - nu*lgamma(z+1);
  Rcpp::NumericVector res = exp(logint-logScaleFac);
  // if (diagnostic) {
  //   return(List::create(Named("z")=z,
  //                       Named("eta_")=eta_,
  //                       Named("logint")=logint,
  //                       Named("beg")=moment*log(z),
  //                       Named("beg2")=z*eta_,
  //                       Named("end")=nu*lgamma(z+1),
  //                       Named("logScaleFac")=logScaleFac,
  //                       Named("raw")=exp(logint-logScaleFac)));
  // } else {
    res = pmin(DBL_MAX, res);
    return(wrap(res));
  // }
}

// [[Rcpp::export(.Rcpp_COMP_Z_asympto)]]
SEXP Rcpp_COMP_Z_asympto(double nu, double pow_lam_nu) {
  double logScaleFac=nu*pow_lam_nu;
  // using Gaunt et al.:
  double nu2 = nu*nu;
  // double nu4 = nu2*nu2;
  // double nu6 = nu4*nu2;
  // double nu8 = nu4*nu4;
  // double nu10 = nu6*nu4;
  double c1 = (nu2-1)/24;
  double c2 = (nu2-1)*(nu2+23)/1152;
  // double c3 = (nu2-1)*(5*nu4-298*nu2+11237)/414720;
  // double c4 = (nu2-1)*(5*nu6-1887*nu4-241041*nu2+2482411)/39813120;
  // double c5 = (nu2-1)*(7*nu8-7420*nu6+1451274*nu4-220083004*nu2+1363929895)/6688604160;
  // double c6 = (nu2-1)*(35*nu10 - 78295*nu8 + 76299326*nu6 + 25171388146*nu4 
  //                   - 915974552561*nu2 + 4175309343349)/4815794995200;
  // double c7 = (nu2-1)*(5*nu6*nu6- 20190*nu10 + 45700491*nu8 - 19956117988*nu6
  //                   +7134232164555*nu4-142838662997982*nu2 + 525035501918789)/115579079884800;
  // (1+c1/logScaleFac+c2/logScaleFac^2+c3/logScaleFac^3+c4/logScaleFac^4+c5/logScaleFac^5+c6/logScaleFac^6+c7/logScaleFac^7)
  double invLSF= 1/logScaleFac;
  //double scaled = 1 + invLSF*(c1+invLSF*(c2+invLSF*(c3+invLSF*(c4+invLSF*(c5+invLSF*(c6+invLSF*c7))))));
  double scaled = 1 + invLSF*(c1+invLSF*c2);
  scaled = scaled/(pow(2*M_PI*pow_lam_nu,(nu-1)/2)*sqrt(nu));
  NumericVector comp_z= NumericVector::create(_["logScaleFac"]=logScaleFac, _["scaled"]=scaled);
  return(comp_z);
}
