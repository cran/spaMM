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
    if (ISNAN(scaled)) { // Writing R xtensions 1.6.4
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


