#include <Rcpp.h>
#include <gsl/gsl_sf_bessel.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export(.nuln_plus_bessel_lnKnu)]]
NumericVector nuln_plus_bessel_lnKnu( Rcpp::NumericVector x, double nu) { //computes log BesselK_nu
  size_t len=x.size();
  NumericVector val(len);
  for(size_t i = 0; i< len ; i++){
    if (std::isinf(x[i])) {
      val[i]= - std::numeric_limits<double>::infinity();
    } else {
      val[i] = gsl_sf_bessel_lnKnu(nu, x[i]) +nu*log(x[i]);
    }
  }
  return(wrap(val));
}
