#include "Ziggurat.h"  
#include <Rcpp.h>
using namespace Rcpp;


static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::export]]
Rcpp::NumericVector zrnorm(int n) {
  Rcpp::NumericVector x(n);
  for (int i=0; i<n; i++) {
    x[i] = zigg.norm();
  }
  return x;
}
