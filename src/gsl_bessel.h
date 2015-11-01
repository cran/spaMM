/*
 This file is part of the spaMM package for R, distributed under 
the terms of the Cecill-2 licence. 

It contains routines dealing with evaluation of Bessel function.

These routines are taken, essentially verbatim, from a circa-2012 version of 
the GNU GSL library, distributed under the terms of the GPL, version 2 or later).
*/


#ifndef H_GSL_BESSEL
#define H_GSL_BESSEL

#include <cmath>
#include <limits>
#include <stdio.h>

//#include "./gsl-1.9/specfunc/chebyshev.h"
//#include "./gsl-1.9/specfunc/cheb_eval.c"

#define GSL_SUCCESS 1
#define GSL_ERROR -1
#define DOMAIN_ERROR -2
/// definition originelle...
#define GSL_ERROR_SELECT_2(a,b) ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
///#define GSL_ERROR_SELECT_2 -3
#define OVERFLOW_ERROR -4


// constantes à vérifier!!
#define GSL_EMAXITER 11 //value as defined in gsl_errno.h
#define GSL_LOG_DBL_MAX 7.0978271289338397e+02 //value as defined in gsl_machine.h
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_DBL_EPSILON std::numeric_limits<double>::epsilon()
#define GSL_SQRT_DBL_EPSILON  sqrt(GSL_DBL_EPSILON)
//#define M_LN10   2.30258509299404568401799145468 /* ln(10) */

struct gsl_sf_result {
  double val;
  double err;
};

struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;

struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

#endif
