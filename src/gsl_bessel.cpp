/*
 This file is part of the spaMM package for R, distributed under 
 the terms of the Cecill-2 licence. 

 It contains routines dealing with evaluation of Bessel function.
 
 These routines are taken, essentially verbatim, from a circa-2012 version of 
 the GNU GSL library, distributed under the terms of the GPL, version 2 or later).
 */


#include "gsl_bessel.h"
#include <Rcpp.h>

using namespace Rcpp;

static inline int cheb_eval_e(const cheb_series * cs,
            const double x,
            gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  {
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
  }

  result->val = d;
  result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

  return GSL_SUCCESS;
}

/* nu = (x+1)/4, -1<x<1, 1/(2nu)(1/Gamma[1-nu]-1/Gamma[1+nu]) */
static double g1_dat[14] = {
  -1.14516408366268311786898152867,
   0.00636085311347084238122955495,
   0.00186245193007206848934643657,
   0.000152833085873453507081227824,
   0.000017017464011802038795324732,
  -6.4597502923347254354668326451e-07,
  -5.1819848432519380894104312968e-08,
   4.5189092894858183051123180797e-10,
   3.2433227371020873043666259180e-11,
   6.8309434024947522875432400828e-13,
   2.8353502755172101513119628130e-14,
  -7.9883905769323592875638087541e-16,
  -3.3726677300771949833341213457e-17,
  -3.6586334809210520744054437104e-20
};
static cheb_series g1_cs = {
  g1_dat,
  13,
  -1, 1,
  7
};

/* nu = (x+1)/4, -1<x<1,  1/2 (1/Gamma[1-nu]+1/Gamma[1+nu]) */
static double g2_dat[15] =
{
  1.882645524949671835019616975350,
 -0.077490658396167518329547945212,
 -0.018256714847324929419579340950,
  0.0006338030209074895795923971731,
  0.0000762290543508729021194461175,
 -9.5501647561720443519853993526e-07,
 -8.8927268107886351912431512955e-08,
 -1.9521334772319613740511880132e-09,
 -9.4003052735885162111769579771e-11,
  4.6875133849532393179290879101e-12,
  2.2658535746925759582447545145e-13,
 -1.1725509698488015111878735251e-15,
 -7.0441338200245222530843155877e-17,
 -2.4377878310107693650659740228e-18,
 -7.5225243218253901727164675011e-20
};
static cheb_series g2_cs = {
  g2_dat,
  14,
  -1, 1,
  8
};

/* based on SLATEC besi0 */

/* chebyshev expansions

 series for bi0        on the interval  0.          to  9.00000d+00
                                        with weighted error   2.46e-18
                                         log weighted error  17.61
                               significant figures required  17.90
                                    decimal places required  18.15

 series for ai0        on the interval  1.25000d-01 to  3.33333d-01
                                        with weighted error   7.87e-17
                                         log weighted error  16.10
                               significant figures required  14.69
                                    decimal places required  16.76


 series for ai02       on the interval  0.          to  1.25000d-01
                                        with weighted error   3.79e-17
                                         log weighted error  16.42
                               significant figures required  14.86
                                    decimal places required  17.09
*/

static double bi0_data[12] = {
  -.07660547252839144951,
  1.92733795399380827000,
   .22826445869203013390,
   .01304891466707290428,
   .00043442709008164874,
   .00000942265768600193,
   .00000014340062895106,
   .00000000161384906966,
   .00000000001396650044,
   .00000000000009579451,
   .00000000000000053339,
   .00000000000000000245
};
static cheb_series bi0_cs = {
  bi0_data,
  11,
  -1, 1,
  11
};

static double ai0_data[21] = {
   .07575994494023796,
   .00759138081082334,
   .00041531313389237,
   .00001070076463439,
  -.00000790117997921,
  -.00000078261435014,
   .00000027838499429,
   .00000000825247260,
  -.00000001204463945,
   .00000000155964859,
   .00000000022925563,
  -.00000000011916228,
   .00000000001757854,
   .00000000000112822,
  -.00000000000114684,
   .00000000000027155,
  -.00000000000002415,
  -.00000000000000608,
   .00000000000000314,
  -.00000000000000071,
   .00000000000000007
};
static cheb_series ai0_cs = {
  ai0_data,
  20,
  -1, 1,
  13
};

static double ai02_data[22] = {
   .05449041101410882,
   .00336911647825569,
   .00006889758346918,
   .00000289137052082,
   .00000020489185893,
   .00000002266668991,
   .00000000339623203,
   .00000000049406022,
   .00000000001188914,
  -.00000000003149915,
  -.00000000001321580,
  -.00000000000179419,
   .00000000000071801,
   .00000000000038529,
   .00000000000001539,
  -.00000000000004151,
  -.00000000000000954,
   .00000000000000382,
   .00000000000000176,
  -.00000000000000034,
  -.00000000000000027,
   .00000000000000003
};
static cheb_series ai02_cs = {
  ai02_data,
  21,
  -1, 1,
  11
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_I0_scaled_e(const double x, gsl_sf_result * result)
{
  double y = fabs(x);

  /* CHECK_POINTER(result) */

  if(y < 2.0 * GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0 - y;
    result->err = 0.5*y*y;
    return GSL_SUCCESS;
  }
  else if(y <= 3.0) {
    const double ey = exp(-y);
    gsl_sf_result c;
    cheb_eval_e(&bi0_cs, y*y/4.5-1.0, &c);
    result->val = ey * (2.75 + c.val);
    result->err = GSL_DBL_EPSILON * fabs(result->val) + ey * c.err;
    return GSL_SUCCESS;
  }
  else if(y <= 8.0) {
    const double sy = sqrt(y);
    gsl_sf_result c;
    cheb_eval_e(&ai0_cs, (48.0/y-11.0)/5.0, &c);
    result->val  = (0.375 + c.val) / sy;
    result->err  = 2.0 * GSL_DBL_EPSILON * (0.375 + fabs(c.val)) / sy;
    result->err += c.err / sy;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double sy = sqrt(y);
    gsl_sf_result c;
    cheb_eval_e(&ai02_cs, 16.0/y-1.0, &c);
    result->val = (0.375 + c.val) / sy;
    result->err  = 2.0 * GSL_DBL_EPSILON * (0.375 + fabs(c.val)) / sy;
    result->err += c.err / sy;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_I0_e(const double x, gsl_sf_result * result)
{
  double y = fabs(x);

  /* CHECK_POINTER(result) */

  if(y < 2.0 * GSL_SQRT_DBL_EPSILON) {
    result->val = 1.0;
    result->err = 0.5*y*y;
    return GSL_SUCCESS;
  }
  else if(y <= 3.0) {
    gsl_sf_result c;
    cheb_eval_e(&bi0_cs, y*y/4.5-1.0, &c);
    result->val  = 2.75 + c.val;
    result->err  = GSL_DBL_EPSILON * (2.75 + fabs(c.val));
    result->err += c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(y < GSL_LOG_DBL_MAX - 1.0) {
    const double ey = exp(y);
    gsl_sf_result b_scaled;
    gsl_sf_bessel_I0_scaled_e(x, &b_scaled);
    result->val  = ey * b_scaled.val;
    result->err  = ey * b_scaled.err + y*GSL_DBL_EPSILON*fabs(result->val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    //OVERFLOW_ERROR(result);
    return OVERFLOW_ERROR;
  }
}

/* based on SLATEC bk0(), bk0e() */

/* chebyshev expansions

 series for bk0        on the interval  0.          to  4.00000d+00
                                        with weighted error   3.57e-19
                                         log weighted error  18.45
                               significant figures required  17.99
                                    decimal places required  18.97

 series for ak0        on the interval  1.25000d-01 to  5.00000d-01
                                        with weighted error   5.34e-17
                                         log weighted error  16.27
                               significant figures required  14.92
                                    decimal places required  16.89

 series for ak02       on the interval  0.          to  1.25000d-01
                                        with weighted error   2.34e-17
                                         log weighted error  16.63
                               significant figures required  14.67
                                    decimal places required  17.20
*/

static double bk0_data[11] = {
  -0.03532739323390276872,
   0.3442898999246284869,
   0.03597993651536150163,
   0.00126461541144692592,
   0.00002286212103119451,
   0.00000025347910790261,
   0.00000000190451637722,
   0.00000000001034969525,
   0.00000000000004259816,
   0.00000000000000013744,
   0.00000000000000000035
};
static cheb_series bk0_cs = {
  bk0_data,
  10,
  -1, 1,
  10
};

static double ak0_data[17] = {
  -0.07643947903327941,
  -0.02235652605699819,
   0.00077341811546938,
  -0.00004281006688886,
   0.00000308170017386,
  -0.00000026393672220,
   0.00000002563713036,
  -0.00000000274270554,
   0.00000000031694296,
  -0.00000000003902353,
   0.00000000000506804,
  -0.00000000000068895,
   0.00000000000009744,
  -0.00000000000001427,
   0.00000000000000215,
  -0.00000000000000033,
   0.00000000000000005
};
static cheb_series ak0_cs = {
  ak0_data,
  16,
  -1, 1,
  10
};

static double ak02_data[14] = {
  -0.01201869826307592,
  -0.00917485269102569,
   0.00014445509317750,
  -0.00000401361417543,
   0.00000015678318108,
  -0.00000000777011043,
   0.00000000046111825,
  -0.00000000003158592,
   0.00000000000243501,
  -0.00000000000020743,
   0.00000000000001925,
  -0.00000000000000192,
   0.00000000000000020,
  -0.00000000000000002
};
static cheb_series ak02_cs = {
  ak02_data,
  13,
  -1, 1,
  8
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K0_scaled_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    //DOMAIN_ERROR(result);
    return DOMAIN_ERROR;
  }
  else if(x <= 2.0) {
    const double lx = log(x);
    const double ex = exp(x);
    int stat_I0;
    gsl_sf_result I0;
    gsl_sf_result c;
    cheb_eval_e(&bk0_cs, 0.5*x*x-1.0, &c);
    stat_I0 = gsl_sf_bessel_I0_e(x, &I0);
    result->val  = ex * ((-lx+M_LN2)*I0.val - 0.25 + c.val);
    result->err  = ex * ((M_LN2+fabs(lx))*I0.err + c.err);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_I0;
  }
  else if(x <= 8.0) {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak0_cs, (16.0/x-5.0)/3.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = c.err / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double sx = sqrt(x);
    gsl_sf_result c;
    cheb_eval_e(&ak02_cs, 16.0/x-1.0, &c);
    result->val  = (1.25 + c.val) / sx;
    result->err  = (c.err + GSL_DBL_EPSILON) / sx;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

/* Evaluate continued fraction CF2, using Thompson-Barnett-Temme method,
 * to obtain values of exp(x)*K_nu and exp(x)*K_{nu+1}.
 *
 * This is unstable for small x; x > 2 is a good cutoff.
 * Also requires |nu| < 1/2.
 */
int
gsl_sf_bessel_K_scaled_steed_temme_CF2(const double nu, const double x,
                                       double * K_nu, double * K_nup1,
                                       double * Kp_nu)
{
  const int maxiter = 10000;

  int i = 1;
  double bi = 2.0*(1.0 + x);
  double di = 1.0/bi;
  double delhi = di;
  double hi    = di;

  double qi   = 0.0;
  double qip1 = 1.0;

  double ai = -(0.25 - nu*nu);
  double a1 = ai;
  double ci = -ai;
  double Qi = -ai;

  double s = 1.0 + Qi*delhi;

  for(i=2; i<=maxiter; i++) {
    double dels;
    double tmp;
    ai -= 2.0*(i-1);
    ci  = -ai*ci/i;
    tmp  = (qi - bi*qip1)/ai;
    qi   = qip1;
    qip1 = tmp;
    Qi += ci*qip1;
    bi += 2.0;
    di  = 1.0/(bi + ai*di);
    delhi = (bi*di - 1.0) * delhi;
    hi += delhi;
    dels = Qi*delhi;
    s += dels;
    if(fabs(dels/s) < GSL_DBL_EPSILON) break;
  }

  hi *= -a1;

  *K_nu   = sqrt(M_PI/(2.0*x)) / s;
  *K_nup1 = *K_nu * (nu + x + 0.5 - hi)/x;
  *Kp_nu  = - *K_nup1 + nu/x * *K_nu;
  if(i == maxiter)
    //GSL_ERROR ("error", GSL_EMAXITER);
    return GSL_ERROR;
  else
    return GSL_SUCCESS;
}

static int gsl_sf_temme_gamma(const double nu, double * g_1pnu, double * g_1mnu, double * g1, double * g2)
{
  const double anu = fabs(nu);    /* functions are even */
  const double x = 4.0*anu - 1.0;
  gsl_sf_result r_g1;
  gsl_sf_result r_g2;
  cheb_eval_e(&g1_cs, x, &r_g1);
  cheb_eval_e(&g2_cs, x, &r_g2);
  *g1 = r_g1.val;
  *g2 = r_g2.val;
  *g_1mnu = 1.0/(r_g2.val + nu * r_g1.val);
  *g_1pnu = 1.0/(r_g2.val - nu * r_g1.val);
  return GSL_SUCCESS;
}


int
gsl_sf_bessel_K_scaled_temme(const double nu, const double x,
                             double * K_nu, double * K_nup1, double * Kp_nu)
{
  const int max_iter = 15000;

  const double half_x    = 0.5 * x;
  const double ln_half_x = log(half_x);
  const double half_x_nu = exp(nu*ln_half_x);
  const double pi_nu   = M_PI * nu;
  const double sigma   = -nu * ln_half_x;
  const double sinrat  = (fabs(pi_nu) < GSL_DBL_EPSILON ? 1.0 : pi_nu/sin(pi_nu));
  const double sinhrat = (fabs(sigma) < GSL_DBL_EPSILON ? 1.0 : sinh(sigma)/sigma);
  const double ex = exp(x);

  double sum0, sum1;
  double fk, pk, qk, hk, ck;
  int k = 0;
  int stat_iter;

  double g_1pnu, g_1mnu, g1, g2;
  int stat_g = gsl_sf_temme_gamma(nu, &g_1pnu, &g_1mnu, &g1, &g2);

  fk = sinrat * (cosh(sigma)*g1 - sinhrat*ln_half_x*g2);
  pk = 0.5/half_x_nu * g_1pnu;
  qk = 0.5*half_x_nu * g_1mnu;
  hk = pk;
  ck = 1.0;
  sum0 = fk;
  sum1 = hk;
  while(k < max_iter) {
    double del0;
    double del1;
    k++;
    fk  = (k*fk + pk + qk)/(k*k-nu*nu);
    ck *= half_x*half_x/k;
    pk /= (k - nu);
    qk /= (k + nu);
    hk  = -k*fk + pk;
    del0 = ck * fk;
    del1 = ck * hk;
    sum0 += del0;
    sum1 += del1;
    if(fabs(del0) < 0.5*fabs(sum0)*GSL_DBL_EPSILON) break;
  }

  *K_nu   = sum0 * ex;
  *K_nup1 = sum1 * 2.0/x * ex;
  *Kp_nu  = - *K_nup1 + nu/x * *K_nu;
#ifdef NO_R_CONSOLE
//    std::cerr<<"bessel_k("<<x<<",nu="<<alpha+nb-1<<"): precision lost in result";
#else
//      REprintf("gsl_sf_bessel_K_scaled_temme : nu=%g,x=%g,K_nu = %g,K_nup1=%g,Kp_nu=%g\n",nu,x,*K_nu,*K_nup1,*Kp_nu);
#endif

  stat_iter = ( k == max_iter ? GSL_EMAXITER : GSL_SUCCESS );
  return GSL_ERROR_SELECT_2(stat_iter, stat_g);
//  return GSL_ERROR_SELECT_2;
}

int
gsl_sf_bessel_Knu_scaled_e10_e(const double nu, const double x, gsl_sf_result_e10 * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0 || nu < 0.0) {
    //DOMAIN_ERROR(result);
    return DOMAIN_ERROR;
  } else {
    int N = (int)(nu + 0.5);
    double mu = nu - N;      /* -1/2 <= mu <= 1/2 */
    double K_mu, K_mup1, Kp_mu;
    double K_nu, K_nup1, K_num1;
    int n, e10 = 0;

    if(x < 2.0) {
      gsl_sf_bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }
    else {
      gsl_sf_bessel_K_scaled_steed_temme_CF2(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }

    /* recurse forward to obtain K_num1, K_nu */
    K_nu   = K_mu;
    K_nup1 = K_mup1;

    for(n=0; n<N; n++) {
      K_num1 = K_nu;
      K_nu   = K_nup1;
      /* rescale the recurrence to avoid overflow */
      if (fabs(K_nu) > GSL_SQRT_DBL_MAX) {
        double p = floor(log(fabs(K_nu))/M_LN10);
        double factor = pow(10.0, p);
        K_num1 /= factor;
        K_nu /= factor;
        e10 += p;
      };
      K_nup1 = 2.0*(mu+n+1)/x * K_nu + K_num1;
    }

    result->val = K_nu;
    result->err = 2.0 * GSL_DBL_EPSILON * (N + 4.0) * fabs(result->val);
    result->e10 = e10;
#ifdef NO_R_CONSOLE
//    std::cerr<<"bessel_k("<<x<<",nu="<<alpha+nb-1<<"): precision lost in result";
#else
//      REprintf("gsl_sf_bessel_Knu_scaled_ez : K_nu = %g,K_mu=%g\n",K_nu, K_mu);
#endif
    return GSL_SUCCESS;
  }
}





int gsl_sf_bessel_lnKnu_e(const double nu, const double x, gsl_sf_result * result) {
  /* CHECK_POINTER(result) */
  if(x <= 0.0 || nu < 0.0) {
      //DOMAIN_ERROR(result);
      return DOMAIN_ERROR;
  } else if(nu == 0.0) {
    gsl_sf_result K_scaled;
    /* This cannot underflow, and
     * it will not throw GSL_EDOM
     * since that is already checked.
     */
    gsl_sf_bessel_K0_scaled_e(x, &K_scaled);
    result->val  = -x + log(fabs(K_scaled.val));
    result->err  = GSL_DBL_EPSILON * fabs(x) + fabs(K_scaled.err/K_scaled.val);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  } else if(x < 2.0 && nu > 1.0) {
    /* Make use of the inequality
     * Knu(x) <= 1/2 (2/x)^nu Gamma(nu),
     * which follows from the integral representation
     * [Abramowitz+Stegun, 9.6.23 (2)]. With this
     * we decide whether or not there is an overflow
     * problem because x is small.
     */
    double ln_bound;
    gsl_sf_result lg_nu;
      ///gsl_sf_lngamma_e(nu, &lg_nu);
      ///vérifier que lgamma de cmath fait la même chose que gsl_sf_lngamma_e
      (&lg_nu)->val = lgamma(nu);
    ln_bound = -M_LN2 - nu*log(0.5*x) + lg_nu.val;
    if(ln_bound > GSL_LOG_DBL_MAX - 20.0) {
      /* x must be very small or nu very large (or both).
       */
      double xi  = 0.25*x*x;
      double sum = 1.0 - xi/(nu-1.0);
      if(nu > 2.0) sum +=  (xi/(nu-1.0)) * (xi/(nu-2.0));
      result->val  = ln_bound + log(sum);
///      result->err  = lg_nu.err; // we have not computed the .err if we have not used gsl_sf_lngamma_e...
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val); ///... hencethis value is not true...
      return GSL_SUCCESS;
    }
#ifdef NO_R_CONSOLE
//    std::cerr<<"bessel_k("<<x<<",nu="<<alpha+nb-1<<"): precision lost in result";
#else
//      REprintf("sous cas x < 2.0 && nu > 1.0 : val = %g,err=%g\n",result->val, result->err);
#endif
    /* can drop-through here */
  }

  {
    /* We passed the above tests, so no problem.
     * Evaluate as usual. Note the possible drop-through
     * in the above code!
     */
    gsl_sf_result_e10 K_scaled;
    //int status =
    gsl_sf_bessel_Knu_scaled_e10_e(nu, x, &K_scaled);
    result->val  = -x + log(fabs(K_scaled.val)) + K_scaled.e10 * M_LN10;
    result->err  = GSL_DBL_EPSILON * fabs(x) + fabs(K_scaled.err/K_scaled.val);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
#ifdef NO_R_CONSOLE
//    std::cerr<<"bessel_k("<<x<<",nu="<<alpha+nb-1<<"): precision lost in result";
#else
//	REprintf("K_scaled.val = %g,K_scaled.err=%g\n",K_scaled.val, K_scaled.err);
#endif

      return GSL_SUCCESS;
  }
}





// [[Rcpp::export(.bessel_lnKnu_e)]]
List bessel_lnKnu_e(Rcpp::NumericVector nu, Rcpp::NumericVector x) { //computes log BesselK_nu
  size_t len=x.size();
  gsl_sf_result result;
  result.val = 0;
  result.err = 0; // both to avoid -Wuninitialized
  IntegerVector status(len);
  NumericVector err(len),val(len);
  for(size_t i = 0; i< len ; i++){
    status[i] = gsl_sf_bessel_lnKnu_e(nu[i], x[i], &result) ;
    val[i] = result.val;
    err[i] = result.err;
  }
  return(List::create(Named("val") = val,Named("err")=err,Named("status")=status));
}
