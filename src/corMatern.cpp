/*
 This file is part of the spaMM package for R, distributed under 
 the terms of the Cecill-2 licence. 
 
 It contains routines dealing with corMatern correlation structures.
 
 It is derived from routines dealing with corSpatial correlation structures, 
 that were part of a circa 2012 version of the nlme package for R,
 made available under the terms of the GNU General Public
 License, version 2, The authors of these routines were described as 
  " Copyright 1997-2005  Douglas M. Bates <bates@stat.wisc.edu>,
			Jose C. Pinheiro,
			Saikat DebRoy " in the source used in this way.

*/

#include <limits>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#define R_NO_REMAP
#include "R.h" // defines USING_R
#include "Rmath.h" // bessel fns

#define longint int

double patched_bessel_k(double x, double alpha, double expo) {
  if (alpha >1000) alpha=1000; /// this is quick patch for addressing problems with glmmPQL -> nlminb -> reaches 'here' with huge nu values -> pb with the calloc /// FR 250113
  /** note vignette of http://cran.r-project.org/web/packages/Bessel/ . But the source code is not appealing, and for large nu in particular **/
  return(bessel_k(x, alpha, expo));
}


/* methods for the virtual class */

extern "C" { //computes Matérn correlation for one distance
  void matern_cor(double *par, double *dist, longint *nug, longint *nusc,double *cor)
{
  double dscale,sc,aux,con,ratio = 1.0;
  if (*nug) ratio = par[2];
  if (*nusc) {
    sc=2.0*sqrt(par[1])/par[0];
    // when nu (par[1]) is large this compares to corGaus
    // with range rho (par[0])
  } else {
    sc=1.0/par[0];
    // initial formulation: when nu is large this compares
    // to corGaus with range par[0]*2*sqrt(nu)
  }
  con = pow(2.0,par[1] - 1) * exp(lgamma(par[1]));
  con = 1.0/con;
  if (*dist < std::numeric_limits<double>::epsilon())  {               // distance ~ 0
    cor[0] = ratio;                                               // correlation is one!
  } else {                                                        // distance > 0
    dscale = *dist * sc;                                           // distance / rho OR 2 sqrt(nu) distance/rho
    aux = con*pow(dscale,par[1])*patched_bessel_k(dscale,par[1],1.0);     // Matérn fonction, with nu = par[1]
    cor[0] = ratio * aux;                                         // multiplies by (1-nugget)
  }
#ifdef NO_R_CONSOLE
  printf("%g\n",*cor);
#endif
 }
}

extern "C" {
  static void matern_mat(double *par, double *dist, longint *n, longint *nug, longint *nusc,
	   double *mat)
{
    longint i, j, np1 = *n + 1;
    double dscale,sc,aux, con, *sdist, ratio = 1.0;
#ifdef NO_R_CONSOLE
    printf("corMatern: computes correlation matrix: rho=%g, nu=%g\n",par[0],par[1]);
#endif
    sdist = dist;
    if (*nug) ratio = par[2];
    if (*nusc) {
      sc=2.0*sqrt(par[1])/par[0];
      // when nu (par[1]) is large this compares to corGaus
      // with range rho (par[0])
	} else {
      sc=1.0/par[0];
      // initial formulation: when nu is large this compares
      // to corGaus with range par[0]*2*sqrt(nu)
    }

    con = pow(2.0,par[1] - 1) * exp(lgamma(par[1]));
    con = 1.0/con;

    for(i = 0; i < *n; i++) {
	mat[i * np1] = 1.0;
	for(j = i+  1; j < *n; j++, sdist++) {
          if (*sdist < std::numeric_limits<double>::epsilon())  {                                  // distance ~ 0 (adjust epsilon value to machine precision??)
	    *(mat + i + j * (*n)) = *(mat + j + i * (*n)) = ratio;      // correlation is one!
	  } else {                                                      // distance > 0
	    dscale = *sdist * sc;                                       // distance / rho OR 2 sqrt(nu) distance/rho
	    aux = con*pow(dscale,par[1])*patched_bessel_k(dscale,par[1],1.0);   // Matérn fonction, with nu = par[1]
	    *(mat + i + j * (*n)) = *(mat + j + i * (*n)) = ratio * aux;// multiplies by (1-nugget)
	  }
	}
    }
#ifdef NO_R_CONSOLE
    printf("done!\n");
#endif
}
}

extern "C" {
static void matern_fact(double *par, double *dist, longint *n, longint *nug, longint *nusc,
	     double *mat,
	     double *logdet)
{
    longint job = 11L, info, i, nsq = *n * (*n), np1 = *n + 1;
    double *work = Calloc(*n, double), *work1 = Calloc(nsq, double);
#ifdef NO_R_CONSOLE
    printf("call to matern_fact rho=%g nu=%g\n",par[0],par[1]);
#endif
#ifndef USING_R
    longint zero = 0L;
#endif
    matern_mat(par, dist, n, nug, nusc, mat);
#ifdef USING_R
    //F77_CALL(chol) (mat, n, n, mat, &info); // R < 3.0.0
    //    dpofa(mat,(*n),(*n));  // cholesky via c++ port of dfopa, chol.cpp and chol.h removed from spaMM >= 1.6.2
    F77_CALL(dpofa) (mat, n, n, &info); // direct call to linpack fn
    // optimization wouldbe required if corMatern class were widely used.
#else
    F77_CALL(chol) (mat, n, work, &zero, &zero, &info); //FR: should no longer work with R >= 3.0.0, but not used
#endif
    for(i = 0; i < *n; i++) {
	work1[i * np1] = 1;
	F77_CALL(dtrsl) (mat, n, n, work1 + i * (*n), &job, &info); // linpack
	*logdet -= log(fabs(mat[i * np1]));
    }
    Memcpy(mat, work1, nsq);
    Free(work); Free(work1);
}
}

extern "C" {
void
matern_matList(double *par, longint *nug, longint *nusc, double *dist, longint *pdims,
		double *minD, double *mat)
{
    longint i, M = pdims[1], *len = pdims + 4,
	*start = len + M;
    double aux;
#ifdef NO_R_CONSOLE
    printf("call to matern_matList\n");
#endif
    /* parameter assumed in unconstrained form */
    par[0] = exp(par[0]); //rho
    par[1] = exp(par[1]); //nu
    if (*nug == 1) {
	aux = exp(par[2]);
	par[2] = 1 / (1.0 + aux);	/* 1 - nugget */
    }

    for(i = 0; i < M;  i++) {
      matern_mat(par, dist + start[i], &len[i], nug, nusc, mat);
      mat += len[i] * len[i];
    }
}
}

extern "C" {
  void matern_factList(double *par, longint *nug, longint *nusc, double *dist, longint *pdims,
		 double *minD, double *FactorL, double *logdet)
  {
    longint i, M = pdims[1], *len = pdims + 4,
	*start = len + M;
    double aux;
#ifdef NO_R_CONSOLE
    printf("call to matern_factList rho= %g, nu=%g\n",par[0],par[1]);
#endif
    /* parameter assumed in unconstrained form */
    par[0] = exp(par[0]);
    par[1] = exp(par[1]);
    if (*nug == 1) {
	aux = exp(par[2]);
	par[2] = 1 / (1.0 + aux);	/* 1 - nugget */
    }
    for(i = 0; i < M;  i++) {
      matern_fact(par, dist + start[i], &len[i], nug, nusc, FactorL, logdet);
      FactorL += len[i] * len[i];
    }
  }
}

void
d_axpy(double *y, double a, double *x, longint n)
{				/* y <- a * x + y  */
  while (n-- > 0) { *y++ += a * *x++; }
}


double *
copy_mat(double *y, longint ldy, double *x, longint ldx,
	 longint nrow, longint ncol)
{				/* y <- x */
  double * ret = y;
  while (ncol-- > 0) { Memcpy(y, x, nrow); y += ldy; x += ldx; }
  return ret;
}

double *
mult_mat(double *z, longint ldz,
	 double *x, longint ldx, longint xrows, longint xcols,
	 double *y, longint ldy, longint ycols)
{				/* z <- x %*% y */
  double *t, *tmp = Calloc((size_t)(xrows * ycols), double);
  int i, j;			/* use tmp so z can be either x or y */

  t = tmp;
  for (i = 0; i < ycols; i++) {
    for (j = 0; j < xcols; j++) {
      d_axpy(t, y[j], x + j * ldx, xrows);
    }
    t += xrows;
    y += ldy;
  }
  copy_mat(z, ldz, tmp, xrows, xrows, ycols);
  Free(tmp);
  return z;
}


extern "C" {
void matern_recalc(double *Xy, longint *pdims, longint *ZXcol, double *par,
	       double *dist, double *minD, longint *nug, longint *nusc, double *logdet)
{
    longint N = pdims[0], M = pdims[1],
	*len = pdims + 4, *start = len + M, i;
    double aux, *sXy;
#ifdef NO_R_CONSOLE
    printf("call to matern_recalc rho=%g, nu=%g\n",par[0],par[1]);
#endif
    /* parameter assumed in unconstrained form */
    par[0] = exp(par[0]);
    par[1] = exp(par[1]);
    if (*nug == 1) {
	aux = exp(par[2]);
	par[2] = 1 / (1.0 + aux); /* 1 - nugget */
    }
    for(i = 0, sXy = Xy; i < M;  i++) {
	double *Factor = Calloc(len[i] * len[i], double);
	matern_fact(par, dist + start[i], &len[i], nug, nusc, Factor, logdet);
	mult_mat(sXy, N, Factor, len[i], len[i], len[i], sXy, N, *ZXcol);
	sXy += len[i];
	Free(Factor);
    }
  }
}


