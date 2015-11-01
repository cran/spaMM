/*
 This file is part of the spaMM package for R, distributed under 
 the terms of the Cecill-2 licence. 
 
 It contains routines dealing with estimation of truncated multivariate normal distribution.
 
 The code is derived from bayesmc.c in version 2.2-5 of the package bayesm for R, 
by  Peter Rossi,Anderson School, UCLA, perossichi@gmail.com;
itself apparently derived from code by r mcculloch 8/04;
and distributed under the GPL-2 licence.

It was modified by F. Rousset in particular to return a log probability, and a standard error
*/


#include <R.h>
#include <Rmath.h>
#include <math.h>


/*
 see bayesm's ghkvec_rcpp.cpp for new version using Halton sequences.
 
 if above=1, then we truncate component i from above at point trunpt[i-1]
        L is lower triangular root of Sigma
	random vector is assumed to have zero mean
    	n is number of draws to use in GHK	
*/

extern "C" { 
 void GHK_oneside(double *L, double* trunpt, int *above, int *dim, int *nrep, double *res, double *se) {
   int i,j,k;
   double mu,tpz,u,QMClog,maxlog=0.0,relL,pa,pb,arg;
   double *z;
   z = (double *)R_alloc(*dim,sizeof(double));
   double *QMClogtable;
   QMClogtable = (double *)R_alloc(*nrep,sizeof(double));
   GetRNGstate();
   for(i=0;i<*nrep;i++) { // loop over random draws
      QMClog=0.0;
      for(j=0;j<*dim;j++) { // loop over dimensions of MVN
        mu=0.0; for(k=0;k<j;k++) mu += L[k*(*dim)+j]*z[k];
        tpz = (trunpt[j]-mu)/L[j*(*dim)+j];
        if(above[j]) {
            pa=0.0; pb = pnorm(tpz,0.0,1.0,1,0);   
        } else {
            pb=1.0; pa = pnorm(tpz,0.0,1.0,1,0);
        }
        QMClog += log(pb-pa);
        u = unif_rand();
        arg=u*pb+(1.-u)*pa;
        if(arg > .999999999) arg=.999999999;
        if(arg < .0000000001) arg=.0000000001;
        z[j] = qnorm(arg,0.0,1.0,1,0);
      } // end loop over dimensions
      QMClogtable[i] = QMClog;
      if (i==0 || QMClog > maxlog) maxlog=QMClog; // computes normalizing factor that will prevent underflow.
   }
   *res = 0.0;
   *se = 0.0;
   /** working on relative values of likelihood 
   In particular, we will compute a SE for logL, taking into account that logL is not a mean of logLi values
   => we compute the SE of relL which is a mean of relLi values, then
   SE(logL)~ SE[L]/L = SE(relL)/relL
   */
    for(i=0;i<*nrep;i++) {
      relL = exp(QMClogtable[i]-maxlog);
      *res += relL; 
      *se += relL*relL; 
    } 
   *res /=  *nrep; // mean relLi
   *se /=  *nrep; // mean relLi^2
   *se -= *res * (*res); // mean relLi-squared - squared(mean relLi)
   *se /=  ((*nrep)-1); // var(relL)
   *se = sqrt(*se); // SE(relL)
   // back to logL
   *se /= *res; // SE(logL) *approximated* by SE[L]/L = SE(relL)/relL
   *res = log(*res); // log(relL)
   *res += maxlog;  // logL
   PutRNGstate();
}
}
