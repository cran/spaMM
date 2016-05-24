#include <R.h>
#include <Rmath.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
SEXP Rcpp_pMVN(NumericMatrix &L, 
                      NumericVector &limits, 
                      IntegerVector &ismax, 
                      NumericMatrix &rand_seq, 
                      int nrep) {
  double mu,tpz,u,QMClog,maxlog=0.0,relL,pa,pb,arg;
  int dimL=L.nrow();
  vector<double> z(dimL);
  bool use_rand_seq=(rand_seq.nrow()>0);
  if (use_rand_seq) {
    if (rand_seq.ncol()!=dimL) {
      Rcerr <<"Invalid 'rand_seq' does not match the dimension of 'L'";
      return R_NilValue;
    } else {
      nrep=rand_seq.nrow();
    }
  }
  vector<double> QMClogtable(nrep);
  GetRNGstate();
  for(int i=0;i< nrep;i++) { // loop over random draws
    QMClog=0.0;
    //Rcout <<  L[0,0] <<" "<< L[1,0] <<" " << L[1,1] << endl;
    for(int lig=0;lig<dimL;lig++) { // loop over dimensions of MVN
      mu=0.;
      for(int col=0;col<lig;col++) {
        mu += L(lig,col)*z[col]; // (,) not [,] which seems to use R-style indices (from 1 !)
        //Rcout << lig <<" "<<col <<" " << L(lig,col) << endl;
      }
      tpz = (limits[lig]-mu)/L(lig,lig);
      if(ismax[lig]) {
        pa=0.0; pb = R::pnorm(tpz,0.0,1.0,1,0);   
      } else {
        pb=1.0; pa = R::pnorm(tpz,0.0,1.0,1,0);
      }
      QMClog += log(pb-pa);
      if (use_rand_seq) {
        u = rand_seq(i,lig);
      } else u = unif_rand(); //GHK
      arg=u*pb+(1.-u)*pa;
      if(arg > .999999999) arg=.999999999;
      if(arg < .0000000001) arg=.0000000001;
      z[lig] = R::qnorm(arg,0.0,1.0,1,0);
    } // end loop over dimensions
    QMClogtable[i] = QMClog;
    if (i==0 || QMClog > maxlog) maxlog=QMClog; // computes normalizing factor that will prevent underflow.
  }
  PutRNGstate();  
  double res = 0.0;
  double se = 0.0;
  /** working on relative values of likelihood 
  In particular, we will compute a SE for logL, taking into account that logL is not a mean of logLi values
  => we compute the SE of relL which is a mean of relLi values, then
  SE(logL)~ SE[L]/L = SE(relL)/relL
  */
  const int nblocks=25;
  if (use_rand_seq && nrep> 2*nblocks) { // suboptimal estimator... 
    int blocksize = nrep / nblocks;
    double bres, sumrelL=0;
    for (int b=0; b < nblocks; b++) {
      bres=0;
      for(int i=b*blocksize;i<(b+1)*blocksize;i++) {
        relL = exp(QMClogtable[i]-maxlog);
        bres += relL; 
      }
      sumrelL += bres; // sum relL
      bres = bres/blocksize; // block mean relL
      res += bres; // sum  block means
      se += bres*bres; // sumsquared block means
    }
    res /=  nblocks; // mean block mean relL
    se /=  nblocks; // mean block mean^2 relL
    se -= res * (res); // 
    se /=  (nblocks-1); // var(block mean)
    se = sqrt(se); // SE(block mean relL)
    for(int i=nblocks*blocksize;i<nrep;i++) {
      relL = exp(QMClogtable[i]-maxlog);
      sumrelL += relL;
    }
    res = sumrelL/nrep; // mean relLi
  } else {
    for(int i=0;i<nrep;i++) {
      relL = exp(QMClogtable[i]-maxlog);
      res += relL; 
      se += relL*relL; 
    } 
    res /=  nrep; // mean relLi
    se /=  nrep; // mean relLi^2
    se -= res * (res); // mean relLi-squared - squared(mean relLi)
    se /=  ((nrep)-1); // var(relL)
    se = sqrt(se); // SE(relL)
  }
  // back to logL
  se /= res; // SE(logL) *approximated* by SE[L]/L = SE(relL)/relL
  res = log(res); // log(relL)
  res += maxlog;  // logL
  return(List::create(_["logInt"]=res, _["seInt"]=se));
}


