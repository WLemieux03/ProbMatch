#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

NumericVector concat(NumericVector ya, NumericVector yb, int m){
  NumericVector ty(m*2);
  for(int i=0; i<m; i++){
    ty[i] = ya[i];
  }
  for(int i=0; i<m; i++){
    ty[m+i] = yb[i];
  }
  return ty;
}

// [[Rcpp::export]]
int mismatches(NumericVector hh, NumericVector y, int m, bool twoSide=true){
  int mm = 0;
  for (int i=0; i<m*2; i++){
    int m1 = 0; int m2 = 0;
    for (int j=0; j<m*2; j++){
      if (y(j)==hh(i)){
        m1=1;
      }
      if (y(i)==hh(j)){
        m2=1;
      }
    }
    if (m1==0){
      mm++;
    }
    if (m2==0 && twoSide){
      mm++;
    }
  }

  return mm;
}

// [[Rcpp::export]]
double matchlk(NumericVector hh, NumericMatrix y, NumericVector Y, int n, int m, int lmt){
  double p = 0;
  if (n==0){return(0);}
  for(int i=0; i<n; i++){
    for (int j=i; j<n; j++){
      int k = 1+(i!=j);
      NumericVector ty = concat(y(i,_), y(j,_), m);
      int w = (mismatches(hh, ty, m)<=lmt);
      p += k*Y(i)*Y(j)*w;
    }
  }

  return p;
}

// [[Rcpp::export]]
double donorlk(NumericMatrix x, NumericVector X, NumericMatrix y, NumericVector Y,int nx,
               int ny, int m, int lmt, double avail, int N){
  double P = 0;
  for(int i=0; i<nx; i++){
    for(int j=i; j<nx; j++){
      int k = 1+(i!=j);
      NumericVector tx = concat(x(i,_), x(j,_), m);
      double p = matchlk(tx, y, Y, ny, m, lmt);
      P += k*X(i)*X(j)*(1-pow((1-p), (avail*N)));
    }
  }
  return P;
}
