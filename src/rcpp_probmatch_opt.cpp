#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

NumericVector concatOpt(NumericVector ya, NumericVector yb, int m){
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
int mismatchesOpt(NumericVector hh, NumericVector y, int m, bool twoSides=true){
  int mm1 = 0; int mm2 = 0;
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
      mm1++;
    }
    if (m2==0 && twoSides){
      mm2++;
    }
  }

  return std::max(mm1,mm2);
}

// [[Rcpp::export]]
int partmmOpt(NumericVector hh, NumericVector y, int m){
  int mm1 = 0; int mm2 = 0;
  for (int i=0; i<m*2; i++){
    int m1 = 0; int m2 = 0;
    for (int j=0; j<m*2; j++){
      if (j<m && y(j)==hh(i)){
        m1=1;
      }
      if (i<m && y(i)==hh(j)){
        m2=1;
      }
    }
    if (m1==0){
      mm1++;
    }
    if (i<m && m2==0){
      mm2++;
    }
  }
  return std::max(mm1-m,mm2);
}

// [[Rcpp::export]]
double matchlkOpt(NumericVector hh, NumericMatrix y, NumericVector Y, int n, int m, int lmt, bool twoSides=true){
  double p = 0;
  if (n==0){return(0);}
  for(int i=0; i<n; i++){
    if (partmmOpt(hh, y(i,_), m)<=lmt){
      for (int j=i; j<n; j++){
        int k = 1+(i!=j);
        NumericVector ty = concatOpt(y(i,_), y(j,_), m);
        int w = (mismatchesOpt(hh, ty, m, twoSides)<=lmt);
        p += k*Y(i)*Y(j)*w;
      }
    }
  }

  return p;
}

// [[Rcpp::export]]
double donorlkOpt(NumericMatrix x, NumericVector X, NumericMatrix y, NumericVector Y,int nx,
               int ny, int m, int lmt, double avail, int N, bool twoSides=true){
  double P = 0;
  for(int i=0; i<nx; i++){
    for(int j=i; j<nx; j++){
      int k = 1+(i!=j);
      NumericVector tx = concatOpt(x(i,_), x(j,_), m);
      double p = matchlkOpt(tx, y, Y, ny, m, lmt, twoSides);
      P += k*X(i)*X(j)*(1-pow((1-p), (avail*N)));
    }
  }
  return P;
}

