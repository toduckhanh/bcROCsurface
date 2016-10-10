/*
 Estimation Volume Under ROC Surface.
 Author: Khanh To Duc
 */

#include <Rcpp.h>
using namespace Rcpp;

static inline double indvus(double a, double b, double c) {
  if((a < b) && (b < c)) {
    return 1.0;
  }else if((a < b) && (b == c)){
    return 0.5;
  }else if((a == b) && (b < c)){
    return 0.5;
  }else if((a == b) && (b == c)){
    return 1.0/6;
  }
  else{
    return 0.0;
  }
}

// [[Rcpp::export]]
double vusC(NumericVector tt, NumericMatrix dd){
  int nn = tt.size();
  double I_ijk = 0.0;
  double den = 0.0, num = 0.0, temp = 0.0;
  for(int i = 0; i < nn; i++){
    for(int j = 0; j < nn; j++){
      if(j != i){
        for(int k = 0; k < nn; k++){
          if((k != j) && (k != i)){
            I_ijk = indvus(tt[i], tt[j], tt[k]);
            temp = dd(i, 0)*dd(j, 1)*dd(k, 2);
            num += I_ijk*temp;
            den += temp;
          }
        }
      }
    }
  }
  double out = num/den;
  return out;
}
