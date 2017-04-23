#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>

NumericVector isPositive(NumericVector x){
  int n = x.size();
  NumericVector out(n);
  int i;
  for(i=0; i<n; i++){
    out[i] = x[i] >= 0;
  }
  return out;
}


//****************************************************************************80
// [[Rcpp::export]]
NumericVector RcppOwenQ1(int nu, double t, NumericVector delta, NumericVector R, Function Tha){
  double a = t/sqrt(nu);
  double b = nu/(nu+t*t);
  double sB = sqrt(b);
  double asB = R::sign(t)*sqrt(t*t/(nu+t*t));
  int J = delta.size();
  if(nu==1){
    NumericVector C = pnorm(R) - isPositive(delta);
    int i;
    for(i=0; i<J; i++){
      NumericVector C1 = Tha(delta[i]*sB, a);
      NumericVector C2 = Tha(R[i], (a*R[i]-delta[i])/R[i]);
      NumericVector C3 = Tha(delta[i]*sB, (delta[i]*a*b-R[i])/b/delta[i]);
      C[i] += 2*(C1[0] - C2[0] - C3[0]);
    }
    return C;
  }
  int n = nu-1;
  NumericMatrix H(n,J);
  NumericMatrix M(n,J);
  H(0,_) = -dnorm(R) * pnorm(a*R-delta);
  M(0,_) = asB*dnorm(delta*sB)*(pnorm(delta*asB)-pnorm((delta*a*b-R)/sB));
  if(nu >= 3){
    H(1,_) = R * H(0,_);
    M(1,_) = b*(delta*a*M(0,_) + a*dnorm(delta*sB)*(dnorm(delta*asB)-dnorm((delta*a*b-R)/sB)));
    if(nu >= 4){
      NumericVector A(n);
      NumericMatrix L(n-2,J);
      A[0] = 1;
      A[1] = 1;
      L(0,_) = a * b * R * dnorm(R) * dnorm(a*R-delta) / 2;
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          L(k,_) = A[k+2] * R * L(k-1,_);
        }
      }
      for(k=2; k<n; k++){
        H(k,_) = A[k] * R * H(k-1,_);
        M(k,_) = (k-1.0)/k * b * (A[k-2] * delta * a * M(k-1,_) + M(k-2,_)) - L(k-2,_);
      }
    }
  }
  if(nu % 2 == 0){
    double sqrt2pi = 2.506628274631000502415765284811;
    NumericVector sum(J);
    int i;
    for(i=0; i<nu-1; i+=2){
      sum += M(i,_)+H(i,_);
    }
    return pnorm(-delta) + sqrt2pi * sum;
  }else{
    NumericVector sum(J);
    int i;
    for(i=1; i<nu-1; i+=2){
      sum += M(i,_)+H(i,_);
    }
    NumericVector C = pnorm(R) - isPositive(delta);
    for(i=0; i<J; i++){
      NumericVector C1 = Tha(delta[i]*sB, a);
      NumericVector C2 = Tha(R[i], (a*R[i]-delta[i])/R[i]);
      NumericVector C3 = Tha(delta[i]*sB, (delta[i]*a*b-R[i])/b/delta[i]);
      C[i] += 2*(C1[0] - C2[0] - C3[0]);
    }
    return C+2*sum;
  }
}
