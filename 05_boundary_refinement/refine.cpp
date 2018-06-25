#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//[[Rcpp::export]]
List refine_step1(arma::mat Yt, int min_len = 5) {
  
  int n = Yt.n_cols;
  arma::mat stat1(n, n); stat1.zeros();
  double max_value = - arma::datum::inf;
  int max_l = 0; int max_r = 0;
  
  arma::mat Mcorr = arma::cor(Yt);
  double sumS = arma::accu(arma::trimatl(Mcorr)) - n;
  int n2 = (1+n)*n/2 - n;
  double tmp = 0;
  for (int i=0; i<=(n-1); i++) {
    int k = i+min_len-1;
    for (int j=k; j<=(n-1); j++) {
      
      arma::mat M = Mcorr.submat(i, i, j, j);
      arma::mat x = arma::trimatl(M);
      int xcol = x.n_cols;
      int n1 = (1+xcol)*xcol/2 - xcol;
      double xsum = arma::accu(x) - xcol;
      double xmean = xsum/n1;
      
      if (n1 == n2){
        tmp = xmean/sqrt(1/double(n1));
      } else {
        
        tmp = (xmean - (sumS - xsum)/double(n2-n1))/sqrt(1/double(n1)+1/double(n2-n1));
      }
      stat1(i,j) = tmp;
      
      if (tmp > max_value) {
        max_value = tmp;
        max_l = i+1;
        max_r = j+1;
      }        
    }
    
  }
  
  return List::create(
    _["max.value"] = max_value,
    _["max.l"] = max_l,
    _["max.r"] = max_r 
    );
}



