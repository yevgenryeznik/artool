#include "stat_tests.hpp"

// t test
int t_test(NumericVector x, NumericVector y, double alpha = 0.05){
  if (x.length() < 2 | y.length() < 2){
    return R_NaReal;
  }
  else {
    int nx = x.length();
    double x_mean = mean(x);
    double sx = sd(x);

    int ny = y.length();
    double y_mean = mean(y);
    double sy = sd(y);

    double nu = (sx*sx/nx+sy*sy/ny)*(sx*sx/nx+sy*sy/ny)/
      ((sx*sx/nx)*(sx*sx/nx)/(nx-1) + (sy*sy/ny)*(sy*sy/ny)/(ny-1));
    double Tn = (x_mean - y_mean)/sqrt(sx*sx/nx + sy*sy/ny);

    bool reject = (std::abs((float)Tn) > Rcpp::qt(NumericVector::create(1-alpha/2), nu))[0];
    return (int)reject;
  }
}


int anova_test(NumericMatrix obs, int K, double alpha = 0.05){
  NumericVector n(K);
  NumericVector SS(K);
  NumericVector Y_mean(K);
  NumericVector Y = obs.column(1);   // responses
  NumericVector trt = obs.column(0); // tretament ids
  for (int k = 1; k <= K; k++){
    Y_mean[k-1] = mean(as<NumericVector>(Y[trt==k]));
    SS[k-1] = sum(Rcpp::pow(as<NumericVector>(Y[trt == k])-Y_mean[k-1], 2));
    n[k-1] = as<NumericVector>(Y[trt == k]).length();
  }
  if (is_true(any(n < K))){
    return R_NaInt;
  }
  else {
    double MSw = sum(SS)/(Y.length()-K);
    double MSb = sum(n*Rcpp::pow(Y_mean - mean(Y), 2))/(K-1);

    double F = MSb/MSw;

    bool reject = (F > Rcpp::qf(NumericVector::create(1-alpha), K-1, Y.length()-K))[0];
    return (int)reject;
  }
}
