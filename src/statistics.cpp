#include "statistics.hpp"


// sample an integer from a set of integers, given probabilities
int sample(IntegerVector range, NumericVector prob) {
  NumericVector cumulative_prob(prob.size()+1);
  double u = runif(1)[0];
  int k;
  for (k = 1; k <= prob.size(); k++) {
    cumulative_prob[k] = cumulative_prob[k-1] + prob[k-1];
    if (cumulative_prob[k-1] < u && u < cumulative_prob[k]) {
      break;
    }
  }
  return range[k-1];
}


// t test
//[[Rcpp::export(.t_test)]]
int t_test(IntegerVector treatment, NumericVector response, double alpha = 0.05){
  NumericVector x = response[treatment == 1];
  NumericVector y = response[treatment == 2];
  
  if (x.length() < 2 or y.length() < 2){
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


//[[Rcpp::export(.anova_test)]]
int anova_test(IntegerVector treatment, NumericVector response, int K, double alpha = 0.05){
  NumericVector n(K);
  NumericVector SS(K);
  NumericVector Y_mean(K);
  IntegerVector k = seq_len(K);
  
  NumericVector Y(response);
  
  for_each(k.begin(), k.end(), [treatment, Y, &n, &SS, &Y_mean](int &k){
    NumericVector Yk = Y[treatment == k];
    Y_mean[k-1] = mean(Yk);
    SS[k-1] = sum(Rcpp::pow(Yk-Y_mean[k-1], 2));
    n[k-1] = Yk.length();
  });
  
  if (is_true(any(n <= 1)) or Y.length() <= K){
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
