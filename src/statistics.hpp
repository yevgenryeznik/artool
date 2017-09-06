// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


// sample an integer from a set of integers, given probabilities
int sample(IntegerVector, NumericVector);
  

// t test
int t_test(NumericVector, NumericVector, double);


// anova test
int anova_test(IntegerVector, NumericVector, int, double);
