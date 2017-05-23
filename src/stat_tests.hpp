// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// t test
int t_test(NumericVector, NumericVector, double);

// anova test
int anova_test(NumericMatrix, int, double);
