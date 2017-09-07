// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <functional>

using namespace Rcpp;
using namespace std;

// bisection method to find a root of a function over interval
double bisection(std::function<double (double)>, double, double, double);

