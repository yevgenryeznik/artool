// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <functional>

using namespace Rcpp;
using namespace std;

// bisection method to find a root of a function over interval
double bisection(std::function<double (double)>, double, double, double);

// Newthon method to find a root of a function over interval
double newthon(std::function<double (double)>, std::function<double (double)>, double, double);

// Secant method to find a root of a function over interval
double secant(std::function<double (double)>, double, double, double);