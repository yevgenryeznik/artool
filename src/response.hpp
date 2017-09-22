// file contains definitions of functions simulating different responses
// given parametric diftribution and treatment assignments

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// binary response
NumericVector response_binary(IntegerVector, List);


// uniform response
NumericVector response_uniform(IntegerVector, List);


// normal response
NumericVector response_normal(IntegerVector, List);


// exponential response
NumericVector response_exp(IntegerVector, List);


// weibull response
NumericVector response_weibull(IntegerVector, List);


// log-logistic response
NumericVector response_loglog(IntegerVector, List);


// log-normal response
NumericVector response_lognormal(IntegerVector, List);


// set response function given 
//    distribution -- response distribution
//    parameter -- List with response distribution parameter(s) values
std::function<NumericVector (IntegerVector)> set_response_function(std::string, List);


// response function, given
//  -- distribution name
//  -- distribution params
//  -- vector of treatment assignments
NumericVector response(std::string, List, IntegerVector);
  