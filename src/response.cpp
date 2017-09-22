// file contains implementations of functions simulating different responses
// given parametric diftribution and treatment assignments

#include "response.hpp"


// binary response
NumericVector response_binary(IntegerVector treatment, List parameter) {
  if (parameter.size() != 1) {
    throw invalid_argument("Input (parameter) of response_binary function must contain 1-item list of parameter values: (prob)!");
  }
  IntegerVector k = seq_len(treatment.size());
  NumericVector prob = as<NumericVector>(parameter["prob"]);
  NumericVector response(treatment.size());
  
  for_each(k.begin(), k.end(), [&response, treatment, prob](int &k){
    response[k-1] = rbinom(1, 1, prob[treatment[k-1]-1])[0];
  });
  
  return response;
}


// uniform response
NumericVector response_uniform(IntegerVector treatment, List parameter) {
  if (parameter.size() != 2) {
    throw invalid_argument("Input (parameter) of response_uniform function must contain 2-item list of parameters values: (min, max)!");
  }
  IntegerVector k = seq_len(treatment.size());
  NumericVector min_ = as<NumericVector>(parameter["min"]);
  NumericVector max_ = as<NumericVector>(parameter["max"]);
  NumericVector response(treatment.size());
  
  for_each(k.begin(), k.end(), [&response, treatment, min_, max_](int &k){
    response[k-1] = runif(1, min_[treatment[k-1]-1], max_[treatment[k-1]-1])[0];
  });
  
  return response;
}


// normal response
NumericVector response_normal(IntegerVector treatment, List parameter) {
  if (parameter.size() != 2) {
    throw invalid_argument("Input (parameter) of response_normal function must contain 2-item list of parameters values: (mean, sd)!");
  }
  IntegerVector k = seq_len(treatment.size());
  NumericVector mean_ = as<NumericVector>(parameter["mean"]);
  NumericVector sd_ = as<NumericVector>(parameter["sd"]);
  NumericVector response(treatment.size());
  
  for_each(k.begin(), k.end(), [&response, treatment, mean_, sd_](int &k){
    response[k-1] = rnorm(1, mean_[treatment[k-1]-1], sd_[treatment[k-1]-1])[0];
  });
  
  return response;
}


// exponential response
NumericVector response_exp(IntegerVector treatment, List parameter) {
  if (parameter.size() != 1) {
    throw invalid_argument("Input (parameter) of response_exp function must contain 1-item list of parameter values: (rate)!");
  }
  IntegerVector k = seq_len(treatment.size());
  NumericVector rate = as<NumericVector>(parameter["rate"]);
  NumericVector response(treatment.size());
  
  for_each(k.begin(), k.end(), [&response, treatment, rate](int &k){
    response[k-1] = rexp(1, rate[treatment[k-1]-1])[0];
  });
  
  return response;
}


// weibull response
NumericVector response_weibull(IntegerVector treatment, List parameter) {
  if (parameter.size() != 2) {
    throw invalid_argument("Input (parameter) of response_weibull function must contain 2-item list of parameters values: (shape, scale)!");
  }
  IntegerVector k = seq_len(treatment.size());
  NumericVector shape = as<NumericVector>(parameter["shape"]);
  NumericVector scale = as<NumericVector>(parameter["scale"]);
  NumericVector response(treatment.size());

  for_each(k.begin(), k.end(), [&response, treatment, shape, scale](int &k){
    double u = runif(1)[0];
    double w = log(-log(1-u));
    response[k-1] = scale[treatment[k-1]-1]*exp(w/shape[treatment[k-1]-1]);
  });
 
  return response;
}



// log-logistic response
NumericVector response_loglog(IntegerVector treatment, List parameter) {
  if (parameter.size() != 2) {
    throw invalid_argument("Input (parameter) of response_weibull function must contain 2-item list of parameters values: (shape, scale)!");
  }
  IntegerVector k = seq_len(treatment.size());
  NumericVector shape = as<NumericVector>(parameter["shape"]);
  NumericVector scale = as<NumericVector>(parameter["scale"]);
  NumericVector response(treatment.size());
  
  for_each(k.begin(), k.end(), [&response, treatment, shape, scale](int &k){
    double u = runif(1)[0];
    double w = log(u/(1-u));
    response[k-1] = scale[treatment[k-1]-1]*exp(w/shape[treatment[k-1]-1]);
  });
  
  return response;
}


// log-logistic response
NumericVector response_lognormal(IntegerVector treatment, List parameter) {
  if (parameter.size() != 2) {
    throw invalid_argument("Input (parameter) of response_weibull function must contain 2-item list of parameters values: (mu, sigma)!");
  }
  IntegerVector k = seq_len(treatment.size());
  NumericVector mu_ = as<NumericVector>(parameter["mu"]);
  NumericVector sigma_ = as<NumericVector>(parameter["sigma"]);
  NumericVector response(treatment.size());
  
  for_each(k.begin(), k.end(), [&response, treatment, mu_, sigma_](int &k){
    double w = qnorm(runif(1))[0];
    response[k-1] = exp(mu_[treatment[k-1]-1]+w*sigma_[treatment[k-1]-1]);
  });
  
  return response;
}


// set response function given 
//    distribution -- response distribution
//    parameter -- List with response distribution parameter(s) values
std::function<NumericVector (IntegerVector)> set_response_function(std::string distribution, List parameter) {
  std::function<NumericVector (IntegerVector)> fcn;
  if (distribution == "binary") { // Binary response
    fcn = [parameter](IntegerVector treatment){
      return response_binary(treatment, parameter);
    };
  }
  else if (distribution == "uniform") { // Uniform response
    fcn = [parameter](IntegerVector treatment){
      return response_uniform(treatment, parameter);
    };
  }
  else if (distribution == "normal") { // Normal response
    fcn = [parameter](IntegerVector treatment){
      return response_normal(treatment, parameter);
    };
  }
  else if (distribution == "exponential") { // Exponential response
    fcn = [parameter](IntegerVector treatment){
      return response_exp(treatment, parameter);
    };
  } 
  else if (distribution == "weibull") { // Weibull response
    fcn = [parameter](IntegerVector treatment){
      return response_weibull(treatment, parameter);
    };
  } 
  else if (distribution == "loglogistic") { // Log-logistic response
    fcn = [parameter](IntegerVector treatment){
      return response_loglog(treatment, parameter);
    };
  } 
  else if (distribution == "lognormal") { // Log-normal response
    fcn = [parameter](IntegerVector treatment){
      return response_lognormal(treatment, parameter);
    };
  } 
  else {
    throw invalid_argument("Inappropriate name of a response distribution");
  }
  return fcn;
}


// response function, given
//  -- distribution name
//  -- distribution params
//  -- vector of treatment assignments
// [[Rcpp::export(.response)]]
NumericVector response(std::string distribution, List parameter, IntegerVector treatment) {
  // response function
  std::function<NumericVector (IntegerVector)> response_function = 
    set_response_function(distribution, parameter);
  
  return response_function(treatment);
}

