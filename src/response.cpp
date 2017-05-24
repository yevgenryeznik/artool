#include "response.hpp"


// Response base class
Response::Response(NumericVector parameters_) {
  parameters = parameters_;
}

// getters of Response class
NumericVector Response::get_parameters() const { return(parameters); }

// setters of Response class
void Response::set_parameters(NumericVector parameters_) { parameters = parameters_; }



// ===== Binary Response =====
BinaryResponse::BinaryResponse(NumericVector parameters_):Response(parameters_) {
  if (parameters_.size() != 1 ) {
    throw invalid_argument("Constructor of a BinaryResponse class must contain 1 parameter: (p)!");
  }
}

NumericVector BinaryResponse::response(int cohort_size) {
  NumericVector parameters = get_parameters();
  return rbinom(cohort_size, 1, parameters[0]);
}
// =========================================



// ===== Normally Distributed Response =====
NormalResponse::NormalResponse(NumericVector parameters_):Response(parameters_) {
  if (parameters_.size() != 2 ) {
    throw invalid_argument("Constructor of a NormalResponse class must contain 2 parameters: (mean, sd)!");
  }
}

NumericVector NormalResponse::response(int cohort_size) {
  NumericVector parameters = get_parameters();
  return rnorm(cohort_size, parameters[0], parameters[1]);
}
// =========================================



// ===== Weibully Distributed Response =====
WeibullResponse::WeibullResponse(NumericVector parameters_):Response(parameters_) {
  if (parameters_.size() != 2 ) {
    throw invalid_argument("Constructor of a NormalResponse class must contain 2 parameters: (mean, sd)!");
  }
}

NumericVector WeibullResponse::response(int cohort_size) {
  NumericVector parameters = get_parameters();
  return rweibull(cohort_size, parameters[0], parameters[1]);
}
// =========================================


