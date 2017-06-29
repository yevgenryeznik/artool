#include "response.hpp"


// Response base class
Response::Response(NumericMatrix parameters_) {
  parameters = parameters_;
}

// getters of Response class
NumericMatrix Response::get_parameters() const { return(parameters); }

// setters of Response class
void Response::set_parameters(NumericMatrix parameters_) { parameters = parameters_; }



// ===== Binary Response =====
BinaryResponse::BinaryResponse(NumericMatrix parameters_):Response(parameters_) {
  if (parameters_.ncol() != 1 ) {
    throw invalid_argument("Constructor of a BinaryResponse class must contain 1-column matrix as parameter: (p)!");
  }
}

double BinaryResponse::response(int treatment) {
  NumericMatrix parameters = get_parameters();
  return rbinom(1, 1, parameters(treatment-1, 0))[0];
}
// =========================================



// ===== Normally Distributed Response =====
NormalResponse::NormalResponse(NumericMatrix parameters_):Response(parameters_) {
  if (parameters_.ncol() != 2 ) {
    throw invalid_argument("Constructor of a NormalResponse class must contain 2-column matrix as parameters: (mean, sd)!");
  }
}

double NormalResponse::response(int treatment) {
  NumericVector parameters = get_parameters();
  return rnorm(1, parameters(treatment-1, 0), parameters(treatment-1, 1))[0];
}
// =========================================



// ===== Weibully Distributed Response =====
WeibullResponse::WeibullResponse(NumericMatrix parameters_):Response(parameters_) {
  if (parameters_.ncol() != 2 ) {
    throw invalid_argument("Constructor of a WeibullResponse class must contain 2-column matrix as parameters: (shape, scale)!");
  }
}

double WeibullResponse::response(int treatment) {
  NumericMatrix parameters = get_parameters();
  return rweibull(1, parameters(treatment-1, 0), parameters(treatment-1, 1))[0];
}
// =========================================


