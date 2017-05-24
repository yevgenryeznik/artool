// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// Response class
class Response {
private:
  NumericVector parameters; // vector of Response parameters

public:

  // constructor of AR class
  Response(NumericVector);


  // getters
  NumericVector get_parameters()  const;

  // setters
  void set_parameters(NumericVector);


  // function which generates response(s):
  //  inputs:
  //    int -- a size of the cohort of patients for which responses are generated
  virtual NumericVector response(int) = 0;
};



// Binary Response
class BinaryResponse: public Response {
public:
  BinaryResponse(NumericVector);
  NumericVector response(int);
};



// Normally Distributed Response
class NormalResponse: public Response {
public:
  NormalResponse(NumericVector);
  NumericVector response(int);
};



// Weibully Distributed Response
class WeibullResponse: public Response {
public:
  WeibullResponse(NumericVector);
  NumericVector response(int);
};
