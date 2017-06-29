// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// Response class
class Response {
private:
  NumericMatrix parameters; // vector of Response parameters

public:

  // constructor of AR class
  Response(NumericMatrix);


  // getters
  NumericMatrix get_parameters()  const;

  // setters
  void set_parameters(NumericMatrix);


  // function which generates response(s):
  //  inputs:
  //    int -- a size of the cohort of patients for which responses are generated
  virtual double response(int) = 0;
};



// Binary Response
class BinaryResponse: public Response {
public:
  BinaryResponse(NumericMatrix);
  double response(int);
};



// Normally Distributed Response
class NormalResponse: public Response {
public:
  NormalResponse(NumericMatrix);
  double response(int);
};



// Weibully Distributed Response
class WeibullResponse: public Response {
public:
  WeibullResponse(NumericMatrix);
  double response(int);
};
