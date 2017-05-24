#include "algorithms.hpp"


// bisection method to find a root of a function over interval [low, upp]
double bisection(std::function<double (double)> fcn, double low, double upp, double tol) {
  double fcn_low = fcn(low);
  double fcn_upp = fcn(upp);
  
  // midpoint of a search interval
  double mid = 0.5*(low + upp);
  double fcn_mid = fcn(mid);
  
  while (abs(fcn_mid) > tol && abs(upp-low) > tol){
    if (fcn_low*fcn_mid > 0) {
      low = mid;
      fcn_low = fcn(low);
    }
    else if (fcn_low*fcn_mid < 0){
      upp = mid;
      fcn_upp = fcn(upp);
    }
    mid = 0.5*(low + upp);
    fcn_mid = fcn(mid);
  }
  
  return mid;
}
