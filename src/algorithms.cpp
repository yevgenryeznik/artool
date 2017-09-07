#include "algorithms.hpp"


// bisection method to find a root of a function over interval [low, upp]
double bisection(std::function<double (double)> fcn, double low, double upp, double tol) {
  if (abs(fcn(low)) < tol) {
    return low;
  }
  else if (abs(fcn(upp)) < tol) {
    return upp;
  }
  else {
    if (fcn(low)*fcn(upp) > 0) {
      throw std::invalid_argument("Function has the same signs at the interval bounds!");
    }
    else {
      double mid = 0.5*(low + upp);
      int low_cond, upp_cond;
      
      while (abs(fcn(mid)) > tol && (upp-low) > tol) {
        low_cond = (int)(fcn(mid)*fcn(low) > 0);
        low = low_cond*mid + (1-low_cond)*low;
        
        upp_cond = (int)(fcn(mid)*fcn(upp) > 0);
        upp = upp_cond*mid + (1-upp_cond)*upp;
        
        mid = 0.5*(low + upp);
      }
      return mid;
    }
  }
}

