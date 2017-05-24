// A file contains description of classes providing Restrict Randomization procedures 
// targeting unequal allocation

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// Restricted Randomization Procedure  (base class)
class RR {
  NumericVector parameters;             // vector of RR procedure parameters
  IntegerVector fixed_allocation_ratio; // vector contains fixed allocation ratio given treatment assignments 
  NumericVector rand_probability;       // vector of randomization probabilities
public:

  // constructor of RR class
  RR(NumericVector, IntegerVector);


  // getters
  NumericVector get_parameters()              const;
  IntegerVector get_fixed_allocation_ratio()  const;
  NumericVector get_rand_probability()        const;
  
  // setters
  void set_parameters(NumericVector);
  void set_rand_probability(NumericVector);
  
  // function which adapts randomization probabilities:
  //  input:
  //    IntegerVector -- vector of an allocation ratio
  virtual void adapt(IntegerVector) =  0;
};



// Below there are classes of Restricted Randomization (RR) Procedures targeting unequal allocation
// Completely Randomized Design
class CRD: public RR {
public:
  CRD(NumericVector, IntegerVector);
  
  // function which adapts randomization probabilities:
  void adapt(IntegerVector) override;
};



// Permuted Block Design
class PBD: public RR {
public:
  PBD(NumericVector, IntegerVector);

  // function which adapts randomization probabilities:
  void adapt(IntegerVector) override;
};



// Block Urn Design
class BUD: public RR {
public:
  BUD(NumericVector, IntegerVector);

  // function which adapts randomization probabilities:
  void adapt(IntegerVector) override;
};



// Mass Weighted Urn Design
class MWUD: public RR {
public:
  MWUD(NumericVector, IntegerVector);

  // function which adapts randomization probabilities:
  void adapt(IntegerVector) override;
};



// Drop-the-Loser rule
class DL: public RR {
  IntegerVector urn;
public:

  DL(NumericVector, IntegerVector);

  // getter
  IntegerVector get_urn() const;

  // setter
  void set_urn(IntegerVector);

  // function to sample a ball from an urn given discrete probability distribution
  int sample_ball();

  // function which adapts randomization probabilities:
  void adapt(IntegerVector) override;
};



// Doubly Adaptive Biased Coin Design
class DBCD: public RR {
public:

  DBCD(NumericVector, IntegerVector);

  // function which adapts randomization probabilities:
  void adapt(IntegerVector) override;
};



// Minimum Quadratic Distance constrained balance optimization
class MinQD: public RR {
public:

  MinQD(NumericVector, IntegerVector);

  // function which adapts randomization probabilities:
  void adapt(IntegerVector) override;
};



// Maximum Entropy constrained balance optimization
class MaxEnt: public RR {
public:

  MaxEnt(NumericVector, IntegerVector);

  // function which adapts randomization probabilities:
  void adapt(IntegerVector) override;
};
// ===== END: Restricted Randomization (RR) Procedures targeting unequal allocation =====




