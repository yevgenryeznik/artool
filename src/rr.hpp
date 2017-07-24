// A file contains description of classes providing Restrict Randomization procedures 
// targeting unequal allocation

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// Restricted Randomization Procedure  (base class)
class RR {

  NumericVector parameters;              // vector of RR procedure parameters
  IntegerVector fixed_allocation_ratio;  // vector contains fixed allocation ratio
  NumericVector target_allocation;       // vector contains fixed allocation ratio
  int number_of_treatments;              // number of treatment groups
  int number_of_subjects;                // number of subjects to randomize
  IntegerVector treatment;               // vector of treatment assignments  
  NumericMatrix rand_probability;        // matrix of randomization probabilities
  
  // operational characteristics
  NumericMatrix alloc_proportion;        // matrix of allocation proportions
  NumericVector forcing_idx;             // vector of forcing indicies
  NumericVector imbalance;               // vector of imbalances
  NumericVector selection_bias;          // vector of selection biases
  
public:

  // constructor of RR class
  RR(NumericVector, IntegerVector, int);


  // getters
  NumericVector get_parameters()              const;
  IntegerVector get_fixed_allocation_ratio()  const;
  NumericVector get_target_allocation()       const;
  int get_number_of_treatments()              const;
  int get_number_of_subjects()                const;
  IntegerVector get_treatment()               const;
  NumericMatrix get_rand_probability()        const;
  NumericMatrix get_alloc_proportion()        const;
  NumericVector get_forcing_idx()             const;
  NumericVector get_imbalance()               const;
  NumericVector get_selection_bias()          const;
  

  // setters
  void set_treatment(int, int);
  void set_rand_probability(int, NumericVector);
  void set_alloc_proportion(int, NumericVector);
  void set_forcing_idx(int, double);
  void set_imbalance(int, double);
  void set_selection_bias(int, double);
  
  
  // sample an integer from a set of integers given probabilities
  int sample(IntegerVector, NumericVector);
  
  // function which randomizes subjects across treatments
  void randomize();

  // function which adapts allocation probabilities
  virtual void adapt(int, IntegerVector) = 0;
};



// Below there are classes of Restricted Randomization (RR) Procedures targeting unequal allocation
// Completely Randomized Design
class CRD: public RR {
public:
  CRD(NumericVector, IntegerVector, int);
  
  // function which adapts allocation probabilities
  void adapt(int, IntegerVector) override;
};



// Permuted Block Design
class PBD: public RR {
public:
  PBD(NumericVector, IntegerVector, int);

  // function which adapts allocation probabilities
  void adapt(int, IntegerVector) override;
};



// Block Urn Design
class BUD: public RR {
public:
  BUD(NumericVector, IntegerVector, int);

  // function which adapts allocation probabilities
  void adapt(int, IntegerVector) override;
};



// Mass Weighted Urn Design
class MWUD: public RR {
public:
  MWUD(NumericVector, IntegerVector, int);

  // function which adapts allocation probabilities
  void adapt(int, IntegerVector) override;
};



// Drop-the-Loser rule
class DL: public RR {
  IntegerVector urn;
public:

  DL(NumericVector, IntegerVector, int);

  // getter
  IntegerVector get_urn() const;

  // setter
  void set_urn(IntegerVector);

  // function which adapts allocation probabilities
  void adapt(int, IntegerVector) override;
};



// Doubly Adaptive Biased Coin Design
class DBCD: public RR {
public:

  DBCD(NumericVector, IntegerVector, int);

  // function which adapts allocation probabilities
  void adapt(int, IntegerVector) override;
};



// Minimum Quadratic Distance constrained balance optimization
class MinQD: public RR {
public:

  MinQD(NumericVector, IntegerVector, int);

  // function which adapts allocation probabilities
  void adapt(int, IntegerVector) override;
};



// Maximum Entropy constrained balance optimization
class MaxEnt: public RR {
public:

  MaxEnt(NumericVector, IntegerVector, int);

  // function which adapts allocation probabilities
  void adapt(int, IntegerVector) override;
};
// ===== END: Restricted Randomization (RR) Procedures targeting unequal allocation =====




