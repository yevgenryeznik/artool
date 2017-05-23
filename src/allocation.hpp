// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// Allocation base class
class Allocation {
private:
  NumericVector parameter_set;            // vector of parameters of a randomization procedure
  NumericVector fixed_allocation_ratio;   // fixed allocation ratio
  IntegerVector allocation_ratio;         // current allocation ratio
  NumericVector target_allocation;        // target allocation proportions

public:

  // constructor of Allocation class
  Allocation(NumericVector parameter_set, NumericVector fixed_allocation_ratio);


  // getters
  NumericVector get_parameter_set()           const;
  NumericVector get_fixed_allocation_ratio()  const;
  NumericVector get_target_allocation()       const;
  IntegerVector get_allocation_ratio()        const;

  // setters
  void set_parameter_set(NumericVector);
  void set_fixed_allocation_ratio(NumericVector);
  void set_target_allocation(NumericVector);
  void set_allocation_ratio(int, int, IntegerVector);


  // function for changing randomization probabilitis adaptively
  virtual NumericVector adapt_rand_probability(int, IntegerVector, NumericVector);
};


// === CRD allocation class: BEGIN ===
// Complete Randomization procedure
class AllocationCRD: public Allocation {
public:

  AllocationCRD(NumericVector parameter_set, NumericVector fixed_allocation_ratio);

};
// === CRD allocation class: END ===


// === PBD allocation class: BEGIN ===
// Permuted Block Design procedure
class AllocationPBD: public Allocation {
public:

  AllocationPBD(NumericVector parameter_set, NumericVector fixed_allocation_ratio);

  // function for changing randomization probabilitis adaptively
  NumericVector adapt_rand_probability(int, IntegerVector, NumericVector) override;
};
// === PBD allocation class: END ===


// === BUD allocation class: BEGIN ===
// Block Urn Design
class AllocationBUD: public Allocation {
public:

  AllocationBUD(NumericVector parameter_set, NumericVector fixed_allocation_ratio);

  // function for changing randomization probabilitis adaptively
  NumericVector adapt_rand_probability(int, IntegerVector, NumericVector) override;
};
// === BUD allocation class: END ===


// === MWUD allocation class: BEGIN ===
// Mass Weighted Urn Design
class AllocationMWUD: public Allocation {
public:

  AllocationMWUD(NumericVector parameter_set, NumericVector fixed_allocation_ratio);

  // function for changing randomization probabilitis adaptively
  NumericVector adapt_rand_probability(int, IntegerVector, NumericVector) override;
};
// === MWUD allocation class: END ===


// === DL allocation class: BEGIN ===
// Drop-the-Loser rule
class AllocationDL: public Allocation {
private:
  NumericVector urn;
public:

  AllocationDL(NumericVector parameter_set, NumericVector fixed_allocation_ratio);

  // getter
  NumericVector get_urn() const;

  // setter
  void set_urn(NumericVector);

  // function to sample a ball from urn given discrete probability distribution
  int sample_ball();

  // function for changing randomization probabilitis adaptively
  NumericVector adapt_rand_probability(int, IntegerVector, NumericVector) override;
};
// === DL allocation class: END ===


// === DBCD allocation class: BEGIN ===
// Doubly Adaptive Biased Coin Design
class AllocationDBCD: public Allocation {
public:

  AllocationDBCD(NumericVector parameter_set, NumericVector fixed_allocation_ratio);

  // function for changing randomization probabilitis adaptively
  NumericVector adapt_rand_probability(int, IntegerVector, NumericVector) override;
};
// === DBCD allocation class: END ===


// === MinQD allocation class: BEGIN ===
// Minimum Quadratic Distance constrained balance optimization
class AllocationMinQD: public Allocation {
public:

  AllocationMinQD(NumericVector parameter_set, NumericVector fixed_allocation_ratio);

  // function for changing randomization probabilitis adaptively
  NumericVector adapt_rand_probability(int, IntegerVector, NumericVector) override;
};
// === MinQD allocation class: END ===


// === MaxEnt allocation class: BEGIN ===
// Maximum Entropy constrained balance optimization
class AllocationMaxEnt: public Allocation {
public:

  AllocationMaxEnt(NumericVector parameter_set, NumericVector fixed_allocation_ratio);

  // function for changing randomization probabilitis adaptively
  NumericVector adapt_rand_probability(int, IntegerVector, NumericVector) override;
};
// === MaxEnt allocation class: END ===

