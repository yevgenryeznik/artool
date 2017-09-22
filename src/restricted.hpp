// A file contains definitions of functions providing Restrict Randomization procedures 
// targeting unequal allocation

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;


// functions below compute randomization probabilities given
//    int (j) -- current subkect's ID
//    IntegerVector (N) -- vector of current allocaton ratios
//    IntegerVector (w) -- vector of target fixed allocation ratios
//    double (p) -- parameter of a randoomization procedure

// Permuted Block Design: PBD(p == b)
List pbd(int, IntegerVector, IntegerVector, double);


// Block Urn Design: BUD(p == lambda)
List bud(int, IntegerVector, IntegerVector, double);


// Mass Weight Urn Design: MWUD(p == alpha)
List mwud(int, IntegerVector, IntegerVector, double);


// Drop-the-Loser: DL(p == a)
// contains additional argument 
//    NumericVector (urn) -- urn state for a design
List dl(int, IntegerVector, IntegerVector, double, IntegerVector);


// Doubly-Adaptive Biased Coind Design: DBCD(p == gamma)
List dbcd(int, IntegerVector, IntegerVector, double);


// Minimum Quadrati Distance: MinQD(p == eta)
List min_qd(int, IntegerVector, IntegerVector, double);


// Maximum Entropy: MaxEnt(p == eta)
List max_ent(int, IntegerVector, IntegerVector, double);


// set randomization procedure given 
//    (IntegerVector) fixed allocation ratio
//    (std::string) randomization procedure
//    (double) parameter of the randomization procedure
std::function<List (int, IntegerVector)> set_rand_procedure(IntegerVector, std::string, double);
  

// run restricted randomization procedure
List restricted(int, IntegerVector, std::string, double, std::string, List);
  

// simulation of restricted randomization procedure
List simulate_restricted(int, int, IntegerVector w, std::string, double, std::string, List);


