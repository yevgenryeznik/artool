#include "allocation.hpp"
#include "stat_tests.hpp"

// Base class
class Trial {
private:
  int number_of_treatments;             // number of treatment groups
  int number_of_subjects;               // number of subjects involved in a trial
  int number_of_simulations;            // number of trial simulations
  int cohort_size;                      // cohort size
  string rand_procedure;                // randomization procedure used
  NumericVector parameter_set;          // vector of parameters of a randomization procedure

  Allocation* allocation;               // allocation class
  List rand_probability;                // list of matrices with randomization probabilities
  NumericMatrix response;               // matrix of subjects' responses
  IntegerMatrix treatment;              // matrix of treatment assignments
  NumericMatrix forcing_idx;            // matrix of forcing indicies
  NumericMatrix imbalance;              // matrix of imbalances
  IntegerMatrix reject;                 // matrix of null hypothesis rejects

public:
  // constructor af a Trial class
  Trial(int number_of_subjects_,
        int number_of_simulations_,
        int cohort_size_,
        string rand_procedure_,
        NumericVector fixed_allocation_ratio_,
        NumericVector parameter_set_) {

    number_of_subjects = number_of_subjects_;
    number_of_simulations = number_of_simulations_;
    cohort_size = cohort_size_;
    rand_procedure = rand_procedure_;

    if (rand_procedure_ == "CRD") {
      allocation = new AllocationCRD(parameter_set_, fixed_allocation_ratio_);
    }
    else if (rand_procedure_ == "PBD") {
      allocation = new AllocationPBD(parameter_set_, fixed_allocation_ratio_);
    }
    else if (rand_procedure_ == "BUD") {
      allocation = new AllocationBUD(parameter_set_, fixed_allocation_ratio_);
    }
    else if (rand_procedure_ == "MWUD") {
      allocation = new AllocationMWUD(parameter_set_, fixed_allocation_ratio_);
    }
    else if (rand_procedure_ == "DL") {
      allocation = new AllocationDL(parameter_set_, fixed_allocation_ratio_);
    }
    else if (rand_procedure_ == "DBCD") {
      allocation = new AllocationDBCD(parameter_set_, fixed_allocation_ratio_);
    }
    else if (rand_procedure_ == "MinQD") {
      allocation = new AllocationMinQD(parameter_set_, fixed_allocation_ratio_);
    }
    else if (rand_procedure_ == "MaxEnt") {
      allocation = new AllocationMaxEnt(parameter_set_, fixed_allocation_ratio_);
    }
    else {
      throw invalid_argument("Inappropriate name of randomization procedure");
    }

    number_of_treatments = allocation->get_fixed_allocation_ratio().size();
    parameter_set = parameter_set_;

    rand_probability = List(number_of_simulations_);
    response = NumericMatrix(number_of_simulations_, number_of_subjects_);
    treatment = IntegerMatrix(number_of_simulations_,number_of_subjects_);
    forcing_idx = NumericMatrix(number_of_simulations_, number_of_subjects_);
    imbalance = NumericMatrix(number_of_simulations_, number_of_subjects_);
    reject = IntegerMatrix(number_of_simulations_, number_of_subjects_);
  }


  // getters
  int get_number_of_treatments()             { return( number_of_treatments ); }
  int get_number_of_subjects()               { return( number_of_subjects ); }
  int get_number_of_simulations()            { return( number_of_simulations ); }
  NumericVector get_fixed_allocation_ratio() { return( allocation->get_fixed_allocation_ratio() ); }
  NumericVector get_target_allocation()      { return( allocation->get_target_allocation() ); }
  List get_rand_probability()                { return( rand_probability ); }
  NumericVector get_parameter_set()          { return( parameter_set ); }
  NumericMatrix get_response()               { return( response ); }
  IntegerMatrix get_treatment()              { return( treatment ); }
  NumericMatrix get_forcing_idx()            { return( forcing_idx); }
  NumericMatrix get_imbalance()              { return( imbalance ); }
  IntegerMatrix get_reject()                 { return( reject ); }

  // setters
  void set_response(int s, int j, double response_) { response(s-1,j-1) = response_; }
  void set_treatment(int s, int j, int treatment_)  { treatment(s-1, j-1) = treatment_; }
  void set_forcing_idx(int s, int j, double forcing_index_) { forcing_idx(s-1, j-1) = forcing_index_; }
  void set_imbalance(int s, int j, double imbalance_)  { imbalance(s-1, j-1) = imbalance_; }
  void set_reject(int s, int j, int reject_)  { reject(s-1, j-1) = reject_; }

  // treatment assignment
  int assign(int j, NumericVector prob) {
    NumericVector cumulative_prob(number_of_treatments+1);
    double u = runif(1)[0];
    int k;

    for (k = 1; k <= number_of_treatments; k++) {
      cumulative_prob[k] = cumulative_prob[k-1] + prob[k-1];
      if (cumulative_prob[k-1] < u && u < cumulative_prob[k]) {
        break;
      }
    }
    return k;
  }


  // simulation of a single trial (with a number s)
  void simulate_trial(int s) {
    IntegerVector treatments = seq_len(number_of_treatments);
    int treatment_;
    double response_;
    double forcing_idx_;
    double imbalance_;
    int reject_;
    NumericMatrix obs(number_of_subjects, 2);

    NumericMatrix prob = NumericMatrix(number_of_subjects, number_of_treatments);
    IntegerVector N(number_of_treatments);
    NumericVector rho = allocation->get_target_allocation();

    for ( int j = 1; j <= number_of_subjects; j++ ) {

      if (rand_procedure == "CRD") {
        allocation->set_allocation_ratio(j, number_of_treatments, treatment.row(s-1));
        prob.row(j-1) = allocation->get_target_allocation();
      }
      else {
        prob.row(j-1) = allocation->adapt_rand_probability(j, treatment.row(s-1), response.row(s-1));
      }


      treatment_ = assign(j, prob.row(j-1));
      response_ = rnorm(1)[0];
      obs.row(j-1) = NumericVector::create(treatment_, response_);

      // imbalance and forcing index
      N = allocation->get_allocation_ratio();
      N[treatment_-1] += 1;
      imbalance_ = sum(Rcpp::pow(as<NumericVector>(N) - j*rho, 2))/j;
      forcing_idx_ = sqrt((float)sum(Rcpp::pow(prob.row(j-1)-rho, 2)));

      set_treatment(s, j, treatment_);
      set_response(s, j, response_);
      // after a single trial simulation forcing index must be transformed into
      // forcing_idx(s-1, j-1) = sum(get_forcing_idx()[seq(0,j-1)])/j;
      set_forcing_idx(s, j, forcing_idx_);
      set_imbalance(s, j, imbalance_);

      // test null hypothesis
      reject_ = anova_test(obs(seq(0,j-1),_), number_of_treatments, 0.05);
      set_reject(s, j, reject_);

    }
    rand_probability[s-1] = prob;
  }


  // simulation of a trial number_of_simulations times
  void simulate() {
    IntegerVector trials = seq_len(number_of_simulations);
    for_each(trials.begin(), trials.end(), [this](int &s){ simulate_trial(s); });
  }

};


RCPP_MODULE(trial) {

  class_<Trial>("Trial")
  .constructor<int,int,int,string,NumericVector,NumericVector>()
  .method("simulate_trial", &Trial::simulate_trial)
  .method("simulate", &Trial::simulate)
  .property("numberOfTreatments", &Trial::get_number_of_treatments)
  .property("numberOfSubjects", &Trial::get_number_of_subjects)
  .property("numberOfSimulations", &Trial::get_number_of_simulations)
  .property("parameter_set", &Trial::get_parameter_set)
  .property("fixed_allocation_ratio", &Trial::get_fixed_allocation_ratio)
  .property("target_allocation", &Trial::get_target_allocation)
  .property("rand_probability", &Trial::get_rand_probability)
  .property("response", &Trial::get_response)
  .property("treatment", &Trial::get_treatment)
  .property("forcing_idx", &Trial::get_forcing_idx)
  .property("imbalance", &Trial::get_imbalance)
  .property("reject", &Trial::get_reject)
  ;
}

