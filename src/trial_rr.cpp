#include "rr.hpp"
#include "response.hpp"
#include "stat_tests.hpp"

// Trial with Restricted Randomization
class TrialRR {
private:
  // setup parameters
  IntegerVector fixed_allocation_ratio; // fixed allocation ratio
  int number_of_treatments;             // number of treatment groups
  int number_of_subjects;               // number of subjects involved in a trial
  int number_of_simulations;            // number of trial simulations
  List resp_distribution;               // List represents response distribution with the followng items:
                                        //  name -- string with distribution name
                                        //  parameters -- list of numeric vectors of response parameters
                                        //                for each treatment
  List rr_procedure;                    // List represents Restricted Randomization procedure used with items:
                                        //  name -- string with RR procedure name
                                        //  parameters -- numeric vector of RR procedure parameters
  double significance_level;            // significance level to test hypothesis H0: mu1 = mu2 = ... = muK
  bool time_drift;                      // if TRUE then this is so called model with a "time drift", i.e. 
                                        // j/n is added to response value, where j is a subject's ID
  
  // agregated classes
  Response* resp;                       // object of Response class to generate responses
  RR* rr;                               // object of RR class to adapt randomization probabilities
  
  // simulations outcome
  List rand_probability;                // list of matrices with randomization probabilities
  List alloc_proportion;                // list of matrices with allocation proportions
  NumericMatrix response;               // matrix of subjects' responses
  IntegerMatrix treatment;              // matrix of treatment assignments
  NumericMatrix forcing_index;          // matrix of forcing indicies
  NumericMatrix imbalance;              // matrix of imbalances
  NumericMatrix mpm;                    // matrix of momentums of probability mass
  IntegerMatrix reject;                 // matrix of null hypothesis rejects

public:
  // constructor af a TrialRR class
  TrialRR(IntegerVector fixed_allocation_ratio_,
          int number_of_subjects_,
          int number_of_simulations_,
          List resp_distribution_,
          List rr_procedure_, 
          double significance_level_):

    fixed_allocation_ratio(fixed_allocation_ratio_),
    number_of_treatments(fixed_allocation_ratio_.size()),
    number_of_subjects(number_of_subjects_),
    number_of_simulations(number_of_simulations_),
    significance_level(significance_level_),
    time_drift(false)
    
    {

    // initialize response object and related fields
    resp_distribution = resp_distribution_;
    string resp_distribution_name = resp_distribution_["name"];
    NumericMatrix resp_parameters = resp_distribution_["parameters"];
    
    if (resp_distribution_name == "Binary") {
      resp = new BinaryResponse(resp_parameters);
    }
    else if (resp_distribution_name == "Normal") {
      resp = new NormalResponse(resp_parameters);
    }
    else if (resp_distribution_name == "Weibull") {
      resp = new WeibullResponse(resp_parameters);
    }
    else {
      throw invalid_argument("Inappropriate name of response distribution");
    }
    
    // initialize rr object and related fields
    rr_procedure = rr_procedure_;
    string rr_procedure_name = rr_procedure_["name"];
    NumericVector rr_procedure_parameters = rr_procedure_["parameters"];
    
    if (rr_procedure_name == "CRD") {
      rr = new CRD(rr_procedure_parameters, fixed_allocation_ratio_, number_of_subjects_);
    }
    else if (rr_procedure_name == "PBD") {
      rr = new PBD(rr_procedure_parameters, fixed_allocation_ratio_, number_of_subjects_);
    }
    else if (rr_procedure_name == "BUD") {
      rr = new BUD(rr_procedure_parameters, fixed_allocation_ratio_, number_of_subjects_);
    }
    else if (rr_procedure_name == "MWUD") {
      rr = new MWUD(rr_procedure_parameters, fixed_allocation_ratio_, number_of_subjects_);
    }
    else if (rr_procedure_name == "DL") {
      rr = new DL(rr_procedure_parameters, fixed_allocation_ratio_, number_of_subjects_);
    }
    else if (rr_procedure_name == "DBCD") {
      rr = new DBCD(rr_procedure_parameters, fixed_allocation_ratio_, number_of_subjects_);
    }
    else if (rr_procedure_name == "MinQD") {
      rr = new MinQD(rr_procedure_parameters, fixed_allocation_ratio_, number_of_subjects_);
    }
    else if (rr_procedure_name == "MaxEnt") {
      rr = new MaxEnt(rr_procedure_parameters, fixed_allocation_ratio_, number_of_subjects_);
    }
    else {
      throw invalid_argument("Inappropriate name of randomization procedure");
    }

    // initialize output characteristics
    rand_probability = List(number_of_simulations_);
    alloc_proportion = List(number_of_simulations_);
    response = NumericMatrix(number_of_simulations_, number_of_subjects_);
    treatment = IntegerMatrix(number_of_simulations_,number_of_subjects_);
    forcing_index = NumericMatrix(number_of_simulations_, number_of_subjects_);
    imbalance = NumericMatrix(number_of_simulations_, number_of_subjects_);
    mpm = NumericMatrix(number_of_simulations_, number_of_subjects_);
    reject = IntegerMatrix(number_of_simulations_, number_of_subjects_);
  }


  // getters
  IntegerVector get_fixed_allocation_ratio() { return( fixed_allocation_ratio ); }
  NumericVector get_target_allocation()      { return( as<NumericVector>(fixed_allocation_ratio)/sum(fixed_allocation_ratio) ); }
  int get_number_of_treatments()             { return( number_of_treatments ); }
  int get_number_of_subjects()               { return( number_of_subjects ); }
  int get_number_of_simulations()            { return( number_of_simulations ); }
  string get_resp_distribution()             { return( resp_distribution["name"] ); }
  NumericMatrix get_resp_parameters()        { return( resp_distribution["parameters"] ); }
  string get_rr_procedure()                  { return( rr_procedure["name"] ); }
  NumericVector get_rr_parameters()          { return( rr_procedure["parameters"] ); }
  double get_significance_level()            { return( significance_level ); }
  bool get_time_drift()                      { return( time_drift ); }
  
  List get_rand_probability()                { return( rand_probability ); }
  List get_proportions()                     { return( alloc_proportion ); }
  NumericMatrix get_response()               { return( response ); }
  IntegerMatrix get_treatment()              { return( treatment ); }
  NumericMatrix get_forcing_index()          { return( forcing_index); }
  NumericMatrix get_imbalance()              { return( imbalance ); }
  NumericMatrix get_mpm()                    { return( mpm ); }
  IntegerMatrix get_reject()                 { return( reject ); }

  // setters
  void set_significance_level(double significance_level_) { significance_level = significance_level_; }
  void set_time_drift(bool time_drift_) { time_drift = time_drift_; }
  void set_response(int s, NumericVector response_) { response.row(s-1) = response_; }
  void set_treatment(int s, IntegerVector treatment_)  { treatment.row(s-1) = treatment_; }
  void set_forcing_index(int s, NumericVector forcing_index_) { forcing_index.row(s-1) = forcing_index_; }
  void set_imbalance(int s, NumericVector imbalance_)  { imbalance.row(s-1) = imbalance_; }
  void set_mpm(int s, NumericVector mpm_)  { mpm.row(s-1) = mpm_; }
  void set_reject(int s, IntegerVector reject_)  { reject.row(s-1) = reject_; }

  // simulation of a single trial (with a number s)
  void simulate_trial(int s) {
    NumericMatrix obs(number_of_subjects, 2);
    NumericVector response_(number_of_subjects);
    IntegerVector reject_(number_of_subjects);
    
    // clean the data before randomization
    rr->set_rand_probability(NumericMatrix(number_of_subjects, number_of_treatments));
    rr->set_alloc_proportion(NumericMatrix(number_of_subjects, number_of_treatments));
    rr->set_treatment(IntegerVector(number_of_subjects));
    rr->set_forcing_index(NumericVector(number_of_subjects));
    rr->set_imbalance(NumericVector(number_of_subjects));
    
    // randomize patients
    rr->run();
    
    rand_probability[s-1] = rr->get_rand_probability();
    alloc_proportion[s-1] = rr->get_alloc_proportion();
    set_treatment(s, rr->get_treatment());
    set_imbalance(s, rr->get_imbalance());
    set_forcing_index(s, rr->get_forcing_index());
    set_mpm(s, rr->get_mpm());
    
    IntegerVector subjects = seq_len(number_of_subjects);
    
    // generate responses given treatment assignments
    for_each(subjects.begin(), subjects.end(), [this, &response_, &obs](int &j){
      response_[j-1] = resp->response(rr->get_treatment()[j-1]);
      if (time_drift) {
        response_[j-1] += (double)j/number_of_subjects;
      }
      obs.row(j-1) = NumericVector::create(rr->get_treatment()[j-1], response_[j-1]);
      });
    set_response(s, response_);

    // test null hypothesis
    for_each(subjects.begin(), subjects.end(), [this, &reject_, &obs](int &j){
      reject_[j-1] = anova_test(obs(seq(0,j-1),_), number_of_treatments, significance_level);
    });
    set_reject(s, reject_);
      
  }


  // simulation of a trial number_of_simulations times
  void simulate() {
    IntegerVector trials = seq_len(number_of_simulations);
    for_each(trials.begin(), trials.end(), [this](int &s){ simulate_trial(s); });
  }
};




RCPP_MODULE(trial) {
  
  class_<TrialRR>("TrialRR")
  .constructor<IntegerVector,int,int,List,List,double>()
  .method("set_significance_level", &TrialRR::set_significance_level)
  .method("set_time_drift", &TrialRR::set_time_drift)
  .method("simulate_trial", &TrialRR::simulate_trial)
  .method("simulate", &TrialRR::simulate)
  .property("fixed_allocation_ratio", &TrialRR::get_fixed_allocation_ratio)
  .property("target_allocation", &TrialRR::get_target_allocation)
  .property("number_of_treatments", &TrialRR::get_number_of_treatments)
  .property("number_of_subjects", &TrialRR::get_number_of_subjects)
  .property("number_of_simulations", &TrialRR::get_number_of_simulations)
  .property("response_distr", &TrialRR::get_resp_distribution)
  .property("response_distr_params", &TrialRR::get_resp_parameters)
  .property("rand_procedure", &TrialRR::get_rr_procedure)
  .property("rand_procedure_params", &TrialRR::get_rr_parameters)
  .property("significance_level", &TrialRR::get_significance_level)
  .property("time_drift", &TrialRR::get_time_drift)
  .property("rand_probability", &TrialRR::get_rand_probability)
  .property("alloc_proportion", &TrialRR::get_proportions)
  .property("response", &TrialRR::get_response)
  .property("treatment", &TrialRR::get_treatment)
  .property("forcing_index", &TrialRR::get_forcing_index)
  .property("imbalance", &TrialRR::get_imbalance)
  .property("mpm", &TrialRR::get_mpm)
  .property("reject", &TrialRR::get_reject)
  ;
  
}


