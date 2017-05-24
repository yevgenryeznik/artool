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
  int cohort_size;                      // cohort size
  List resp_distribution;               // List represents response distribution with the followng items:
                                        //  name -- string with distribution name
                                        //  parameters -- list of numeric vectors of response parameters
                                        //                for each treatment
  List rr_procedure;                    // List represents Restricted Randomization procedure used with items:
                                        //  name -- string with RR procedure name
                                        //  parameters -- numeric vector of RR procedure parameters

  // agregated classes
  Response* resp;                       // object of Response class to generate responses
  RR* rr;                               // object of RR class to adapt randomization probabilities
  
  // simulations outcome
  List rand_probability;                // list of matrices with randomization probabilities
  NumericMatrix response;               // matrix of subjects' responses
  IntegerMatrix treatment;              // matrix of treatment assignments
  NumericMatrix forcing_idx;            // matrix of forcing indicies
  NumericMatrix imbalance;              // matrix of imbalances
  IntegerMatrix reject;                 // matrix of null hypothesis rejects

public:
  // constructor af a TrialRR class
  TrialRR(IntegerVector fixed_allocation_ratio_,
          int number_of_subjects_,
          int number_of_simulations_,
          int cohort_size_,
          List resp_distribution_,
          List rr_procedure_):

    fixed_allocation_ratio(fixed_allocation_ratio_),
    number_of_treatments(fixed_allocation_ratio_.size()),
    number_of_subjects(number_of_subjects_),
    number_of_simulations(number_of_simulations_),
    cohort_size(cohort_size_)
    
    {

    // initialize resp object and related fields
    resp_distribution = resp_distribution_;
    string resp_distribution_name = resp_distribution_["name"];
    List resp_parameters = resp_distribution_["parameters"];
    
    if (resp_distribution_name == "Binary") {
      resp = new BinaryResponse(NumericVector::create(0.5));
    }
    else if (resp_distribution_name == "Normal") {
      resp = new NormalResponse(NumericVector::create(0, 1));
    }
    else if (resp_distribution_name == "Weibull") {
      resp = new WeibullResponse(NumericVector::create(1, 1));
    }
    else {
      throw invalid_argument("Inappropriate name of response distribution");
    }
    
    // initialize rr object and related fields
    rr_procedure = rr_procedure_;
    string rr_procedure_name = rr_procedure_["name"];
    NumericVector rr_procedure_parameters = rr_procedure_["parameters"];
    
    if (rr_procedure_name == "CRD") {
      rr = new CRD(rr_procedure_parameters, fixed_allocation_ratio_);
    }
    else if (rr_procedure_name == "PBD") {
      rr = new PBD(rr_procedure_parameters, fixed_allocation_ratio_);
    }
    else if (rr_procedure_name == "BUD") {
      rr = new BUD(rr_procedure_parameters, fixed_allocation_ratio_);
    }
    else if (rr_procedure_name == "MWUD") {
      rr = new MWUD(rr_procedure_parameters, fixed_allocation_ratio_);
    }
    else if (rr_procedure_name == "DL") {
      rr = new DL(rr_procedure_parameters, fixed_allocation_ratio_);
    }
    else if (rr_procedure_name == "DBCD") {
      rr = new DBCD(rr_procedure_parameters, fixed_allocation_ratio_);
    }
    else if (rr_procedure_name == "MinQD") {
      rr = new MinQD(rr_procedure_parameters, fixed_allocation_ratio_);
    }
    else if (rr_procedure_name == "MaxEnt") {
      rr = new MaxEnt(rr_procedure_parameters, fixed_allocation_ratio_);
    }
    else {
      throw invalid_argument("Inappropriate name of randomization procedure");
    }

    // initialize output characteristics
    rand_probability = List(number_of_simulations_);
    response = NumericMatrix(number_of_simulations_, number_of_subjects_);
    treatment = IntegerMatrix(number_of_simulations_,number_of_subjects_);
    forcing_idx = NumericMatrix(number_of_simulations_, number_of_subjects_);
    imbalance = NumericMatrix(number_of_simulations_, number_of_subjects_);
    reject = IntegerMatrix(number_of_simulations_, number_of_subjects_);
  }


  // getters
  IntegerVector get_fixed_allocation_ratio() { return( fixed_allocation_ratio ); }
  NumericVector get_target_allocation()      { return( as<NumericVector>(fixed_allocation_ratio)/sum(fixed_allocation_ratio) ); }
  int get_number_of_treatments()             { return( number_of_treatments ); }
  int get_number_of_subjects()               { return( number_of_subjects ); }
  int get_number_of_simulations()            { return( number_of_simulations ); }
  int get_cohort_size()                      { return( cohort_size ); }
  string get_resp_distribution()             { return( resp_distribution["name"] ); }
  List get_resp_parameters()                 { return( resp_distribution["parameters"] ); }
  string get_rr_procedure()                  { return( rr_procedure["name"] ); }
  NumericVector get_rr_parameters()          { return( rr_procedure["parameters"] ); }
  
  List get_rand_probability()                { return( rand_probability ); }
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
    NumericVector rho = get_target_allocation();

    for ( int j = 1; j <= number_of_subjects; j++ ) {
      // adaptation
      rr->adapt(N);
      prob.row(j-1) = rr->get_rand_probability();
      // teratment assignment
      treatment_ = assign(j, prob.row(j-1));
      set_treatment(s, j, treatment_);
      
      // generate response for a treatment assigned
      resp->set_parameters(as<List>(resp_distribution["parameters"])[treatment_-1]);
      response_ = resp->response(1)[0];
      set_response(s, j, response_);
      
      
      obs.row(j-1) = NumericVector::create(treatment_, response_);
      
      // imbalance and forcing index
      N[treatment_-1] += 1;
      imbalance_ = sum(Rcpp::pow(as<NumericVector>(N) - j*rho, 2))/j;
      forcing_idx_ = sqrt((float)sum(Rcpp::pow(prob.row(j-1)-rho, 2)));

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

  class_<TrialRR>("TrialRR")
  .constructor<IntegerVector,int,int,int,List,List>()
  .method("simulate_trial", &TrialRR::simulate_trial)
  .method("simulate", &TrialRR::simulate)
  .property("fixedAllocationRatio", &TrialRR::get_fixed_allocation_ratio)
  .property("targetAllocation", &TrialRR::get_target_allocation)
  .property("numberOfTreatments", &TrialRR::get_number_of_treatments)
  .property("numberOfSubjects", &TrialRR::get_number_of_subjects)
  .property("numberOfSimulations", &TrialRR::get_number_of_simulations)
  .property("cohortSize", &TrialRR::get_cohort_size)
  .property("responseDistribution", &TrialRR::get_resp_distribution)
  .property("responseDistributionParams", &TrialRR::get_resp_parameters)
  .property("randomizationProcedure", &TrialRR::get_rr_procedure)
  .property("randomizationProcedureParams", &TrialRR::get_rr_parameters)
  .property("randomizationProbability", &TrialRR::get_rand_probability)
  .property("response", &TrialRR::get_response)
  .property("treatment", &TrialRR::get_treatment)
  .property("forcing_idx", &TrialRR::get_forcing_idx)
  .property("imbalance", &TrialRR::get_imbalance)
  .property("reject", &TrialRR::get_reject)
  ;
}


