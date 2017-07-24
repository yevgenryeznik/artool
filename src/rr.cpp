// A file contains implementation of classes providing Restrict Randomization procedures 
// targeting unequal allocation

#include "rr.hpp"
#include "algorithms.hpp"

// Restricted Randomization Procedure (base class)
RR::RR(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_): 
  parameters(parameters_),
  fixed_allocation_ratio(fixed_allocation_ratio_),
  number_of_subjects(number_of_subjects_) 
  
  {
    target_allocation = as<NumericVector>(fixed_allocation_ratio_)/sum(fixed_allocation_ratio_);
    number_of_treatments = fixed_allocation_ratio_.size();
    treatment = IntegerVector(number_of_subjects_);
    rand_probability = NumericMatrix(number_of_subjects, number_of_treatments);
  }


// getters of Restricted Randomization Procedure (base class)
NumericVector RR::get_parameters()              const { return parameters; }
IntegerVector RR::get_fixed_allocation_ratio()  const { return fixed_allocation_ratio; }
NumericVector RR::get_target_allocation()       const { return target_allocation; }
int RR::get_number_of_treatments()              const { return number_of_treatments ; }
int RR::get_number_of_subjects()                const { return number_of_subjects; }
IntegerVector RR::get_treatment()               const { return treatment; }
NumericMatrix RR::get_rand_probability()        const { return rand_probability; }
NumericMatrix RR::get_alloc_proportion()        const { return alloc_proportion; }
NumericVector RR::get_forcing_idx()             const { return forcing_idx; }
NumericVector RR::get_imbalance()               const { return imbalance; }
NumericVector RR::get_selection_bias()          const { return selection_bias; }


// setters of Restricted Randomization Procedure (base class)
void RR::set_treatment(int j, int treatment_)              { treatment[j-1] = treatment_; }  
void RR::set_rand_probability(int j, NumericVector prob)   { rand_probability.row(j-1) = prob; }
void RR::set_alloc_proportion(int j, NumericVector prop)   { alloc_proportion.row(j-1) = prop; }
void RR::set_forcing_idx(int j, double forcing_idx_)       { forcing_idx[j-1] = forcing_idx_; }
void RR::set_imbalance(int j, double imbalance_)           { imbalance[j-1] = imbalance_; }
void RR::set_selection_bias(int j, double selection_bias_) { selection_bias[j-1] = selection_bias_; }


// sample an integer from a set of integers
int RR::sample(IntegerVector range, NumericVector prob) {
  NumericVector cumulative_prob(prob.size()+1);
  double u = runif(1)[0];
  int k;
  for (k = 1; k <= prob.size(); k++) {
    cumulative_prob[k] = cumulative_prob[k-1] + prob[k-1];
    if (cumulative_prob[k-1] < u && u < cumulative_prob[k]) {
      break;
    }
  }
  return range[k-1];
}


// function which randomizes subjects across treatments
void RR::randomize() {
  NumericVector rho = get_target_allocation();
  IntegerVector N(number_of_treatments);
  int treatment_;
  double forcing_idx_;
  double imbalance_;
  double selection_bias_;
  
  NumericVector fi(get_number_of_subjects());
  IntegerVector guess(get_number_of_subjects());
  
  for ( int j = 1; j <= number_of_subjects; j++ ) {
    // adaptation
    adapt(j, N);
    
    // teratment assignment
    treatment_ = sample(seq_len(number_of_treatments), rand_probability.row(j-1));
    set_treatment(j, treatment_);
    N[treatment_-1] += 1;
    
    // operational characteristics
    set_alloc_proportion(j, as<NumericVector>(N)/sum(N));
    
    // imbalance 
    imbalance_ = sqrt((float)sum(Rcpp::pow(as<NumericVector>(N) - j*rho, 2)))/j;
    set_imbalance(j, imbalance_);
    
    // forcing index
    fi[j-1] = sum(Rcpp::pow(get_rand_probability().row(j-1)-rho, 2));
    forcing_idx_ = sum(fi[seq(0, j-1)])/j;
    set_forcing_idx(j, forcing_idx_);
    
    // sequential bias
    
    if (get_rand_probability()(j-1, treatment_-1) == Rcpp::max(get_rand_probability().row(j-1))) {
      guess(j-1) = 1;
    }
    selection_bias_ = mean(as<NumericVector>(guess)[seq(0, j-1)]);  
    set_selection_bias(j, selection_bias_);
  }   
}


// implementation of CRD class
CRD::CRD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void CRD::adapt(int j, IntegerVector N) {
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());
  set_rand_probability(j, w/sum(w));
}



// implementation of PBD class
PBD::PBD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void PBD::adapt(int j, IntegerVector N){
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  double b = get_parameters()[0];                                   // PBD has one parameter (b)
  double bsize = b*sum(w);                                          // block size

  // set randomization probabilities based on current allocation ratio (N)
  int k = std::floor((float)((j-1)/bsize));
  set_rand_probability(j, (w*b*(1+k)-as<NumericVector>(N))/(bsize*(1+k)-(j-1)));

};



// implementation of BUD class
BUD::BUD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void BUD::adapt(int j, IntegerVector N){
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  double lambda = get_parameters()[0] ;                             // BUD has one parameter (lambda)

  // set randomization probabilities based on current allocation ratio (N)
  int k = min(Rcpp::floor(as<NumericVector>(N)/w));
  set_rand_probability(j, (w*(lambda+k)-as<NumericVector>(N))/(sum(w)*(lambda+k)-(j-1)));

};



// implementation of MWUD class
MWUD::MWUD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void MWUD::adapt(int j, IntegerVector N){
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  double alpha = get_parameters()[0] ;                              // MWUD has one parameter (alpha)

  // set randomization probabilities based on current allocation ratio (N)
  NumericVector num(number_of_treatments);
  for (int k = 1; k <= number_of_treatments; k++) {
    num[k-1] = Rcpp::max(NumericVector::create(alpha*w[k-1] - N[k-1] + (j-1)*w[k-1], 0));
  }
  set_rand_probability(j, num/sum(num));
  
};



// implementation of DL class
DL::DL(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

IntegerVector DL::get_urn() const { return urn; }

void DL::set_urn(IntegerVector urn_) { urn = urn_; }

// int DL::sample_ball() {
//   NumericVector prob = as<NumericVector>(urn)/sum(urn);
//   NumericVector cumulative_prob(urn.size()+1);
//   double u = runif(1)[0];
//   int ball;
// 
//   for (ball = 1; ball <= urn.size(); ball++) {
//     cumulative_prob[ball] = cumulative_prob[ball-1] + prob[ball-1];
//     if (cumulative_prob[ball-1] < u && u < cumulative_prob[ball]) {
//       break;
//     }
//   }
//   return ball;
// }

void DL::adapt(int j, IntegerVector N){
  IntegerVector w = get_fixed_allocation_ratio();                   // fixed allocation ratio
  double a = get_parameters()[0] ;                                  // DL has one parameter (a)
  IntegerVector urn_(get_number_of_treatments()+1);                 // urn state

  NumericVector prob(get_number_of_treatments());
  
  // set randomization probabilities based on current urn state
  bool flag = true;
  int ball;
  if (j == 1) {
    urn_[0] = 1;
    urn_[seq(1, get_number_of_treatments())] = w;
    set_urn(urn_);
  }
  else{
    urn_ = get_urn();
  }
  while(flag){
    ball = sample(seq_len(urn.size()), as<NumericVector>(urn)/sum(urn))-1;
    if (ball == 0) {
      urn_[seq(1, get_number_of_treatments())] = urn_[seq(1, get_number_of_treatments())]+a*w;
      set_urn(urn_);
    }
    else {
      prob[ball-1] = 1;
      urn_[ball] = urn_[ball]-1;
      set_urn(urn_);
      flag = false;
    }
  }
  set_rand_probability(j, prob);
  
};


// implementation of DBCD class
DBCD::DBCD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void DBCD::adapt(int j, IntegerVector N) {
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  int number_of_treatments = get_fixed_allocation_ratio().size();   // number of treatments
  double gm = get_parameters()[0] ;                                 // DBCD has one parameter (gamma)
  NumericVector rho = w/sum(w);                                     // target allocation 

  NumericVector prob(number_of_treatments);

  // set randomization probabilities based on current allocation ratio (N)
  // if all treatments have at least one assignment
  if (  is_true(all(N > 0)) ) {
    for ( int k = 1; k <= number_of_treatments; k++) {
      prob[k-1] = rho[k-1]*pow(rho[k-1]/(N[k-1]/(double)(j-1)), gm);
    }
    prob = prob/sum(prob);
  }
  else {
    prob = rho;
  }
  set_rand_probability(j, prob);
  
}



// implementation of MinQD class
MinQD::MinQD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void MinQD::adapt(int j, IntegerVector N){
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  double eta = get_parameters()[0] ;                                // MinQD has one parameter (eta)
  NumericVector rho = w/sum(w);                                     // target allocation 

  double mu;
  NumericVector prob(number_of_treatments);

  // the hypothetical "lack of balance"
  NumericVector B(number_of_treatments);

  // set randomization probabilities based on current allocation ratio (N)
  for (int k = 1; k <= number_of_treatments; k++) {
    // hypothetical imbalance if a subject is assigned to a treatment k
    IntegerVector N1(N.size());
    N1[k-1] += 1;
    N1 += N;
    B[k-1] = max(Rcpp::abs(as<NumericVector>(N1)/(double)(j)-rho));
  }
  if (var(B) <= 1e-16) {
    prob = rho;
  }
  else {
    mu = 2/((number_of_treatments-1)*var(B))*eta*(sum(B*rho)-min(B));
    prob = rho-0.5*mu*(B-mean(B));
  }

  set_rand_probability(j, prob);
};



// implementation of MaxEnt class
MaxEnt::MaxEnt(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void MaxEnt::adapt(int j, IntegerVector N){
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  double eta = get_parameters()[0] ;                                // MinQD has one parameter (eta)
  NumericVector rho = w/sum(w);                                     // target allocation 

  double mu;
  NumericVector prob(number_of_treatments);
  
  // the hypothetical "lack of balance"
  NumericVector B(number_of_treatments);

  // compute randomization probabilities based on current allocation ratio (N)
  NumericVector num(number_of_treatments);
  for (int k = 1; k <= number_of_treatments; k++) {
    // hypothetical imbalance if a subject is assigned to a treatment k
    IntegerVector N1(N.size());
    N1[k-1] += 1;
    N1 += N;
    B[k-1] = max(Rcpp::abs(as<NumericVector>(N1)/(double)(j)-rho));
  }
  if (var(B) <= 1e-16) {
    prob = rho;
  }
  else {
    // we have to find a zero of a one-variable (mu) function.
    // bisection method is used -- implemeted in the file algorithms.cpp

    // function to find zero of
    std::function<double (double)> fcn = [eta, B, rho](double mu) {
      return min(B)*eta + (1-eta)*sum(rho*B) - sum(B*rho*exp(-mu*B))/sum(rho*exp(-mu*B));
    };
    
    mu = bisection(fcn, 0, 20/max(B), 1e-5);
      
    prob = rho*exp(-mu*B)/sum(rho*exp(-mu*B));
  }

  set_rand_probability(j, prob);
};



