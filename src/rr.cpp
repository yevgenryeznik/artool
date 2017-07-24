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
    rand_probability = NumericMatrix(number_of_subjects_, number_of_treatments);
    alloc_proportion = NumericMatrix(number_of_subjects_, number_of_treatments);
    forcing_index = NumericVector(number_of_subjects_);
    imbalance = NumericVector(number_of_subjects_);
    mpm = NumericVector(number_of_subjects_);
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
NumericVector RR::get_forcing_index()           const { return forcing_index; }
NumericVector RR::get_imbalance()               const { return imbalance; }
NumericVector RR::get_mpm()                     const { return mpm; }

// setters of Restricted Randomization Procedure (base class)
void RR::set_treatment(int j, int treatment_)              { treatment[j-1] = treatment_; }  
void RR::set_treatment(IntegerVector treatment_)           { treatment = treatment_; }  
void RR::set_rand_probability(int j, NumericVector prob)   { rand_probability.row(j-1) = prob; }
void RR::set_rand_probability(NumericMatrix prob)          { rand_probability = prob; }
void RR::set_alloc_proportion(int j, NumericVector prop)   { alloc_proportion.row(j-1) = prop; }
void RR::set_alloc_proportion(NumericMatrix prop)          { alloc_proportion = prop; }
void RR::set_forcing_index(int j, double forcing_index_)   { forcing_index[j-1] = forcing_index_; }
void RR::set_forcing_index(NumericVector forcing_index_)   { forcing_index = forcing_index_; }
void RR::set_imbalance(int j, double imbalance_)           { imbalance[j-1] = imbalance_; }
void RR::set_imbalance(NumericVector imbalance_)           { imbalance = imbalance_; }
void RR::set_mpm(int j, double mpm_)                       { mpm[j-1] = mpm_; }
void RR::set_mpm(NumericVector mpm_)                       { mpm = mpm_; }


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
void RR::run() {
  NumericVector rho = get_target_allocation();
  IntegerVector N(number_of_treatments);
  IntegerVector N1(get_number_of_subjects());
  double forcing_index_;
  double imbalance_;
  double mpm_;

  NumericVector fi(get_number_of_subjects());
  
  for ( int j = 1; j <= number_of_subjects; j++ ) {
    N1 = N;
    // adaptation and assignment
    randomize(j, N);
    
    // update current allocation ratio
    N[treatment[j-1]-1] += 1;
    
    // allocation proportion
    set_alloc_proportion(j, as<NumericVector>(N)/sum(N));
    
    // imbalance 
    imbalance_ = sqrt((float)sum(Rcpp::pow(as<NumericVector>(N) - j*rho, 2)))/j;
    set_imbalance(j, imbalance_);
    
    // forcing index
    fi[j-1] = sum(Rcpp::pow(get_rand_probability().row(j-1)-rho, 2));
    forcing_index_ = mean(fi[seq(0, j-1)]);
    set_forcing_index(j, forcing_index_);
    
    // momentum of probability mass
    mpm_ = 0;
    for (int k = 1; k <= number_of_treatments; k++) {
      N1[k-1] += 1;
      mpm_ += rand_probability.row(j-1)[k-1]*sqrt((float)sum(Rcpp::pow(as<NumericVector>(N1) - j*rho, 2)));
      N1[k-1] -= 1;
    }
    set_mpm(j, mpm_);
  }   
}


// implementation of CRD class
CRD::CRD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void CRD::randomize(int j, IntegerVector N) {
  // fixed allocation ratio
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());
  
  // current alloction probabilities
  NumericVector prob(get_number_of_treatments());
  
  int treatment_;
  
  prob = w/sum(w);
  treatment_ = sample(seq_len(get_number_of_treatments()), prob);
  
  set_treatment(j, treatment_);
  set_rand_probability(j, prob);
}



// implementation of PBD class
PBD::PBD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void PBD::randomize(int j, IntegerVector N){
  // fixed allocation ratio
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());
  
  // PBD has one parameter (b)
  double b = get_parameters()[0];          
  
  // block size
  double bsize = b*sum(w);                                          
  
  // current alloction probabilities
  NumericVector prob(get_number_of_treatments());
  
  int treatment_;
  
  // set randomization probabilities based on current allocation ratio (N)
  int k = std::floor((float)((j-1)/bsize));
  prob = (w*b*(1+k)-as<NumericVector>(N))/(bsize*(1+k)-(j-1));
  treatment_ = sample(seq_len(get_number_of_treatments()), prob);
  
  set_treatment(j, treatment_);
  set_rand_probability(j, prob);
};



// implementation of BUD class
BUD::BUD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void BUD::randomize(int j, IntegerVector N){
  // fixed allocation ratio
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());
  
  // BUD has one parameter (lambda)
  double lambda = get_parameters()[0] ;                             

  // current alloction probabilities
  NumericVector prob(get_number_of_treatments());
  
  int treatment_;
  
  // set randomization probabilities based on current allocation ratio (N)
  int k = min(Rcpp::floor(as<NumericVector>(N)/w));
  prob = (w*(lambda+k)-as<NumericVector>(N))/(sum(w)*(lambda+k)-(j-1));
  treatment_ = sample(seq_len(get_number_of_treatments()), prob);
  
  set_treatment(j, treatment_);
  set_rand_probability(j, prob);
};



// implementation of MWUD class
MWUD::MWUD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void MWUD::randomize(int j, IntegerVector N){
  // fixed allocation ratio
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());
  
  // MWUD has one parameter (alpha)
  double alpha = get_parameters()[0] ;

  // current alloction probabilities
  NumericVector prob(get_number_of_treatments());
  
  int treatment_;

  // set randomization probabilities based on current allocation ratio (N)
  NumericVector num(get_number_of_treatments());
  for (int k = 1; k <= get_number_of_treatments(); k++) {
    num[k-1] = Rcpp::max(NumericVector::create(alpha*w[k-1] - N[k-1] + (j-1)*w[k-1], 0));
  }
  prob = num/sum(num);
  treatment_ = sample(seq_len(get_number_of_treatments()), prob);
  
  set_treatment(j, treatment_);
  set_rand_probability(j, prob);
};



// implementation of DL class
DL::DL(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {
  
  // initialize urn 
  urn = IntegerVector(get_number_of_treatments()+1);
}


void DL::randomize(int j, IntegerVector N){
  // fixed allocation ratio
  IntegerVector w = get_fixed_allocation_ratio();      
  
  // DL has one parameter (a)
  double a = get_parameters()[0] ;                      
  
  // current alloction probabilities
  NumericVector prob(get_number_of_treatments());
  
  int treatment_;
  
  // set randomization probabilities based on current urn state
  bool flag = true;
  if (j == 1) {
    urn[0] = 1;
    urn[seq(1, get_number_of_treatments())] = w;
  }
  while(flag){
    treatment_ = sample(seq_len(urn.size()), as<NumericVector>(urn)/sum(urn))-1;
    if (treatment_ == 0) {
      urn[seq(1, get_number_of_treatments())] = urn[seq(1, get_number_of_treatments())]+a*w;
    }
    else {
      prob = as<NumericVector>(urn)[seq(1, get_number_of_treatments())]/
        sum(urn[seq(1, get_number_of_treatments())]);
      urn[treatment_] = urn[treatment_]-1;
      flag = false;
    }
  }
  
  set_treatment(j, treatment_);
  set_rand_probability(j, prob);
};


// implementation of DBCD class
DBCD::DBCD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void DBCD::randomize(int j, IntegerVector N) {
  // fixed allocation ratio
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());
  
  // DBCD has one parameter (gamma)
  double gm = get_parameters()[0] ;                                 
  
  // current alloction probabilities
  NumericVector prob(get_number_of_treatments());
  
  int treatment_;

  // target allocation proportion
  NumericVector rho = w/sum(w);                                     

  // set randomization probabilities based on current allocation ratio (N)
  // if all treatments have at least one assignment
  if (  is_true(all(N > 0)) ) {
    for ( int k = 1; k <= get_number_of_treatments(); k++) {
      prob[k-1] = rho[k-1]*pow(rho[k-1]/(N[k-1]/(double)(j-1)), gm);
    }
    prob = prob/sum(prob);
  }
  else {
    prob = rho;
  }
  treatment_ = sample(seq_len(get_number_of_treatments()), prob);
  
  set_treatment(j, treatment_);
  set_rand_probability(j, prob);
}



// implementation of MinQD class
MinQD::MinQD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void MinQD::randomize(int j, IntegerVector N){
  // fixed allocation ratio
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());
  
  // MinQD has one parameter (eta)
  double eta = get_parameters()[0] ;
  
  // current alloction probabilities
  NumericVector prob(get_number_of_treatments());
  
  int treatment_;
  
  // target allocation 
  NumericVector rho = w/sum(w);                                     

  double mu;

  // the hypothetical "lack of balance"
  NumericVector B(get_number_of_treatments());

  // set randomization probabilities based on current allocation ratio (N)
  for (int k = 1; k <= get_number_of_treatments(); k++) {
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
    mu = 2/((get_number_of_treatments()-1)*var(B))*eta*(sum(B*rho)-min(B));
    prob = rho-0.5*mu*(B-mean(B));
  }
  treatment_ = sample(seq_len(get_number_of_treatments()), prob);
  
  set_treatment(j, treatment_);
  set_rand_probability(j, prob);
};



// implementation of MaxEnt class
MaxEnt::MaxEnt(NumericVector parameters_, IntegerVector fixed_allocation_ratio_, int number_of_subjects_):
  RR(parameters_, fixed_allocation_ratio_, number_of_subjects_) {}

void MaxEnt::randomize(int j, IntegerVector N){
  // fixed allocation ratio
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());
  
  // MaxEnt has one parameter (eta)
  double eta = get_parameters()[0] ;           
  
  // current alloction probabilities
  NumericVector prob(get_number_of_treatments());
  
  int treatment_;
  
  // target allocation 
  NumericVector rho = w/sum(w);                                     

  double mu;

  // the hypothetical "lack of balance"
  NumericVector B(get_number_of_treatments());

  // compute randomization probabilities based on current allocation ratio (N)
  NumericVector num(get_number_of_treatments());
  for (int k = 1; k <= get_number_of_treatments(); k++) {
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
    
    //mu = bisection(fcn, 0, 100/max(B), 1e-5);
    mu = secant(fcn, 0, 1, 1e-5);
    prob = rho*exp(-mu*B)/sum(rho*exp(-mu*B));
  }
  treatment_ = sample(seq_len(get_number_of_treatments()), prob);
  
  set_treatment(j, treatment_);
  set_rand_probability(j, prob);
};



