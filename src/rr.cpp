// A file contains implementation of classes providing Restrict Randomization procedures 
// targeting unequal allocation

#include "rr.hpp"
#include "algorithms.hpp"

// Restricted Randomization Procedure (base class)
RR::RR(NumericVector parameters_, IntegerVector fixed_allocation_ratio_) {
  parameters = parameters_;
  fixed_allocation_ratio = fixed_allocation_ratio_;
  rand_probability = NumericVector(fixed_allocation_ratio_.size());
}


// getters of Restricted Randomization Procedure (base class)
NumericVector RR::get_parameters()             const { return parameters; }
IntegerVector RR::get_fixed_allocation_ratio() const { return fixed_allocation_ratio; }
NumericVector RR::get_rand_probability()       const { return rand_probability; }

// setters of Restricted Randomization Procedure (base class)
void RR::set_parameters(NumericVector parameters_) { parameters = parameters_; }
void RR::set_rand_probability(NumericVector rand_probability_) { rand_probability = rand_probability_; }



// implementation of CRD class
CRD::CRD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_):
  RR(parameters_, fixed_allocation_ratio_) {}

void CRD::adapt(IntegerVector allocation_ratio) {
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());
  set_rand_probability(w/sum(w));
}



// implementation of PBD class
PBD::PBD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_):
  RR(parameters_, fixed_allocation_ratio_) {}

void PBD::adapt(IntegerVector allocation_ratio){
  IntegerVector N = allocation_ratio;
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  double b = get_parameters()[0];                                   // PBD has one parameter (b)
  double bsize = b*sum(w);                                          // block size
  int j = sum(N) + 1;                                               // current subjects' ID 

  // set randomization probabilities based on current allocation ratio (N)
  int k = std::floor((float)((j-1)/bsize));
  set_rand_probability((w*b*(1+k)-as<NumericVector>(N))/(bsize*(1+k)-(j-1)));

};



// implementation of BUD class
BUD::BUD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_):
  RR(parameters_, fixed_allocation_ratio_) {}

void BUD::adapt(IntegerVector allocation_ratio){
  IntegerVector N = allocation_ratio;
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  double lambda = get_parameters()[0] ;                             // BUD has one parameter (lambda)
  int j = sum(N) + 1;                                               // current subjects' ID 
  
  // set randomization probabilities based on current allocation ratio (N)
  int k = min(Rcpp::floor(as<NumericVector>(N)/w));
  set_rand_probability((w*(lambda+k)-as<NumericVector>(N))/(sum(w)*(lambda+k)-(j-1)));

};



// implementation of MWUD class
MWUD::MWUD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_):
  RR(parameters_, fixed_allocation_ratio_) {}

void MWUD::adapt(IntegerVector allocation_ratio){
  IntegerVector N = allocation_ratio;
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  double alpha = get_parameters()[0] ;                              // MWUD has one parameter (alpha)
  int j = sum(N) + 1;                                               // current subjects' ID 
  
  // set randomization probabilities based on current allocation ratio (N)
  NumericVector num(number_of_treatments);
  for (int k = 1; k <= number_of_treatments; k++) {
    num[k-1] = Rcpp::max(NumericVector::create(alpha*w[k-1] - N[k-1] + (j-1)*w[k-1], 0));
  }
  set_rand_probability(num/sum(num));
  
};



// implementation of DL class
DL::DL(NumericVector parameters_, IntegerVector fixed_allocation_ratio_):
  RR(parameters_, fixed_allocation_ratio_) {}

IntegerVector DL::get_urn() const { return urn; }

void DL::set_urn(IntegerVector urn_) { urn = urn_; }

int DL::sample_ball() {
  NumericVector prob = as<NumericVector>(urn)/sum(urn);
  NumericVector cumulative_prob(urn.size()+1);
  double u = runif(1)[0];
  int ball;

  for (ball = 1; ball <= urn.size(); ball++) {
    cumulative_prob[ball] = cumulative_prob[ball-1] + prob[ball-1];
    if (cumulative_prob[ball-1] < u && u < cumulative_prob[ball]) {
      break;
    }
  }
  return ball;
}

void DL::adapt(IntegerVector allocation_ratio){
  IntegerVector N = allocation_ratio;
  IntegerVector w = get_fixed_allocation_ratio();                   // fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  double a = get_parameters()[0] ;                                  // DL has one parameter (a)
  IntegerVector urn_(number_of_treatments+1);                       // urn state
  int j = sum(N) + 1;                                               // current subjects' ID 
  
  NumericVector prob(number_of_treatments);
  
  // set randomization probabilities based on current urn state
  bool flag = true;
  int ball;
  if (j == 1) {
    urn_[0] = 1;
    urn_[seq(1, number_of_treatments)] = w;
    set_urn(urn_);
  }
  else{
    urn_ = get_urn();
  }
  while(flag){
    ball = sample_ball()-1;
    if (ball == 0) {
      urn_[seq(1, number_of_treatments)] = urn_[seq(1, number_of_treatments)]+a*w;
      set_urn(urn_);
    }
    else {
      prob[ball-1] = 1;
      urn_[ball] = urn_[ball]-1;
      set_urn(urn_);
      flag = false;
    }
  }
  set_rand_probability(prob);
  
};


// implementation of DBCD class
DBCD::DBCD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_):
  RR(parameters_, fixed_allocation_ratio_) {}

void DBCD::adapt(IntegerVector allocation_ratio) {
  IntegerVector N = allocation_ratio;
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  int number_of_treatments = get_fixed_allocation_ratio().size();   // number of treatments
  double gm = get_parameters()[0] ;                                 // DBCD has one parameter (gamma)
  NumericVector rho = w/sum(w);                                     // target allocation 
  int j = sum(N) + 1;                                               // current subjects' ID 
  
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
  set_rand_probability(prob);
  
}



// implementation of MinQD class
MinQD::MinQD(NumericVector parameters_, IntegerVector fixed_allocation_ratio_):
  RR(parameters_, fixed_allocation_ratio_) {}

void MinQD::adapt(IntegerVector allocation_ratio){
  IntegerVector N = allocation_ratio;
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  double eta = get_parameters()[0] ;                                // MinQD has one parameter (eta)
  NumericVector rho = w/sum(w);                                     // target allocation 
  int j = sum(N) + 1;                                               // current subjects' ID 
  
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

  set_rand_probability(prob);
};



// implementation of MaxEnt class
MaxEnt::MaxEnt(NumericVector parameters_, IntegerVector fixed_allocation_ratio_):
  RR(parameters_, fixed_allocation_ratio_) {}

void MaxEnt::adapt(IntegerVector allocation_ratio){
  IntegerVector N = allocation_ratio;
  NumericVector w = as<NumericVector>(get_fixed_allocation_ratio());// fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  double eta = get_parameters()[0] ;                                // MinQD has one parameter (eta)
  NumericVector rho = w/sum(w);                                     // target allocation 
  int j = sum(N) + 1;                                               // current subjects' ID 
  
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

  set_rand_probability(prob);
};



