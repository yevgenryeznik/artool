#include "allocation.hpp"


// Allocation base class
Allocation::Allocation(NumericVector parameter_set_, NumericVector fixed_allocation_ratio_):
  parameter_set(parameter_set_),
  fixed_allocation_ratio(fixed_allocation_ratio_),
  target_allocation(fixed_allocation_ratio_/sum(fixed_allocation_ratio_))
{}


// getters of Allocation class
NumericVector Allocation::get_parameter_set()           const { return(parameter_set); }
NumericVector Allocation::get_fixed_allocation_ratio()  const { return(fixed_allocation_ratio); }
NumericVector Allocation::get_target_allocation()       const { return(target_allocation); }
IntegerVector Allocation::get_allocation_ratio()        const { return(allocation_ratio); }


// compute alocation ratio before j-th subject (vector N)
void Allocation::set_allocation_ratio(int j, int number_of_treatments, IntegerVector treatment) {
  IntegerVector N(number_of_treatments);
  if ( j > 1) {
    for ( int k = 1; k <= number_of_treatments; k++) {
      N[k-1] = count(treatment[seq(0, j-2)].begin(), treatment[seq(0,j-2)].end(), k);
    }
  }
  allocation_ratio = N;
}


// setters of Allocation class
void Allocation::set_parameter_set(NumericVector parameter_set_) {
  parameter_set = parameter_set_;
}
void Allocation::set_fixed_allocation_ratio(NumericVector fixed_allocation_ratio_)  {
  fixed_allocation_ratio = fixed_allocation_ratio_;
}
void Allocation::set_target_allocation(NumericVector target_allocation_) {
  target_allocation = target_allocation_;
}


// function for changing randomization probabilitis adaptively
NumericVector Allocation::adapt_rand_probability(int j, IntegerVector treatment, NumericVector response) {
  return target_allocation;
}


// === CRD allocation class: BEGIN ===
// constructor of CRD Allocation class
AllocationCRD::AllocationCRD(NumericVector parameter_set_, NumericVector fixed_allocation_ratio_):
  Allocation(parameter_set_, fixed_allocation_ratio_) {}
// === CRD allocation class: END ===


// === PBD allocation class: BEGIN ===
// constructor of PBD Allocation class
AllocationPBD::AllocationPBD(NumericVector parameter_set_, NumericVector fixed_allocation_ratio_):
  Allocation(parameter_set_, fixed_allocation_ratio_) {}


// function for changing randomization probabilitis adaptively according to PBD procedure
NumericVector AllocationPBD::adapt_rand_probability(int j, IntegerVector treatment, NumericVector response){
  NumericVector w = get_fixed_allocation_ratio();                   // fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  NumericVector prob(number_of_treatments);                         // treatment randomization probabilities
  double b = get_parameter_set()[0] ;                               // PBD has one parameter (b)
  double bsize = b*sum(w);                                          // block size


  IntegerVector N(number_of_treatments);
  // compute alocation ratio before j-th subject (vector N)
  set_allocation_ratio(j, number_of_treatments, treatment);
  N = get_allocation_ratio();

  // compute randomization probabilities based on current allocation ratio (N)
  int k = std::floor((float)((j-1)/bsize));
  prob = (w*b*(1+k)-as<NumericVector>(N))/(bsize*(1+k)-(j-1));

  return prob;
};
// === PBD allocation class: END ===


// === BUD allocation class: BEGIN ===
// constructor of BUD Allocation class
AllocationBUD::AllocationBUD(NumericVector parameter_set_, NumericVector fixed_allocation_ratio_):
  Allocation(parameter_set_, fixed_allocation_ratio_) {}


// function for changing randomization probabilitis adaptively according to BUD procedure
NumericVector AllocationBUD::adapt_rand_probability(int j, IntegerVector treatment, NumericVector response){
  NumericVector w = get_fixed_allocation_ratio();                   // fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  NumericVector prob(number_of_treatments);                         // treatment randomization probabilities
  double lambda = get_parameter_set()[0] ;                          // BUD has one parameter (lambda)


  IntegerVector N(number_of_treatments);
  // compute alocation ratio before j-th subject (vector N)
  set_allocation_ratio(j, number_of_treatments, treatment);
  N = get_allocation_ratio();

  // compute randomization probabilities based on current allocation ratio (N)
  int k = min(Rcpp::floor(as<NumericVector>(N)/w));
  prob = (w*(lambda+k)-as<NumericVector>(N))/(sum(w)*(lambda+k)-(j-1));

  return prob;
};
// === BUD allocation class: END ===


// === MWUD allocation class: BEGIN ===
// constructor of MWUD Allocation class
AllocationMWUD::AllocationMWUD(NumericVector parameter_set_, NumericVector fixed_allocation_ratio_):
  Allocation(parameter_set_, fixed_allocation_ratio_) {}


// function for changing randomization probabilitis adaptively according to MWUD procedure
NumericVector AllocationMWUD::adapt_rand_probability(int j, IntegerVector treatment, NumericVector response){
  NumericVector w = get_fixed_allocation_ratio();                   // fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  NumericVector prob(number_of_treatments);                         // treatment randomization probabilities
  double alpha = get_parameter_set()[0] ;                           // MWUD has one parameter (alpha)


  IntegerVector N(number_of_treatments);
  // compute alocation ratio before j-th subject (vector N)
  set_allocation_ratio(j, number_of_treatments, treatment);
  N = get_allocation_ratio();

  // compute randomization probabilities based on current allocation ratio (N)
  NumericVector num(number_of_treatments);
  for (int k = 1; k <= number_of_treatments; k++) {
    num[k-1] = Rcpp::max(NumericVector::create(alpha*w[k-1] - N[k-1] + (j-1)*w[k-1], 0));
  }
  prob = num/sum(num);

  return prob;
};
// === MWUD allocation class: END ===


// === DL allocation class: BEGIN ===
// constructor of DL Allocation class
AllocationDL::AllocationDL(NumericVector parameter_set_, NumericVector fixed_allocation_ratio_):
  Allocation(parameter_set_, fixed_allocation_ratio_) {}


// getter
NumericVector AllocationDL::get_urn() const { return urn; }


// setter
void AllocationDL::set_urn(NumericVector urn_) { urn = urn_; }


// function to sample a ball from urn given discrete probability distribution
int AllocationDL::sample_ball() {
  NumericVector prob = urn/sum(urn);
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


// function for changing randomization probabilitis adaptively according to MWUD procedure
NumericVector AllocationDL::adapt_rand_probability(int j, IntegerVector treatment, NumericVector response){
  NumericVector w = get_fixed_allocation_ratio();                   // fixed allocation ratio
  int number_of_treatments = w.size();                              // number of treatments
  NumericVector prob(number_of_treatments);                         // treatment randomization probabilities
  double a = get_parameter_set()[0] ;                               // DL has one parameter (a)
  NumericVector urn_(number_of_treatments+1);                        // urn state

  IntegerVector N(number_of_treatments);
  // compute alocation ratio before j-th subject (vector N)
  set_allocation_ratio(j, number_of_treatments, treatment);

  // compute randomization probabilities based on current urn state
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
      urn[seq(1, number_of_treatments)] = urn[seq(1, number_of_treatments)]+a*w;
      set_urn(urn);
    }
    else {
      prob[ball-1] = 1;
      urn[ball] = urn[ball]-1;
      set_urn(urn);
      flag = false;
    }
  }

  return prob;
};
// === DL allocation class: END ===


// === DBCD allocation class: BEGIN ===
// constructor of DBCD Allocation class
AllocationDBCD::AllocationDBCD(NumericVector parameter_set_, NumericVector fixed_allocation_ratio_):
  Allocation(parameter_set_, fixed_allocation_ratio_) {}


// function for changing randomization probabilitis adaptively according to DBCD procedure
NumericVector AllocationDBCD::adapt_rand_probability(int j, IntegerVector treatment, NumericVector response) {

  int number_of_treatments = get_fixed_allocation_ratio().size();   // number of treatments
  NumericVector prob(number_of_treatments);                         // treatment randomization probabilities
  double gm = get_parameter_set()[0] ;                              // DBCD has one parameter (gamma)
  NumericVector rho = get_target_allocation();                      // target allocation


  IntegerVector N(number_of_treatments);

  set_allocation_ratio(j, number_of_treatments, treatment);
  N = get_allocation_ratio();

  if ( j > number_of_treatments) {
    // compute alocation ratio before j-th subject (vector N)

    // compute randomization probabilities based on current allocation ratio (N)
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
  }
  else {
    prob = rho;
  }
  return prob;
}
// === DBCD allocation class: END ===


// === MinQD allocation class: BEGIN ===
// constructor of MinQD Allocation class
AllocationMinQD::AllocationMinQD(NumericVector parameter_set_, NumericVector fixed_allocation_ratio_):
  Allocation(parameter_set_, fixed_allocation_ratio_) {}


// function for changing randomization probabilitis adaptively according to MinQD procedure
NumericVector AllocationMinQD::adapt_rand_probability(int j, IntegerVector treatment, NumericVector response){
  NumericVector w = get_fixed_allocation_ratio();                   // fixed allocation ratio
  NumericVector rho = get_target_allocation();                      // target allocation
  int number_of_treatments = w.size();                              // number of treatments
  NumericVector prob(number_of_treatments);                         // treatment randomization probabilities
  double eta = get_parameter_set()[0] ;                             // MinQD has one parameter (eta)
  double mu;

  IntegerVector N(number_of_treatments);

  // the hypothetical "lack of balance"
  NumericVector B(number_of_treatments);

  // compute randomization probabilities based on current allocation ratio (N)
  NumericVector num(number_of_treatments);
  for (int k = 1; k <= number_of_treatments; k++) {
    // compute alocation ratio before j-th subject (vector N)
    set_allocation_ratio(j, number_of_treatments, treatment);
    N = get_allocation_ratio();
    N[k-1] += 1;
    B[k-1] = max(Rcpp::abs(as<NumericVector>(N)/(double)(j)-rho));
  }
  if (var(B) <= 1e-16) {
    prob = rho;
  }
  else {
    mu = 2/((number_of_treatments-1)*var(B))*eta*(sum(B*rho)-min(B));
    prob = rho-0.5*mu*(B-mean(B));
  }

  return prob;
};
// === MinQD allocation class: END ===


// === MaxEnt allocation class: BEGIN ===
// constructor of MaxEnt Allocation class
AllocationMaxEnt::AllocationMaxEnt(NumericVector parameter_set_, NumericVector fixed_allocation_ratio_):
  Allocation(parameter_set_, fixed_allocation_ratio_) {}


// function for changing randomization probabilitis adaptively according to MaxEnt procedure
NumericVector AllocationMaxEnt::adapt_rand_probability(int j, IntegerVector treatment, NumericVector response){
  NumericVector w = get_fixed_allocation_ratio();                   // fixed allocation ratio
  NumericVector rho = get_target_allocation();                      // target allocation
  int number_of_treatments = w.size();                              // number of treatments
  NumericVector prob(number_of_treatments);                         // treatment randomization probabilities
  double eta = get_parameter_set()[0] ;                             // MaxEnt has one parameter (eta)

  IntegerVector N(number_of_treatments);

  // the hypothetical "lack of balance"
  NumericVector B(number_of_treatments);

  // compute randomization probabilities based on current allocation ratio (N)
  NumericVector num(number_of_treatments);
  for (int k = 1; k <= number_of_treatments; k++) {
    // compute alocation ratio before j-th subject (vector N)
    set_allocation_ratio(j, number_of_treatments, treatment);
    N = get_allocation_ratio();
    N[k-1] += 1;
    B[k-1] = max(Rcpp::abs(as<NumericVector>(N)/(double)(j)-rho));
  }
  if (var(B) <= 1e-16) {
    prob = rho;
  }
  else {
    // we have to find a zero of a one-variable (mu) function.
    // bisection method is used.

    // bisection starts
    // function to find zero of
    auto fcn = [](double mu, double eta, NumericVector B, NumericVector rho) {
      return min(B)*eta + (1-eta)*sum(rho*B) - sum(B*rho*exp(-mu*B))/sum(rho*exp(-mu*B));
    };

    // left bound of a search interval
    double lb = 0, fcn_lb = fcn(lb, eta, B, rho);

    // right bound of a search interval
    double rb = 20/max(B), fcn_rb = fcn(rb, eta, B, rho);

    // midpoint of a search interval
    double mu = 0.5*(lb+rb), fcn_mu = fcn(mu, eta, B, rho);

    while (abs(fcn_mu) > 1e-5 && abs(rb-lb) > 1e-5){
      if (fcn_lb*fcn_mu > 0) {
        lb = mu;
        fcn_lb = fcn(lb, eta, B, rho);
      }
      else if (fcn_lb*fcn_mu < 0){
        rb = mu;
        fcn_rb = fcn(rb, eta, B, rho);
      }
      mu = 0.5*(lb+rb);
      fcn_mu = fcn(mu, eta, B, rho);
    }
    // bisection ends

    prob = rho*exp(-mu*B)/sum(rho*exp(-mu*B));
  }

  return prob;
};
// === MaxEnt allocation class: END ===



