// A file contains implementation of functions providing Restricted Randomization procedures 
// targeting unequal allocation

#include "restricted.hpp"
#include "response.hpp"
#include "algorithms.hpp"
#include "statistics.hpp"


// functions below compute randomization probabilities given
//    j -- current subkect's ID
//    N -- vector of current allocaton ratios
//    w -- vector of target fixed allocation ratios
//    p -- parameter of a randoomization procedure

// Completely Randomized Design: CRD
List crd(int j, IntegerVector N, IntegerVector w, double p) {
  int number_of_treatments = w.size();
  
  NumericVector prob = as<NumericVector>(w)/sum(w);
  int trt = sample(seq_len(number_of_treatments), prob);
  
  return List::create(_["treatement"] = trt, 
                      _["probability"] = prob);
  
}


// Permuted Block Design: PBD(p == b)
List pbd(int j, IntegerVector N, IntegerVector w, double p) {
  int number_of_treatments = w.size();
  
  // block size
  double bsize = p*sum(w);
  
  int k = std::floor((float)((j-1)/bsize));
  NumericVector prob = (as<NumericVector>(w)*p*(1+k)-as<NumericVector>(N))/
    (bsize*(1+k)-(j-1));
  int trt = sample(seq_len(number_of_treatments), prob);
  
  return List::create(_["treatement"] = trt, 
                      _["probability"] = prob);
  
}



// Block Urn Design: BUD(p == lambda)
List bud(int j, IntegerVector N, IntegerVector w, double p) {
  int number_of_treatments = w.size();
  
  int k = min(Rcpp::floor(as<NumericVector>(N)/as<NumericVector>(w)));
  NumericVector prob = (as<NumericVector>(w)*(p+k)-as<NumericVector>(N))/
    (sum(w)*(p+k)-(j-1));
  int trt = sample(seq_len(number_of_treatments), prob);
  
  return List::create(_["treatement"] = trt, 
                      _["probability"] = prob);
  
}


// Mass Weight Urn Design: MWUD(p == alpha)
List mwud(int j, IntegerVector N, IntegerVector w, double p) {
  int number_of_treatments = w.size();
  
  NumericVector prob(number_of_treatments);
  IntegerVector k = seq_len(number_of_treatments);
  
  for_each(k.begin(), k.end(), [&j, &N, &w, &p, &prob](int &k){
    prob[k-1] = Rcpp::max(NumericVector::create(p*w[k-1] - N[k-1] + (j-1)*w[k-1], 0));
  });
  prob = prob/sum(prob);
  int trt = sample(k, prob);

  return List::create(_["treatement"] = trt, 
                      _["probability"] = prob);
}


// Drop-the-Loser: DL(p == a)
List dl(int j, IntegerVector N, IntegerVector w, double p, IntegerVector urn){
  int number_of_treatments = w.size();
  
  int trt;
  NumericVector prob(number_of_treatments);
  bool flag = true;
  while(flag){
    trt = sample(seq_len(urn.size()), as<NumericVector>(urn)/sum(urn))-1;
    if (trt == 0) {
      urn[seq(1, number_of_treatments)] = urn[seq(1, number_of_treatments)]+(int)p*w;
    }
    else {
      prob = as<NumericVector>(urn)[seq(1, number_of_treatments)]/
        sum(urn[seq(1, number_of_treatments)]);
      urn[trt] = urn[trt]-1;
      flag = false;
    }
  }

  return List::create(_["treatement"] = trt, 
                      _["probability"] = prob);
  
}


// Doubly-Adaptive Biased Coind Design: DBCD(p == gamma)
List dbcd(int j, IntegerVector N, IntegerVector w, double p) {
  int number_of_treatments = w.size();
  
  // target allocation proportion
  NumericVector rho = as<NumericVector>(w)/sum(w);                                     
  
  NumericVector prob(number_of_treatments);
  IntegerVector k = seq_len(number_of_treatments);
  
  if (  is_true(all(N > 0)) ) {
    for_each(k.begin(), k.end(), [&j, &N, &rho, &p, &prob](int &k){
      prob[k-1] = rho[k-1]*pow(rho[k-1]/(N[k-1]/(double)(j-1)), p);
    });
    prob = prob/sum(prob);
  }
  else {
    prob = rho;
  }
  int trt = sample(k, prob);
  
  return List::create(_["treatement"] = trt, 
                      _["probability"] = prob);
}


// Minimum Quadrati Distance: MinQD(p == eta)
List min_qd(int j, IntegerVector N, IntegerVector w, double p){
  int number_of_treatments = w.size();
  
  // target allocation 
  NumericVector rho = as<NumericVector>(w)/sum(w);                                     
  
  NumericVector prob(number_of_treatments);
  IntegerVector k = seq_len(number_of_treatments);

  // the hypothetical "lack of balance"
  NumericVector B(number_of_treatments);
  
  for_each(k.begin(), k.end(), [j, N, rho, &B](int &k){
    IntegerVector N1 = N;
    N1[k-1] += 1;
    B[k-1] = max(Rcpp::abs(as<NumericVector>(N1)/(double)(j)-rho));
    N1[k-1] -= 1;
  }); 
  
  if (var(B) <= 1e-16) {
    prob = rho;
  }
  else {
    double mu = 2/((w.size()-1)*var(B))*p*(sum(B*rho)-min(B));
    prob = rho-0.5*mu*(B-mean(B));
  }
  
  int trt = sample(k, prob);
  
  return List::create(_["treatement"] = trt, 
                      _["probability"] = prob);
}


// Maximum Entropy: MaxEnt(p == eta)
List max_ent(int j, IntegerVector N, IntegerVector w, double p){
  int number_of_treatments = w.size();
  
  // target allocation 
  NumericVector rho = as<NumericVector>(w)/sum(w);                                     
  
  NumericVector prob(number_of_treatments);
  IntegerVector k = seq_len(number_of_treatments);
  
  // the hypothetical "lack of balance"
  NumericVector B(number_of_treatments);
  
  for_each(k.begin(), k.end(), [j, N, rho, &B](int &k){
    IntegerVector N1 = N;
    N1[k-1] += 1;
    B[k-1] = max(Rcpp::abs(as<NumericVector>(N1)/(double)(j)-rho));
    N1[k-1] -= 1;
  }); 
  
  if (p != 1) {
    if (var(B) <= 1e-16) {
      prob = rho;
    }
    else {
      // we have to find a zero of a one-variable (mu) function.
      // bisection method is used -- implemeted in the file algorithms.cpp
      
      // function to find zero of
      std::function<double (double)> fcn = [p, B, rho](double mu) {
        return min(B)*p + (1-p)*sum(rho*B) - sum(B*rho*exp(-mu*B))/sum(rho*exp(-mu*B));
      };
      
      double mu = bisection(fcn, 0, 100/max(B), 1e-5);
      prob = rho*exp(-mu*B)/sum(rho*exp(-mu*B));
    }
  }
  else {
    prob[B == min(B)] = rho[B == min(B)];
    prob = prob/(double)sum(prob);
  }
  int trt = sample(k, prob);
  
  return List::create(_["treatement"] = trt, 
                      _["probability"] = prob);
}


// set randomization procedure given 
//    w -- fixed allocation ratio
//    procedure -- randomization procedure
//    p -- parameter of the randomization procedure
std::function<List (int, IntegerVector)> set_rand_procedure(IntegerVector w, std::string procedure, double p) {
  std::function<List (int, IntegerVector)> fcn;
  if (procedure == "CRD") {      // Completely Randomized Design
    fcn = [w, p](int j, IntegerVector N){
      return crd(j, N, w, p);
    };
  }
  else if (procedure == "PBD") { // Permuted Block Design
    fcn = [w, p](int j, IntegerVector N){
      return pbd(j, N, w, p);
    };
  }
  else if (procedure == "BUD") { // Block Urn Design
    fcn = [w, p](int j, IntegerVector N){
      return bud(j, N, w, p);
    };
  }
  else if (procedure == "MWUD") { // Mass Weighted Urn Design
    fcn = [w, p](int j, IntegerVector N){
      return mwud(j, N, w, p);
    };
  }
  else if (procedure == "DL") {   // Drop-the-Loser
    IntegerVector urn(w.size()+1);
    urn[0] = 1;
    urn[seq(1, w.size())] = w;
    fcn = [w, p, urn](int j, IntegerVector N){
      return dl(j, N, w, p, urn);
    };
  }
  else if (procedure == "DBCD") { // Doubly-Adaptive Biased Coin
    fcn = [w, p](int j, IntegerVector N){
      return dbcd(j, N, w, p);
    };
  }
  else if (procedure == "MinQD") { // Minimum Quadratic Distance
    fcn = [w, p](int j, IntegerVector N){
      return min_qd(j, N, w, p);
    };
  }
  else if (procedure == "MaxEnt") { // Maximum Entropy
    fcn = [w, p](int j, IntegerVector N){
      return max_ent(j, N, w, p);
    };
  }
  else {
    throw invalid_argument("Inappropriate name of a randomization procedure");
  }
  return fcn;
}


// restricted randomization procedure
//[[Rcpp::export(.restricted)]]
List restricted(int number_of_subjects, IntegerVector w, std::string procedure, double p, 
            std::string distribution, List parameter){
  int number_of_treatments = w.size();
  IntegerVector k = seq_len(number_of_treatments);

  // target allocation proportions
  NumericVector rho = as<NumericVector>(w)/sum(w);

  // current allocation ratio
  IntegerVector N(number_of_treatments);
  IntegerVector N1(number_of_treatments);

  // matrix of randomization probabilities
  NumericMatrix probability(number_of_subjects, number_of_treatments);

  // matrix of allocation proportions
  NumericMatrix allocation(number_of_subjects, number_of_treatments);

  // treatment assignments
  IntegerVector treatment(number_of_subjects);

  // responses
  NumericVector response(number_of_subjects);
  
  // imbalance
  NumericVector imbalance(number_of_subjects);

  // forcing index
  NumericVector fi(number_of_subjects);
  NumericVector forcing_index(number_of_subjects);

  // momentum of probability mass
  NumericVector mpm1_(number_of_subjects);
  NumericVector mpm1(number_of_subjects);
  
  NumericVector mpm2_(number_of_subjects);
  NumericVector mpm2(number_of_subjects);

  // List with current assignment
  List assignment;

  // randomization procedure
  std::function<List (int, IntegerVector)> rand_procedure = 
    set_rand_procedure(w, procedure, p);

  // response function
  std::function<NumericVector (IntegerVector)> response_function = 
    set_response_function(distribution, parameter);
  
  
  for(int j = 1; j <= number_of_subjects; j++) {
    N1 = N;
    assignment = rand_procedure(j, N);
    treatment[j-1] = assignment[0];
    probability.row(j-1) = as<NumericVector>(assignment[1]);

    N[treatment[j-1]-1] += 1;

    // allocation proportion
    allocation.row(j-1) = as<NumericVector>(N)/sum(N);

    // imbalance
    imbalance[j-1] = sqrt((float)sum(Rcpp::pow(as<NumericVector>(N) - j*rho, 2)))/j;
    mpm2_[j-1] = imbalance[j-1]*j;
    mpm2[j-1] = mean(mpm2_[seq(0, j-1)]);
    
    // forcing index
    fi[j-1] = sum(Rcpp::pow(probability.row(j-1)-rho, 2));
    forcing_index[j-1] = mean(fi[seq(0, j-1)]);

    // momentum of probability mass
    for (int k = 1; k <= number_of_treatments; k++) {
      N1[k-1] += 1;
      mpm1_[j-1] += probability.row(j-1)[k-1]*sqrt((float)sum(Rcpp::pow(as<NumericVector>(N1) - j*rho, 2)));
      N1[k-1] -= 1;
    }
    mpm1[j-1] = mean(mpm1_[seq(0, j-1)]);
    
  }  
  
  response = response_function(treatment);
  
  return List::create(_["treatment"] = treatment,
                      _["response"] = response,
                      _["Imb"] = imbalance,
                      _["FI"] = forcing_index,
                      _["MPM1"] = mpm1, 
                      _["MPM2"] = mpm2, 
                      _["probability"] = probability, 
                      _["allocation"] = allocation);
}


// simulation of restricted randomization procedure
//[[Rcpp::export(.simulate_restricted)]]
List simulate_restricted(int number_of_simulations, 
                         int number_of_subjects, 
                         IntegerVector w, 
                         std::string procedure, 
                         double p, 
                         std::string distribution, 
                         List parameter) {
  
  IntegerMatrix treatment(number_of_simulations, number_of_subjects);
  NumericMatrix response(number_of_simulations, number_of_subjects);
  NumericMatrix Imb(number_of_simulations, number_of_subjects);
  NumericMatrix FI(number_of_simulations, number_of_subjects);
  NumericMatrix MPM1(number_of_simulations, number_of_subjects);
  NumericMatrix MPM2(number_of_simulations, number_of_subjects);
  NumericMatrix Prob(number_of_subjects, w.size());
  List PropList(number_of_simulations);
  
  
  
  // run trial number_of_simulations times
  IntegerVector simulations = seq_len(number_of_simulations)-1;
  for_each(simulations.begin(), simulations.end(), 
           [&treatment, &response, &Imb, &FI, &MPM1, &MPM2, &Prob, &PropList,
            number_of_subjects, w, procedure, p, distribution, parameter](int &s){
              
              List trial = restricted(number_of_subjects, w, procedure, p, distribution, parameter);
              
              // treatment asignments within s-th simulation
              treatment.row(s) = as<IntegerVector>(trial[0]); 
              
              // responses within s-th simulation
              response.row(s) = as<NumericVector>(trial[1]); 
              
              // imbalances vs subject within s simulation
              Imb.row(s) = as<NumericVector>(trial[2]);
              
              // forcing indecies vs subject within s simulation
              FI.row(s) = as<NumericVector>(trial[3]);
              
              // momentum of probbaility mass vs subject within s-th simulation
              MPM1.row(s) = as<NumericVector>(trial[4]); // Alex' approach
              MPM2.row(s) = as<NumericVector>(trial[5]); // Olga's approach
              
              // unconditional allocation probability
              IntegerVector trt = seq_len(w.size());
              for_each(trt.begin(), trt.end(), [&Prob, trial](int k){
                Prob.column(k-1) = Prob.column(k-1)+as<NumericMatrix>(trial[6]).column(k-1);
              }); 
              
              
              // allocation proportions
              PropList[s] = as<NumericMatrix>(trial[7]);
              
            });
  
  // operating characteristics of the simulations
  NumericVector MI(number_of_subjects);
  NumericVector AFI(number_of_subjects);
  NumericVector AMPM1(number_of_subjects);
  NumericVector AMPM2(number_of_subjects);
  
  IntegerVector subject = seq_len(number_of_subjects);
  for_each(subject.begin(), subject.end(), 
           [Imb, FI, MPM1, MPM2, PropList, w, 
            &MI, &AFI, &AMPM1, &AMPM2](int &j) {
              MI[j-1] = max(Imb.column(j-1));          // maximum imbalance
              AFI[j-1] = mean(FI.column(j-1));         // average FI
              AMPM1[j-1] = mean(MPM1.column(j-1));     // average MPM1 (Alex' approach)
              AMPM2[j-1] = mean(MPM2.column(j-1));     // average MPM2 (Olga's approach)
            });
  
  
  return List::create(_["treatment"] = treatment,
                      _["response"] = response,
                      _["MI"] = MI,
                      _["AFI"] = AFI,
                      _["AMPM1"] = AMPM1, 
                      _["AMPM2"] = AMPM2,
                      _["Probability"] = Prob/number_of_simulations, 
                      _["Allocation"] = PropList);
  
}






