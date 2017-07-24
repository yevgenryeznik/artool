# setup and simulate clinical trial

Rcpp::loadModule("trial", TRUE)

# function to initialize a trial with Restricted Randomization procedure
initialize_rr_trial <- function(w, nsbj, nsim, cohort, 
                                resp, resp_params, 
                                proc, proc_params, 
                                alpha, time_drift = FALSE) {
  rr_trial <- new(TrialRR, w, nsbj, nsim, cohort, 
                  list(name = resp, parameters = resp_params), 
                  list(name = proc, parameters = proc_params),
                  alpha)
  rr_trial$set_time_drift(time_drift)
  
  return( rr_trial )
}


# function to simulate a trial with Restricted Randomization procedure

simulate_rr_tiral <- function(rr_trial) {
  # call method simulate()
  rr_trial$simulate()
  design <- ifelse(is.na(rr_trial$randomizationProcedureParams), 
                   rr_trial$randomizationProcedure, 
                   paste0(rr_trial$randomizationProcedure, " (", 
                          rr_trial$randomizationProcedureParams, ")"))
  
  # get all the operation characteristics
  # 1) Maximal Imbalance vs. Number of Subjects: MI(n)
  MI <- apply(rr_trial$imbalance, 2, max)
  
  # 2) Average Forcing Index vs. Number of Subjects: AFI(n)
  AFI <- apply(rr_trial$forcing_idx, 2, mean)
  
  # 3) Average Selection Bias vs. Number of Subjects: ASB(n)
  ASB <- apply(rr_trial$selection_bias, 2, mean)
  
  # 4) Variability of allocation proportions vs. Number of Subjects: ASD(n)
  ASD <- do.call('cbind',
            lapply(seq_len(rr_trial$numberOfTreatments),
                   function(treatment) {
                     apply(do.call('cbind',
                                   lapply(rr_trial$allocationProportion,
                                          function(prop){
                                            prop[,treatment]
                                          })
                     ), 1, sd)
                   })
    )
  ASD <- sqrt(seq_len(rr_trial$numberOfSubjects)*apply(ASD^2, 1, sum))

  # 5) Allocation proportions in the end of a trial vs. Number of Simulations: AP(s)
  AP <- as.data.frame(
    do.call('rbind',
            lapply(rr_trial$allocationProportion,
                   function(prop) {
                     prop[rr_trial$numberOfSubjects, ]
                   })
            )
  )
  
  col_names <- unlist(
    lapply(seq_len(rr_trial$numberOfTreatments),
           function(k){
             sprintf("treatment%d", k)
           })
  )
  
  names(AP) <- col_names
                   
  # 6) Unconditional allocation probability vs. Number of Subjects: Pi(n)
  Pi <- as.data.frame(
    do.call('cbind',
            lapply(seq_len(rr_trial$numberOfTreatments),
                   function(treatment) {
                     apply(do.call('cbind',
                                   lapply(rr_trial$randomizationProbability,
                                          function(prob){
                                            prob[,treatment]
                                          })
                     ), 1, mean)
                   })
    )
  )
  
  col_names <- unlist(
    lapply(seq_len(rr_trial$numberOfTreatments),
           function(k){
             sprintf("Pi%d", k)
           })
  )
  
  names(Pi) <- col_names
  
  # 7) Type I error/Power + 95% CI
  q <- qnorm(0.975)
  reject_mean <- apply(rr_trial$reject, 2, mean, na.rm = TRUE)
  reject_se <- apply(rr_trial$reject, 2, sd, na.rm = TRUE)/sqrt(rr_trial$numberOfSimulations)
  reject_low <- reject_mean-q*reject_se
  reject_upp <- reject_mean+q*reject_se
  
  reject <- data.frame(reject_low, reject_mean, reject_upp)
  
  op <- data.frame(design, 
                   w = paste0("w = (", paste(rr_trial$fixedAllocationRatio,collapse=" "), ")"),
                   nsbj = seq_len(rr_trial$numberOfSubjects),
                   MI, AFI, ASB, ASD, Pi, reject)
  
  prop <- data.frame(design, 
                     w = paste0("w = (", paste(rr_trial$fixedAllocationRatio,collapse=" "), ")"),
                     nsim = seq_len(rr_trial$numberOfSimulations),
                     AP)
                     
  return(list(op = op, prop = prop))
}