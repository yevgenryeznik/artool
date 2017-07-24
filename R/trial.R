# setup and simulate clinical trial

Rcpp::loadModule("trial", TRUE)

# function to initialize a trial with Restricted Randomization procedure
initialize_rr_trial <- function(w, nsbj, nsim, 
                                resp, resp_params, 
                                proc, proc_params, 
                                alpha, time_drift = FALSE) {
  rr_trial <- new(TrialRR, w, nsbj, nsim,
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

  # summary of simulated data
  target <- paste0("w = (", paste(rr_trial$fixed_allocation_ratio, collapse = ':'),")")
  
  design <- rr_trial$rand_procedure
  
  procedure <- ifelse(design == "CRD", 
                      paste0(design),
                      paste0(design, " (", rr_trial$rand_procedure_params[1], ")"))
  
  response_type <- ifelse(min(rr_trial$response_distr_params[, 1]) == 
                          max(rr_trial$response_distr_params[, 1]), 
                          "flat dose-response",
                          "monotone dose-response")
  
  response_model <- ifelse(rr_trial$time_drift, 
                           "model M2 (with \"time drift\")",
                           "model M1 (without \"time drift\")")
  
  treatments <- map_chr(seq_len(rr_trial$number_of_treatments), 
                        function(k) { paste("treatment", k) })
  
  # 1) Maximal Imbalance vs. Number of Subjects: MI(n)
  MI <- rr_trial$imbalance %>%
    apply(2, max) %>%
    as_data_frame() %>%
    set_names(c("MI")) %>%
    bind_cols(subject = seq_len(nrow(.)))

  # 2) Average Forcing Index vs. Number of Subjects: AFI(n)
  AFI <- rr_trial$forcing_index %>%
    apply(2, mean) %>%
    as_data_frame() %>%
    set_names(c("AFI")) %>%
    bind_cols(subject = seq_len(nrow(.)))
  
  # 3) Variability of allocation proportions vs. Number of Subjects: ASD(n)
  ASD <- rr_trial$alloc_proportion %>%
    map(function(ap) {
      as_data_frame(ap) %>% 
        bind_cols(subject = seq_len(nrow(ap)))
      }) %>%
    bind_rows() %>% 
    gather(key, value, -subject) %>%
    group_by(subject, key) %>%
    summarise(VAR = var(value)) %>%
    group_by(subject) %>%
    summarise(ASDj = sum(VAR)) %>%
    mutate(ASD = sqrt(.$subject*.$ASDj)) %>%
    .[, c("subject", "ASD")]

  # 4) Allocation proportions in the end of a trial vs. Number of Simulations: AP(s)
  AP <- rr_trial$alloc_proportion %>%
    map(function(ap) {
      as_data_frame(matrix(ap[nrow(ap), ], ncol = length(ap[nrow(ap),])))
    }) %>%
    bind_rows() %>%
    set_names(treatments) %>%
    mutate(target = target, design = design, procedure = procedure, 
           response_type = response_type, response_model = response_model)
  
  # 5) Unconditional allocation probability vs. Number of Subjects: Pi(n)
  Pi <- rr_trial$rand_probability %>%
    map(function(prob) {
      as_data_frame(prob) %>% 
        bind_cols(subject = seq_len(nrow(prob)))
    }) %>%
    bind_rows() %>% 
    gather(key, value, -subject) %>%
    group_by(subject, key) %>%
    summarise(probability = mean(value)) %>%
    spread(key, probability) %>%
    set_names(c("subject", treatments)) %>%
    mutate(target = target, design = design, procedure = procedure, 
           response_type = response_type, response_model = response_model)
  
  # 6) Momentum of Probability Mass vs. Number of Subjects: MPM(n)
  MPM <- rr_trial$mpm %>%
    apply(2, mean) %>%
    as_data_frame() %>%
    set_names(c("MPM")) %>%
    bind_cols(subject = seq_len(nrow(.)))
  
  # 7) Type I error/Power +- SE vs. Number of Subjects: alpha(n), alphaSE(n)
  alpha <- t(rr_trial$reject) %>%
    as_data_frame() %>%
    bind_cols(subject = seq_len(nrow(.))) %>%
    gather(key, value, -subject) %>%
    group_by(subject) %>%
    summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), n = n(), nobs = sum(!is.na(value))) %>%
    filter(nobs > 1) %>%
    mutate(min = map_dbl(mean-1.96*sd/sqrt(n), ~max(., 0)), 
           max = map_dbl(mean+1.96*sd/sqrt(n), ~min(., 1)),
           target = target, design = design, procedure = procedure, 
           response_type = response_type, response_model = response_model) %>%
    .[, c("subject", "mean", "min", "max", "target", "design", "procedure", "response_type", "response_model")]
    
  
  op <- bind_cols(MI, AFI, ASD, MPM) %>% 
    .[, c("subject", "MI", "AFI", "ASD", "MPM")] %>%
    mutate(target = target, design = design, procedure = procedure, 
           response_type = response_type, response_model = response_model)
  
  return( list(op = op, AP = AP, Pi = Pi, alpha = alpha) )

}

