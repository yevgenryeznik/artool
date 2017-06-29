# setup and simulate clinical trial

Rcpp::loadModule("trial", TRUE)

initialize_rr_trial <- function(w, nsbj, nsim, cohort, 
                                resp, resp_params, 
                                proc, proc_params, 
                                alpha, time_drift = FALSE) {
  
  return( new(TrialRR, w, nsbj, nsim, cohort, 
              list(name = resp, parameters = resp_params), 
              list(name = proc, parameters = proc_params),
              alpha, time_drift ))
}
