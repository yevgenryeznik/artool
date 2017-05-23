# setup and simulate clinical trial

Rcpp::loadModule("trial", TRUE)

setup_trial <- function(w, nsbj, nsim, cohort, proc, proc_param) {
  return( new(Trial, nsbj, nsim, cohort, proc, w, proc_param) )
}
