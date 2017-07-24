# fixed allocation ratio
w <- rbind(
  c( 1,  1,  1,  1),
  c( 2,  1,  1,  2),
  c( 4,  3,  2,  1),
  c(37, 21, 21, 21)
)

# number of subjects
nsbj <- 200

# number of simulations
nsim <- 1000

# cohort size
cohort <- 1

# response
response <- list(
   list(resp_type = "flat", resp = "Normal", resp_params = rbind(c(0,1), c(0,1), c(0,1), c(0,1))),
   list(resp_type = "monotone", resp = "Normal", resp_params = rbind(c(0,1), c(0.1,1), c(0.2,1), c(0.25,1)))
)

# randomization procedure
design <- list(
  list(proc = "CRD",    proc_params = NA),
  list(proc = "PBD",    proc_params = 1),
  list(proc = "PBD",    proc_params = 2),
  list(proc = "PBD",    proc_params = 3),
  list(proc = "PBD",    proc_params = 4),
  list(proc = "PBD",    proc_params = 5),
  list(proc = "BUD",    proc_params = 2),
  list(proc = "BUD",    proc_params = 3),
  list(proc = "BUD",    proc_params = 4),
  list(proc = "BUD",    proc_params = 5),
  list(proc = "MWUD",   proc_params = 2),
  list(proc = "MWUD",   proc_params = 4),
  list(proc = "MWUD",   proc_params = 6),
  list(proc = "MWUD",   proc_params = 8),
  list(proc = "DL",     proc_params = 2),
  list(proc = "DL",     proc_params = 4),
  list(proc = "DL",     proc_params = 6),
  list(proc = "DL",     proc_params = 8),
  list(proc = "DBCD",   proc_params = 1),
  list(proc = "DBCD",   proc_params = 2),
  list(proc = "DBCD",   proc_params = 4),
  list(proc = "DBCD",   proc_params = 5),
  list(proc = "DBCD",   proc_params = 10),
  list(proc = "MaxEnt", proc_params = 0.05),
  list(proc = "MaxEnt", proc_params = 0.1),
  list(proc = "MaxEnt", proc_params = 0.25),
  list(proc = "MaxEnt", proc_params = 0.05),
  list(proc = "MaxEnt", proc_params = 1)
)

# significance level
alpha <- 0.05

# is a time drift included
time_drift <- list(
  list(model = "M1", value = FALSE), 
  list(model = "M2", value = TRUE)
)


# combine data to create scenario
scenario <- list()
for(w_row in seq_len(nrow(w))) {
  for(design_id in seq_along(design)) {
    for(response_id in seq_along(response)) {
      for(time_drift_id in seq_along(time_drift)) {
        scenario <- append(scenario,
                           list(list(
                             w = w[w_row,],
                             nsbj = nsbj,
                             nsim = nsim,
                             cohort = cohort,
                             resp_type = response[[response_id]]$resp_type,
                             resp = response[[response_id]]$resp,
                             resp_params = response[[response_id]]$resp_params,
                             proc = design[[design_id]]$proc,
                             proc_params = design[[design_id]]$proc_params,
                             alpha = alpha, 
                             time_drift_model = time_drift[[time_drift_id]]$model,
                             time_drift = time_drift[[time_drift_id]]$value
                           ))
        )
      }
    }
  }
}

# run simulations
simulated_data <- parallel::mcMap(
  function(sc) {
    rr_trial <- initialize_rr_trial(
      sc$w, 
      sc$nsbj, 
      sc$nsim, 
      sc$cohort,
      sc$resp, 
      sc$resp_params,
      sc$proc,
      sc$proc_params,
      sc$alpha, 
      sc$time_drift
    )
    
    return(simulate_rr_tiral(rr_trial))
  }, scenario, mc.cores = parallel::detectCores())


op <- do.call('rbind',
  Map(function(sc, simulated_data) {
  cbind(simulated_data$op, resp_type = sc$resp_type, resp_model = sc$time_drift_model)
}, scenario, simulated_data))


write.table(op, "./data/example1_op.csv2", dec = ".", sep = ";")


prop <- do.call('rbind',
              Map(function(sc, simulated_data) {
                cbind(simulated_data$prop, resp_type = sc$resp_type, resp_model = sc$time_drift_model)
              }, scenario, simulated_data))


write.table(prop, "./data/example1_prop.csv2", dec = ".", sep = ";")

