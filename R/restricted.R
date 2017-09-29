#' Restricted randomization
#' 
#' Simulates a clinical trial one time with retricted randomizaion procedure
#' 
#' @author Yevgen Ryeznik (\email{yevgen.ryeznik@gmail.com}), Oleksandr Sverdlov
#' 
#' @param nsbj number of subjects to randomize.
#' @param w vector of fixed allocation ratio (integers with GCD = 1).
#' @param proc name of randomization procedure \code{c("CRD", "PBD", "BUD", "MWUD", 
#'    "DBCD", "DL", "MinQD", "MaxEnt")}.
#' @param proc_param tuning parameter of randomization procedure (NA for \code{"CRD"}).
#' @param distr distribution of response \code{c("binary", "uniform", "normal", 
#'    "exponential", "weibull", "loglogistic", "lognormal")}. See \link[artool]{response}
#'    for details.
#' @param distr_param parameter(s) of the response distribution (a list with numeric vector(s) with the
#'    same length as \code{w}). See \link[artool]{response}.
#'
#' @return a list with three items: 
#'    \itemize{
#'       \item \code{treatement} -- a vector of treatment assignments,
#'       \item \code{response} -- a vector of responses,
#'       \item \code{Imb} -- a vector of imbalances at each step,
#'       \item \code{FI} -- a vector of forcing indecies at each step,
#'       \item \code{MPM1} -- a vector of momentums of probability mass based on expected
#'             \emph{predicted} imbalance,
#'       \item \code{CMPM1} -- cumulative MPM1
#'       \item \code{MPM2} -- a vector of momentums of probability mass based on the expected
#'              imbalance,
#'       \item \code{CMPM2} -- cumulative MPM2
#'       \item \code{probability} -- a matrix with allocation probabilities at each step,
#'       \item \code{allocation} -- a matrix of allocation proportions at each step.
#'    }
#' @examples 
#'    rr(200, c(4, 3, 2, 1), "DL", 2, "normal", list(mean = rep(0, 4), sd = rep(1, 4)))
#'    rr(200, c(31, 13), "MaxEnt", 0.5, "normal", list(mean = c(0, 0.2), sd = c(1, 1))) 
#'       
#' @useDynLib artool   
#' 
#' @export
rr <- function(nsbj, w, proc, proc_param, distr, distr_param) {
  artool:::.restricted(nsbj, w, proc, proc_param, distr, distr_param)
} 


#' Restricted randomization
#' 
#' Makes multiple simulates of a clinical trial with retricted randomizaion procedures
#' 
#' @author Yevgen Ryeznik (\email{yevgen.ryeznik@gmail.com}), Oleksandr Sverdlov
#' 
#' @param nsim number of simulation runs.
#' @param nsbj number of subjects to randomize.
#' @param w vector of fixed allocation ratio (integers with GCD = 1).
#' @param proc name of randomization procedure \code{c("CRD", "PBD", "BUD", "MWUD", 
#'    "DBCD", "DL", "MinQD", "MaxEnt")}.
#' @param proc_param tuning parameter of randomization procedure (NA for \code{"CRD"}).
#' @param distr distribution of response \code{c("binary", "uniform", "normal", 
#'    "exponential", "weibull", "loglogistic", "lognormal")}. See \link[artool]{response}
#'    for details.
#' @param distr_param parameter(s) of the response distribution (a list with numeric vector(s) with the
#'    same length as \code{w}). See \link[artool]{response}.
#'
#' @return a list with three items: 
#'    \itemize{
#'       \item \code{treatment} -- a matrix of treatment assignments with \code{nsim} rows
#'       and \code{nsbj} columns, i.e. r-th row contains tretament assignments of r-th 
#'       simulation,
#'       \item \code{response} -- a matrix of responses with \code{nsim} rows
#'       and \code{nsbj} columns, i.e. r-th row contains responses of r-th 
#'       simulation, 
#'       \item{\code{op} -- a data frame which contains operating characteristics:
#'           \itemize{
#'              \item{\code{MI} -- maximum imbalance,}
#'              \item{\code{AFI} -- average forcing index,}
#'              \item{\code{AMPM1} -- average momentum of probability mass based on expected \emph{predicted} imbalance,}
#'              \item{\code{ACMPM1} -- average cumulative momentum of probability mass based on expected \emph{predicted} imbalance,}
#'              \item{\code{AMPM2} -- average momentum of probability mass based on the expected imbalance,}
#'              \item{\code{ACMPM2} -- average cumulative momentum of probability mass based on expected \emph{predicted} imbalance,}
#'              \item{\code{ASD} -- average standard deviation of allocation proportions,}
#'          }
#'       }
#'       \item{\code{probability} -- a data frame with unconditional allocation probabilities vs number of subjects
#'    for each treatment;}
#'       \item{\code{allocation} -- a data frame with allocation proportions in the end of a trial
#'    vs number of simulations for each treatment (could be usefull to make a boxplot of allocation proportion distribution).}
#'    } 
#' @examples 
#'    simulate_rr(1000, 200, c(4, 3, 2, 1), "DL", 2, "normal", list(mean = rep(0, 4), sd = rep(1, 4)))
#'    simulate_rr(1000, 200, c(31, 13), "MaxEnt", 0.5, "normal", list(mean = c(0, 0.2), sd = c(1, 1)))    
#' @useDynLib artool   
#' 
#' @export

simulate_rr <- function(nsim, nsbj, w, proc, proc_param, distr, distr_param) {
  target <- paste0("w = (", paste(w, collapse=","), ")")
  design <- proc
  procedure <- if_else(design == "CRD", "CRD", paste0(design, " (", proc_param, ")"))
  subject <- seq_len(nsbj)

  # simulate wiith Rcpp-based function "simulate_restricted"
  trial <- artool:::.simulate_restricted(nsim, nsbj, w, proc, proc_param, distr, distr_param)
  
  # extract alocation proportions
  allocation <- trial$Allocation %>%
    map(~ {
      bind_cols(
        data_frame(target, design, procedure, subject), 
        as_data_frame(.)
      ) %>%
        set_names(c("target", "design", "procedure", "subject", paste("treatment", seq_along(w))))
    }) %>%
    bind_rows()
  
  # ASD
  ASD <- allocation %>%
    gather(variable, value, -target, -design, -procedure, -subject) %>%
    group_by(variable, target, design, procedure, subject) %>%
    summarise(value_var = var(value)) %>%
    group_by(target, design, procedure, subject) %>%
    summarise(ASD = sum(value_var)) %>%
    mutate(ASD = sqrt(subject*ASD))
  
  # operational characteritsics
  op <- data_frame(target, 
                   design, 
                   procedure, 
                   subject, 
                   MI = trial$MI, 
                   AFI = trial$AFI, 
                   AMPM1 = trial$AMPM1,
                   ACMPM1 = trial$ACMPM1,
                   AMPM2 = trial$AMPM2, 
                   ACMPM2 = trial$ACMPM2) %>%
    inner_join(ASD, by = c("target", "design", "procedure", "subject"))
  
  # unconditional allocation probabilities
  probability <- bind_cols(
    data_frame(target, design, procedure, subject),
    as_data_frame(trial$Probability) %>% 
      set_names(paste("treatment", seq_along(w)))
  )
  
  # final allocation proportions
  allocation <- allocation %>% 
    filter(subject == max(subject)) %>%
    select(-subject)
  
  return ( list(
      treatment = trial$treatment, 
       response = trial$response,
             op = op, 
    probability = probability, 
     allocation = allocation
    ) )
}
