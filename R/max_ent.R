#library(tidyverse)

# bisection procedure to find a zero of a function
bisection <- function(fun, low, upp, tol = 1e-5) {
  if (abs(fun(low)) < tol) {
    return(low)
  }
  else if (abs(fun(upp)) < tol) {
    return(upp)
  }
  else {
    if (fun(low)*fun(upp) > 0) {
      stop("Function has the same signs at the interval bounds!")
    }
    else {
      mid <- 0.5*(low + upp)
      while (abs(fun(mid)) > tol && (upp-low) > tol) {
        low_cond <- as.numeric(fun(mid)*fun(low) > 0)
        low <- low_cond*mid + (1-low_cond)*low
        
        upp_cond <- as.numeric(fun(mid)*fun(upp) > 0)
        upp <- upp_cond*mid + (1-upp_cond)*upp
        
        mid <- 0.5*(low + upp)
      }
      return(mid)
    }
  }
}

# Newthon method
newthon <- function(fcn, fcn1, x0, tol = 1e-5){
  x1 <- x0
  while (abs(fcn(x1)) > tol) {
    
    x1 <- x1 - fcn(x1)/fcn1(x0)
    print(c(x1, fcn(x1), fcn1(x1), fcn(x1)/fcn1(x1)))
  }
  print(c("stop", paste(", x1 = ", x1)))
  return(x1)
}



# secant method
secant <- function(fcn, x1, x2, tol = 1e-5){
  x3 <- x2-fcn(x2)*(x2-x1)/(fcn(x2)-fcn(x1))
  while (abs(fcn(x3)) > tol && abs(x2-x1) > tol) {
    
    x1 <- x2
    x2 <- x3
    x3 <- x2-fcn(x2)*(x2-x1)/(fcn(x2)-fcn(x1))
    print(c(x1, x2, x3, fcn(x3)))
  }
  print(c("stop", paste(", x = ", x3)))
  return(x3)
}


#' Maximum Entropy Constrained Balance Design (MaxEnt(eta))
# @param w fixed allocation ration (target)
# @param eta paraemters of the procedure (0 <= eta <= 1)
# @param nsbj number of subjects to randomize
max_ent <- function(w, eta, nsbj) {
  # number of treatments
  ntrt <- length(w)

  # vector of treatment assignments
  trt <- rep(0, nsbj)

  # target allocation proportion
  rho <- w/sum(w)

  # allocation probabilities per subject
  prob <- matrix(0, nrow = nsbj, ncol = ntrt)

  # current allocation ration
  N <- rep(0, ntrt)

  # the hypothetical "lack of balance"
  B <- rep(0, ntrt)

  for(j in seq_len(nsbj)) {
    for (k in seq_len(ntrt)) {
      # hypothetical imbalance if a subject is assigned to a treatment k
      N1 <- rep(0, length(N))
      N1[k] <- 1
      N1  <- N1 + N
      B[k] <- max(abs(N1/j-rho))
    }
    if (var(B) <= 1e-16) {
      prob[j,] <- rho
    }
    else {
      # we have to find a zero of a one-variable (mu) function.
      # bisection method is used

      fcn <- function(mu) {
        min(B)*eta + (1-eta)*sum(rho*B) - sum(B*rho*exp(-mu*B))/sum(rho*exp(-mu*B))
      }
      
      fcn1 <- function(mu) {
        (sum(B*B*rho*exp(-mu*B))*sum(rho*exp(-mu*B)) - sum(B*rho*exp(-mu*B))^2)/
          sum(rho*exp(-mu*B))^2
      }
      # mu <- bisection(fcn, 0, 1000/max(B))
      # mu <- newthon(fcn, fcn1, 0)
      mu <- secant(fcn, 0, 1)
      prob[j,] <- rho*exp(-mu*B)/sum(rho*exp(-mu*B))
    }
    trt[j] <- sample(seq_len(ntrt), 1, TRUE, prob[j,])
    N[trt[j]] <- N[trt[j]] + 1
  }
  return(prob)
}


#' function to plot allocation probabilities vs number of subjects
# @param prob matrix of the allocation probabilities
visualize_prob <- function(prob) {
  prob_df <- bind_cols(sbj = seq_len(nrow(prob)), as_data_frame(prob))
  prob_df_names <- c("subject", map_chr(seq_len(ncol(prob)),
                                        function(k){
                                          paste("treatment", k)
                                        }))
  prob_df %>%
    set_names(prob_df_names) %>%
    gather(key, probability, -subject) %>%
    group_by(subject, probability) %>%
    ggplot(aes(x = subject, y = probability))+
      geom_point(size = 0.75) +
      geom_line(size = 0.5)+
      facet_wrap(~ key, ncol = 1, scales = "free_y")

}
