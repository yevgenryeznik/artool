# modification of MaxEnt design in order to make it ARP  
max_ent2 <- function(number_of_subjects, w, distribution, parameter) {
  w2 <- rep(1, sum(w))
  rho <- w/sum(w)
  parameter1 <- map(seq_along(parameter), function(i) {
    unlist(map(seq_along(w), ~ rep(parameter[[i]][.], w[.]), w))
  })
  trial <- restricted(number_of_subjects, w2, "MaxEnt", 1, distribution, parameter1)
  
  # transform treatments' ids back 
  trial$treatment <- map_dbl(trial$treatment, function(trt, w) {
    seq_along(w)[
    map_lgl(seq_along(w), function(id, w) {
      vtrt <- cumsum(c(0, w))
      vtrt[id] < trt && trt <= vtrt[id+1]
    }, w)]
  }, w)
  
  # transform probabilities
  trial$probability <- invoke(rbind, map(seq_len(nrow(trial$probability)), ~ {
    map_dbl(seq_along(w), function(id, w){
      vtrt <- cumsum(c(0, w))
      sum(trial$probability[., seq(vtrt[id]+1, vtrt[id+1])])
    }, w)
  }))
  
  # transform proportions
  trial$proportion <- invoke(rbind, map(seq_len(nrow(trial$proportion)), ~ {
    map_dbl(seq_along(w), function(id, w){
      vtrt <- cumsum(c(0, w))
      sum(trial$proportion[., seq(vtrt[id]+1, vtrt[id+1])])
    }, w)
  }))
  
  # imbalance
  trial$imbalance <- map_dbl(seq_len(nrow(trial$proportion)), ~ {
    sqrt(sum((as.numeric(.*trial$proportion[.,])-.*rho)^2))/.
  })
  
  # forcing index
  fi <- map_dbl(seq_len(nrow(trial$probability)), ~ {
    sqrt(sum((as.numeric(.*trial$proportion[.,])-.*rho)^2))/.
  })
  trial$forcing_index <- map_dbl(seq_along(fi), ~ mean(fi[1:.]))
  
  # mpm
  trial$mpm2 <- trial$imbalance*seq_along(trial$imbalance)
  
  return(trial)
} 
  
