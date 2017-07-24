
overall_performance <- function(op) {
  imb_vs_fi <- simulations$imb_vs_fi
  crd <- subset(imb_vs_fi, name == "CRD")
  pbd1 <- subset(imb_vs_fi, name == "PBD (1.00)")
  UI <- (imb_vs_fi$max_imb - pbd1$max_imb)/(crd$max_imb - pbd1$max_imb)
  UR <- (imb_vs_fi$mean_fi - crd$mean_fi)/(pbd1$mean_fi - crd$mean_fi)
  wI <- wR <- 1
  G <- sqrt(((wI*UI)^2+(wR*UR)^2)/(wI^2+wR^2))
  ovp <- data.frame(design = imb_vs_fi$name, G = G, max_imb = imb_vs_fi$max_imb, mean_fi = imb_vs_fi$mean_fi)
  simulations$ovp <- cbind(rank = seq_len(nrow(ovp)), ovp[with(ovp, order(G)), ])
  names(simulations$ovp) <- c("Rank", "Design", "G", "Maximal Imbalance", "Mean FI")
  
}

# boxplots of allocation proportions
boxplot_allocation_proportion <- function(prop, save_plot = TRUE) {

  prop_melted <- reshape2::melt(prop,id.vars = c("nsim", "design", "w", "resp_type", "resp_model"))
  
  w <- unique(as.character(prop_melted$w))
  
  for (r in seq_along(w)) {
    ggplot(subset(prop_melted, w == w[r])) +
      geom_boxplot(aes(x = factor(design), y = value, fill=design))+
      facet_grid(variable ~ ., scales = "free_y")
    
    if (save_plot) {
      ggsave(sprintf("./plots/prop_boxplot_%s.pdf", w[r]), width=6, height=8)
    }
  }
  
}