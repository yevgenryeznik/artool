run_example1 <- function() {
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
  
  # response
  response <- list(
    list(resp_type = "flat", resp = "Normal", resp_params = rbind(c(0,1), c(0,1), c(0,1), c(0,1))),
    list(resp_type = "monotone", resp = "Normal", resp_params = rbind(c(0,1), c(0.1,1), c(0.2,1), c(2,1)))
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
    list(proc = "MaxEnt", proc_params = 0.5),
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
  
  if (!file.exists("./data/sim_ex1.Rda")) {
    # run simulations
    sim_ex1 <- parallel::mcMap(
      function(sc) {
        rr_trial <- initialize_rr_trial(
          sc$w,
          sc$nsbj,
          sc$nsim,
          sc$resp,
          sc$resp_params,
          sc$proc,
          sc$proc_params,
          sc$alpha,
          sc$time_drift
        )
        
        return(simulate_rr_tiral(rr_trial))
      }, scenario, mc.cores = parallel::detectCores())
    
    save(sim_ex1, file = "./data/sim_ex1.Rda")
  }
}


summary_example1 <- function(all_designs = TRUE) {
  load("./data/sim_ex1.Rda")
  
  # fixed allocation ratio
  w <- rbind(
    c( 1,  1,  1,  1),
    c( 2,  1,  1,  2),
    c(37, 21, 21, 21),
    c( 4,  3,  2,  1)
  )
  # colors to use in plots
  color <- c("red", "steelblue", "seagreen", "orange", "purple")

  # shapes to use in plots
  shape <- c(21, 22, 23, 24, 3)

  # collect op data from simulations into a single data frame (all designs)
  selected_designs <- c("BUD", "CRD", "DBCD", "DL", "MaxEnt", "MWUD","PBD")
  selected <- "[all_designs]"
  
  op <- sim_ex1 %>%
    map(~ .$op) %>%
    bind_rows() %>%
    filter(design %in% selected_designs) %>%
    select(target, design, procedure, subject, MI, AFI, ASD, MPM) %>%
    group_by(target, design, procedure, subject) %>%
    summarise(MI = max(MI), AFI = mean(AFI), ASD = mean(ASD), MPM = mean(MPM))
  
  design <- unique(op$design)
  procedure <- unique(op$procedure)
  
  id <- unlist(map(design, 
      function(design_elem, procedure){
        designs <- map_chr(procedure, 
                function(procedure_elem){
                  strsplit(procedure_elem, " ")[[1]][1]
                })
        seq_len(sum(design_elem == designs))
      }, procedure))
  
  color_values <- color[id]
  names(color_values) <- procedure
  
  shape_values <- shape[id]
  names(shape_values) <- procedure
  
  design_values <- c(shape, 4, 8)
  names(design_values) <- design
  
  target_ <- unique(op$target)
  
  # visualize MI vs subject, AFI vs subject, MPM vs subject
  # MI vs AFI, MI vs MPM, AFI vs MPM
  map(seq_along(target_), function(i, target_, op){
      sub_op <- op %>% 
        filter(target == target_[i]) 
      
      # MPM vs subject
      sub_op %>%
        ggplot(aes(x = subject, y = MPM,  group = procedure))+
        geom_point(aes(shape = procedure, fill = procedure, color = procedure), size = 2)+
        geom_line(aes(color = procedure), size = 0.1)+
        scale_shape_manual(values = shape_values)+
        scale_color_manual(values = color_values)+
        scale_fill_manual(values = color_values)+
        facet_wrap(~ design, scales = "free_y", ncol = 1)+
        xlab("number of subjects")+
        ggtitle(paste0("Momentum of Probability Mass (MPM) vs. Number of Subjects: ", target_[i]))+
        theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
              axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
              axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
              axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
              title = element_text(family = "Helvetica", face = "bold", size = 14),
              strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
              strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
              legend.position = "right",
              legend.title = element_blank(),
              legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    ggsave(paste0("./data/plots/", selected, " mpm_vs_subject (", target_[i],").pdf"), width = 16, height = 9, units = "in")
    
    # MI vs subject
    sub_op %>%
      ggplot(aes(x = subject, y = MI,  group = procedure))+
      geom_point(aes(shape = procedure, fill = procedure, color = procedure), size = 2)+
      geom_line(aes(color = procedure), size = 0.1)+
      scale_shape_manual(values = shape_values)+
      scale_color_manual(values = color_values)+
      scale_fill_manual(values = color_values)+
      facet_wrap(~ design, scales = "free_y", ncol = 1)+
      xlab("number of subjects")+
      ggtitle(paste0("Maximum Imbalance (MI) vs. Number of Subjects: ", target_[i]))+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    ggsave(paste0("./data/plots/", selected, " mi_vs_subject (", target_[i],").pdf"), width = 16, height = 9, units = "in")

    # AFI vs subject
    sub_op %>%
      ggplot(aes(x = subject, y = AFI,  group = procedure))+
      geom_point(aes(shape = procedure, fill = procedure, color = procedure), size = 2)+
      geom_line(aes(color = procedure), size = 0.1)+
      scale_shape_manual(values = shape_values)+
      scale_color_manual(values = color_values)+
      scale_fill_manual(values = color_values)+
      facet_wrap(~ design, scales = "free_y", ncol = 1)+
      xlab("number of subjects")+
      ggtitle(paste0("Average Forcing Index (AFI) vs. Number of Subjects: ", target_[i]))+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    ggsave(paste0("./data/plots/", selected, " afi_vs_subject (", target_[i],").pdf"), width = 16, height = 9, units = "in")
    
    
    # MPM vs MI
    sub_op %>%
      ggplot(aes(x = MI, y = MPM,  group = procedure))+
      geom_point(aes(shape = procedure, fill = procedure, color = procedure), size = 2)+
      geom_line(aes(color = procedure), size = 0.1)+
      scale_shape_manual(values = shape_values)+
      scale_color_manual(values = color_values)+
      scale_fill_manual(values = color_values)+
      facet_wrap(~ design, scales = "free_y", ncol = 1)+
      ggtitle(paste0("Momentum of Probability Mass (MPM) vs. Maximum Imbalance: ", target_[i]))+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    ggsave(paste0("./data/plots/", selected, " mpm_vs_mi (", target_[i],").pdf"), width = 16, height = 9, units = "in")
    
    # AFI vs MI
    sub_op %>%
      ggplot(aes(x = MI, y = AFI,  group = procedure))+
      geom_point(aes(shape = procedure, fill = procedure, color = procedure), size = 2)+
      geom_line(aes(color = procedure), size = 0.1)+
      scale_shape_manual(values = shape_values)+
      scale_color_manual(values = color_values)+
      scale_fill_manual(values = color_values)+
      facet_wrap(~ design, scales = "free_y", ncol = 1)+
      ggtitle(paste0("Average Forcing Index (AFI) vs. Maximum Imbalance: ", target_[i]))+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    ggsave(paste0("./data/plots/", selected, " afi_vs_mi (", target_[i],").pdf"), width = 16, height = 9, units = "in")
    
    # AFI vs MPM
    sub_op %>%
      ggplot(aes(x = AFI, y = MPM,  group = procedure))+
      geom_point(aes(shape = procedure, fill = procedure, color = procedure), size = 2)+
      geom_line(aes(color = procedure), size = 0.1)+
      scale_shape_manual(values = shape_values)+
      scale_color_manual(values = color_values)+
      scale_fill_manual(values = color_values)+
      facet_wrap(~ design, scales = "free", ncol = 1)+
      ggtitle(paste0("Momentum of Probability Mass (MPM) vs. Average Forcing Index (AFI): ", target_[i]))+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    ggsave(paste0("./data/plots/", selected, " mpm_vs_afi (", target_[i],").pdf"), width = 16, height = 9, units = "in")
    
  }, target_, op) 
  
  # ===============================================
  
  # overall performance
  ovp <- op %>%
    group_by(target, subject)
  
  mi <- ovp %>%
    filter(procedure %in% c("CRD", "PBD (1)")) %>%
    .[, c("target", "procedure", "subject", "MI")] %>%
    spread(procedure, MI) %>% 
    rename(MI_CRD = CRD, MI_PBD1 = `PBD (1)`) %>%
    group_by(target, subject)
  
  afi <- ovp %>%
    filter(procedure %in% c("CRD", "PBD (1)")) %>%
    .[, c("target", "procedure", "subject", "AFI")] %>%
    spread(procedure, AFI) %>% 
    rename(AFI_CRD = CRD, AFI_PBD1 = `PBD (1)`) %>%
    group_by(target, subject)
  
  wI <- wR <- 1
  ovpG <- ovp %>%
    inner_join(mi, by = c("target", "subject")) %>%
    inner_join(afi, by = c("target", "subject")) %>%
    mutate(UI = MI/(MI_CRD-MI_PBD1)-MI_PBD1/(MI_CRD-MI_PBD1),
           UR = AFI/(AFI_PBD1-AFI_CRD)-AFI_CRD/(AFI_PBD1-AFI_CRD),
           G = sqrt(((wI*UI)^2 +(wR*UR)^2)/(wI^2 + wR^2))) %>%
    .[, c("target", "design", "procedure", "subject", "G")]
  
  target_ <- unique(ovpG$target)
  
  # plot ovp heatmaps G vs. (procedure, number of subjects)
  map(seq_along(target_), function(i, target_, ovpG){
    sub_ovpG <- ovpG %>% 
      filter(!is.nan(G)  & subject > 9 & target == target_[i]) %>% 
      .[, c("G", "target", "design", "procedure", "subject")] 
    
    sub_ovpG %>%
      ggplot(aes(subject, procedure)) +
      geom_tile(aes(fill=G))+
      scale_fill_gradientn(colors = rainbow(7))+
      scale_x_continuous(breaks = c(10, seq(25, 200, by = 25)))+
      ggtitle(paste0("Overall performance: ", target_))+
      xlab("number of subjects")+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    ggsave(paste0("./data/plots/", selected, " ovp_heatmap (", target_[i],").pdf"), width = 16, height = 9, units = "in")
    
  }, target_, ovpG)
  
  # create *.xlxs file with ovp(s) for a final sample size
  wb <- loadWorkbook("./data/plots/ovp4.xlsx", create = TRUE)
  map(seq_along(target_), function(i, target_, ovpG, wb) {
    createSheet(wb, name =gsub(":", " ", target_))
    sub_ovpG <- ovpG %>% 
      filter(subject == max(.$subject) & target == target_[i]) %>% 
      .[, c("G", "target", "subject", "procedure")] %>%
      arrange(G) %>%
      .[, c("procedure", "G")] %>%
      mutate(rank = seq_along(.$G)) 

    writeWorksheet(wb, sub_ovpG, sheet = gsub(":", " ", target_[i]))
    
  }, target_, ovpG, wb)
  saveWorkbook(wb) 
  
  # ===============================================

  # AFI vs MI scatter plot
  op %>%
    filter(subject %in% c(25, 50, 100, 200)) %>%
    ggplot(aes(x = MI, y = AFI,  group = procedure))+
    geom_point(aes(shape = design, fill = procedure, color = procedure), size = 4)+
    scale_shape_manual(values = design_values)+
    scale_color_manual(values = color_values)+
    scale_fill_manual(values = color_values)+
    facet_grid(target ~ subject, scales = "free")+
    ggtitle("Average Forcing Index (AFI) vs. Maximum Imbalance (MI)")+
    theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
          axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
          title = element_text(family = "Helvetica", face = "bold", size = 14),
          strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
          strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
  
  ggsave(paste0("./data/plots/", selected, " afi_vs_mi_selected_subjects", ".pdf"), width = 16, height = 9, units = "in")
  
  # MPM vs MI scatter plot
  op %>%
    filter(subject %in% c(25, 50, 100, 200)) %>%
    ggplot(aes(x = MI, y = MPM,  group = procedure))+
    geom_point(aes(shape = design, fill = procedure, color = procedure), size = 4)+
    scale_shape_manual(values = design_values)+
    scale_color_manual(values = color_values)+
    scale_fill_manual(values = color_values)+
    facet_grid(target ~ subject, scales = "free")+
    ggtitle("Momentum of Probability Mass (MPM) vs. Maximum Imbalance (MI)")+
    theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
          axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
          title = element_text(family = "Helvetica", face = "bold", size = 14),
          strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
          strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
  
  ggsave(paste0("./data/plots/", selected, " mpm_vs_mi_selected_subjects", ".pdf"), width = 16, height = 9, units = "in")


  # AFI vs MPM scatter plot
  op %>%
    filter(subject %in% c(25, 50, 100, 200)) %>%
    ggplot(aes(x = MPM, y = AFI,  group = procedure))+
    geom_point(aes(shape = design, fill = procedure, color = procedure), size = 4)+
    scale_shape_manual(values = design_values)+
    scale_color_manual(values = color_values)+
    scale_fill_manual(values = color_values)+
    facet_grid(target ~ subject, scales = "free")+
    ggtitle("Average Forcing Index (AFI) vs. Momentum of Probability Mass (MPM)")+
    theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
          axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
          title = element_text(family = "Helvetica", face = "bold", size = 14),
          strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
          strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
  
  ggsave(paste0("./data/plots/", selected, " afi_vs_mpm_selected_subjects", ".pdf"), width = 16, height = 9, units = "in")
  
  # ======================================
  
  # allocation proportion boxplots
  ap <- sim_ex1 %>%
    map(~ .$AP) %>%
    bind_rows() %>%
    filter(design %in% selected_designs) %>%
    gather(treatment, proportion, -target, -design, -procedure, -response_type, -response_model) %>%
    group_by(target, design, procedure, treatment)
  
  map(seq_along(target_), function(i, target_, ap) {
    ap %>%
      filter(target == target_[i]) %>%
      ggplot(aes(x = procedure, y = proportion))+
      geom_boxplot(aes(fill = procedure))+
      ggtitle(paste0("Allocation proportion: ", target_[i]))+
      ylab("allocation proportion")+
      facet_wrap(~ treatment, ncol = 1)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 6),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    
    ggsave(paste0("./data/plots/", selected, " ap_boxplots (", target_[i],").pdf"), width = 16, height = 9, units = "in")
    
  }, target_, ap)
  
  # unconditional allocation probability plots
  Pi <- sim_ex1 %>%
    map(~ .$Pi) %>%
    bind_rows() %>%
    filter(design %in% selected_designs) %>%
    gather(treatment, probability, -subject, -target, -design, -procedure, -response_type, -response_model) %>%
    group_by(target, design, procedure, subject, treatment) %>%
    summarise(probability = mean(probability))
  
  design_ <- unique(Pi$design)
  
  map(seq_along(target_), function(i, target_, design_, Pi) {
    sub_Pi <- Pi %>%
      filter(target == target_[i])
    map(seq_along(design_), function(j, target_i, design_, sub_Pi) {
      sub_Pi %>%
        filter(design == design_[j]) %>%
        ggplot(aes(x = subject, y = probability, group = treatment))+
        geom_point(aes(color = treatment), size = 0.75)+
        geom_line(aes(color = treatment), size = 0.25)+
        scale_y_continuous(limits = c(0, 1), breaks = w[i,]/sum(w[i,]))+
        ggtitle(paste0("Unconditional allocation probability: ", target_i, ", ", design_[j]))+
        xlab("number of subjects")+
        ylab("unconditional allocation probability")+
        facet_wrap(~ procedure, nrow = 1)+
        theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
              axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
              axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
              axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
              title = element_text(family = "Helvetica", face = "bold", size = 14),
              strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
              strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
      
      
      ggsave(paste0("./data/plots/", "[arp_property]", " pi_plots (", target_i, ", ", design[j], ").pdf"), width = 16, height = 9, units = "in")
    }, target_[i], design_, sub_Pi)
  }, target_, design_, Pi)

  # ======================================
  
  # collect op data from simulations into a single data frame (selcted designs:
  # BUD, CRD, DBCD, DL)
  
  selected_designs <- c("BUD", "CRD", "DBCD", "DL")
  selected <- "[4_designs]"
  
  # allocation proportion boxplots
  ap <- sim_ex1 %>%
    map(~ .$AP) %>%
    bind_rows() %>%
    filter(design %in% selected_designs) %>%
    gather(treatment, proportion, -target, -design, -procedure, -response_type, -response_model) %>%
    group_by(target, design, procedure, treatment)
  
  map(seq_along(target_), function(i, target_, ap) {
    ap %>%
      filter(target == target_[i]) %>%
      ggplot(aes(x = procedure, y = proportion))+
      geom_boxplot(aes(fill = procedure))+
      scale_y_continuous(breaks = w[i,]/sum(w[i,]))+
      ggtitle(paste0("Allocation proportion: ", target_[i]))+
      ylab("allocation proportion")+
      facet_wrap(~ treatment, ncol = 1)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    
    ggsave(paste0("./data/plots/", selected, " ap_boxplots (", target_[i],").pdf"), width = 16, height = 9, units = "in")
    
  }, target_, ap)

    
  op <- sim_ex1 %>%
    map(~ .$op) %>%
    bind_rows() %>%
    filter(design %in% selected_designs) %>%
    select(target, design, procedure, subject, MI, AFI, ASD, MPM) %>%
    group_by(target, design, procedure, subject) %>%
    summarise(MI = max(MI), AFI = mean(AFI), ASD = mean(ASD), MPM = mean(MPM))
  
  design <- unique(op$design)
  procedure <- unique(op$procedure)
  
  id <- unlist(map(design, 
                   function(design_elem, procedure){
                     designs <- map_chr(procedure, 
                                        function(procedure_elem){
                                          strsplit(procedure_elem, " ")[[1]][1]
                                        })
                     seq_len(sum(design_elem == designs))
                   }, procedure))
  
  color_values <- color[id]
  names(color_values) <- procedure
  
  shape_values <- shape[id]
  names(shape_values) <- procedure
  
  design_values <- c(shape, 4, 8)
  names(design_values) <- design
  
  target_ <- unique(op$target)
  
  # ASD vs AFI for selected subjects
  op %>%
    filter(subject %in% c(50, 100, 200)) %>%
    ggplot(aes(x = AFI, y = ASD))+
    geom_point(aes(shape = factor(subject), color = procedure), size = 3)+
    geom_line(aes(linetype = design, color = procedure), size = 1)+
    labs(shape="# of subjects")+
    ggtitle("\"Average Standard Deviation\" of Allocation Proportions (ASD) vs. Average Forcing Index (AFI)")+
    ylab("allocation proportion")+
    facet_wrap(~ target, ncol = 1)+
    theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
          axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
          title = element_text(family = "Helvetica", face = "bold", size = 14),
          strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
          strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
          legend.position = "right",
          #legend.title = element_blank(c),
          legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
  
   
  ggsave(paste0("./data/plots/", selected, " asd_vs_afi.pdf"), width = 16, height = 9, units = "in")
  
  
  # type I error/power
  alpha <- sim_ex1 %>%
    map(~ .$alpha) %>%
    bind_rows() %>%
    filter(design %in% selected_designs)

  resp_type_ <- unique(alpha$response_type)
  resp_model_ <- unique(alpha$response_model)
  
  map(seq_along(target_), function(k, target_, resp_type_, resp_model_, aplha) {
    sub_alpha <- alpha %>%
      filter(target == target_[k])
    
    map(seq_along(resp_type_), function(i, target_k, resp_type_, resp_model_, sub_alpha) {
      sub_sub_alpha <- sub_alpha %>%
        filter(response_type == resp_type_[i])
      
      map(seq_along(resp_model_), function(j, target_k, resp_type_i, resp_model_, sub_sub_alpha){
        sub_sub_alpha %>%
          filter(response_model == resp_model_[j]) %>%
          ggplot(aes(x = subject, y = mean, ymin = min, ymax = max))+
            geom_ribbon(fill = "lightblue")+
            geom_line(size = 0.5)+
            xlab("number of subjects")+
            ylab("")+
            scale_x_continuous(breaks = c(1, seq(25, 200, by = 25)))+
            ggtitle(paste0(resp_type_i, ", ", resp_model_[j], ", ", target_k))+
            facet_wrap( ~ procedure, ncol = 5, scales = "free_y")+
            theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
                  axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
                  axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
                  axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
                  title = element_text(family = "Helvetica", face = "bold", size = 14),
                  strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
                  strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
                  legend.position = "none",
                  legend.title = element_blank(),
                  legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
        
        ggsave(paste0("./data/plots/", selected, "tIerror-power (", target_k, ", ", resp_type[i], ", ",  j,  ").pdf"), width = 16, height = 9, units = "in")    
      }, target_[k], resp_type_[i], resp_model_, sub_sub_alpha)
    }, target_[k], resp_type_, resp_model_, sub_alpha)
  }, target_, resp_type_, resp_model_, alpha)
  
}


