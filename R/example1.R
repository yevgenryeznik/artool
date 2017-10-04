#' simulation example for a paper in "Statistics in Medicine"
#' @param nsim number of simulations
#'
#' @export

run_example1 <- function(nsim) {
  # number of simulations
  nsim <- 10000
  
  # number of subjects
  nsbj <- 200
  
  # fixed allocation ratio
  w <- rbind(
    c( 1,  1,  1,  1),
    c( 2,  1,  1,  2),
    c( 4,  3,  2,  1),
    c(37, 21, 21, 21)
  )
  
  # number of treatments
  ntrt <- dim(w)[2] 

  # randomization procedure(s)
  proc <- c("CRD",  rep("PBD", 5), rep("BUD", 4), rep("MWUD", 4), rep("DL", 4), rep("DBCD", 5), rep("MaxEnt", 5))
  proc_param <- c(NA, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 4, 6, 8, 2, 4, 6, 8, 1, 2, 4, 5, 10, 0.05, 0.1, 0.25, 0.5, 1)
  
  # response distribution
  distr <- "normal"
  distr_param_flat <- list(mean = rep(0, 4), sd = rep(1, 4))
  # distr_param_monotone <- list(mean = c(0, 0.25, 0.5, 0.75), sd = rep(1, 4))  #example2
  # distr_param_monotone <- list(mean = c(0, 0.1, 0.25, 0.5), sd = rep(1, 4))  #example3
  # distr_param_monotone <- list(mean = c(0, 0.1, 0.25, 0.75), sd = rep(1, 4)) #example4
  # distr_param_monotone <- list(mean = c(0, 0.1, 0.25, 0.6), sd = rep(1, 4))  #example5
   distr_param_monotone <- list(mean = c(0, 0.1, 0.3, 0.7), sd = rep(1, 4))  #example5
  
  # significance level
  alpha <- 0.05
  
  # combine data to create scenario
  scenario <- invoke(rbind, 
                     map(seq_len(nrow(w)), 
                         ~ cbind(w_row = ., proc_id = seq_along(proc))
                     )
  )
  
  example1 <- list()
  for(sc in seq_len(nrow(scenario))){
    target <- as.numeric(w[scenario[sc,"w_row"],])
    proc_ <- proc[scenario[sc, "proc_id"]]
    proc_param_ <- proc_param[scenario[sc, "proc_id"]]
    procedure <- if_else(proc_ == "CRD", "CRD", paste0(proc_, " (", proc_param_, ")"))
    cat(paste0("scenario # ", sc, " is being simulated: target = ", 
               paste0("w = (", paste(target, collapse=","), ")"), 
               ", procedure = ", procedure, "\n"))
    
    trial <- simulate_rr(nsim, nsbj, target, proc_, proc_param_, distr, distr_param_flat)

    # Type I error/power 
    tIerror <- map(seq_len(nrow(trial$treatment)), ~ {
        treatment <- trial$treatment[., ]
        subject <- seq_along(treatment)
        
        # response1 is an original flat response (M1 model)
        response1 <- trial$response[., ] 
        flat_m1 <- map_dbl(subject, ~ artool:::.anova_test(treatment[1:.], response1[1:.], ntrt, alpha))
        
        # response2 is an original flat response with time drift added (M2 model)
        response2 <- response1 + subject/nsbj 
        flat_m2 <- map_dbl(subject, ~ artool:::.anova_test(treatment[1:.], response2[1:.], ntrt, alpha))
        
        # reponse3 is a monotone response (M1 model)
        response3 <- response(distr, distr_param_monotone, treatment)
        monotone_m1 <- map_dbl(subject, ~ artool:::.anova_test(treatment[1:.], response3[1:.], ntrt, alpha))
        
        # response4 is a monotone response response with time drift added (M2 model)
        response4 <- response3 + subject/nsbj
        monotone_m2 <- map_dbl(subject, ~ artool:::.anova_test(treatment[1:.], response4[1:.], ntrt, alpha))
        
        data_frame(subject, flat_m1, flat_m2, monotone_m1, monotone_m2)
        }) %>% bind_rows() %>%
      gather(variable, reject, -subject) %>%
      filter(!is.na(reject)) %>%
      group_by(variable, subject) %>%
      summarise(reject_mean = mean(reject), 
                reject_se = sd(reject)/sqrt(nsim))
    
    tIerror <- trial$op %>% 
      select(target, design, procedure, subject) %>% 
      inner_join(tIerror, by = "subject")
    
    example1[[sc]] <- list(
      op = trial$op,
      probability = trial$probability,
      allocation = trial$allocation, 
      tIerror = tIerror
    )
  } 
  save(example1, file = "./data/example1.Rda")
}


summary_example1 <- function(data_file = "./data/example1.Rda", summary_folder = "summary_example1_eps", ext = ".eps") {
  load(file = data_file)
  
  # colors to use in plots
  color <- c("red", "steelblue", "seagreen", "orange", "purple")

  # shapes to use in plots
  shape <- c(21, 22, 23, 24, 3)

  # consider all designs first
  selected_designs <- c("BUD", "CRD", "DBCD", "DL", "MaxEnt", "MWUD","PBD")
  selected <- "[all_designs]"
  
  # operational characteristics
  op <- example1 %>%
    map(~ {.$op}) %>%
      bind_rows() %>%
      gather(variable, value, -target, -design, -procedure, -subject) %>%
      filter(variable %in% c("MI", "AFI", "AMPM2", "ACMPM2")) 
  
  target_list <- unique(op$target)
  op$target <- factor(op$target, levels = target_list, ordered = TRUE)
    
  variable_list <- unique(op$variable)
  variable_labels <- c("maximum imbalance", "average forcing index", "average momentum of probability mass", "average cumulative momentum of probability mass")
  op$variable <- factor(op$variable, levels = variable_list, labels = variable_labels, ordered = TRUE)
  #op$variable1 <- factor(op$variable1, levels = variable_list, ordered = TRUE)
  
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
  
  # with the command below we create a set of plots:
  #   - MI vs subjcet 
  #   - AFI vs subject
  #   - AMPM2 vs subject
  #   - ACMPM2 vs subject
  # given target
  
  op_plots <- op %>% 
    group_by(target, variable, variable) %>% 
    do(plot = ggplot(data = ., aes(x = subject, y = value, group = procedure))+
         geom_point(aes(shape = procedure, fill = procedure, color = procedure), size = 2)+
         geom_line(aes(color = procedure), size = 0.1)+
         scale_shape_manual(values = shape_values)+
         scale_color_manual(values = color_values)+
         scale_fill_manual(values = color_values)+
         facet_wrap(~ design, scales = "free_y", ncol = 1)+
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
               legend.text  = element_text(family = "Helvetica", face = "bold", size = 10)))
  
  for (i in seq_len(nrow(op_plots))) {
    op_plots$plot[[i]]$labels$y <- as.character(op_plots$variable[i])
  }
  
  op_plot_files <- map2(op_plots$variable, op_plots$target, ~ paste0("./data/", summary_folder, "/", selected, " ", .x, "_vs_subject, ", .y, ext))
  map2(op_plot_files, op_plots$plot , ggsave, width = 16, height = 9, units = "in", dpi = 300)
  
  
  # overall performance
  ovp <- op %>% 
    spread(variable, value) %>%
    group_by(target, subject)
  
  mi <- ovp %>%
    filter(procedure %in% c("CRD", "MaxEnt (1)")) %>%
    select(target, procedure, subject, ACMPM = `average cumulative momentum of probability mass`) %>%
    spread(procedure, ACMPM) %>% 
    rename(ACMPM_CRD = CRD, ACMPM_MaxEnt1 = `MaxEnt (1)`) %>%
    group_by(target, subject)
  
  afi <- ovp %>%
    filter(procedure %in% c("CRD", "MaxEnt (1)")) %>%
    select(target, procedure, subject, AFI = `average forcing index`) %>%
    spread(procedure, AFI) %>% 
    rename(AFI_CRD = CRD, AFI_MaxEnt1 = `MaxEnt (1)`) %>%
    group_by(target, subject)
  
  wI <- wR <- 1
  ovpG <- ovp %>%
    select(target, 
           design, 
           procedure, 
           subject, 
           ACMPM = `average cumulative momentum of probability mass`,
           AFI = `average forcing index`) %>%
    inner_join(mi, by = c("target", "subject")) %>%
    inner_join(afi, by = c("target", "subject")) %>%
    mutate(UI = ACMPM/(ACMPM_CRD-ACMPM_MaxEnt1)-ACMPM_MaxEnt1/(ACMPM_CRD-ACMPM_MaxEnt1),
           UR = AFI/(AFI_MaxEnt1-AFI_CRD)-AFI_CRD/(AFI_MaxEnt1-AFI_CRD),
           G = sqrt(((wI*UI)^2 +(wR*UR)^2)/(wI^2 + wR^2))) %>%
    select(target, design, procedure, subject, G) %>%
    filter(!is.nan(G))
  
  # heat map plots
  ovpG_plots <- ovpG %>%
    group_by(target) %>%
    do(plot = ggplot(data = ., aes(subject, procedure)) +
         geom_tile(aes(fill=G))+
         scale_fill_gradientn(colors = rainbow(7))+
         scale_x_continuous(breaks = c(1, seq(25, 200, by = 25)))+
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
               legend.text  = element_text(family = "Helvetica", face = "bold", size = 10)))
  
  ovpG_plot_files <- map(ovpG_plots$target, ~ paste0("./data/", summary_folder, "/", selected, " ovpG_heatmap, ", ., ext))
  map2(ovpG_plot_files, ovpG_plots$plot , ggsave, width = 16, height = 9, units = "in", dpi = 300)
  
  
  # create *.xlxs file with ovp(s) for a final sample size
  wb <- loadWorkbook(paste0("./data/", summary_folder, "/ovp4.xlsx"), create = TRUE)
  ovpG_xlsx <- ovpG %>% 
    select(target, procedure, subject, G) %>%
    group_by(target) %>% 
    filter(subject == max(subject)) %>%
    arrange(target, G) %>%
    select(target, procedure, G) %>%
    group_by(target) %>%
    mutate(rank = seq_along(G))
  
  map(seq_along(target_list), function(i, target_, ovpG_xlsx, wb) {
    createSheet(wb, name = target_[i])
    ovpG_xlsx_sheet <- ovpG_xlsx %>% 
      ungroup() %>%
      filter(target == target_[i]) %>% 
      select(rank, procedure, G)

    writeWorksheet(wb, ovpG_xlsx_sheet, sheet = target_[i])
    
  }, target_list, ovpG_xlsx, wb)
  saveWorkbook(wb) 
  
  # ===============================================

  # AFI vs AMPM2 scatter plot
  op_sc <-  op %>% 
    filter(subject %in% c(25, 50, 100, 200) & 
           variable %in% c("average cumulative momentum of probability mass",
                           "average momentum of probability mass",
                           "average forcing index")) %>%
    spread(variable, value) %>%
    select(target, design, procedure, subject, 
           AFI = `average forcing index`,
           AMPM = `average momentum of probability mass`,
           ACMPM = `average cumulative momentum of probability mass`)
    #group_by(target) %>% 
  op_sc %>%
    ggplot(aes(x = ACMPM, y = AFI, group = procedure))+
         geom_point(aes(shape = design, fill = procedure, color = procedure), size = 4)+
         scale_shape_manual(values = design_values)+
         scale_color_manual(values = color_values)+
         scale_fill_manual(values = color_values)+
         facet_grid(target ~ subject, scales = "free")+
         #ggtitle("Average Forcing Index (AFI) vs. Average Momentum of Probability Mass (AMPM)")+
         xlab("average cumulative momentum of probability mass (ACMPM)")+
         ylab("average forcing index (AFI)")+
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
  
  ggsave(paste0("./data/", summary_folder, "/", selected, " afi_vs_acmpm2", ext), width = 16, height = 9, units = "in", dpi = 300)

  # ======================================
  
  # allocation proportion boxplots
  ap <- example1 %>%
    map(~ {.$allocation}) %>%
    bind_rows() %>%
    filter(design %in% selected_designs) %>%
    gather(treatment, proportion, -target, -design, -procedure) %>%
    group_by(target, design, procedure, treatment)

  map(seq_along(target_list), function(i, target_, ap) {
    w <- as.numeric(unlist(regmatches(target_[i], gregexpr("[[:digit:]]+", target_[i]))))
    ap %>%
      filter(target == target_[i]) %>%
      ggplot(aes(x = procedure, y = proportion))+
      geom_boxplot(aes(fill = procedure))+
      scale_y_continuous(breaks = w/sum(w))+
      #ggtitle(paste0("Allocation proportion: ", target_[i]))+
      ylab("allocation proportion")+
      facet_wrap(~ treatment, ncol = 1)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 6),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    
    ggsave(paste0("./data/", summary_folder, "/", selected, " ap_boxplots (", target_[i],")", ext), width = 16, height = 9, units = "in", dpi = 300)
    
  }, target_list, ap)
  
  # unconditional allocation probability plots
  Pi <- example1 %>%
    map(~ {.$probability}) %>%
    bind_rows() %>%
    gather(treatment, value, -target, -design, -procedure, -subject) %>%
    group_by(target, design, procedure, treatment, subject)
  

  design_ <- unique(Pi$design)
  
  map(seq_along(target_list), function(i, target_, design_, Pi) {
    sub_Pi <- Pi %>%
      filter(target == target_[i])
    map(seq_along(design_), function(j, target_i, design_, sub_Pi) {
      w <- as.numeric(unlist(regmatches(target_i, gregexpr("[[:digit:]]+", target_i))))
      sub_Pi %>%
        filter(design == design_[j]) %>%
        ggplot(aes(x = subject, y = value, group = treatment))+
        geom_point(aes(color = treatment), size = 1.25)+
        geom_line(aes(color = treatment), size = 0.5)+
        scale_y_continuous(limits = c(0, 1), breaks = w/sum(w))+
        scale_color_manual(values = c("darkblue", "darkgreen", "black", "red"), labels = c("treatement 1", "treatement 2", "treatement 3", "treatement 4"))+
        #ggtitle(paste0("Unconditional allocation probability: ", target_i, ", ", design_[j]))+
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
              legend.text  = element_text(family = "Helvetica", face = "bold", size = 16))
      
      
      ggsave(paste0("./data/", summary_folder, "/", "[arp_property]", " pi_plots (", target_i, ", ", design[j], ")", ext), width = 16, height = 9, units = "in", dpi = 300)
    }, target_[i], design_, sub_Pi)
  }, target_list, design_, Pi)

  # plot for the paper in SIM
  w <- c(4,3,2,1)
  Pi %>%
    filter(target == "w = (4,3,2,1)" & design %in% c("DBCD", "MaxEnt", "MWUD")) %>%
    ggplot(aes(x = subject, y = value, group = treatment))+
    geom_point(aes(color = treatment), size = 1.25)+
    geom_line(aes(color = treatment), size = 0.5)+
    scale_y_continuous(limits = c(0, 1), breaks = w/sum(w))+
    scale_color_manual(values = c("darkblue", "darkgreen", "black", "red"), labels = c("treatement 1", "treatement 2", "treatement 3", "treatement 4"))+
    #ggtitle(paste0("Unconditional allocation probability: ", target_i, ", ", design_[j]))+
    xlab("number of subjects")+
    ylab("unconditional allocation probability")+
    facet_wrap(~ procedure, ncol = 5)+
    theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
          axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
          axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
          title = element_text(family = "Helvetica", face = "bold", size = 14),
          strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
          strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text  = element_text(family = "Helvetica", face = "bold", size = 16))
  
  ggsave(paste0("./data/", summary_folder, "/", "[arp_property]", " pi_plots (w = (4,3,2,1), DBCD,MaxEnt, MWUD)", ext), width = 16, height = 9, units = "in", dpi = 300)
  # ======================================
  
  # cselect 4 designs to continuou:
  # BUD, CRD, DBCD, DL)
  
  selected_designs <- c("BUD", "CRD", "DBCD", "DL")
  selected <- "[4_designs]"
  
  # allocation proportion boxplots
  map(seq_along(target_list), function(i, target_, ap) {
    w <- as.numeric(unlist(regmatches(target_[i], gregexpr("[[:digit:]]+", target_[i]))))
    ap %>%
      filter(design %in% selected_designs) %>%
      filter(target == target_[i]) %>%
      ggplot(aes(x = procedure, y = proportion))+
      geom_boxplot(aes(fill = procedure))+
      scale_y_continuous(breaks = w/sum(w))+
      #ggtitle(paste0("Allocation proportion: ", target_[i]))+
      ylab("allocation proportion")+
      facet_wrap(~ treatment, ncol = 1)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
            title = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 14),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
    
    
    ggsave(paste0("./data/", summary_folder, "/", selected, " ap_boxplots (", target_[i],")", ext), width = 16, height = 9, units = "in", dpi = 300)
    
  }, target_list, ap)

    
  op1 <- example1 %>%
    map(~ {.$op}) %>%
    bind_rows() %>%
    filter(design %in% selected_designs) %>%
    select(-MI, -AMPM1, -ACMPM1, -AMPM2, -ACMPM2)
  op1$target <- factor(op1$target, levels = target_list, ordered = TRUE)
  
  design <- unique(op1$design)
  procedure <- unique(op1$procedure)
  
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
  
  design_values <- shape[1:4]
  names(design_values) <- design
  
  # ASD vs AFI for selected subjects
  sub_op11 <- op1 %>%
    filter(subject %in% c(50, 100, 200)) 
  sub_op12 <- op1 %>%
    filter(subject == 100 & target == "w = (1,1,1,1)") 
  
    ggplot()+
    geom_point(data = sub_op11, aes(x = AFI, y = ASD, shape = factor(subject), color = procedure), size = 3)+
    geom_line(data = sub_op11, aes(x = AFI, y = ASD, linetype = design, color = procedure), size = 1)+
    geom_point(data = sub_op12, aes(x = AFI, y = ASD, shape = factor(subject), color = procedure), size = 3)+
    geom_text(data = sub_op12, aes(x = AFI, y = ASD, label=procedure, color = procedure),
              size = 4, angle = 60 , family = "Helvetica", check_overlap = FALSE, nudge_y = 0.15)+
    labs(shape="# of subjects")+
    #ggtitle("\"Average Standard Deviation\" of Allocation Proportions (ASD) vs. Average Forcing Index (AFI)")+
    ylab("ASD")+
    facet_wrap(~ target, ncol = 1, scales = "free_y")+
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
  
   
  ggsave(paste0("./data/", summary_folder, "/", selected, " asd_vs_afi", ext), width = 16, height = 9, units = "in", dpi = 300)
  
  
  # type I error/power
  # tIerror <- example1 %>%
  #   map(~ .$tIerror) %>%
  #   bind_rows() %>%
  #   filter(design %in% selected_designs)
  # 
  # resp_type_ <- unique(alpha$response_type)
  # resp_model_ <- unique(alpha$response_model)
  # 
  # map(seq_along(target_), function(k, target_, resp_type_, resp_model_, aplha) {
  #   sub_alpha <- alpha %>%
  #     filter(target == target_[k])
  #   
  #   map(seq_along(resp_type_), function(i, target_k, resp_type_, resp_model_, sub_alpha) {
  #     sub_sub_alpha <- sub_alpha %>%
  #       filter(response_type == resp_type_[i])
  #     
  #     map(seq_along(resp_model_), function(j, target_k, resp_type_i, resp_model_, sub_sub_alpha){
  #       sub_sub_alpha %>%
  #         filter(response_model == resp_model_[j]) %>%
  #         ggplot(aes(x = subject, y = mean, ymin = min, ymax = max))+
  #           geom_ribbon(fill = "lightblue")+
  #           geom_line(size = 0.5)+
  #           xlab("number of subjects")+
  #           ylab("")+
  #           scale_x_continuous(breaks = c(1, seq(25, 200, by = 25)))+
  #           ggtitle(paste0(resp_type_i, ", ", resp_model_[j], ", ", target_k))+
  #           facet_wrap( ~ procedure, ncol = 5, scales = "free_y")+
  #           theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
  #                 axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
  #                 axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
  #                 axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
  #                 title = element_text(family = "Helvetica", face = "bold", size = 14),
  #                 strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
  #                 strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
  #                 legend.position = "none",
  #                 legend.title = element_blank(),
  #                 legend.text  = element_text(family = "Helvetica", face = "bold", size = 10))
  #       
  #       ggsave(paste0("./data/", summary_folder, "/", selected, "tIerror-power (", target_k, ", ", resp_type[i], ", ",  j,  ").pdf"), width = 16, height = 9, units = "in")    
  #     }, target_[k], resp_type_[i], resp_model_, sub_sub_alpha)
  #   }, target_[k], resp_type_, resp_model_, sub_alpha)
  # }, target_, resp_type_, resp_model_, alpha)
  
  # type I error/power table
  alpha <- example1 %>%
    map(~ {.$tIerror}) %>%
    bind_rows() %>%
    filter(subject %in% c(50, 100, 150, 200) ) %>%
    select(target, design, procedure, subject, 
           variable, alpha = reject_mean) %>%
    arrange(target, design, procedure, subject) %>%
    spread(variable, alpha)
  
  # create *.xlxs file with type I error/power values
  
  wb <- loadWorkbook(paste0("./data/", summary_folder, "/alpha.xlsx"), create = TRUE)
  map(seq_along(target_list), function(i, target_, alpha, wb) {
    createSheet(wb, name =gsub(":", " ", target_))
    sub_alpha1 <- alpha %>% 
      filter(target == target_[i] & subject == 50) %>%
      select(procedure, 
             `n = 50, Flat Response (No Drift)` = flat_m1, 
             `n = 50, Flat Response (With Drift)` = flat_m2, 
             `n = 50, Monotone Response (No Drift)` = monotone_m1, 
             `n = 50, Monotone Response (With Drift)` = monotone_m2)
    
    sub_alpha2 <- alpha %>% 
      filter(target == target_[i] & subject == 100) %>%
      select(procedure, 
             `n = 100, Flat Response (No Drift)` = flat_m1, 
             `n = 100, Flat Response (With Drift)` = flat_m2, 
             `n = 100, Monotone Response (No Drift)` = monotone_m1, 
             `n = 100, Monotone Response (With Drift)` = monotone_m2)

    sub_alpha3 <- alpha %>% 
      filter(target == target_[i] & subject == 150) %>%
      select(procedure, 
             `n = 150, Flat Response (No Drift)` = flat_m1, 
             `n = 150, Flat Response (With Drift)` = flat_m2, 
             `n = 150, Monotone Response (No Drift)` = monotone_m1, 
             `n = 150, Monotone Response (With Drift)` = monotone_m2)

    sub_alpha4 <- alpha %>% 
      filter(target == target_[i] & subject == 200) %>%
      select(procedure, 
             `n = 200, Flat Response (No Drift)` = flat_m1, 
             `n = 200, Flat Response (With Drift)` = flat_m2, 
             `n = 200, Monotone Response (No Drift)` = monotone_m1, 
             `n = 200, Monotone Response (With Drift)` = monotone_m2)
    sub_alpha <- inner_join(
      inner_join(sub_alpha1, sub_alpha2, by = c("procedure")),
      inner_join(sub_alpha3, sub_alpha4, by = c("procedure")),
      by = c("procedure")
    )
    
    writeWorksheet(wb, sub_alpha, sheet = target_[i])
    
  }, target_list, alpha, wb)
  saveWorkbook(wb) 
}


