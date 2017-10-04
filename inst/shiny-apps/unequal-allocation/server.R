#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(artool)

shinyServer(function(input, output) {

  sim_data <- reactiveValues(
    op = NULL,
    ovp = NULL, # overall performance
    probability = NULL,
    allocation = NULL,
    tIerror = NULL
  )
  
  observeEvent(input$simulate, {

    proc_selected <- c(input$proc1_check,
                       input$proc2_check,
                       input$proc3_check,
                       input$proc4_check,
                       input$proc5_check,
                       input$proc6_check,
                       input$proc7_check,
                       input$proc8_check)

    proc_name <- c("CRD",
                   input$proc2_name,
                   input$proc3_name,
                   input$proc4_name,
                   input$proc5_name,
                   input$proc6_name,
                   input$proc7_name,
                   input$proc8_name)

    proc_param <- as.numeric(c("NA",
                               input$proc2_param,
                               input$proc3_param,
                               input$proc4_param,
                               input$proc5_param,
                               input$proc6_param,
                               input$proc7_param,
                               input$proc8_param))

    proc <- proc_name[proc_selected]
    proc_param <- proc_param[proc_selected]
    
    
    distr <- "normal"
    distr_param <- list(
      mean = as.numeric(unlist(strsplit(input$resp_mean,","))),
      sd = as.numeric(unlist(strsplit(input$resp_sd,",")))
    )
      
    w <- as.numeric(unlist(strsplit(input$w,",")))
    sim_data$w <- w
    ntrt <- length(w)
    sim_data$shape_values <- c(3, 17, 8, 15, 16, 25, 18, 10)[proc_selected]
    sim_data$color_values <- c("red", "darkorange", "gold3", "darkgreen", 
                               "blue", "darkblue", "darkviolet", "black")[proc_selected]
    
    nsbj <- as.numeric(input$nsbj)
    nsim <- as.numeric(input$nsim)
    alpha <- as.numeric(input$alpha)
    
    trials <- Map(function(proc, proc_param){
      simulate_rr(nsim, nsbj, w, proc, proc_param, distr, distr_param)
    }, proc, proc_param)


    sim_data$op <- trials %>%
      map(~ {.$op}) %>%
      bind_rows() %>%
      select(target, design, procedure, subject, MI, AFI, AMPM = AMPM1) %>%
      gather(variable, value, -target, -design, -procedure, -subject)

          
    sim_data$allocation <- trials %>%
      map(~ {.$allocation}) %>%
      bind_rows() %>%
      gather(variable, value, -target, -design, -procedure)
    
    sim_data$probability <- trials %>%
      map(~ {.$probability}) %>%
      bind_rows() %>%
      gather(variable, value, -target, -design, -procedure, -subject)
    
    sim_data$tIerror <- trials %>%
      map(~ {
        mtreatment <- .$treatment
        mresponse <- .$response 
        procedure <- .$op$procedure
        tIerror <- map(seq_len(nrow(mtreatment)), ~ {
          treatment <- as.numeric(mtreatment[., ])
          response <- as.numeric(mresponse[., ])
          subject <- seq_along(treatment)
          reject <- map_dbl(subject, ~ artool:::.anova_test(treatment[1:.], response[1:.], ntrt, alpha))
          data_frame(procedure, subject, reject)
        })  %>% bind_rows() %>%
          filter(!is.na(reject)) %>%
          group_by(procedure, subject) %>%
          summarise(reject_mean = mean(reject),
                    reject_se = sd(reject)/sqrt(nsim))
      }) %>%
      bind_rows()
    
    
    ovp <- sim_data$op %>%
      filter(subject == max(.$subject) & variable %in% c("ACMPM2", "AFI")) %>%
      group_by(target, subject)
    
    mi <- ovp %>%
      filter(procedure %in% c("CRD", "MaxEnt (1)") & variable == "ACMPM2") %>%
      select(target, subject, procedure, value) %>% 
      spread(procedure, value) %>% 
      rename(ACMPM_CRD = CRD, ACMPM_MaxEnt1 = `MaxEnt (1)`) %>%
      group_by(target, subject)
    
    afi <- ovp %>%
      filter(procedure %in% c("CRD", "MaxEnt (1)") & variable == "AFI") %>%
      select(target, subject, procedure, value) %>% 
      spread(procedure, value) %>% 
      rename(AFI_CRD = CRD, AFI_MaxEnt1 = `MaxEnt (1)`) %>%
      group_by(target, subject)
    
    wI <- wR <- 1
    sim_data$ovp <- ovp %>%
      select(target, subject, procedure, variable, value) %>% 
      spread(variable, value) %>% 
      inner_join(mi, by = c("target", "subject")) %>%
      inner_join(afi, by = c("target", "subject")) %>%
      mutate(UI = ACMPM2/(ACMPM_CRD-ACMPM_MaxEnt1)-ACMPM_MaxEnt1/(ACMPM_CRD-ACMPM_MaxEnt1),
             UR = AFI/(AFI_MaxEnt1-AFI_CRD)-AFI_CRD/(AFI_MaxEnt1-AFI_CRD),
             G = sqrt(((wI*UI)^2 +(wR*UR)^2)/(wI^2 + wR^2))) %>%
      arrange(G) %>%
      mutate(Rank = seq_len(nrow(.)))%>%
      .[c("Rank", "procedure", "G")]
  })

  output$allocation_boxplot <- renderPlot({
    if (is.null(sim_data$allocation)) return()
    sim_data$allocation %>%
      ggplot()+
        geom_boxplot(aes(x = factor(procedure), y = value, fill = factor(procedure)))+
        xlab("randomization procedure")+
        ylab("allocation proportion")+
      facet_grid(variable ~ ., scales = "free_y")+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 12),
            legend.title = element_blank(),
            legend.text = element_text(family = "Helvetica", face = "bold", size = 12))
  })

  output$probability_plot <- renderPlot({
    if (is.null(sim_data$probability)) return()
    sim_data$probability %>%
    ggplot()+
      geom_line(aes(x = subject, y = value), size = 0.5)+
      geom_point(aes(x = subject, y = value), size = 1.5)+
      scale_y_continuous(limits = c(0, 1), breaks = sim_data$w/sum(sim_data$w))+
      xlab("number of subjects")+
      ylab("uncond. allocation probability")+
      facet_grid(procedure~variable)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 12))
  })

  output$op_plot <- renderPlot({
    if (is.null(sim_data$op)) return()
    sim_data$op %>%
    ggplot()+
      geom_line(aes(x = subject, y = value, color = procedure), size = 0.5)+
      geom_point(aes(x = subject, y = value, shape = procedure, color = procedure), size = 2.5)+
      scale_shape_manual(values = sim_data$shape_values)+
      scale_color_manual(values = sim_data$color_values)+
      xlab("number of subjects")+
      ylab("")+
      facet_wrap(~ variable, scales = "free_y", ncol = 1)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 12),
            legend.title = element_blank(),
            legend.text = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text = element_text(family = "Helvetica", face = "bold", size = 12))
  })

  output$tIerror_plot <- renderPlot({
    if (is.null(sim_data$tIerror)) return()
    sim_data$tIerror %>%
    ggplot()+
      geom_ribbon(aes(x=subject, ymin = reject_mean-reject_se, ymax = reject_mean+reject_se), fill = "lightblue")+
      geom_line(aes(x = subject, y = reject_mean), size = 1.2)+
      xlab("number of subjects")+
      ylab("")+
      facet_wrap(~ procedure, scales = "free_y", ncol = 1)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 12),
            legend.title = element_blank(),
            legend.text = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text = element_text(family = "Helvetica", face = "bold", size = 12))
  })

  output$ovp_table <- renderTable({
    if (is.null(sim_data$ovp)) return()
    sim_data$ovp
    }, spacing = "m", digits = 3)

})
  
  