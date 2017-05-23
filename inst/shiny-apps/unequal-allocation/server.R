#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rartool)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  simulations <- reactiveValues(props = NULL,
                                pi_mean = NULL,
                                imb = NULL,
                                imb_vs_fi = NULL,
                                ovp = NULL, # overall performance
                                typeIerror = NULL)

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

    proc_name <- proc_name[proc_selected]
    proc_param <- proc_param[proc_selected]

    q <- qnorm(1-as.numeric(input$alpha)/2)

    w <- as.numeric(unlist(strsplit(input$w,",")))
    simulations$w <- w
    nsbj <- as.numeric(input$nsbj)
    nsim <- as.numeric(input$nsim)

    trials <- Map(function(proc, param){
      setup_trial(w, nsbj, nsim, 1, proc, param)
    }, proc_name, proc_param)

    sim_data <- Map(function(trial, proc, param){
                      trial$simulate()
                      col_names <- unlist(
                        lapply(seq_len(trial$numberOfTreatments),
                          function(k){
                            sprintf("treatment #%d", k)
                          })
                      )
                      design <- ifelse(is.na(param), sprintf("%s", proc), sprintf("%s (%4.2f)", proc, param))
                      nsbj <- seq_len(trial$numberOfSubjects)

                      # compute proportions
                      props <- as.data.frame(
                        do.call('rbind',
                                lapply(seq_len(nrow(trial$treatment)),
                                       function(r) {
                                         mapply(function(k){
                                                  sum(trial$treatment[r,]==k)/trial$numberOfSubjects
                                                }, seq_len(trial$numberOfTreatments))
                                        })
                              )
                      )
                      names(props) <- col_names
                      props <- cbind(name = design, props)
                      props <- reshape2::melt(props, id.vars = c("name"), vars = names(props))

                      # compute unconditional probabilities
                      pi_mean <- as.data.frame(
                        do.call('cbind',
                                lapply(seq_len(trial$numberOfTreatments),
                                       function(i) {
                                         apply(do.call('cbind',
                                                       lapply(trial$rand_probability,
                                                              function(prob){
                                                                prob[,i]
                                                              })
                                                       ), 1, mean)
                                       })
                                )
                      )
                      names(pi_mean) <- col_names
                      pi_mean <- cbind(name = design, nsbj = nsbj, pi_mean)
                      pi_mean <- reshape2::melt(pi_mean, id.vars = c("name", "nsbj"), vars = names(pi_mean))

                      # compute maximum imbalance
                      max_imb <- apply(trial$imbalance, 2, max)

                      # compute median and mean FI
                      median_fi <- apply(
                        mapply(
                          function(j) {
                            if (j == 1) { return(trial$forcing_idx[,1]) }
                            else { return(apply(trial$forcing_idx[,1:j], 1, mean)) }
                          }, seq_len(trial$numberOfSubjects)
                        ), 2, median)

                      mean_fi <- apply(
                        mapply(
                          function(j) {
                            if (j == 1) { return(trial$forcing_idx[,1]) }
                            else { return(apply(trial$forcing_idx[,1:j], 1, mean)) }
                          }, seq_len(trial$numberOfSubjects)
                        ), 2, mean)
                      imb <- data.frame(median_fi = median_fi, mean_fi,  max_imb = max_imb)

                      imb <- cbind(name = design, nsbj = nsbj, imb)
                      imb_vs_fi <- imb[nrow(imb), ]
                      imb <- reshape2::melt(imb, id.vars = c("name", "nsbj"), vars = names(imb))
                      imb$variable <- factor(imb$variable,
                                             levels = c("max_imb", "median_fi", "mean_fi"),
                                             labels = c("Maximal Imbalance", "Median Forcing Index", "Mean Forcing Index"))


                      # typeIerror
                      typeIerror_mean <- apply(trial$reject, 2, mean, na.rm = TRUE)
                      typeIerror_se <- apply(trial$reject, 2, sd, na.rm = TRUE)/sqrt(trial$numberOfSimulations)
                      typeIerror <- data.frame(name = design, nsbj = nsbj, mean = typeIerror_mean, error = q*typeIerror_se)

                      return(list(props = props,
                                  pi_mean = pi_mean,
                                  imb = imb,
                                  imb_vs_fi = imb_vs_fi,
                                  typeIerror = typeIerror
                                  ))
                    }, trials, proc_name, proc_param)

    simulations$props <- do.call('rbind', lapply(seq_along(sim_data), function(i){sim_data[[i]]$props}))
    simulations$pi_mean <- do.call('rbind', lapply(seq_along(sim_data), function(i){sim_data[[i]]$pi_mean}))
    simulations$imb <- do.call('rbind', lapply(seq_along(sim_data), function(i){sim_data[[i]]$imb}))
    simulations$imb_vs_fi <- do.call('rbind', lapply(seq_along(sim_data), function(i){sim_data[[i]]$imb_vs_fi}))
    simulations$typeIerror <- do.call('rbind', lapply(seq_along(sim_data), function(i){sim_data[[i]]$typeIerror}))

    # compute overall performance
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
  })

  output$prop_boxplot <- renderPlot({
    if (is.null(simulations$props)) return()
    ggplot(simulations$props)+
      geom_boxplot(aes(x = factor(name), y = value, fill = factor(name)))+
      #scale_y_continuous(limits = c(0, 1), breaks = c(0, prop$w/sum(prop$w), 1))+
      xlab("randomization procedure")+
      ylab("proportion")+
      facet_wrap(~ variable, ncol=1, scales = "free_y")+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 12),
            legend.title = element_blank(),
            legend.text = element_text(family = "Helvetica", face = "bold", size = 12))
  })

  output$pi_mean_plot <- renderPlot({
    if (is.null(simulations$pi_mean)) return()
    ggplot(simulations$pi_mean)+
      geom_line(aes(x = nsbj, y = value), size = 0.5)+
      geom_point(aes(x = nsbj, y = value), size = 1.5)+
      scale_y_continuous(limits = c(0, 1), breaks = simulations$w/sum(simulations$w))+
      xlab("number of subjects")+
      ylab("uncond. allocation probability")+
      facet_grid(name~variable)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 12))
  })

  output$imbalance_plot <- renderPlot({
    if (is.null(simulations$imb)) return()
    ggplot(simulations$imb)+
      geom_line(aes(x = nsbj, y = value, color = name), size = 0.5)+
      geom_point(aes(x = nsbj, y = value, color = name), size = 1.5)+
       #scale_y_continuous(limits = c(0, 1))+
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

  output$typeIerror_plot <- renderPlot({
    if (is.null(simulations$typeIerror)) return()
    ggplot(simulations$typeIerror)+
      geom_ribbon(aes(x=nsbj, ymin = mean-error, ymax = mean+error), fill = "lightblue")+
      geom_line(aes(x = nsbj, y = mean), size = 1.2)+
      #scale_y_continuous(limits = c(0, 1))+
      xlab("number of subjects")+
      ylab("")+
      facet_wrap(~ name, scales = "free_y", ncol = 1)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 12),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 12),
            legend.title = element_blank(),
            legend.text = element_text(family = "Helvetica", face = "bold", size = 12),
            strip.text = element_text(family = "Helvetica", face = "bold", size = 12))
  })

  output$performance_table <- renderTable({
    if (is.null(simulations$ovp)) return()
    simulations$ovp
    }, spacing = "m", digits = 3)
})
