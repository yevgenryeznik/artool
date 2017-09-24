#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(artool)

# Define UI for application
choices <- c("CRD", "PBD", "BUD", "MWUD", "DL", "DBCD", "MinQD", "MaxEnt")
shinyUI(fluidPage(

  # Application title
  titlePanel("Simulation of clinical trials targeting unequal allocation"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      textInput('w', 'Fixed allocation ratio (comma delimited)', "4, 3, 2, 1"),
      textInput('nsbj', 'Number of subjects', "200"),
      textInput('nsim', 'Number of simulations', "1000"),
      tags$b("Response: Normal"),
      wellPanel(
        textInput('resp_mean', 'Means (comma delimited)', "0, 0, 0, 0"),
        textInput('resp_sd', 'SDs (comma delimited)', "1, 1, 1, 1")
      ),
      textInput('alpha', 'Significance level', "0.05"),
      tags$b("Randomization procedures"),
      wellPanel(
      fluidRow(
        div(style="display:inline-block",checkboxInput("proc1_check", "CRD", 1))
      ),
      fluidRow(
        div(style="display:inline-block",checkboxInput("proc2_check", "", 1)),
        div(style="display:inline-block",selectInput("proc2_name", "", choices, "PBD",width='100px'),
            tags$head(tags$style(type="text/css", "#proc2_name {max-width: 150px}"))),
        div(style="display:inline-block",textInput("proc2_param", "", "1", width = '110px'),
            tags$head(tags$style(type="text/css", "#proc2_param {max-width: 150px}")))
      ),
      fluidRow(
        div(style="display:inline-block",checkboxInput("proc3_check", "", 1)),
        div(style="display:inline-block",selectInput("proc3_name", "", choices, "BUD",width='100px'),
            tags$head(tags$style(type="text/css", "#proc3_name {max-width: 150px}"))),
        div(style="display:inline-block",textInput("proc3_param", "", "2", width = '110px'),
            tags$head(tags$style(type="text/css", "#proc3_param {max-width: 150px}")))
      ),
      fluidRow(
        div(style="display:inline-block",checkboxInput("proc4_check", "", 1)),
        div(style="display:inline-block",selectInput("proc4_name", "", choices, "MWUD",width='100px'),
            tags$head(tags$style(type="text/css", "#proc4_name {max-width: 150px}"))),
        div(style="display:inline-block",textInput("proc4_param", "", "2", width = '110px'),
            tags$head(tags$style(type="text/css", "#proc4_param {max-width: 150px}")))
      ),
      fluidRow(
        div(style="display:inline-block",checkboxInput("proc5_check", "", 1)),
        div(style="display:inline-block",selectInput("proc5_name", "", choices, "DL",width='100px'),
            tags$head(tags$style(type="text/css", "#proc5_name {max-width: 150px}"))),
        div(style="display:inline-block",textInput("proc5_param", "", "2", width = '110px'),
            tags$head(tags$style(type="text/css", "#proc5_param {max-width: 150px}")))
      ),
      fluidRow(
        div(style="display:inline-block",checkboxInput("proc6_check", "", 1)),
        div(style="display:inline-block",selectInput("proc6_name", "", choices, "DBCD",width='100px'),
            tags$head(tags$style(type="text/css", "#proc6_name {max-width: 150px}"))),
        div(style="display:inline-block",textInput("proc6_param", "", "2", width = '110px'),
            tags$head(tags$style(type="text/css", "#proc6_param {max-width: 150px}")))
      ),
      fluidRow(
        div(style="display:inline-block",checkboxInput("proc7_check", "", 1, width = "10%")),
        div(style="display:inline-block",selectInput("proc7_name", "", choices, "MinQD",width='100px'),
            tags$head(tags$style(type="text/css", "#proc7_name {max-width: 150px}"))),
        div(style="display:inline-block",textInput("proc7_param", "", "0.5", width = '110px'),
            tags$head(tags$style(type="text/css", "#proc7_param {max-width: 150px}")))
      ),
      fluidRow(
        div(style="display:inline-block",checkboxInput("proc8_check", "", 1)),
        div(style="display:inline-block",selectInput("proc8_name", "", choices, "MaxEnt",width='100px'),
            tags$head(tags$style(type="text/css", "#proc8_name {max-width: 150px}"))),
        div(style="display:inline-block",textInput("proc8_param", "", "0.5", width = '110px'),
            tags$head(tags$style(type="text/css", "#proc8_param {max-width: 150px}")))
      )),
      actionButton('simulate', 'Simulate')
    ),

    # Show a plot of the generated distribution
    mainPanel(
       #plotOutput("prop_boxplot")
      tabsetPanel(
        tabPanel("Proportions Boxplots",
                 p(),
                 p("The plot below shows distributions of subjects' allocation proportions per treatment arm."),
                 p(),
                 plotOutput("allocation_boxplot", width = "100%", height = "1000px")),

        tabPanel("Unconditional Allocation Probability",
                 p(),
                 p("The plot below shows unconditional allocation probability per treatment for given designs."),
                 p(),
                 plotOutput("probability_plot", width = "100%", height = "1000px")),

        tabPanel("Imbalance/Forcing Index/Momentum of Probability Mass",
                 p(),
                 p("The plots below show maximum imbalance (MI), average FI (AFI), and average momentum of probability mass (AMPM) for given designs."),
                 p(),
                 plotOutput("op_plot", width = "100%", height = "1000px")),
        tabPanel("Overall Performance",
                 p(),
                 p(),
                 tableOutput("ovp_table")),
        tabPanel("Type I error",
                 p(),
                 p("The plots below show the type I error/power control through the simulations."),
                 p(),
                 plotOutput("tIerror_plot", width = "100%", height = "1000px")),
        tabPanel("Read Me (Description)",
                 p(),
                 includeHTML("description.html"))
      )

    )
  )
))
