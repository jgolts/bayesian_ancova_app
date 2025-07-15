library(tidyverse)
library(glue)
library(readxl)
library(parallel)
library(lmerTest)
library(parameters)
library(brms)
library(bayestestR)
library(ggdist)
library(ggfortify)
library(shiny)
library(bslib)

plotPosterior <- function(df, param) {
  
  df <- df %>%
    filter(Parameter == param)
  
  ggplot(data = df, aes(x = Value, y = Parameter)) +
    stat_halfeye(aes(fill = after_stat(level)), .width = c(.95, 1),
                 point_interval = mean_hdi) +
    scale_fill_brewer(na.translate = F, palette = "YlOrRd") +
    labs(title = "Posterior PDF with HDPI",
         x = "Estimate",
         y = NULL,
         fill = "HDPI")
}

ui <- page_fluid(
  titlePanel("Bayesian ANCOVA App for Continuous Outcome Data"),
  navset_tab( 
    nav_panel("Dataset Upload",
              pageWithSidebar(
                headerPanel(""),
                sidebarPanel(
                  fileInput("upload", "Please upload a dataset with baseline 
                  and follow-up measurements, treatment group variable, as well
                  as a patient identifier, in wide format
                            variable, in .xls or .xlsx format:",
                            buttonLabel = "Upload dataset",
                            accept = c(".xls", ".xlsx"),
                            multiple = F),
                  tableOutput("files"),
                  numericInput(inputId = "sheet",
                               label = "Select relevant sheet of Excel file",
                               value = 1, min = 1),
                  selectInput(inputId = "id", 
                              label = "Select identifier variable:",
                              choices = NULL),
                  selectInput(inputId = "group",
                              label = "Select treatment group variable:",
                              choices = NULL),
                  conditionalPanel(condition = "input.group != null &&
                                    input.group != ''",
                                   selectInput(inputId = "ref",
                                               label = "Which treatment group is the reference 
                                group?",
                                               choices = NULL),
                                   helpText("Note: This app currently performs analysis
                              only on trials with two groups")
                  ),
                  selectInput(inputId = "baseline", 
                              label = "Select baseline variable:",
                              choices = NULL),
                  selectInput(inputId = "fwup", 
                              label = "Select follow-up variable:",
                              choices = NULL),
                  helpText("Note: Adjustment for baseline covariates is not
                  currently supported. This app performs complete case analysis
                           only. Rows with missing baseline or follow-up
                           will be deleted")
                ),
                mainPanel(
                  titlePanel("Dataset preview"),
                  DT::DTOutput("contents")
                )
              )
    ),
    nav_panel("Bayesian Settings and Visualisation", 
              pageWithSidebar(
                headerPanel(""),
                sidebarPanel(
                  accordion(
                    accordion_panel(
                      "Set seed",
                      numericInput(inputId = "seed",
                                   label = "Seed for MCMC sampling",
                                   value = 157,
                                   min = 0)
                    ),
                    accordion_panel(
                      "Model intercept prior",
                      numericInput(inputId = "interceptMean",
                                   label = "Mean of normal prior distribution
                                  for the model intercept",
                                   value = 0),
                      numericInput(inputId = "interceptVariance",
                                   label = "Variance of normal prior
                                  distribution for model intercept",
                                   value = 10000,
                                   min = 0)
                    ),
                    accordion_panel(
                      "Baseline effect prior",
                      numericInput(inputId = "baselineMean",
                                   label = "Mean of normal prior distribution
                                   for model baseline effect",
                                   value = 0),
                      numericInput(inputId = "baselineVariance",
                                   label = "Variance of normal prior 
                                   distribution for model baseline effect",
                                   value = 10000,
                                   min = 0)
                    ),
                    accordion_panel(
                      "Treatment effect prior",
                      numericInput(inputId = "treatmentMean",
                                   label = "Mean of normal prior distribution
                                   for model treatment effect",
                                   value = 0),
                      numericInput(inputId = "treatmentVariance",
                                   label = "Variance of normal prior 
                                   distribution for model treatment effect",
                                   value = 10000,
                                   min = 0)
                    ),
                    accordion_panel(
                      "SD of residuals prior",
                      numericInput(inputId = "sigmaScale",
                                   label = "Scale of half-Cauchy prior
                                   distribution for standard deviation of
                                   outcome variable residuals",
                                   value = 25,
                                   min = 0)
                    ),
                    accordion_panel(
                      "MCMC iterations",
                      numericInput(inputId = "iter",
                                   label = "Set number of total MCMC
                                   iterations; half will be used for burn-in",
                                   value = 20000,
                                   min = 2)
                    ),
                    accordion_panel(
                      "Treatment effect direction",
                      selectInput(inputId = "nulleffect",
                                  label = "Is the treatment effect of interest 
                              above or below zero?",
                                  list(
                                    "TE > 0" = ">",
                                    "TE < 0" = "<")
                      )
                    ),
                    input_task_button("activate", "Fit models"),
                    div(),
                    helpText("Please allow a few minutes for Bayesian
                             model fitting")
                  )
                ),
                mainPanel(
                  accordion(
                    accordion_panel(
                      "Model intercept prior visualisation",
                      plotOutput("prior1plot")
                    ),
                    accordion_panel(
                      "Baseline effect prior visualisation",
                      plotOutput("prior2plot")
                    ),
                    accordion_panel(
                      "Treatment effect prior visualisation",
                      plotOutput("prior3plot")
                    ),
                    accordion_panel(
                      "SD of residuals prior visualisation",
                      plotOutput("prior4plot")
                    )
                  )
                )
              )),
    nav_panel("Bayesian Diagnostics",
              accordion(
                accordion_panel(
                  "MCMC trace plot",
                  plotOutput("traceplot"),
                  helpText("Please allow a moment for trace plot rendering")
                ),
                accordion_panel(
                  "Rhat",
                  DT::DTOutput("rhat")
                ),
                accordion_panel(
                  "MCMC trace plot (x2 iterations)",
                  plotOutput("traceplotx2"),
                  helpText("Please allow a moment for trace plot rendering")
                ),
                accordion_panel(
                  "Rhat (x2 iterations)",
                  DT::DTOutput("rhat2")
                ),
                accordion_panel(
                  "Autocorrelation plot",
                  plotOutput("acfplot")
                ),
                accordion_panel(
                  "Posterior draws histogram",
                  plotOutput("histplot")
                )
              )
    ),
    nav_panel("Results",
              accordion(
                accordion_panel(
                  "Frequentist results",
                  accordion(
                    accordion_panel(
                      "ANOVA-CHANGE",
                      DT::DTOutput("anova")
                    ),
                    accordion_panel(
                      "ANCOVA",
                      DT::DTOutput("ancova")
                    ),
                    accordion_panel(
                      "LMM",
                      DT::DTOutput("lmm")
                    )
                  )
                ),
                accordion_panel(
                  "Bayesian results",
                  accordion(
                    accordion_panel(
                      "Coefficient estimates and HDPIs",
                      DT::DTOutput("brm")
                    ),
                    accordion_panel(
                      "Posterior probability",
                      textOutput("postprob")
                    ),
                    accordion_panel(
                      "Posterior distributions of coefficents with HDPIs",
                      accordion(
                        accordion_panel(
                          "Model intercept",
                          plotOutput("interplot")
                        ),
                        accordion_panel(
                          "Baseline effect",
                          plotOutput("baseplot")
                        ),
                        accordion_panel(
                          "Treatment effect",
                          plotOutput("teplot")
                        ),
                        accordion_panel(
                          "SD of residuals",
                          plotOutput("sdplot")
                        )
                      )
                    )
                  )
                )
              ))
  ),
  id = "tab" 
)

server <- function(input, output, session) {
  data <- reactive({
    req(input$upload)
    validate(
      need(input$sheet <= length(excel_sheets(input$upload$datapath)),
           "Selected sheet not in dataset")
    )
    read_excel(input$upload$datapath, sheet = input$sheet)
  })
  
  output$contents <- DT::renderDT({
    DT::datatable(data())
  })  
  
  observe({
    req(input$upload)
    updateNumericInput(session, "sheet", max = length(
      excel_sheets(input$upload$datapath)
    ))
  })
  
  observe({
    data <- data()
    updateSelectInput(session, "id", choices = names(data))
    updateSelectInput(session, "group", choices = names(data))
    updateSelectInput(session, "baseline", choices = names(data))
    updateSelectInput(session, "fwup", choices = names(data))
  })
  
  var <- reactive({
    req(input$group)
    as.factor(data()[[input$group]])
  })
  
  observeEvent(input$group, {
    req(var())
    updateSelectInput(session, "ref", choices = levels(var()))
  })
  
  output$prior1plot <- renderPlot({
    mean <- input$interceptMean
    sd <- sqrt(input$interceptVariance)
    
    x_vals <- seq(mean - 4 * sd, mean + 4 * sd, length.out = 1000)
    df <- tibble(x = x_vals, y = dnorm(x_vals, mean = mean, sd = sd))
    
    ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) + 
      labs(x = "Model intercept",
           y = "Probability density")
  })
  
  output$prior2plot <- renderPlot({
    mean <- input$baselineMean
    sd <- sqrt(input$baselineVariance)
    
    x_vals <- seq(mean - 4 * sd, mean + 4 * sd, length.out = 1000)
    df <- tibble(x = x_vals, y = dnorm(x_vals, mean = mean, sd = sd))
    
    ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) + 
      labs(x = "Baseline effect",
           y = "Probability density")
  })
  
  output$prior3plot <- renderPlot({
    mean <- input$treatmentMean
    sd <- sqrt(input$treatmentVariance)
    
    x_vals <- seq(mean - 4 * sd, mean + 4 * sd, length.out = 1000)
    df <- tibble(x = x_vals, y = dnorm(x_vals, mean = mean, sd = sd))
    
    ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) + 
      labs(x = "Treatment effect",
           y = "Probability density")
  })
  
  output$prior4plot <- renderPlot({
    scale <- input$sigmaScale
    
    x_vals <- seq(0, 10 * scale, length.out = 1000)
    df <- tibble(x = x_vals, y = dcauchy(x_vals, 0, scale))
    
    ggplot(df, aes(x = x, y = y)) +
      geom_line(linewidth = 1.2) + 
      labs(x = "SD of model outcome",
           y = "Probability density")
  })
  
  sel_data <- reactive({
    
    validate(
      need(input$id, "Select an identifier variable"),
      need(input$group, "Select a treatment group variable"),
      need(input$baseline, "Select baseline variable"),
      need(input$fwup, "Select follow-up variable"),
      need(input$ref, "Please select a reference group."),
      need(length(unique(na.omit(data()[[input$group]]))) == 2, 
           "Only two-group datasets are currently supported."),
      need(is.numeric(data()[[input$baseline]]),
           "Baseline variable must be numeric."),
      need(is.numeric(data()[[input$fwup]]),
           "Follow-up variable must be numeric."),
      need(length(unique(c(input$id, input$group, input$baseline,
                           input$fwup))) == 4,
           "One of the dataset variables was selected more than once."),
      need(
        nrow(data()) == n_distinct(data()[[input$id]]),
        "Multiple rows per patient ID detected. Please ensure dataset is wide."
      )
      
    )
    
    df <- data() %>%
      select(input$id, input$group, input$baseline, input$fwup) %>% 
      drop_na() %>% 
      rename(old_group = !!sym(input$group)) %>% 
      mutate(change = .[[input$fwup]] - .[[input$baseline]]) %>% 
      mutate(!!input$group := case_when(
        old_group == input$ref ~ 0,
        .default = 1)) %>% 
      select(input$id, input$group, input$baseline, input$fwup, change)
  })
  
  longdata <- reactive({
    
    sd <- sel_data()
    validate(
      need(nrow(sd) > 0, "Not enough complete cases available for analysis.")
    )
    
    sd %>%
      select(input$id, input$group, input$baseline, input$fwup) %>% 
      pivot_longer(cols = c(input$baseline, input$fwup),
                   names_to = "time", values_to = "score") %>% 
      mutate(time = factor(time, levels = c(input$baseline, input$fwup)))
  })
  
  lmmformula <- reactive({
    
    req(sel_data(), longdata())
    
    as.formula(glue("score ~ time * `{input$group}` + (1| `{input$id}`)"))
  })
  
  anovaformula <- reactive({
    
    req(sel_data(), longdata())
    
    as.formula(glue("change ~ `{input$group}`"))
  })
  
  ancovaformula <- reactive({
    
    req(sel_data(), longdata())
    
    as.formula(glue("`{input$fwup}` ~ `{input$baseline}` + `{input$group}`"))
  })
  
  ffit <- eventReactive(
    input$activate,
    {
      req(lmmformula(), anovaformula(), ancovaformula())
      
      tryCatch({
        list(
          lmm = lmer(lmmformula(), data = longdata()),
          anova = lm(anovaformula(), data = sel_data()),
          ancova = lm(ancovaformula(), data = sel_data())
        )
      }, error = function(e) {
        showNotification("Error fitting frequentist models: check inputs",
                         type = "error")
        NULL
      })
    }
  )
  
  bfit <- eventReactive(
    input$activate,
    {
      req(ancovaformula())
      
      tryCatch({
        stanparams <- stanvar(input$interceptMean, "intercept_mean") +
          stanvar(sqrt(input$interceptVariance), "intercept_sd") +
          stanvar(input$sigmaScale, "sigma_scale")
        
        priors <- c(
          prior(normal(intercept_mean, intercept_sd), class = "Intercept"),
          prior_string(paste0("normal(", input$baselineMean, ", ",
                              sqrt(input$baselineVariance), ")"), 
                       class = "b", coef = paste0(input$baseline)),
          prior_string(paste0("normal(", input$treatmentMean, ", ",
                              sqrt(input$treatmentVariance), ")"), 
                       class = "b", coef = paste0(input$group)),
          prior(cauchy(0, sigma_scale), class = "sigma")
        )
        
        list(
          bayes = brm(ancovaformula(), data = sel_data(), silent = 1,
                      thin = 1, iter = input$iter,
                      warmup = input$iter / 2, prior = priors,
                      stanvars = stanparams, seed = input$seed,
                      cores = detectCores(), chains = 4),
          
          bayesx2 = brm(ancovaformula(), data = sel_data(), silent = 1,
                        thin = 1, iter = input$iter * 2,
                        warmup = input$iter, prior = priors,
                        stanvars = stanparams, seed = input$seed,
                        cores = detectCores(), chains = 4)
        )
      }, error = function(e) {
        showNotification("Error fitting Bayesian models: check priors and
                         inputs",
                         type = "error")
        NULL
      })
    })
  
  observeEvent(input$activate, {
    bfit()
    ffit()
  })
  
  draws <- reactive({
    
    as_draws_df(bfit()[["bayes"]])
  })
  
  output$traceplot <- renderPlot({
    
    mcmc_plot(bfit()[["bayes"]], type = "trace")
  })
  
  output$traceplotx2 <- renderPlot({
    
    mcmc_plot(bfit()[["bayesx2"]], type = "trace")
  })
  
  output$rhat <- DT::renderDT({
    
    req(bfit())
    
    sum <- summary(bfit()[["bayes"]])
    rbind(sum$fixed, sum$spec_pars) %>%
      tibble %>% 
      mutate(Coefficient = c("Intercept", input$baseline, input$group,
                             "sigma"),
             Rhat = round(Rhat, 2)) %>%
      select(Coefficient, Rhat)
  })
  
  output$rhat2 <- DT::renderDT({
    
    req(bfit())
    
    sum <- summary(bfit()[["bayesx2"]])
    rbind(sum$fixed, sum$spec_pars) %>%
      tibble %>% 
      mutate(Coefficient = c("Intercept", input$baseline, input$group,
                             "sigma"),
             Rhat = round(Rhat, 2)) %>%
      select(Coefficient, Rhat)
  })
  
  output$acfplot <- renderPlot({
    
    mcmc_plot(bfit()[["bayes"]], type = "acf")
  })
  
  output$histplot <- renderPlot({
    
    mcmc_plot(bfit()[["bayes"]], type = "hist")
  })
  
  output$anova <- DT::renderDT({
    
    req(ffit())
    
    sum_anova <- summary(ffit()[["anova"]])
    ci_anova <- confint(ffit()[["anova"]])
    bind_cols(sum_anova$coefficients, ci_anova) %>% 
      tibble() %>% 
      rename("P-Value" = "Pr(>|t|)") %>%
      mutate(Coefficient = c("Intercept", input$group),
             `95% CI` = paste0("(", round(`2.5 %`, 2), ", ",
                               round(`97.5 %`, 2), ")"),
             Estimate = round(Estimate, 2),
             `P-Value` = round(`P-Value`, 4)) %>%
      select(Coefficient, Estimate, `P-Value`, `95% CI`)
  })
  
  output$ancova <- DT::renderDT({
    
    req(ffit())
    
    sum_ancova <- summary(ffit()[["ancova"]])
    ci_ancova <- confint(ffit()[["ancova"]])
    bind_cols(sum_ancova$coefficients, ci_ancova) %>% 
      tibble() %>% 
      rename("P-Value" = "Pr(>|t|)") %>%
      mutate(Coefficient = c("Intercept", input$baseline, input$group),
             `95% CI` = paste0("(", round(`2.5 %`, 2), ", ",
                               round(`97.5 %`, 2), ")"),
             Estimate = round(Estimate, 2),
             `P-Value` = round(`P-Value`, 4)) %>% 
      select(Coefficient, Estimate, `P-Value`, `95% CI`)
  })
  
  output$lmm <- DT::renderDT({
    
    req(ffit())
    
    sum_lmm <- summary(ffit()[["lmm"]], ddf = "Kenward-Roger")
    ci_lmm <- model_parameters(ffit()[["lmm"]], ci_method = "kenward")
    bind_cols(sum_lmm$coefficients, drop_na(as_tibble(ci_lmm$CI_low)),
              drop_na(as_tibble(ci_lmm$CI_high))) %>% 
      rename("P-Value" = "Pr(>|t|)") %>%
      mutate(Coefficient = c("Intercept", "Time", input$group, 
                             paste0("Time:", input$group)),
             `95% CI` = paste0("(", round(`value...6`, 2), ", ",
                               round(`value...7`, 2), ")"),
             Estimate = round(Estimate, 2),
             `P-Value` = round(`P-Value`, 4)) %>% 
      select(Coefficient, Estimate, `P-Value`, `95% CI`)
  })
  
  output$brm <- DT::renderDT({
    
    req(bfit(), draws())
    
    sum_brm <- summary(bfit()[["bayes"]])
    sigma_ci_brm <- bayestestR::hdi(draws()$sigma) %>% 
      mutate(Parameter = "sigma") %>% 
      select(Parameter, CI, CI_low, CI_high)
    ci_brm <- bayestestR::hdi(bfit()[["bayes"]]) %>% 
      select(Parameter, CI, CI_low, CI_high) %>% 
      rbind(sigma_ci_brm)
    
    rbind(sum_brm$fixed, sum_brm$spec_pars) %>% 
      cbind(ci_brm) %>% 
      tibble() %>% 
      mutate(Coefficient = c("Intercept", input$baseline, input$group, "sigma"),
             Estimate = round(Estimate, 1),
             `95% HDPI` = paste0("(", round(CI_low, 1), ", ",
                                 round(CI_high, 1), ")")) %>% 
      select(Coefficient, Estimate, `95% HDPI`)
  })
  
  output$postprob <- renderPrint({
    
    req(draws())
    
    ops <- list(">" = `>`, "<" = `<`)
    prob <- mean(ops[[input$nulleffect]](draws()[,3], 0))
    if (input$nulleffect == "<"){
      print(paste0(round(prob*100 ,2), "% of posterior draws lie below the",
                   " line of null effect"))
    } else if (input$nulleffect == ">"){
      print(paste0(round(prob*100,2), "% of posterior draws lie above the line",
                   " of null effect"))
    }
  })
  
  draws_long <- reactive({
    
    req(draws())
    
    draws()[,1:4] %>% 
      pivot_longer(cols = everything(),
                   names_to = "Parameter", values_to = "Value")
  })
  
  output$interplot <- renderPlot({
    
    plotPosterior(draws_long(), "b_Intercept")
  })
  
  output$baseplot <- renderPlot({
    
    plotPosterior(draws_long(), paste0("b_", make.names(input$baseline)))
  })
  
  output$teplot <- renderPlot({
    
    plotPosterior(draws_long(), paste0("b_", make.names(input$group)))
  })
  
  output$sdplot <- renderPlot({
    
    plotPosterior(draws_long(), "sigma")
  })
  
}

shinyApp(ui = ui, server = server)