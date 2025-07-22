library(tidyverse)
library(glue)
library(readxl)
library(parallel)
library(lmerTest)
library(parameters)
library(brms)
library(bayestestR)
library(DescTools)
library(ggdist)
library(ggfortify)
library(shiny)
library(bslib)

plotPosterior <- function(df, param, size, type, ct) {
  
  dsize <- size/100
  
  interval <- switch(
    tolower(type),
    "hdi" = switch(tolower(ct),
                   "mean"   = mean_hdi,
                   "median" = median_hdi,
                   "mode"   = mode_hdi),
    "ci" = switch(tolower(ct),
                  "mean"   = mean_qi,
                  "median" = median_qi,
                  "mode"   = mode_qi)
  )
  
  df <- df %>%
    filter(Parameter == param)
  
  ggplot(data = df, aes(x = Value, y = Parameter)) +
    stat_halfeye(aes(fill = after_stat(level)), .width = c(dsize, 1),
                 point_interval = interval) +
    scale_fill_brewer(na.translate = F, palette = "YlOrRd") +
    labs(title = paste("Posterior PDF with", size, "%", toupper(type),
                       "and", ct),
         x = "Estimate",
         y = NULL,
         fill = toupper(type)
    )
}

plotPosteriorCDF <- function(df, param, size, type, ct) {
  
  dsize <- size/100
  
  interval <- switch(
    tolower(type),
    "hdi" = switch(tolower(ct),
                   "mean"   = mean_hdi,
                   "median" = median_hdi,
                   "mode"   = mode_hdi),
    "ci" = switch(tolower(ct),
                  "mean"   = mean_qi,
                  "median" = median_qi,
                  "mode"   = mode_qi)
  )
  
  df <- df %>%
    filter(Parameter == param)
  
  ggplot(data = df, aes(x = Parameter, y = Value)) +
    stat_ccdfinterval(aes(fill = after_stat(level)), .width = c(dsize, 1),
                 point_interval = interval, justification = 1) +
    scale_fill_brewer(na.translate = F, palette = "YlOrRd") +
    labs(title = paste("Posterior CDF with", size, "%", toupper(type),
                       "and", ct),
         x = "Estimate",
         y = NULL,
         fill = toupper(type)
    )
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
                                               label = "Which treatment group 
                                               is the reference 
                                group?",
                                               choices = NULL),
                                   helpText("Note: This app currently performs 
                                   analysis
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
                                   value = 10,
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
                              below or above zero?",
                                  list(
                                    "TE < 0" = "<",
                                    "TE > 0" = ">")
                      )
                    ),
                    accordion_panel(
                      "Central tendency measure for plots",
                      selectInput(inputId = "plotCT",
                                  label = "Which measure of central tendency 
                                  should be plotted for posterior 
                                  visualisation?",
                                  list(
                                    "Mean" = "mean",
                                    "Median" = "median",
                                    "Mode" = "mode"
                                    )
                                  )
                    ),
                    accordion_panel(
                      "Intervals for plots",
                      selectInput(inputId = "plotCI",
                                  label = "Credible intervals or HDIs for 
                                  posterior visualisation?",
                                  list(
                                    "Credible Intervals" = "ci",
                                    "HDIs" = "hdi"
                                    )
                                  )
                    ),
                    accordion_panel(
                      "Interval size",
                      numericInput(inputId = "sizeCI",
                                   label = "What width should be used for 
                                   CIs/HDIs?",
                                   value = 95,
                                   min = 50,
                                   max = 99
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
                  helpText("Please allow a moment for trace plot rendering"),
                  downloadButton(outputId = "dl_traceplot", label = "Download")
                ),
                accordion_panel(
                  "Rhat",
                  DT::DTOutput("rhat")
                ),
                accordion_panel(
                  "MCMC trace plot (x2 iterations)",
                  plotOutput("traceplotx2"),
                  helpText("Please allow a moment for trace plot rendering"),
                  downloadButton(outputId = "dl_traceplotx2",
                                 label = "Download")
                ),
                accordion_panel(
                  "Rhat (x2 iterations)",
                  DT::DTOutput("rhat2")
                ),
                accordion_panel(
                  "Autocorrelation plot",
                  plotOutput("acfplot"),
                  downloadButton(outputId = "dl_acfplot", label = "Download")
                ),
                accordion_panel(
                  "Posterior draws histogram",
                  plotOutput("histplot"),
                  downloadButton(outputId = "dl_histplot", label = "Download")
                )
              )
    ),
    nav_panel("Results",
              accordion(
                accordion_panel(
                  "Frequentist results",
                  DT::DTOutput("ancova"),
                  downloadButton(outputId = "dl_ancova", label = "Download")
                ),
                accordion_panel(
                  "Bayesian results",
                  accordion(
                    accordion_panel(
                      "Coefficient estimates and intervals",
                      DT::DTOutput("brm"),
                      downloadButton(outputId = "dl_brm", label = "Download")
                    ),
                    accordion_panel(
                      "Posterior probability",
                      textOutput("postprob")
                    ),
                    accordion_panel(
                      "Plots of posterior distributions of coefficents with 
                      selected intervals",
                      accordion(
                        accordion_panel(
                          "PDFs",
                          accordion(
                            accordion_panel(
                            "Model intercept",
                            plotOutput("interplotpdf"),
                            downloadButton(outputId = "dl_interplotpdf",
                                           label = "Download")
                            ),
                            accordion_panel(
                            "Baseline effect",
                            plotOutput("baseplotpdf"),
                            downloadButton(outputId = "dl_baseplotpdf",
                                           label = "Download")
                            ),
                            accordion_panel(
                            "Treatment effect",
                            plotOutput("teplotpdf"),
                            downloadButton(outputId = "dl_teplotpdf",
                                           label = "Download")
                            ),
                            accordion_panel(
                            "SD of residuals",
                            plotOutput("sdplotpdf"),
                            downloadButton(outputId = "dl_sdplotpdf",
                                           label = "Download")
                            )
                            )
                          ),
                        accordion_panel(
                          "CCDFs",
                          accordion(
                            accordion_panel(
                              "Model intercept",
                              plotOutput("interplotcdf"),
                              downloadButton(outputId = "dl_interplotcdf",
                                             label = "Download")
                            ),
                            accordion_panel(
                              "Baseline effect",
                              plotOutput("baseplotcdf"),
                              downloadButton(outputId = "dl_baseplotcdf",
                                             label = "Download")
                            ),
                            accordion_panel(
                              "Treatment effect",
                              plotOutput("teplotcdf"),
                              downloadButton(outputId = "dl_teplotcdf",
                                             label = "Download")
                            ),
                            accordion_panel(
                              "SD of residuals",
                              plotOutput("sdplotcdf"),
                              downloadButton(outputId = "dl_sdplotcdf",
                                             label = "Download")
                            )
                          )
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

  ancovaformula <- reactive({
    
    req(sel_data())
    
    sd <- sel_data()
    validate(
      need(nrow(sd) > 0, "Not enough complete cases available for analysis.")
    )
    
    as.formula(glue("`{input$fwup}` ~ `{input$baseline}` + `{input$group}`"))
  })
  
  ffit <- eventReactive(
    input$activate,
    {
      req(ancovaformula())
      
      tryCatch({
        list(
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
                      cores = detectCores(), chains = min(4, detectCores())),
          
          bayesx2 = brm(ancovaformula(), data = sel_data(), silent = 1,
                        thin = 1, iter = input$iter * 2,
                        warmup = input$iter, prior = priors,
                        stanvars = stanparams, seed = input$seed,
                        cores = detectCores(), chains = min(4, detectCores()))
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
    bind_rows(sum$fixed, sum$spec_pars) %>%
      tibble %>% 
      mutate(Coefficient = c("Intercept", input$baseline, input$group,
                             "sigma"),
             Rhat = round(Rhat, 2)) %>%
      select(Coefficient, Rhat)
  })
  
  output$rhat2 <- DT::renderDT({
    
    req(bfit())
    
    sum <- summary(bfit()[["bayesx2"]])
    bind_rows(sum$fixed, sum$spec_pars) %>%
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
  
  output$dl_traceplot <- downloadHandler(
    filename = function() { paste0("traceplot_", Sys.Date(), ".jpg") },
    content = function(file) {
      ggsave(file, plot = mcmc_plot(bfit()[["bayes"]], type = "trace"),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_traceplotx2 <- downloadHandler(
    filename = function() { paste0("traceplot_x2_", Sys.Date(), ".jpg") },
    content = function(file) {
      ggsave(file, plot = mcmc_plot(bfit()[["bayesx2"]], type = "trace"),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_acfplot <- downloadHandler(
    filename = function() { paste0("acfplot_", Sys.Date(), ".jpg") },
    content = function(file) {
      ggsave(file, plot = mcmc_plot(bfit()[["bayes"]], type = "acf"),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_histplot <- downloadHandler(
    filename = function() { paste0("histplot_", Sys.Date(), ".jpg") },
    content = function(file) {
      ggsave(file, plot = mcmc_plot(bfit()[["bayes"]], type = "hist"),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )

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
  
  output$dl_ancova <- downloadHandler(
    filename = function() { paste0("frequentist_results_", Sys.Date(),
                                   ".csv") },
    content = function(file) {
      req(ffit())
      sum_ancova <- summary(ffit()[["ancova"]])
      ci_ancova <- confint(ffit()[["ancova"]])
      table_out <- bind_cols(sum_ancova$coefficients, ci_ancova) %>%
        tibble() %>%
        rename("P-Value" = "Pr(>|t|)") %>%
        mutate(Coefficient = c("Intercept", input$baseline, input$group),
               `95% CI` = paste0("(", round(`2.5 %`, 2), ", ",
                                 round(`97.5 %`, 2), ")"),
               Estimate = round(Estimate, 2),
               `P-Value` = round(`P-Value`, 4)) %>%
        select(Coefficient, Estimate, `P-Value`, `95% CI`)
      write.csv(table_out, file, row.names = FALSE)
    }
  )
  
  output$brm <- DT::renderDT({
    
    req(bfit(), draws())
    
    draws <- as_tibble(draws())
    
    size <- input$sizeCI / 100
    size_label <- paste0(input$sizeCI, "%")
    
    sd_brm <- draws %>% 
      select(b_Intercept:sigma) %>% 
      sapply(sd) %>% 
      as_tibble() %>% 
      rename("SD" = "value")
    
    mode_brm <- draws %>% 
      select(b_Intercept:sigma) %>% 
      sapply(DescTools::Mode) %>% 
      as_tibble()
    
    median_brm <- draws %>% 
      select(b_Intercept:sigma) %>% 
      sapply(median) %>% 
      as_tibble() %>% 
      rename("Median" = "value")
    
    sum_brm <- summary(bfit()[["bayes"]])
    
    hdi_brm <- bayestestR::hdi(bfit()[["bayes"]], ci = size) %>% 
      select(Parameter, CI_low, CI_high)
    sigma_hdi <- bayestestR::hdi(draws$sigma, ci = size) %>% 
      mutate(Parameter = "sigma") %>% 
      select(Parameter, CI_low, CI_high)
    hdi_brm <- bind_rows(hdi_brm, sigma_hdi) %>% 
      rename(HDI_low = CI_low, HDI_high = CI_high)
    
    ci_brm <- bayestestR::ci(bfit()[["bayes"]], ci = size) %>% 
      select(Parameter, CI_low, CI_high)
    sigma_ci <- bayestestR::ci(draws$sigma, ci = size) %>% 
      mutate(Parameter = "sigma") %>% 
      select(Parameter, CI_low, CI_high)
    ci_brm <- bind_rows(ci_brm, sigma_ci)
    
    result_table <- tibble(
      Coefficient = c("Intercept", input$baseline, input$group, "sigma"),
      Mean = round(c(sum_brm$fixed$Estimate, sum_brm$spec_pars$Estimate), 2),
      Mode = round(mode_brm, 2),
      Median = round(median_brm$Median, 2),
      SD = round(sd_brm$SD, 2),
      !!paste0(size_label, " HDPI") := paste0("(", round(hdi_brm$HDI_low, 2),
                                              ", ", round(hdi_brm$HDI_high, 2),
                                              ")"),
      !!paste0(size_label, " CI") := paste0("(", round(ci_brm$CI_low, 2), ", ",
                                            round(ci_brm$CI_high, 2), ")")
    )
    
    DT::datatable(result_table)
  })
  
  output$dl_brmtable <- downloadHandler(
    filename = function() {
      paste0("bayesian_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(bfit(), draws())
      
      size <- input$sizeCI / 100
      size_label <- paste0(input$sizeCI, "%")
      
      sd_brm <- draws()[,1:4] %>%  
        map_dbl(sd) %>% 
        as_tibble() %>% 
        rename("SD" = "value")
      
      mode_brm <- draws()[,1:4] %>% 
        map_dbl(mode) %>% 
        as_tibble() %>% 
        rename("Mode" = "value")
      
      median_brm <- draws()[,1:4] %>% 
        map_dbl(median) %>% 
        as_tibble() %>% 
        rename("Median" = "value")
      
      sum_brm <- summary(bfit()[["bayes"]])
      
      hdi_brm <- bayestestR::hdi(bfit()[["bayes"]], ci = size) %>% 
        select(Parameter, CI_low, CI_high)
      sigma_hdi <- bayestestR::hdi(draws()$sigma, ci = size) %>% 
        mutate(Parameter = "sigma") %>% 
        select(Parameter, CI_low, CI_high)
      hdi_brm <- bind_rows(hdi_brm, sigma_hdi) %>% 
        rename(HDI_low = CI_low, HDI_high = CI_high)
      
      ci_brm <- bayestestR::ci(bfit()[["bayes"]], ci = size) %>% 
        select(Parameter, CI_low, CI_high)
      sigma_ci <- bayestestR::ci(draws()$sigma, ci = size) %>% 
        mutate(Parameter = "sigma") %>% 
        select(Parameter, CI_low, CI_high)
      ci_brm <- bind_rows(ci_brm, sigma_ci)
      
      result_table <- tibble(
        Coefficient = c("Intercept", input$baseline, input$group, "sigma"),
        Mean = round(c(sum_brm$fixed$Estimate, sum_brm$spec_pars$Estimate), 2),
        Mode = round(mode_brm$Mode, 2),
        Median = round(median_brm$Median, 2),
        SD = round(sd_brm$SD, 2),
        !!paste0(size_label, " HDPI") := paste0("(", round(hdi_brm$HDI_low, 2),
                                                ", ", round(hdi_brm$HDI_high, 2),
                                                ")"),
        !!paste0(size_label, " CI") := paste0("(", round(ci_brm$CI_low, 2), ", ",
                                              round(ci_brm$CI_high, 2), ")")
      )
      
      write.csv(result_table, file, row.names = FALSE)
    }
  )
  
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
  
  output$interplotpdf <- renderPlot({
    
    plotPosterior(draws_long(), "b_Intercept", input$sizeCI, input$plotCI,
                  input$plotCT)
  })
  
  output$baseplotpdf <- renderPlot({
    
    plotPosterior(draws_long(), paste0("b_", make.names(input$baseline)),
                  input$sizeCI, input$plotCI, input$plotCT)
  })
  
  output$teplotpdf <- renderPlot({
    
    plotPosterior(draws_long(), paste0("b_", make.names(input$group)),
                  input$sizeCI, input$plotCI, input$plotCT)
  })
  
  output$sdplotpdf <- renderPlot({
    
    plotPosterior(draws_long(), "sigma", input$sizeCI, input$plotCI,
                  input$plotCT)
  })
  
  output$interplotcdf <- renderPlot({
    
    plotPosteriorCDF(draws_long(), "b_Intercept", input$sizeCI, input$plotCI,
                     input$plotCT)
  })
  
  output$baseplotcdf <- renderPlot({
    
    plotPosteriorCDF(draws_long(), paste0("b_", make.names(input$baseline)),
                  input$sizeCI, input$plotCI, input$plotCT)
  })
  
  output$teplotcdf <- renderPlot({
    
    plotPosteriorCDF(draws_long(), paste0("b_", make.names(input$group)),
                  input$sizeCI, input$plotCI, input$plotCT)
  })
  
  output$sdplotcdf <- renderPlot({
    
    plotPosteriorCDF(draws_long(), "sigma", input$sizeCI, input$plotCI,
                     input$plotCT)
  })
  
  output$dl_interplotpdf <- downloadHandler(
    filename = function() { paste0("intercept_posterior_pdf_", Sys.Date(),
                                   ".png") },
    content = function(file) {
      ggsave(file, plot = plotPosterior(draws_long(), "b_Intercept",
                                        input$sizeCI, input$plotCI,
                                        input$plotCT),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_baseplotpdf <- downloadHandler(
    filename = function() { paste0("baseline_posterior_pdf_", Sys.Date(),
                                   ".png") },
    content = function(file) {
      ggsave(file, plot = plotPosterior(draws_long(),
                                        paste0("b_",
                                               make.names(input$baseline)),
                                        input$sizeCI, input$plotCI,
                                        input$plotCT),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_teplotpdf <- downloadHandler(
    filename = function() { paste0("treatment_posterior_pdf_", Sys.Date(),
                                   ".png") },
    content = function(file) {
      ggsave(file, plot = plotPosterior(draws_long(),
                                        paste0("b_", make.names(input$group)),
                                        input$sizeCI, input$plotCI,
                                        input$plotCT),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_sdplotpdf <- downloadHandler(
    filename = function() { paste0("sigma_posterior_pdf_", Sys.Date(),
                                   ".png") },
    content = function(file) {
      ggsave(file, plot = plotPosterior(draws_long(), "sigma",
                                        input$sizeCI, input$plotCI,
                                        input$plotCT),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_interplotcdf <- downloadHandler(
    filename = function() { paste0("intercept_posterior_cdf_", Sys.Date(),
                                   ".png") },
    content = function(file) {
      ggsave(file, plot = plotPosteriorCDF(draws_long(), "b_Intercept",
                                           input$sizeCI, input$plotCI,
                                           input$plotCT),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_baseplotcdf <- downloadHandler(
    filename = function() { paste0("baseline_posterior_cdf_", Sys.Date(),
                                   ".png") },
    content = function(file) {
      ggsave(file, plot = plotPosteriorCDF(draws_long(),
                                           paste0("b_",
                                                  make.names(input$baseline)),
                                           input$sizeCI, input$plotCI,
                                           input$plotCT),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_teplotcdf <- downloadHandler(
    filename = function() { paste0("treatment_posterior_cdf_", Sys.Date(),
                                   ".png") },
    content = function(file) {
      ggsave(file, plot = plotPosteriorCDF(draws_long(),
                                           paste0("b_",
                                                  make.names(input$group)),
                                           input$sizeCI, input$plotCI,
                                           input$plotCT),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  output$dl_sdplotcdf <- downloadHandler(
    filename = function() { paste0("sigma_posterior_cdf_", Sys.Date(),
                                   ".png") },
    content = function(file) {
      ggsave(file, plot = plotPosteriorCDF(draws_long(), "sigma",
                                           input$sizeCI, input$plotCI,
                                           input$plotCT),
             width = 10, height = 5, units = "in", dpi = 300)
    }
  )
  
  
}

shinyApp(ui = ui, server = server)