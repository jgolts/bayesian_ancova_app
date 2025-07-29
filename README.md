# bayesian_ancova_app
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
![R Version](https://img.shields.io/badge/R-4.5.0-blue)

This is an R Shiny app that performs Bayesian ANCOVA on pre/post measurement datasets. It was created originally as part of my research project on Bayesian ANCOVA for pre/post data, code for which can be found in [my bancova_analysis repo](https://github.com/jgolts/bancova_analysis). Please feel free to fork it and modify to your liking. To run it without downloading the files, install the R package `shiny`, and enter `runGitHub( "bayesian_ancova_app", "jgolts")` into the console.

## Features
+ Custom dataset upload in .xls/.xlsx format, currently supporting only wide format data
+ Analysis currently restricted to two groups only, with covariate adjustment only for pre-test value
+ Performs frequentist and Bayesian ANCOVA for pre/post data
+ Allows user-defined priors for all model parameters - normal priors for coefficients, and half-Cauchy for SD of residuals
+ Visualises priors as density plots
+ Custom MCMC iterations and seed
+ Auto-detects number of cores and sets number of chains accordingly, to a minimum of 4 chains per STAN default
+ Frequentist model results (coefficients, p-values, confidence intervals)
+ ETI and HDI credible intervals supported at user-specified width for each parameter
+ Displays mean, median, and SD of posterior distribution for each parameter
+ MCMC diagnostics: Trace plots and Rhat, repeated at double the iterations, autocorrelation graphs, posterior histograms
+ Posterior distribution plots with custom width ETI/HDI and mean/median - PDF and CCDF
+ Download options for all plots and tables

## Planned Features
+ Expansion to 3+ group analysis
+ Additional covariate adjustment
