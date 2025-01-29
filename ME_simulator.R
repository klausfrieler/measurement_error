require(univOutl)
require(lavaan)
require(TSI) # install from github: install.packages("remotes") /n remotes::install_github("mmansolf/TSI")
require(simex)
require(truncnorm)
library(tidyverse)

source("scenario1.R")
source("scenario2.R")
source("scenario3.R")

messagef <- function(...) message(sprintf(...))
scale_ <- function(x) as.numeric(scale(x))

#'@export
ME_simulator <- R6::R6Class("ME_simulator",
                              
                              private = list(
                              ),
                              
                              public = list(
                                coefficients = list(b0 = 1, b1 = .4, b2 = .2),
                                error_types = c("classical", "systematic", "heteroscedastic", "differential"),
                                measurement_errors = c(.05, 0.33, 0.67, 1, 1.52),
                                me_diffs = c(0, 0.1, 0.5),
                                methods = c("no_correction",
                                             "outlier_exclusion",
                                             #"Weighting",
                                             "LV",
                                             "MI",
                                             "simex"),
                                n_sample = 50,
                                n_batch = 50,
                                normal = F,
                                verbose = F,
                                seed = NULL,
                                scenario = 1,
                                results = NULL,
                                result_stats = NULL,
                                simulated_data = NULL,
                                
                                initialize = function(b0 = 1,
                                                      b1 = 0.4,
                                                      b2 = 0.2,
                                                      error_types = c("classical", "systematic", "heteroscedastic", "differential"),
                                                      measurement_errors = c(.05, 0.33, 0.67, 1, 1.52),
                                                      me_diffs = c(0, 0.1, 0.5),
                                                      n_sample = 50,
                                                      n_batch = 50,
                                                      normal = FALSE,
                                                      verbose = FALSE,
                                                      scenario = c(1, 2, 3),
                                                      methods = c("no_correction",
                                                                  "outlier_exclusion",
                                                                  #"Weighting",
                                                                  "LV",
                                                                  "MI",
                                                                  "simex"),
                                                      seed = NULL) {
                                  self$coefficients <- list(b0 = b0, b1 = b1, b2 = b2)
                                  self$error_types <- error_types
                                  self$methods <- methods
                                  self$scenario <- scenario[1]
                                  self$n_sample <- n_sample
                                  self$n_batch <- n_batch
                                  self$measurement_errors <- measurement_errors
                                  self$me_diffs <- me_diffs
                                  self$verbose <- verbose
                                  self$seed <- seed
                                  self$normal <- normal
                                  self$results <- NULL
                                  self$result_stats <- NULL
                                  self$simulated_data <- NULL
                                  #if(scenario == 2){
                                  #  self$methods <- setdiff(self$methods, "MI")
                                  #}
                                  if(scenario == 3){
                                    if(length(error_types) > 1 || error_types[[1]] != "heteroscedastic"){
                                      warning("For scenario 3 only heteroscedastic error available. Fixing.")
                                    }
                                    self$error_types <- "heteroscedastic"
                                    if(any(!is.character(measurement_errors))){
                                      browser()
                                      stop("Scenario 3 uses categorial error levels (very_high, high, medium and low)")
                                    }
                                    bad_levels <- setdiff(measurement_errors, scenario3_me_levels)
                                    if(length(bad_levels) > 1){
                                      browser()
                                      stop("Scenario 3 uses categorial error levels (very_high, high, medium and low)")
                                    }
                                    self$measurement_errors <- measurement_errors
                                    self$me_diffs <- 0
                                  }
                                  invisible(self)
                                },
                                
                                print = function(...) {
                                  cat("ME Simulator: \n")
                                  cat("  Coefficients: \n", 
                                      sprintf("    %s = %.2f\n", names(self$coefficients), self$coefficients), 
                                      "\n", 
                                      sep = "")
                                  cat("  Metadata:  ", "\n", 
                                      sprintf("    %s: %s\n", 
                                              c("scenario", 
                                                "N (sample)", 
                                                "N (batch)", 
                                                "Seed",
                                                "Measurement Errors", 
                                                "ME differences", 
                                                "Error types", 
                                                "Corrections"), 
                                              c(self$scenario, 
                                                self$n_sample, 
                                                self$n_batch, 
                                                ifelse(is.null(self$seed), "--", self$seed),
                                                paste(self$measurement_errors, collapse = ", "),
                                                paste(self$me_diffs, collapse = ", "),
                                                paste(self$error_types, collapse = ", "),
                                                paste(self$methods, collapse = ", ")
                                                )), "\n", sep = "")
                                  cat("  State:  ", "\n", 
                                      sprintf("    %s: %s\n", "No. Results", length(self$results)), sep = "")
                                  invisible(self)
                                },
                                
                                generate_data = function(error_type, me_raw){
                                  if(self$scenario == 1){
                                    generator_func <- generate_data_scenario1
                                  }
                                  else if(self$scenario == 2){
                                    generator_func <- generate_data_scenario2
                                  }
                                  else if(self$scenario == 3){
                                    generator_func <- generate_data_scenario3
                                  }
                                  else{
                                    stop(sprintf("Invalid scenario: %s", scenario))
                                  }
                                  generator_func(
                                    b0 = self$coefficients$b0,
                                    b1 = self$coefficients$b1,
                                    b2 = self$coefficients$b2,
                                    error_type = error_type,
                                    me = me_raw,
                                    n_sample = self$n_sample,
                                    n_batch = self$n_batch,
                                    normal = self$normal
                                  )
                                },
                                
                                get_coefficients = function(df, method, me, scenario){
                                  if(self$scenario == 1){
                                    coef_func <- get_coefs_scenario1
                                  }
                                  else if(self$scenario == 2){
                                    coef_func <- get_coefs_scenario2
                                  }
                                  else if(self$scenario == 3){
                                    coef_func <- get_coefs_scenario3
                                  }
                                  else{
                                    stop(sprintf("Invalid scenario: %s", scenario))
                                  }
                                  coef_func(df, method, me)
                                },
                                
                                eval  = function(coefs) {
                                  coefs <- coefs %>% 
                                    mutate(true_coefs = unlist(self$coefficients),  
                                           error = value - true_coefs, 
                                           rel_error = error/true_coefs,
                                           abs_error = abs(error),
                                           rel_abs_error = abs_error/abs(true_coefs)
                                    )
                                  coefs %>% select(name = term, true_coefs, coefs = value, everything())
                                },
                                
                                run  = function(fname = NULL, seed = NULL) {
                                  #sample size (n) is a parameter, possible values are 50, 500, 5000; but no separate loop, instead run simulator several times

                                  if(!is.null(seed)){
                                    set.seed(seed)
                                  }
                                  if(!is.null(self$seed)){
                                    set.seed(seed)
                                  }
                                  self$results <- 
                                    map_dfr(self$measurement_errors, function(me_raw) {
                                      map_dfr(self$error_types, function(et) {
                                      if(self$scenario != 3){
                                        messagef("Simulating data with '%s' error (me_raw = %.2f)...", et, me_raw)
                                      }
                                      else{
                                        messagef("Simulating data with '%s' error (me_raw = %s)...", et, me_raw)
                                      }
                                        simulated_data <- self$generate_data(error_type = et, me = me_raw)
                                      
                                      map_dfr(self$me_diffs, function(md) {
                                        if(self$scenario != 3){
                                          me <- me_raw * (1 +  md)
                                          messagef("...  %d batches with me_diff = %.2f -> me_total = %.2f.", self$n_batch, md, me)
                                        }
                                        else{
                                          me <- me_raw
                                        }
                                        
                                        map_dfr(1:self$n_batch, function(i) {
                                          #simulation
                                          #rowser()
                                          messagef("************* BATCH: %d ***********", i)
                                          #analysis methods and evaluation, i.e. comparison of coefficients to ground truth from simulation
                                          df <- simulated_data[simulated_data$batch == i, ] %>% select(-batch)
                                          #browser()
                                          map_dfr(self$methods, function(method) {
                                            if (self$verbose)
                                              messagef("...using correction method '%s' [%d] ", method, i)
                                            coefs <- self$get_coefficients(df, method, me)
                                            #browser()
                                            ret <- self$eval(coefs) %>%
                                              mutate(
                                                method = !!method,
                                                batch = i,
                                                error_type = et,
                                                measurement_error = me,
                                                measurement_error_diff = md
                                              )
                                            #browser()
                                            ret
                                          })
                                        })
                                      })
                                    })  %>%
                                      mutate(n_sample = self$n_sample,
                                             n_batch = self$n_batch,
                                             measurement_error_raw = me_raw)
                                  })
                                  if(self$scenario == 3){
                                    self$results <- self$results %>% 
                                      mutate(measurement_error = factor(measurement_error, 
                                                                        levels = scenario3_me_levels))
                                  }
                                  self$get_stats()
                                  #print(self$result_stats)
                                  if(!is.null(fname)){
                                    self$save(fname)
                                  }
                                  invisible(self)
                                  
                                },
                                
                                save = function(fname){
                                  saveRDS(self, fname)
                                },
                                
                                load = function(fname){
                                  readRDS(fname)
                                },
                                
                                get_stats = function(){
                                  if(is.null(self$results)){
                                    warning("No results available, run simulation first.")
                                  }
                                  
                                  metrics <- c("error", 
                                              "rel_error",
                                              "abs_error",
                                              "rel_abs_error")
                                  
                                  self$result_stats <- self$results %>%
                                    group_by(me_raw = measurement_error_raw, 
                                             me_diff = measurement_error_diff, 
                                             error_type, method, 
                                             name) %>%
                                    summarise(across(all_of(metrics), mean), .groups = "drop") 
                                  
                                },
                                #' @param metric to display on y-axis
                                #' @param coef use all or a single coefficient
                                #' @param x_var Use this on x-axis
                                #' @param x_label Label on x-aixs 
                                #' @param color_var Use this for color code
                                #' @param color_label Label for color_var (legend title)
                                diagnostics = function(metric = "rel_error",
                                                       coef = c("all", "b0", "b1", "b2"),
                                                       x_var = "name",
                                                       x_label = "Coefficient",
                                                       color_var = "measurement_error",
                                                       color_label = "ME",
                                                       max_se = NULL,
                                                       alpha = .1,
                                                       with_plot = T) {
                                  coef <- match.arg(coef)
                                  metric <- metric[[1]]
                                  if(is.null(self$results)){
                                    warning("No results available, run simulation first.")
                                    return(invisible(self))
                                  }
                                  simul_data <- self$results
                                  if(coef != "all"){
                                    simul_data <- simul_data %>% filter(name == coef)
                                  }
                                  if(!is.null(max_se)){
                                    simul_data <- simul_data %>% filter(se <= max_se)
                                    if(nrow(simul_data) == 0){
                                      warning("No results left after SE filtering.")
                                      return(invisible(self))
                                    }
                                  }
                                  if(with_plot){
                                    q <- simul_data %>% 
                                      mutate(!!sym(color_label) := factor(!!sym(color_var)),
                                             !!sym(x_label) := factor(!!sym(x_var))) %>%
                                      ggplot(aes(x = !!sym(x_label), 
                                                 y = !!sym(metric), 
                                                 color = !!sym(color_label)))
                                    q <- q + geom_boxplot(width = .2, outliers = FALSE)
                                    q <- q + geom_jitter(alpha = alpha, position = position_dodge(width = .2))
                                    q <- q + facet_grid(error_type ~ method)
                                    q <- q  + theme_bw()
                                    if(coef != "all"){
                                      q <- q + ggtitle(sprintf("coeficient = %s", coef))
                                    }
                                    if(str_detect(metric, "rel")){
                                      q <- q + scale_y_continuous(labels = scales::percent)
                                    }
                                    print(q)
                                  }
                                  if(is.null(self$result_stats)){
                                    self$get_stats()
                                  }
                                  ret <- self$result_stats %>% filter(metric == !!metric)
                                  if(coef != "all"){
                                    ret <- ret %>% filter(name == coef)
                                  }
                                  ret
                                }
                                
                                ),

                              # End public
                              #####################
                              active = list(

                              )
)

