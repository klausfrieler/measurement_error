library(rempsyc)
library(flextable)
diagnostics_1 <- function(simul_data,
                          metric = "rel_error",
                          coef = c("all", "b0", "b1", "b2"),
                          x_var = "name",
                          x_label = "Coefficient",
                          color_var = "measurement_error",
                          color_label = "ME",
                          with_plot = T) {
  coef <- match.arg(coef)
  if (coef != "all") {
    simul_data <- simul_data %>% filter(name == coef)
  }
  stats <- simul_data %>%
    group_by(me_raw = measurement_error_raw,
             me_diff = measurement_error_diff,
             error_type,
             method,
             name) %>%
    summarise(!!sym(metric) := mean(!!sym(metric)), .groups = "drop") %>%
    arrange(metric)
  if (with_plot) {
    q <- simul_data %>%
      mutate(!!sym(color_label) := factor(!!sym(color_var)),!!sym(x_label) := factor(!!sym(x_var))) %>%
      ggplot(aes(
        x = !!sym(x_label),
        y = !!sym(metric),
        color = !!sym(color_label)
      ))
    q <- q + geom_boxplot(width = .2, outliers = FALSE)
    q <- q + geom_jitter(alpha = .1, position = position_dodge(width = .2))
    q <- q + facet_grid(error_type ~ method)
    q <- q  + theme_bw()
    if (coef != "all") {
      q <- q + ggtitle(sprintf("coeficient = %s", coef))
    }
    if (str_detect(metric, "rel")) {
      q <- q + scale_y_continuous(labels = scales::percent)
    }
    print(q)
  }
  stats
}


diagnostics_2 <- function(simul_data,
                          metric = "rel_error",
                          coef = c("all", "b0", "b1", "b2"),
                          x_var = "name",
                          x_label = "Coefficient",
                          color_var = "measurement_error",
                          color_label = "ME",
                          with_plot = T) {
  coef <- match.arg(coef)
  
  if (coef != "all") {
    simul_data <- simul_data %>% filter(name == coef)
  }
  
  stats <- simul_data %>%
    group_by(me_raw = measurement_error_raw,
             me_diff = measurement_error_diff,
             error_type,
             method,
             name) %>%
    summarise(!!sym(metric) := mean(!!sym(metric)), .groups = "drop") %>%
    arrange(metric)
  
  if (with_plot) {
    q <- simul_data %>%
      mutate(!!sym(color_label) := factor(!!sym(color_var)),!!sym(x_label) := factor(!!sym(x_var))) %>%
      ggplot(aes(
        x = !!sym(x_label),
        y = !!sym(metric),
        color = !!sym(color_label)
      ))
    q <- q + geom_boxplot(width = .2, outliers = FALSE)
    q <- q + geom_jitter(alpha = .1, position = position_dodge(width = .2))
    q <- q + facet_grid(error_type ~ method)
    q <- q  + theme_bw()
    
    if (coef != "all") {
      q <- q + ggtitle(sprintf("coeficient = %s", coef))
    }
    
    if (str_detect(metric, "rel")) {
      q <- q + scale_y_continuous(labels = scales::percent)
    }
    
    print(q)
  }
  
  stats
}
#error_types = c("classical", "systematic", "heteroscedastic", "differential"),

simu_def_full <- expand_grid(
  n_samples = c(50, 500, 5000),
  scenario = c(1, 2),
  n_batch = 50,
  error_types = c("classical", "systematic", "heteroscedastic", "differential"),
  measurement_errors = c(.05, 0.33, 0.67, 1, 1.52),
  me_diffs = c(0, 0.1, 0.5),
  methods = c("no_correction", "outlier_exclusion", "LV", "MI", "simex")
)

simu_def_test_12 <- expand_grid(
  n_samples = c(50, 500),
  scenario = c(1, 2),
  n_batch = 1,
  error_types = c("classical"),
  measurement_errors = as.character(c(.05, .5)),
  me_diffs = c(0),
  methods = c("no_correction")
)

simu_def_test_3 <- expand_grid(
  n_samples = c(50, 500),
  scenario = c(3),
  n_batch = 1,
  error_types = c("heteroscedastic"),
  measurement_errors = c("low", "medium"),
  me_diffs = c(0),
  methods = c("no_correction")
)

simu_def_test <- bind_rows(simu_def_test_12, simu_def_test_3)

ground_truth <- list(
  scenario1 = c(b0 = 1, b1 = -0.5, b2 = 0.25),
  scenario2 = c(b0 = 1, b1 = 0.4, b2 = 0.2),
  scenario3 = c(b0 = -1, b1 = 0.5, b2 = 0.3)
) 


### Scenario 1
run_simulations <- function(scenarios = c(1, 2, 3), id = "scenario", config = simu_def){
  ret <- map_dfr(scenarios, function(sc) {
    
    gt <- ground_truth[[sprintf("scenario%d", sc)]]
    def <- config %>%  filter(scenario == sc)
    n_batch <- def$n_batch[1]
    
    me <- def %>% pull(measurement_errors) %>% unique()    
    if(sc != 3){
      me <- as.numeric(me)
    }
    map_dfr(def$n_samples, function(n) {
      if(sc == 3){
        browser()
      }
      fname <- sprintf("%s%d_n=%d.rds", id, sc, n)
      simulator <- ME_simulator$new(
        scenario = sc,
        n_sample = n,
        n_batch = n_batch,
        b0 = gt[["b0"]],
        b1 = gt[["b1"]],
        b2 = gt[["b2"]],
        me_diffs = def %>% pull(me_diffs) %>% unique(),
        methods = def %>% pull(methods) %>% unique(),
        error_types = def %>% pull(error_types) %>% unique(),
        measurement_errors = me,
      )
      simulator$run()
      simulator$save(fname)
      #browser()
      simulator$results %>% mutate(N = n_sample,
                                   measurement_error = as.character(simulator$results$measurement_error),
                                   measurement_error_raw = as.character(simulator$results$measurement_error_raw),
                                   scenario = sc,
                                   id = id)
    })
  })
  cutoff <- ret %>% filter(method == "no_correction") %>% pull(se) %>% max()
  ret <- ret %>% mutate(status = factor(se <= cutoff, levels = c(FALSE, TRUE), labels = c("bad", "clean")))
  saveRDS( ret, file = sprintf("%s_all.rds", id))
  ret
}

make_tables <- function(simulations, with_save = F){
  T_A1 <- simulations %>% filter(measurement_error_diff == 0) %>%
    group_by(method, N, measurement_error) %>%
    summarize(
      mean_se = mean(se, na.rm = T),
      min_se = min(se, na.rm = T),
      max_se = max(se, na.rm = T)
    )
  nT_A1 <- nice_table(T_A1)
  if(with_save){
    flextable::save_as_docx(nT_A1, path = "S1_SEs_by_MExNxMethod_table.docx")
  }
  nT_A1
  
}

make_plots <- function(simulations, with_save = T, id = "test"){
  dat_clean <- simulations %>% filter(status == "clean")
  
  diagnostics_1(dat_clean %>% filter(scenario == 1, measurement_error_diff == 0), metric = "rel_error")
  if(with_save)ggsave(sprintf("%s1_all_clean.png", id))
  
  #Figures 3a-c
  sample_sizes <- unique(simulations %>% filter(scenario == 1) %>% pull(N))
  for(i in sample_sizes){
    browser()
    diagnostics_2(
      dat_clean %>% filter(N == i),
      coef = "b1",
      x_var = "measurement_error_diff",
      x_label = "ME Difference",
      color_var = "measurement_error_raw",
      color_label = "ME raw"
    )
    if(with_save)ggsave(sprintf("%s1_N=%d_b1_clean.png", id, i))

  }
  

}