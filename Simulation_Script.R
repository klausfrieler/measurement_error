## Commands to run the simulations and analyses
library(rempsyc)
library(flextable)
library(tidyverse)

#Diagnostics and plotting funcitons

diagnostics_1 <- function(simul_data, metric = "rel_error", 
                        coef = c("all", "b0", "b1", "b2"),
                        x_var = "name", x_label = "Coefficient",
                        color_var = "measurement_error", color_label = "ME",
                        with_plot = T){
  coef <- match.arg(coef)
  if(coef != "all"){
    simul_data <- simul_data %>% filter(name == coef)
  }
  stats <- simul_data %>%
    group_by(me_raw = measurement_error_raw, me_diff = measurement_error_diff, error_type, method, name) %>%
    summarise(!!sym(metric) := mean(!!sym(metric)), .groups = "drop") %>%
    arrange(metric)
  if(with_plot){
    q <- simul_data %>% 
      mutate(!!sym(color_label) := factor(!!sym(color_var)),
             !!sym(x_label) := factor(!!sym(x_var))) %>%
      ggplot(aes(x = !!sym(x_label), y = !!sym(metric), color = !!sym(color_label)))
    q <- q + geom_boxplot(width = .2, outliers = FALSE)
    q <- q + geom_jitter(alpha = .1, position = position_dodge(width = .2))
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
  stats
}


diagnostics_2 <- function(simul_data, metric = "rel_error", 
                         coef = c("all", "b0", "b1", "b2"),
                         x_var = "name", x_label = "Coefficient",
                         color_var = "measurement_error", color_label = "ME",
                         with_plot = T){
  coef <- match.arg(coef)
  if(coef != "all"){
    simul_data <- simul_data %>% filter(name == coef)
  }
  stats <- simul_data %>%
    group_by(me_raw = measurement_error_raw, me_diff = measurement_error_diff, error_type, method, name) %>%
    summarise(!!sym(metric) := mean(!!sym(metric)), .groups = "drop") %>%
    arrange(metric)
  if(with_plot){
    q <- simul_data %>% 
      mutate(!!sym(color_label) := factor(!!sym(color_var)),
             !!sym(x_label) := factor(!!sym(x_var))) %>%
      ggplot(aes(x = !!sym(x_label), y = !!sym(metric), color = !!sym(color_label)))
    q <- q + geom_boxplot(width = .2, outliers = FALSE)
    q <- q + geom_jitter(alpha = .1, position = position_dodge(width = .2))
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
  stats
}

filter_solutions <- function(d, N, method){
  bad_solution <- NULL
  for(i in seq(along=N)){
    tmp_d <- d[d$N==N[i],]
    se_max <- max(tmp_d[tmp_d$method=="no_correction" ,"se"])
    se_cutoff <- se_max * sqrt(1 + 1.52^2)
    for(j in seq(along=method)){
      bad_solution <- c(bad_solution, tmp_d[tmp_d$method==method[j] & tmp_d$se > se_cutoff, "solution_index"])
    }
  }
  bad_solution <- unique(bad_solution)
  d <- subset(d, !(solution_index %in% bad_solution))
  d <- d[!is.na(d$se),]
}

### Generation of plots and tables
### Assumes that data has been generated already and is available in ME_sim/simulations

### Scenario 1

s1_all <- readRDS("simulations/scenario1_all.rds")

s1_unq_sol <- unique(s1_all[,c("N","batch","error_type","measurement_error","measurement_error_diff","method")])
s1_unq_sol$solution_index <- 1:nrow(s1_unq_sol)
s1_all <- merge(s1_all, s1_unq_sol)

s1_clean <- filter_solutions(s1_all, N=c(50,500,5000),method=c("LV","MI","simex","outlier_exclusion"))

#Figure 1
png("figures/s1_all_allCoefs.png")
diagnostics_1(s1_all[s1_all$measurement_error_diff==0,],metric="rel_error")
dev.off()

#Figure 2
png("figures/s1_all_allCoefs_SE.png")
diagnostics_1(s1_all[s1_all$measurement_error_diff==0,],metric="se")
dev.off()

#Table A1: Coefficient SE by level of ME, sample size and correction method 
T_A1 <- s1_all[s1_all$measurement_error_diff==0,] %>% 
  group_by(method,N,measurement_error) %>%
  dplyr::summarize(mean_se = mean(se, na.rm=T), min_se = min(se, na.rm=T),max_se = max(se, na.rm=T))
nT_A1 <- nice_table(T_A1)
#flextable::save_as_docx(nT_A1, path = "S1_SEs_by_MExNxMethod_table.docx")

#Compute SE threshold for filtering bad results
#se_cutoff <- max(s1_all[s1_all$method=="no_correction","se"])
#s1_clean <- na.omit(s1_all[s1_all$se<=se_cutoff,])
#s1_bad <- na.omit(s1_all[s1_all$se>se_cutoff,])
#table(s1_bad$method)/sum(table(s1_bad$method))

#Figure 3
png("figures/s1_all_allCoefs_clean.png")
diagnostics_1(s1_all[s1_all$measurement_error_diff==0,],metric="rel_error")
dev.off()

#Figures 4a-c
png("figures/s1_n=50_b1_clean.png")
diagnostics_2(s1_clean[s1_clean$N==50,], coef = "b1", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME")
dev.off()

png("figures/s1_n=500_b1_clean.png")
diagnostics_2(s1_clean[s1_clean$N==500,], coef = "b1", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME")
dev.off()

png("figures/s1_n=5000_b1_clean.png")
diagnostics_2(s1_clean[s1_clean$N==5000,], coef = "b1", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME")
dev.off()

T_A2 <- s1_clean[s1_clean$name=="b1",] %>% 
  group_by(method,n_sample,measurement_error_raw, error_type, measurement_error_diff) %>%
  dplyr::summarize(mean_b1 = mean(rel_error, na.rm=T), 
            CI_LB = mean(rel_error, na.rm=T) - 2 * sd(rel_error, na.rm = TRUE),
            CI_UB = mean(rel_error, na.rm=T) + 2 * sd(rel_error, na.rm = TRUE)
              )
nT_A2 <- nice_table(T_A2)

### Scenario 2

s2_all <- readRDS("simulations/scenario2_all.rds")

s2_unq_sol <- unique(s2_all[,c("N","batch","error_type","measurement_error","measurement_error_diff","method")])
s2_unq_sol$solution_index <- 1:nrow(s2_unq_sol)
s2_all <- merge(s2_all, s2_unq_sol)

s2_clean <- filter_solutions(s2_all, N=c(50,500,5000),method=c("LV","MI","simex","outlier_exclusion"))

#Figures 5a-f
png("figures/s2_n=50_b1_clean.png")
diagnostics_2(s2_clean[s2_clean$N==50,], coef = "b1", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME")
dev.off()

png("figures/s2_n=50_b2_clean.png")
diagnostics_2(s2_clean[s2_clean$N==50,], coef = "b2", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME")
dev.off()

png("figures/s2_n=500_b1_clean.png")
diagnostics_2(s2_clean[s2_clean$N==500,], coef = "b1", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME")
dev.off()

png("figures/s2_n=500_b2_clean.png")
diagnostics_2(s2_clean[s2_clean$N==500,], coef = "b2", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME")
dev.off()

png("figures/s2_n=5000_b1_clean.png")
diagnostics_2(s2_clean[s2_clean$N==5000,], coef = "b1", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME")
dev.off()

png("figures/s2_n=5000_b2_clean.png")
diagnostics_2(s2_clean[s2_clean$N==5000,], coef = "b2", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME")
dev.off()


### Scenario 3

###### Data simulations

# Scenario 1, sample size N = 50
#s1_n50 <- ME_simulator$new(scenario = 1, n_sample = 50, b0 = 1, b1 = -0.5, b2 = 0.25)
#s1_n50$run()
#s1_n50$save("s1_n50.rds")

# sample size N = 500
#s1_n500 <- ME_simulator$new(scenario = 1, n_sample = 500, b0 = 1, b1 = -0.5, b2 = 0.25)
#s1_n500$run()
#s1_n500$save("s1_n500.rds")

# sample size N = 5000
#s1_n5000 <- ME_simulator$new(scenario = 1, n_sample = 5000, b0 = 1, b1 = -0.5, b2 = 0.25)
#s1_n5000$run()
#s1_n5000$save("s1_n5000.rds")

# aggregate the 3 sample sizes and compute overview of distribution of coefficient bias and standard error
#dat50 <- data.frame(s1_n50$results,N=rep(50,nrow(s1_n50$results)))
#dat500 <- data.frame(s1_n500$results,N=rep(500,nrow(s1_n500$results)))
#dat5000 <- data.frame(s1_n5000$results,N=rep(5000,nrow(s1_n5000$results)))
#dat <- rbind(dat50,dat500,dat5000)
#save(dat,file="simulations/scenario1_all.rds")

### Scenario 2
# sample size N = 50
#s2_n50 <- ME_simulator$new(scenario = 2, n_sample = 50, n_batch=50, b0 = 1, b1 = 0.4, b2 = 0.2)
#s2_n50$run()
#s2_n50$save("s2_n50.rds")

# sample size N = 500
#s2_n500 <- ME_simulator$new(scenario = 2, n_sample = 500, n_batch=50, b0 = 1, b1 = 0.4, b2 = 0.2)
#s2_n500$run()
#s2_n500$save("s2_n500.rds")

# sample size N = 5000
#s2_n5000 <- ME_simulator$new(scenario = 2, n_sample = 5000, n_batch=50, b0 = 1, b1 = 0.4, b2 = 0.2)
#s2_n5000$run()
#s2_n5000$save("s2_n5000.rds")

#dat50 <- data.frame(s2_n50$results,N=rep(50,nrow(s2_n50$results)))
#dat500 <- data.frame(s2_n500$results,N=rep(500,nrow(s2_n500$results)))
#dat5000 <- data.frame(s2_n5000$results,N=rep(5000,nrow(s2_n5000$results)))
#dat <- rbind(dat50,dat500,dat5000)
#save(dat,file="simulations/scenario2_all.rds")