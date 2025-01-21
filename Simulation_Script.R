## Commands to run the simulations and analyses
library(rempsyc)
library(flextable)
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

### Scenario 1

# sample size N = 50
s1_n50 <- ME_simulator$new(scenario = 1, n_sample = 50, b0 = 1, b1 = -0.5, b2 = 0.25)
s1_n50$run()
s1_n50$save("s1_n50.rds")

# sample size N = 500
s1_n500 <- ME_simulator$new(scenario = 1, n_sample = 500, b0 = 1, b1 = -0.5, b2 = 0.25)
s1_n500$run()
s1_n500$save("s1_n500.rds")

# sample size N = 5000
s1_n5000 <- ME_simulator$new(scenario = 1, n_sample = 5000, b0 = 1, b1 = -0.5, b2 = 0.25)
s1_n5000$run()
s1_n5000$save("s1_n5000.rds")

# aggregate the 3 sample sizes and compute overview of distribution of coefficient bias and standard error
dat50 <- data.frame(s1_n50$results,N=rep(50,nrow(s1_n50$results)))
dat500 <- data.frame(s1_n500$results,N=rep(500,nrow(s1_n500$results)))
dat5000 <- data.frame(s1_n5000$results,N=rep(5000,nrow(s1_n5000$results)))
dat <- rbind(dat50,dat500,dat5000)
save(dat,file="S1_Results_AllN.rds")

#Figure 1
png("S1_AllN_AllCoefs.png")
diagnostics_1(dat[dat$measurement_error_diff==0,],metric="rel_error")
devo.off()

#Figure 2
png("S1_AllN_AllCoefs_SE.png")
diagnostics_1(dat[dat$measurement_error_diff==0,],metric="se")
dev.off()

#Table A1
T_A1 <- dat[dat$measurement_error_diff==0,] %>% 
  group_by(method,N,measurement_error) %>%
  summarize(mean_se = mean(se, na.rm=T), min_se = min(se, na.rm=T),max_se = max(se, na.rm=T))
nT_A1 <- nice_table(T_A1)
#flextable::save_as_docx(nT_A1, path = "S1_SEs_by_MExNxMethod_table.docx")

#Compute SE threshold for filtering bad results
se_cutoff <- max(dat[dat$method=="no_correction","se"])
dat_clean <- na.omit(dat[dat$se<=se_cutoff,])
dat_bad <- na.omit(dat[dat$se>se_cutoff,])
table(dat_bad$method)/sum(table(dat_bad$method))

#Figure 1b
png("S1_AllN_AllCoefs_AfterFiltering.png")
diagnostics_1(dat_clean[dat_clean$measurement_error_diff==0,],metric="rel_error")
dev.off()

#Figures 3a-c
png("S1_N50_b1_cleandat.png")
diagnostics_2(dat_clean[dat_clean$N==50,], coef = "b1", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME raw")
dev.off()

png("S1_N500_b1_cleandat.png")
diagnostics_2(dat_clean[dat_clean$N==500,], coef = "b1", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME raw")
dev.off()

png("S1_N5000_b1_cleandat.png")
diagnostics_2(dat_clean[dat_clean$N==5000,], coef = "b1", x_var = "measurement_error_diff", x_label = "ME Difference", color_var = "measurement_error_raw", color_label = "ME raw")
dev.off()

### Scenario 2
# sample size N = 50
s2_n50 <- ME_simulator$new(scenario = 2, n_sample = 50, n_batch=50, b0 = 1, b1 = 0.4, b2 = 0.2)
s2_n50$run()
s2_n50$save("s2_n50.rds")

# sample size N = 500
s2_n500 <- ME_simulator$new(scenario = 2, n_sample = 500, n_batch=50, b0 = 1, b1 = 0.4, b2 = 0.2)
s2_n500$run()
s2_n500$save("s2_n500.rds")

# sample size N = 5000
s2_n5000 <- ME_simulator$new(scenario = 2, n_sample = 5000, n_batch=50, b0 = 1, b1 = 0.4, b2 = 0.2)
s2_n5000$run()
s2_n5000$save("s2_n5000.rds")

dat50 <- data.frame(s2_n50$results,N=rep(50,nrow(s2_n50$results)))
dat500 <- data.frame(s2_n500$results,N=rep(500,nrow(s2_n500$results)))
dat5000 <- data.frame(s2_n5000$results,N=rep(5000,nrow(s2_n5000$results)))
dat <- rbind(dat50,dat500,dat5000)
save(dat,file="S2_Results_AllN.rds")