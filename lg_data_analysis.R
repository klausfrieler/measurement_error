### Study 2: Anaysis empirical data from LongGold
library(tidyverse)
source("my_TSI.R")

analyze_lg_data <- function(){
  lg_dat <- read.csv("data/lg_data.csv", colClasses=c("p_id"="character")) %>% as_tibble()
  #lg_dat <- read.csv("c:/Users/klaus.frieler/Downloads/select_vars.csv",colClasses=c("p_id"="character"))
  ### Pre-process data:
  ### use only cases where BAT, MIQ, GMS.MT are not missing
  lg_dat <- lg_dat %>% filter(!is.na(BAT.score), !is.na(MIQ.score), !is.na(GMS.musical_training))
  
  ### remove all cases without a valid ID (p_id not starting with "0")
  lg_dat <- lg_dat %>% filter(substr(p_id, 1, 1) == "0")
  
  ### use each p_id only once, selecting the observation from the latest test wave or highest age
  
  lg_dat <- lg_dat %>% 
    group_by(p_id) %>% 
    mutate(last = test_wave == max(test_wave)) %>% 
    ungroup() %>% 
    filter(last) %>% select(-last)
  
  ### demographics
  
  # mean(lg_dat$age)
  # sd(lg_dat$age)
  # table(lg_dat$gender) / sum(table(lg_dat$gender))
  
  ### Data analysis with and without correction; only MI correction for now
  
  m_naive <- lm(BAT.score ~ MIQ.score + GMS.musical_training, data = lg_dat)
  coefs_naive <- broom::tidy(m_naive) %>% 
    select(term, value = estimate, se = std.error) %>% 
    mutate(term = c("b0", "b1", "b2"), 
           method = "no_correction")
  
  ##outlier exclusion
  browser()
  ol_y <- suppressMessages(boxB(
    x = lg_dat$BAT.error,
    k = 1.5,
    method = 'resistant'
  )$outliers)
  
  ol_x <- suppressMessages(boxB(
    x = lg_dat$MIQ.error,
    k = 1.5,
    method = 'resistant'
  )$outliers)
  
  ol_z <- suppressMessages(boxB(
    x = lg_dat$GMS.musical_training,
    k = 1.5,
    method = 'resistant'
  )$outliers)
  
  
  ol <- c(ol_y, ol_x, ol_z)
  tmp_data <- lg_dat
  if (length(ol) >= 1) {
    #coefs <- coef(lm(y ~ x + z, data = df[-ol, ]))
    tmp_data <- lg_dat[-ol, ]
  }
  coefs_outlier <- broom::tidy(lm(BAT.score ~ MIQ.score + GMS.musical_training, data = tmp_data)) %>%
    select(term, value = estimate, se = std.error) %>%
    mutate(term = c("b0", "b1", "b2"), 
           method = "outlier_exclusion")
  
  
  ###weighting
  inv_err <- (1 /sqrt(lg_dat$MIQ.error ^ 2 + lg_dat$BAT.error ^ 2)) #* (1 / lg_dat$GMS.Musical_training.error ^ 2)
  
  coefs_w0 <- broom::tidy(lm(BAT.score ~ MIQ.score + GMS.musical_training, weights = inv_err, data = lg_dat)) %>% 
    select(term, value = estimate, se = std.error) %>% 
    mutate(term = c("b0", "b1", "b2"),  
           method = "weighting0")
  inv_err <- (1 / (lg_dat$MIQ.error ^ 2 + lg_dat$BAT.error ^ 2)) #* (1 / lg_dat$GMS.Musical_training.error ^ 2)
  
  coefs_w1 <- broom::tidy(lm(BAT.score ~ MIQ.score + GMS.musical_training, weights = inv_err, data = lg_dat)) %>% 
    select(term, value = estimate, se = std.error) %>% 
    mutate(term = c("b0", "b1", "b2"),  
           method = "weighting1")

  inv_err <- (1 / (lg_dat$MIQ.error ^ 2) * (1 / lg_dat$BAT.error ^ 2)) #* (1 / lg_dat$GMS.Musical_training.error ^ 2)
  
  coefs_w2 <- broom::tidy(lm(BAT.score ~ MIQ.score + GMS.musical_training, weights = inv_err, data = lg_dat)) %>% 
    select(term, value = estimate, se = std.error) %>% 
    mutate(term = c("b0", "b1", "b2"),  
           method = "weighting2")

  inv_err <- (1 / (lg_dat$MIQ.error ^ 2) + (1 / lg_dat$BAT.error ^ 2)) #* (1 / lg_dat$GMS.Musical_training.error ^ 2)
  
  coefs_w3 <- broom::tidy(lm(BAT.score ~ MIQ.score + GMS.musical_training, weights = inv_err, data = lg_dat)) %>% 
    select(term, value = estimate, se = std.error) %>% 
    mutate(term = c("b0", "b1", "b2"),  
           method = "weighting3")
  
  inv_err <- 1 / (lg_dat$MIQ.error ^ 2) #* (1 / lg_dat$GMS.Musical_training.error ^ 2)
  
  coefs_w4 <- broom::tidy(lm(BAT.score ~ MIQ.score + GMS.musical_training, weights = inv_err, data = lg_dat)) %>% 
    select(term, value = estimate, se = std.error) %>% 
    mutate(term = c("b0", "b1", "b2"),  
           method = "weighting4")
  
  browser()
  ### LAVAAN  
  
  av_se_x <- mean(lg_dat$MIQ.error)
  av_se_y <- mean(lg_dat$BAT.error)
  av_se_z <- mean(.01)
  
  r.lx <- 1 / (1 + av_se_x^2) #compute reliability of latent variable (i.e., common variance) from measurement error
  r.x <- 1 - r.lx    #compute unique variance from measurement error
  
  r.ly <- 1 / (1 + av_se_y^2) #compute reliability of latent variable (i.e., common variance) from measurement error
  r.y <- 1 - r.ly    #compute unique variance from measurement error
  
  r.lz <- 1 / (1 + av_se_z^2) #compute reliability of latent variable (i.e., common variance) from measurement error
  r.z <- 1 - r.lz    #compute unique variance from measurement error
  
  m <- '
        lx =~ 1 * x
        lx ~~ psi * lx
        x ~~ theta * x
        theta == psi * %s/%s

        ly =~ 1 * y
        ly ~~ tau * ly
        y ~~ eta * y
        eta == tau * %s/%s

        lz =~ 1 * z
        lz ~~ phi * lz
        z ~~ zeta * z
        zeta == phi * %s/%s

        y ~ lx + lz
        x ~ 0*1
        z ~ 0*1
        '
  
  m <- sprintf(m, r.x, r.lx, r.y, r.ly, r.z, r.lz)
  
  fit <- suppressWarnings(sem(
    m,
    data = lg_dat %>% select(x = MIQ.score, y = BAT.score, z = GMS.musical_training),
    estimator = "MLR",
    meanstructure = T,
    rstarts = 5,
    optim.dx.tol = 0.001
  ))
  #browser()
  #print(parameterEstimates(fit))
  coefs_lv <- suppressWarnings(parameterEstimates(fit)[c(17, 10, 11), c("label", "est", "se")] %>% as_tibble())
  coefs_lv <- coefs_lv %>%
    select(term = label, value = est, se) %>%
    mutate(term = c("b0", "b1", "b2"), 
           method = "LV")
  
  
  model_simex <- simex(
    m_naive,
    SIMEXvariable = c("MIQ.score", "BAT.score", "GMS.musical_training"),
    B = 200,
    asymptotic = F,
    measurement.error = cbind(lg_dat$MIQ.error, lg_dat$BAT.error, rnorm(nrow(lg_dat), 0, .1)^2)
  )
  sum <- summary(model_simex)$coefficients$jackknife
  coefs_simex <- as.data.frame(sum) %>%
    as_tibble() %>%
    select(value = Estimate, se = `Std. Error`) %>%
    mutate(term = c("b0", "b1", "b2"),
           method = "simex")
  
  mice.dfs = my_TSI(
    lg_dat %>% select(where(is.numeric)) %>% select(-c(test_year, test_wave, age.months)),
    os_names = c("BAT.score", "MIQ.score"),
    se_names = c("BAT.error", "MIQ.error"),
    metrics = "z",
    score_types = "ML",
    separated = T,
    mice_args = c(
      m = 2,
      maxit = 10,
      printFlag = T
    )
  )
  coefs_MI <- pool(with(mice.dfs, lm(true_BAT.score ~ true_MIQ.score + GMS.musical_training))) %>%
    summary() %>%
    select(term, value = estimate, se = std.error) %>%
    mutate(term = c("b0", "b1", "b2"), 
           method = "MI")
  
  coefs_MI

  
  bind_rows(coefs_naive, coefs_w, coefs_MI, coefs_simex, coefs_lv, coefs_outlier)
}
