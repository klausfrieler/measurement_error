library(tidyverse)

mv_poison <-function(n = 100, r = .2, lambda = 1){
  normal_mu <- c(0, 0)
  sigma <- matrix(c(1, r, r, 1), byrow = T, ncol = 2)
  normal <- MASS::mvrnorm(n, normal_mu, sigma)
  unif <- pnorm(normal)
  #browser()
  pois <- t(qpois(t(unif), lambda))
  pois1 <- t(qpois(t(unif[,1]), lambda))
  normal2 <- t(qnorm(t(unif[, 2]))) * 25 + 100
  return(cbind(pois1, normal2))
  #pois  
}

simulate_wy <- function(n = 100, sd_error = 1,  a = 1,  r = -.7, relia = .4, lambda = 1){
  predictors <- MASS::mvrnorm(n = n, 
                               mu = c(0, 0), 
                               Sigma = matrix(c(1, r, r, 1), 
                                              byrow = T, ncol =2)) %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    set_names(c("practice_hours", "vitamin_d"))
  
  predictors <- mv_poison(n = n, r = r, lambda = lambda)%>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    set_names(c("practice_hours", "vitamin_d"))
  
  error <- rnorm(n, 0, sd_error)
  practice_hours_measured <- predictors$practice_hours + rnorm(n, sd = sqrt(1/relia - 1) * lambda)
  
  vitamin_d_measured <- predictors$vitamin_d + rnorm(n, sd = sqrt(1/relia - 1) * lambda)
  #browser()
  #print(var(predictors$practice_hours)/var(practice_hours_measured))
  #browser()
  musical_performance <- a *predictors$practice_hours  + error 
  
  predictors %>% bind_cols(tibble(musical_performance, error)) %>% 
    mutate(practice_hours_measured = practice_hours_measured, 
           vitamin_d_measured = vitamin_d,
           musical_performance_normed = 25 * as.numeric(scale(musical_performance)) + 100)
}
correct_wy_models <- function(n = 1000){
  set.seed(666)  
  simu <- simulate_wy(n = n) 
  browser()
  var_map <- c("y" = "musical_performance", "x" = "practice_hours", "z" = "vitamin_d")
  simu <- simu %>% 
    mutate(y = musical_performance_normed, 
           x = practice_hours_measured, 
           z = vitamin_d_measured) %>% 
    mutate(x_se = .4, y_se = 1, z_se = .4)
  coefs <- 
    map_dfr(setdiff(wy_methods, "MI"), function(m){
      get_coefs_wy(simu, method = m) %>% mutate(method = m)
    
  })
  true <- simu %>% lm(y ~ practice_hours + vitamin_d, data = . )
  true_coefs <- 
    broom::tidy(true) %>%
    select(term, value = estimate, se = std.error) %>%
    mutate(term = c("b0", "b1", "b2"), method = "true_model")
  browser()
  raw_coefs <- coefs
  coefs[coefs$term == "b0",]$value <- (true_coefs[true_coefs$term == "b0",]$value - coefs[coefs$term == "b0",]$value)/true_coefs[true_coefs$term == "b0",]$value  
  coefs[coefs$term == "b1",]$value <- (true_coefs[true_coefs$term == "b1",]$value - coefs[coefs$term == "b1",]$value)/true_coefs[true_coefs$term == "b1",]$value  
  coefs[coefs$term == "b2",]$value <- (true_coefs[true_coefs$term == "b2",]$value - coefs[coefs$term == "b2",]$value)/true_coefs[true_coefs$term == "b2",]$value  
  coefs <- coefs %>% mutate(se = 0)
  coefs 
}

wy_methods <- c("no_correction", "outlier_exclusion", "weighting", "LV", "MI", "simex")

get_coefs_wy <- function(df, method){
  messagef("Westfall-Yarkoni: correcting with method '%s'", method)
  
  if (method == "no_correction") {
    coefs <- broom::tidy(lm(y ~ x + z, data = df)) %>% 
      select(term, value = estimate, se = std.error) %>% 
      mutate(term = c("b0", "b1", "b2"))
  }
  else if (method == "outlier_exclusion") {
    ol_y <- suppressMessages(boxB(
      x = df$y,
      k = 1.5,
      method = 'resistant'
    )$outliers)
    
    ol_x <- suppressMessages(boxB(
      x = df$x,
      k = 1.5,
      method = 'resistant'
    )$outliers)
    
    ol_z <- suppressMessages(boxB(
      x = df$z,
      k = 1.5,
      method = 'resistant'
    )$outliers)
    
    
    ol <- c(ol_y, ol_x, ol_z)
    tmp_data <- df
    if (length(ol) >= 1) {
      #coefs <- coef(lm(y ~ x + z, data = df[-ol, ]))
      tmp_data <- df[-ol, ]
    }
    coefs <- broom::tidy(lm(y ~ x + z, data = tmp_data)) %>%
      select(term, value = estimate, se = std.error) %>%
      mutate(term = c("b0", "b1", "b2"))
  }
  
  # else if (method == "weighting_sd") {
  #   #inv_err <- mean(1 / df$y_se^2, 1 / df$x_se^2, 1 / df$z_se^2)
  #   #inv_err <- rep(mean(1 / df$z_se ^ 2), nrow(df))
  #   #inv_err <- 1 / df$x_se ^ 2
  #   inv_err <- 1 / sqrt(df$x_se ^ 2 + df$z_se ^ 2)
  #   coefs <- broom::tidy(lm(y ~ x + z, weights = inv_err, data = df)) %>% 
  #     select(term, value = estimate, se = std.error) %>% 
  #     mutate(term = c("b0", "b1", "b2"))
  # }
  else if (method == "weighting") {
    #inv_err <- mean(1 / df$y_se^2, 1 / df$x_se^2, 1 / df$z_se^2)
    #inv_err <- rep(mean(1 / df$z_se ^ 2), nrow(df))
    #inv_err <- 1 / df$x_se ^ 2
    inv_err <- 1 / (df$x_se ^ 2 + df$z_se ^ 2)
    coefs <- broom::tidy(lm(y ~ x + z, weights = inv_err, data = df)) %>% 
      select(term, value = estimate, se = std.error) %>% 
      mutate(term = c("b0", "b1", "b2"))
  }
  
  else if (method == "LV") {
    av_se_x <- mean(df$x_se)
    av_se_y <- mean(df$y_se)
    av_se_z <- mean(df$z_se)
    
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
      data = df,
      estimator = "MLR",
      meanstructure = T,
      rstarts = 5,
      optim.dx.tol = 0.001
    ))
    #browser()
    #print(parameterEstimates(fit))
    coefs <- suppressWarnings(parameterEstimates(fit)[c(17, 10, 11), c("label", "est", "se")] %>% as_tibble())
    coefs <- coefs %>%
      select(term = label, value = est, se) %>%
      mutate(term = c("b0", "b1", "b2"))
  }
  
  else if (method == "MI") {
    tryCatch({
      mice.dfs = my_TSI(
        df %>% select(-c(xt, yt, zt)),
        os_names = c("y", "x", "z"),
        se_names = c("y_se", "x_se", "z_se"),
        metrics = "z",
        score_types = "ML",
        separated = T,
        mice_args = c(
          m = 20,
          maxit = 10,
          printFlag = F
        )
      )
      
      coefs <- pool(with(mice.dfs, lm(true_y ~ true_x + true_z))) %>%
        summary() %>%
        select(term, value = estimate, se = std.error) %>%
        mutate(term = c("b0", "b1", "b2")) %>% 
        as_tibble()
    }, error = function(e) {
      warning(sprintf("MSI error for error level"))
      coefs <<- tibble(
        term = c("b0", "b1", "b2"),
        value = c(NA, NA, NA),
        se = c(NA, NA, NA)
      )
    })
  }
  else if (method  == "simex") {
    #browser()
    model_naive <- lm(y ~ x + z, data = df)
    model_simex <- simex(
      model_naive,
      SIMEXvariable = c("x", "y", "z"),
      B = 200,
      asymptotic = F,
      measurement.error = cbind(df$x_se, df$y_se, df$z_se)
    )
    browser()
    sum <- summary(model_simex)$coefficients$jackknife
    coefs <- as.data.frame(sum) %>%
      as_tibble() %>%
      select(value = Estimate, se = `Std. Error`) %>%
      mutate(term = c("b0", "b1", "b2"))
  }
  else{
    stop(sprintf("Unknown method: %s", method))
  }
  #browser()
  coefs 
}


get_wy_plots <- function(n = 1000){
  set.seed(666)  
  simu <- simulate_wy(n = n) 
  mod1 <- simu %>% lm(musical_performance_normed ~ practice_hours_measured, data = . ) 
  mod2 <- simu %>% lm(musical_performance_normed ~ practice_hours, data = . ) 
  ph_mod1 <- simu %>% lm(vitamin_d_measured ~ practice_hours_measured, data = . ) 
  ph_mod2 <- simu %>% lm(vitamin_d_measured ~ practice_hours, data = . ) 
  browser()
  simu <- simu %>% 
    mutate(musical_performance_normed_control_1 = residuals(mod1),
          musical_performance_normed_control_2 = residuals(mod2),
          vitamin_d_control_1 = residuals(ph_mod1),
          vitamin_d_control_2 = residuals(ph_mod2)
    )
  r <- cor(simu$vitamin_d, simu$musical_performance_normed)
  r1 <- cor(simu$vitamin_d, simu$musical_performance_normed_control_1)
  r2 <- cor(simu$vitamin_d, simu$musical_performance_normed_control_2)
  
  simu <- simu %>% 
    pivot_longer(cols =  
                   c("musical_performance_normed", 
                     "musical_performance_normed_control_1", 
                     "musical_performance_normed_control_2")) 
  #browser()  
  q <- simu %>% 
    ggplot(aes(x = vitamin_d_measured, y = value)) 
  name_labs <- c(sprintf("Simple correlation: r = %.2f**", r), 
                 sprintf("Controlled with error: r = %.2f**", r1),
                 sprintf("Controlled without error: r = %.2f (n.s.)", r2))
  names(name_labs) <- c("musical_performance_normed", 
                       "musical_performance_normed_control_1", 
                       "musical_performance_normed_control_2")
  q <- q + geom_point(alpha = .2) 
  q <- q + facet_wrap(~name, scale = "free_y", labeller =  labeller(name = name_labs))  
  q <- q + geom_smooth(method = "lm", color = "indianred", se = F, linewidth = 1.2)
  q <- q + theme_bw(base_size = 12)
  q <- q + theme(strip.background = element_rect(fill = "white"))
  q <- q + labs(x = "Vitamin D [mg]", y = "Musical Performance [a.u.]")
  
  ggsave(q, file = "docs/westfall_yarkoni_remake.png", dpi = 500)
  q  
}

# NOT RUN
# simulate_wy(
#   n = 1000,
#   a = 1,
#   relia = .4,
#   lambda = 1
# ) %>%
#   select(-c(error, practice_hours_measured)) %>%
#   correlation::correlation(partial = T, p_adjust = "none")
