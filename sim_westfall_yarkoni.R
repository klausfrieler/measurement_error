library(tidyverse)

get_coefs_wy <- function(df, method, measurement_error){
  if (method == "no_correction") {
    coefs <- broom::tidy(lm(y ~ x + zt, data = df)) %>% 
      select(term, value = estimate, se = std.error) %>% 
      mutate(term = c("b0", "b1", "b2"))
  }
  else if (method == "outlier_exclusion") {
    ol <- suppressMessages(
      boxB(x = df$x,
           k = 1.5,
           method = 'resistant')$outliers)
    tmp_data <- df
    if (length(ol) >= 1) {
      #coefs <- coef(lm(y ~ x + z, data = df[-ol, ]))
      tmp_data <- df[-ol,]
    }
    coefs <- broom::tidy(lm(y ~ x + zt, data = tmp_data)) %>% 
      select(term, value = estimate, se = std.error) %>% 
      mutate(term = c("b0", "b1", "b2"))
  }
  
  # if (methods[j] == "weighting") {
  #   inv_err <- rep(mean(1 / df$xe ^ 2), nrow(df))
  #   coefs <- coef(lm(yt ~ x + zt, weights = inv_err, data = df))
  #   eval_w <- rbind(eval_w, eval(
  #     coefs = coefs,
  #     b0 = b0,
  #     b1 = b1,
  #     b2 = b2
  #   ))
  # }
  
  else if (method == "LV") {
    r.lx <- 1 / (1 + measurement_error ^ 2) #compute reliability of latent variable (i.e. common variance) from measurement error
    r.x <- 1 - r.lx    #compute unique variance from measurement error
    
    m <- '
        lx =~ 1 * x
        lx ~~ psi * lx
        x ~~ theta * x
        theta == psi * %s/%s
        y ~ lx + zt
        '
    # theta == psi * r.x/r.lx => error variance / factor variance
    m <- sprintf(m, r.x, r.lx)
    
    fit <- sem(
      m,
      data = df,
      estimator = "MLR",
      meanstructure = T
    )
    
    coefs <- parameterEstimates(fit)[c(9, 4, 5), c("label", "est", "se")] %>% as_tibble()
    coefs <- coefs %>% 
      select(term = label, value = est, se) %>% 
      mutate(term = c("b0", "b1", "b2"))
    
  }
  else if (method == "MI") {
    #browser()
    mice.dfs = TSI(
      df,
      os_names = c("x"),
      #se_names = c("xe"),
      metrics = "z",
      score_types = "CTT",
      reliability = c("x" = 1 / (1 + measurement_error ^ 2)),
      separated = F,
      mice_args = c(
        m = 20,
        maxit = 10,
        printFlag = F
      )
    )
    #coefs <- pool(with(mice.dfs, lm(y ~ true_x + zt))) %>% pluck("pooled") %>% pluck("estimate")
    coefs <- pool(with(mice.dfs, lm(y ~ true_x + zt))) %>%
      summary() %>%
      select(term, value = estimate, se = std.error) %>%
      mutate(term = c("b0", "b1", "b2"))
  }
  
  else if(method == "simex") {
    model_naive <- lm(y ~ x + zt, data = df)
    model_simex <- simex(model_naive, SIMEXvariable = "x", B = 200, 
                         asymptotic = F, 
                         measurement.error = measurement_error)
    sum <- summary(model_simex)$coefficients$jackknife
    coefs <- as.data.frame(sum) %>%  
      as_tibble() %>% 
      select(value = Estimate, se = `Std. Error`) %>% 
      mutate(term = c("b0", "b1", "b2"))
  }
  else{
    stop(sprintf("Unknown method: %s", method))
  }
  coefs
}

MTS2 <- c(
  +1.000, +0.681, +0.609, +0.557, -0.098, -0.125,				
  +0.681,	+1.000, +0.667, +0.501, -0.132, -0.174,  			
  +0.609, +0.667,	+1.000, +0.417, -0.331, -0.394,		
  +0.557,	+0.501,	+0.417,	+1.000, +0.003, -0.021,
  -0.098,	-0.132,	-0.331,	+0.003,	+1.000, +0.927,
  -0.125,	-0.174,	-0.394,	-0.021,	+0.927,	+1.000)

mat_from_lower_tri <- function(x, with_diag = T, names = NULL){
  n <- floor(sqrt(2 * length(x)))
  mat <- matrix(rep(0, n*n), ncol = n)
  upper_indices <- which(upper.tri(mat, diag = TRUE), arr.ind = TRUE)  
  mat[upper_indices] <- x
  browser()
  ret <- mat %>%  t() %>% as.dist() %>% as.matrix()
  diag(ret) <- 1
  if(!is.null(names)){
    colnames(ret) <- names
    row.names(ret) <- names
  }
  ret
}
mts2_names <- c("OC",
"MU",
"JB",
"F",
"A",
"PR")
 
mts3_names <-c("A", "F",  "MU", "JB", "OC", "PR")
mts4_names <-c("MU", "JB", "F", "PR", "OC",  "A" )

MTS3 <- c(1,				
0.35, 1,				
0.33, 0.017, 1, 		
0.073, -0.172,	-0.303, 1, 	
-0.155, -0.337, -0.119, 0.102, 1,
-0.796, -0.343, -0.224, -0.086,	0.263,	1)

MTS4 <- c(
  1, 				
  0.643,	1,			
  0.386,	0.251,	1, 		
  0.145,	0.115,	0.007,	1, 	
  -0.192,	-0.042,	-0.294,	-0.005,	1, 
  -0.32,	-0.337,	-0.264,	-0.297,	0.062,	1
  )
  
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
  #var_map <- c("y" = "musical_performance", "x" = "practice_hours", "z" = "vitamin_d")
  simu <- simu %>% 
    mutate(y = musical_performance_normed, 
           x = practice_hours_measured, 
           zt = vitamin_d) %>% 
    mutate(x_se = .4, y_se = 1)
  coefs <- 
    map_dfr(setdiff(wy_methods, "MI"), function(m){
      get_coefs_wy(simu, method = m, measurement_error = .4) %>% 
        mutate(method = m)
    
  })
  true <- simu %>% lm(musical_performance ~ practice_hours + vitamin_d, data = . )
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

wy_methods <- c("no_correction", "outlier_exclusion",  "LV", "MI", "simex")

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

coef_plot <- function(coefs_data){
  coefs_data <- coefs_data %>% 
    mutate(method = factor(method, 
                           levels = c("no_correction", 
                                      "weighting", 
                                      "outlier_exclusion",
                                      "simex",
                                      "brms", 
                                      "MI", 
                                      "LV"))) 
  q <- coefs_data %>% ggplot(aes(x = term, y = value, fill = method))
  q <- q + geom_col(position = position_dodge(width = 1))
  q <- q + geom_point(position = position_dodge(width = 1))
  q <- q + geom_errorbar(aes(ymin = value - 2*se, ymax = value + 2*se), 
                         position = position_dodge(width = 1), width = .7 )
  q <- q + scale_fill_brewer(palette = "Set1")
  q <- q + labs(x = "Cofficient", y = "Value", fill = "")
  q <- q + theme_bw()
  q <- q + theme(legend.position = c(.6, .25), legend.background = element_blank())
  q <- q + theme(legend.position = "bottom")
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
#
# Correction:
# correct_wy_models() %>% coef_plot() 