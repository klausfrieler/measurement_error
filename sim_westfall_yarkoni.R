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
  #print(var(predictors$practice_hours)/var(practice_hours_measured))
  #browser()
  musical_performance <- a *predictors$practice_hours  + error 
  
  predictors %>% bind_cols(tibble(musical_performance, error)) %>% 
    mutate(practice_hours_measured = practice_hours_measured, 
           musical_performance_normed = 25 * as.numeric(scale(musical_performance)) + 100)
}


get_wy_plots <- function(n = 1000){
  simu <- simulate_wy(n = n) 
  mod1 <- simu %>% lm(musical_performance_normed ~ practice_hours_measured, data = . ) 
  mod2 <- simu %>% lm(musical_performance_normed ~ practice_hours, data = . ) 
  ph_mod1 <- simu %>% lm(vitamin_d ~ practice_hours_measured, data = . ) 
  ph_mod2 <- simu %>% lm(vitamin_d ~ practice_hours, data = . ) 
  
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
    ggplot(aes(x = vitamin_d, y = value)) 
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
