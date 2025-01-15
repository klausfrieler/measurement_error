generate_data_scenario1 <- function(b0, b1, b2, error_type, me = 1, n_sample = 50, n_batch = 50, n_groups = 10, normal = T){
  n <- n_sample * n_batch
  if(!normal){
    xt <- scale_(runif(
      n = n,
      min = -1,
      max = 1)) 
  } else{
    xt <- rnorm(n, mean = 0, sd = 1)
  }
  
  if(!normal){
    zt <- scale_(18 + rpois(n = n, lambda = 5))
  }
  else{
    stopifnot(n %% n_groups == 0)
    zt <- unlist(replicate(n = n/n_groups, sample(1:n_groups), simplify = FALSE))
    zt <- scale_(zt)
  }
  me <- rtruncnorm(
    n,
    a = 0,
    b = 2 * me,
    mean = me,
    sd = me / 10
  )
  #browser()
  xe <- rnorm(n = n,
              mean = 0,
              sd = me * sd(xt))
  
  if (error_type == "classical") {
    x <- xt + xe
    #print(colnames(df))
  }
  else if (error_type == "systematic") {
    x <- me + (1 + me) * xt + xe
    #browser()
  }
  else if (error_type == "heteroscedastic") {
    xe <- rnorm(n = n,
                mean = 0,
                sd = sqrt(me * abs(xt)^1.3))
    if(any(is.na(xe) | !is.finite(xe))){
      #browser()
    }
    x <- xt + xe
  }
  else if (error_type == "differential") {
    yt <- b0 + b1 * xt + b2 * zt
    xe <- rnorm(n = n,
                mean = 0,
                sd = sqrt(me * abs(yt)^1.3))
    x <- xt + xe
  }
  else{
    stop(sprinf("Unknow error type: %s", error_type))
  }
  #ye <- rnorm(n, 0, 2*sqrt(me))
  yt <- b0 + b1 * xt + b2 * zt
  ye <- rnorm(n, 0, sd = 0.05)
  y <- yt + ye
  tibble(y = y, x = x, xe = xe, xt = xt, zt = zt, batch = rep(1:n_batch, n/n_batch))
  
}

get_coefs_scenario1 <- function(df, method, measurement_error){
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
      df %>% select(-c(xt, xe)),
      os_names = c("x"),
      #se_names = c("xe"),
      metrics = "z",
      score_types = "CTT",
      reliability = c("x" = 1 / (1 + measurement_error ^ 2)),
      separated = F,
      mice_args = c(
        m = 50,
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

