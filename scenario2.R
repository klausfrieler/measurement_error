generate_data_scenario2 <- function(b0, b1, b2, error_type, me = 1, n_sample = 50, n_batch = 50, normal = F){
  n <- n_sample * n_batch
  #browser()
  if(!normal){
    xt <- scale_(5+rpois(
      n = n,
      lambda = 2))
    xt <- scale_(xt)
  } else{
    xt <- rnorm(n, mean = 0, sd = 1)
  }
  
  if(!normal){
    zt <- xt + runif(n = n, min = 1, max = 10) 
    zt <- scale_(zt)
  }
  else{
    zt <- rnorm(n, 0, 5) #xt + 
    zt <- scale_(zt)
  }
  #me <- rtruncnorm(
  #  n,
  #  a = 0,
  #  b = 2 * me,
  #  mean = me,
  #  sd = me / 10
  #)
  
  xe <- rnorm(n = n,
              mean = 0,
              sd = me * sd(xt))
  
  ze <- rnorm(n = n,
              mean = 0,
              sd = me * sd(zt))
  
  if (error_type == "classical") {
    x <- xt + xe
    z <- zt + ze
  }
  else if (error_type == "systematic") {
    x <- me + (1 + me) * xt + xe
    x <- scale_(x)
    z <- me + (1 + me) * zt + ze
    z <- scale_(z)
  }
  else if (error_type == "heteroscedastic") {
    xe <- rnorm(n = n,
                mean = 0,
                sd = sqrt(me * abs(xt)^1.3))
    if(any(is.na(xe) | !is.finite(xe))){
      #browser()
    }
    x <- xt + xe
    
    ze <- rnorm(n = n,
                mean = 0,
                sd = sqrt(me * rev(abs(zt))^1.3))
    if(any(is.na(ze) | !is.finite(ze))){
      #browser()
    }
    z <- zt + ze
  }
  else if (error_type == "differential") {
    yt <- b0 + b1 * xt + b2 * zt
    xe <- rnorm(n = n,
                mean = 0,
                sd = sqrt(me * abs(yt)^1.3))
    x <- xt + xe
    
    ze <- rnorm(n = n,
                mean = 0,
                sd = sqrt(me * rev(abs(yt))^1.3))
    z <- zt + ze
  }
  else{
    stop(sprinf("Unknow error type: %s", error_type))
  }
  #ye <- rnorm(n, 0, 2*sqrt(me))
  yt <- b0 + b1 * xt + b2 * zt
  ye <- rnorm(n, 0, sd = 0.05)
  y <- yt + ye
  tibble(y = y, x = x, xe = xe, xt = xt, zt = zt, ze = ze, z = z, batch = rep(1:n_batch, n/n_batch))
  
}


get_coefs_scenario2 <- function(df, method, measurement_error){
  if (method == "no_correction") {
    coefs <- broom::tidy(lm(y ~ x + z, data = df)) %>% 
      select(term, value = estimate, se = std.error) %>% 
      mutate(term = c("b0", "b1", "b2"))
    #coefs <- coef(lm(y ~ x + z, data = df))
  }
  else if (method == "outlier_exclusion") {
    ol_x <- suppressMessages(
      boxB(x = df$x,
           k = 1.5,
           method = 'resistant')$outliers)
    
    ol_z <- suppressMessages(
      boxB(x = df$z,
           k = 1.5,
           method = 'resistant')$outliers)
    ol <- c(ol_x, ol_z)
    tmp_data <- df    
    if (length(ol) >= 1) {
      #coefs <- coef(lm(y ~ x + z, data = df[-ol, ]))
      tmp_data <- df[-ol,]
    }
    coefs <- broom::tidy(lm(y ~ x + z, data = tmp_data)) %>% 
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
    
    r.lz <- 1 / (1 + measurement_error ^ 2) #compute reliability of latent variable (i.e. common variance) from measurement error
    r.z <- 1 - r.lz    #compute unique variance from measurement error
    
    m <- '
        lx =~ 1 * x
        lx ~~ psi * lx
        x ~~ theta * x
        theta == psi * %s/%s
        
        lz =~ 1 * z
        lz ~~ phi * lz
        z ~~ zeta * z
        zeta == phi * %s/%s
        
        y ~ lx + lz
        x ~ 0*1
        z ~ 0*1
        '
    # theta == psi * r.x/r.lx => error variance / factor variance
    m <- sprintf(m, r.x, r.lx, r.z, r.lz)
    
    fit <- sem(
      m,
      data = df,
      estimator = "MLR",
      meanstructure = T,
      rstarts=5,
      optim.dx.tol = 0.001
    )
    #browser()
    #print(parameterEstimates(fit))
    coefs <- parameterEstimates(fit)[c(13, 7, 8), c("label", "est", "se")] %>% as_tibble()
    coefs <- coefs %>% 
      select(term = label, value = est, se) %>% 
      mutate(term = c("b0", "b1", "b2"))
  }
  else if (method == "MI") {
    #browser()
    tryCatch( 
      {
      mice.dfs = TSI(
        df %>% select(-c(xt, xe, zt, ze)),
        os_names = c("x","z"),
        #se_names = c("xe"),
        metrics = "z",
        score_types = "CTT",
        reliability = c(
          "x" = 1 / (1 + measurement_error ^ 2),
          "z" = 1 / (1 + measurement_error ^ 2)
        ),
        separated = F,
        mice_args = c(
          m = 50,
          maxit = 10,
          printFlag = F
        )
      )
      #coefs <- pool(with(mice.dfs, lm(y ~ true_x + true_z))) %>% pluck("pooled") %>% pluck("estimate")
      coefs <- pool(with(mice.dfs, lm(y ~ true_x + true_z))) %>% 
        summary() %>% 
        select(term, value = estimate, se = std.error) %>% 
        mutate(term = c("b0", "b1", "b2"))
      #browser()
      },
      error = function(e) {
        coefs <<- tibble(term = c("b0", "b1", "b2"), value = c(NA, NA, NA), se = c(NA, NA, NA))
      }
   )
  }
  
  else if(method == "simex") {
    model_naive <- lm(y ~ x + z, data = df)
    model_simex <- simex(model_naive, SIMEXvariable = c("x","z"), 
                         B = 200, asymptotic = F, 
                         measurement.error = cbind(measurement_error, measurement_error))
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
