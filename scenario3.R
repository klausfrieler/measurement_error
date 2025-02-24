## Scenario 3
source("my_TSI.R")

read_item_banks <- function() {
  bat_ib <- read.csv("data/BAT_item_bank.csv")
  bat_ib <- bat_ib %>% rename(difficulty = difficulty_without_track_effect)
  bat_ib <<- bat_ib %>% select(discrimination, difficulty, guessing, inattention) %>% as.matrix()
  
  miq_ib <- read.csv("data/MIQ_item_bank.csv", sep = ";")
  miq_ib <<- miq_ib %>% select(discrimination, difficulty, guessing, inattention) %>% as.matrix()
  
}

##BAT: items, me and reliability found through trial and error with n = 50
### items me    reliability
### 5     1.52  0.3
### 10    1     0.5
### 25    0.65  0.7
### 41    0.49  0.8
### 95    0.33  0.9
### MIQ does not achieve reliability beyond 0.8 due to limited number of items


# simulate_respondents <- function() {
#   foo_25 <- catR::simulateRespondents(
#     thetas = rnorm(50),
#     itemBank = bat_ib,
#     start = start,
#     stop = stop  <- list(
#       rule = "length",
#       thr = 25,
#       alpha = alpha
#     ),
#     test = test,
#     final = final
#   )
#   
#   miq_25 <- catR::simulateRespondents(
#     thetas = rnorm(50),
#     itemBank = bat_ib,
#     start = start,
#     stop = stop  <- list(
#       rule = "length",
#       thr = 25,
#       alpha = alpha
#     ),
#     test = test,
#     final = final
#   )
#   
# }

scenario3_me_levels <- c("low", "medium", "high", "very_high")
my_simulate_respondents <- function(thetas, item_bank, start, stop, test, final) {
  map_dfr(1:length(thetas), function(i) {
    ret <- catR::randomCAT(
      trueTheta = thetas[i],
      itemBank = item_bank,
      start = start,
      stop = stop,
      test = test,
      final = final
    )
    tibble(estimated_theta = ret$thFinal, final_se = ret$seFinal)
  })
}

generate_data_scenario3 <- function(b0 = -1,
                                    b1 = 0.5,
                                    b2 = 0.3,
                                    error_type = c("heteroscedastic"),
                                    me = scenario3_me_levels,
                                    n_sample = 50,
                                    n_batch = 50,
                                    normal = NULL) {
  #the amount of measure error (me) is a parameter, possible values are: 0, 0.33, 0.65, 1, 1.52
  #sample size (n) is a parameter, possible values are 50, 500, 5000
  #include parameter for level of ME info: individual-level vs. variable-level?
  #include parameter for true ME info vs. aggregated/approximate ME info?
  me <- match.arg(me)
  read_item_banks()
  
  bat_me_tab <- tibble(
    no_items = c(5, 10, 25, 41),
    me = c(1.52, 1, 0.65, 0.49),
    rel = c(0.3, 0.5, 0.7, 0.8),
    level = c("very_high", "high", "medium", "low")
  )
  
  miq_me_tab <- tibble(
    no_items = c(2, 5, 25, 95),
    me = c(1.60, 0.98, 0.66, 0.5),
    rel = c(0.28, 0.51, 0.69, 0.8),
    level = c("very_high", "high", "medium", "low")
  )

  ###catR parameters
  alpha <- 0.05
  start <- list(nrItems = 1,
                theta = 0,
                startSelect = "MFI")
  stop  <- list(rule = "length", thr = 10, alpha = alpha)
  final <- list(method = "WL", alpha = alpha)
  test  <- list(
    method = "BM",
    itemSelect = "bOpt",
    priorDist = "norm",
    priorPar = c(0, 1)
  )
  
  n <- n_sample * n_batch
  bat_n_it <- bat_me_tab[bat_me_tab$level == me, ]$no_items
  miq_n_it <- miq_me_tab[miq_me_tab$level == me, ]$no_items
  if (error_type == "heteroscedastic") {
    #BAT scores
    
    #GMS.MT scores
    zt <- scale(runif(n = n, min = 1, max = 7)) %>% as.numeric()
    z_se <- rtruncnorm(
      n = n,
      a = 1,
      b = 7,
      mean = 4,
      #@daniel what to use here?!
      sd = sd(zt) #bat_me_tab$me[i] * sd(xt)
    ) %>% as.numeric()
    z <- zt + z_se
    
    #MIQ scores
    xt <- rnorm(n = n, 0, 1)
    tictoc::tic()
    x_tmp <- my_simulate_respondents(
      thetas = xt ,
      item_bank = miq_ib,
      start = start,
      stop = list(
        rule = "length",
        thr = miq_n_it,
        alpha = alpha
      ),
      test = test,
      final = final
    )
    tictoc::toc(func.toc= function(tic, toc, msg){
      sprintf("MIQ: %.3f elapsed", toc - tic)
    })
    x <- x_tmp$estimated_theta
    x_se <- x_tmp$final_se
    
    yt <- b0 + b1 * xt + b2 * zt
    #yt <- rnorm(n = n_sample, 0, 1)
    tictoc::tic()
    y_tmp <- my_simulate_respondents(
      thetas = yt ,
      item_bank = bat_ib,
      start = start,
      stop = list(
        rule = "length",
        thr = bat_n_it,
        alpha = alpha
      ),
      test = test,
      final = final
    )
    tictoc::toc(func.toc= function(tic, toc, msg){
      sprintf("BAt: %.3f elapsed", toc - tic)
    })
    
    y <- y_tmp$estimated_theta
    y_se <- y_tmp$final_se
    
    df <- tibble(
      yt = yt,
      y = y,
      y_se = y_se,
      xt = xt,
      x = x,
      x_se = x_se,
      zt = zt,
      z = z,
      z_se = z_se
    )
    #colnames(df) <- c("yt", "y", "y_se", "xt", "x", "x_se", "zt", "z", "z_se")
  }
  else{
    stop(sprintf("Invalid error type '%' for CAT", error_type))
  }
  df %>% mutate(batch = rep(1:n_batch, n/n_batch))
}

scenario3_methods <- c("no_correction", "outlier_exclusion", "weighting", "LV", "MI", "simex")

get_coefs_scenario3 <- function(df, method,  measurement_error_level){
  messagef("Scenario3: correcting with method '%s' for error level '%s'", method, measurement_error_level)
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
  
  else if (method == "weighting") {
    #inv_err <- mean(1 / df$y_se^2, 1 / df$x_se^2, 1 / df$z_se^2)
    #inv_err <- rep(mean(1 / df$z_se ^ 2), nrow(df))
    inv_err <- 1 / df$x_se ^ 2
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
          m = 5,
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
  
