## Scenario 3

bat_ib <- read.csv("../BAT_item_bank.csv")
bat_ib <- bat_ib %>% rename(difficulty = difficulty_without_track_effect)
bat_ib <- bat_ib %>% select(discrimination, 
                                  difficulty, 
                                  guessing, 
                                  inattention) %>% as.matrix()

miq_ib <- read.csv("../MIQ_item_bank.csv",sep=";")
miq_ib <- miq_ib %>% select(discrimination, 
                            difficulty, 
                            guessing, 
                            inattention) %>% as.matrix()

alpha <- 0.05
start <- list(nrItems = 1, theta = 0, startSelect = "MFI")
stop  <- list(rule = "length", thr = 10, alpha = alpha)
test  <- list(method = "BM", itemSelect = "bOpt", priorDist = "norm", priorPar = c(0, 1))
final <- list(method = "WL", alpha = alpha)

##BAT: items, me and reliability found through trial and error with n=50
### items me    reliability
### 5     1.52  0.3
### 10    1     0.5
### 25    0.65  0.7
### 41    0.49  0.8
### 95    0.33  0.9
### MIQ does not achieve reliability beyond 0.8 due to limited number of items

bat_me_tab <- data.frame(no_items = c(5,10,25,41), me = c(1.52,1,0.65,0.49), rel = c(0.3,0.5,0.7,0.8), level = c("very_high","high","medium","low"))
miq_me_tab <- data.frame(no_items = c(2,5,25,95), me = c(1.60,0.98,0.66,0.5), rel = c(0.28,0.51,0.69,0.8), level = c("very_high","high","medium","low"))

foo_25 <- simulateRespondents(thetas=rnorm(50), itemBank = bat_ib, start = start, stop = stop  <- list(rule = "length", thr = 25, alpha = alpha), test = test, final = final)
miq_25 <- simulateRespondents(thetas=rnorm(50), itemBank = bat_ib, start = start, stop = stop  <- list(rule = "length", thr = 25, alpha = alpha), test = test, final = final)



sim_s3 <- function(b0 = -1, b1 = 0.5, b2 = 0.3, error_type = c("heteroscedastic"), me = bat_me_tab$level, n_sample = 50, n_batch = 50, methods=c("No_correction", "Outlier_Exclusion","Weighting", "LV","MI")){
  #the amount of measure error (me) is a parameter, possible values are: 0, 0.33, 0.65, 1, 1.52
  #sample size (n) is a parameter, possible values are 50, 500, 5000
  #include parameter for level of ME info: individual-level vs. variable-level?
  #include parameter for true ME info vs. aggregated/approximate ME info?

  for(i in seq(along=me)){
    bat_n_it <- bat_me_tab[bat_me_tab$me==me[i], "no_items"]
    miq_n_it <- miq_me_tab[miq_me_tab$me==me[i], "no_items"] 
      
    if(error_type=="heteroscedastic"){
      #BAT scores
      yt <- rnorm(n = n_sample,0,1)
      y_tmp <- simulateRespondents(thetas= bat_t , itemBank = bat_ib, start = start, stop = list(rule = "length", thr = bat_n_it, alpha = alpha), test = test, final = final)
      y <- bat_tmp$final.values.df$estimated.theta
      y_se <- bat_tmp$final.values.df$final.SE
      
      #GMS.MT scores
      xt <- scale(runif(n = n_sample, min = 1, max = 7))
      xe <- scale(rtruncnorm(n = n_sample, a = 1, b = 7, mean = 4, sd = me*sd(xt)))
      x <- xt + xe
      
      #MIQ scores
      zt <- rnorm(n = n_sample,0,1)
      z_tmp <- simulateRespondents(thetas= miq_t , itemBank = miq_ib, start = start, stop = list(rule = "length", thr = miq_n_it, alpha = alpha), test = test, final = final)
      z <- miq_tmp$final.values.df$estimated.theta
      z_se <- miq_tmp$final.values.df$final.SE
      
      df <- as.data.frame(cbind(yt,y,y_se,xt,x,x_se,zt,z,z_se))
      colnames(df) <- c("yt","y","y_se","xt","x","x_se","zt","z","z_se")
    }
    
    
    
    #analysis methods and evaluation, i.e. comparison of coefficients to ground truth from simulation  
    for(j in seq(along=methods)){
      
      if(methods[j] == "No_correction"){
        coefs <- coef(lm(y ~ x + z, data = df))
        eval_nc <- rbind(eval_nc, eval(coefs=coefs, b0=b0, b1=b1, b2=b2))
      }
      if(methods[j] == "Outlier_Exclusion"){
        ol_y <- suppressMessages(
          boxB(x = df$y,
               k = 1.5,
               method = 'resistant')$outliers)
        
        ol_x <- suppressMessages(
          boxB(x = df$x,
               k = 1.5,
               method = 'resistant')$outliers)
        
        ol_z <- suppressMessages(
          boxB(x = df$z,
               k = 1.5,
               method = 'resistant')$outliers)
        
        
        ol <- c(ol_y, ol_x, ol_z)
        tmp_data <- df    
        if (length(ol) >= 1) {
          #coefs <- coef(lm(y ~ x + z, data = df[-ol, ]))
          tmp_data <- df[-ol,]
        }
        coefs <- broom::tidy(lm(y ~ x + z, data = tmp_data)) %>% 
          select(term, value = estimate, se = std.error) %>% 
          mutate(term = c("b0", "b1", "b2"))
      }
      
      if(methods[j] =="Weighting"){
        inv_err <- mean(1/df$y_se^2, 1/df$x_se^2,1/df$z_se^2)
        coefs <- coef(lm(y ~ x + z, weights=inv_err, data = df))
        eval_w <- rbind(eval_w, eval(coefs=coefs, b0=b0, b1=b1, b2=b2))
      }
      
      if(methods[j] =="LV"){
        av_se_x <- mean(df$x_se)
        av_se_y <- mean(df$y_se)
        av_se_z <- mean(df$z_se)
        
        r.lx <- 1/(1+av_se_x^2) #compute reliability of latent variable (i.e. commn variance) from measurement error
        r.x <- 1 - r.lx    #compute unique variance from measurement error
        
        r.ly <- 1/(1+av_se_y^2) #compute reliability of latent variable (i.e. commn variance) from measurement error
        r.y <- 1 - r.ly    #compute unique variance from measurement error
        
        r.lz <- 1/(1+av_se_z^2) #compute reliability of latent variable (i.e. commn variance) from measurement error
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
      
      if(methods[j] =="MI"){
        tryCatch( 
          {
            mice.dfs=TSI(df %>% select(-c(xt, yt, zt)),
                     os_names=c("y","x","z"),
                     se_names=c("y_se","x_se","z_se"),
                     metrics="z",
                     score_types="ML",
                     separated=T,
                     mice_args=c(
                         m = 50,
                         maxit = 10,
                         printFlag = F))
        
            coefs <- pool(with(mice.dfs, lm(true_y ~ true_x + true_z))) %>% 
            summary() %>% 
            select(term, value = estimate, se = std.error) %>% 
            mutate(term = c("b0", "b1", "b2"))
          },
          error = function(e) {
            coefs <<- tibble(term = c("b0", "b1", "b2"), value = c(NA, NA, NA), se = c(NA, NA, NA))
          }
        )
      }
      
      else if(method == "simex") {
        model_naive <- lm(y ~ x + z, data = df)
        model_simex <- simex(model_naive, SIMEXvariable = c("x","y","z"), 
                             B = 200, asymptotic = F, 
                             measurement.error = cbind(x_se, y_se, z_se))
        sum <- summary(model_simex)$coefficients$jackknife
        coefs <- as.data.frame(sum) %>%  
          as_tibble() %>% 
          select(value = Estimate, se = `Std. Error`) %>% 
          mutate(term = c("b0", "b1", "b2"))
      }
      else{
        stop(sprintf("Unknown method: %s", method))
      }
    }
  }
}




