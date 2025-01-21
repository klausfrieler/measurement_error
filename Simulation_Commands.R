## Commands to run the simulations and analyses

# sample size N = 50
s1_n50 <- ME_simulator$new(scenario = 1, n_sample = 50, b0 = 1, b1 = -0.5, b2 = 0.25)
s1_n50$run()

# sample size N = 500
s1_n500 <- ME_simulator$new(scenario = 1, n_sample = 500, b0 = 1, b1 = -0.5, b2 = 0.25)
s1_n500$run()

# aggregate the 3 sample sizes an create