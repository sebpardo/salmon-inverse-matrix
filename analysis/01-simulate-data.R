### Simulating smolt and return abudances with specified S1, S2, and Pr values
### Includes observation error (as the CV of the true obs value) in the simulation process

library(tidyverse)

simulate_salmon <- function(years, smolt_n, S1, S2, Pr, obs_cv) {
  nyears <- length(years)
  if (nyears != length(smolt_n)) stop("Smolt abundance and number of years differ in length")
  if (nyears != length(S1)) stop("S1 and number of years differ in length")
  if (nyears != length(S2)) stop("S2 and number of years differ in length")
  if (nyears != length(Pr)) stop("Pr and number of years differ in length")
 
  # Adding smolt observation error 
  smolt_obs <-  smolt_n + rnorm(nyears, 0, smolt_n * obs_cv)
  
  ### equations
  SW1_true <- smolt_obs * S1 * Pr # grilse abundance year i + 1
  SW2_true <- smolt_obs * S1 * (1 - Pr) * S2 # 2sw abundance year i + 2
  
  # adding observation error
  error_SW1 <- rnorm(nyears, 0, SW1_true * obs_cv)
  SW1_obs <- SW1_true + error_SW1
  SW2_obs <-   SW2_true + rnorm(nyears, 0, SW2_true * obs_cv)
  
  # RS1 <- SW1 * S2 # RS abundance from year i + 1 grilse in year i + 2
  # RS2 <- SW2 * lead(S2) # RS abundance from 2SW fish in year i + 3
  # RS <- RS1 + lag(RS2) # RS total is in year i + 2

  simuldat <- tibble(year = seq_len(length(years) + 2),
                     smolt_true = c(smolt_n, NA, NA),
                     smolt_obs = c(smolt_obs, NA, NA),
                     SW1_obs = c(NA, SW1_obs, NA),
                     SW2_obs = c(NA, NA, SW2_obs),
                     SW1_true = c(NA, SW1_true, NA),
                     SW2_true = c(NA, NA, SW2_true),
                   # RS = c(NA, NA, RS),
                      S1_true = c(NA, S1, NA),
                      S2_true = c(NA, NA, S2),
                      Pr_true = c(NA, Pr, NA))
  simuldat
}

# Setting seed for reproducibility
set.seed(42)

years <- 1:20
nyears <- length(years)

smolt_n <- rep(10000, nyears) # smolt abundance in year i
S1 <- seq(0.02, 0.2, length.out = nyears) # marine surival 1st winter year i + 1
S2 <- rep(0.4, nyears) # marine surival 1st winter year i + 2
#Pr <- seq(0.8, 0.4, length.out = nyears) # proportion returning as grilse year i + 1
Pr <- rep(0.95, nyears) # proportion returning as grilse year i + 1
obs_cv <- 0.05 # coefficient of variation of observation error 

# Need to create 6 scenarios, with all combinations of:
# fixed and variable Pr, 1SW dominated, mixed 1SW-2SW, and 2SW dominated
# and assess how estimates are biased

simuldat1 <- simulate_salmon(years = years, smolt_n = smolt_n, 
                             S1 = S1, S2 = S2, Pr = Pr, obs_cv = obs_cv) # fixed Pr, grilse-dominated
simuldat2 <- simulate_salmon(years = years, smolt_n = smolt_n, 
                             S1 = S1, S2 = S2, Pr = runif(nyears, 0.6, 0.95), obs_cv = obs_cv) # variable Pr, grilse-dominated

simuldat3 <- simulate_salmon(years = years, smolt_n = smolt_n, 
                             S1 = S1, S2 = S2, Pr = rep(0.4, nyears), obs_cv = obs_cv) # fixed Pr, mixed 1SW-2SW
simuldat4 <- simulate_salmon(years = years, smolt_n = smolt_n, 
                             S1 = S1, S2 = S2, Pr = runif(nyears, 0.2, 0.7), obs_cv = obs_cv) 

# Two additional scenarios with 2SW-dominated populations
simuldat5 <-  simulate_salmon(years = years, smolt_n = smolt_n, 
                              S1 = S1, S2 = S2, Pr = rep(0.15, nyears), obs_cv = obs_cv) # fixed Pr, mixed 2SW dominated
simuldat6 <- simulate_salmon(years = years, smolt_n = smolt_n, 
                              S1 = S1, S2 = S2, Pr = runif(nyears, 0.05, 0.3), obs_cv = obs_cv) # variable Pr, mixed 2SW dominated

save(simuldat1, simuldat2, simuldat3, simuldat4, simuldat5, simuldat6, file = "data/simulated-data.rda")
