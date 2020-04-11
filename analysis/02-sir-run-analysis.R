### Running SIR for simulated datasets
### Nov 1st 2019

library(tidyverse)

source("analysis/99-sir-functions.R")
source("analysis/99-plot-sir-output.R")
source("analysis/01-simulate-data.R")

iter <- 30000
iter_samples <- 10000
# Estimate S1 and S2, leave everything else constant

#rlnorm(iter, 0.6, 0.3) %>% -. %>% exp %>% hist

# Priors from Bayesian model
# S1 <- rlnorm(iter, 0.6, 0.3) %>% -. %>% exp # Converting the Z1 prior used in Stan to a distribution of S1s
# S2 <- rlnorm(iter, 0.2, 0.3) %>% -. %>% exp 
# Prg <- plogis(rlogis(iter, 0, 0.8))
# 
# # Prior for 1SW-dominated rivers (e.g. Conne)
# Prg <- plogis(rlogis(iter, 2.4, 0.35))

# # Old uninformative priors
# S1 <- runif(iter, 0, 0.9)
# S2 <- runif(iter, 0, 0.9)
# Prg <-  runif(iter, 0.2, 0.99)

# simuldat1 = fixed Pr, grilse-dominated
# simuldat2 = variable Pr, grilse-dominated
# simuldat3 = fixed Pr, 1SW-2SW mix
# simuldat4 = variable Pr, 1SW-2SW mix

# add simulation data
popvecs <- simuldat1 %>%
  #rename(small = SW1_obs) %>%
  mutate(smolt = lag(smolt_obs)) %>%
  filter(year > 3) %>% # & year < 18) %>%
  select(year, smolt, SW1_obs, SW2_obs)

# years for simulated data
siryears <- head(popvecs$year, -2) # remove the last year


##################
### RUNNING SIR
##################


# simuldat1 = fixed Pr, grilse-dominated
# simuldat2 = variable Pr, grilse-dominated
# simuldat3 = fixed Pr, 1SW-2SW mix
# simuldat4 = variable Pr, 1SW-2SW mix
# simuldat5 = fixed Pr, 2SW-dominated
# simuldat6 = variable Pr, 2SW-dominated

# Priors from Bayesian model
S1 <- rlnorm(iter, 0.6, 0.3) %>% -. %>% exp # Converting the Z1 prior used in Stan to a distribution of S1s
S2 <- rlnorm(iter, 0.2, 0.3) %>% -. %>% exp 
# Prior for 1SW-dominated rivers (e.g. Conne)
Prg1 <- plogis(rlogis(iter, 2.4, 0.35))
# Wider prior for Pr (non 1SW-dominated rivers)
# Prg <- plogis(rlogis(iter, 0, 0.8))
Prg12 <- plogis(rlogis(iter, -0.5, 0.5))

Prg <- plogis(rlogis(iter, -1, 0.5))
hist(plogis(rlogis(iter, -0.5, 0.5)))


# simuldat1 = fixed Pr, grilse-dominated
out1 <- run_sir(simuldat = simuldat1, S1 = S1, S2 = S2, Prg = Prg1, constant = 30)
plot_sir_output(siroutput = out1$siroutput, siryears = siryears, posteriors = out1$posteriors, 
                simuldat = simuldat1, S1 = S1, S2 = S2, Prg = Prg1) 
sir_convergence(out1)

# simuldat2 = variable Pr, grilse-dominated
out2 <- run_sir(simuldat = simuldat2, S1 = S1, S2 = S2, Prg = Prg1, constant = 70)
plot_sir_output(siroutput = out2$siroutput, siryears = siryears, posteriors = out2$posteriors, 
                simuldat = simuldat2, S1 = S1, S2 = S2, Prg = Prg1) 
sir_convergence(out2)

# simuldat3 = fixed Pr, 1SW-2SW mix
out3 <- run_sir(simuldat = simuldat3, S1 = S1, S2 = S2, Prg = Prg, constant = 100)
plot_sir_output(siroutput = out3$siroutput, siryears = siryears, posteriors = out3$posteriors, 
                simuldat = simuldat3, S1 = S1, S2 = S2, Prg = Prg) 
sir_convergence(out3)


# simuldat4 = variable Pr, 1SW-2SW mix
out4 <- run_sir(simuldat = simuldat4, S1 = S1, S2 = S2, Prg = Prg12, constant = 70)
plot_sir_output(siroutput = out4$siroutput, siryears = siryears, posteriors = out4$posteriors, 
                simuldat = simuldat4, S1 = S1, S2 = S2, Prg = Prg12) 
sir_convergence(out4)


# simuldat5 = fixed Pr, 2SW-dominated
out5 <- run_sir(simuldat = simuldat5, S1 = S1, S2 = S2, Prg = Prg, constant = 160)
plot_sir_output(siroutput = out5$siroutput, siryears = siryears, posteriors = out5$posteriors, 
                simuldat = simuldat5, S1 = S1, S2 = S2, Prg = Prg) 
sir_convergence(out5)

# simuldat6 = variable Pr, 2SW-dominated
out6 <- run_sir(simuldat = simuldat6, S1 = S1, S2 = S2, Prg = Prg, constant = 200)
plot_sir_output(siroutput = out6$siroutput, siryears = siryears, posteriors = out6$posteriors, 
                simuldat = simuldat6, S1 = S1, S2 = S2, Prg = Prg) 
sir_convergence(out6)

bind_cols(
  sir_convergence(out1),
  sir_convergence(out2),
  sir_convergence(out3),
  sir_convergence(out4),
  sir_convergence(out5),
  sir_convergence(out6)
)


save(out1, out2, out3, out4, out5, out6, file = "data/SIR-simul-out-informative.rda")
#save(out1, out2, out3, out4, out5, out6, file = "~/Projects/SIR-simul-out-informative.rda")


# save(siroutput1, siroutput2, siroutput3, siroutput4,
#      posteriors1, posteriors2, posteriors3, posteriors4,
#      simuldat1, simuldat2, simuldat3, simuldat4, file = "~/Projects/SIR-simul-out.rda")
# 
# save(siroutput5, siroutput6,
#      posteriors5, posteriors6,
#      simuldat5, simuldat6, file = "~/Projects/SIR-simul-out-2SW-dom.rda")

