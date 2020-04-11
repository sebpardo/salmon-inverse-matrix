### supplementary material
### estimation of marine survival for declining populations
library(tidyverse)
library(cowplot)

source("analysis/01-simulate-data.R")

set.seed(12345)

years <- 1:20
nyears <- length(years)
smolt_n <- rep(10000, nyears) # smolt abundance in year i
S1 <- seq(0.2, 0.02, length.out = nyears) # marine surival 1st winter year i + 1
S2 <- rep(0.4, nyears) # marine surival 1st winter year i + 2
#Pr <- seq(0.8, 0.4, length.out = nyears) # proportion returning as grilse year i + 1
Pr1 <- rep(0.95, nyears) # proportion returning as grilse year i + 1
obs_cv <- 0.05 # coefficient of variation of observation error 
Pr1var <- runif(nyears, 0.6, 0.95)


Prm <- rep(0.4, nyears) # fixed Pr, mixed 1SW-2SW
Prmvar <-  runif(nyears, 0.2, 0.7) # Two additional scenarios with 2SW-dominated populations
Pr2 <- rep(0.15, nyears) # fixed Pr, mixed 2SW dominated
Pr2var <- runif(nyears, 0.05, 0.3)

simuldat1d <- simulate_salmon(years = years, smolt_n = smolt_n, 
                             S1 = S1, S2 = S2, Pr = Pr1, obs_cv = obs_cv) 
simuldat2d <- simulate_salmon(years = years, smolt_n = smolt_n, 
                             S1 = S1, S2 = S2, Pr = Pr1var, obs_cv = obs_cv) 
simuldat3d <- simulate_salmon(years = years, smolt_n = smolt_n, 
                              S1 = S1, S2 = S2, Pr = Prm, obs_cv = obs_cv) 
simuldat4d <- simulate_salmon(years = years, smolt_n = smolt_n, 
                              S1 = S1, S2 = S2, Pr = Prmvar, obs_cv = obs_cv) 
simuldat5d <- simulate_salmon(years = years, smolt_n = smolt_n, 
                              S1 = S1, S2 = S2, Pr = Pr2, obs_cv = obs_cv) 
simuldat6d <- simulate_salmon(years = years, smolt_n = smolt_n, 
                              S1 = S1, S2 = S2, Pr = Pr2var, obs_cv = obs_cv) 

simuldatsd <- bind_rows(mutate(simuldat1d, scenario = 1, pr_type = "Fixed Pr", return_type = "1SW-dominated"),
                       mutate(simuldat2d, scenario = 2, pr_type = "Variable Pr", return_type = "1SW-dominated"),
                       mutate(simuldat3d, scenario = 3, pr_type = "Fixed Pr", return_type = "Mixed 1SW-2SW"),
                       mutate(simuldat4d, scenario = 4, pr_type = "Variable Pr", return_type = "Mixed 1SW-2SW"),
                       mutate(simuldat5d, scenario = 5, pr_type = "Fixed Pr", return_type = "2SW-dominated"),
                       mutate(simuldat6d, scenario = 6, pr_type = "Variable Pr", return_type = "2SW-dominated"))

ann_text <- data.frame(scenario = 1:6, year = 4, lab = letters[1:6],
                       S1_true = 0.003, S1_med = 0.47,
                       pr_type = rep(c("Fixed Pr", "Variable Pr"), 3),
                       return_type = rep(c("1SW-dominated", "Mixed 1SW-2SW", "2SW-dominated"), each = 2))

allsimullongd <- simuldatsd %>%
  select(year, smolt_obs, SW1_obs, SW2_obs, smolt_true, SW1_true, SW2_true, scenario) %>%
  gather(key = lifestage, value = number, -year, -scenario) %>%
  left_join(select(ann_text, scenario, pr_type, return_type), by = "scenario")


alltrued <-  allsimullongd %>%
  filter(lifestage %in% c("smolt_true", "SW1_true", "SW2_true") & year %in% 4:20) %>%
  mutate(lifestage = fct_recode(lifestage, smolt_obs = "smolt_true", SW1_obs = "SW1_true", SW2_obs = "SW2_true")) %>%
  mutate(number = ifelse(lifestage == "smolt_obs", number/10, number))

allsimullongd %>%
  filter(lifestage %in% c("smolt_obs", "SW1_obs", "SW2_obs") & year %in% 4:20) %>%
  mutate(number = ifelse(lifestage == "smolt_obs", number/10, number)) %>%
  ggplot(aes(year, number)) + 
  geom_line(aes(year, number, linetype = lifestage), color = "gray60", data = alltrued, inherit.aes = FALSE) +
  geom_point(aes(year, number, shape = lifestage), size = 1.2, alpha = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Smolt numbers")) +
  facet_wrap(return_type ~ pr_type, ncol = 2, scale = "free_y") + 
  #scale_y_log10() + 
  labs(x = "Year", y = "Return numbers", linetype = "") + 
  scale_shape_discrete(labels  = c("Smolts", "One-sea-winter", "Two-sea-winter")) +
  scale_linetype_discrete(labels  = c("Smolts", "One-sea-winter", "Two-sea-winter")) +
  theme(legend.position = "bottom") +
  labs(color  = "Guide name", linetype = "Guide name", shape = "Guide name")


ggsave("figures/simul_timeseries_declining.pdf", width = 6, height = 7)

### Running SIR
source("analysis/99-sir-functions.R")
source("analysis/99-plot-sir-output.R")

iter <- 30000
iter_samples <- 10000


# add simulation data
popvecs <- simuldat1d %>%
  #rename(small = SW1_obs) %>%
  mutate(smolt = lag(smolt_obs)) %>%
  filter(year > 3) %>% # & year < 18) %>%
  select(year, smolt, SW1_obs, SW2_obs)

# years for simulated data
siryears <- head(popvecs$year, -2) # remove the last year


# Priors from Bayesian model
S1 <- rlnorm(iter, 0.6, 0.3) %>% -. %>% exp # Converting the Z1 prior used in Stan to a distribution of S1s
S2 <- rlnorm(iter, 0.2, 0.3) %>% -. %>% exp 
# Prior for 1SW-dominated rivers (e.g. Conne)
Prg1 <- plogis(rlogis(iter, 2.4, 0.35))
# Wider prior for Pr (non 1SW-dominated rivers)
# Prg <- plogis(rlogis(iter, 0, 0.8))
Prg12 <- plogis(rlogis(iter, -0.5, 0.5))


Prg <- plogis(rlogis(iter, -1, 0.5))


# simuldat1 = fixed Pr, grilse-dominated
out1d <- run_sir(simuldat = simuldat1d, S1 = S1, S2 = S2, Prg = Prg1, constant = 30)
plot_sir_output(siroutput = out1d$siroutput, siryears = siryears, posteriors = out1d$posteriors, 
                simuldat = simuldat1d, S1 = S1, S2 = S2, Prg = Prg1) 
sir_convergence(out1d)

# simuldat2 = variable Pr, grilse-dominated
out2d <- run_sir(simuldat = simuldat2d, S1 = S1, S2 = S2, Prg = Prg1, constant = 70)
plot_sir_output(siroutput = out2d$siroutput, siryears = siryears, posteriors = out2d$posteriors, 
                simuldat = simuldat2d, S1 = S1, S2 = S2, Prg = Prg1) 
sir_convergence(out2d)

# simuldat3 = fixed Pr, 1SW-2SW mix
out3d <- run_sir(simuldat = simuldat3d, S1 = S1, S2 = S2, Prg = Prg, constant = 100)
plot_sir_output(siroutput = out3d$siroutput, siryears = siryears, posteriors = out3d$posteriors, 
                simuldat = simuldat3d, S1 = S1, S2 = S2, Prg = Prg) 
sir_convergence(out3d)


# simuldat4 = variable Pr, 1SW-2SW mix
out4d <- run_sir(simuldat = simuldat4d, S1 = S1, S2 = S2, Prg = Prg12, constant = 70)
plot_sir_output(siroutput = out4d$siroutput, siryears = siryears, posteriors = out4d$posteriors, 
                simuldat = simuldat4d, S1 = S1, S2 = S2, Prg = Prg12) 
sir_convergence(out4d)


# simuldat5 = fixed Pr, 2SW-dominated
out5d <- run_sir(simuldat = simuldat5d, S1 = S1, S2 = S2, Prg = Prg, constant = 160)
plot_sir_output(siroutput = out5d$siroutput, siryears = siryears, posteriors = out5d$posteriors, 
                simuldat = simuldat5d, S1 = S1, S2 = S2, Prg = Prg) 
sir_convergence(out5d)

# simuldat6 = variable Pr, 2SW-dominated
out6d <- run_sir(simuldat = simuldat6d, S1 = S1, S2 = S2, Prg = Prg, constant = 200)
plot_sir_output(siroutput = out6d$siroutput, siryears = siryears, posteriors = out6d$posteriors, 
                simuldat = simuldat6d, S1 = S1, S2 = S2, Prg = Prg) 
sir_convergence(out6d)


bind_cols(
  sir_convergence(out1d),
  sir_convergence(out2d),
  sir_convergence(out3d),
  sir_convergence(out4d),
  sir_convergence(out5d),
  sir_convergence(out6d)
)

save(out1d, out2d, out3d, out4d, out5d, out6d, file = "data/SIR-simul-out-informative-declining.rda")
#save(out1d, out2d, out3d, out4d, out5d, out6d, file = "~/Projects/SIR-simul-out-informative-declining.rda")


allposteriorsd <- bind_rows(mutate(reduce(out1d$posteriors, bind_rows), scenario = 1),
                           mutate(reduce(out2d$posteriors, bind_rows), scenario = 2),
                           mutate(reduce(out3d$posteriors, bind_rows), scenario = 3),
                           mutate(reduce(out4d$posteriors, bind_rows), scenario = 4),
                           mutate(reduce(out5d$posteriors, bind_rows), scenario = 5),
                           mutate(reduce(out6d$posteriors, bind_rows), scenario = 6)) %>% 
  as_tibble()

simuldatsd <- bind_rows(mutate(out1d$simuldat, scenario = 1, pr_type = "Fixed Pr", return_type = "1SW-dominated"),
                       mutate(out2d$simuldat, scenario = 2, pr_type = "Variable Pr", return_type = "1SW-dominated"),
                       mutate(out3d$simuldat, scenario = 3, pr_type = "Fixed Pr", return_type = "Mixed 1SW-2SW"),
                       mutate(out4d$simuldat, scenario = 4, pr_type = "Variable Pr", return_type = "Mixed 1SW-2SW"),
                       mutate(out5d$simuldat, scenario = 5, pr_type = "Fixed Pr", return_type = "2SW-dominated"),
                       mutate(out6d$simuldat, scenario = 6, pr_type = "Variable Pr", return_type = "2SW-dominated"))


median_est_sird <- allposteriorsd %>%
  group_by(scenario, year) %>%
  summarise(S1_med = median(S1),
            S1_1q = quantile(S1, 0.25),
            S1_3q = quantile(S1, 0.75),
            S1_mean = mean(S1),
            S2_med = median(S2),
            S2_1q = quantile(S2, 0.25),
            S2_3q = quantile(S2, 0.75),
            Prg_med = median(Prg),
            Prg_1q = quantile(Prg, 0.25),
            Prg_3q = quantile(Prg, 0.75)) 

# specify annotation text
ann_text <- data.frame(scenario = 1:6, year = 4, lab = letters[1:6],
                       S1_true = 0.003, S1_med = 0.47,
                       pr_type = rep(c("Fixed Pr", "Variable Pr"), 3),
                       return_type = rep(c("1SW-dominated", "Mixed 1SW-2SW", "2SW-dominated"), each = 2))

median_est_trued <- left_join(median_est_sird, simuldatsd, by = c("year", "scenario")) %>%
  ungroup() %>%
  mutate(S1_ratio = S1_med/S1_true,
         scenario = as.character(scenario)) 

median_est_trued %>%
  #mutate(return_type = as.factor(return_type)) %>%
  mutate(return_type = factor(return_type, levels = c("1SW-dominated", "Mixed 1SW-2SW", "2SW-dominated"))) %>%
  ggplot(aes(S1_true, S1_med)) + 
  xlim(0, 0.18) + ylim(0, 0.25) +
  xlab("Simulated S1 value") + ylab("Estimated S1 value") +
  geom_abline(intercept = 0, slope = 1, color = "grey80", linetype = 2) + 
  geom_errorbar(aes(ymin = S1_1q, ymax = S1_3q, width = 0)) +
  # geom_line(aes(group = scenario)) + 
  geom_point(size = 2) + 
  geom_smooth(aes(group = scenario), method = "lm", se = FALSE) +
  facet_grid(return_type ~ pr_type) +
  cowplot::theme_cowplot() + 
  #theme_classic() + 
  cowplot::panel_border(colour = "gray80", size = 0.5, linetype = 1,
                        remove = FALSE) + 
  geom_text(data = ann_text, label = paste0("(",letters[1:6],")"), size = 5)


ggsave("figures/simul_plots_6_declining.png", width = 5.8, height = 6.5)
ggsave("figures/simul_plots_6_declining.pdf", width = 5.8, height = 6.5)


med_est_longd <- gather(allposteriorsd, parameter, value, -year, -weights, -scenario) %>%
  group_by(scenario, year, parameter) %>%
  summarise(med = median(value),
            mean = mean(value),
            sd = sd(value),
            cv = sd/mean,
            q1 = quantile(value, 0.25),
            q3 = quantile(value, 0.75))

simuldats_longd <- simuldatsd %>%
  select(year, S1_true:return_type) %>%
  rename(Prg_true = Pr_true) %>%
  gather(parameter, true_value, -year, -scenario, -pr_type, -return_type) %>%
  mutate(parameter = str_replace(parameter, "_true", ""))


med_est_true_longd <- left_join(med_est_longd, simuldats_longd, by = c("year", "scenario", "parameter")) %>%
  ungroup() %>%
  mutate(scenario = as.character(scenario),
         scenario2 = paste0(pr_type, ",\n", return_type),
         perc_diff = (med - true_value) * 100 / true_value
  )

med_est_true_longd %>%
  mutate(return_type = factor(return_type, levels = c("2SW-dominated", "Mixed 1SW-2SW", "1SW-dominated"))) %>%
  mutate(scenario2 = factor(scenario2)) %>%
  mutate(pr_type = factor(pr_type)) %>%
  mutate(scenario2 = fct_reorder2(scenario2, pr_type, return_type)) %>%
  mutate(parameter = str_replace(parameter, "g", "\\[g\\]"),
         parameter = str_replace(parameter, "([0-9])", "\\[\\1\\]")) %>%
  mutate(parameter = factor(parameter, levels = c("S[1]", "S[2]", "Pr[g]"))) %>%
  ggplot(aes(year, med)) + geom_point() +
  #geom_line(aes(year, true_value), color = "blue") +
  geom_point(aes(year, true_value), color = "blue", size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = q1, ymax = q3, width = 0)) +
  facet_grid(parameter ~ scenario2, scales = "free_y", 
             labeller = labeller(.rows = label_parsed)) +
  theme(strip.text.y = element_text(angle = 0)) +
  xlab("Year") + ylab("Parameter value")

ggsave("figures/sir_out_18_declining.png", width = 10, height = 5.5)
ggsave("figures/sir_out_18_declining.pdf", width = 10, height = 5.5)
