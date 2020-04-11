library(tidyverse)
library(cowplot)
library(xtable)

load("data/SIR-simul-out-informative.rda")  

old.theme <- theme_get()

#theme_set(theme_cowplot())

allposteriors <- bind_rows(mutate(reduce(out1$posteriors, bind_rows), scenario = 1),
                           mutate(reduce(out2$posteriors, bind_rows), scenario = 2),
                           mutate(reduce(out3$posteriors, bind_rows), scenario = 3),
                           mutate(reduce(out4$posteriors, bind_rows), scenario = 4),
                           mutate(reduce(out5$posteriors, bind_rows), scenario = 5),
                           mutate(reduce(out6$posteriors, bind_rows), scenario = 6)) %>% 
  as_tibble()

simuldats <- bind_rows(mutate(out1$simuldat, scenario = 1, pr_type = "Fixed Pr", return_type = "1SW-dominated"),
                       mutate(out2$simuldat, scenario = 2, pr_type = "Variable Pr", return_type = "1SW-dominated"),
                       mutate(out3$simuldat, scenario = 3, pr_type = "Fixed Pr", return_type = "Mixed 1SW-2SW"),
                       mutate(out4$simuldat, scenario = 4, pr_type = "Variable Pr", return_type = "Mixed 1SW-2SW"),
                       mutate(out5$simuldat, scenario = 5, pr_type = "Fixed Pr", return_type = "2SW-dominated"),
                       mutate(out6$simuldat, scenario = 6, pr_type = "Variable Pr", return_type = "2SW-dominated"))


median_est_sir <- allposteriors %>%
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
saveRDS(median_est_sir, file = "data/median_est_sir_6scenarios.rds")

# specify annotation text
ann_text <- data.frame(scenario = 1:6, year = 4, lab = letters[1:6],
                       S1_true = 0.003, S1_med = 0.47,
                       pr_type = rep(c("Fixed Pr", "Variable Pr"), 3),
                       return_type = rep(c("1SW-dominated", "Mixed 1SW-2SW", "2SW-dominated"), each = 2))

median_est_true <- left_join(median_est_sir, simuldats, by = c("year", "scenario")) %>%
  ungroup() %>%
  mutate(S1_ratio = S1_med/S1_true,
         scenario = as.character(scenario)) 

median_est_true %>%
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


ggsave("figures/simul_plots_6.png", width = 5.8, height = 6.5)
ggsave("figures/simul_plots_6.pdf", width = 5.8, height = 6.5)

#####

med_est_long <- gather(allposteriors, parameter, value, -year, -weights, -scenario) %>%
  group_by(scenario, year, parameter) %>%
  summarise(med = median(value),
            mean = mean(value),
            sd = sd(value),
            cv = sd/mean,
            q1 = quantile(value, 0.25),
            q3 = quantile(value, 0.75))
  
simuldats_long <- simuldats %>%
  select(year, S1_true:return_type) %>%
  rename(Prg_true = Pr_true) %>%
  gather(parameter, true_value, -year, -scenario, -pr_type, -return_type) %>%
  mutate(parameter = str_replace(parameter, "_true", ""))


med_est_true_long <- left_join(med_est_long, simuldats_long, by = c("year", "scenario", "parameter")) %>%
  ungroup() %>%
  mutate(scenario = as.character(scenario),
         scenario2 = paste0(pr_type, ",\n", return_type),
         perc_diff = (med - true_value) * 100 / true_value
         )

med_est_true_long %>%
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
  
ggsave("figures/sir_out_18.png", width = 10, height = 5.5)
ggsave("figures/sir_out_18.pdf", width = 10, height = 5.5)

### Create LaTeX table with biases of estimated parameters
outperc <- med_est_true_long %>%
  select(scenario2, everything(), scenario, -q1, -q3, 
         -pr_type, -return_type, -mean, -med, -sd, -true_value, -cv) %>%
  mutate(scn = scenario) %>%
  mutate(scenario = str_replace(scenario2, "\n", " ")) %>%
  select(scenario, parameter, everything()) %>%
  select(-scenario2) %>%
  mutate(scenario = str_replace(scenario, "Pr", "\\Pg")) %>%
  spread(key = parameter, value = perc_diff) %>%
  select(scenario, scn, year, S1, S2, Prg)


perc1 <- filter(outperc, scn == 1) %>% select(-scn)
perc2 <- filter(outperc, scn == 2) %>% select(-scn)
perc3 <- filter(outperc, scn == 3) %>% select(-scn)
perc4 <- filter(outperc, scn == 4) %>% select(-scn)
perc5 <- filter(outperc, scn == 5) %>% select(-scn)
perc6 <- filter(outperc, scn == 6) %>% select(-scn)

allpercs <- bind_cols(select(perc1, -scenario),
          select(perc2, -scenario, -year),
          select(perc3, -scenario, -year),
          select(perc4, -scenario, -year),
          select(perc5, -scenario, -year),
          select(perc6, -scenario, -year))


print.xtable(xtable(allpercs, digits = c(0, 0, rep(1, 18))),
             file = "perc-diff-6scenarios.tex",
             include.rownames = FALSE,
             include.colnames = FALSE,
             only.contents = TRUE)

outperc

outperc %>%
  group_by(scenario) %>%
  summarise(scn = unique(scn)) %>%
  arrange(scn)



# Average percentage difference in each scenario
meanperc <- outperc %>%
  group_by(scn) %>%
  summarise(s1diffmean = mean(S1),
            s2diffmean = mean(S2),
            prdiffmean = mean(Prg))

meanpercrow <- as.matrix(meanperc[,-1]) %>% 
  t %>% as.vector() %>%
  matrix(nrow = 1) %>% as.tibble %>%
  mutate(year = "mean") %>%
  select(year, everything())

print.xtable(xtable(meanpercrow, digits = c(0, 0, rep(1, 18))),
             file = "perc-diff-mean.tex",
             include.rownames = FALSE,
             include.colnames = FALSE,
             only.contents = TRUE,
             hline.after = NULL)

### Create LaTeX table with simulated return abundances in all 6 scenarios
simulreturns <- bind_cols(select(out1$simuldat, year, smolt_obs, SW1_obs, SW2_obs),
          select(out2$simuldat, smolt_obs, SW1_obs, SW2_obs),
          select(out3$simuldat, smolt_obs, SW1_obs, SW2_obs),
          select(out4$simuldat, smolt_obs, SW1_obs, SW2_obs),
          select(out5$simuldat, smolt_obs, SW1_obs, SW2_obs),
          select(out6$simuldat, smolt_obs, SW1_obs, SW2_obs))
print.xtable(xtable(simulreturns[1:20, ], digits = c(0, 0, 0, rep(1, 17))),
             file = "simul-returns.tex",
             include.rownames = FALSE,
             include.colnames = FALSE,
             NA.string = "-",
             only.contents = TRUE)

#### Plot of time series of simulated return abundances in all six scenarios

allsimullong <- simuldats %>%
  select(year, smolt_obs, SW1_obs, SW2_obs, smolt_true, SW1_true, SW2_true, scenario) %>%
  gather(key = lifestage, value = number, -year, -scenario) %>%
  left_join(select(ann_text, scenario, pr_type, return_type), by = "scenario")

#allsimullong$lifestage <- factor(allsimullong$lifestage, levels = c("SW1", "SW2", "smolt"))

alltrue <-  allsimullong %>%
  filter(lifestage %in% c("smolt_true", "SW1_true", "SW2_true") & year %in% 4:20) %>%
  mutate(lifestage = fct_recode(lifestage, smolt_obs = "smolt_true", SW1_obs = "SW1_true", SW2_obs = "SW2_true")) %>%
  mutate(number = ifelse(lifestage == "smolt_obs", number/10, number))
  
allsimullong %>%
  filter(lifestage %in% c("smolt_obs", "SW1_obs", "SW2_obs") & year %in% 4:20) %>%
  mutate(number = ifelse(lifestage == "smolt_obs", number/10, number)) %>%
ggplot(aes(year, number)) + 
  geom_line(aes(year, number, linetype = lifestage), color = "gray60", data = alltrue, inherit.aes = FALSE) +
  geom_point(aes(year, number, shape = lifestage), size = 1.2, alpha = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Smolt numbers")) +
  facet_wrap(return_type ~ pr_type, ncol = 2, scale = "free_y") + 
  #scale_y_log10() + 
  labs(x = "Year", y = "Return numbers", linetype = "") + 
  scale_shape_discrete(labels  = c("Smolts", "One-sea-winter", "Two-sea-winter")) +
  scale_linetype_discrete(labels  = c("Smolts", "One-sea-winter", "Two-sea-winter")) +
  theme(legend.position = "bottom") +
  labs(color  = "Guide name", linetype = "Guide name", shape = "Guide name")

  ggsave("figures/simul_timeseries.pdf", width = 6, height = 7)
  
