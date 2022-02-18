
# Libraries ---------------------------------------------------------------
rm(list = ls())

# for data cleaning and plotting
library(tidyverse)
library(lubridate)


# for analysis
library(brms)
library(cmdstanr) # backend for brms
library(broom)
library(loo)      # model selection 
library(ggeffects)

library(survival)
library(survminer)

# load in data ------------------------------------------------------------

# read "Combined_Flight_final.csv" from "Final_Values_Folder. Be careful about
# Working directories and where new files are being written.
flight<- read.csv(file.choose(), stringsAsFactors = FALSE, header = T)

# tidying up some factors etc...
flight <- flight %>%
  mutate(Family = factor(Family)) %>%
  mutate(ID = factor(ID)) %>%
  mutate(Treatment = factor(Treatment)) %>%
  mutate(Sex = factor(Sex)) %>%
  mutate(Activity = factor(Activity)) %>%
  mutate(Block = factor(Block)) %>%
  mutate(Assay = factor(Assay)) %>%
  mutate(Chamber = factor(Chamber)) %>%
  mutate(Session = factor(Session)) %>%
  mutate(Size = as.numeric(Size)) %>%
  mutate(total_peak_speed = as.numeric(total_peak_speed)) %>%
  mutate(sel_disappearance = case_when(
    Lifespan < 60 ~ "1",
    Lifespan < 90 ~ "2",
    Lifespan >= 90 ~ "3")) %>%
  mutate(sel_disappearance = factor(sel_disappearance)) %>%
  mutate(rotations = (total_distance/0.88))


# Creating a subset to work with only beetles that flew. Removing missing values.
active0 <- flight %>%
  filter(Activity == "Flight") %>%
  filter(!is.na(Size)) %>%
  filter(age_at_flight > 0) %>%
  mutate(fly.binary = ifelse(total_distance > 0, 1,0))

active0[is.na(active0)] <- 0

active <- active0 %>%
  filter(total_distance > 0)

# checking the proportion of zeros in measure of distance
100*sum(active0$total_distance == 0)/nrow(active0)

# brms Models ------------------------------------------------------------

distance_brm0<- brm(bf(
  total_distance ~ sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = hurdle_gamma(),
  warmup = 1000, 
  iter = 7500, 
  thin = 4,
  control = list(adapt_delta = 0.95),
  chains = 6,
  cores = 6,
  seed = 42,
  save_all_pars = TRUE)


distance_brm1a<- brm(bf(
  total_distance ~ Assay + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Assay + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = hurdle_gamma(),
  warmup = 1000, 
  iter = 7500, 
  thin = 4,
  control = list(adapt_delta = 0.95),
  chains = 6,
  cores = 6,
  seed = 42,
  save_all_pars = TRUE)


distance_brm1b<- brm(bf(
  total_distance ~ Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Assay + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = hurdle_gamma(),
  warmup = 1000, 
  iter = 7500, 
  thin = 4,
  control = list(adapt_delta = 0.95),
  chains = 6,
  cores = 6,
  seed = 42,
  save_all_pars = TRUE)


distance_brm2<- brm(bf(
  total_distance ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = hurdle_gamma(),
  warmup = 1000, 
  iter = 7500, 
  thin = 4,
  control = list(adapt_delta = 0.95),
  chains = 6,
  cores = 6,
  seed = 42,
  save_all_pars = TRUE)


distance_brm3<- brm(bf(
  total_distance ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = hurdle_gamma(),
  warmup = 1000, 
  iter = 7500, 
  thin = 4,
  control = list(adapt_delta = 0.95),
  chains = 6,
  cores = 6,
  seed = 42,
  save_all_pars = TRUE)


d_brm0_loo <- loo(distance_brm0, cores = 6)
d_brm1a_loo <- loo(distance_brm1a, cores = 6)
d_brm1b_loo <- loo(distance_brm1b, cores = 6)
d_brm2_loo <- loo(distance_brm2, cores = 6)
d_brm3_loo <- loo(distance_brm3, cores = 6)

plot(distance_brm0)
pp_check(distance_brm3, ndraws = 100) + ggplot2::xlim(0,10000)

plot(d_brm3_loo)

n <- loo_compare(d_brm0_loo,d_brm1a_loo,d_brm1b_loo,d_brm2_loo,d_brm3_loo)
print(n,moment_match = T, simplify = F)

# Survival Models ---------------------------------------------------------


# Be careful about Working directories and where new files are being written.
lifedata<- read.csv(file.choose(), stringsAsFactors = FALSE, header = T)

# Create more informative ID variable. Remove junk/pilot blocks. Create ordered 
# factor of Treatment for models. Remove unnecessary columns
lifedata <- lifedata %>%
  mutate(ID = paste0(Family,ID)) %>%
  filter(Block != "A" & Block != "E") %>%
  mutate(oTreatment = ordered(Treatment, levels = c('A', 'B'))) %>%
  mutate(Family = factor(Family)) %>%
  mutate(ID = factor(ID)) %>%
  mutate(Block = factor(Block)) %>%
  select(Family,ID,Sex,oTreatment,Size,Censored,Lifespan,Activity,Block)

# Convert individual measuers to long-form "state per day/age" format. MUST use 
# Lifespan as response variable in models. I made the mistake with timegroup before.
longlife <- survSplit(Surv(Lifespan,Censored) ~., 
                      lifedata,
                      cut=c(unique(lifedata$Lifespan)), 
                      episode ="timegroup")

# look only at individuals who flew.
longlife_flyers <- longlife[longlife$Activity == "Flight",]

system.time(surv0 <- brm(Censored ~ s(Lifespan, k = 20) + (1|ID) + (1|Family) + (1|Block),
                             family = bernoulli(link = "cloglog"), 
                             data = longlife_flyers, 
                             warmup = 1000, 
                             iter = 7500, 
                             thin = 4,
                             control = list(adapt_delta = 0.95),
                             chains = 4,
                             cores = 4,
                             seed = 95, 
                             backend = "cmdstanr", 
                             threads = threading(2))
)

system.time(surv1 <- brm(Censored ~ oTreatment + s(Lifespan, k = 20) + (1|ID) + (1|Family) + (1|Block),
                             family = bernoulli(link = "cloglog"), 
                             data = longlife_flyers, 
                             warmup = 1000, 
                             iter = 7500, 
                             thin = 4,
                             control = list(adapt_delta = 0.95),
                             chains = 4,
                             cores = 4,
                             seed = 95, 
                             backend = "cmdstanr", 
                             threads = threading(2))
)


system.time(surv2 <- brm(Censored ~ oTreatment + s(Lifespan, k = 20) + s(Lifespan, by = oTreatment) + (1|ID) + (1|Family) + (1|Block),
  family = bernoulli(link = "cloglog"), 
  data = longlife_flyers, 
  warmup = 1000, 
  iter = 7500, 
  thin = 4,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 95, 
  backend = "cmdstanr", 
  threads = threading(2))
)

saveRDS(surv0, file = "surv0.rds")
saveRDS(surv1, file = "surv0.rds")
saveRDS(surv2, file = "surv0.rds")

plot(conditional_smooths(surv0))
plot(msms)
pp_check(test_surv)
brms::pp_check(test_surv, type = "bars", ndraws = 100)

y <- ggemmeans(
  test_surv,
  terms = c("Lifespan", "oTreatment"),
  ci.lvl = 0.95,
  type = "fe",
  back.transform = FALSE)
plot(y)

ggplot(y, aes(x = x, y = log(predicted), col = group)) + 
  geom_smooth(method = "loess", se = F) +
  geom_ribbon(data = y, aes(ymin = log(conf.low), ymax = log(conf.high), alpha = 0.01, fill = group)) +
  facet_wrap(~group)






