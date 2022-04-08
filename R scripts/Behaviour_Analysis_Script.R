# Script for analysing flight behaviours of ageing beetles.

# this script focussed on longitudinal measures of flight performance using bayesian
# models. SUrvival analysis was shifted to another script to avoid working with
# multiple dataframe formats

# Libraries ---------------------------------------------------------------
rm(list = ls())

# for data cleaning and plotting
library(plyr)
library(dplyr)
library(lubridate)


# for analysis
library(brms)
library(cmdstanr) # backend for brms
library(broom)
library(loo)      # model selection 
library(ggeffects)
library(glmmTMB)

library(DHARMa)

# load in data ------------------------------------------------------------

#setwd("C:/Users/s1945485/OneDrive - University of Edinburgh/Documents/Edinburgh/My Data/testing")

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

flight <- flight %>%
  filter(Family != "xU") %>%
  filter(Family != "xAP") %>%
  filter(Family != "xAQ") %>%
  filter(Family != "xX") %>%
  filter(Family != "xAS") %>%
  filter(Family != "xAR") %>%
  filter(Family != "xW") %>%
  filter(Family != "xAT") %>%
  filter(Family != "xB") %>%
  filter(Family != "xAM") %>%
  filter(Family != "xAE") %>%
  filter(Family != "xY") %>%
  filter(Family != "xAN") %>%
  filter(Family != "xAB") %>%
  filter(Family != "xAC") %>%
  filter(Family != "xA") %>%
  filter(Family != "xC") %>%
  filter(Family != "xV") %>%
  filter(Family != "xR") %>%
  filter(Family != "xAL")


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
active0$revolution <- round((active0$total_distance/0.88),0)

# Exploring Priors ------------------------------------------------------------
get_prior(bf(
  total_distance ~ sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
  hu ~ sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)),
data = active0,
family = "hurdle_gamma")



distance_brm0_logprior<- brm(
  bf(
  total_distance ~ sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
  hu ~ sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)
  ),
  data = active0,
  family = "hurdle_gamma",
  warmup = 200, 
  iter = 2000, 
  thin = 2,
  control = list(adapt_delta = 0.95),
  prior = c(prior(normal(0,2),class = Intercept),
            prior(normal(0,2),class = b),
            prior(normal(0,2),class = sd),
            prior(logistic(0,2),class = Intercept, dpar = hu),
            prior(normal(0,2),class = b, dpar = hu),
            prior(normal(0,2),class = sd, dpar = hu)),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm0_logprior, "distance_brm0_logprior.RDS")
d_logprior <- loo(distance_brm0_logprior, cores = 4,moment_match = TRUE)

distance_brm0_normprior<- brm(
  bf(
    total_distance ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
    hu ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)
  ),
  data = active0,
  family = "hurdle_gamma",
  warmup = 200, 
  iter = 2000, 
  thin = 2,
  control = list(adapt_delta = 0.95),
  prior = c(prior(normal(0,2),class = Intercept),
            prior(normal(0,2),class = b),
            prior(normal(0,2),class = sd),
            prior(normal(0,2),class = Intercept, dpar = hu),
            prior(normal(0,2),class = b, dpar = hu),
            prior(normal(0,2),class = sd, dpar = hu)),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm0_normprior, "distance_brm0_normprior.RDS")
d_normprior <- loo(distance_brm0_normprior, cores = 4,moment_match = TRUE)
save

distance_brm0_defprior<- brm(
  bf(
    total_distance ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
    hu ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)
  ),
  data = active0,
  family = "hurdle_gamma",
  warmup = 200, 
  iter = 2000, 
  thin = 2,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm0_defprior, "distance_brm0_defprior.RDS")
d_defprior <- loo(distance_brm0_defprior, cores = 4,moment_match = TRUE)

n <- loo_compare(d_defprior,d_normprior,d_logprior)
print(n,moment_match = T, simplify = F)


# BRMS Models for Distance ------------------------------------------------

distance_brm0<- brm(
  bf(
    total_distance ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
    hu ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm0, "distance_brm0.RDS")

distance_brm1a<- brm(bf(
  total_distance ~ Assay + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
  hu ~ Assay + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm1a, "distance_brm1a.RDS")


distance_brm1b<- brm(bf(
  total_distance ~ Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
  hu ~ Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm1b, "distance_brm1b.RDS")


distance_brm2<- brm(bf(
  total_distance ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
  hu ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm2, "distance_brm2.RDS")


distance_brm3<- brm(bf(
  total_distance ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
  hu ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm3, "distance_brm3.RDS")


d_brm0_loo_less_session <- loo(distance_brm0_less_session, cores = 4,moment_match = TRUE)
d_brm1a_loo_less_session <- loo(distance_brm1a_less_session, cores = 4,moment_match = TRUE)
d_brm1b_loo_less_session <- loo(distance_brm1b_less_session, cores = 4,moment_match = TRUE)
d_brm2_loo_less_session <- loo(distance_brm2_less_session, cores = 4,moment_match = TRUE)
d_brm3_loo_less_session <- loo(distance_brm3_less_session, cores = 4,moment_match = TRUE)


n <- loo_compare(d_brm0_loo,d_brm1a_loo,d_brm1b_loo,d_brm2_loo,d_brm3_loo)
print(n,simplify = F)

plot(distance_brm3)

bayesplot::ppc_dens_overlay(y = active0$total_distance, 
                            yrep = posterior_predict(distance_brm3, ndraws = 2000))  +
  ggplot2::xlim(0,20000)

sjPlot::plot_model(distance_brm3, type = "pred", terms = c("Assay","Treatment")) + coord_cartesian(ylim = c(0,8000))

conditional_effects(distance_brm3, robust = F,
                         effects = c("Assay:Treatment"),
                         method = "fitted")
summary(distance_brm3)



# BRMS Models for Peak Speed ----------------------------------------------



peak_brm0<- brm(
  bf(
    total_peak_speed ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
    hu ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber))+ (1|Session),
  data = active0,
  family = "hurdle_gamma",
  warmup = 100, 
  iter = 10000, 
  thin = 10,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm0_less_session, "distance_brm0.RDS")

distance_brm1a_less_session<- brm(bf(
  total_distance ~ Assay + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Assay + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm1a_less_session, "distance_brm1a_less_session.RDS")


distance_brm1b_less_session<- brm(bf(
  total_distance ~ Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm1b_less_session, "distance_brm1b_less_session.RDS")


distance_brm2_less_session<- brm(bf(
  total_distance ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm2_less_session, "distance_brm2_less_session.RDS")


peak_brm3<- brm(bf(
  total_peak_speed ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
  hu ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 200, 
  iter = 2000, 
  thin = 2,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(peak_brm, "peak_brm3.RDS")


d_brm0_loo_less_session <- loo(distance_brm0_less_session, cores = 4,moment_match = TRUE)
d_brm1a_loo_less_session <- loo(distance_brm1a_less_session, cores = 4,moment_match = TRUE)
d_brm1b_loo_less_session <- loo(distance_brm1b_less_session, cores = 4,moment_match = TRUE)
d_brm2_loo_less_session <- loo(distance_brm2_less_session, cores = 4,moment_match = TRUE)
d_brm3_loo_less_session <- loo(distance_brm3_less_session, cores = 4,moment_match = TRUE)



n <- loo_compare(d_brm0_loo,d_brm1a_loo,d_brm1b_loo,d_brm2_loo,d_brm3_loo)
print(n,simplify = F)

plot(peak_brm3)

bayesplot::ppc_dens_overlay(y = active0$total_distance, 
                            yrep = posterior_predict(distance_brm3, ndraws = 2000))  +
  ggplot2::xlim(0,20000)

sjPlot::plot_model(peak_brm3, type = "pred", terms = c("Assay","Treatment"), show.zeroinf = T) 

conditional_effects(peak_brm3, robust = F,
                    effects = c("Assay:Treatment"),
                    method = "fitted")
summary(peak_brm3)



# BRMS Models for Duration ------------------------------------------------


peak_brm0<- brm(
  bf(
    total_peak_speed ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
    hu ~ sel_disappearance + offset(log(trial_length)) +
      (1|ID) + (1|Family) + (1|Block) + (1|Chamber))+ (1|Session),
  data = active0,
  family = "hurdle_gamma",
  warmup = 100, 
  iter = 10000, 
  thin = 10,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm0_less_session, "distance_brm0.RDS")

distance_brm1a_less_session<- brm(bf(
  total_distance ~ Assay + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Assay + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm1a_less_session, "distance_brm1a_less_session.RDS")


distance_brm1b_less_session<- brm(bf(
  total_distance ~ Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm1b_less_session, "distance_brm1b_less_session.RDS")


distance_brm2_less_session<- brm(bf(
  total_distance ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
  hu ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 2500, 
  iter = 25000, 
  thin = 15,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(distance_brm2_less_session, "distance_brm2_less_session.RDS")


duration_brm3<- brm(bf(
  total_duration ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session),
  hu ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
    (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session)),
  data = active0,
  family = "hurdle_gamma",
  warmup = 200, 
  iter = 2000, 
  thin = 2,
  control = list(adapt_delta = 0.95),
  chains = 4,
  cores = 4,
  seed = 123,
  save_pars = save_pars(all = TRUE))
saveRDS(peak_brm, "peak_brm3.RDS")


d_brm0_loo_less_session <- loo(distance_brm0_less_session, cores = 4,moment_match = TRUE)
d_brm1a_loo_less_session <- loo(distance_brm1a_less_session, cores = 4,moment_match = TRUE)
d_brm1b_loo_less_session <- loo(distance_brm1b_less_session, cores = 4,moment_match = TRUE)
d_brm2_loo_less_session <- loo(distance_brm2_less_session, cores = 4,moment_match = TRUE)
d_brm3_loo_less_session <- loo(distance_brm3_less_session, cores = 4,moment_match = TRUE)



n <- loo_compare(d_brm0_loo,d_brm1a_loo,d_brm1b_loo,d_brm2_loo,d_brm3_loo)
print(n,simplify = F)

plot(peak_brm3)

bayesplot::ppc_dens_overlay(y = active0$total_distance, 
                            yrep = posterior_predict(distance_brm3, ndraws = 2000))  +
  ggplot2::xlim(0,20000)

sjPlot::plot_model(duration_brm3, type = "pred", terms = c("Assay","Treatment"), show.zeroinf = T) 

conditional_effects(duration_brm3, robust = F,
                    effects = c("Assay:Treatment"),
                    method = "fitted")
summary(duration_brm3)


# glmmtmb distance models ----------------------------------------------------------
null_model <- glmmTMB(total_distance ~ sel_disappearance + offset(log(trial_length)) +
                        (1|Family/ID) + (1|Block) + (1|Chamber), 
                      zi=~sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                      active0, family=ziGamma(link = "log"))

assay_model_0 <- glmmTMB(total_distance ~ Assay + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         zi=~sel_disappearance +(1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
assay_model_1 <- glmmTMB(total_distance ~ Assay + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         zi=~ Assay + sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))

treatment_model_0 <- glmmTMB(total_distance ~ Treatment + sel_disappearance + offset(log(trial_length)) +
                               (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                             zi=~sel_disappearance +(1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                             active0, family=ziGamma(link = "log"))
treatment_model_1 <- glmmTMB(total_distance ~ Treatment + sel_disappearance + offset(log(trial_length)) +
                               (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                             zi=~Treatment +sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                             active0, family=ziGamma(link = "log"))

add_model_0 <- glmmTMB(total_distance ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
                         (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       zi=~sel_disappearance +(1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       active0, family=ziGamma(link = "log"))
add_model_1 <- glmmTMB(total_distance ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
                         (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       zi=~Assay +sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       active0, family=ziGamma(link = "log"))
add_model_2 <- glmmTMB(total_distance ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
                         (1|ID) +sel_disappearance + (1|Family) + (1|Block) + (1|Chamber), 
                       zi=~Treatment + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       active0, family=ziGamma(link = "log"))
add_model_3 <- glmmTMB(total_distance ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
                         (1|ID) +sel_disappearance + (1|Family) + (1|Block) + (1|Chamber), 
                       zi=~Assay+Treatment + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       active0, family=ziGamma(link = "log"))

inter_model_0 <- glmmTMB(total_distance ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         zi=~sel_disappearance +(1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
inter_model_1 <- glmmTMB(total_distance ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         zi=~Assay +sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
inter_model_2 <- glmmTMB(total_distance ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + sel_disappearance +(1|Family) + (1|Block) + (1|Chamber), 
                         zi=~Treatment + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
inter_model_3 <- glmmTMB(total_distance ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + sel_disappearance +(1|Family) + (1|Block) + (1|Chamber), 
                         zi=~ Assay + Treatment + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
inter_model_4 <- glmmTMB(total_distance ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session), 
                         zi=~ Assay*Treatment + sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session) + offset(log(trial_length)), 
                         active0, family=ziGamma(link = "log"))

anova(null_model, 
      assay_model_0, assay_model_1,
      treatment_model_0, treatment_model_1,
      add_model_0, add_model_1, add_model_2, add_model_3,
      inter_model_0,inter_model_1,inter_model_2,inter_model_3,inter_model_4)

sjPlot::plot_model(inter_model_4, type = "pred", terms = c("Assay","Treatment"), show.zeroinf = T) + coord_cartesian(ylim = c(0,8000))
plot(ggpredict(inter_model_4, c("Assay","Treatment"), type = "zero_inflated")) + coord_cartesian(ylim = c(0,8000))

car::Anova(inter_model_4, type = "III",component = "cond")
car::Anova(inter_model_4, type = "III",component = "zi")

# glmmtmb models ----------------------------------------------------------
null_model <- glmmTMB(total_peak_speed ~ sel_disappearance + offset(log(trial_length)) +
                        (1|Family/ID) + (1|Block) + (1|Chamber), 
                      zi=~sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                      active0, family=ziGamma(link = "log"))

assay_model_0 <- glmmTMB(total_peak_speed ~ Assay + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         zi=~sel_disappearance +(1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
assay_model_1 <- glmmTMB(total_peak_speed ~ Assay + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         zi=~ Assay + sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))

treatment_model_0 <- glmmTMB(total_peak_speed ~ Treatment + sel_disappearance + offset(log(trial_length)) +
                               (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                             zi=~sel_disappearance +(1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                             active0, family=ziGamma(link = "log"))
treatment_model_1 <- glmmTMB(total_peak_speed ~ Treatment + sel_disappearance + offset(log(trial_length)) +
                               (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                             zi=~Treatment +sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                             active0, family=ziGamma(link = "log"))

add_model_0 <- glmmTMB(total_peak_speed ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
                         (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       zi=~sel_disappearance +(1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       active0, family=ziGamma(link = "log"))
add_model_1 <- glmmTMB(total_peak_speed ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
                         (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       zi=~Assay +sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       active0, family=ziGamma(link = "log"))
add_model_2 <- glmmTMB(total_peak_speed ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
                         (1|ID) +sel_disappearance + (1|Family) + (1|Block) + (1|Chamber), 
                       zi=~Treatment + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       active0, family=ziGamma(link = "log"))
add_model_3 <- glmmTMB(total_peak_speed ~ Assay + Treatment + sel_disappearance + offset(log(trial_length)) +
                         (1|ID) +sel_disappearance + (1|Family) + (1|Block) + (1|Chamber), 
                       zi=~Assay+Treatment + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                       active0, family=ziGamma(link = "log"))

inter_model_0 <- glmmTMB(total_peak_speed ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         zi=~sel_disappearance +(1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
inter_model_1 <- glmmTMB(total_peak_speed ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         zi=~Assay +sel_disappearance + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
inter_model_2 <- glmmTMB(total_peak_speed ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + sel_disappearance +(1|Family) + (1|Block) + (1|Chamber), 
                         zi=~Treatment + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
inter_model_3 <- glmmTMB(total_peak_speed ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + sel_disappearance +(1|Family) + (1|Block) + (1|Chamber), 
                         zi=~ Assay + Treatment + (1|ID) + (1|Family) + (1|Block) + (1|Chamber), 
                         active0, family=ziGamma(link = "log"))
inter_model_4 <- glmmTMB(total_peak_speed ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session), 
                         zi=~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) + (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session), 
                         active0, family=ziGamma(link = "log"))


anova(null_model, 
      assay_model_0, assay_model_1,
      treatment_model_0, treatment_model_1,
      add_model_0, add_model_1, add_model_2, add_model_3,
      inter_model_0,inter_model_1,inter_model_2,inter_model_3,inter_model_4)

sjPlot::plot_model(inter_model_4, type = "pred", terms = c("Assay","Treatment"), show.zeroinf = T) +coord_cartesian(ylim = c(0,2.5))
plot(ggpredict(inter_model_4, c("Assay","Treatment"), type = "zero_inflated")) +coord_cartesian(ylim = c(0,2.5))
res<-simulateResiduals(inter_model_4)
plot(res)

car::Anova(inter_model_4, type = "III",component = "cond")
car::Anova(inter_model_4, type = "III",component = "zi")


inter_model_4 <- glmmTMB(total_duration ~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) +
                           (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session), 
                         zi=~ Assay*Treatment + sel_disappearance + offset(log(trial_length)) + (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + (1|Session), 
                         active0, family=ziGamma(link = "log"))
