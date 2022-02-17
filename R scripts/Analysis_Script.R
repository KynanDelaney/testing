
# Libraries ---------------------------------------------------------------
rm(list = ls())

# for data cleaning and plotting
library(tidyverse)
library(lubridate)


# for analysis
library(brms)
library(cmdstanr) # backend for brms
library(loo)      # model selection 
library(glmmTMB)
library(gamm4)

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
    Lifespan > 90 ~ "3")) %>%
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



# glmmTMB Models ------------------------------------------------------------

# checking the proportion of zeros in measure of distance
100*sum(active0$total_distance == 0)/nrow(active0)


hurdle_distance_0 <- glmmTMB(total_distance ~ sel_disappearance + 
                               (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + offset(log(trial_length)),
                         data=active0,
                         ziformula=~1,
                         family = ziGamma(link = "log"))

summary(hurdle_distance_0)

hurdle_distance_1a <- glmmTMB(total_distance ~ Treatment + sel_disappearance + 
                                (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + offset(log(trial_length)),
                             data=active0,
                             ziformula=~1,
                             family = ziGamma(link = "log"))

summary(hurdle_distance_1a)


hurdle_distance_1b <- glmmTMB(total_distance ~ Assay + sel_disappearance + 
                                (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + offset(log(trial_length)),
                              data=active0,
                              ziformula=~1,
                              family = ziGamma(link = "log"))

summary(hurdle_distance_1b)

hurdle_distance_2 <- glmmTMB(total_distance ~ Assay + Treatment + sel_disappearance + 
                               (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + offset(log(trial_length)),
                             data=active0,
                             ziformula=~1,
                             family = ziGamma(link = "log"))

summary(hurdle_distance_2)

hurdle_distance_3_0 <- glmmTMB(total_distance ~ Assay*Treatment + sel_disappearance + 
                               (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + offset(log(trial_length)),
                             data=active0,
                             ziformula=~1,
                             family = ziGamma(link = "log"))

hurdle_distance_3_1 <- glmmTMB(total_distance ~ Assay*Treatment + sel_disappearance + 
                                 (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + offset(log(trial_length)),
                               data=active0,
                               ziformula=~Assay + Treatment + sel_disappearance + 
                                 (1|ID) + (1|Block) + (1|Chamber),
                               family = ziGamma(link = "log"))

summary(hurdle_distance_3)

anova(hurdle_distance_3_0, hurdle_distance_3_1)

distance_brm0<- brm(bf(total_distance ~ sel_disappearance + 
                         (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + offset(log(trial_length)),
                       hu ~ sel_disappearance + 
                         (1|ID) + (1|Family) + (1|Block) + (1|Chamber) + offset(log(trial_length))),
                    data = active0,
                    family = hurdle_gamma(),
                    warmup = 500, 
                    iter   = 2500, 
                    thin = 2,
                    control = list(adapt_delta = 0.95, max_treedepth = 15),
                    chains = 4, 
                    cores  = 4,
                    backend = "cmdstanr", 
                    threads = threading(2),
                    seed = 95)

pp_check(distance_brm0, ndraws = 100) + ggplot2::xlim(0,10000)

# Survival Models ---------------------------------------------------------


# read "Combined_Flight_final.csv" from "Final_Values_Folder. Be careful about
# Working directories and where new files are being written.
lifedata<- read.csv(file.choose(), stringsAsFactors = FALSE, header = T)

lifedata <- lifedata %>%
  mutate(ID = paste0(Family,ID)) %>%
  filter(Block != "A" & Block != "E") %>%
  mutate(oTreatment = ordered(Treatment, levels = c('A', 'B'))) %>%
  mutate(Family = factor(Family)) %>%
  mutate(ID = factor(ID)) %>%
  mutate(Block = factor(Block)) %>%
  select(Family,ID,Sex,oTreatment,Size,Censored,Lifespan,Activity,Block)


longlife <- survSplit(Surv(Lifespan,Censored) ~., 
                      lifedata,
                      cut=c(unique(lifedata$Lifespan)), 
                      episode ="timegroup")

longlife_flyers <- longlife[longlife$Activity == "Flight",]

system.time(gam.1 <- gamm4(Censored ~ oTreatment + s(Lifespan) + s(Lifespan, by = oTreatment), random = ~(1|ID) + (1|Family) + (1|Block),
                                data = longlife_flyers, family = binomial(link="cloglog")))


longlife_flyers <- gammit::predict_gamm(gam.1$gam, newdata = longlife_flyers, re.form = NULL, se = T, keep_prediction_data = T)
summary(gam.1$gam)

mort_1_fit <- ggplot(longlife_flyers, aes(x = Lifespan, y = prediction, col = oTreatment, group = oTreatment)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp")) +
  geom_ribbon(data = longlife_flyers, aes(ymin = prediction - 1.96*se, ymax = prediction + 1.96*se), alpha = 0.3, color = NA) +
  theme_classic() +
  labs(x = "Age (days)", y = "Log(mortality)") +
  xlim(0, 180)







