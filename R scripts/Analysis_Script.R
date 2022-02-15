
# Libraries ---------------------------------------------------------------
rm(list = ls())

# for data cleaning and plotting
library(tidyverse)
library(lubridate)
library(gridExtra)
library(gghalves)
library(ggbeeswarm)

# for analysis
library(brms)
library(cmdstanr) # backend for brms
library(loo)      # model selection 

library(survival)
library(survminer)





# load in data ------------------------------------------------------------

# Read-in main beetle measurements from "Final_Values_Folder".
beetles<- read.csv(file.choose(),
                   stringsAsFactors = FALSE, header = T)

# Make a more informative individual ID. Remove extraneous columns.
beetles <- beetles %>%
  mutate(ID = paste0(Family,ID)) %>%
  select(-Weight,-LastObserved,-Death.Removal,-Notes)

beetles$ID <- factor(beetles$ID)

# Converts a wide format of three trials to repeated-measures long format.
beetles <- reshape(beetles, 
              varying=list(Trial= c("X1st.Trial","X2nd.Trial", "X3rd.Trial"), 
                           Chamber= c("X1st.Chamber","X2nd.Chamber","X3rd.Chamber"), 
                           StartTime= c("X1st.StartTime","X2nd.StartTime","X3rd.StartTime"), 
                           Session= c("X1st.Session","X2nd.Session","X3rd.Session"), 
                           Comments= c("X1st.Comments","X2nd.Comments","X3rd.Comments") ), 
              v.names=c("Trial", "Chamber", "StartTime", "Session","Comments"), 
              direction="long",  
              times=1:3,
              timevar="Assay")

# removing extra columns. Labelling Assays more informatively. "Waiting" beetles
# can easily be filtered out - "Assay != 0,". Distinct() and filter() remove 
# duplications or empty rows created by reshape() above. age_at_flight is calculated 
# last, for each trial.
beetles <- beetles %>%
  select(-Comments, -id)%>%
  mutate(Assay = ifelse(Activity == "Waiting", 0, Assay)) %>%
  distinct(ID,Assay, .keep_all = T) %>%
  filter(!is.na(Session)) %>%
  mutate(age_at_flight = as.numeric(difftime(dmy(Trial), dmy(Eclosion))))

# Read-in un-paired flight  measurements from "Final_Values_Folder".
flight<- read.csv(file.choose(),
                  stringsAsFactors = FALSE, header = T)


# combining beetle id and measurements with trial, session and chamber info.
flight <- list(beetles,flight) %>% reduce(full_join, by = c("Session", "Chamber" = "Channel"))

# family shouldn't be missing from any rows.
flight<- flight %>%
  filter(!is.na(Family))

# outputting all this cleaned/paired data to a useable csv.
#write.csv(flight, file = "Combined_Flight_final.csv")

# read "Combined_Flight_final.csv" from "Final_Values_Folder. Be careful about
# Working directories and where new files are being written.
flight<- read.csv(file.choose(), stringsAsFactors = FALSE, header = T)


flight$Sex <- ifelse(flight$Sex == "", "M", flight$Sex)
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
  mutate(total_peak_speed = as.numeric(total_peak_speed))

# Removing data that either came from the pilot, or failed week.
flight <- flight %>%
  filter(Block != "A" & Block != "E")

# Creating a subset to work with only beetles that flew. Removing missing values.
active0 <- flight %>%
  filter(Activity == "Flight") %>%
  filter(!is.na(Size)) %>%
  filter(age_at_flight > 0) %>%
  mutate(fly.binary = ifelse(total_distance > 0, 1,0))

active0[is.na(active0)] <- 0

active <- active0 %>%
  filter(total_distance > 0)




# Plotting and data checks ----------------------------------------


# A function to add sample size to plots. "y = 6.5" just hard-codes the height of 
# the label, not the value.
n_fun <- function(x){
  return(data.frame(y = 6.5, label = length(x)))
}

# Visualising the effect of treatment on body size across blocks
ggplot(active0, aes(x = factor(Treatment), y = as.numeric(Size), colour = Block)) +
  geom_boxplot() +
  theme(legend.position = "none")+
  theme_classic() +
  stat_summary(fun.data = n_fun, geom = "text",
               aes(group=interaction(Block,Treatment)),
               hjust = 0.5, position = position_dodge(0.6)) 
  

# Visualisation of the effect of treatment across sex
ggplot(active0, aes(x = factor(Treatment), y = as.numeric(Size), colour = Sex)) +
  geom_boxplot() +
  theme(legend.position = "none")+
  theme_classic() +
  stat_summary(fun.data = n_fun, geom = "text",
               aes(group=interaction(Sex,Treatment)),
               hjust = 0.5, position = position_dodge(0.6))


# visualisation of the effect of treatment across sex
ggplot(active0, aes(x = factor(Treatment), y = as.numeric(Size))) +
  geom_violin() +
  theme(legend.position = "none")+
  theme_classic() +
  stat_summary(fun.data = n_fun, geom = "text",
               aes(group=interaction(Treatment)),
               hjust = 0.5, position = position_dodge(0.6))


n_fun_y_50 <- function(x){
  return(data.frame(y = c(50), label = length(x)))
}

# Checking the age distributions in each assay
ggplot(active0, aes(x = Assay, y = age_at_flight)) +
  geom_violin() +
  stat_summary(fun.data = n_fun_y_50, geom = "text")



# Plotting data trends with 0s included -----------------------------------
library(gridExtra)
library(ggplot2)
library(patchwork)

# changing some labels
active0 <- active0 %>%
  mutate(Treatment = ifelse(Treatment == "A", "Impoverished", "Enhanced")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Impoverished", "Enhanced"))) %>%
  mutate(Assay = recode_factor(Assay, "1" = "10 - 20 days", "2" = "60 - 70 days", "3" = "90 - 100 days"))
active <- active %>%
  mutate(Treatment = ifelse(Treatment == "A", "Impoverished", "Enhanced")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Impoverished", "Enhanced"))) %>%
  mutate(Assay = recode_factor(Assay, "1" = "10 - 20 days", "2" = "60 - 70 days", "3" = "90 - 100 days"))



n_fun_y_17.5 <- function(x){
  return(data.frame(y = 17.5, label = length(x)))
}

# plot raw distance, in kilometres, uncorrected for "trial_length"
distance_inc_0 <- ggplot(active0, aes(x = Assay, y = total_distance/1000, col = Treatment)) +
  geom_half_boxplot(outlier.color = NA) + geom_half_violin(side = "r") + 
  stat_summary(fun.data = n_fun_y_17.5, geom = "text",
               aes(group=interaction(Treatment)),
               hjust = 0.5, position = position_dodge(0.6)) + 
  labs(title = "Including non-flyers") + ylab("Kilometres flown") + theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

distance_ex_0 <- ggplot(active, aes(x = Assay, y = total_distance/1000, col = Treatment)) +
  geom_half_boxplot(outlier.color = NA) + geom_half_violin(side = "r") +
  stat_summary(fun.data = n_fun_y_17.5, geom = "text",
               aes(group=interaction(Treatment)),
               hjust = 0.5, position = position_dodge(0.6)) +
  labs(title = "Excluding non-flyers") + theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



n_fun_y_1.5 <- function(x){
  return(data.frame(y = 1.5, label = length(x)))
}

# Average speed in metres per second, calculated by distance/duration of actual flight
avg0 <- ggplot(active0, aes(x = Assay, y = total_avg_speed, col = Treatment)) +
  geom_half_boxplot(outlier.color = NA) + geom_half_violin(side = "r") +
  stat_summary(fun.data = n_fun_y_1.5, geom = "text",
               aes(group=interaction(Treatment)),
               hjust = 0.5, position = position_dodge(0.6)) +
  ylab("Average speed (m/s)")  + ylim(0,1.5) + theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

avg1 <- ggplot(active, aes(x = Assay, y = total_avg_speed, col = Treatment)) +
  geom_half_boxplot(outlier.color = NA) + geom_half_violin(side = "r") +
  stat_summary(fun.data = n_fun_y_1.5, geom = "text",
               aes(group=interaction(Treatment)),
               hjust = 0.5, position = position_dodge(0.6)) + ylim(0,1.5) + theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




n_fun_y_3.1 <- function(x){
  return(data.frame(y = 3.1, label = length(x)))
}

# peak speed recorded during assay. uncorrected for abbreviated trials. Capped at
# 3 for fast flyers, 0.1 for slow flyers and 0 for non-flyers.
peak0 <- ggplot(active0, aes(x = Assay, y = total_peak_speed, col = Treatment)) +
  geom_half_boxplot(outlier.color = NA) + geom_half_violin(side = "r") +
  stat_summary(fun.data = n_fun_y_3.1, geom = "text",
               aes(group=interaction(Treatment)),
               hjust = 0.5, position = position_dodge(0.6)) + ylim(0,3.2) +
  ylab("Peak speed (m/s)") + theme_classic()

peak1 <- ggplot(active, aes(x = Assay, y = total_peak_speed, col = Treatment)) +
  geom_half_boxplot(outlier.color = NA) + geom_half_violin(side = "r") +
  stat_summary(fun.data = n_fun_y_3.1, geom = "text",
               aes(group=interaction(Treatment)),
               hjust = 0.5, position = position_dodge(0.6)) + ylim(0,3.2) + theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


active0 <- active0 %>%
  mutate(fly.binary = ifelse(total_distance > 0, 1,0))

fly_prop <- ggplot(active0, aes(x = Assay, y = fly.binary, col = Treatment)) +
  stat_summary(fun=mean, geom="point", position = position_dodge(0.6)) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.2, position = position_dodge(0.6)) +
  ylim(0,1) + ylab("Proportion of beetles flying") + labs(title = "Proportion of beetles flying") +
  theme_classic() + theme(legend.title = element_blank())

ggplot(active0, aes(x = Assay, y = fly.binary, col = Treatment)) +
  geom_smooth() +
  ylim(0,1) + ylab("Proportion of beetles flying") + labs(title = "Proportion of beetles flying") +
  theme_classic() + theme(legend.title = element_blank())

design <- "
  123
  145
  167
"
combined <- fly_prop + distance_inc_0 + distance_ex_0 + avg0 + avg1 + peak0 + peak1 + plot_layout(design=design) & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")



# Converts a wide format of three trials to repeated-measures long format.
quali <- gather(active0, key = "flight.length", value = "flight.count",
                marathon_flights, long_flights, medium_flights, short_flights, skip_flights)

quali <- quali %>%
  filter(flight.count > 0) %>%
  mutate(flight.length = factor(flight.length,
                                levels = c("skip_flights",
                                           "short_flights", 
                                           "medium_flights", 
                                           "long_flights", 
                                           "marathon_flights"))) %>%
  filter(flight.length != "skip_flights")  

n_fun_y_50 <- function(x){
  return(data.frame(y = 50, label = length(x)))
}
ggplot(quali, aes(x = Assay, y = flight.count, col = Treatment)) +
  geom_boxplot() + 
  stat_summary(fun.data = n_fun_y_50, geom = "text",
               aes(group=interaction(Treatment, flight.length)),
               hjust = 0.5, position = position_dodge(0.6)) +
  facet_wrap(~flight.length)

# reworking data-frame to look at qualitative changes in flight.

# Converts a wide format of three trials to repeated-measures long format.
quali <- gather(active, key = "flight.length", value = "flight.count",
                marathon_flights, long_flights, medium_flights, short_flights, skip_flights)

quali <- quali %>%
  filter(flight.count > 0) %>%
  mutate(flight.length = factor(flight.length,
                                levels = c("skip_flights",
                                           "short_flights", 
                                           "medium_flights", 
                                           "long_flights", 
                                           "marathon_flights"))) %>%
  filter(flight.length != "skip_flights") %>%
  filter(flight.length != "marathon_flights") 

ggplot(quali, aes(x = flight.length, y = flight.count, col = Treatment)) +
  geom_boxplot() +
  facet_wrap(~Assay)




# Early models ------------------------------------------------------------

active0$total_distance <- round(active0$total_distance)
active0$Lifespan <- scale(active0$Lifespan)
null_distance_model <- brm(total_distance ~ 1 + Lifespan + (1|ID) + (1|Family) + (1|Block),
                  data = active0,
                  family = zero_inflated_negbinomial(),
                  warmup = 500, 
                  iter   = 2500, 
                  thin = 2,
                  prior = c(prior("normal(0, 5)", class = "b"),
                            prior("normal(0, 5)", class = "Intercept") +
                            prior("beta(2,2)", class = "zi")),
                  control = list(adapt_delta = 0.95, max_treedepth = 15),
                  chains = 4, 
                  cores  = 4,
                  backend = "cmdstanr", 
                  threads = threading(2),
                  seed = 95)

pp_check(null_distance_model, ndraws = 100) + ggplot2::xlim(0,10000)
saveRDS(null_distance_model, file.choose())

med_distance_model <- brm(total_distance ~ Assay+Treatment + Lifespan + (1|ID) + (1|Family) + (1|Block),
                          data = active0,
                          family = zero_inflated_negbinomial(),
                          warmup = 500, 
                          iter   = 2500, 
                          thin = 2,
                          control = list(adapt_delta = 0.95, max_treedepth = 15),
                          chains = 4, 
                          cores  = 4,
                          backend = "cmdstanr", 
                          threads = threading(2),
                          seed = 95)
saveRDS(med_distance_model, file.choose())

full_distance_model <- brm(total_distance ~ Assay*Treatment + Lifespan + (1|ID) + (1|Family) + (1|Block),
                              data = active0,
                              family = zero_inflated_negbinomial(),
                              warmup = 500, 
                              iter   = 2500, 
                              thin = 2,
                              control = list(adapt_delta = 0.95, max_treedepth = 15),
                              chains = 4, 
                              cores  = 4,
                              backend = "cmdstanr", 
                              threads = threading(2),
                              seed = 95)

saveRDS(full_distance_model, file.choose())
pp_check(full_distance_model)


null_prop_model <- brm(fly.binary ~ 1 + Lifespan + 
                         (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
                           data = active0,
                           family = bernoulli(),
                           warmup = 500, 
                           iter   = 2500, 
                           thin = 2,
                           control = list(adapt_delta = 0.95, max_treedepth = 15),
                           chains = 4, 
                           cores  = 4,
                           backend = "cmdstanr", 
                           save_pars = save_pars(all = T),
                           threads = threading(2),
                           seed = 95)


med_prop_model <- brm(fly.binary ~ Treatment + Assay + Lifespan + 
                        (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
                       data = active0,
                       family = bernoulli(),
                       warmup = 500, 
                       iter   = 2500, 
                       thin = 2,
                       control = list(adapt_delta = 0.95, max_treedepth = 15),
                       chains = 4, 
                       cores  = 4,
                       backend = "cmdstanr",
                       save_pars = save_pars(all = T),
                       threads = threading(2),
                       seed = 95)


full_prop_model <- brm(fly.binary ~ Treatment*Assay + Lifespan + 
                         (1|ID) + (1|Family) + (1|Block) + (1|Chamber),
                      data = active0,
                      family = bernoulli(),
                      warmup = 500, 
                      iter   = 2500, 
                      thin = 2,
                      control = list(adapt_delta = 0.95, max_treedepth = 15),
                      chains = 4, 
                      cores  = 4,
                      backend = "cmdstanr",
                      save_pars = save_pars(all = T),
                      threads = threading(2),
                      seed = 95)

pp_check(null_prop_model, ndraws = 100)
pp_check(med_prop_model, ndraws = 100)
pp_check(full_prop_model, ndraws = 100)


loo_null_prop <- loo::loo(null_prop_model)
loo_med_prop <- loo::loo(med_prop_model)
loo_full_prop <- loo::loo(full_prop_model)

n <- loo_compare(x = list(loo_null_prop,loo_med_prop, loo_full_prop))

print(n, simplify = F)

loo_med_prop_model_lin <- loo::loo(med_prop_model_lin)
loo_med_prop_model_quad <- loo::loo(med_prop_model_quad)

