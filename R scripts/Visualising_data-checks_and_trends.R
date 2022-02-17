# Script for visualising flight data

# Author: Kynan Delaney
# Date: 14/02/2022

# Raw data plots of flight behaviours expressed by adult beetles at up to three
# time-points/ages (10-20 days old, 60-70 days old, 90+ days old).

# Libraries ---------------------------------------------------------------
rm(list = ls())

# for data cleaning and plotting
library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(gghalves)
library(ggbeeswarm)

# load in data ------------------------------------------------------------

# Be careful about Working directories and where new files are being written.
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
  mutate(total_peak_speed = as.numeric(total_peak_speed))

# Creating a subset to work with only beetles that flew. Removing missing values.
# Active0 includes all beetles, and makes a new factor for flight/no flight. Renaming
# some factors for better plotting.
# Active only includes assays in which the individual flew.
active0 <- flight %>%
  filter(Activity == "Flight") %>%
  filter(!is.na(Size)) %>%
  filter(age_at_flight > 0) %>%
  mutate(Treatment = ifelse(Treatment == "A", "Impoverished", "Enhanced")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Impoverished", "Enhanced"))) %>%
  mutate(Assay = recode_factor(Assay, "1" = "10 - 20 days", "2" = "60 - 70 days", "3" = "90 - 100 days")) %>%
  mutate(fly.binary = ifelse(total_distance > 0, 1,0))

active0[is.na(active0)] <- 0

active <- active0 %>%
  filter(total_distance > 0)


# Plotting and data checks ----------------------------------------

# A function to add sample size to plots. "y = 6.5" just hard-codes the height of 
# the label, not the value.
n_fun_y_6.5 <- function(x){
  return(data.frame(y = 6.5, label = length(x)))
}

# Visualising the effect of treatment on body size across blocks
ggplot(active0, aes(x = factor(Treatment), y = as.numeric(Size), colour = Block)) +
  geom_boxplot() +
  theme(legend.position = "none")+
  theme_classic() +
  stat_summary(fun.data = n_fun_y_6.5, geom = "text",
               aes(group=interaction(Block,Treatment)),
               hjust = 0.5, position = position_dodge(0.6)) 


# Visualisation of the effect of treatment across sex
ggplot(active0, aes(x = factor(Treatment), y = as.numeric(Size), colour = Sex)) +
  geom_boxplot() +
  theme(legend.position = "none")+
  theme_classic() +
  stat_summary(fun.data = n_fun_y_6.5, geom = "text",
               aes(group=interaction(Sex,Treatment)),
               hjust = 0.5, position = position_dodge(0.6))


# A function to add sample size to plots. "y = 50" just hard-codes the height of 
# the label, not the value.
n_fun_y_50 <- function(x){
  return(data.frame(y = c(50), label = length(x)))
}

# Checking the age distributions in each assay and block
ggplot(active0, aes(x = Assay, y = age_at_flight, col = Block)) +
  geom_violin() +
  stat_summary(fun.data = n_fun_y_50, geom = "text", aes(group=interaction(Block,Treatment)))



# Plotting data trends with 0s included -----------------------------------


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


ggplot(data = active0,aes(x = Assay, y = total_distance, group = ID)) +
  geom_point() +
  geom_line() 
