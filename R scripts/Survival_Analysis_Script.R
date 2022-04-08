# Script for analysing Age-specific mortality in the beetles

# This script is mainly about GAMs and other spline-based approaches to estimating
# differences between groups, as opposed to estimating exact "rates" of ageing.


#  Clear environment ------------------------------------------------------


rm(list = ls())


# Libraries ---------------------------------------------------------------

# for data cleaning and plotting.
library(plyr)
library(dplyr)
library(lubridate)

# for analysis.
library(brms)
library(broom)
library(loo)  # model selection. 
library(mgcv)

library(mgcViz)
library(DHARMa)

# for visualisation.
library(ggeffects)
library(survival)
library(survminer)
library(patchwork)


# Read in and tidy raw data ---------------------------------------------------

# Be careful about Working directories and where new files are being written.
lifedata<- read.csv(file.choose(), stringsAsFactors = FALSE, header = T)

lifedata <- lifedata %>%
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

# Create more informative ID variable. Remove junk/pilot blocks. Create ordered 
# factor of Treatment for models. Remove unnecessary columns.
lifedata <- lifedata %>%
  mutate(ID = paste0(Family,ID)) %>%
  filter(Block != "A") %>%
  filter(Lifespan >= 10) %>%
  mutate(oTreatment = ordered(Treatment, levels = c('A', 'B'))) %>%
  mutate(group = ifelse(Treatment == "A", "impoverished", "enhanced")) %>%
  mutate(oActivity = ordered(Activity, levels = c('Flight', 'Waiting'))) %>%
  mutate(Family = factor(Family)) %>%
  mutate(ID = factor(ID)) %>%
  mutate(Block = factor(Block)) %>%
  mutate(Size = as.numeric(Size)) %>%
  select(Family,ID,Sex,oTreatment,oActivity,Size,Censored,Lifespan,Treatment,group,Activity,Block)

#write.csv(lifedata, "lifespan.csv", row.names = F)

# read in and work with cleaned data --------------------------------------

# Be careful about Working directories and where new files are being written.
#lifedata<- read.csv("CleanData/lifespan.csv", stringsAsFactors = FALSE, header = T)

# summary statistic for difference in size between treatment groups.
new <- subset(lifedata, select = c(Size, group))
new <- na.omit(new)
new %>% 
  mutate(group = factor(group)) %>%
  group_by(group) %>% 
  summarize(mean = mean(Size))

# Convert individual measuers to long-form "state per day/age" format. MUST use 
# Lifespan as response variable in models. I made the mistake with timegroup before.
longlife <- survSplit(Surv(Lifespan,Censored) ~., 
                      lifedata,
                      cut=c(unique(lifedata$Lifespan)), 
                      episode ="timegroup")

# look only at individuals who flew.
longlife_flyers <- longlife[longlife$Activity == "Flight",]


# GAMMs of age-specific mortality fit with mgcv -----------------------------------------

# GAMM fit with bam() to Massively speed up run-time. Only looking at adult 
# post-10 days/sexual maturity lifespan.
system.time(mgcv_surv_0 <- bam(Censored ~ s(Lifespan) +
                                 s(Family, bs = "re") + s(Block, bs = "re"),
                               family = binomial(link = "cloglog"), 
                               data = subset(longlife_flyers, Lifespan >= 10), 
                               select = T,
                               discrete = T,
                               nthreads=6, 
                               control =  gam.control(trace = TRUE)))

saveRDS(mgcv_surv_0, file = "mgcv_surv_0.rds")

system.time(mgcv_surv_1 <- bam(Censored ~ oTreatment + s(Lifespan) +
                                 s(Family, bs = "re") + s(Block, bs = "re"),
                               family = binomial(link = "cloglog"), 
                               data = subset(longlife_flyers, Lifespan >= 10), 
                               select = T,
                               discrete = T,
                               nthreads=6, 
                               control =  gam.control(trace = TRUE)))
saveRDS(mgcv_surv_1, file = "mgcv_surv_1.rds")



system.time(mgcv_surv_2k9 <- bam(Censored ~ oTreatment + s(Lifespan) + s(Lifespan, by = oTreatment) + 
                                 s(Family, bs = "re") + s(Block, bs = "re"),
                               family = binomial(link = "cloglog"), 
                               data = subset(longlife_flyers, Lifespan >= 10), 
                               select = T,
                               discrete = T,
                               nthreads=6, 
                               control =  gam.control(trace = TRUE)))
saveRDS(mgcv_surv_2, file = "mgcv_surv_2.rds")

gam.check(mgcv_surv_2,pch=19,cex=.3)

simulationOutput <- simulateResiduals(mgcv_surv_2, refit = T)
plot(simulationOutput,quantreg = T)
testDispersion(simulationOutput)


anova(mgcv_surv_0,mgcv_surv_1,mgcv_surv_2, test = "LRT")
anova(mgcv_surv_2)

# extract model estimates for plotting
y <- ggemmeans(
  mgcv_surv_2,
  terms = c("Lifespan", "oTreatment"),
  ci.lvl = 0.95,
  type = "fixed",
  rg.limit = 100000,
  back.transform = F)
plot(y)

# convert death rate to mortality. set confidence intervals to appropriate scale
y$predicted <- log(y$predicted)
y$conf.low <- log(y$conf.low)
y$conf.high <- log(y$conf.high)
y$group <- ifelse(y$group == "A", "Impoverished", "Enhanced")

mort_plot0 <- ggplot(y, aes(x = x, y = predicted, col = group)) +
  geom_line() +
  geom_ribbon(data = y, aes(ymin = conf.low, ymax = conf.high, alpha = 0.9, fill = group)) +
  ylab("Log(mortality)") +
  xlab("Age (days)") + 
  ggtitle("Age-specific mortality of flying beetles") +
  theme_classic() +
  xlim(10,180) +
  scale_x_continuous(breaks = c(10,50,100,150,180)) +
  theme(plot.title = element_text(size = 24)) 
  

reduced <- lifedata[lifedata$Lifespan >=10,]
reduced$group <- factor(reduced$group)
surv_obj <- survfit(Surv(Lifespan,Censored) ~ group, data = reduced)
surv_plot <- ggsurvplot(surv_obj, title = "Survivorship of flying beetles", conf.int = T, font.title = 24, xlim = c(10,180)) 
surv_plot <- surv_plot$plot + theme(legend.position = "none") + scale_x_continuous(breaks = c(10,50,100,150,180)) + xlab("Age (days)")  

combined <- surv_plot + mort_plot0 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")


# GAMMs of age-specific mortality fit with brms ---------------------------

system.time(brms_surv_2 <- brm(Censored ~ oTreatment + s(Lifespan, k = 4) + s(Lifespan, by = oTreatment, k = 4) + 
                      (1|ID) + (1|Family) + (1|Block),
                      family = bernoulli(link = "logit"),
                      data = subset(longlife_flyers, Lifespan >= 10),
                      warmup = 500,
                      iter = 2500,
                      thin = 2,
                      control = list(adapt_delta = 0.95),
                      chains = 4,
                      cores = 4,
                      seed = 123,
                      save_pars = save_pars(all = TRUE)))
                   

saveRDS(brms_surv_2, "Survival models/brms_surv_2.rds")
plot(brms_surv_2)

z <- ggemmeans(
  brms_surv_2,
  terms = c("Lifespan", "oTreatment"),
  ci.lvl = 0.95,
  type = "fixed",
  rg.limit = 100000,
  back.transform = F)
plot(z)

# convert death rate to mortality. set confidence intervals to appropriate scale
z$predicted <- log(z$predicted)
z$conf.low <- log(z$conf.low)
z$conf.high <- log(z$conf.high)
z$group <- ifelse(z$group == "A", "Impoverished", "Enhanced")

mort_plot1 <- ggplot(z, aes(x = x, y = predicted, col = group)) +
  geom_line() +
  geom_ribbon(data = z, aes(ymin = conf.low, ymax = conf.high, alpha = 0.9, fill = group)) +
  ylab("Log(mortality)") +
  xlab("Age (days)") + 
  ggtitle("Age-specific mortality of flying beetles") +
  theme_classic() +
  xlim(10,180) +
  scale_x_continuous(breaks = c(10,50,100,150,180)) +
  theme(plot.title = element_text(size = 24)) 



combined1 <- mort_plot0 + mort_plot1 & theme(legend.position = "bottom")
combined1 + plot_layout(guides = "collect")
