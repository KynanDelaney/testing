# Script for analysing Age-specific mortality in the beetles

# this script is mainly about GAMs and other spline-based approaches to estimating
# differences between groups, as opposed to estimating exact "rates" of ageing.

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
library(mgcv)

library(survival)
library(survminer)

# Survival Models ---------------------------------------------------------


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
# factor of Treatment for models. Remove unnecessary columns
lifedata <- lifedata %>%
  mutate(ID = paste0(Family,ID)) %>%
  filter(Block != "A") %>%
  filter(Lifespan >= 10) %>%
  mutate(oTreatment = ordered(Treatment, levels = c('A', 'B'))) %>%
  mutate(Family = factor(Family)) %>%
  mutate(ID = factor(ID)) %>%
  mutate(Block = factor(Block)) %>%
  select(Family,ID,Sex,oTreatment,Size,Censored,Lifespan,Activity,Block)
new$Size <- as.numeric(new$Size)

b<- lm(Size ~ group, data = new)
lifedata$group <- ifelse(lifedata$oTreatment == "A", "impoversidhed", "enhanced")
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

#write.csv(x=longlife_flyers, file="Lifespan.csv", row.names = F)

# look only at individuals who flew.
longlife_flyers <- longlife[longlife$Activity == "Flight",]

system.time(surv0 <- brm(Censored ~ s(Lifespan, k = 20) + (1|ID) + (1|Family) + (1|Block),
                         family = bernoulli(link = "cloglog"), 
                         data = longlife_flyers, 
                         warmup = 500, 
                         iter = 2500, 
                         thin = 2,
                         control = list(adapt_delta = 0.95),
                         chains = 3,
                         cores = 3,
                         seed = 95, 
                         backend = "cmdstanr", 
                         threads = threading(2))
)
saveRDS(surv0, file = "surv0.rds")


system.time(surv1 <- brm(Censored ~ oTreatment + s(Lifespan, k = 20) + (1|ID) + (1|Family) + (1|Block),
                         family = bernoulli(link = "cloglog"), 
                         data = longlife_flyers, 
                         warmup = 500, 
                         iter = 2500, 
                         thin = 2,
                         control = list(adapt_delta = 0.95),
                         chains = 3,
                         cores = 3,
                         seed = 95, 
                         backend = "cmdstanr", 
                         threads = threading(2))
)
saveRDS(surv1, file = "surv1.rds")


system.time(surv2 <- brm(Censored ~ oTreatment + s(Lifespan, k = 20) + s(Lifespan, by = oTreatment) + (1|ID) + (1|Family) + (1|Block),
                         family = bernoulli(link = "cloglog"), 
                         data = longlife_flyers, 
                         warmup = 500, 
                         iter = 2500, 
                         thin = 2,
                         control = list(adapt_delta = 0.95),
                         chains = 3,
                         cores = 3,
                         seed = 95, 
                         backend = "cmdstanr", 
                         threads = threading(2))
)


saveRDS(surv2, file = "surv2.rds")

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


system.time(mgcv_surv_0 <- gam(Censored ~ s(Lifespan) +
                                 s(ID, bs = "re") + s(Family, bs = "re") + s(Block, bs = "re"),
                               family = binomial(link = "cloglog"), data = longlife_flyers, method = 'REML',
                               select = T))
saveRDS(mgcv_surv_0, file = "mgcv_surv_0.rds")


system.time(mgcv_surv_1 <- gam(Censored ~ oTreatment + s(Lifespan) +
                                 s(ID, bs = "re") + s(Family, bs = "re") + s(Block, bs = "re"),
                               family = binomial(link = "cloglog"), data = longlife_flyers, method = 'REML',
                               select = T))
saveRDS(mgcv_surv_1, file = "mgcv_surv_1.rds")


system.time(mgcv_surv_2 <- gam(Censored ~ oTreatment + s(Lifespan) + s(Lifespan, by = oTreatment) + 
                                 s(ID, bs = "re") + s(Family, bs = "re") + s(Block, bs = "re"),
                               family = binomial(link = "cloglog"), data = longlife_flyers, method = 'REML',
                               select = T))
saveRDS(mgcv_surv_2, file = "mgcv_surv_2.rds")

anova(mgcv_surv_1,mgcv_surv_2, test = "LRT")



y <- ggpredict(
  mgcv_surv_1,
  terms = c("Lifespan", "oTreatment"),
  ci.lvl = 0.95,
  type = "fe",
  back.transform = F)
plot(y)
y$predicted <- -log(1-y$predicted)
y$conf.low <- -log(1-y$conf.low)
y$conf.high <- -log(1-y$conf.high)
y$group <- ifelse(y$group == "A", "impoverished", "Enhanced")

ggplot(y, aes(x = x, y = log(predicted), col = group)) + 
  geom_smooth(method = "loess", se = F) +
  geom_ribbon(data = y, aes(ymin = log(conf.low), ymax = log(conf.high), alpha = 0.9, fill = group)) +
  ylab("Log(mortality)") +
  xlab("Age (days)") + 
  ggtitle("Age-specific mortality of flying beetles") +
  theme_classic() +
  theme(plot.title = element_text(size = 24))

reduced <- lifedata[lifedata$Activity == "Flight",]
reduced$group <- ifelse(reduced$oTreatment == "A", "Impoverished", "Enhanced")
kap.1.fit <- survfit(Surv(Lifespan,Censored) ~ group, data = reduced)
kap_1_fit <- autoplot(kap.1.fit,
                      xlim = c(0, 175),
                      conf.int = T,
                      fun = "pct",
                      xlab = "Age (days)", ylab = "Survivorship")
n <- ggsurvplot(kap.1.fit, title = "Survivorship of flying beetles", font.title = 24)
n$plot + theme(legend.position = "none")
