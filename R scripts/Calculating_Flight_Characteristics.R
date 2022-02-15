# Template Script for processing flightmill data ---------------

# Author: Kynan Delaney
# Date: 10/11/2021

# A script to streamline calculating flight characteristics of beetles. I tried 
# to implement thresholds and data cleaning consistent with: 
# https://www.jove.com/t/53377/a-simple-flight-mill-for-the-study-of-tethered-flight-in-insects





rm(list = ls())

# Libraries and options ---------------------------------------------------

library(dplyr)
library(readr)
library(plyr)
library(purrr)
library(lubridate)
library(data.table)
library(chron)

options("digits.secs"=6)

# Loading in extracted CSVs -----------------------------------------------
# Basic data cleaning and eliminating duplicated consecutive values in the dataframe --------------

#setwd("~/Edinburgh/My Data/Flight Data/New Extracted Data")

rm(list = ls())

# Reading in and looking at data.
df <- read.csv(file.choose())

# automatically updates session number from each new file read-in.
session <- as.numeric(df$session[1])
str(df)

# Data is inserted into the database in batches. Th is arranges data by "Primary Key"
# of row number (id). Prevents errors in time-elapsed measures.
df <- df %>%
  arrange(id)

ddply(df, "session", summarize,
      start = min((ymd_hms(recorded_time))),
      end = max((ymd_hms(recorded_time))),
      length = (difftime(end,start, units = c("secs"))))

# Remove observations from mill set-up time. Based on written recorded start time.
df <- df %>% 
  filter(recorded_time >= "2021-09-01 12:25:00")
#session <- "6b"

# Changing time from generic time-stamp to time elapsed from session start. Removing
# unnecessary columns.
# Converting channel values from raw voltage steps to binary event/non-event.
df <- df %>%
  mutate(time = difftime(ymd_hms(recorded_time),first(ymd_hms(df$recorded_time)), units = c("secs"))) %>%
  select(id, session, time, channel_0:channel_6) %>%
  mutate_each(funs(ifelse(.<300, 0,1)), -id, - session, - time)

df %>%
  summarize(start = first(time),
            end = max(time),
            length = end-start)


# Discounting observations outside of the 4-hour window.
df <- df %>%
  filter(time <= 14400)



trial_length <- as.numeric(max(df$time) - first(df$time))

# A legacy of the old script. Discounts the first 2400 observations. Imprecise!
#flight <- tail(flight, -2400)

# A data.table command. allows each duplicated consecutive values be removed
# WITHOUT having to split up dataframe to preserve timestamp information.
# An unlooped [:(] method for removing duplicated "1" values. Change "df" 
# placemarker for real dataframe.
setDT(df)

df[, channel_0 := ifelse(duplicated(rleid(channel_0)), 0, channel_0), by = "session"]
df[, channel_1 := ifelse(duplicated(rleid(channel_1)), 0, channel_1), by = "session"]
df[, channel_2 := ifelse(duplicated(rleid(channel_2)), 0, channel_2), by = "session"]
df[, channel_3 := ifelse(duplicated(rleid(channel_3)), 0, channel_3), by = "session"]
df[, channel_4 := ifelse(duplicated(rleid(channel_4)), 0, channel_4), by = "session"]
df[, channel_5 := ifelse(duplicated(rleid(channel_5)), 0, channel_5), by = "session"]
df[, channel_6 := ifelse(duplicated(rleid(channel_6)), 0, channel_6), by = "session"]



# Calculating speeds and correcting for error/noise -----------------------
chambers <- c("channel_0",
              "channel_1",
              "channel_2",
              "channel_3",
              "channel_4",
              "channel_5",
              "channel_6")

# A loop to calculate the "time since last event" within each chamber.
for (i in chambers){
  df[,paste0(i, "_intervals") :=c(NA,diff(time)) , by = i]
}

# Unlooped method example.
#df[,channel_0_intervals :=c(NA,diff(time)) , by = channel_0]


# Converting data from 'data.table' syntax to standard dataframe.
df <- as.data.frame(df)

# Turning time since last rotation into a measure of speed, in metres per second,
# while removing excess time since last "0" values.
df <- df %>%
  mutate(channel_0_intervals = ifelse(channel_0 == 0, 0, 0.88/channel_0_intervals)) %>%
  mutate(channel_1_intervals = ifelse(channel_1 == 0, 0, 0.88/channel_1_intervals)) %>%
  mutate(channel_2_intervals = ifelse(channel_2 == 0, 0, 0.88/channel_2_intervals)) %>%
  mutate(channel_3_intervals = ifelse(channel_3 == 0, 0, 0.88/channel_3_intervals)) %>%
  mutate(channel_4_intervals = ifelse(channel_4 == 0, 0, 0.88/channel_4_intervals)) %>%
  mutate(channel_5_intervals = ifelse(channel_5 == 0, 0, 0.88/channel_5_intervals)) %>%
  mutate(channel_6_intervals = ifelse(channel_6 == 0, 0, 0.88/channel_6_intervals))

# Removing rows of all zeros that bloat out the dataframe. Doesn't work with column
# names for some reason. Makes future calculating and plotting quicker.
df <- df %>% 
  filter(rowSums(.[4:17]) > 0)


# First diagnostic plots --------------------------------------------------
# Setting plot parameters to compare thresholded and unthresholded data.
#par(mfrow = c(1,1))

# "Diagnostic" plot of flightmill runs. Layering individuals over each other to 
# help set reasonable threshold values below.
plot(channel_0_intervals ~ time, data = df, ylim=c(0,6))
points(df$channel_1_intervals, col = "red")
points(df$channel_2_intervals, col = "green")
points(df$channel_3_intervals, col = "aquamarine4")
points(df$channel_4_intervals, col = "blue")
points(df$channel_5_intervals, col = "purple")
points(df$channel_6_intervals, col = "orange")
abline(h=3)


# Thresholds continued ----------------------------------------------------

# Setting upper threshold to account for "noise" from sensor readings/obstructions.
# This happens separate from the lower threshold because high thresholds are purely 
# noise, whereas low values are more likely true rotations that may be the aftermath
# of beetle activity, or poor flight starts. 3m/s equates to ~3.4 revolutions 
# every second.
df <- df %>%
  mutate(channel_0_intervals = ifelse(channel_0_intervals > 3, 0, channel_0_intervals)) %>%
  mutate(channel_1_intervals = ifelse(channel_1_intervals > 3, 0, channel_1_intervals)) %>%
  mutate(channel_2_intervals = ifelse(channel_2_intervals > 3, 0, channel_2_intervals)) %>%
  mutate(channel_3_intervals = ifelse(channel_3_intervals > 3, 0, channel_3_intervals)) %>%
  mutate(channel_4_intervals = ifelse(channel_4_intervals > 3, 0, channel_4_intervals)) %>%
  mutate(channel_5_intervals = ifelse(channel_5_intervals > 3, 0, channel_5_intervals)) %>%
  mutate(channel_6_intervals = ifelse(channel_6_intervals > 3, 0, channel_6_intervals))

# Potentially unnecessary wrangling of noise correction ---------------------------------------

# Removing events (1) associated with those high values. These values are being
# removed to recalculate speeds in the absence of sensor-noise. This may be necessary
# if there is a lot of noise, but may also be unnecessary wrangling of the data.
#
# An alternative to this step is simply bound the speed values in upper and lower 
# thresholds and treat noise the same as drift/slow spins. Measures of distance 
# will be the same, but some calculations of average speed will be affected.
df <- df %>%
  mutate(channel_0 = ifelse(channel_0_intervals == 0, 0, channel_0)) %>%
  mutate(channel_1 = ifelse(channel_1_intervals == 0, 0, channel_1)) %>%
  mutate(channel_2 = ifelse(channel_2_intervals == 0, 0, channel_2)) %>%
  mutate(channel_3 = ifelse(channel_3_intervals == 0, 0, channel_3)) %>%
  mutate(channel_4 = ifelse(channel_4_intervals == 0, 0, channel_4)) %>%
  mutate(channel_5 = ifelse(channel_5_intervals == 0, 0, channel_5)) %>%
  mutate(channel_6 = ifelse(channel_6_intervals == 0, 0, channel_6)) 

# A loop to calculate the "time since last event" within each chamber. Repeated
# after raw data corrected for sensor-noise. setDT and as.data.frame 
setDT(df)
for (i in chambers){
  df[,paste0(i, "_intervals") :=c(NA,diff(time)) , by = i]
}
df <- as.data.frame(df)

# Turning time since last rotation into a measure of speed, while removing excess
# time since last "0" values. Repeated after raw data corrected for sensor-noise.
df <- df %>%
  mutate(channel_0_intervals = ifelse(channel_0 == 0, 0, 0.88/channel_0_intervals)) %>%
  mutate(channel_1_intervals = ifelse(channel_1 == 0, 0, 0.88/channel_1_intervals)) %>%
  mutate(channel_2_intervals = ifelse(channel_2 == 0, 0, 0.88/channel_2_intervals)) %>%
  mutate(channel_3_intervals = ifelse(channel_3 == 0, 0, 0.88/channel_3_intervals)) %>%
  mutate(channel_4_intervals = ifelse(channel_4 == 0, 0, 0.88/channel_4_intervals)) %>%
  mutate(channel_5_intervals = ifelse(channel_5 == 0, 0, 0.88/channel_5_intervals)) %>%
  mutate(channel_6_intervals = ifelse(channel_6 == 0, 0, 0.88/channel_6_intervals))

# Implementing lower threshold to account for drifting after flying events, or 
# poor flight starts. 0.1m/s equates to 1 revolution every ~9 seconds. Very slow!
df <- df %>%
  mutate(across(10:16, ~replace(., . <0.1, 0)))

# Removing rows of all zeros that bloat out the dataframe. Repeated to account for
# thresholding eliminations.
df <- df %>% 
  filter(rowSums(.[4:17]) > 0)
  

# "Diagnostic" plot of flightmill runs. Layering individuals over each other to 
# help set reasonable threshold values below. Par(mfrow = c(2,1)) above lets us
# look at both thresholded and unthresholded data.
#plot(channel_0_intervals ~ time, data = df, ylim=c(0,6))
#points(df$channel_1_intervals, col = "red")
#points(df$channel_2_intervals, col = "green")
#points(df$channel_3_intervals, col = "aquamarine4")
#points(df$channel_4_intervals, col = "blue")
#points(df$channel_5_intervals, col = "purple")
#points(df$channel_6_intervals, col = "orange")
#abline(h=3)



# Splitting up dataframe into separate channels and timestamps --------

# A massive vector of letters A >>> ZZZ for labelling individual flight events.
let_vec <- c(LETTERS,
             do.call("paste0",CJ(LETTERS,LETTERS)),
             do.call("paste0",CJ(LETTERS,LETTERS,LETTERS)))

# Splitting up channels and removing all 0 values. Length of each df * 0.88 is 
# equal to total distance travelled.
df$time <- as.numeric(df$time)
df0 <- df %>% 
  filter(channel_0_intervals != 0) %>% select(time, channel_0_intervals)

# Calculating the time difference between each element in the list. Large gaps
# differentiate one flying event from another.
setDT(df0)
df0[,interval_difference :=c(0,diff(time))]
df0 <- as.data.frame(df0)

# Loop for labelling each flight event separated by >14 seconds. First cluster of 
# activity is labelled "A", increasing to B,C,D,...etc for each distinct cluster 
# of flying activity, or lone/paired data points.
dat_vec <- df0$interval_difference
group_vec <- matrix(nrow = length(dat_vec), ncol = 1)
x <- 1

for (i in 1:length(dat_vec)){
  if(dat_vec[i] > 14){
    x <- x+1
    group_vec[i,1] <- let_vec[x]
  } else {
    group_vec[i,1] <- let_vec[x]
  }
}
group_vec[length(dat_vec)] <- let_vec[x]
df0$grouping <- group_vec


# Same process as above and repeated for each subsequent channel.
df1 <- df %>%
  filter(channel_1_intervals != 0) %>% select(time, channel_1_intervals)
setDT(df1)
df1[,interval_difference :=c(0,diff(time))]
df1 <- as.data.frame(df1)

dat_vec <- df1$interval_difference
group_vec <- matrix(nrow = length(dat_vec), ncol = 1)
x <- 1

for (i in 1:length(dat_vec)){
  if(dat_vec[i] > 14){
    x <- x+1
    group_vec[i,1] <- let_vec[x]
  } else {
    group_vec[i,1] <- let_vec[x]
  }
}
group_vec[length(dat_vec)] <- let_vec[x]
df1$grouping <- group_vec


df2 <- df %>%
  filter(channel_2_intervals != 0) %>% select(time, channel_2_intervals)
setDT(df2)
df2[,interval_difference :=c(0,diff(time))]
df2 <- as.data.frame(df2)

dat_vec <- df2$interval_difference
group_vec <- matrix(nrow = length(dat_vec), ncol = 1)
x <- 1

for (i in 1:length(dat_vec)){
  if(dat_vec[i] > 14){
    x <- x+1
    group_vec[i,1] <- let_vec[x]
  } else {
    group_vec[i,1] <- let_vec[x]
  }
}
group_vec[length(dat_vec)] <- let_vec[x]
df2$grouping <- group_vec


# Dear god how could this be looped?
df3 <- df %>%
  filter(channel_3_intervals != 0) %>% select(time, channel_3_intervals)
setDT(df3)
df3[,interval_difference :=c(0,diff(time))]
df3 <- as.data.frame(df3)

dat_vec <- df3$interval_difference
group_vec <- matrix(nrow = length(dat_vec), ncol = 1)
x <- 1

for (i in 1:length(dat_vec)){
  if(dat_vec[i] > 14){
    x <- x+1
    group_vec[i,1] <- let_vec[x]
  } else {
    group_vec[i,1] <- let_vec[x]
  }
}
group_vec[length(dat_vec)] <- let_vec[x]
df3$grouping <- group_vec


df4 <- df %>%
  filter(channel_4_intervals != 0) %>% select(time, channel_4_intervals)
setDT(df4)
df4[,interval_difference :=c(0,diff(time))]
df4 <- as.data.frame(df4)

dat_vec <- df4$interval_difference
group_vec <- matrix(nrow = length(dat_vec), ncol = 1)
x <- 1

for (i in 1:length(dat_vec)){
  if(dat_vec[i] > 14){
    x <- x+1
    group_vec[i,1] <- let_vec[x]
  } else {
    group_vec[i,1] <- let_vec[x]
  }
}
group_vec[length(dat_vec)] <- let_vec[x]
df4$grouping <- group_vec


df5 <- df %>%
  filter(channel_5_intervals != 0) %>% select(time, channel_5_intervals)
setDT(df5)
df5[,interval_difference :=c(0,diff(time))]
df5 <- as.data.frame(df5)

dat_vec <- df5$interval_difference
group_vec <- matrix(nrow = length(dat_vec), ncol = 1)
x <- 1

for (i in 1:length(dat_vec)){
  if(dat_vec[i] > 14){
    x <- x+1
    group_vec[i,1] <- let_vec[x]
  } else {
    group_vec[i,1] <- let_vec[x]
  }
}
group_vec[length(dat_vec)] <- let_vec[x]
df5$grouping <- group_vec


df6 <- df %>%
  filter(channel_6_intervals != 0) %>% select(time, channel_6_intervals)
setDT(df6)
df6[,interval_difference :=c(0,diff(time))]
df6 <- as.data.frame(df6)

dat_vec <- df6$interval_difference
group_vec <- matrix(nrow = length(dat_vec), ncol = 1)
x <- 1

for (i in 1:length(dat_vec)){
  if(dat_vec[i] > 14){
    x <- x+1
    group_vec[i,1] <- let_vec[x]
  } else {
    group_vec[i,1] <- let_vec[x]
  }
}
group_vec[length(dat_vec)] <- let_vec[x]
df6$grouping <- group_vec


# Flight description and summary tables -----------------------------------

# Summary table for channel 0. Gives a breakdown of each individual flight event.
# defined by periods of activity more than 14 seconds apart. Dist = distance 
# travelled during that flight. Number of elements in vector (i.e. rotations) times
# the circumference of the mill. Duration is the time elapsed since start of each
# flight stint and the end, in seconds. Peak = fastest speed recorded within each 
# flight stint. General_avg = distance/duration. Speed_avg = the mean of recorded
# speed values that make up the elements in this column. These two values should
# strongly correlate, but are differently sensitive to threshold values. 
# filter() removes lone and paired recorded events that are some activity (false
# start, swinging or other action) that might not be considered a true signal of 
# flight.

# if there was no activity in a chamber, summary will fail, and an empty chamber list 
# will need to be made. This empty dataframe prevents errors in that case.
channel_0_summary <- data.frame(dist=c(0),duration=c(0),peak=c(0),general_avg=c(0),speed_avg=c(0))

channel_0_summary <- df0 %>%
  dplyr::group_by(grouping) %>%
  dplyr::filter(length(grouping) > 2) %>%
  dplyr::summarise(dist = length(channel_0_intervals)*0.88, 
            duration = last(time)-first(time), 
            peak = max(channel_0_intervals), 
            general_avg = dist/duration, 
            speed_avg = mean(channel_0_intervals))

c0_chann <- 0
c0_total_distance <- sum(channel_0_summary$dist)
c0_total_duration <- sum(channel_0_summary$duration)
c0_total_average_speed <- c0_total_distance/c0_total_duration
c0_total_peak_speed <- max(channel_0_summary$peak)
c0_marathon_flights <- nrow(subset(channel_0_summary, duration >= 7200))
c0_long_flights <- nrow(subset(channel_0_summary, duration <= 7200 & duration >= 3600))
c0_medium_flights <- nrow(subset(channel_0_summary, duration <= 3600 & duration >= 600))
c0_short_flights <- nrow(subset(channel_0_summary, duration <= 600 & duration >= 60))
c0_skip_flights <- nrow(subset(channel_0_summary, duration <= 60))

# if there was no activity in a chamber, summary will fail, and an empty chamber list 
# will need to be made. This empty dataframe prevents errors in that case.
channel_1_summary <- data.frame(dist=c(0),duration=c(0),peak=c(0),general_avg=c(0),speed_avg=c(0))

channel_1_summary <- df1 %>%
  dplyr::group_by(grouping) %>%
  dplyr::filter(length(grouping) > 2) %>%
  dplyr::summarise(dist = length(channel_1_intervals)*0.88, 
            duration = last(time)-first(time), 
            peak = max(channel_1_intervals), 
            general_avg = dist/duration, 
            speed_avg = mean(channel_1_intervals))

c1_chann <- 1
c1_total_distance <- sum(channel_1_summary$dist)
c1_total_duration <- sum(channel_1_summary$duration)
c1_total_average_speed <- c1_total_distance/c1_total_duration
c1_total_peak_speed <- max(channel_1_summary$peak)
c1_marathon_flights <- nrow(subset(channel_1_summary, duration >= 7200))
c1_long_flights <- nrow(subset(channel_1_summary, duration <= 7200 & duration >= 3600))
c1_medium_flights <- nrow(subset(channel_1_summary, duration <= 3600 & duration >= 600))
c1_short_flights <- nrow(subset(channel_1_summary, duration <= 600 & duration >= 60))
c1_skip_flights <- nrow(subset(channel_1_summary, duration <= 60))

# if there was no activity in a chamber, summary will fail, and an empty c2 list 
# will need to be made. This empty dataframe prevents errors in that case.
channel_2_summary <- data.frame(dist=c(0),duration=c(0),peak=c(0),general_avg=c(0),speed_avg=c(0))

channel_2_summary <- df2 %>%
  dplyr::group_by(grouping) %>%
  dplyr::filter(length(grouping) > 2) %>%
  dplyr::summarise(dist = length(channel_2_intervals)*0.88, 
            duration = last(time)-first(time), 
            peak = max(channel_2_intervals), 
            general_avg = dist/duration, 
            speed_avg = mean(channel_2_intervals))


c2_chann <- 2
c2_total_distance <- sum(channel_2_summary$dist)
c2_total_duration <- sum(channel_2_summary$duration)
c2_total_average_speed <- c2_total_distance/c2_total_duration
c2_total_peak_speed <- max(channel_2_summary$peak)
c2_marathon_flights <- nrow(subset(channel_2_summary, duration >= 7200))
c2_long_flights <- nrow(subset(channel_2_summary, duration <= 7200 & duration >= 3600))
c2_medium_flights <- nrow(subset(channel_2_summary, duration <= 3600 & duration >= 600))
c2_short_flights <- nrow(subset(channel_2_summary, duration <= 600 & duration >= 60))
c2_skip_flights <- nrow(subset(channel_2_summary, duration <= 60))


# if there was no activity in a chamber, summary will fail, and an empty chamber list 
# will need to be made. This empty dataframe prevents errors in that case.
channel_3_summary <- data.frame(dist=c(0),duration=c(0),peak=c(0),general_avg=c(0),speed_avg=c(0))

channel_3_summary <- df3 %>%
  dplyr::group_by(grouping) %>%
  dplyr::filter(length(grouping) > 2) %>%
  dplyr::summarise(dist = length(channel_3_intervals)*0.88, 
            duration = last(time)-first(time), 
            peak = max(channel_3_intervals), 
            general_avg = dist/duration, 
            speed_avg = mean(channel_3_intervals))

c3_chann <- 3
c3_total_distance <- sum(channel_3_summary$dist)
c3_total_duration <- sum(channel_3_summary$duration)
c3_total_average_speed <- c3_total_distance/c3_total_duration
c3_total_peak_speed <- max(channel_3_summary$peak)
c3_marathon_flights <- nrow(subset(channel_3_summary, duration >= 7200))
c3_long_flights <- nrow(subset(channel_3_summary, duration <= 7200 & duration >= 3600))
c3_medium_flights <- nrow(subset(channel_3_summary, duration <= 3600 & duration >= 600))
c3_short_flights <- nrow(subset(channel_3_summary, duration <= 600 & duration >= 60))
c3_skip_flights <- nrow(subset(channel_3_summary, duration <= 60))


# if there was no activity in a chamber, summary will fail, and an empty chamber list 
# will need to be made. This empty dataframe prevents errors in that case.
channel_4_summary <- data.frame(dist=c(0),duration=c(0),peak=c(0),general_avg=c(0),speed_avg=c(0))

channel_4_summary <- df4 %>%
  dplyr::group_by(grouping) %>%
  dplyr::filter(length(grouping) > 2) %>%
  dplyr::summarise(dist = length(channel_4_intervals)*0.88, 
            duration = last(time)-first(time), 
            peak = max(channel_4_intervals), 
            general_avg = dist/duration, 
            speed_avg = mean(channel_4_intervals))

c4_chann <- 4
c4_total_distance <- sum(channel_4_summary$dist)
c4_total_duration <- sum(channel_4_summary$duration)
c4_total_average_speed <- c4_total_distance/c4_total_duration
c4_total_peak_speed <- max(channel_4_summary$peak)
c4_marathon_flights <- nrow(subset(channel_4_summary, duration >= 7200))
c4_long_flights <- nrow(subset(channel_4_summary, duration <= 7200 & duration >= 3600))
c4_medium_flights <- nrow(subset(channel_4_summary, duration <= 3600 & duration >= 600))
c4_short_flights <- nrow(subset(channel_4_summary, duration <= 600 & duration >= 60))
c4_skip_flights <- nrow(subset(channel_4_summary, duration <= 60))


# if there was no activity in a chamber, summary will fail, and an empty chamber list 
# will need to be made. This empty dataframe prevents errors in that case.
channel_5_summary <- data.frame(dist=c(0),duration=c(0),peak=c(0),general_avg=c(0),speed_avg=c(0))

channel_5_summary <- df5 %>%
  dplyr::group_by(grouping) %>%
  dplyr::filter(length(grouping) > 2) %>%
  dplyr::summarise(dist = length(channel_5_intervals)*0.88, 
            duration = last(time)-first(time), 
            peak = max(channel_5_intervals), 
            general_avg = dist/duration, 
            speed_avg = mean(channel_5_intervals))

c5_chann <- 5
c5_total_distance <- sum(channel_5_summary$dist)
c5_total_duration <- sum(channel_5_summary$duration)
c5_total_average_speed <- c5_total_distance/c5_total_duration
c5_total_peak_speed <- max(channel_5_summary$peak)
c5_marathon_flights <- nrow(subset(channel_5_summary, duration >= 7200))
c5_long_flights <- nrow(subset(channel_5_summary, duration <= 7200 & duration >= 3600))
c5_medium_flights <- nrow(subset(channel_5_summary, duration <= 3600 & duration >= 600))
c5_short_flights <- nrow(subset(channel_5_summary, duration <= 600 & duration >= 60))
c5_skip_flights <- nrow(subset(channel_5_summary, duration <= 60))


# if there was no activity in a chamber, summary will fail, and an empty chamber list 
# will need to be made. This empty dataframe prevents errors in that case.
channel_6_summary <- data.frame(dist=c(0),duration=c(0),peak=c(0),general_avg=c(0),speed_avg=c(0))

channel_6_summary <- df6 %>%
  dplyr::group_by(grouping) %>%
  dplyr::filter(length(grouping) > 2) %>%
  dplyr::summarise(dist = length(channel_6_intervals)*0.88, 
            duration = last(time)-first(time), 
            peak = max(channel_6_intervals), 
            general_avg = dist/duration, 
            speed_avg = mean(channel_6_intervals))

c6_chann <- 6
c6_total_distance <- sum(channel_6_summary$dist)
c6_total_duration <- sum(channel_6_summary$duration)
c6_total_average_speed <- c6_total_distance/c6_total_duration
c6_total_peak_speed <- max(channel_6_summary$peak)
c6_marathon_flights <- nrow(subset(channel_6_summary, duration >= 7200))
c6_long_flights <- nrow(subset(channel_6_summary, duration <= 7200 & duration >= 3600))
c6_medium_flights <- nrow(subset(channel_6_summary, duration <= 3600 & duration >= 600))
c6_short_flights <- nrow(subset(channel_6_summary, duration <= 600 & duration >= 60))
c6_skip_flights <- nrow(subset(channel_6_summary, duration <= 60))


# Combining all descriptions of flight into dataframe to append to main csv.
chann <- c(c0_chann,
           c1_chann,
           c2_chann,
           c3_chann,
           c4_chann,
           c5_chann,
           c6_chann)
total_distance <- c(c0_total_distance,
                    c1_total_distance,
                    c2_total_distance,
                    c3_total_distance,
                    c4_total_distance,
                    c5_total_distance,
                    c6_total_distance)
total_duration <- c(c0_total_duration,
                    c1_total_duration,
                    c2_total_duration,
                    c3_total_duration,
                    c4_total_duration,
                    c5_total_duration,
                    c6_total_duration)
total_avg_speed <- c(c0_total_average_speed,
                     c1_total_average_speed,
                     c2_total_average_speed,
                     c3_total_average_speed,
                     c4_total_average_speed,
                     c5_total_average_speed,
                     c6_total_average_speed)
total_peak_speed <- c(c0_total_peak_speed,
                      c1_total_peak_speed,
                      c2_total_peak_speed,
                      c3_total_peak_speed,
                      c4_total_peak_speed,
                      c5_total_peak_speed,
                      c6_total_peak_speed)
marathon_flights <- c(c0_marathon_flights,
                      c1_marathon_flights,
                      c2_marathon_flights,
                      c3_marathon_flights,
                      c4_marathon_flights,
                      c5_marathon_flights,
                      c6_marathon_flights)
long_flights <- c(c0_long_flights,
                  c1_long_flights,
                  c2_long_flights,
                  c3_long_flights,
                  c4_long_flights,
                  c5_long_flights,
                  c6_long_flights)
medium_flights <- c(c0_medium_flights,
                    c1_medium_flights,
                    c2_medium_flights,
                    c3_medium_flights,
                    c4_medium_flights,
                    c5_medium_flights,
                    c6_medium_flights)
short_flights <- c(c0_short_flights,
                   c1_short_flights,
                   c2_short_flights,
                   c3_short_flights,
                   c4_short_flights,
                   c5_short_flights,
                   c6_short_flights)
skip_flights <- c(c0_skip_flights,
                 c1_skip_flights,
                 c2_skip_flights,
                 c3_skip_flights,
                 c4_skip_flights,
                 c5_skip_flights,
                 c6_skip_flights)
flight <- data.frame(session,
                     chann,
                     total_distance, 
                     total_duration, 
                     total_avg_speed,
                     total_peak_speed,
                     marathon_flights,
                     long_flights,
                     medium_flights,
                     short_flights,
                     skip_flights,
                     trial_length)

 #flight
write.table(flight, "flight.csv", sep = ",", col.names = !file.exists("flight.csv"), append = T, row.names = F)
#write.csv(flight,"flight.csv", row.names = T)
#plot(speed_avg ~ log(dist), data = channel_5_summary, col = "blue")
#points(general_avg ~ log(dist), data = channel_5_summary, col = "red")
#abline(0,1)
#summary(lm(speed_avg~general_avg, data = channel_5_summary))
