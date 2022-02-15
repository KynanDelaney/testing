# Script for combining flightmill and beetle data ---------------

# A script to merge anonymous chamber information and beetle id and measurements.
# Reshapes wide-format individual data to long-format.
# Author: Kynan Delaney
# Date: 14/02/2022


# Libraries ---------------------------------------------------------------
rm(list = ls())

# for data cleaning and plotting
library(tidyverse)
library(lubridate)


# load in data ------------------------------------------------------------

# Read-in main beetle measurements from "Size_Manipulation_Final.csv".
beetles<- read.csv(file.choose(),
                   stringsAsFactors = FALSE, header = T)

# Make a more informative individual ID. Remove extraneous columns.
beetles <- beetles %>%
  mutate(ID = paste0(Family,ID)) %>%
  select(-Weight,-LastObserved,-Death.Removal,-Notes)

beetles$ID <- factor(beetles$ID)

# Converts a wide format of three trials to repeated-measures long format.
beetles <- reshape(beetles, 
                   varying=list(Trial.Date= c("X1st.Trial","X2nd.Trial", "X3rd.Trial"), 
                                Chamber= c("X1st.Chamber","X2nd.Chamber","X3rd.Chamber"), 
                                StartTime= c("X1st.StartTime","X2nd.StartTime","X3rd.StartTime"), 
                                Session= c("X1st.Session","X2nd.Session","X3rd.Session"), 
                                Comments= c("X1st.Comments","X2nd.Comments","X3rd.Comments") ), 
                   v.names=c("Trial.Date", "Chamber", "StartTime", "Session","Comments"), 
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
  mutate(age_at_flight = as.numeric(difftime(dmy(Trial.Date), dmy(Eclosion))))

# Read-in un-paired flight  measurements from "Chamber_information_0.1_to_3_thresholds".
flight<- read.csv(file.choose(),
                  stringsAsFactors = FALSE, header = T)


# combining beetle id and measurements with trial, session and chamber info.
flight <- list(beetles,flight) %>% reduce(full_join, by = c("Session", "Chamber" = "Channel"))

# family shouldn't be missing from any rows.
flight<- flight %>%
  filter(!is.na(Family))

# Removing data that either came from the pilot, or failed week.
flight <- flight %>%
  filter(Block != "A" & Block != "E")
flight$Sex <- ifelse(flight$Sex == "", "M", flight$Sex)

# outputting all this cleaned/paired data to a useable csv.
write.csv(flight, file = "Combined_Flight_final.csv")
