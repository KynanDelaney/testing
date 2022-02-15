# Template Script for extracting/processing flightmill data ---------------

# Author: Kynan Delaney
# Date: 10/11/2021

# A script to streamline extracting data from the Pis for initial cleaning.
# I define the database connection, simulate data, insert it into the database
# and recover it in a processable format.

# IMPORTANT NOTE: DEVELOPMENT DATABASE CAN BE OVERWRITTEN. PRODUCTION DATABASE
# SHOULD BE RESTRICTED FOR READ-ONLY REMOTE ACCESS.




rm(list = ls())

# Libraries and options ---------------------------------------------------

library(RPostgres) 
library(RPostgreSQL)
library(dplyr)
library(readr)
library(plyr)
library(purrr)
library(lubridate)
library(data.table)
library(chron)

options("digits.secs"=6)


# Defining connection to database -----------------------------------------

name <- "blackmill"
ip <- "172.20.42.199"

# Details for connecting to database. "Host" may change if ip-address changes.
con <- dbConnect(RPostgres::Postgres(), dbname = name, host=ip, port='5432', user="pi", password = "raspberry")

# Returns list of Tables stored in this database. Should only return ONE.
mill_colour <- dbListTables(con) 

# Print out a sample of the table. Basically "head(mydata)". Takes ages to run.
#dbReadTable(con, mill_colour)


# Accessing, Extracting and Exporting data --------------------------------

#-------------------------------------------------------#
# DOUBLE-CHECK the destination folder for all the data. #
#-------------------------------------------------------#

# Dplyr approach to using database data.
working <- tbl(con, mill_colour)

# Create a list of unique session IDs to later filter and export.
session_list <- working %>% distinct(session) %>% pull()

# A loop for writing CSV for each unique session ID. 
for (i in session_list){
  working %>%
    filter(session == i) %>%
    write.csv(paste0(mill_colour,"_session_", i, ".csv"))}

# Disconnecting from database when finished.
dbDisconnect(con) 


