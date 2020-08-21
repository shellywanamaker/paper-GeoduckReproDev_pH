#Title: Apex Tables - 2019 Point Whitney broodstock experiment
#Project: FFAR
#Author: Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20200716
#See Readme file for details

rm(list=ls())
# Install packages if not already in your library
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("ggplot2" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggplot2') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
# Load packages and pacage version/date/import/depends info
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(lubridate)


#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/paper-GeoduckReproDev_pH/data/water_chem") #set working

### CONICAL Seawater chemistry Data - Analysis, Graphs, Tables (APEX DATA) ####

#CONTINUOUS EXPERIMENTAL APEX DATA 
#Load Apex Data 
APEX_1<-read.csv("Apex/Apex_data_20181026-20181107.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_1.dat <- APEX_1 %>% dplyr::select(c("date","probe.value","probe.name.1", "probe.value.1", "probe.name.4", "probe.value.4",
                                         "probe.name.5", "probe.value.5", "probe.name.23", "probe.value.23", "probe.name.24", "probe.value.24",
                                         "probe.name.25", "probe.value.25", "probe.name.26", "probe.value.26", "probe.name.27", "probe.value.27")) %>% na.omit()

APEX_2<-read.csv("Apex/Apex_data_20181107-20181112.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_2.dat <- APEX_2 %>% dplyr::select(c("date","probe.value","probe.name.1", "probe.value.1", "probe.name.4", "probe.value.4",
                                         "probe.name.5", "probe.value.5", "probe.name.23", "probe.value.23", "probe.name.24", "probe.value.24",
                                         "probe.name.25", "probe.value.25", "probe.name.26", "probe.value.26", "probe.name.27", "probe.value.27")) %>% na.omit()

APEX_3<-read.csv("Apex/Apex_data_20181107-20181128.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_3.dat <- APEX_3 %>% dplyr::select(c("date","probe.value","probe.name.1", "probe.value.1", "probe.name.4", "probe.value.4",
                                         "probe.name.5", "probe.value.5", "probe.name.23", "probe.value.23", "probe.name.24", "probe.value.24",
                                         "probe.name.25", "probe.value.25", "probe.name.26", "probe.value.26", "probe.name.27", "probe.value.27")) %>% na.omit()


APEX_4<-read.csv("Apex/Apex_data_20181115-20181212.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_4.dat <- APEX_4 %>% dplyr::select(c("date","probe.value","probe.name.1", "probe.value.1", "probe.name.4", "probe.value.4",
                                         "probe.name.5", "probe.value.5", "probe.name.23", "probe.value.23", "probe.name.24", "probe.value.24",
                                         "probe.name.25", "probe.value.25", "probe.name.26", "probe.value.26", "probe.name.27", "probe.value.27")) %>% na.omit()


APEX_5<-read.csv("Apex/Apex_data_20181201-20181219.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_5.dat <- APEX_5 %>% dplyr::select(c("date","probe.value","probe.name.1", "probe.value.1", "probe.name.4", "probe.value.4",
                                         "probe.name.5", "probe.value.5", "probe.name.23", "probe.value.23", "probe.name.24", "probe.value.24",
                                         "probe.name.25", "probe.value.25", "probe.name.26", "probe.value.26", "probe.name.27", "probe.value.27")) %>% na.omit()


APEX_6<-read.csv("Apex/Apex_data_20181215-20190103.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_6.dat <- APEX_6 %>% dplyr::select(c("date","probe.value","probe.name.1", "probe.value.1", "probe.name.4", "probe.value.4",
                                         "probe.name.5", "probe.value.5", "probe.name.23", "probe.value.23", "probe.name.24", "probe.value.24",
                                         "probe.name.25", "probe.value.25", "probe.name.26", "probe.value.26", "probe.name.27", "probe.value.27")) %>% na.omit()
                          


APEX_7<-read.csv("Apex/Apex_data_20190101-20190114.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_7.dat <- APEX_7 %>% na.omit() %>% dplyr::select(c("date","probe.value","probe.name.1", "probe.value.1", "probe.name.4", "probe.value.4",
                                         "probe.name.5", "probe.value.5", "probe.name.23", "probe.value.23", "probe.name.24", "probe.value.24",
                                         "probe.name.25", "probe.value.25", "probe.name.26", "probe.value.26", "probe.name.27", "probe.value.27")) %>% na.omit()


APEX_8<-read.csv("Apex/Apex_data_20190101-20190201.csv", header=T, sep=",", na.string="NA", as.is=T) 
APEX_8.dat <- APEX_8 %>% dplyr::select(c("date","probe.value","probe.name.1", "probe.value.1", "probe.name.4", "probe.value.4",
                                         "probe.name.5", "probe.value.5", "probe.name.23", "probe.value.23", "probe.name.24", "probe.value.24",
                                         "probe.name.25", "probe.value.25", "probe.name.26", "probe.value.26", "probe.name.27", "probe.value.27")) %>% na.omit()


colnames(APEX_1)
colnames(APEX_2)
colnames(APEX_3)
colnames(APEX_8)

APEX_data<- do.call("rbind", list(APEX_1.dat, APEX_3.dat, APEX_4.dat, APEX_5.dat, APEX_6.dat, APEX_7.dat, APEX_8.dat)) # bind all data together
APEX_data$date_form <-as.POSIXct(APEX_data$date, format="%m/%d/%Y %H:%M:%S") #convert date format
APEX_data$date_form 
APEX_data$X <- 1:nrow(APEX_data) # make new column for cumulative number of rows

colnames(APEX_data)
plot(APEX_data$date_form, APEX_data$probe.value.1) #pH
plot(APEX_data$date_form, APEX_data$probe.value.4) #temp 2 -- not used
plot(APEX_data$date_form, APEX_data$probe.value.5) # pH 2  -- not used
plot(APEX_data$date_form, APEX_data$probe.value.23) # temp 3
plot(APEX_data$date_form, APEX_data$probe.value.24) # pH 3
plot(APEX_data$date_form, APEX_data$probe.value.25) # temp 4
plot(APEX_data$date_form, APEX_data$probe.value.26) # pH 4
plot(APEX_data$date_form, APEX_data$probe.value.27) # pH 6 -- looks like temperature
