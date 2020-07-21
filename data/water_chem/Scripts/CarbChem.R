#Title: SewaterChem-Analysis - 2019 Point Whitney broodstock experiment
#Project: FFAR
#Author: Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20200720
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
library(dplyr)
library(seacarb) #version: 3.2 Date/Publication: 2017-06-19 Depends: R (>= 2.10), oce, gsw Imports: NA

#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/paper-GeoduckReproDev_pH/data/water_chem/") #set working

#pH calibration (Tris)--------------------------
path <-("C:/Users/samjg/Documents/My_Projects/paper-GeoduckReproDev_pH/data/water_chem/Titrator/pH_Calibration_Files")
file.names<-list.files(path = path, pattern = "csv$")
file.names

pH.cals <- data.frame(matrix(NA, nrow=length(file.names), ncol=3, dimnames=list(file.names,c("Date", "Intercept", "Slope")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Calib.Data <-read.table(file.path(path,file.names[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
  model <-lm(mVTris ~ TTris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
  coe <- coef(model) #extracts the coeffecients
  summary(model)$r.squared
  plot(Calib.Data$mVTris, Calib.Data$TTris)
  pH.cals[i,2:3] <- coe #inserts them in the dataframe
  pH.cals[i,1] <- substr(file.names[i],1,8) #stores the file name in the Date column
}
colnames(pH.cals) <- c("Calib.Date",  "Intercept",  "Slope") #rename columns
pH.cals #view data

#constants for use in pH calculation 
R <- 8.31447215 #gas constant in J mol-1 K-1 
F <-96485.339924 #Faraday constant in coulombs mol-1

#read in probe measurements of pH, temperature, and salinity from tanks
daily <- read.csv("Titrator/Daily_Temp_pH_Sal.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

#merge with Seawater chemistry file
SW.chem <- merge(pH.cals, daily, by="Calib.Date")

#calculate pH.Total
mvTris <- SW.chem$Temperature*SW.chem$Slope+SW.chem$Intercept #calculate the mV of the tris standard using the temperature mv relationships in the measured standard curves 
STris<-34.5 #salinity of the Tris
phTris<- (11911.08-18.2499*STris-0.039336*STris^2)*(1/(SW.chem$Temperature+273.15))-366.27059+ 0.53993607*STris+0.00016329*STris^2+(64.52243-0.084041*STris)*log(SW.chem$Temperature+273.15)-0.11149858*(SW.chem$Temperature+273.15) #calculate the pH of the tris (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
SW.chem$pH.Total<-phTris+(mvTris/1000-SW.chem$pH_mV/1000)/(R*(SW.chem$Temperature+273.15)*log(10)/F) #calculate the pH on the total scale (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)

#output figures
pdf("C:/Users/samjg/Documents/My_Projects/paper-GeoduckReproDev_pH/data/water_chem/Output/Daily_Treatment_Measures_1.pdf")
par(mfrow=c(3,1))
plot(SW.chem$Treatment, SW.chem$Temperature, xlab="Treatment", ylab="Temperature°C", ylim=c(5,14), las=2)
plot(SW.chem$Treatment, SW.chem$pH.Total, xlab="Treatment", ylab="pH Total Scale", ylim=c(6.5,8.2), las=2)
plot(SW.chem$Treatment, SW.chem$Salinity, xlab="Treatment", ylab="Salinity psu", ylim=c(26,32), las=2)
dev.off()

pdf("C:/Users/samjg/Documents/My_Projects/paper-GeoduckReproDev_pH/data/water_chem/Output/Daily_Treatment_Measures_2.pdf")
par(mfrow=c(3,1))
plot(SW.chem$Sample.ID, SW.chem$Temperature, xlab="Tank", ylab="Temperature°C", ylim=c(5,14),las=2)
plot(SW.chem$Sample.ID, SW.chem$pH.Total, xlab="Tank", ylab="pH Total Scale", ylim=c(6.5,8.2),las=2)
plot(SW.chem$Sample.ID, SW.chem$Salinity, xlab="Tank", ylab="Salinity psu", ylim=c(26,32),las=2)
dev.off()

##### DISCRETE TA CALCULATIONS #####
TA <- read.csv("Titrator/Cumulative_TA_Output.csv", header=TRUE, sep=",", na.strings="NA")  #read in  TA results

##### SEAWATER CHEMISTRY ANALYSIS FOR DISCRETE MEASUREMENTS#####
#Seawater chemistry table from simultaneous TA, pH, temperature and salinity measurements
#merge calculated pH and daily measures with TA data and run seacarb
SW.chem$Date_Sample.ID <- paste(SW.chem$Date, SW.chem$Sample.ID, sep='_') #generate new row with concatenated sample id
TA$Calib.Date <- substr(TA$Sample.ID, 1,8)
TA$ID <- as.numeric(substr(TA$Sample.ID, 10,10)) # call substring as numeric
TA <- na.omit(TA) # ommit the TA of CRM data
TA$Sample.ID_2 <- paste('TANK',TA$ID, sep='')
TA$Date_Sample.ID <- paste((substr(TA$Sample.ID,1,8)), TA$Sample.ID_2, sep='_') #generate new row with concatenated sample id

SW.chem <- merge(SW.chem,TA, by="Date_Sample.ID", all = TRUE, sort = T) #merge seawater chemistry with total alkalinity

#Calculate CO2 parameters using seacarb
carb.output <- carb(flag=8, var1=SW.chem$pH.Total, var2=SW.chem$TA/1000000, S= SW.chem$Salinity, T=SW.chem$Temperature, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d") #calculate seawater chemistry parameters using seacarb
carb.output$ALK <- carb.output$ALK*1000000 #convert to µmol kg-1
carb.output$CO2 <- carb.output$CO2*1000000 #convert to µmol kg-1
carb.output$HCO3 <- carb.output$HCO3*1000000 #convert to µmol kg-1
carb.output$CO3 <- carb.output$CO3*1000000 #convert to µmol kg-1
carb.output$DIC <- carb.output$DIC*1000000 #convert to µmol kg-1
carb.output <- carb.output[,-c(1,4,5,8,10:13,19)] #subset variables of interest
carb.output <- cbind(SW.chem$Date,  SW.chem$Date_Sample.ID,  SW.chem$Treatment, carb.output) #combine the sample information with the seacarb output
colnames(carb.output) <- c("Date",  "Tank",  "Treatment",	"Salinity",	"Temperature", "pH",	"CO2",	"pCO2","HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat") #Rename columns to describe contents

write.table(carb.output, "C:/Users/samjg/Documents/My_Projects/paper-GeoduckReproDev_pH/data/water_chem/Output/Seawater_chemistry_table_Output_All.csv", sep=",", row.names = FALSE) #save data




