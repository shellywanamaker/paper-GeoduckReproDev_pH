#Title: ChemTables
#Project: FFAR
#Author: Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20200720
#See Readme file for details

rm(list=ls())

# Install packages if not already in your library
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("ggplot2" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggplot2') 
if ("ggpubr" %in% rownames(installed.packages()) == 'FALSE') install_github('ggpubr') 
if ("Rmisc" %in% rownames(installed.packages()) == 'FALSE') install.packages('Rmisc') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("lsmeans" %in% rownames(installed.packages()) == 'FALSE') install.packages('lsmeans') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("reshape" %in% rownames(installed.packages()) == 'FALSE') install.packages('reshape') 
if ("multcompView" %in% rownames(installed.packages()) == 'FALSE') install.packages('multcompView') 
if ("lmtest" %in% rownames(installed.packages()) == 'FALSE') install.packages('lmtest') 
if ("car" %in% rownames(installed.packages()) == 'FALSE') install.packages('car') 

# Load packages and pacage version/date/import/depends info
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(ggpubr)         # Version: 0.1.8 Date: 2018-08-30, Depends: R (>= 3.1.0), ggplot2, magrittrImports: ggrepel, grid, ggsci, stats, utils, tidyr, purrr, dplyr(>=0.7.1), cowplot, ggsignif, scales, gridExtra, glue, polynom
library(Rmisc)          # Version: 1.5 Packaged: 2013-10-21, Depends: lattice, plyr
library(plotrix)        # Version: 3.7-4, Date/Publication: 2018-10-03
library(lsmeans)        # Version: 2.27-62, Date/Publication: 2018-05-11, Depends: methods, R (>= 3.2)
library(gridExtra)      # Version: 2.3, Date/Publication: 2017-09-09, Imports: gtable, grid, grDevices, graphics, utils
library(reshape)        # Version: 0.8.7, Date/Publication: 2017-08-06, Depends: R (>= 2.6.1) Imports: plyr
library(multcompView)   # Version: 0.1-7, Date/Publication: 2015-07-31, Imports: grid
library(Rmisc)
library(lmtest)
library(car)

#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/Github_repositories/paper-GeoduckReproDev_pH/data/water_chem/") #set working

#read in CarbChem.R output Seawater chemistry csv file
chem <- read.csv("Output/Seawater_chemistry_table_Output_All.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
chem # view data
chem.om <- na.omit(chem)
chem.om$Tank <- substr(chem.om$Tank,14,14)

#####################   make long table ################################################  #
chem.treat <- chem.om %>% select(-("Tank")) # drop "Tank" to create long table via meltin next line
chem.treat.LONG <- melt(chem.treat, id.vars=c("Date", "Treatment")) # uses tidyr to make a long table from wide
chem.tank <- chem.om %>% select(-("Treatment")) # drop "Tank" to create long table via meltin next line
chem.tank.LONG <- melt(chem.tank, id.vars=c("Date", "Tank")) # uses tidyr to make a long table from wide

################### mean and SEM tables for TREATMENT ################### # 
TABLE.treat <- ddply(chem.treat.LONG, c("Treatment","variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                   N = length(na.omit(value)), #count the sample size removing NA
                   mean = mean(value), #calculate average 
                   sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean
TABLE.treat$mean <- signif(TABLE.treat$mean,digits=3) # reduce number of sig digits to three
TABLE.treat$sem <- signif(TABLE.treat$sem,digits=3) # reduce number of sig digits to three
TABLE.treat_wide <- reshape(TABLE.treat, idvar = c("Treatment", "N"), timevar = "variable", direction = "wide")
TABLE.treat_wide # view table

################### mean and SEM tables for TREATMENT ################### # 
TABLE.tank<- ddply(chem.tank.LONG, c("Tank","variable"), summarise, #Calculate descriptive stats by Treatment and exposure 
                     N = length(na.omit(value)), #count the sample size removing NA
                     mean = mean(value), #calculate average 
                     sem = sd(value)/sqrt(N)) # calcualte the standard error of the mean
TABLE.tank$mean <- signif(TABLE.tank$mean,digits=3) # reduce number of sig digits to three
TABLE.tank$sem <- signif(TABLE.tank$sem,digits=3) # reduce number of sig digits to three
TABLE.tank_wide <- reshape(TABLE.tank, idvar = c("Tank", "N"), timevar = "variable", direction = "wide")
TABLE.tank_wide # view table


############################ FINAL TABLES WITH ST.ERROR ###################################### #
# final table TREATMENT
TABLE.treatment_FINAL <- data.frame(matrix(nrow = 3, ncol = 1))
TABLE.treatment_FINAL$N <- TABLE.treat_wide$N
TABLE.treatment_FINAL$Treatment <- TABLE.treat_wide$Treatment
TABLE.treatment_FINAL$Salinity <- paste(TABLE.treat_wide$mean.Salinity, TABLE.treat_wide$sem.Salinity, sep=" ? ")
TABLE.treatment_FINAL$Temperature <- paste(TABLE.treat_wide$mean.Temperature, TABLE.treat_wide$sem.Temperature, sep=" ? ")
TABLE.treatment_FINAL$pH <- paste(TABLE.treat_wide$mean.pH, TABLE.treat_wide$sem.pH, sep=" ? ")
TABLE.treatment_FINAL$CO2 <- paste(TABLE.treat_wide$mean.CO2, TABLE.treat_wide$sem.CO2, sep=" ? ")
TABLE.treatment_FINAL$pCO2 <- paste(TABLE.treat_wide$mean.pCO2, TABLE.treat_wide$sem.pCO2, sep=" ? ")
TABLE.treatment_FINAL$HCO3 <- paste(TABLE.treat_wide$mean.HCO3, TABLE.treat_wide$sem.HCO3, sep=" ? ")
TABLE.treatment_FINAL$CO3 <- paste(TABLE.treat_wide$mean.CO3, TABLE.treat_wide$sem.CO3, sep=" ? ")
TABLE.treatment_FINAL$DIC <- paste(TABLE.treat_wide$mean.DIC, TABLE.treat_wide$sem.DIC, sep=" ? ")
TABLE.treatment_FINAL$TA <- paste(TABLE.treat_wide$mean.TA, TABLE.treat_wide$sem.TA, sep=" ? ")
TABLE.treatment_FINAL$Aragonite.Sat <- paste(TABLE.treat_wide$mean.Aragonite.Sat, TABLE.treat_wide$sem.Aragonite.Sat, sep=" ? ")
TABLE.treatment_FINAL <- TABLE.treatment_FINAL[,-1] # view table

# final table TANK
TABLE.tank_FINAL <- data.frame(matrix(nrow = 6, ncol = 1))
TABLE.tank_FINAL$N <- TABLE.tank_wide$N
TABLE.tank_FINAL$Tank <- TABLE.tank_wide$Tank
TABLE.tank_FINAL$Salinity <- paste(TABLE.tank_wide$mean.Salinity, TABLE.tank_wide$sem.Salinity, sep=" ? ")
TABLE.tank_FINAL$Temperature <- paste(TABLE.tank_wide$mean.Temperature, TABLE.tank_wide$sem.Temperature, sep=" ? ")
TABLE.tank_FINAL$pH <- paste(TABLE.tank_wide$mean.pH, TABLE.tank_wide$sem.pH, sep=" ? ")
TABLE.tank_FINAL$CO2 <- paste(TABLE.tank_wide$mean.CO2, TABLE.tank_wide$sem.CO2, sep=" ? ")
TABLE.tank_FINAL$pCO2 <- paste(TABLE.tank_wide$mean.pCO2, TABLE.tank_wide$sem.pCO2, sep=" ? ")
TABLE.tank_FINAL$HCO3 <- paste(TABLE.tank_wide$mean.HCO3, TABLE.tank_wide$sem.HCO3, sep=" ? ")
TABLE.tank_FINAL$CO3 <- paste(TABLE.tank_wide$mean.CO3, TABLE.tank_wide$sem.CO3, sep=" ? ")
TABLE.tank_FINAL$DIC <- paste(TABLE.tank_wide$mean.DIC, TABLE.tank_wide$sem.DIC, sep=" ? ")
TABLE.tank_FINAL$TA <- paste(TABLE.tank_wide$mean.TA, TABLE.tank_wide$sem.TA, sep=" ? ")
TABLE.tank_FINAL$Aragonite.Sat <- paste(TABLE.tank_wide$mean.Aragonite.Sat, TABLE.tank_wide$sem.Aragonite.Sat, sep=" ? ")
TABLE.tank_FINAL <- TABLE.tank_FINAL[,-1] # view table

# save output table
write.table(TABLE.treatment_FINAL,"C:/Users/samjg/Documents/My_Projects/paper-GeoduckReproDev_pH/data/water_chem/Output/Chem.Table.TREATMENT.csv",sep=",", row.names=FALSE)  # write table to output folder
write.table(TABLE.tank_FINAL,"C:/Users/samjg/Documents/My_Projects/paper-GeoduckReproDev_pH/data/water_chem/Output/Chem.Table.TANK.csv",sep=",", row.names=FALSE)  # write table to output folder
