#Title: Pgenerosa histology staging analysis
#Project: FFAR
#Author: Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20200716
#See Readme file for details

rm(list=ls())

# Load packages and pacage version/date/import/depends info
library(dplyr)          # Version 0.7.6, Packaged: 2018-06-27, Depends: R (>= 3.1.2)Imports: assertthat (>= 0.2.0), bindrcpp (>= 0.2.0.9000), glue (>=1.1.1), magrittr (>= 1.5), methods, pkgconfig (>= 2.0.1), R6(>= 2.2.2), Rcpp (>= 0.12.15), rlang (>= 0.2.0), tibble (>=1.3.1), tidyselect (>= 0.2.3), utils
library(ggplot2)        # Version 2.2.1, Packaged: 2016-12-30, Depends: R (>= 3.1)Imports: digest, grid, gtable (>= 0.1.1), MASS, plyr (>= 1.7.1),reshape2, scales (>= 0.4.1), stats, tibble, lazyeval
library(tidyr)
library(ggpubr)
library(viridis)
library(hrbrthemes)
library(gridExtra)

#set working directory--------------------------------------------------------------------------------------------------
setwd("C:/Users/samjg/Documents/My_Projects/Pgenerosa_histology/RAnalysis/") #set working

# upload data
histology<-read.csv("Data/Staging/Hist_Staging_Trigg_Crandall.csv", header=T, sep=",", na.string="NA", as.is=T) 
histology # view data
histology$Stage <- ifelse(is.na(histology$Stage_Crandall), histology$Stage_Trigg, histology$Stage_Crandall) # call Grace's Stages (Stage_Crandall) - however when NA (photos analyzed by Shelly only) replace with Stage_Trigg

# RE: 'Stage' and 'Staging_number' columns
# Immature: very early active (1),
# early active (2)) 
# Mature: late active (3)
# ripe (4))
# Spent (5)

hist_table <- histology %>% dplyr::select(c('Tank','Treatment','Date', 'Sex', 'Slide.Number', 'Geoduck.ID', 'Stage'))
hist_table.OM <- na.omit(hist_table) # ommit the single NA value
hist_table.OM$Stage <- as.character(hist_table.OM$Stage)
# calcaulate the proportions of staging number by date and treatment 
prop_hist_stage <-
  hist_table.OM %>% 
  group_by(Date, Treatment, Sex, Stage) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
prop_hist_stage$prop <- prop_hist_stage$freq*100

prop_Jan23_Female <- prop_hist_stage %>% dplyr::filter('123' %in% Date) %>%  dplyr::filter('F' %in% Sex)
prop_Jan23_Male <- prop_hist_stage %>% dplyr::filter('123' %in% Date) %>%  dplyr::filter('M' %in% Sex)
prop_Feb21_Female <- prop_hist_stage %>% dplyr::filter('221' %in% Date) %>%  dplyr::filter('F' %in% Sex)
prop_Feb21_Male <- prop_hist_stage %>% dplyr::filter('221' %in% Date) %>%  dplyr::filter('M' %in% Sex)

### separate datasets for six, date, and treatment 
# female
prop_Jan23_Female_Amb <- prop_Jan23_Female %>%
  dplyr::filter('Ambient' %in% Treatment)

prop_Jan23_Female_Low <- prop_Jan23_Female %>%
  dplyr::filter('Low' %in% Treatment)

prop_Feb21_Female_Amb <- prop_Feb21_Female %>%
  dplyr::filter('Ambient' %in% Treatment)

prop_Feb21_Female_Low <- prop_Feb21_Female %>%
  dplyr::filter('Low' %in% Treatment)
# Male
prop_Jan23_Male_Amb <- prop_Jan23_Male %>%
  dplyr::filter('Ambient' %in% Treatment)

prop_Jan23_Male_Low <- prop_Jan23_Male %>%
  dplyr::filter('Low' %in% Treatment)

prop_Feb21_Male_Amb <- prop_Feb21_Male %>%
  dplyr::filter('Ambient' %in% Treatment)

prop_Feb21_Male_Low <- prop_Feb21_Male %>%
  dplyr::filter('Low' %in% Treatment)


# STACKED BAR PLOT 

# Stacked
prop_hist_stage$Date_Treatment <- paste(prop_hist_stage$Date, prop_hist_stage$Treatment, sep ="_") # new column merges date and treat for stacked bar plots
prop_female <- prop_hist_stage %>% dplyr::filter(Sex %in% 'F') # new dataset to plot JUST female histology
levels(prop_female$Stage) <- c("7", "6", "5")

prop_male <- prop_hist_stage %>% dplyr::filter(Sex %in% 'M')  # new dataset to plot JUST male histology

StackedBar_FEMALE <- ggplot(prop_female, aes(fill=Treatment, y= prop, x=Treatment, alpha = Stage)) + 
                      scale_fill_manual(values = c("#00BFC4","#F8766D")) +
                      geom_bar(position="stack", stat="identity",colour="white") +
                      ggtitle("Female histology staging") +
                      # theme_ipsum() +
                      theme_classic() +
                      xlab("Date_Treatment") +
                      ylab("number of samples") +
                       facet_wrap(~ Date+Sex)
StackedBar_FEMALE# view plot

StackedBar_MALE <- ggplot(prop_male, aes(fill=Treatment, y= prop, x=Treatment, alpha = Stage)) + 
                      scale_fill_manual(values = c("#00BFC4","#F8766D")) +
                      geom_bar(position="stack", stat="identity",colour="white") +
                      ggtitle("Male histology staging") +
                      # theme_ipsum() +
                      theme_classic() +
                      xlab("Date_Treatment") +
                      ylab("number of samples") +
                      facet_wrap(~ Date)
StackedBar_MALE# view plot

StackedBar_ALL <- ggplot(prop_hist_stage, aes(fill=Stage, y= prop, x=Treatment)) + 
                    geom_bar(position="stack", stat="identity",colour="white") +
                    #scale_fill_viridis(discrete = T) +
                    scale_fill_grey(start=0.8, end=0.2) +
                    ggtitle("Female histology staging") +
                    # theme_ipsum() +
                    theme_classic() +
                    xlab("Date_Treatment") +
                    ylab("number of samples") +
                    facet_wrap(~ Sex + Date)
StackedBar_ALL# view plot

#  SAVE
stacked_pop <- ggarrange(StackedBar_FEMALE, StackedBar_MALE,  ncol = 1, nrow = 2)
ggsave(file="Stacked_barplots_Staging_PROPORTION.pdf", stacked_pop, width = 12, height = 24, units = c("in")) 
ggsave(file="Stacked_barplots_Staging.pdf", StackedBar_ALL, width = 24, height = 12, units = c("in")) 


# DONUT PLOTS
F_amb_123 <- ggdonutchart(prop_Jan23_Female_Amb, "prop", label = "Stage",
             lab.pos = "in", lab.size = 10, lab.font = "white",
             title = "Female_Ambient_0123", 
             fill = "Stage", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00"))

F_low_123 <- ggdonutchart(prop_Jan23_Female_Low, "prop", label = "Stage",
             lab.pos = "in", lab.size = 10, lab.font = "white",
             title = "Female_Low_0123", 
             fill = "Stage", color = "white", palette = c("#56B4E9", "#0072B2"))
             
F_amb_221 <- ggdonutchart(prop_Feb21_Female_Amb, "prop", label = "Stage",
             lab.pos = "in", lab.size = 10, lab.font = "white",
             title = "Female_Ambient_0221", 
             fill = "Stage", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00", "#999999"))

F_low_221 <- ggdonutchart(prop_Feb21_Female_Low, "prop", label = "Stage",
             lab.pos = "in", lab.size = 10, lab.font = "white",
             title = "Female_Low_0221", 
             fill = "Stage", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00"))


M_amb_123 <- ggdonutchart(prop_Jan23_Male_Amb, "prop", label = "Stage",
                          lab.pos = "in", lab.size = 10, lab.font = "white",
                          title = "Male_Ambient_0123", 
                          fill = "Stage", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00", "#999999"))

M_low_123 <- ggdonutchart(prop_Jan23_Male_Low, "prop", label = "Stage",
                          lab.pos = "in", lab.size = 10, lab.font = "white",
                          title = "Male_Low_0123", 
                          fill = "Stage", color = "white", palette = c("#56B4E9","#CC79A7"))

M_amb_221 <- ggdonutchart(prop_Feb21_Male_Amb, "prop", label = "Stage",
                          lab.pos = "in", lab.size = 10, lab.font = "white",
                          title = "Male_Ambient_0221",                
                          fill = "Stage", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00", "#999999"))

M_low_221 <- ggdonutchart(prop_Feb21_Male_Low, "prop", label = "Stage",
                          lab.pos = "in", lab.size = 10, lab.font = "white",
                          title = "Male_Low_0221", 
                          fill = "Stage", color = "white", palette = c("#56B4E9", "#0072B2", "#D55E00", "#999999"))


ggarrange(M_amb_123, M_low_123,F_amb_123, F_low_123,
          M_amb_221, M_low_221, F_amb_221, F_low_221,  ncol = 4, nrow = 2)
