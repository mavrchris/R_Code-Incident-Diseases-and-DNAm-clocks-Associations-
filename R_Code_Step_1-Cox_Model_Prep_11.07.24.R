
# Step 1 #

### RENV to ensure reproducibility ###
install.packages("renv")
library(renv)

# create lock file #
renv::init()

# Update locked file with any new packages used #
renv::snapshot()

# Download +install all used packages in lock file (with correct version) #
renv::restore()


### Set Directory ###

# UNI PC #
#setwd("Y:/marioni-lab/Christos/GS data") #

# MY PC# 

setwd ("Y:/Christos/GS data") 


### load in : ###

### Packages ###
library(dplyr)

### Epigenetic Clocks ###
clock_data<- read.csv("2024-03-29_clock_output.csv")

### DuneDin PoAm  ###
poam_clock_data <- read.csv("2024-03-29_poam.csv")

### load in disease, Covariate and Mortality data  ###
disease_data<- read.csv("2024-03-29_diseases.csv")
cov_data<- read.csv("2024-04-19_covariates.csv")
mortality_data<- read.csv("2024-04-19_deaths.csv")

### Load in file linking ids to meth dis ###
GSK_MAP <- readRDS("GS20K_Targets_18869.rds")

### load in appt records ###
appt_rds <- readRDS("202-04-26_gs_appt.rds")



##1##

## Filter Clocks to relevant predictors ##
filt_clock_data<- select(clock_data, c(SampleID,DNAmAgeHannum, DNAmPhenoAge, DNAmGrimAge,
                                       DNAmAgeHannumAdjAge, DNAmPhenoAgeAdjAge, DNAmGrimAgeAdjAge))

### Now, Merge clock data (filteres + poam + anchor point) with anchor point = sample ID , by.x= "SampleID" in clock_data, by.y= "x" (name of id)in poam ###
merged_clock_data <- merge(filt_clock_data,poam_clock_data, by.x= "SampleID", by.y = "X", all.x = TRUE)

### merge appt data to GSK ###
appt_GSK <- merge(appt_rds, GSK_MAP[,c("X","Sample_Name")], by.x= "id", by.y= "Sample_Name")

merged_clock_data_1 <- merge(merged_clock_data,appt_GSK, by.x= "SampleID", by.y= "X" )

### merge merged_clock_data with mortality data ### # if columns are same, no need for by.x or by.y just by = id ##
merged_clock_data_2<- merge(merged_clock_data_1, mortality_data, by.x= "id", by.y= "id", all.x= TRUE )

#Check NA= alive (not in mortality data)# 
merged_clock_data_2_alive <- merged_clock_data_2 %>% 
  filter(!is.na(dod_ym))

### merge merged_clocks_2 with cov, USING GSK as anchors ###
merged_clock_data_3 <- merge(merged_clock_data_2, cov_data, by.x= "id", by.y= "id")

### Calculate Time to censor: ###

### first ### 

### Substr APPT dates ###
merged_clock_data_3$GS_yoa <- substr(merged_clock_data_3$ym, 1,4) %>% as.numeric()
merged_clock_data_3$GS_moa <- substr(merged_clock_data_3$ym, 5,6) %>% as.numeric()

### substr YOD dates ### 
merged_clock_data_3$GS_yod <- substr(merged_clock_data_3$dod_ym, 1,4) %>% as.numeric()
merged_clock_data_3$GS_mod <- substr(merged_clock_data_3$dod_ym, 5,6) %>% as.numeric()

# recheck Check NA= as people have not died# 
merged_clock_data_3_alive <- merged_clock_data_3 %>% 
  filter(!is.na(GS_yod))

### new column : 202204- GS APPT ###
merged_clock_data_3$cutoff_minus_GSAPPT <- (2022 - merged_clock_data_3$GS_yoa) + ((04- merged_clock_data_3$GS_moa)/12) %>% as.numeric()

### new column : dod- GSappt ### 
merged_clock_data_3$DOD_minus_GSAPPT <- (merged_clock_data_3$GS_yod - merged_clock_data_3$GS_yoa) + ((merged_clock_data_3$GS_mod - merged_clock_data_3$GS_moa)/12) %>% as.numeric()

## New column T_Censor ##
merged_clock_data_3$t_censor <- ifelse(!is.na(merged_clock_data_3$dod_ym), merged_clock_data_3$DOD_minus_GSAPPT,  merged_clock_data_3$cutoff_minus_GSAPPT)

##  Create new DF trimmed  ###

merged_clock_data_4 <- subset(merged_clock_data_3, select = - c(ym, dod_ym, GS_yoa, GS_moa, GS_yod, GS_mod, cutoff_minus_GSAPPT, DOD_minus_GSAPPT))
# UP TO HERE SAVE AS SRIPT #

saveRDS(merged_clock_data_4, file="Merged_clock_data_13May2024.RDS")
