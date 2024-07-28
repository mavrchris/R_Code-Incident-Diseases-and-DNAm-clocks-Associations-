#RENV#


### RENV to ensure reproducibility ###
install.packages("renv")
library(renv)

# create lock file #
renv::init()

# Update locked file with any new packages used #
renv::snapshot()

# Download +install all used packages in lock file (with correct version) #
renv::restore()


#STEP 1#
### Set Directory ###

setwd ("Y:/Christos/RENV")


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

# STEP 2#
library(dplyr)
# Install and Load Packages # 
install.packages("survival", lib="Y:/Christos/GS data")
library (survival)


# Read in step_1_R RDS object <- Merged_clock_data_4

merged_clock_data_4 <- readRDS("Merged_clock_data_13May2024.RDS")

###  Read in Diseases  ###

diseases<- read.csv("2024-03-29_diseases.csv", stringsAsFactors = FALSE)

### Unique() , list of all individual diseases in disease data  ###

dis_list<- unique(diseases$Disease)

## check length/number of unique diseases ##

length(dis_list)
dis_list

## create list objects for cox models ## 

res.cox.Grim = list()
res.cox.Han = list()
res.cox.Pheno = list()
res.cox.Dunedin = list()

# list for F prop, Age % and event #

female_proportions = list()
age_percentiles = list()
n_event = list()
median_age_per_disease <- list()

# list for PROP. Hz. Asuumption test #
res.cox.zph.Han <- list()
res.cox.zph.Pheno <- list()
res.cox.zph.Grim <- list()
res.cox.zph.Dunedin <- list()



### FOR- Loop intiation ###


for (x in 1:length(dis_list)){ 
  disease = dis_list[x]
  disease_data <- diseases %>% filter(Disease== disease)
  print(paste0('disease ', x, ' of ' , length(dis_list)))
  
  ### Calculate time to incident  : ###
  
  ###  First : new column GP diagnosis= DiseaSED  = time of incident  time year and month  ###
  disease_data$yod <- substr(disease_data$dt1_ym, 1,4) %>% as.numeric 
  disease_data$mod <- substr(disease_data$dt1_ym, 5,6) %>% as.numeric 
  
  ### then new column of time  of GS appt ### 
  disease_data$yoa <- substr(disease_data$gs_app,1,4) %>% as.numeric 
  disease_data$moa <- substr(disease_data$gs_app,5,6) %>% as.numeric 
  
  ### new column time to disease, same principle with substr, reminder, /12== to give result in years ###
  disease_data$t_disease <- (disease_data$yod - disease_data$yoa) + ((disease_data$mod - disease_data$moa)/12)
  
  ### Remove excess columns ###
  filt_disease_data <- subset(disease_data, select = - c(mod,yoa,moa))
  
  #3#
  
  ### merge filtered disease data  with clock data, utilizing anchor points ### 
  merged_data <- merge(filt_disease_data, merged_clock_data_4, by.x= "id", by.y= "id", all.y= TRUE)
  
  
  #4#
  
  ## Remove individuals with incident flag of 0 = had disease from before (prevalent)= not allowing for predictions, whilst keeping na ##
  
  #keep na # 
  merged_data_2 <- merged_data[merged_data$incident != 0 | is.na(merged_data$incident), ]
  
  #5#
  
  ## New column : incident disease = 1, prevalent before gs appt= 0, so t_disease not negative =>DISEASE PREDICTION (DP)## 
  merged_data_2$DP <- ifelse(merged_data_2$t_disease >0,1,0)
  
  ## 5.A 1980 filter, removal of wrong data, filter for year of disease onset <1980 whilst keeping NA ##
  merged_data_3 <- merged_data_2 %>% filter (is.na(yod) | yod> 1980)
  
  # 5.B # 
  
  ## neg id , removal of t_disease less than zero (negative) ## 
  id_neg <- merged_data_3$id[which(is.na(merged_data_3$t_disease) | merged_data_3$t_disease < 0)]
  
  
  ## ordered by date, therefore always takes first event ##
  merged_data_ord <- merged_data_3[(order(merged_data_3$dt1_ym)),]
  
  ## Acount for Duplicates, due to primary and secondary appt, => considers only first ## 
  
  merged_data_ord <- merged_data_ord[ !duplicated(na.omit(merged_data_ord$id)),]
  
  #6#
  
  
  ### New column, time to event, event is either death or cutoff for non-diseased (t_censor) or time to disease for diseased ###
  merged_data_ord$event <- ifelse(is.na(merged_data_ord$incident),0,1)
  merged_data_ord$t_event <- ifelse(merged_data_ord$event == 1, merged_data_ord$t_disease,merged_data_ord$t_censor) 
  
  #trim#
  final_data<-  subset(merged_data_ord, select = - c(yod,X.x,X.y,dt1_ym,gs_appt))
  
  
  
  
  
  ### Move to next disease if no events for this disease or else to avoid break of loop ###
  
  if(length(which(final_data$event==1)) == 0){
    next
  }
  
  # Sex filtering #
  
  prop_F <-  length(which(final_data$sex == "F" & final_data$event ==1))/ length(which(final_data$event ==1 ))
  if (prop_F > 0.9 ) { 
    
    final_data = final_data %>% filter(sex == "F") 
    
  }
  
  if (prop_F < 0.1) {
    final_data = final_data %>% filter(sex == "M") 
    
  }
  
  female_proportions[[disease]] = prop_F
  
  # AGE filter  # 
  
  # first new column age of onset #
  
  final_data$onset_age <- final_data$age + final_data$t_disease
  
  #Calculate percentiles in vector #
  
  onset_age_boundaries <- final_data %>%
    summarise(
      P2.5 = quantile(onset_age, 0.025, na.rm = TRUE),
      P97.5 = quantile(onset_age, 0.975, na.rm = TRUE)
      
      
      
      #### Create columns for upper and loweR and filter in cox , did not work whe using columns from vector ####
      
      final_data$Upper_onset<- onset_age_boundaries$P97.5
      final_data$Lower_onset<- onset_age_boundaries$P2.5
      
      age_percentiles[[disease]] = onset_age_boundaries
      
      
      
      # return a vector of ncases for each disease # 
      
      
      
      # Calculate the sum of events where event == 1 for the current disease
      n_event <- sum(final_data$event == 1)
      
      final_data$n_event <- n_event
      
      # FILTER ABOVE N EVENT 30 #
      
      if (any(final_data$n_event < 30)) {
        next
      }
      
      ### for appendix file 3: (disease:ncase:ncontrol:mean_age:SD_age:prof_F)  ###
      
      # appnedix table 3 data generation #
      mean_age <- mean(final_data$age, na.rm = TRUE)
      sd_age <- sd(final_data$age, na.rm = TRUE)
      ncase <- n_event
      ncontrol <- 17855 - ncase
      percent_female <- prop_F * 100
      
      # Append summary statistics to the list# 
      disease_summary[[disease]] <- data.frame(
        disease = disease,
        ncase = ncase,
        ncontrol = ncontrol,
        mean_age = mean_age,
        sd_age = sd_age,
        percent_female = percent_female
      )
      
      #7#
      
      # Run Cox Regression Model #
      
      # Construct 4 models for 4 clocks, event= disease/death/cut-off, therefore use t_event == including non-diseased as controls #
      # if/else accounts for some diseases having only F/M therefore not including sex in covariate #
      if(prop_F > 0.1 | prop_F < 0.9) {
        res.cox.Han[[disease]] <- coxph(Surv(t_event,event)~ scale(DNAmAgeHannumAdjAge) + age + bmi + years  + ever_smoke + pack_years + units  + rank , data = final_data[final_data$age + final_data$t_event>final_data$Lower_onset & final_data$age + final_data$t_event<final_data$Upper_onset,])
        
        res.cox.Pheno[[disease]] <- coxph(Surv(t_event,event)~ scale(DNAmPhenoAgeAdjAge) + age + bmi + years + ever_smoke + pack_years + units  + rank , data = final_data[final_data$age + final_data$t_event>final_data$Lower_onset & final_data$age + final_data$t_event<final_data$Upper_onset,])
        
        res.cox.Grim[[disease]] <- coxph(Surv(t_event,event)~ scale(DNAmGrimAgeAdjAge) + age + bmi + years   + ever_smoke + pack_years + units + rank , data = final_data[final_data$age + final_data$t_event>final_data$Lower_onset & final_data$age + final_data$t_event<final_data$Upper_onset,])
        
        res.cox.Dunedin[[disease]] <- coxph(Surv(t_event,event)~ scale(DunedinPACE) + age + bmi + years   + ever_smoke + pack_years + units  + rank , data = final_data[final_data$age + final_data$t_event>final_data$Lower_onset & final_data$age + final_data$t_event<final_data$Upper_onset,])
        
      } else{
        
        res.cox.Han[[disease]] <- coxph(Surv(t_event,event)~ scale(DNAmAgeHannumAdjAge) + age + bmi + years  + ever_smoke + pack_years + units  + rank + sex , data = final_data[final_data$age + final_data$t_event>final_data$Lower_onset & final_data$age + final_data$t_event<final_data$Upper_onset,])
        
        res.cox.Pheno[[disease]] <- coxph(Surv(t_event,event)~ scale(DNAmPhenoAgeAdjAge) + age + bmi + years + ever_smoke + pack_years + units  + rank + sex , data = final_data[final_data$age + final_data$t_event>final_data$Lower_onset & final_data$age + final_data$t_event<final_data$Upper_onset,])
        
        res.cox.Grim[[disease]] <- coxph(Surv(t_event,event)~ scale(DNAmGrimAgeAdjAge) + age + bmi + years   + ever_smoke + pack_years + units + rank + sex , data = final_data[final_data$age + final_data$t_event>final_data$Lower_onset & final_data$age + final_data$t_event<final_data$Upper_onset,])
        
        res.cox.Dunedin[[disease]] <- coxph(Surv(t_event,event)~ scale(DunedinPACE) + age + bmi + years   + ever_smoke + pack_years + units  + rank + sex , data = final_data[final_data$age + final_data$t_event>final_data$Lower_onset & final_data$age + final_data$t_event<final_data$Upper_onset,])
      }
      
      
      #Proportional Hazards Assumption test#
      #HZR#
      print(paste0('running zph for ', disease))
      res.cox.zph.Han[[disease]] = cox.zph(res.cox.Han[[disease]])$table[1,]
      res.cox.zph.Pheno[[disease]] = cox.zph(res.cox.Pheno[[disease]])$table[1,]
      res.cox.zph.Grim[[disease]] = cox.zph(res.cox.Grim[[disease]])$table[1,]
      res.cox.zph.Dunedin[[disease]] = cox.zph(res.cox.Dunedin[[disease]])$table[1,]
      
      
      #less load/data kept#
      rm(final_data)
      gc() 
} # close for-loop #

#HZR: tabulate them essentially # 
res.cox.zph.all.Dunedin = do.call("rbind", res.cox.zph.Dunedin)
res.cox.zph.all.Grim = do.call("rbind", res.cox.zph.Grim)
res.cox.zph.all.Pheno = do.call("rbind", res.cox.zph.Pheno)
res.cox.zph.all.Han = do.call("rbind", res.cox.zph.Han)


# AS.DATA.FRAME <- converts to allow filter #
res.cox.zph.all.Dunedin <- as.data.frame(res.cox.zph.all.Dunedin)
res.cox.zph.all.Grim <- as.data.frame(res.cox.zph.all.Grim)
res.cox.zph.all.Pheno<- as.data.frame(res.cox.zph.all.Pheno)
res.cox.zph.all.Han<- as.data.frame(res.cox.zph.all.Han)

# see significant ones =  null assumption is not met #
significant_results_Dun_all <- res.cox.zph.all.Dunedin[res.cox.zph.all.Dunedin$p < 0.05/299, ]
significant_results_Grim_all <- res.cox.zph.all.Grim[res.cox.zph.all.Grim$p < 0.05/299, ]
significant_results_Pheno_all <- res.cox.zph.all.Pheno[res.cox.zph.all.Pheno$p < 0.05/299, ]
significant_results_Han_all <- res.cox.zph.all.Han[res.cox.zph.all.Han$p < 0.05/299, ]


# EXPORT AS CSV- appendix file #
write.csv(res.cox.zph.all.Dunedin, "res.cox.zph.all.Dunedin.csv")
write.csv(res.cox.zph.all.Grim, "res.cox.zph.all.Grim.csv")
write.csv(res.cox.zph.all.Pheno, "res.cox.zph.all.Pheno.csv")
write.csv(res.cox.zph.all.Han, "res.cox.zph.all.Han.csv")


# SAVE DF AS RDS #
saveRDS(res.cox.Han, file = "res.cox.Han_10_JUNE_AGE_SEX_2024.RDS")
saveRDS(res.cox.Pheno,  file = "res.cox.Pheno_10_JUNE_AGE_SEX_2024.RDS")
saveRDS(res.cox.Grim,  file = "res.cox.Grim_10_JUNE_AGE_SEX_2024.RDS")
saveRDS(res.cox.Dunedin,  file = "res.cox.Dunedin_10_JUNE_AGE_SEX_2024.RDS")


# NEW APPENDIX TABLE 3##
#Define the file path to save the CSV file#
file_path_2 <- "Y:/Christos/MScR II/Final Results/appendix_table_3.csv"  

# tabulate essentially # 
disease_summary = do.call("rbind", disease_summary)

# Write the dataframe to a CSV file
write.csv(disease_summary, file = file_path_2, row.names = FALSE)



# STEP 3 #
install.packages("dplyr")
library(dplyr)
install.packages("ggplot2")
library(ggplot2)
install.packages("survival")
library(survival)
install.packages("patchwork")
library(patchwork)
install.packages("cowplot")
library(cowplot)


# First READ IN COX MODELS # 

res.cox.Han <- readRDS("res.cox.Han_10_JUNE_AGE_SEX_2024.RDS")
res.cox.Grim <- readRDS("res.cox.Grim_10_JUNE_AGE_SEX_2024.RDS")
res.cox.Pheno <- readRDS("res.cox.Pheno_10_JUNE_AGE_SEX_2024.RDS")
res.cox.Dunedin <- readRDS("res.cox.Dunedin_10_JUNE_AGE_SEX_2024.RDS")



# Next, create a data frame containing cox model outputs # 


### required columns: disease, Ncase, N, clock, HazardRatio, Lower CI, Upper CI, P-value ### 

#1. GrimAge #

## return a vector of N for each model ## sapply, (s-> vector), so streamlining, remove , and make sapply #

ntotal_Grim = sapply(res.cox.Grim, function(x){summary(x)$n})


# return a vector of ncases for each disease # 

ncases_Grim = sapply(res.cox.Grim, function(x){summary(x)$nevent})

# return a vector of hazard ratios # 

hr_Grim = sapply(res.cox.Grim, function(x){summary(x)$coefficient["scale(DNAmGrimAgeAdjAge)",2]})
pval_Grim = sapply(res.cox.Grim, function(x){summary(x)$coefficient["scale(DNAmGrimAgeAdjAge)",5]})

lower_ci_Grim = sapply(res.cox.Grim, function(x){exp(confint(x))["scale(DNAmGrimAgeAdjAge)",1]})
upper_ci_Grim = sapply(res.cox.Grim, function(x){exp(confint(x))["scale(DNAmGrimAgeAdjAge)",2]})

## Plotting DF  ##

plot_df_grim = data.frame(dis = names(res.cox.Grim),
                          ncase = ncases_Grim,
                          n = ntotal_Grim,
                          clock = 'GrimAge',
                          hr = hr_Grim,
                          lci = lower_ci_Grim,
                          uci = upper_ci_Grim, 
                          p = pval_Grim)


#2. Han #

# return a vector of N for each model #

ntotal_Han = unlist(lapply(res.cox.Han, function(x){summary(x)$n}))

# return a vector of ncases for each disease # 

ncases_Han = unlist(lapply(res.cox.Han, function(x){summary(x)$nevent}))

# return a vector of hazard ratios # 

hr_Han = sapply(res.cox.Han, function(x){summary(x)$coefficient["scale(DNAmAgeHannumAdjAge)",2]})
pval_Han = sapply(res.cox.Han, function(x){summary(x)$coefficient["scale(DNAmAgeHannumAdjAge)",5]})

lower_ci_Han = sapply(res.cox.Han, function(x){exp(confint(x))["scale(DNAmAgeHannumAdjAge)",1]})
upper_ci_Han = sapply(res.cox.Han, function(x){exp(confint(x))["scale(DNAmAgeHannumAdjAge)",2]})

## Plotting DF  ##

plot_df_han = data.frame(dis = names(res.cox.Han),
                         ncase = ncases_Han,
                         n = ntotal_Han,
                         clock = 'Hannum',
                         hr = hr_Han,
                         lci = lower_ci_Han,
                         uci = upper_ci_Han, 
                         p = pval_Han)

# Dunedin #

# return a vector of N for each model #

ntotal_Dun = sapply(res.cox.Dunedin, function(x){summary(x)$n})

# return a vector of ncases for each disease # 

ncases_Dun = sapply(res.cox.Dunedin, function(x){summary(x)$nevent})

# return a vector of hazard ratios # 

hr_Dun = sapply(res.cox.Dunedin, function(x){summary(x)$coefficient["scale(DunedinPACE)",2]})
pval_Dun = sapply(res.cox.Dunedin, function(x){summary(x)$coefficient["scale(DunedinPACE)",5]})

lower_ci_Dun = sapply(res.cox.Dunedin, function(x){exp(confint(x))["scale(DunedinPACE)",1]})
upper_ci_Dun= sapply(res.cox.Dunedin, function(x){exp(confint(x))["scale(DunedinPACE)",2]})

## Plotting DF  ##

plot_df_dun = data.frame(dis = names(res.cox.Dunedin),
                         ncase = ncases_Dun,
                         n = ntotal_Dun,
                         clock = 'Dunedin',
                         hr = hr_Dun,
                         lci = lower_ci_Dun,
                         uci = upper_ci_Dun, 
                         p = pval_Dun)

# Pheno #

# return a vector of N for each model #

ntotal_Pheno = sapply(res.cox.Pheno, function(x){summary(x)$n})

# return a vector of ncases for each disease # 

ncases_Pheno = sapply(res.cox.Pheno, function(x){summary(x)$nevent})

# return a vector of hazard ratios # 

hr_Pheno = sapply(res.cox.Pheno, function(x){summary(x)$coefficient["scale(DNAmPhenoAgeAdjAge)",2]})
pval_Pheno = sapply(res.cox.Pheno, function(x){summary(x)$coefficient["scale(DNAmPhenoAgeAdjAge)",5]})

lower_ci_Pheno = sapply(res.cox.Pheno, function(x){exp(confint(x))["scale(DNAmPhenoAgeAdjAge)",1]})
upper_ci_Pheno = sapply(res.cox.Pheno, function(x){exp(confint(x))["scale(DNAmPhenoAgeAdjAge)",2]})

## Plotting DF  ##

plot_df_pheno = data.frame(dis = names(res.cox.Pheno),
                           ncase = ncases_Pheno,
                           n = ntotal_Pheno,
                           clock = 'PhenoAge',
                           hr = hr_Pheno,
                           lci = lower_ci_Pheno,
                           uci = upper_ci_Pheno, 
                           p = pval_Pheno)



### Filter Before Plots ### 

# Filter 1 already filtered, sanity check: ncount > 30 #

plot_df_dun_1 <- plot_df_dun %>% filter (ncase > 30)

plot_df_grim_1 <- plot_df_grim %>% filter (ncase > 30)

plot_df_han_1 <- plot_df_han %>% filter (ncase > 30)

plot_df_pheno_1<- plot_df_pheno %>% filter (ncase > 30)


# Filter 2 : P value/ Bonferroni Corr., P/ ndisease x nclocks # Precise -> 299 #



plot_df_dun_P_filt <- plot_df_dun_1 %>% filter (p < 0.05/299)

plot_df_grim_P_filt<- plot_df_grim_1 %>% filter (p < 0.05/299)

plot_df_han_P_filt  <- plot_df_han_1 %>% filter (p < 0.05/299)

plot_df_pheno_P_filt <- plot_df_pheno_1 %>% filter (p < 0.05/299)

### For graphs with : 1 disease and all 4 clocks  ###

plot_df_all <- rbind(plot_df_dun_1,plot_df_grim_1, plot_df_han_1, plot_df_pheno_1)
# as 4 clocks, bonferoni correction filter= 4 x299 = 1196 #
sig_diseases = unique(plot_df_all$dis[plot_df_all$p < 0.05/1196])
plot_df_all_2 = plot_df_all %>% filter(dis%in% sig_diseases)

### Unique() , list of all individual diseases in disease data  ###
diseases<- read.csv("2024-03-29_diseases.csv")

dis_list_unique<- unique(diseases$Disease)

## check length/number of unique diseases ##

length(dis_list_unique)







### Order By Highest HR ###

plot_df_dun_3 <- plot_df_dun_P_filt %>%
  mutate(dis = reorder(dis, hr))

plot_df_grim_3 <- plot_df_grim_P_filt %>%
  mutate(dis = reorder(dis, hr))

plot_df_han_3 <- plot_df_han_P_filt %>%
  mutate(dis = reorder(dis, hr))

plot_df_pheno_3 <- plot_df_pheno_P_filt %>%
  mutate(dis = reorder(dis, hr))

### Plotting ###
# total : #
plot_df <- rbind(plot_df_grim_P_filt,plot_df_dun_P_filt, plot_df_han_P_filt, plot_df_pheno_P_filt)


colnames(plot_df_all_2)[colnames(plot_df_all_2) == "clock"] <- "Clock"

plot_df_all_2$Significant <- ifelse(plot_df_all_2$p < 0.05/1196, 'Significant', 'Non-significant')

plot_df_all_2

# Define the file path to save the CSV file
file_path <- "Y:/Christos/MScR II/Final Results/plot_df_all_2.csv"  # Change this to the desired file path


# Write the dataframe to a CSV file
write.csv(plot_df_all_2, file = file_path, row.names = FALSE)

#for shapes of HR#

plot_df_all_3 <- plot_df_all_2 %>%
  mutate(shape = case_when(
    p > 0.05 ~ "non-significant",
    p < 0.05 & p >= 0.05/299 ~ "significant",
    p < 0.05/299 ~ "Bonf.-significant"
  ))

#incorrect data#
disease_to_remove <- "covid"

# Filter out the rows with the specified disease
plot_df_all_4 <- plot_df_all_3 %>%
  filter(dis != disease_to_remove)


####################
library(ggplot2)
library(dplyr)
library(cowplot)

# Plot creation, 1 pdf, multiple plots within #

# Create a list to store plots
plot_list <- list()

# Create plots and store them in the list #
for (i in unique(as.character(plot_df_all_4$dis))) {
  # Subset data for the current disease #
  tmp_plot <- plot_df_all_4 %>% filter(dis == i)
  
  # Extract ncase for the current disease from tmp_plot #
  ncase <- unique(tmp_plot$ncase)
  
  # Create the forrest plot #
  plot_out <- ggplot(tmp_plot, aes(x = reorder(Clock, +hr), y = as.numeric(hr), col = Clock, fill = Clock)) +
    geom_point(aes(shape = shape), size = 2, stroke = 0.3, position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = as.numeric(lci), ymax = as.numeric(uci)), width = 0.1, position = position_dodge(width = 0.8)) +
    ylab('Hazard Ratio (95% CI)') +
    # Remove X-axis title #
    xlab('') +  
    geom_hline(yintercept = 1) +
    coord_flip() +
    # Add disease name and n_case to the title #
    ggtitle(paste(i, ",  n_case:", ncase)) +  
    scale_shape_manual(values = c("non-significant" = 16, "significant" = 22, "Bonf.-significant" = 24),
                       labels = c("non-significant" = expression("P-Non-Significant  ">=" 0.05"), "significant" = "P-Nominal < 0.05", "Bonf.-significant" = expression("P-Bonferroni < " ~ 1.67 %*% 10^-24))) +
    theme(
      plot.title = element_text(hjust = 0.08, size= 11),
      axis.text.y = element_blank(),   # Remove y-axis text labels (which are x-axis labels after CORD-FLIP) #
      axis.title.y = element_blank(),  # Remove y-axis title (which is x-axis title after CORD-FLIP) #
      axis.ticks.y = element_blank()   # Remove y-axis ticks (which are x-axis ticks after CORD-FLIP) #
    )
  
  # Add plot to the list
  plot_list[[i]] <- plot_out
}

# Choose the number of plots per page  #
plots_per_page <- 9  

# Calculate the number of pages needed #
num_pages <- ceiling(length(plot_list) / plots_per_page)

# Create a PDF file with multiple pages #


pdf("Combined_Plots_JULLY_28.pdf", width = 11, height = 8.5)  

# Extract plots for the current page #
for (page in 1:num_pages) {
  
  start_idx <- (page - 1) * plots_per_page + 1
  end_idx <- min(page * plots_per_page, length(plot_list))
  plots_for_page <- plot_list[start_idx:end_idx]
  
  # Combine plots into a grid , 3 columns per page #
  combined_plot <- plot_grid(plotlist = plots_for_page, ncol = 3)  
  
  # Print the combined plot to the PDF #
  print(combined_plot)
}

# Close the PDF #
print (dev.off())

