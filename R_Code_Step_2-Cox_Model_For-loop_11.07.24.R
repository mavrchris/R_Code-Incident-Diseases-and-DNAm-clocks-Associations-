# Step 2 #

#RENV to ensure reproducibility #
install.packages("renv")
library(renv)
# create lock file #
renv::init()
# Update locked file with any new packages used #
renv::snapshot()
# Download +install all used packages in lock file (with correct version) #
renv::restore()

## Data prep, before loop ##

### Set Directory ###
# Uni- PC #
#setwd("Y:/marioni-lab/Christos/GS data")

# MY PC# 

setwd ("Y:/Christos/GS data") 

###  Load Packages ###
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



