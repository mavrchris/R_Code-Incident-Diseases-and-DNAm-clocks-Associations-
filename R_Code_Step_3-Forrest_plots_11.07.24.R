# Step # : data Visualization #

#RENV to ensure reproducibility #
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
#setwd("Y:/marioni-lab/Christos/GS data")

# MY PC# 

setwd ("Y:/Christos/GS data") 

###  Load Packages ###
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
                       labels = c("non-significant" = "Non-significant", "significant" = "Significant", "Bonf.-significant" = "Bonf.-significant")) +
    theme(
      plot.title = element_text(hjust = 0.08),
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

pdf("Combined_Plots_JULY_08_third.pdf", width = 11, height = 8.5)  

# Extract plots for the current page #
for (page in 1:num_pages) {
  
  start_idx <- (page - 1) * plots_per_page + 1
  end_idx <- min(page * plots_per_page, length(plot_list))
  plots_for_page <- plot_list[start_idx:end_idx]
  
  # Combine plots into a grid , 3 columns per page #
  combined_plot <- plot_grid(plotlist = plots_for_page, ncol = 3)  
  
  # Print the combined plot to the PDF #
  print(combined_plot)
  
  # Close the PDF #
  print (dev.off())
}


