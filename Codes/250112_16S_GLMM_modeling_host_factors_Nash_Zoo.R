#Researcher: Chloe Cummins
#Date: January 12, 2025
#GLMM modeling of host-associated factors at the Nashville Zoo (and maybe Chatt. Zoo)


#Set working directory
setwd("G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Hellbender Gut Microbiome Manuscript Stuff/Manuscript_Draft_Code_and_Figs")

#Previous code originally stored in "G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Hellbender Gut Microbiome Manuscript Stuff/Manuscript_Code_Check"
#and in "G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Database/Data/Thesis_16S_Analyses"

#Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(ggpubr)
library(RColorBrewer)
library(tibble)
library(data.table)
library(plotly)
library(phyloseq)
library(forcats)
library(hillR)
library(car)
library(DHARMa)
library(performance)
library(lmodel2)  #for SMI calculation
library(glmmTMB)
library(emmeans)
library(buildmer)
library(ggh4x)
library(DspikeIn)
library(MuMIn)
library(AICcmodavg)


#Importing the decontam and rarefied data
otu_count <- fread("240820_HB_16S_OTU_count_table_decontam.csv")
tax_table <- fread("240820_HB_16S_taxonomy_decontam.csv")
metadata <- fread("240820_HB_16S_metadata_decontam.csv")

##### Preparing the data for GLMM modeling #####

#Converting all of the data into data frames
otu_table_2 <- as.data.frame(otu_count)

tax_table_2 <- as.data.frame(tax_table)

metadata_2 <- as.data.frame(metadata)

#Formatting OTU table for Hill number calculations
otu_table_2 <- otu_table_2 %>% 
  column_to_rownames("index")

##### Calculating Hill numbers for alpha diversity and adding to the metadata file #####

metadata_2$h0 <- hill_taxa(otu_table_2, q = 0, MARGIN = 1)
#q = 0 is essentially a richness metric, placing more emphasis on rare taxa

#metadata_2$h1 <- hill_taxa(otu_table_2, q = 1, MARGIN = 1)
#q = 1 is Hill-Shannon diversity, where abundant and rare taxa are equally weighted
#and evenness is accounted for

#metadata_2$h2 <- hill_taxa(otu_table_2, q = 2, MARGIN = 1)
#q = 2 is inverse Simpson, which places more emphasis on common/abundant species

#metadata_2$h3 <- hill_taxa(otu_table_2, q = 3, MARGIN = 1)
#q = 3 is where dominant, abundant species are given the most weight



##### Creating Nashville only data frame #####

#Filtering for only cloacal swab samples from the Nashville Zoo
meta_Nash <- filter(metadata_2, site == "Nashville Zoo" & env_medium == "cloaca swab")

#Making sure sample count is the same after filtering out samples
meta_otu_inter <- intersect(meta_Nash$index, otu_count$index)
otu_table_Nash <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
otu_table_Nash <- otu_table_Nash[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
otu_table_Nash <- otu_table_Nash %>% 
  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_Nash))
tax_table_Nash <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_Nash <- otu_table_Nash %>% select(all_of(intersected))



##### Creating Chattanooga only data frame #####

#meta_Chatt <- filter(metadata_2, site == "Chattanooga Zoo")

#Making sure sample count is the same after filtering out samples
#meta_otu_inter <- intersect(meta_Chatt$index, otu_count$index)
#otu_table_Chatt <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
#otu_table_Chatt <- otu_table_Chatt[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_Chatt <- otu_table_Chatt %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected <- intersect(tax_table_2$otu, colnames(otu_table_Chatt))
#tax_table_Chatt <- tax_table_2 %>% filter(otu %in% intersected)
#otu_table_Chatt <- otu_table_Chatt %>% select(all_of(intersected))



##### Creating Nashville and Chattanooga data frame #####

#meta_zoo <- filter(metadata_2, site == "Nashville Zoo" & env_medium == "cloaca swab" |
#                       site == "Chattanooga Zoo")

#Making sure sample count is the same after filtering out samples
#meta_otu_inter <- intersect(meta_zoo$index, otu_count$index)
#otu_table_zoo <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
#otu_table_zoo <- otu_table_zoo[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_zoo <- otu_table_zoo %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected <- intersect(tax_table_2$otu, colnames(otu_table_zoo))
#tax_table_zoo <- tax_table_2 %>% filter(otu %in% intersected)
#otu_table_zoo <- otu_table_zoo %>% select(all_of(intersected))




#####Preparing Nashville Zoo data for modeling#####

#Adding a column in the metadata to calculate age of animal at collection (in years)
meta_Nash$hatch_date <- as.Date(meta_Nash$hatch_date, "%m/%d/%Y")

meta_Nash$collection_date <- as.Date(meta_Nash$collection_date, "%m/%d/%Y")

meta_Nash$age_at_collection <- difftime(meta_Nash$collection_date, 
                                        meta_Nash$hatch_date, units = "days") / 365
#converting the days result into years because difftime does not offer years as an argument

meta_Nash <- mutate(meta_Nash,
                    age_at_collection = as.integer(age_at_collection), 
                    age_at_collection = str_replace_all(age_at_collection, "days", ""), 
                    #age_at_collection = as.factor(age_at_collection)
                    age_at_collection = as.numeric(age_at_collection)
)
#converting the decimals to whole numbers and removing the word "days" from the column

#Coding sex as a factor
meta_Nash$sex <- as.factor(meta_Nash$sex)



##### Scaled Mass Index (SMI) Calculation - Nashville Zoo only #####

#Natural log transforming the mass and TL data in order to perform standardized major axis (SMA) regression
#per the steps outlined in Peig and Green 2009
meta_Nash$ln_mass <- log(meta_Nash$weight)

meta_Nash$ln_TL <- log(meta_Nash$total_length)

#Performing a model II SMA regression
sma <- lmodel2(ln_mass ~ ln_TL, data = meta_Nash, nperm = 99)
#SMA slope =  2.840927 - accessed using sma[["regression.results"]][["Slope"]] with the third number listed
#this slope is the scaling exponent for the SMI calculation

#Calculating the scaling exponent for the SMI equation
bsma <- sma[["regression.results"]][["Slope"]][[3]]
bsma

#Calculating the mean total length for the SMI equation
l0 <- mean(meta_Nash$total_length, na.rm = TRUE)
l0
#mean TL = 31.96296 cm

#Calculating SMI based on the equation
meta_Nash$SMI <- (meta_Nash$weight * (l0/meta_Nash$total_length) ^ bsma)

#meta_Nash$scaled_SMI <- scale(meta_Nash$SMI)



##### Scaled Mass Index (SMI) Calculation - Chattanooga Zoo only #####

#Natural log transforming the mass and TL data in order to perform standardized major axis (SMA) regression
#per the steps outlined in Fig. 1 of Peig and Green 2009
#meta_Chatt$ln_mass <- log(meta_Chatt$weight)

#meta_Chatt$ln_TL <- log(meta_Chatt$total_length)

#Performing a model II SMA regression
#sma <- lmodel2(ln_mass ~ ln_TL, data = meta_Chatt, nperm = 99)
#SMA slope =  3.863747
#this slope is the scaling exponent for the SMI calculation

#bsma <- sma[["regression.results"]][["Slope"]][[3]]
#bsma

#l0 <- mean(meta_Chatt$total_length, na.rm = TRUE)
#l0
#mean TL = 32.71667 cm

#meta_Chatt$SMI <- (meta_Chatt$weight * (l0/meta_Chatt$total_length) ^ bsma)

#meta_Chatt$scaled_SMI <- scale(meta_Chatt$SMI)



##### Scaled Mass Index (SMI) Calculation - Nashville and Chattanooga Zoo #####

#Natural log transforming the mass and TL data in order to perform standardized major axis (SMA) regression
#per the steps outlined in Fig. 1 of Peig and Green 2009
#meta_zoo$ln_mass <- log(meta_zoo$weight)

#meta_zoo$ln_TL <- log(meta_zoo$total_length)

#Performing a model II SMA regression
#sma <- lmodel2(ln_mass ~ ln_TL, data = meta_zoo, nperm = 99)
#SMA slope =  3.116265
#this slope is the scaling exponent for the SMI calculation

#bsma <- sma[["regression.results"]][["Slope"]][[3]]
#bsma

#l0 <- mean(meta_zoo$total_length, na.rm = TRUE)
#l0
#mean TL = 32.02449 cm

#meta_zoo$SMI <- (meta_zoo$weight * (l0/meta_zoo$total_length) ^ bsma)

#meta_zoo$scaled_SMI <- scale(meta_zoo$SMI)



#1. Evaluating whether richness is affected by host-associated factors at sample collection - Nashville Zoo only----

#As a note, this model may be dropped or may change based on the conversation that we had during my defense
#about the SMI results. However, if the model seems ok after double checking and we can feel confident about the 
#results, then it may be included in the manuscript.


#Assessing normality of richness (response variable) at Nash Zoo
shapiro.test(meta_Nash$h0) #p-value < 2.2e-16, not normally distributed
#Data is not normally distributed, indicating that a GLMM may produce a better fit than an LMM

#Looking at the distribution of the response variable to double check that the stats. match
densityplot(meta_Nash$h0, main = "Distribution of Richness")
#right skew

#Histogram of richness
hist(meta_Nash$h0)

#Q-Q plot of richness
qqnorm(meta_Nash$h0)
qqline(meta_Nash$h0, col = "red")
#there seem to be some points near the top of the graph

#Q-Q plot to see which age group is causing the points near the top in the previous plot
h0_qq_age <- ggplot(meta_Nash, aes(sample = h0, color = age_at_collection)) + 
  stat_qq() +
  stat_qq_line() +
  theme_classic()
print(h0_qq_age)
#points seem to be across all age groups

#Q-Q plot to see which sex is causing the points near the top in the previous plot
h0_qq_sex <- ggplot(meta_Nash, aes(sample = h0, color = sex)) + 
  stat_qq() +
  stat_qq_line() +
  theme_classic()
print(h0_qq_sex)
#points seem to be across both sexes


## GLMM of host factors on richness ##

#Evaluating whether I should use a tweedie distribution or a poisson distribution
#h0_1 <- glmmTMB(h0 ~ age_at_collection + (1 | sample_name), data = meta_Nash, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
h0_2 <- glmmTMB(h0 ~ age_at_collection + (1 | sample_name), data = meta_Nash, family = poisson())

AIC(h0_2)
AICc(h0_2)
summary(h0_2)
#The inclusion of tweedie is mainly an artifact in the code since tweedie is mainly used for continuous, zero-inflated data
#so it may not be appropriate in this case and was kept to maintain the numbering of the remaining models


#Evaluating the addition of more fixed effects that are host associated
h0_3 <- glmmTMB(h0 ~ age_at_collection + sex + (1 | sample_name), data = meta_Nash, family = poisson())

h0_4 <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 | sample_name), data = meta_Nash, family = poisson())

AIC(h0_2, h0_3, h0_4)
AICc(h0_2)
AICc(h0_3)
AICc(h0_4)
#h0_4 has the lowest AIC and AICc

#Evaluating the addition of more random effects
h0_5 <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 | sample_name) + (1 | dna_experiment_date), 
                data = meta_Nash, family = poisson())

h0_6 <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 | sample_name) + (1 | collection_date), 
                data = meta_Nash, family = poisson())

h0_7 <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 | sample_name) + (1 | plate_ID), 
                        data = meta_Nash, family = poisson())

h0_8 <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 | sample_name) + (1 | well_location), 
                data = meta_Nash, family = poisson())

h0_9 <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 | sample_name) + (1 | dna_experiment_date) + 
                     (1 | collection_date) + (1 | plate_ID), data = meta_Nash, family = poisson())

h0_10 <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 | sample_name) + (1 | dna_experiment_date) + 
                  (1 | collection_date) + (1 | plate_ID) + (1 | well_location), data = meta_Nash, family = poisson())

#All models returned a warning message, but from what I understand based on the glmmTMB documentation, 
#this warning may be ok as long as it is not accompanied by a convergence error/warning. 

#Link to documentation (https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html)
#"This warning occurs when the optimizer visits a region of parameter space that is invalid. 
#It is not a problem as long as the optimizer has left that region of parameter space upon convergence, 
#which is indicated by an absence of the model convergence warnings described above."

AIC(h0_4, h0_5, h0_6, h0_7, h0_8, h0_9, h0_10)
AICc(h0_4)
AICc(h0_5)
AICc(h0_6)
AICc(h0_7)
AICc(h0_8)
AICc(h0_9)
AICc(h0_10)
#h0_10 has the lowest AIC and AICc

#Testing overdispersion (ratio greater than 1)
overdisp_ratio = deviance(h0_10)/df.residual(h0_10)
overdisp_ratio
#no overdispersion based on ratio

summary(h0_10)
#sex and age not significant but SMI is significant

group_means <- emmeans(h0_10, ~ age_at_collection + sex + SMI)
group_means

plot(residuals(h0_10))
#residuals look decent

sim_res <- simulateResiduals(h0_10)
plot(sim_res)
#dispersion and outlier tests are not significant but KS is significant

qqnorm(residuals(h0_10))
qqline(residuals(h0_10))
#qqplot looks ok - there is some deviation in both tails

Anova(h0_10, type = 2)
#SMI/body condition has a significant effect on richness but age and sex do not

r.squaredGLMM(h0_10)

#Adding a random slope to the model to see if that provides better fit after changing age at collection to numeric instead of factor
h0_10a <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 + age_at_collection | sample_name) + (1 | dna_experiment_date) + 
                    (1 | collection_date) + (1 | plate_ID) + (1 | well_location), data = meta_Nash, family = poisson())

AIC(h0_10, h0_10a)
AICc(h0_10)
AICc(h0_10a)
#model with the random slope has lower AIC and AICc

#Testing overdispersion (ratio greater than 1)
overdisp_ratio = deviance(h0_10a)/df.residual(h0_10a)
overdisp_ratio
#no overdispersion based on ratio

summary(h0_10a)
#sex and age not significant but SMI is slightly significant

group_means <- emmeans(h0_10a, ~ age_at_collection + sex + SMI)
group_means

plot(residuals(h0_10a))
#residuals look decent

sim_res <- simulateResiduals(h0_10a)
plot(sim_res)
#dispersion and outlier tests are not significant but KS is significant

qqnorm(residuals(h0_10a))
qqline(residuals(h0_10a))
#qqplot looks ok - there is some deviation in both tails

Anova(h0_10a, type = 2)
#SMI/body condition has a significant effect on richness but age and sex do not

r.squaredGLMM(h0_10a)


#Testing to see if interactions between the fixed effects provide a better fit to the model
#h0_11 <- glmmTMB(h0 ~ age_at_collection * sex * SMI + (1 | sample_name) + (1 | dna_experiment_date) + 
#                   (1 | collection_date) + (1 | plate_ID) + (1 | well_location), data = meta_Nash, family = poisson())
#Since there are no females in the age 6 group, those variables got dropped and no values were calculated for those comparisons
#I am unsure if this model should be used since not all of the comparisons from the interactions can be calculated

#AIC(h0_11)

#summary(h0_11)

#group_means <- emmeans(h0_11, ~ age_at_collection * sex * SMI)
#group_means

#plot(residuals(h0_11))

#qqnorm(residuals(h0_10))
#qqline(residuals(h0_10))

#sim_res <- simulateResiduals(h0_11)
#plot(sim_res)

#Anova(h0_11, type = 2)


#Conclusion: proceed with h0_10 because it is more interpretable and appears to have a fairly decent fit based on DHARMa and residuals
#and because including the interaction did not obviously/significantly improve fit and meant that some comparisons were not calculated


#testOutliers(h0_10)
#testOutliers(sim_res)
#check_outliers(sim_res, type = "default")
#check_outliers(sim_res, type = "binomial")
#check_outliers(sim_res, type = "bootstrap")



#2. Evaluating whether richness is affected by body condition - Chattanooga Zoo only----

#Assessing normality of richness (response variable) at Nash Zoo
#shapiro.test(meta_Chatt$h0) #p-value < 0.04513, not normally distributed
#Data is not normally distributed, indicating that a GLMM may produce a better fit than an LMM

#Looking at the distribution of the response variable to double check that the stats. match
#densityplot(meta_Chatt$h0, main = "Distribution of Richness")
#right skew

#Histogram of richness
#hist(meta_Chatt$h0)

#Q-Q plot of richness
#qqnorm(meta_Chatt$h0)
#qqline(meta_Chatt$h0, col = "red")



## GLMM of SMI on richness ##

#Evaluating whether I should use a tweedie distribution or a poisson distribution
#h0_1 <- glmmTMB(h0 ~ SMI + (1 | sample_name), data = meta_Chatt, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
#h0_2 <- glmmTMB(h0 ~ SMI + (1 | sample_name), data = meta_Chatt, family = poisson())

#AIC(h0_1, h0_2)
#summary(h0_2)
#Poisson had lower AIC 

#Evaluating the addition of more random effects
#h0_3 <- glmmTMB(h0 ~ SMI + (1 | sample_name) + (1 | dna_experiment_date), 
#                data = meta_Chatt, family = poisson())

#h0_4 <- glmmTMB(h0 ~ SMI + (1 | sample_name) + (1 | collection_date), 
#                data = meta_Chatt, family = poisson())

#h0_5 <- glmmTMB(h0 ~ SMI + (1 | sample_name) + (1 | plate_ID), 
#                data = meta_Chatt, family = poisson())

#h0_6 <- glmmTMB(h0 ~ SMI + (1 | sample_name) + (1 | well_location), 
#                data = meta_Chatt, family = poisson())

#h0_7 <- glmmTMB(h0 ~ SMI + (1 | sample_name) + (1 | dna_experiment_date) + 
#                  (1 | collection_date) + (1 | plate_ID), data = meta_Chatt, family = poisson())

#h0_8 <- glmmTMB(h0 ~ SMI + (1 | sample_name) + (1 | dna_experiment_date) + 
#                   (1 | collection_date) + (1 | plate_ID) + (1 | well_location), data = meta_Chatt, family = poisson())

#AIC(h0_2, h0_3, h0_4, h0_5, h0_6, h0_7, h0_8)
#Addition of random effects aside from sample name did not improve model fit
#Proceed with h0_2 since it is a simpler model and has lower AIC

#overdisp_ratio = deviance(h0_2)/df.residual(h0_2)
#overdisp_ratio

#summary(h0_2)

#group_means <- emmeans(h0_2, ~ SMI)
#group_means

#plot(residuals(h0_2))

#sim_res <- simulateResiduals(h0_2)
#plot(sim_res)

#Anova(h0_2, type = 2)
#SMI does not have a sig. effect on richness when Chatt is analyzed by itself



##### Random GLMM test with Nash and Chatt Zoo samples analyzed together for SMI effect on richness #####

#h0_test <- glmmTMB(h0 ~ SMI + (1 | sample_name) + (1 | dna_experiment_date) + 
#                   (1 | collection_date) + (1 | plate_ID) + (1 | well_location), data = meta_zoo, family = poisson())

#summary(h0_test)
#plot(residuals(h0_test))

#sim_res <- simulateResiduals(h0_test)
#plot(sim_res)

#Anova(h0_test, type = 2)
#SMI has a sig. effect on richness with Nash and Chatt are analyzed together



##### GLMM modelling without individuals that have high richness #####

#Filter out individuals that have high richness
#high_rich <- c(11973, 12120, 12187, 12105, 11959, 11953, 12055, 12144, 12079, 22964)

#meta_Nash_filter <- filter(meta_Nash, !animal_biosample %in% high_rich)


#h0_10 <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 | sample_name) + (1 | dna_experiment_date) + 
#                   (1 | collection_date) + (1 | plate_ID) + (1 | well_location), data = meta_Nash, family = poisson())

#h0_10_filter <- glmmTMB(h0 ~ age_at_collection + sex + SMI + (1 | sample_name) + (1 | dna_experiment_date) + 
#                          (1 | collection_date) + (1 | plate_ID) + (1 | well_location), data = meta_Nash_filter, family = poisson())

#AIC(h0_10, h0_10_filter)

#overdisp_ratio = deviance(h0_10_filter)/df.residual(h0_10_filter)
#overdisp_ratio

#summary(h0_10_filter)

#group_means <- emmeans(h0_10_filter, ~ age_at_collection + sex + SMI)
#group_means

#plot(residuals(h0_10_filter))

#sim_res <- simulateResiduals(h0_10_filter)
#plot(sim_res)

#Anova(h0_10_filter, type = 2)
#SMI/body condition has a significant effect on richness but age and sex do not

#r.squaredGLMM(h0_10_filter)
