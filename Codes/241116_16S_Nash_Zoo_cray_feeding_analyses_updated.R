#Researcher: Chloe Cummins
#Start Date: November 16, 2024
#Updating 16S Nashville Zoo crayfish feeding thesis analyses based on discussions 
#with lab members during November 2024


#Set working directory to manuscript code check folder
setwd("G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Hellbender Gut Microbiome Manuscript Stuff/Manuscript_Draft_Code_and_Figs")

#Previous code stored in "G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Hellbender Gut Microbiome Manuscript Stuff/Manuscript_Code_Check"
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
library(gapminder)
library(ape)
library(plotly)
library(phyloseq)
library(forcats)
library(hillR)
library(betapart)
library(geomtextpath)
library(reshape2)
library(ggbreak)
library(car)
library(DHARMa)
library(performance)
library(lmodel2)  #for SMI calculation
library(glmmTMB)
library(emmeans)
library(buildmer)
library(MuMIn)
library(AICcmodavg)

#Note: this is a fairly large library/package - needed for LEfSe analysis
library(microbiomeMarker)


#Importing the decontam data
otu_count <- fread("240820_HB_16S_OTU_count_table_decontam.csv")
tax_table <- fread("240820_HB_16S_taxonomy_decontam.csv")
metadata <- fread("240820_HB_16S_metadata_decontam.csv")

#Setting seed for reproducibility
set.seed(012219)



##### Preparing the data for analysis #####

#Converting all of the data into data frames
otu_table_2 <- as.data.frame(otu_count)

tax_table_2 <- as.data.frame(tax_table)

metadata_2 <- as.data.frame(metadata)

#Filtering the metadata to only include Nashville Zoo cloacal swab samples 
#from animals with complete crayfish feeding data - hatch year 2018
meta_Nash_cray <- filter(metadata_2, site == "Nashville Zoo" & env_medium == "cloaca swab" &
                           hatch_date == "10/5/2018")

#Making sure sample count is the same for the taxonomy and OTU tables after filtering out samples
#If the sizes of the metadata, OTU, and taxonomy tables are not the same after you filter samples out of the metadata
#some analyses will throw errors
meta_otu_inter <- intersect(meta_Nash_cray$index, otu_count$index)
otu_table_2 <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table since it is likely a contaminant (present even in neg. controls)
otu_table_2 <- otu_table_2[,-2] #filtering out the OTU00001 column which represents E. Shigella (based on tax_table_2)

#Format the OTU table for the better plotting, joining, etc.  
otu_table_2 <- otu_table_2 %>% 
  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_2))
tax_table_2 <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_2 <- otu_table_2 %>% select(all_of(intersected))

#Setting certain variables as factors for analysis
meta_Nash_cray <- mutate(meta_Nash_cray, 
                         diet_crayfish = factor(diet_crayfish),
                         sex = factor(sex, levels = c("M", "F")), 
                         collection_date = factor(collection_date), 
                         last_crayfish_feeding = factor(last_crayfish_feeding), 
                         sample_name = factor(sample_name))

#Adding a column to the data frame to allow for better plotting of collection date (since there were two days in Feb. 2023)
meta_Nash_cray <- meta_Nash_cray %>% 
  mutate(collect_date_graph = case_when(collection_date == "2/22/2023" ~ "Feb. 2023", 
                                        collection_date == "2/21/2023" ~ "Feb. 2023", 
                                        collection_date == "5/1/2023" ~ "May 2023", 
                                        collection_date == "11/7/2023" ~ "Nov. 2023",
                                        collection_date == "1/23/2024" ~ "Jan. 2024",
                                        collection_date == "4/23/2024" ~ "April 2024"),
         collect_date_graph = factor(collect_date_graph, order = T, 
                                     levels = c("Feb. 2023", "May 2023", "Nov. 2023", "Jan. 2024", "April 2024")))
#As a note, we ended up deciding to focus on the middle three time points for most of these analyses due to the crayfish feeding experimental design



##### Calculating Hill numbers for alpha diversity and adding to the metadata file #####

#All of the available Hill numbers were calculated here in case they were needed for anything
#but all of my analyses mainly use richness (h0)

meta_Nash_cray$h0 <- hill_taxa(otu_table_2, q = 0, MARGIN = 1)
#q = 0 is essentially a richness metric, placing more emphasis on rare taxa
#Margin is set to the default 1 in this case since the function expects sites, or in this case samples, to be rows

#meta_Nash_cray$h1 <- hill_taxa(otu_table_2, q = 1, MARGIN = 1)
#q = 1 is Hill-Shannon diversity, where abundant and rare taxa are equally weighted
#and evenness is accounted for

#meta_Nash_cray$h2 <- hill_taxa(otu_table_2, q = 2, MARGIN = 1)
#q = 2 is inverse Simpson, which places more emphasis on common/abundant species

#meta_Nash_cray$h3 <- hill_taxa(otu_table_2, q = 3, MARGIN = 1)
#q = 3 is where dominant, abundant species are given the most weight



##### Scaled Mass Index (SMI) Analyses #####

#Natural log transforming the mass and TL data in order to perform standardized major axis (SMA) regression
#per the equation and steps outlined in Peig and Green 2009
#Link to paper: https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0706.2009.17643.x
meta_Nash_cray$ln_mass <- log(meta_Nash_cray$weight)

meta_Nash_cray$ln_TL <- log(meta_Nash_cray$total_length)

#Performing a model II SMA regression
sma <- lmodel2(ln_mass ~ ln_TL, data = meta_Nash_cray, nperm = 99)

bsma <- sma[["regression.results"]][["Slope"]][[3]]
bsma
#SMA slope =  3.147391
#this slope is the scaling exponent for the SMI calculation

l0 <- mean(meta_Nash_cray$total_length, na.rm = TRUE)
l0
#mean TL = 30.49381 cm

#Calculating SMI based on the equation
meta_Nash_cray$SMI <- (meta_Nash_cray$weight * (l0/meta_Nash_cray$total_length) ^ bsma)

#meta_Nash_cray$scaled_SMI <- scale(meta_Nash_cray$SMI)



##### Adding a column to the metadata for crayfish feeding groups #####

#Specifying Group 1 animals - indv. fed crayfish by Nov. 2023 (2nd time point)
group1 <- c("UHM759", "UHM748", "UHM779", "UHM818", "UHM820", "UHM811", "UHM810",
            "UHM819", "UHM806", "UHM827", "UHM784", "UHM829")

#Filtering the metadata and adding a column for Group 1 individuals
meta_Nash_cray_g1 <- meta_Nash_cray %>% 
  subset(sample_name %in% group1) %>% 
  mutate(cray_group = "Group 1")

#Double checking that the numbers for the group are what we would expect
group_broad_summary <- meta_Nash_cray_g1 %>% 
  group_by(diet_crayfish, collection_date, sex, last_crayfish_feeding) %>% 
  tally()
group_broad_summary


#Specifying Group 2 animals - indv. fed crayfish by last full sampling event in Jan. 2024
group2 <- c("UHM747", "UHM746", "UHM749", "UHM782", "UHM783", "UHM813", "UHM776",
            "UHM832", "UHM777", "UHM775")

#Filtering the metadata and adding a column for Group 2 individuals
meta_Nash_cray_g2 <- meta_Nash_cray %>% 
  subset(sample_name %in% group2) %>% 
  mutate(cray_group = "Group 2")

#Double checking that the numbers for the group are what we would expect
group_broad_summary <- meta_Nash_cray_g2 %>% 
  group_by(diet_crayfish, collection_date, sex, last_crayfish_feeding) %>% 
  tally()
group_broad_summary


#Specifying the few animals in the female group, which were fed very early on, unlike the males
#female_group <- c("UHM755", "UHM763", "UHM754", "UHM762", "UHM764", "UHM761", "UHM756")

#Filtering the metadata for Female Group individuals
#meta_Nash_cray_female <- meta_Nash_cray %>% 
#  subset(sample_name %in% female_group) %>% 
#  mutate(cray_group = "Female Group")

#Double checking that the numbers for the group are what we would expect
#group_broad_summary <- meta_Nash_cray_female %>% 
#  group_by(diet_crayfish, collection_date, sex, last_crayfish_feeding) %>% 
#  tally()
#group_broad_summary
#As a note, some of the females were inconsistently fed and went long timespans
#from last cray feeding to the sampling date, unlike many of the males

#Joining all of the group data frames together
#meta_cray_feeding_full <- rbind(meta_Nash_cray_g1, meta_Nash_cray_g2, meta_Nash_cray_female)
#meta_cray_feeding_full$cray_group <- as.factor(meta_cray_feeding_full$cray_group)

#Creating a just male data frame for modeling 
meta_cray_feedingM <- rbind(meta_Nash_cray_g1, meta_Nash_cray_g2)
meta_cray_feedingM$cray_group <- as.factor(meta_cray_feedingM$cray_group)



##### Creating time point data frames for modeling and analysis #####

#NOTE: I made multiple variations of dataframes for the crayfish feeding groups since we were still
#in the process of solidifying the overall design we wanted to proceed with, many of these can be ignored


#Data frame for crayfish feeding groups for 2 time points
#meta_cray_2t <- filter(meta_cray_feeding_full, collection_date == "11/7/2023" |
#                                collection_date == "1/23/2024")

  #all males
#meta_cray_2tM <- filter(meta_cray_feedingM, collection_date == "11/7/2023" |
#                                collection_date == "1/23/2024")

#Data frame for crayfish feeding groups for 3 time points - full data, including females
#meta_cray_3t <- filter(meta_cray_feeding_full, collection_date == "5/1/2023" |
#                                collection_date == "11/7/2023" | collection_date == "1/23/2024")

  #all males
meta_cray_3tM <- filter(meta_cray_feedingM, collection_date == "5/1/2023" |
                                collection_date == "11/7/2023" | collection_date == "1/23/2024")

#Data frame for crayfish feeding groups for 3 time points - reduced data
#The reduced data simply follows the same set of indv. across time (they are not missing any time points)
#set.seed(012219) #for reproducibility
#group1_reduced <- c("UHM748", "UHM779", "UHM818", "UHM827", "UHM829")

#Deleting UHM775 from the vector because it was missing a time point
#group2_new <- group2[-10]
#group2_reduced <- sample(group2_new, size = 5)

#female_new <- female_group[-5]
#female_reduced <- sample(female_new, size = 5)
  #UHM764 had to be excluded because it was missing a time point

#meta_cray_3t_reduced <- filter(meta_cray_feeding_full, sample_name %in% group1_reduced |
#                                 sample_name %in% group2_reduced | sample_name %in% female_reduced)

#meta_cray_3t_reduced <- filter(meta_cray_3t_reduced, collection_date == "5/1/2023" |
#                                 collection_date == "11/7/2023" | collection_date == "1/23/2024")

  #all males
#meta_cray_3tM_reduced <- filter(meta_cray_feeding_full, sample_name %in% group1_reduced |
#                                 sample_name %in% group2_reduced)

#meta_cray_3tM_reduced <- filter(meta_cray_3tM_reduced, collection_date == "5/1/2023" |
#                                 collection_date == "11/7/2023" | collection_date == "1/23/2024")


#Creating all male data frame for 1st and last time points (May 2023 and Jan 2024)
#this is for the before and after crayfish comparison (including betapart and LEfSe analysis)
meta_cray_2tM_bva <- filter(meta_cray_feedingM, collection_date == "5/1/2023" |
                          collection_date == "1/23/2024")


##### Statistical Modeling - Alpha Diversity #####

#1. Evaluating whether crayfish feeding across 2 time points has an impact on richness via Hill number (h0)-----

#Assessing normality of richness (response variable) among the crayfish feeding groups at Nash Zoo
#shapiro.test(meta_cray_2tM$h0) #p-value < 4.645e-12, not normally distributed
#Data is not normally distributed, indicating that a GLMM may produce a better fit than an LMM

#Looking at the distribution of the response variable to double check that the stats. match
#densityplot(meta_cray_2tM$h0, main = "Distribution of Richness")
#a few outliers seem to be forming a right skew, otherwise potentially seems normal

#Histogram of richness
#hist(meta_cray_2tM$h0)

#Q-Q plot of richness
#qqnorm(meta_cray_2tM$h0)
#qqline(meta_cray_2tM$h0, col = "red")
#there seem to be some outliers near the top of the graph

#Q-Q plot to see which feeding group is causing the outliers in the previous plot
#h0_qq <- ggplot(meta_cray_2tM, aes(sample = h0, color = cray_group)) + 
#  stat_qq() +
#  stat_qq_line() +
#  theme_classic()
#print(h0_qq)
#outliers seem to be in all groups


  ## GLMM of diet's effect on richness ##

#Evaluating whether I should use a tweedie distribution or a poisson distribution
#h0_1 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name), data = meta_cray_2tM, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
#h0_2 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name), data = meta_cray_2tM, family = poisson())

#AIC(h0_1, h0_2)
#summary(h0_2)
#Anova(h0_2, type = 2)
#the tweedie distribution produced a false convergence issue so I added the control argument
#but based on what I understand tweedie is better suited to continuous data with lots of zeros
#and may not be most appropriate in this case since the data is not overdispersed

#plot(residuals(h0_2))

#sim_res <- simulateResiduals(h0_2)
#plot(sim_res)


#Evaluating the addition of more random effects according to potential importance
#h0_3 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + (1 | dna_experiment_date),
#                data = meta_cray_2tM, family = poisson())

#h0_4 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_cray_2tM, family = poisson())

#h0_5 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID) + (1 | SMI), data = meta_cray_2tM, family = poisson())

#h0_6 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + 
#                  (1 + dna_experiment_date | plate_ID) + (1 | SMI), data = meta_cray_2tM, family = poisson())

#Testing to the tweedie model
#h0_7 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + (1 | dna_experiment_date), data = meta_cray_2tM, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

#AIC(h0_3, h0_4, h0_5, h0_6, h0_7)

#overdisp_ratio = deviance(h0_5)/df.residual(h0_5)
#overdisp_ratio

#summary(h0_5)
#Anova(h0_5, type = 2)

#group_means <- emmeans(h0_5, ~ cray_group * collection_date)
#group_means

#emmeans(h0_5, pairwise ~ cray_group | collection_date)

#plot(residuals(h0_5))

#sim_res <- simulateResiduals(h0_5)
#plot(sim_res)


#2. Evaluating whether crayfish feeding across 3 time points with reduced data has an impact on richness via Hill number (h0)-----

#Assessing normality of richness (response variable) among the crayfish feeding groups at Nash Zoo
#shapiro.test(meta_cray_3tM_reduced$h0) #p-value < 7.72e-09, not normally distributed
#Data is not normally distributed, indicating that a GLMM may produce a better fit than an LMM

#Looking at the distribution of the response variable to double check that the stats. match
#densityplot(meta_cray_3tM_reduced$h0, main = "Distribution of Richness")
#a few outliers seem to be forming a right skew, otherwise potentially seems normal

#Histogram of richness
#hist(meta_cray_3tM_reduced$h0)

#Q-Q plot of richness
#qqnorm(meta_cray_3tM_reduced$h0)
#qqline(meta_cray_3tM_reduced$h0, col = "red")
#there seem to be some outliers near the top of the graph

#Q-Q plot to see which feeding group is causing the outliers in the previous plot
#h0_qq <- ggplot(meta_cray_3tM_reduced, aes(sample = h0, color = cray_group)) + 
#  stat_qq() +
#  stat_qq_line() +
#  theme_classic()
#print(h0_qq)
#outliers seem to be in all groups


## GLMM of diet's effect on richness ##

#Evaluating whether I should use a tweedie distribution or a poisson distribution
#h0_1 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name), data = meta_cray_3tM_reduced, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
#h0_2 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name), data = meta_cray_3tM_reduced, family = poisson())

#AIC(h0_1, h0_2)
#summary(h0_2)
#the tweedie distribution produced a false convergence issue so I added the control argument
#but based on what I understand tweedie is better suited to continuous data with lots of zeros
#and may not be most appropriate in this case since the data is not overdispersed

#plot(residuals(h0_2))

#sim_res <- simulateResiduals(h0_2)
#plot(sim_res)


#Evaluating the addition of more random effects according to potential importance
#h0_3 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + (1 | dna_experiment_date),
#                data = meta_cray_3tM_reduced, family = poisson())

#h0_4 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_cray_3tM_reduced, family = poisson())

#h0_5 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID) + (1 | SMI), data = meta_cray_3tM_reduced, family = poisson())

#h0_6 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + 
#                  (1 + dna_experiment_date | plate_ID) + (1 | SMI), data = meta_cray_3tM_reduced, family = poisson())

#Testing to the tweedie model
#h0_7 <- glmmTMB(h0 ~ cray_group * collection_date + (1 | sample_name) + (1 | dna_experiment_date), data = meta_cray_3tM_reduced, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

#AIC(h0_3, h0_4, h0_5, h0_6, h0_7)

#overdisp_ratio = deviance(h0_5)/df.residual(h0_5)
#overdisp_ratio

#summary(h0_5)

#group_means <- emmeans(h0_5, ~ cray_group * collection_date)
#group_means

#emmeans(h0_5, pairwise ~ cray_group | collection_date)

#plot(residuals(h0_5))

#sim_res <- simulateResiduals(h0_5)
#plot(sim_res)



#3. Evaluating whether crayfish feeding across 3 time points has an impact on richness via Hill number (h0)-----

meta_cray_3tM <- meta_cray_3tM %>% 
  mutate(collect_point = case_when(collection_date == "5/1/2023" ~ 1, 
                                   collection_date == "11/7/2023" ~ 2, 
                                   collection_date == "1/23/2024" ~ 3))

#meta_cray_3tM$collect_point <- as.Date(meta_cray_3tM$collection_date, "%m/%d/%Y")
#meta_cray_3tM$collect_point <- as.POSIXct.Date(meta_cray_3tM$collect_point, format = "%Y-%m-%d")


#Assessing normality of richness (response variable) among the crayfish feeding groups at Nash Zoo
shapiro.test(meta_cray_3tM$h0) #p-value < 4.94-13, not normally distributed
#Data is not normally distributed, indicating that a GLMM may produce a better fit than an LMM

#Looking at the distribution of the response variable to double check that the stats. match
densityplot(meta_cray_3tM$h0, main = "Distribution of Richness")
#a few points seem to be forming a right skew, otherwise potentially seems somewhat normal

#Histogram of richness
hist(meta_cray_3tM$h0)

#Q-Q plot of richness
qqnorm(meta_cray_3tM$h0)
qqline(meta_cray_3tM$h0, col = "red")
#there seem to be some outliers in the response variable near the top of the graph

#Q-Q plot to see which points correspond to which feeding group
h0_qq <- ggplot(meta_cray_3tM, aes(sample = h0, color = cray_group)) + 
  stat_qq() +
  stat_qq_line() +
  theme_classic()
print(h0_qq)
#outliers seem to be in both groups


## GLMM of diet's effect on richness ##

#Evaluating whether I should use a tweedie distribution or a poisson distribution
#testing a tweedie distribution was suggested by one of our postdocs, Mitra
#h0_1 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name), data = meta_cray_3tM, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
h0_2 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name), data = meta_cray_3tM, family = poisson())

AIC(h0_2)
AICc(h0_2) #Adding AICc from here forward to account for small sample sizes
summary(h0_2)
#The inclusion of tweedie here is mainly an artifact in the code where I tested to see if it may fit well - based on the recommendation of one of my lab members
#However, based on what I understand tweedie is better suited to continuous data with lots of zeros
#and therefore may not be the most appropriate in this case since richness is discrete and should not be zero-inflated
#So I commented tweedie out but kept it in the code to avoid messing up the numbering of the rest of the models

#Quick test of the residuals before specifying the model further
plot(residuals(h0_2))

sim_res <- simulateResiduals(h0_2)
plot(sim_res)


#Evaluating the addition of more random effects according to potential importance
h0_3 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date),
                data = meta_cray_3tM, family = poisson())

h0_4 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
                  (1 | plate_ID), data = meta_cray_3tM, family = poisson())

h0_5 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
                  (1 | plate_ID) + (1 | SMI), data = meta_cray_3tM, family = poisson())
#I originally included SMI as a random effect in these models but changed it to a fixed effect (see additional code below)
#based on a conversation with my committee members,
#and it was suggested to test SMI as a fixed effect (potentially through stepwise selection) since it is a
#continuous variable to see if it should be included as fixed or random.

#trying a nested random effect for h0_6
h0_6 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + 
                  (1 + dna_experiment_date | plate_ID) + (1 | SMI), data = meta_cray_3tM, family = poisson())
#h0_6 gave a model convergence warning - do not use

AIC(h0_3, h0_4, h0_5, h0_6) #I assessed the AIC of h0_6 to double check the model convergence issue since the AIC would be NA in that case
AICc(h0_3)
AICc(h0_4)
AICc(h0_5)
#h0_5 has the lowest AIC and lowest AICc

#Calculating whether the model may be overdispersed (ratio greater than 1)
overdisp_ratio = deviance(h0_5)/df.residual(h0_5)
overdisp_ratio
#no overdispersion

summary(h0_5)

group_means <- emmeans(h0_5, ~ cray_group * collect_point)
group_means

emmeans(h0_5, pairwise ~ cray_group | collect_point)

plot(residuals(h0_5))

sim_res <- simulateResiduals(h0_5)
plot(sim_res)

Anova(h0_5, type = 2)

#Stepwise selection with SMI as a fixed effect - based on suggestions by committee
#I basically took the h0_5 model structure and changed SMI to a fixed effect
h0_1_step <- buildglmmTMB(h0 ~ cray_group * collect_point + SMI + (1 | sample_name) + 
                  (1 + dna_experiment_date | plate_ID) + (1 | dna_experiment_date) + (1 | plate_ID), 
                  data = meta_cray_3tM, family = poisson(), buildmerControl = 
                    buildmerControl(crit = "AIC"))  #I had to use AIC here since AICc is not available for this function
(f <- formula(h0_1_step@model)) #returning the model suggested by stepwise selection

#Stepwise selection with SMI as a random effect
h0_2_step <- buildglmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + 
                            (1 + dna_experiment_date | plate_ID) + (1 | dna_experiment_date) + 
                            (1 | plate_ID) + (1 | SMI), 
                          data = meta_cray_3tM, family = poisson(), buildmerControl = 
                            buildmerControl(crit = "AIC"))
(f <- formula(h0_2_step@model))
#Crayfish feeding group was dropped during model selection for this stepwise - making it not an ideal model 
#since there is theoretical justification to keep that fixed effect in the model because that is one of the variables that we want to test

#Testing the model identified by stepwise selection that includes SMI as a fixed effect
h0_7 <- glmmTMB(h0 ~ collect_point + cray_group + collect_point * cray_group + 
                  SMI + (1 | sample_name) + (1 | dna_experiment_date), data = meta_cray_3tM, family = poisson())

#Including plate ID as an additional random effect similar to h0_5
h0_8 <- glmmTMB(h0 ~ cray_group * collect_point + SMI + (1 | sample_name) + (1 | dna_experiment_date) +
                  (1 | plate_ID), data = meta_cray_3tM, family = poisson())
#Both models returned a warning message, but from what I understand based on the glmmTMB documentation, 
#this warning may be ok as long as it is not accompanied by a convergence error/warning. 

#Link to documentation (https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html)
#"This warning occurs when the optimizer visits a region of parameter space that is invalid. 
#It is not a problem as long as the optimizer has left that region of parameter space upon convergence, 
#which is indicated by an absence of the model convergence warnings described above."

AIC(h0_7, h0_8)
AICc(h0_7)
AICc(h0_8)
#h0_7 has a lower AICc

summary(h0_7)

plot(residuals(h0_7))

#Creating a normal QQ plot
qqnorm(residuals(h0_7))
qqline(residuals(h0_7))

#Simulating and testing DHARMa residuals
sim_res <- simulateResiduals(h0_7)
plot(sim_res) #I mainly paid attention to the QQ plot here
testDispersion(h0_7) #DHARMa indicates that there is no significant overdispersion
testUniformity(h0_7) #KS, dispersion, and outlier tests are all non-significant

overdisp_ratio = deviance(h0_7)/df.residual(h0_7)
overdisp_ratio
#ratio = 5.9 -  potentially overdispersed??? but DHARMa calculated a different dispersion value and found that dispersion was not significant
#I will attempt to fit a negative binomial distribution as a just in case later in the code to see if that improves anything
#but overall I feel like the DHARMa results may be more reliable than my hand calculation above since DHARMa uses a simulation-based approach

group_means <- emmeans(h0_7, ~ cray_group * collect_point + SMI)
group_means

r.squaredGLMM(h0_7)

#testOutliers(h0_7)
#testOutliers(sim_res)
#check_outliers(sim_res, type = "default")
#check_outliers(sim_res, type = "binomial")
#check_outliers(sim_res, type = "bootstrap")

#Attempting to fit negative binomial to see if it is any better
h0_9 <- glmmTMB(h0 ~ collect_point * cray_group + 
                  SMI + (1 | sample_name) + (1 | dna_experiment_date), data = meta_cray_3tM, family = nbinom1())
#nbinom1 uses a linear parameterization based on the glmmTMB family distribution help page 

h0_10 <- glmmTMB(h0 ~ collect_point * cray_group + 
                  SMI + (1 | sample_name) + (1 | dna_experiment_date), data = meta_cray_3tM, family = nbinom2())
#nbinom1 uses a quadratic parameterization based on the glmmTMB family distribution help page 

AIC(h0_9, h0_10)
AICc(h0_9)
AICc(h0_10)
#AIC lower for nbinom2 for h0_10

summary(h0_10)

group_means <- emmeans(h0_10, ~ cray_group * collect_point)
group_means

emmeans(h0_10, pairwise ~ cray_group | collect_point)

plot(residuals(h0_10))
#residuals look a little worse, there are more scattered points near the upper end of the residuals

#Creating a normal QQ plot
qqnorm(residuals(h0_10))
qqline(residuals(h0_10))
#QQ plot also looks worse and has more deviation in one of the tails

sim_res <- simulateResiduals(h0_10)
plot(sim_res)
testDispersion(h0_10) #DHARMa indicates that there is significant overdispersion
testUniformity(h0_10) #KS and outlier test are not significant but the dispersion test is significant
 

Anova(h0_10, type = 2)

#Overall, the negative binomial model (h0_10) has a lower AIC, but it does not
#seem to improve model fit that well based on the residuals compared to Poisson distribution
#Also, this is somewhat confusing since I was receiving conflicting information about the dispersion in the Poisson model (h0_7)


#Conclusion: I decided to proceed with h0_7 since it seems to have better residuals and QQ plots 
#as well as no issues with DHARMa residuals/QQ plot (no significant KS, outlier, or dispersion). 




#4. Evaluating whether crayfish feeding across 3 time points has an impact on Shannon diversity via Hill number (h1)-----

#Assessing normality of Shannon (response variable) among the crayfish feeding groups at Nash Zoo
#shapiro.test(meta_cray_3tM$h1) #p-value < 7.935e-12, not normally distributed
#Data is not normally distributed, indicating that a GLMM may produce a better fit than an LMM

#Looking at the distribution of the response variable to double check that the stats. match
#densityplot(meta_cray_3tM$h1, main = "Distribution of Shannon Diversity")
#a few outliers seem to be forming a right skew, otherwise potentially seems normal

#Histogram of Shannon
#hist(meta_cray_3tM$h1)

#Q-Q plot of richness
#qqnorm(meta_cray_3tM$h1)
#qqline(meta_cray_3tM$h1, col = "red")
#there seem to be some outliers near the top of the graph

#Q-Q plot to see which feeding group is causing the outliers in the previous plot
#h1_qq <- ggplot(meta_cray_3tM, aes(sample = h1, color = cray_group)) + 
#  stat_qq() +
#  stat_qq_line() +
#  theme_classic()
#print(h1_qq)
#outliers seem to be in all groups

#Stepwise selection with SMI as a fixed effect
#h1_1_step <- buildglmmTMB(h1 ~ cray_group * collect_point + SMI + (1 | sample_name) + 
#                            (1 + dna_experiment_date | plate_ID) + (1 | dna_experiment_date) + (1 | plate_ID), 
#                          data = meta_cray_3tM, family = Gamma(), buildmerControl = 
#                            buildmerControl(crit = "AIC"))
#(f <- formula(h1_1_step@model))

#Stepwise selection with SMI as a random effect
#h1_2_step <- buildglmmTMB(h1 ~ cray_group * collect_point + (1 | sample_name) + 
#                            (1 + dna_experiment_date | plate_ID) + (1 | dna_experiment_date) + 
#                            (1 | plate_ID) + (1 | SMI), 
#                          data = meta_cray_3tM, family = Gamma(), buildmerControl = 
#                            buildmerControl(crit = "AIC"))
#(f <- formula(h1_2_step@model))

#Comparing gamma and tweedie distributions
#h1_1 <- glmmTMB(h1 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID) + (1 | SMI), data = meta_cray_3tM, family = Gamma())

#h1_2 <- glmmTMB(h1 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date), data = meta_cray_3tM, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

#AIC(h1_1, h1_2)

 
#h1_3 <- glmmTMB(h1 ~ cray_group * collect_point + SMI + (1 | sample_name) + (1 | dna_experiment_date), data = meta_cray_3tM, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

#AIC(h1_2, h1_3)

#summary(h1_2)

#group_means <- emmeans(h1_2, ~ cray_group * collect_point)
#group_means

#emmeans(h1_2, pairwise ~ cray_group | collect_point)

#plot(residuals(h1_2))

#sim_res <- simulateResiduals(h1_2)
#plot(sim_res)


#5. Evaluating whether crayfish feeding across 2 time points (before vs. after) has an impact on richness via Hill number (h0)-----
#meta_cray_2tM_bva <- meta_cray_2tM_bva %>% 
#  mutate(collect_point = case_when(collection_date == "5/1/2023" ~ 1, 
#                                   collection_date == "11/7/2023" ~ 2, 
#                                   collection_date == "1/23/2024" ~ 3))

#meta_cray_2tM_bva$collect_point <- as.Date(meta_cray_2tM_bva$collection_date, "%m/%d/%Y")
#meta_cray_2tM_bva$collect_point <- as.POSIXct.Date(meta_cray_2tM_bva$collect_point, format = "%Y-%m-%d")


#Assessing normality of richness (response variable) among the crayfish feeding groups at Nash Zoo
#shapiro.test(meta_cray_2tM_bva$h0) #p-value < 2.421e-10, not normally distributed
#Data is not normally distributed, indicating that a GLMM may produce a better fit than an LMM

#Looking at the distribution of the response variable to double check that the stats. match
#densityplot(meta_cray_2tM_bva$h0, main = "Distribution of Richness")
#a few outliers seem to be forming a right skew, otherwise potentially seems normal

#Histogram of richness
#hist(meta_cray_2tM_bva$h0)

#Q-Q plot of richness
#qqnorm(meta_cray_2tM_bva$h0)
#qqline(meta_cray_2tM_bva$h0, col = "red")
#there seem to be some outliers near the top of the graph

#Q-Q plot to see which feeding group is causing the outliers in the previous plot
#h0_qq <- ggplot(meta_cray_2tM_bva, aes(sample = h0, color = cray_group)) + 
#  stat_qq() +
#  stat_qq_line() +
#  theme_classic()
#print(h0_qq)
#outliers seem to be in group 1


## GLMM of diet's effect on richness ##

#Evaluating whether I should use a tweedie distribution or a poisson distribution
#h0_1 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name), data = meta_cray_2tM_bva, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
#h0_2 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name), data = meta_cray_2tM_bva, family = poisson())

#AIC(h0_1, h0_2)
#summary(h0_2)


#plot(residuals(h0_2))

#sim_res <- simulateResiduals(h0_2)
#plot(sim_res)


#Evaluating the addition of more random effects according to potential importance
#h0_3 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date),
#                data = meta_cray_2tM_bva, family = poisson())

#h0_4 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_cray_2tM_bva, family = poisson())

#h0_5 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID) + (1 | SMI), data = meta_cray_2tM_bva, family = poisson())

#h0_6 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + 
#                  (1 + dna_experiment_date | plate_ID) + (1 | SMI), data = meta_cray_2tM_bva, family = poisson())

#AIC(h0_3, h0_4, h0_5, h0_6)

#overdisp_ratio = deviance(h0_5)/df.residual(h0_5)
#overdisp_ratio
#no overdispersion based on ratio

#summary(h0_5)

#group_means <- emmeans(h0_5, ~ cray_group * collect_point)
#group_means

#emmeans(h0_5, pairwise ~ cray_group | collect_point)

#plot(residuals(h0_5))

#sim_res <- simulateResiduals(h0_5)
#plot(sim_res)
#testDispersion(h0_5)#DHARMa says there is overdispersion

#Anova(h0_5, type = 2)

#Stepwise selection with SMI as a fixed effect
#h0_1_step <- buildglmmTMB(h0 ~ cray_group * collect_point + SMI + (1 | sample_name) + 
#                            (1 + dna_experiment_date | plate_ID) + (1 | dna_experiment_date) + (1 | plate_ID), 
#                          data = meta_cray_2tM_bva, family = poisson(), buildmerControl = 
#                            buildmerControl(crit = "AIC"))
#(f <- formula(h0_1_step@model))

#Stepwise selection with SMI as a random effect
#h0_2_step <- buildglmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + 
#                            (1 + dna_experiment_date | plate_ID) + (1 | dna_experiment_date) + 
#                            (1 | plate_ID) + (1 | SMI), 
#                          data = meta_cray_2tM_bva, family = poisson(), buildmerControl = 
#                            buildmerControl(crit = "AIC"))
#(f <- formula(h0_2_step@model))
#Crayfish feeding group was dropped during model selection for this step - making it not an ideal model

#Using the model identified by stepwise selection
#h0_7 <- glmmTMB(h0 ~ collect_point + cray_group + collect_point * cray_group +
#                  (1 | sample_name), data = meta_cray_2tM_bva, family = poisson())

#Using SMI as a fixed effect as a just in case
#h0_8 <- glmmTMB(h0 ~ cray_group * collect_point + SMI + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_cray_2tM_bva, family = poisson())

#AIC(h0_5, h0_7, h0_8)
#h0_5 has lower AIC

#summary(h0_7)

#plot(residuals(h0_7))

#sim_res <- simulateResiduals(h0_7)
#plot(sim_res)
#testDispersion(h0_7) #DHARMa indicates that there is overdispersion

#overdisp_ratio = deviance(h0_7)/df.residual(h0_7)
#overdisp_ratio
#ratio = 3.1 -  potentially overdispersed
#Both h0_5 and h0_7 may be overdispersed but h0_5 has a better fit based on residuals and AIC

#Trying to fit a Tweedie distribution and negative binomial distribution to see if they work better
#h0_9 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_cray_2tM_bva, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

#h0_10 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                   (1 | plate_ID), data = meta_cray_2tM_bva, family = nbinom1())

#h0_11 <- glmmTMB(h0 ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                   (1 | plate_ID), data = meta_cray_2tM_bva, family = nbinom2())

#AIC(h0_9, h0_10, h0_11, h0_5)
#AIC the same for both - will use nbinom1

#summary(h0_9)
 

#group_means <- emmeans(h0_9, ~ cray_group * collect_point)
#group_means

#emmeans(h0_9, pairwise ~ cray_group | collect_point)

#plot(residuals(h0_11))

#sim_res <- simulateResiduals(h0_11)
#plot(sim_res)
#Tweedie and nbinom2 (tested separately) seem to produce similar results to Poisson

#Anova(h0_9, type = 2)


#Testing a negative binomial distribution with SMI included as a fixed effect
#similar to the model used for 3 time points (just diff. distribution based on overdispersion for only 2 time points)
#h0_12 <- glmmTMB(h0 ~ cray_group * collect_point + SMI + (1 | sample_name) + (1 | dna_experiment_date), 
#                 data = meta_cray_2tM_bva, family = poisson())

#h0_13 <- glmmTMB(h0 ~ cray_group * collect_point + SMI + (1 | sample_name) + (1 | dna_experiment_date), 
#                 data = meta_cray_2tM_bva, family = nbinom2())


#AIC(h0_5, h0_9, h0_11, h0_12, h0_13)

#summary(h0_13)

#plot(residuals(h0_13))

#sim_res <- simulateResiduals(h0_13)
#plot(sim_res)

#overdisp_ratio = deviance(h0_13)/df.residual(h0_13)
#overdisp_ratio

#Anova(h0_13, type = 2)



##### Bray-Curtis NMDS for Males vs. Females #####

#Making sure sample count is the same after filtering out samples
#meta_otu_inter <- intersect(meta_cray_3t$index, otu_count$index)
#otu_table_3t <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
#otu_table_3t <- otu_table_3t[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_3t <- otu_table_3t %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected <- intersect(tax_table_2$otu, colnames(otu_table_3t))
#tax_table_3t <- tax_table_2 %>% filter(otu %in% intersected)
#otu_table_3t <- otu_table_3t %>% select(all_of(intersected))

#Preparing the NMDS
#set.seed(012219) #set seed for reproducibility
#bray_dist_MvF <- vegdist(otu_table_3t, method = "bray", na.rm = T)

#Checking for homogeneity of dispersion
#bray_beta_diet_MvF <- betadisper(bray_dist_MvF, meta_cray_3t$cray_group)
#permutest(bray_beta_diet_MvF)
#TukeyHSD(bray_beta_diet_MvF)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.438 :)

#Checking for homogeneity of dispersion
#bray_beta_time_MvF <- betadisper(bray_dist_MvF, meta_cray_3t$collection_date)
#permutest(bray_beta_time_MvF)
#TukeyHSD(bray_beta_time_MvF)
#no significant difference in dispersion between sampling dates, p = 0.421 :)


#Checking for homogeneity of dispersion
#bray_beta_MvF <- betadisper(bray_dist_MvF, meta_cray_3t$sex)
#permutest(bray_beta_MvF)
#TukeyHSD(bray_beta_MvF)
#no significant difference in dispersion between M/F, p = 0.721 :)

#Performing the PERMANOVA
#bray_adonis_MvF <- adonis2(bray_dist_MvF ~ h0 * sex + collection_date + SMI, data = meta_cray_3t,
#                            permutations = 999, by = "terms")
#bray_adonis_MvF
#no significant effect of any of the variables

#bray_nmds_MvF <- metaMDS(bray_dist_MvF, k = 3) #conduct non-metric multidimensional scaling

#bray_nmds_MvF$stress
#stress value is 0.0981499, which is fairly decent

#bray_nmds_df_MvF <- as.data.frame(bray_nmds_MvF[["points"]]) %>% #generate dataframe with NMDS points 
#  rownames_to_column("index") %>% 
#  left_join(meta_cray_3t, by = join_by("index"))

#Plotting the Bray-Curtis NMDS
#bray_nmds_MvF_plot <- ggplot(bray_nmds_df_MvF, aes(x = MDS1, y = MDS2, color = sex, shape = cray_group, fill = sex))+
#  stat_ellipse(aes(group = sex, fill = sex), type = "t", geom = "polygon", 
#               alpha = 0.25, linetype = "twodash") +
#  scale_shape_manual(values = c(21, 24, 23)) +
#  geom_point(size = 2.5, alpha = 0.85, color = "black") +
#  scale_fill_brewer(palette = "Accent", labels = c("Male", "Female")) +
#  scale_color_brewer(palette = "Accent", labels = c("Male", "Female")) + 
#  labs(title = "16S Bray-Curtis", color = "Sex", shape = "Crayfish Feeding Group", fill = "Sex") +
#  xlab("NMDS1") +
#  ylab("NMDS2") +
#  annotate(geom = "text", x = max(bray_nmds_df_MvF$MDS1)*-2.7, y = max(bray_nmds_df_MvF$MDS2)*0.65, 
#           label = "Richness: p = 0.999, R2 = 0.004\nSex: p = 0.441, R2 = 0.014\nCollection Date: p = 0.477, R2 = 0.027\nRichness:Sex = p = 0.830, R2 = 0.008",
#           size = 3, hjust = 0, vjust = 0.1) +
#  geom_text(data = subset(bray_nmds_df_MvF, collect_date_graph == "May 2023"), aes(-3, -2), 
#                          label = "Richness: p = 0.999, R2 = 0.004\nSex: p = 0.441, R2 = 0.014\nCollection Date: p = 0.477, R2 = 0.027\nRichness:Sex = p = 0.830, R2 = 0.008", 
#            check.overlap = T) +
#  geom_text(data = ann_text, mapping = aes(x = x, y = y, label = ann_text$lab)) +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
#  facet_wrap(~ collect_date_graph)
#print(bray_nmds_MvF_plot)



#ann_text <- data.frame(x = -3,
#                       y = -2,
#                       lab = "Richness: p = 0.999, R2 = 0.004\nSex: p = 0.441, R2 = 0.014\nCollection Date: p = 0.477, R2 = 0.027\nRichness:Sex = p = 0.830, R2 = 0.008",
#                       collect_date_graph = "May 2023", 
#                       sex = "M", 
#                       cray_group = "Group1")



##### Raup-Crick NMDS for Males vs. Females ######

#I will be using the OTU table made for Bray-Curtis in the previous step for the 3 time points 

#Preparing the NMDS
#set.seed(012219) #set seed for reproducibility
#raup_dist_MvF <- vegdist(otu_table_3t, method = "raup", na.rm = T)

#raup_dist_MvF2 <- raupcrick(otu_table_3t, null = "r1", nsim = 999)

#Checking for homogeneity of dispersion
#raup_beta_diet_MvF <- betadisper(raup_dist_MvF2, meta_cray_3t$cray_group)
#permutest(raup_beta_diet_MvF)
#TukeyHSD(raup_beta_diet_MvF)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.147 :)

#Checking for homogeneity of dispersion
#raup_beta_time_MvF <- betadisper(raup_dist_MvF2, meta_cray_3t$collection_date)
#permutest(raup_beta_time_MvF)
#TukeyHSD(raup_beta_time_MvF)
#no significant difference in dispersion between sampling dates, p = 0.541 :)

#Checking for homogeneity of dispersion
#raup_beta_MvF <- betadisper(raup_dist_MvF2, meta_cray_3t$sex)
#permutest(raup_beta_MvF)
#TukeyHSD(raup_beta_MvF)
#no significant difference in dispersion between M/F, p = 0.125 :)

#Performing the PERMANOVA
#raup_adonis_MvF <- adonis2(raup_dist_MvF2 ~ h0 * sex + collection_date + SMI, data = meta_cray_3t,
#                           permutations = 999, by = "terms")
#raup_adonis_MvF
#no significant effect of any of the variables

#raup_nmds_MvF <- metaMDS(raup_dist_MvF2, k = 3, trymax = 150, previous.best) #conduct non-metric multidimensional scaling

#raup_nmds_MvF$stress
#stress value is 0.1029734, which is fairly decent

#raup_nmds_df_MvF <- as.data.frame(raup_nmds_MvF[["points"]]) %>% #generate dataframe with NMDS points 
#  rownames_to_column("index") %>% 
#  left_join(meta_cray_3t, by = join_by("index"))

#Plotting the Raup-Crick NMDS
#raup_nmds_MvF_plot <- ggplot(raup_nmds_df_MvF, aes(x = MDS1, y = MDS2, color = sex, shape = cray_group, fill = sex))+
#  stat_ellipse(aes(group = sex, fill = sex), type = "t", geom = "polygon", 
#               alpha = 0.25, linetype = "twodash") +
#  scale_shape_manual(values = c(21, 24, 23)) +
#  geom_point(size = 2.5, alpha = 0.85, color = "black") +
#  scale_fill_brewer(palette = "Accent", labels = c("Male", "Female")) +
#  scale_color_brewer(palette = "Accent", labels = c("Male", "Female")) + 
#  labs(title = "16S Raup-Crick", color = "Sex", shape = "Crayfish Feeding Group", fill = "Sex") +
#  xlab("NMDS1") +
#  ylab("NMDS2") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
#  facet_wrap(~ collect_date_graph)
#print(raup_nmds_MvF_plot)



##### Bray-Curtis NMDS for Male Feeding Groups - 3 Time Points #####


#Creating new column for better plotting of the NMDS
meta_cray_3tM <- meta_cray_3tM %>% 
  mutate(cray_group_time = case_when(cray_group == "Group 1" & collection_date == "5/1/2023" ~ "Group 1 T1", 
                                     cray_group == "Group 1" & collection_date == "11/7/2023" ~ "Group 1 T2",
                                     cray_group == "Group 1" & collection_date == "1/23/2024" ~ "Group 1 T3",
                                     cray_group == "Group 2" & collection_date == "5/1/2023" ~ "Group 2 T1", 
                                     cray_group == "Group 2" & collection_date == "11/7/2023" ~ "Group 2 T2",
                                     cray_group == "Group 2" & collection_date == "1/23/2024" ~ "Group 2 T3"))

#Making sure sample count is the same after filtering out samples
meta_otu_inter <- intersect(meta_cray_3tM$index, otu_count$index)
otu_table_3tM <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
otu_table_3tM <- otu_table_3tM[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
otu_table_3tM <- otu_table_3tM %>% 
  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_3tM))
tax_table_3tM <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_3tM <- otu_table_3tM %>% select(all_of(intersected))

#Preparing the NMDS
set.seed(012219) #set seed for reproducibility
bray_dist_M <- vegdist(otu_table_3tM, method = "bray", na.rm = T)

#Checking for homogeneity of dispersion for feeding groups
bray_beta_diet_M <- betadisper(bray_dist_M, meta_cray_3tM$cray_group)
permutest(bray_beta_diet_M)
TukeyHSD(bray_beta_diet_M)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.744 :)

#Checking for homogeneity of dispersion for sample collection date
bray_beta_time_M <- betadisper(bray_dist_M, meta_cray_3tM$collection_date)
permutest(bray_beta_time_M)
TukeyHSD(bray_beta_time_M)
#no significant difference in dispersion between sampling dates, p = 0.665 :)

#Checking for homogeneity of dispersion for group and sample collection date
bray_beta_gtime <- betadisper(bray_dist_M, meta_cray_3tM$cray_group_time)
permutest(bray_beta_gtime)
TukeyHSD(bray_beta_gtime)
#no significant difference in dispersion between group + sampling dates, p = 0.55 :)

#Performing the PERMANOVA
bray_adonis_M <- adonis2(bray_dist_M ~ h0 * cray_group + collection_date, data = meta_cray_3tM,
                           permutations = 999, by = "terms")
bray_adonis_M
#significant effect of cray_group, collection_date, and h0:cray_group


bray_nmds_M <- metaMDS(bray_dist_M, k = 3) #conduct non-metric multidimensional scaling

bray_nmds_M$stress
#stress value is 0.09330998, which is fairly decent

#Creating a dataframe with the Bray-Curtis distances and metadata
bray_nmds_df_M <- as.data.frame(bray_nmds_M[["points"]]) %>% #generate dataframe with NMDS points 
  rownames_to_column("index") %>% 
  left_join(meta_cray_3tM, by = join_by("index"))


#Plotting the Bray-Curtis NMDS - facet wrapped according to time point
bray_nmds_M_plot <- ggplot(bray_nmds_df_M, aes(x = MDS1, y = MDS2, color = cray_group, shape = diet_crayfish, fill = cray_group))+
  stat_ellipse(aes(group = cray_group, fill = cray_group), type = "t", geom = "polygon", 
               alpha = 0.2, linetype = "twodash") +
  scale_shape_manual(values = c(21, 24), labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding")) +
  geom_point(size = 2.5, alpha = 0.85, color = "black") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_color_brewer(palette = "Pastel1") + 
  labs(title = "16S Bray-Curtis", color = "Crayfish Feeding Group", shape = "Diet in Captivity", fill = "Crayfish Feeding Group") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  facet_wrap(~ collect_date_graph)
print(bray_nmds_M_plot)


#Getting the ColorBrewer hash codes for the color palettes I want to use
brewer.pal(n = 9, "Oranges")

brewer.pal(n = 9, "Blues")

#Plotting the Bray-Curtis NMDS - all time points on one plot
bray_nmds_M_plot <- ggplot(bray_nmds_df_M, aes(x = MDS1, y = MDS2, color = cray_group_time, shape = diet_crayfish, fill = cray_group_time))+
  stat_ellipse(aes(group = cray_group_time, fill = cray_group_time), type = "t", geom = "polygon", 
               alpha = 0.2, linetype = "twodash") +
  scale_shape_manual(values = c(21, 24), labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding")) +
  geom_point(size = 2.5, alpha = 0.85, color = "black") +
  scale_fill_manual(values = c("#FEE6CE", "#FDAE6B", "#F16913", "#C6DBEF", "#6BAED6", "#2171B5")) +
  scale_color_manual(values = c("#FEE6CE", "#FDAE6B", "#F16913", "#C6DBEF", "#6BAED6", "#2171B5")) + 
  labs(color = "Crayfish Feeding Group and\nSampling Time Point", shape = "Zoo Diet", fill = "Crayfish Feeding Group and\nSampling Time Point") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate("text", x = -Inf, y = -Inf, label = "Stress = 0.09", hjust = -0.15, vjust = -1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
#facet_wrap(~ collect_date_graph)
print(bray_nmds_M_plot)

#Adding centroid positions to the NMDS plot
bray_nmds_M_plot_cent <- bray_nmds_M_plot +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = cray_group_time), size = 5, shape = 21)
print(bray_nmds_M_plot_cent)


#Plotting betadisper results

  ##For collection dates
#bray.betadisper.ds <- betadisper(bray_dist_M, meta_cray_3tM$collect_date_graph)
#par(mfrow=c(1,2))
#plot(bray.betadisper.ds, label = T, label.cex = 0.5)
#boxplot(bray.betadisper.ds)
#par(mfrow=c(1,1))

  ##For feeding groups
#bray.betadisper.ds <- betadisper(bray_dist_M, meta_cray_3tM$cray_group)
#par(mfrow=c(1,2))
#plot(bray.betadisper.ds, label = T, label.cex = 0.5)
#boxplot(bray.betadisper.ds)
#par(mfrow=c(1,1))

  ##Creating a variable that combines feeding group and time
#meta_cray_3tM$cray_time <- interaction(meta_cray_3tM$cray_group, meta_cray_3tM$collect_date_graph)

  ##For feeding groups across time
#bray.betadisper.ds <- betadisper(bray_dist_M, meta_cray_3tM$cray_time)
#par(mfrow=c(1,2))
#plot(bray.betadisper.ds, label = T, label.cex = 0.4)
#boxplot(bray.betadisper.ds)
#par(mfrow=c(1,1))

#permutest(bray.betadisper.ds)
#TukeyHSD(bray.betadisper.ds)



##### Raup-Crick NMDS for Male Feeding Groups - 3 Time Points #####


#Preparing the NMDS
set.seed(012219) #set seed for reproducibility

#Creating random sub-sample based on betadisper results
#Random sub-sample test
meta_sub_raup <- meta_cray_3tM %>% 
  group_by(cray_group) %>% 
  slice_sample(n = 26, replace = F) %>%
  ungroup()
#subsampling to 26 samples for the smallest group size (Group 2) based on betadisper results

#Making sure sample count is the same after filtering out samples
meta_otu_inter <- intersect(meta_sub_raup$index, otu_count$index)
otu_table_3tM_sub <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
otu_table_3tM_sub <- otu_table_3tM_sub[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
otu_table_3tM_sub <- otu_table_3tM_sub %>% 
  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_3tM_sub))
tax_table_2_sub <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_3tM_sub <- otu_table_3tM_sub %>% select(all_of(intersected))

#Creating the Raup-Crick distance matrix
raup_dist_M <- vegdist(otu_table_3tM_sub, method = "raup", na.rm = T)

#raup_dist_M2 <- raupcrick(otu_table_3tM, null = "r1", nsim = 999)
#saveRDS(raup_dist_M2, file = "250109_raupcrick_dist_3tM.Rds")
#raup_dist_M2 <- readRDS("250109_raupcrick_dist_3tM.Rds")

#Checking for homogeneity of dispersion
raup_beta_diet_M <- betadisper(raup_dist_M, meta_sub_raup$cray_group)
permutest(raup_beta_diet_M)
TukeyHSD(raup_beta_diet_M)
#sig. diff. between feeding groups with all of the data included (before subsampling), p = 0.001
#sig. diff. between feeding groups after equal subsampling, p = 0.001, but PERMANOVA is robust to diff. in dispersion for balanced sampling designs
#see Anderson and Walsh 2013 for reference (https://esajournals.onlinelibrary.wiley.com/doi/10.1890/12-2010.1)

boxplot(raup_beta_diet_M)
#group 2 has higher dispersion both before and after subsampling

##IGNORE THIS PART##
#Original test for subsampling
#group_broad_summary = meta_cray_3tM%>% 
#  group_by(cray_group) %>% 
#  tally()
#group_broad_summary
#group 1: n = 28, group 2: n = 26
#group 1 and group 2 only differ by two samples - I tried subsampling to 26 samples 
#but dispersion was still significant between the two groups


#Checking for homogeneity of dispersion
raup_beta_time_M <- betadisper(raup_dist_M, meta_sub_raup$collection_date)
permutest(raup_beta_time_M)
TukeyHSD(raup_beta_time_M)
#sig. difference between sample collection dates with all of the data included (before subsampling), p = 0.011
#no sig. differences between collection dates after subsampling, p = 0.146

boxplot(raup_beta_time_M)
#May 2023 has the least dispersion before subsampling

#Checking for homogeneity of dispersion for group and sample collection date
raup_beta_gtime <- betadisper(raup_dist_M, meta_sub_raup$cray_group_time)
permutest(raup_beta_gtime)
TukeyHSD(raup_beta_gtime)
#sig. difference between group and collection date with all of the data included (before subsampling), p = 
#sig. diff. between feeding groups after equal subsampling, p = 0.03, but PERMANOVA is robust to diff. in dispersion for balanced sampling designs

#Performing the PERMANOVA
raup_adonis_M <- adonis2(raup_dist_M ~ h0 + cray_group + collection_date, data = meta_sub_raup,
                           permutations = 999, by = "terms")
raup_adonis_M
#h0 and collection_date should be removed compared to Bray-Curtis because it resulted in a neg. R2 value, indicating poor fit

#Removing h0 and collection date due to negative R2 values and poor fit for these variables
raup_adonis_M <- adonis2(raup_dist_M ~ cray_group, data = meta_sub_raup,
                         permutations = 999, by = "terms")
raup_adonis_M
#feeding group significant

raup_nmds_M <- metaMDS(raup_dist_M, k = 3) #conduct non-metric multidimensional scaling

raup_nmds_M$stress
#stress = 0.0834537


#Creating a dataframe with Raup-Crick distances and metadata
raup_nmds_df_M <- as.data.frame(raup_nmds_M[["points"]]) %>% #generate dataframe with NMDS points 
  rownames_to_column("index") %>% 
  left_join(meta_cray_3tM, by = join_by("index"))


#Plotting the Raup-Crick NMDS - facet wrapped according to time point
raup_nmds_M_plot <- ggplot(raup_nmds_df_M, aes(x = MDS1, y = MDS2, color = cray_group, shape = diet_crayfish, fill = cray_group))+
  stat_ellipse(aes(group = cray_group, fill = cray_group), type = "t", geom = "polygon", 
               alpha = 0.2, linetype = "twodash") +
  scale_shape_manual(values = c(21, 24), labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding")) +
  geom_point(size = 2.5, alpha = 0.85, color = "black") +
  scale_fill_brewer(palette = "Pastel1") +
  scale_color_brewer(palette = "Pastel1") + 
  labs(title = "16S Raup-Crick", color = "Crayfish Feeding Group", shape = "Diet in Captivity", fill = "Crayfish Feeding Group") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  facet_wrap(~ collect_date_graph)
print(raup_nmds_M_plot)


#Plotting the Raup-Crick NMDS - all time points on one plot
raup_nmds_M_plot <- ggplot(raup_nmds_df_M, aes(x = MDS1, y = MDS2, color = cray_group_time, shape = diet_crayfish, fill = cray_group_time))+
  stat_ellipse(aes(group = cray_group_time, fill = cray_group_time), type = "t", geom = "polygon", 
               alpha = 0.2, linetype = "twodash") +
  scale_shape_manual(values = c(21, 24), labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding")) +
  geom_point(size = 2.5, alpha = 0.85, color = "black") +
  scale_fill_manual(values = c("#FEE6CE", "#FDAE6B", "#F16913", "#C6DBEF", "#6BAED6", "#2171B5")) +
  scale_color_manual(values = c("#FEE6CE", "#FDAE6B", "#F16913", "#C6DBEF", "#6BAED6", "#2171B5")) + 
  labs(color = "Crayfish Feeding Group and\nSampling Time Point", shape = "Zoo Diet", fill = "Crayfish Feeding Group and\nSampling Time Point") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate("text", x = -Inf, y = -Inf, label = "Stress = 0.08", hjust = -0.15, vjust = -1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
#facet_wrap(~ collect_date_graph)
print(raup_nmds_M_plot)


#Adding centroid positions to the NMDS
raup_nmds_M_plot_cent <- raup_nmds_M_plot +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = cray_group_time), size = 5, shape = 21)
print(raup_nmds_M_plot_cent)

#Creating a multipanel figure for Bray-Curtis and Raup-Crick with the centroids
bray_raup_multi <- ggarrange(bray_nmds_M_plot_cent, raup_nmds_M_plot_cent, ncol = 1, 
                             nrow = 2, common.legend = TRUE, legend = "right", 
                             labels = c("a)", "b)"))
print(bray_raup_multi)


#3D Raup-Crick NMDS for wonky result with raupcrick() matrix
#xaxis <- list(title = "NMDS1")
#yaxis <- list(title = "NMDS2")
#zaxis <- list(title = "NMDS3")

#raup_3d_nmds <- plot_ly(raup_nmds_df_M, x = ~MDS1, y = ~MDS2, z = ~MDS3,
#                        type = "scatter3d", mode = "markers",
#                        color = ~ cray_group, 
#                        size = I(80)) %>% 
#  layout(title = "16S Raup-Crick NMDS with raupcrick() matrix", scene = list(xaxis = xaxis, yaxis = yaxis, zaxis = zaxis))
#raup_3d_nmds



##### LEFSe Analysis for Nov. Sampling for Male Feeding Groups #####

#As a note, the LEFSe analysis could barely identify differential taxa between the two groups;
#it only identified 5 markers and they were all unclassified taxa. - Not sure what went wrong 
#so I proceeded with analyzing according to time point and doing a before and after feeding instead of by groups since the groups may become too similar after feeding.

#meta_cray_NovM <- filter(meta_cray_3tM, collection_date == "11/7/2023")

#meta_cray_NovM %>% 
#  group_by(cray_group) %>% 
#  tally()
#group 1: n = 10, group 2: n = 9

#Making sure sample count is the same after filtering out samples from the other time points
#meta_otu_inter_lefse_NovM <- intersect(meta_cray_NovM$index, otu_count$index)
#otu_table_NovM <- otu_count %>% filter(index %in% meta_otu_inter_lefse_NovM)

#Filtering Escherichia Shigella from the OTU table
#otu_table_NovM <- otu_table_NovM[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_NovM <- otu_table_NovM %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected_cray_lefse <- intersect(tax_table_2$otu, colnames(otu_table_NovM))
#tax_table_NovM <- tax_table_2 %>% filter(otu %in% intersected_cray_lefse)
#otu_table_2NovM <- otu_table_NovM %>% select(all_of(intersected_cray_lefse))

#Converting the OTU and taxonomy tables to a matrix for phyloseq object
#tax_table_NovM_mat <- tax_table_NovM

#tax_table_NovM_mat [, "Species"] <- NA

#tax_table_NovM_mat <- tax_table_NovM_mat %>%
#  column_to_rownames("otu") %>%
#  rename_with(str_to_title) %>%
#  subset(select = -c(Size)) %>% 
#  as.matrix()

#otu_table_NovM_mat <- as.matrix(otu_table_NovM)

#meta_cray_NovM_mat <- meta_cray_NovM %>% 
#  column_to_rownames("index")

#Creating a phyloseq object with the subsampled data - since it is required for the LEFSE function
#ASV = otu_table(otu_table_NovM_mat, taxa_are_rows = F)
#TAX = tax_table(tax_table_NovM_mat)
#metadata.phyloseq = sample_data(meta_cray_NovM_mat)

#physeq_NovM = phyloseq(ASV, TAX, metadata.phyloseq)
#physeq_NovM

#Conducting the LEFSe analysis
#lefse_NovM <- run_lefse(physeq_NovM, group = "cray_group", taxa_rank = "all", 
#                        norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2)
#lefse_NovM

#To get ColorBrewer palette hash codes
#brewer.pal(n = 6, "Pastel1")

#Plotting the LEFse as a cladogram
#lefse_clad_NovM <- plot_cladogram(lefse_NovM, color = c("#FBB4AE", "#B3CDE3"))
#print(lefse_clad_NovM)

#lefse_clad +
#  scale_fill_brewer(palette = "Accent", labels = c("NOT Crayfish Fed", "Crayfish Fed"))

#Attempting to plot the analysis as a bar chart/histogram
#lefse_hist <- plot_ef_bar(lefse_NovM)
#to show just genus level markers, you must first run the lefse at the genus level
#print(lefse_hist)

#Cleaning up the lefse plot via ggplot
#lefse_hist +
#  scale_fill_brewer(palette = "Accent", labels = c("NOT Crayfish Fed", "Crayfish Fed")) +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



##### LEFSe Analysis for Both Male Feeding Groups across 3 Time Points - before vs. after #####

#meta_cray_3tM %>% 
#  group_by(diet_crayfish) %>% 
#  tally()
#before crayfish: n = 23, after crayfish: n = 31

#Converting the OTU and taxonomy tables to a matrix for phyloseq object
#tax_table_3tM_mat <- tax_table_3tM

#tax_table_3tM_mat [, "Species"] <- NA

#tax_table_3tM_mat <- tax_table_3tM_mat %>%
#  column_to_rownames("otu") %>%
#  rename_with(str_to_title) %>%
#  subset(select = -c(Size)) %>% 
#  as.matrix()

#otu_table_3tM_mat <- as.matrix(otu_table_3tM)

#meta_cray_3tM_mat <- meta_cray_3tM %>% 
#  column_to_rownames("index")

#Creating a phyloseq object with the subsampled data - since it is required for the LEFSE function
#ASV = otu_table(otu_table_3tM_mat, taxa_are_rows = F)
#TAX = tax_table(tax_table_3tM_mat)
#metadata.phyloseq = sample_data(meta_cray_3tM_mat)

#physeq_3tM = phyloseq(ASV, TAX, metadata.phyloseq)
#physeq_3tM

#Conducting the LEFSe analysis
#lefse_3tM <- run_lefse(physeq_3tM, group = "diet_crayfish", taxa_rank = "all", 
#                        norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2)
#lefse_3tM

#To get ColorBrewer palette hash codes
#brewer.pal(n = 6, "Pastel2")

#Plotting the LEFse as a cladogram
#lefse_clad_3tM <- plot_cladogram(lefse_3tM, color = c("#B3E2CD", "#FDCDAC"), 
#                                 clade_label_level = 3)
#print(lefse_clad_3tM)

#lefse_clad_3tM +
#  scale_fill_brewer(palette = "Pastel2", labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding"))

#Rerunning the lefse at the genus level 
#lefse_3tM_gen <- run_lefse(physeq_3tM, group = "diet_crayfish", taxa_rank = "Genus", 
#                       norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2)

#Attempting to plot the analysis as a bar chart/histogram
#lefse_hist <- plot_ef_bar(lefse_3tM_gen)
#to show just genus level markers, you must first run the lefse at the genus level
#print(lefse_hist)

#Cleaning up the lefse histogram plot via ggplot
#lefse_hist +
#  scale_fill_brewer(palette = "Pastel2", labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding")) +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



##### LEFSe Analysis Comparing Each Male Feeding Group across 3 Time Points #####

#meta_cray_3tM %>% 
#  group_by(diet_crayfish) %>% 
#  tally()
#group 1: n = 28, group 2: n = 26

#Converting the OTU and taxonomy tables to a matrix for phyloseq object
#tax_table_3tM_mat <- tax_table_3tM

#tax_table_3tM_mat [, "Species"] <- NA

#tax_table_3tM_mat <- tax_table_3tM_mat %>%
#  column_to_rownames("otu") %>%
#  rename_with(str_to_title) %>%
#  subset(select = -c(Size)) %>% 
#  as.matrix()

#otu_table_3tM_mat <- as.matrix(otu_table_3tM)

#meta_cray_3tM_mat <- meta_cray_3tM %>% 
#  column_to_rownames("index")

#Creating a phyloseq object with the subsampled data - since it is required for the LEFSE function
#ASV = otu_table(otu_table_3tM_mat, taxa_are_rows = F)
#TAX = tax_table(tax_table_3tM_mat)
#metadata.phyloseq = sample_data(meta_cray_3tM_mat)

#physeq_3tM = phyloseq(ASV, TAX, metadata.phyloseq)
#physeq_3tM

#Conducting the LEFSe analysis
#lefse_3tM <- run_lefse(physeq_3tM, group = "cray_group", taxa_rank = "all", 
#                       norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2)
#lefse_3tM

#To get ColorBrewer palette hash codes
#brewer.pal(n = 6, "Pastel1")

#Plotting the LEFse as a cladogram
#lefse_clad_3tM <- plot_cladogram(lefse_3tM, color = c("#FBB4AE", "#B3CDE3"), 
#                                 clade_label_level = 3)
#print(lefse_clad_3tM)

#lefse_clad_3tM +
#  scale_fill_brewer(palette = "Pastel1", labels = c("Group 1", "Group 2"))

#Rerunning the lefse at the genus level 
#lefse_3tM_gen <- run_lefse(physeq_3tM, group = "cray_group", taxa_rank = "Genus", 
#                           norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2)

#Attempting to plot the analysis as a bar chart/histogram
#lefse_hist <- plot_ef_bar(lefse_3tM_gen)
#to show just genus level markers, you must first run the lefse at the genus level
#print(lefse_hist)

#Cleaning up the lefse histogram plot via ggplot
#lefse_hist +
#  scale_fill_brewer(palette = "Pastel1", labels = c("Group 1", "Group 2")) +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



##### Betapart Density Plot for Both Male Feeding Groups across 3 Time Points - before vs. after #####

#Filtering the metadata to separate cray fed and not cray fed groups
#meta_cray_3tM0 <- filter(meta_cray_3tM, diet_crayfish == "0")

#meta_cray_3tM1 <- filter(meta_cray_3tM, diet_crayfish == "1")

#Making sure sample count is the same after filtering out samples for not cray fed group
#meta_otu_inter_cray <- intersect(meta_cray_3tM0$index, otu_count$index)
#otu_table_cray_betapart0 <- otu_count %>% filter(index %in% meta_otu_inter_cray)

#otu_table_cray_betapart0 <- otu_table_cray_betapart0[,-2] #filtering out the OTU00001 column for E. Shigella

#otu_table_cray_betapart0 <- otu_table_cray_betapart0 %>% 
#  column_to_rownames("index")

#otu_table_cray_betapart0[otu_table_cray_betapart0 > 0] <- 1

#Making sure sample count is the same after filtering out samples for cray fed group
#meta_otu_inter_cray <- intersect(meta_cray_3tM1$index, otu_count$index)
#otu_table_cray_betapart1 <- otu_count %>% filter(index %in% meta_otu_inter_cray)

#otu_table_cray_betapart1 <- otu_table_cray_betapart1[,-2] #filtering out the OTU00001 column for E. Shigella

#otu_table_cray_betapart1 <- otu_table_cray_betapart1 %>% 
#  column_to_rownames("index")

#otu_table_cray_betapart1[otu_table_cray_betapart1 > 0] <- 1

#Converting the OTU table into a betapart object
#cray_core0 <- betapart.core(otu_table_cray_betapart0)
#cray_core1 <- betapart.core(otu_table_cray_betapart1)

#Multiple site dissimilarity - just for GP not being used for distribution plot
#cray_mult <- beta.multi(cray_core, index.family = "sorensen")

#Pairwise dissimilarity - just for GP not being used for distribution plot
#cray_pair <- beta.pair(cray_core, index.family = "sorensen")

#Resample 100 times across 20 samples for both groups
#cray_samp0 <- beta.sample(cray_core0, index.family = "sorensen", sites = 20, sample = 100)
#cray_samp1 <- beta.sample(cray_core1, index.family = "sorensen", sites = 20, sample = 100)

#Extracting the values needed for plotting the distribution of SNE, SIM, and SOR
#Specifically, I am extracting the data frame containing the values for each resample
#SNE = nestedness, SIM = simpson dissimilarity/turnover, SOR = total dissimilarity based on Sorensen
#dist0 <- cray_samp0$sampled.values
#dist1 <- cray_samp1$sampled.values

#Melting to long format for plotting
#reshape0 <- melt(dist0)
#reshape1 <- melt(dist1)

#Adding a column to indicate crayfish feeding before joining the two df together to plot using ggplot
#reshape0$feeding_group <- "BEFORE Crayfish Feeding"
#reshape1$feeding_group <- "AFTER Crayfish Feeding"

#Joining the two data frames for plotting
#dist_complete <- rbind(reshape0, reshape1)

#Plotting the distribution of beta diversity metrics across both groups
#beta_dist_plot <- ggplot(dist_complete, aes(x = value, linetype = variable, fill = feeding_group)) + 
#  geom_density(alpha = 0.5, linewidth = 0.8) +
#  scale_x_break(c(0.1, 0.80)) +
#  scale_fill_manual(values = c("#FDCDAC","#B3E2CD")) +
#  scale_linetype_manual(values = c("twodash", "dotted", "solid")) +
#  labs(fill = "Diet in Captivity", linetype = "Beta Diversity Measures") +
#  xlab("Beta Diversity") +
#  ylab("Density") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#print(beta_dist_plot)



##### Betapart Density Plot for Each Male Feeding Group across 3 Time Points #####

#Filtering the metadata to separate cray fed and not cray fed groups
#meta_cray_3tM_g1 <- filter(meta_cray_3tM, cray_group == "Group 1")

#meta_cray_3tM_g2 <- filter(meta_cray_3tM, cray_group == "Group 2")

#Making sure sample count is the same after filtering out samples for not cray fed group
#meta_otu_inter_cray <- intersect(meta_cray_3tM_g1$index, otu_count$index)
#otu_table_cray_betapartg1 <- otu_count %>% filter(index %in% meta_otu_inter_cray)

#otu_table_cray_betapartg1 <- otu_table_cray_betapartg1[,-2] #filtering out the OTU00001 column for E. Shigella

#otu_table_cray_betapartg1 <- otu_table_cray_betapartg1 %>% 
#  column_to_rownames("index")

#otu_table_cray_betapartg1[otu_table_cray_betapartg1 > 0] <- 1

#Making sure sample count is the same after filtering out samples for cray fed group
#meta_otu_inter_cray <- intersect(meta_cray_3tM_g2$index, otu_count$index)
#otu_table_cray_betapartg2 <- otu_count %>% filter(index %in% meta_otu_inter_cray)

#otu_table_cray_betapartg2 <- otu_table_cray_betapartg2[,-2] #filtering out the OTU00001 column for E. Shigella

#otu_table_cray_betapartg2 <- otu_table_cray_betapartg2 %>% 
#  column_to_rownames("index")

#otu_table_cray_betapartg2[otu_table_cray_betapartg2 > 0] <- 1

#Converting the OTU table into a betapart object
#cray_coreg1 <- betapart.core(otu_table_cray_betapartg1)
#cray_coreg2 <- betapart.core(otu_table_cray_betapartg2)

#Multiple site dissimilarity - just for GP not being used for distribution plot
#cray_mult <- beta.multi(cray_core, index.family = "sorensen")

#Pairwise dissimilarity - just for GP not being used for distribution plot
#cray_pair <- beta.pair(cray_core, index.family = "sorensen")

#Resample 100 times across 25 samples for both groups
#cray_sampg1 <- beta.sample(cray_coreg1, index.family = "sorensen", sites = 25, sample = 100)
#cray_sampg2 <- beta.sample(cray_coreg2, index.family = "sorensen", sites = 25, sample = 100)

#Extracting the values needed for plotting the distribution of SNE, SIM, and SOR
#Specifically, I am extracting the data frame containing the values for each resample
#SNE = nestedness, SIM = simpson dissimilarity/turnover, SOR = total dissimilarity based on Sorensen
#distg1 <- cray_sampg1$sampled.values
#distg2 <- cray_sampg2$sampled.values

#Melting to long format for plotting
#reshapeg1 <- melt(distg1)
#reshapeg2 <- melt(distg2)

#Adding a column to indicate crayfish feeding before joining the two df together to plot using ggplot
#reshapeg1$feeding_group <- "Group 1"
#reshapeg2$feeding_group <- "Group 2"

#Joining the two data frames for plotting
#dist_complete <- rbind(reshapeg1, reshapeg2)

#Plotting the distribution of beta diversity metrics across both groups
#beta_dist_plot <- ggplot(dist_complete, aes(x = value, linetype = variable, fill = feeding_group)) + 
#  geom_density(alpha = 0.5, linewidth = 0.8) +
#  scale_x_break(c(0.07, 0.85)) +
#  scale_fill_manual(values = c("#FBB4AE", "#B3CDE3")) +
#  scale_linetype_manual(values = c("twodash", "dotted", "solid")) +
#  labs(fill = "Diet in Captivity", linetype = "Beta Diversity Measures") +
#  xlab("Beta Diversity") +
#  ylab("Density") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#print(beta_dist_plot)



##### Plotting Richness According to 3 Time Points for Male Feeding Groups #####

#h0_3tM_plot <- ggplot(meta_cray_3tM, aes(x = collect_date_graph, y = h0, fill = cray_group)) + 
#  geom_boxplot(position = "dodge", alpha = 0.8) +
#  scale_fill_brewer(palette = "Pastel1", name = "Crayfish Feeding Group") +
#  labs(x = "", y = "Effective Number of OTUs (q = 0)") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#        legend.position = "right")
#print(h0_3tM_plot)



##### Bray-Curtis NMDS for Both Male Feeding Groups Across 3 Time Points - before vs. after #####

#Making sure sample count is the same after filtering out samples
#meta_otu_inter <- intersect(meta_cray_3tM$index, otu_count$index)
#otu_table_3tM <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
#otu_table_3tM <- otu_table_3tM[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_3tM <- otu_table_3tM %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected <- intersect(tax_table_2$otu, colnames(otu_table_3tM))
#tax_table_3tM <- tax_table_2 %>% filter(otu %in% intersected)
#otu_table_3tM <- otu_table_3tM %>% select(all_of(intersected))

#Preparing the NMDS
#set.seed(012219) #set seed for reproducibility
#bray_dist_M <- vegdist(otu_table_3tM, method = "bray", na.rm = T)

#Checking for homogeneity of dispersion
#bray_beta_diet_M <- betadisper(bray_dist_M, meta_cray_3tM$diet_crayfish)
#permutest(bray_beta_diet_M)
#TukeyHSD(bray_beta_diet_M)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.857 :)

#Performing the PERMANOVA
#bray_adonis_M_bva <- adonis2(bray_dist_M ~ diet_crayfish, data = meta_cray_3tM,
#                         permutations = 999, by = "terms")
#bray_adonis_M_bva
#significant effect of diet_crayfish

#bray_nmds_M <- metaMDS(bray_dist_M, k = 3) #conduct non-metric multidimensional scaling

#bray_nmds_M$stress
#stress value is 0.09330998, which is fairly decent

#bray_nmds_df_M_bva <- as.data.frame(bray_nmds_M[["points"]]) %>% #generate dataframe with NMDS points 
#  rownames_to_column("index") %>% 
#  left_join(meta_cray_3tM, by = join_by("index"))

#Plotting the Bray-Curtis NMDS
#bray_nmds_M_bva_plot <- ggplot(bray_nmds_df_M_bva, aes(x = MDS1, y = MDS2, color = diet_crayfish, shape = diet_crayfish, fill = diet_crayfish))+
#  stat_ellipse(aes(group = diet_crayfish, fill = diet_crayfish), type = "t", geom = "polygon", 
#               alpha = 0.25, linetype = "twodash") +
#  scale_shape_manual(values = c(21, 24), labels = c("BEFORE Crayfish Feeding (captive diet)", "AFTER Crayfish Feeding (wild diet)")) +
#  geom_point(size = 2.5, alpha = 0.85, color = "black") +
#  scale_fill_brewer(palette = "Pastel2", labels = c("BEFORE Crayfish Feeding (captive diet)", "AFTER Crayfish Feeding (wild diet)")) +
#  scale_color_brewer(palette = "Pastel2", labels = c("BEFORE Crayfish Feeding (captive diet)", "AFTER Crayfish Feeding (wild diet)")) + 
#  labs(title = "16S Bray-Curtis", color = "Diet in Captivity", shape = "Diet in Captivity", fill = "Diet in Captivity") +
#  xlab("NMDS1") +
#  ylab("NMDS2") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
#print(bray_nmds_M_bva_plot)



##### Raup-Crick NMDS for Both Male Feeding Groups Across 3 Time Points - before vs. after #####

#Preparing the NMDS
#set.seed(012219) #set seed for reproducibility
#raup_dist_M <- vegdist(otu_table_3tM, method = "raup", na.rm = T)

#raup_dist_M2 <- raupcrick(otu_table_3tM, null = "r1", nsim = 999)

#Checking for homogeneity of dispersion
#raup_beta_diet_M <- betadisper(raup_dist_M2, meta_cray_3tM$diet_crayfish)
#permutest(raup_beta_diet_M)
#TukeyHSD(raup_beta_diet_M)
#significant difference in dispersion between the crayfish feeding groups, p = 0.887


#Performing the PERMANOVA
#raup_adonis_M_bva <- adonis2(raup_dist_M2 ~ h0 * collection_date + diet_crayfish, data = meta_cray_3tM,
#                         permutations = 999, by = "terms")
#raup_adonis_M_bva
#significant difference in cray_group and collection_date with slight sig. of h0:collection_date

#raup_nmds_M <- metaMDS(raup_dist_M2, k = 3) #conduct non-metric multidimensional scaling

#raup_nmds_M$stress
#received following warning:
#Warning message:
#In metaMDS(raup_dist_M2, k = 3) :
#  stress is (nearly) zero: you may have insufficient data
#stress value is 7.47236e-05, which is odd
#stressplot(raup_nmds_M)


#raup_nmds_M <- metaMDS(raup_dist_M, k = 3, trymax = 150, previous.best) #conduct non-metric multidimensional scaling

#raup_nmds_M$stress
#stress = 0.08605129 when using vegdist distance matrix


#raup_nmds_df_M <- as.data.frame(raup_nmds_M[["points"]]) %>% #generate dataframe with NMDS points 
#  rownames_to_column("index") %>% 
#  left_join(meta_cray_3tM, by = join_by("index"))

#Plotting the Raup-Crick NMDS
#raup_nmds_M_plot <- ggplot(raup_nmds_df_M, aes(x = MDS1, y = MDS2, color = cray_group, shape = diet_crayfish, fill = cray_group))+
#  stat_ellipse(aes(group = cray_group, fill = cray_group), type = "t", geom = "polygon", 
#               alpha = 0.25, linetype = "twodash") +
#  scale_shape_manual(values = c(21, 24), labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding")) +
#  geom_point(size = 2.5, alpha = 0.85, color = "black") +
#  scale_fill_brewer(palette = "Pastel1") +
#  scale_color_brewer(palette = "Pastel1") + 
#  labs(title = "16S Raup-Crick - raupcrick() matrix", color = "Crayfish Feeding Group", shape = "Diet in Captivity", fill = "Crayfish Feeding Group") +
#  xlab("NMDS1") +
#  ylab("NMDS2") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
#  facet_wrap(~ collect_date_graph)
#print(raup_nmds_M_plot)



##### Bray-Curtis NMDS for Both Male Feeding Groups Across 2 Time Points - before vs. after #####

#Creating new column for better plotting of the NMDS
#meta_cray_2tM_bva <- meta_cray_2tM_bva %>% 
#  mutate(cray_group_time = case_when(cray_group == "Group 1" & collection_date == "5/1/2023" ~ "Group 1 Time Point 1", 
#                                     cray_group == "Group 1" & collection_date == "1/23/2024" ~ "Group 1 Time Point 2", 
#                                     cray_group == "Group 2" & collection_date == "5/1/2023" ~ "Group 2 Time Point 1", 
#                                     cray_group == "Group 2" & collection_date == "1/23/2024" ~ "Group 2 Time Point 2"))

#Making sure sample count is the same after filtering out samples
#meta_otu_inter <- intersect(meta_cray_2tM_bva$index, otu_count$index)
#otu_table_2tM_bva <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
#otu_table_2tM_bva <- otu_table_2tM_bva[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_2tM_bva <- otu_table_2tM_bva %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected <- intersect(tax_table_2$otu, colnames(otu_table_2tM_bva))
#tax_table_2tM_bva <- tax_table_2 %>% filter(otu %in% intersected)
#otu_table_2tM_bva <- otu_table_2tM_bva %>% select(all_of(intersected))

#Preparing the NMDS
#set.seed(012219) #set seed for reproducibility
#bray_dist_M_bva <- vegdist(otu_table_2tM_bva, method = "bray", na.rm = T)

#Checking for homogeneity of dispersion
#bray_beta_gtime_M <- betadisper(bray_dist_M_bva, meta_cray_2tM_bva$cray_group_time)
#permutest(bray_beta_gtime_M)
#TukeyHSD(bray_beta_gtime_M)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.761 :)

#Checking for homogeneity of dispersion
#bray_beta_diet_M <- betadisper(bray_dist_M_bva, meta_cray_2tM_bva$cray_group)
#permutest(bray_beta_diet_M)
#TukeyHSD(bray_beta_diet_M)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.618 :)

#Checking for homogeneity of dispersion
#bray_beta_date_M <- betadisper(bray_dist_M_bva, meta_cray_2tM_bva$collection_date)
#permutest(bray_beta_date_M)
#TukeyHSD(bray_beta_date_M)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.815 :)

#Performing the PERMANOVA
#bray_adonis_M_bva <- adonis2(bray_dist_M_bva ~ cray_group * collection_date, data = meta_cray_2tM_bva,
#                             permutations = 999, by = "terms")
#bray_adonis_M_bva
#no sig. effect of anything

#bray_nmds_M_bva <- metaMDS(bray_dist_M_bva, k = 3) #conduct non-metric multidimensional scaling

#bray_nmds_M_bva$stress
#stress value is 0.08785398, which is fairly decent

#bray_nmds_df_M_bva <- as.data.frame(bray_nmds_M_bva[["points"]]) %>% #generate dataframe with NMDS points 
#  rownames_to_column("index") %>% 
#  left_join(meta_cray_2tM_bva, by = join_by("index"))

#Plotting the Bray-Curtis NMDS
#bray_nmds_M_bva_plot <- ggplot(bray_nmds_df_M_bva, aes(x = MDS1, y = MDS2, color = cray_group_time, shape = collect_date_graph, fill = cray_group_time))+
#  stat_ellipse(aes(group = cray_group_time, fill = cray_group_time), type = "t", geom = "polygon", 
#               alpha = 0.25, linetype = "twodash") +
#  scale_shape_manual(values = c(21, 24)) +
#  geom_point(size = 2.5, alpha = 0.85, color = "black") +
#  scale_fill_brewer(palette = "Paired") +
#  scale_color_brewer(palette = "Paired") + 
#  labs(title = "16S Bray-Curtis", color = "Crayfish Feeding Group", shape = "Collection Date", fill = "Crayfish Feeding Group") +
#  xlab("NMDS1") +
#  ylab("NMDS2") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
#print(bray_nmds_M_bva_plot)


##### Raup-Crick NMDS for Both Male Feeding Groups Across 2 Time Points - before vs. after #####

#Preparing the NMDS
#set.seed(012219) #set seed for reproducibility
#raup_dist_M_bva <- vegdist(otu_table_2tM_bva, method = "raup", na.rm = T)

#raup_dist_M2_bva <- raupcrick(otu_table_2tM_bva, null = "r1", nsim = 999)
#saveRDS(raup_dist_M2_bva, file = "250109_raupcrick_dist_2tM_bva.Rds")
#raup_dist_M2_bva <- readRDS("250109_raupcrick_dist_2tM_bva.Rds")

#Random sub-sample test
#meta_sub_raup <- meta_cray_2tM_bva %>% 
#  group_by(cray_group_time) %>% 
#  slice_sample(n = 6, replace = F) %>%
#  ungroup()

#Making sure sample count is the same after filtering out samples
#meta_otu_inter <- intersect(meta_sub_raup$index, otu_count$index)
#otu_table_2_sub <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
#otu_table_2_sub <- otu_table_2_sub[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_2_sub <- otu_table_2_sub %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected <- intersect(tax_table_2$otu, colnames(otu_table_2_sub))
#tax_table_2_sub <- tax_table_2 %>% filter(otu %in% intersected)
#otu_table_2_sub <- otu_table_2_sub %>% select(all_of(intersected))

#raup_dist_M_bva <- vegdist(otu_table_2_sub, method = "raup", na.rm = T)

#Checking for homogeneity of dispersion
#raup_beta_gtime_M <- betadisper(raup_dist_M_bva, meta_cray_2tM_bva$cray_group_time)
#permutest(raup_beta_gtime_M)
#TukeyHSD(raup_beta_gtime_M)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.28 :)
#sig. diff. using vegdist between G2 T1 (n = 8) and G1 T2 (n = 12) with G2 T1 showing more dispersion

#boxplot(raup_beta_gtime_M)

#Checking for homogeneity of dispersion
#raup_beta_diet_M <- betadisper(raup_dist_M_bva, meta_cray_2tM_bva$cray_group)
#permutest(raup_beta_diet_M)
#TukeyHSD(raup_beta_diet_M)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.082 :)
#sig. diff. using vegdist with group 2 showing more dispersion

#boxplot(raup_beta_diet_M)

#Checking for homogeneity of dispersion
#raup_beta_date_M <- betadisper(raup_dist_M_bva, meta_cray_2tM_bva$collection_date)
#permutest(raup_beta_date_M)
#TukeyHSD(raup_beta_date_M)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.916 :)
#sig. diff. using vegdist with T1 showing more dispersion

#boxplot(raup_beta_date_M)

#Performing the PERMANOVA
#raup_adonis_M_bva <- adonis2(raup_dist_M_bva ~ cray_group * collection_date, data = meta_cray_2tM_bva,
#                             permutations = 999, by = "terms")
#raup_adonis_M_bva
#no sig. effect of anything with raupcrick and vegdist

#raup_nmds_M_bva <- metaMDS(raup_dist_M_bva, k = 3) #conduct non-metric multidimensional scaling

#raup_nmds_M_bva$stress
#stress value for raupcrick is odd again = 8.242804e-05
#stress value is 0.05522396 for vegdist, which is fairly decent

#raup_nmds_df_M_bva <- as.data.frame(raup_nmds_M_bva[["points"]]) %>% #generate dataframe with NMDS points 
#  rownames_to_column("index") %>% 
#  left_join(meta_cray_2tM_bva, by = join_by("index"))

#Filtering outliers to see if it improves the raupcrick plot
#raup_nmds_df_M_bva <- filter(raup_nmds_df_M_bva, index != "UHM810_34458" & index != "UHM818_34455")

#Plotting the Raup-Crick NMDS
#raup_nmds_M_bva_plot <- ggplot(raup_nmds_df_M_bva, aes(x = MDS1, y = MDS2, color = cray_group_time, shape = collect_date_graph, fill = cray_group_time)) +
#  stat_ellipse(aes(group = cray_group_time, fill = cray_group_time), type = "t", geom = "polygon", 
#               alpha = 0.25, linetype = "twodash") +
#  scale_shape_manual(values = c(21, 24)) +
#  geom_point(size = 2.5, alpha = 0.85, color = "black") +
#  scale_fill_brewer(palette = "Paired") +
#  scale_color_brewer(palette = "Paired") + 
#  labs(title = "16S Raup-Crick", color = "Crayfish Feeding Group", shape = "Collection Date", fill = "Crayfish Feeding Group") +
#  xlab("NMDS1") +
#  ylab("NMDS2") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
#print(raup_nmds_M_bva_plot)


##### LEfSe Analysis for Both Male Feeding Groups across 2 Time Points - before vs. after #####

#Preparing the data and making sure sample count is the same after filtering out samples for the 2tM_bva data frame
meta_otu_inter <- intersect(meta_cray_2tM_bva$index, otu_count$index)
otu_table_2tM_bva <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
otu_table_2tM_bva <- otu_table_2tM_bva[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
otu_table_2tM_bva <- otu_table_2tM_bva %>% 
  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_2tM_bva))
tax_table_2tM_bva <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_2tM_bva <- otu_table_2tM_bva %>% select(all_of(intersected))

#Checking to make sure sample numbers are as we would expect after removing the middle Nov. time point
meta_cray_2tM_bva %>% 
  group_by(collection_date) %>% 
  tally()
#before crayfish (5/1/2023): n = 14, after crayfish (1/23/2024): n = 21

#Converting the OTU and taxonomy tables to a matrix for phyloseq object
tax_table_2tM_bva_mat <- tax_table_2tM_bva

#tax_table_2tM_bva_mat [, "Species"] <- NA

#Cleaning up the taxonomy matrix to add to the phyloseq object for the analysis
tax_table_2tM_bva_mat <- tax_table_2tM_bva_mat %>%
  column_to_rownames("otu") %>%
  rename_with(str_to_title) %>%
  subset(select = -c(Size)) %>% 
  as.matrix()

#Converting the OTU table and metadata to matrices for the phyloseq object
otu_table_2tM_bva_mat <- as.matrix(otu_table_2tM_bva)

meta_cray_2tM_bva_mat <- meta_cray_2tM_bva %>% 
  column_to_rownames("index")

#Creating a phyloseq object with the data - since it is required for the LEfSe function
ASV = otu_table(otu_table_2tM_bva_mat, taxa_are_rows = F)
TAX = tax_table(tax_table_2tM_bva_mat)
metadata.phyloseq = sample_data(meta_cray_2tM_bva_mat)

physeq_2tM_bva = phyloseq(ASV, TAX, metadata.phyloseq)
physeq_2tM_bva

#Conducting the LEfSe analysis
lefse_2tM_bva <- run_lefse(physeq_2tM_bva, group = "collect_date_graph", taxa_rank = "all", 
                       norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2) 
#grouping based on time point since that represents before vs. after crayfish feeding
#most of the other settings here are the default
lefse_2tM_bva

#To get ColorBrewer palette hash codes for plotting the cladogram and histogram
brewer.pal(n = 6, "Pastel2")

#Plotting the LEfSe as a cladogram
lefse_clad_2tM_bva <- plot_cladogram(lefse_2tM_bva, color = c("#B3E2CD", "#FDCDAC"), 
                                 clade_label_level = 3)
print(lefse_clad_2tM_bva)
#There was a warning message here but the plot was still produced

#Adding labels to the plot to clarify before and after crayfish feeding instead of date
lefse_clad_2tM_bva +
  scale_fill_brewer(palette = "Pastel2", labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding"))


##Attempting to clean up the look of the cladogram to avoid text overlap##
lefse_clad_2tM_bva2 <- plot_cladogram(lefse_2tM_bva, color = c("#B3E2CD", "#FDCDAC"), 
                                    clade_label_level = 7, clade_label_font_size = 2.5, 
                                    marker_legend_param = list(ncol = 3))
print(lefse_clad_2tM_bva2) #this moves all of the taxonomic labels to the legend

#Adding labels for before and after instead of date
lefse_clad_2tM_bva2 +
  scale_fill_brewer(palette = "Pastel2", labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding"))


#Testing an alternative method to clean up the look of the cladogram - changing the default only_marker argument
lefse_clad_2tM_bva3 <- plot_cladogram(lefse_2tM_bva, color = c("#B3E2CD", "#FDCDAC"), 
                                     clade_label_level = 6, only_marker = T, 
                                     clade_label_font_size = 2.5, marker_legend_param = list(ncol = 3))
print(lefse_clad_2tM_bva3)

lefse_clad_2tM_bva3 <- lefse_clad_2tM_bva3 +
  scale_fill_brewer(palette = "Pastel2", labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding"))


#Rerunning the LEfSe at the genus level in order to cleanly plot the histogram
lefse_2tM_bva_gen <- run_lefse(physeq_2tM_bva, group = "collect_date_graph", taxa_rank = "Genus", 
                           norm = "CPM", kw_cutoff = 0.05, lda_cutoff = 2)
#As a note, the same differentially abundant genera are identified with this analysis as if you use taxa_rank = "all" like above 
#so this is just to whittle down the genera so that the histogram is not too messy

#Attempting to plot the analysis as a bar chart/histogram
lefse_hist <- plot_ef_bar(lefse_2tM_bva_gen)
#to show just genus level markers, you must first run the LEfSe at the genus level
print(lefse_hist)

#Cleaning up the LEfSe histogram plot via ggplot
lefse_hist <- lefse_hist +
  scale_fill_brewer(palette = "Pastel2", labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(lefse_hist)

#Joining the cladogram and histogram into a multipanel figure
clad_hist_multi <- ggarrange(lefse_clad_2tM_bva3, lefse_hist, ncol = 1, 
                             nrow = 2, labels = c("a)", "b)"), heights = c(1.75, 1))
print(clad_hist_multi)


##### Betapart Distribution for Both Male Feeding Groups across 2 Time Points - before vs. after #####

#This analysis is based on steps/code outlined in Baselga and Orme 2012 
#(https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00224.x)

#Filtering the metadata to separate before and after crayfish feeding
meta_cray_2tM_bva_may <- filter(meta_cray_2tM_bva, collection_date == "5/1/2023")

meta_cray_2tM_bva_jan <- filter(meta_cray_2tM_bva, collection_date == "1/23/2024")

#Making sure sample count is the same after filtering out samples for before cray feeding
meta_otu_inter_cray <- intersect(meta_cray_2tM_bva_may$index, otu_count$index)
otu_table_cray_betapart0 <- otu_count %>% filter(index %in% meta_otu_inter_cray)

otu_table_cray_betapart0 <- otu_table_cray_betapart0[,-2] #filtering out the OTU00001 column for E. Shigella

otu_table_cray_betapart0 <- otu_table_cray_betapart0 %>% 
  column_to_rownames("index")

otu_table_cray_betapart0[otu_table_cray_betapart0 > 0] <- 1

#Making sure sample count is the same after filtering out samples for after cray feeding
meta_otu_inter_cray <- intersect(meta_cray_2tM_bva_jan$index, otu_count$index)
otu_table_cray_betapart1 <- otu_count %>% filter(index %in% meta_otu_inter_cray)

otu_table_cray_betapart1 <- otu_table_cray_betapart1[,-2] #filtering out the OTU00001 column for E. Shigella

otu_table_cray_betapart1 <- otu_table_cray_betapart1 %>% 
  column_to_rownames("index")

otu_table_cray_betapart1[otu_table_cray_betapart1 > 0] <- 1

#Converting the OTU table into a betapart object
cray_core0 <- betapart.core(otu_table_cray_betapart0)
cray_core1 <- betapart.core(otu_table_cray_betapart1)

#Betapart object with all samples from both groups
#otu_table_2tM_bva_betapart <- otu_table_2tM_bva
#otu_table_2tM_bva_betapart[otu_table_2tM_bva_betapart > 0] <- 1
#cray_core <- betapart.core(otu_table_2tM_bva_betapart)

#Multiple site dissimilarity - not used for distribution plot
#cray_mult <- beta.multi(cray_core, index.family = "sorensen")

#Pairwise dissimilarity - not used for distribution plot
#cray_pair <- beta.pair(cray_core, index.family = "sorensen")

#Resample 100 times across 12 samples for both groups
set.seed(012219)  #setting seed for reproducibility
cray_samp0 <- beta.sample(cray_core0, index.family = "sorensen", sites = 12, sample = 100)
cray_samp1 <- beta.sample(cray_core1, index.family = "sorensen", sites = 12, sample = 100)

#Extracting the values needed for plotting the distribution of SNE, SIM, and SOR
#Specifically, I am extracting the data frame containing the values for each resample
#SNE = nestedness, SIM = Simpson dissimilarity/turnover, SOR = total dissimilarity based on Sorensen
dist0 <- cray_samp0$sampled.values
dist1 <- cray_samp1$sampled.values

#Melting to long format for plotting
reshape0 <- melt(dist0)
reshape1 <- melt(dist1)

#Adding a column to indicate crayfish feeding before joining the two df together to plot using ggplot
reshape0$feeding_group <- "BEFORE Crayfish Feeding"
reshape1$feeding_group <- "AFTER Crayfish Feeding"

#Joining the two data frames for plotting
dist_complete <- rbind(reshape0, reshape1)

#Changing the order to before and after for plotting, otherwise it plots alphabetically (listing after cray feeding first in the plot)
dist_complete$feeding_group <- factor(dist_complete$feeding_group, levels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding"))

#Plotting the distribution of beta diversity metrics across both groups
beta_dist_plot <- ggplot(dist_complete, aes(x = value, linetype = variable, fill = feeding_group)) + 
  geom_density(alpha = 0.5, linewidth = 0.8) +
  scale_x_break(c(0.15, 0.75)) +  #adding an x-axis break since there is a large gap between the SNE and SIM + SOR
  scale_fill_manual(values = c("#B3E2CD","#FDCDAC")) +
  scale_linetype_manual(values = c("twodash", "dotted", "solid")) +
  labs(linetype = "Beta Diversity Measures", fill = "Zoo Diet") +
  xlab("Beta Diversity") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(beta_dist_plot)

#Reordering the legend so that fill comes first
beta_dist_plot <- beta_dist_plot + guides(fill = guide_legend(order = 1))

print(beta_dist_plot)


##### Statistical Modeling - Beta Diversity #####


#1. Evaluating whether crayfish feeding across 3 time points has an impact on Bray-Curtis beta dispersion-----

#Extracting the Bray-Curtis distances generated by betadisper
#meta_3tM_bray <- cbind(meta_cray_3tM, bray_dis = bray_beta_diet_M$distances)

#Comparing an ordered beta distribution and a gamma distribution
#bc_1 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name), data = meta_3tM_bray, family = Gamma())

#bc_2 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name), data = meta_3tM_bray, family = ordbeta())
#Choosing to test an ordered beta distribution since its values are bound between 0 and 1 (like proportions or probabilities)
#like the distances produced by betadisper whereas gamma can model continuous positive values with larger numbers

#AIC(bc_1, bc_2)
#AIC for gamma is lower (-74 versus -67) but ordered beta seems more appropriate since the distances are bounded 
#so I will proceed with bc_2 but can change back to bc_1 if need be

#Determining whether to include SMI as a fixed effect or a random effect
#bc_3 <- glmmTMB(bray_dis ~ cray_group * collect_point + SMI + (1 | sample_name), data = meta_3tM_bray, family = ordbeta())

#bc_4 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | SMI), data = meta_3tM_bray, family = ordbeta())
#bc_4 produced singular fit warning - DO NOT USE

#AIC(bc_3, bc_4, bc_2)
#Including SMI as a fixed effect did not improve model fit - drop from the model and use bc_2

#Determining whether to include richness as a fixed effect
#bc_5 <- glmmTMB(bray_dis ~ cray_group * collect_point + h0 + (1 | sample_name), data = meta_3tM_bray, family = ordbeta())

#AIC(bc_5, bc_2)
#Including richness did not significantly improve fit - proceed with simpler model (bc_2)

#Testing the inclusion of additional random effects
#bc_6 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date), 
#                data = meta_3tM_bray, family = ordbeta())

#bc_7 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_3tM_bray, family = ordbeta())

#AIC(bc_6, bc_7, bc_2)
#Inclusion of additional random effects did not improve model fit - bc_2 is still the best model

#summary(bc_2)

#group_means <- emmeans(bc_2, ~ cray_group * collect_point)
#group_means

#emmeans(bc_2, pairwise ~ cray_group | collect_point)

#plot(residuals(bc_2))

#sim_res <- simulateResiduals(bc_2)
#plot(sim_res)

#Anova(bc_2, type = 2)


#2. Evaluating whether crayfish feeding across 2 time points (before vs after) has an impact on Bray-Curtis beta dispersion-----

#Renaming betadisper object from bray NMDS above
#bray_beta_diet_2tM <- betadisper(bray_dist_M_bva, meta_cray_2tM_bva$cray_group)
#If changed to collection date, similar results still occur

#Extracting the Bray-Curtis distances generated by betadisper
#meta_2tM_bva_bray <- cbind(meta_cray_2tM_bva, bray_dis = bray_beta_diet_2tM$distances)

#Comparing an ordered beta distribution and a gamma distribution
#bc_1 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name), data = meta_2tM_bva_bray, family = Gamma())

#bc_2 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name), data = meta_2tM_bva_bray, family = ordbeta())
#Choosing to test an ordered beta distribution since its values are bound between 0 and 1 (like proportions or probabilities)
#like the distances produced by betadisper whereas gamma can model continuous positive values with larger numbers

#AIC(bc_1, bc_2)


#Determining whether to include SMI as a fixed effect or a random effect
#bc_3 <- glmmTMB(bray_dis ~ cray_group * collect_point + SMI + (1 | sample_name), data = meta_2tM_bva_bray, family = ordbeta())

#bc_4 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | SMI), data = meta_2tM_bva_bray, family = ordbeta())
#bc_4 produced singular fit warning - DO NOT USE

#AIC(bc_3, bc_4, bc_2)
#Including SMI as a fixed effect did not improve model fit - drop from the model and use bc_2

#Determining whether to include richness as a fixed effect
#bc_5 <- glmmTMB(bray_dis ~ cray_group * collect_point + h0 + (1 | sample_name), data = meta_2tM_bva_bray, family = ordbeta())

#AIC(bc_5, bc_2)
#Including richness did not significantly improve fit - proceed with simpler model (bc_2)

#Testing the inclusion of additional random effects
#bc_6 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date), 
#                data = meta_2tM_bva_bray, family = ordbeta())

#bc_7 <- glmmTMB(bray_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_2tM_bva_bray, family = ordbeta())

#AIC(bc_6, bc_7, bc_2)
#Inclusion of additional random effects did not improve model fit - bc_2 is still the best model

#summary(bc_2)

#group_means <- emmeans(bc_2, ~ cray_group * collect_point)
#group_means

#emmeans(bc_2, pairwise ~ cray_group | collect_point)

#plot(residuals(bc_2))

#sim_res <- simulateResiduals(bc_2)
#plot(sim_res)

#Anova(bc_2, type = 2)



#3. Evaluating whether crayfish feeding across 3 time points has an impact on Raup-Crick beta dispersion-----

#Using output from raupcrick function - can remove 2 to make vegdist matrix
#raup_beta_diet_M <- betadisper(raup_dist_M2, meta_cray_3tM$diet_crayfish)

#Extracting the Raup-Crick distances generated by betadisper
#meta_3tM_raup <- cbind(meta_cray_3tM, raup_dis = raup_beta_diet_M$distances)

#Comparing an ordered beta distribution and a gamma distribution
#rc_1 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name), data = meta_3tM_raup, family = Gamma())

#rc_2 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name), data = meta_3tM_raup, family = ordbeta())
#Choosing to test an ordered beta distribution since its values are bound between 0 and 1 (like proportions or probabilities)
#like the distances produced by betadisper whereas gamma can model continuous positive values with larger numbers

#AIC(rc_1, rc_2)
#AIC(rc_2)
#AIC for gamma is NA because gamma distribution cannot accept 0 values which are generated with Raup-Crick 
#but ordered beta seems to work

#Determining whether to include SMI as a fixed effect or a random effect
#rc_3 <- glmmTMB(raup_dis ~ cray_group * collect_point + SMI + (1 | sample_name), data = meta_3tM_raup, family = ordbeta())

#rc_4 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | SMI), data = meta_3tM_raup, family = ordbeta())

#AIC(rc_3, rc_4, rc_2)
#Including SMI as a fixed effect did not improve model fit - drop from the model and use rc_2

#Determining whether to include richness as a fixed effect
#rc_5 <- glmmTMB(raup_dis ~ cray_group * collect_point + h0 + (1 | sample_name), data = meta_3tM_raup, family = ordbeta())

#AIC(rc_5, rc_2)
#Including richness did not significantly improve fit - proceed with simpler model (rc_2)

#Testing the inclusion of additional random effects
#rc_6 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date), 
#                data = meta_3tM_raup, family = ordbeta())

#rc_7 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_3tM_raup, family = ordbeta())

#AIC(rc_6, rc_7, rc_2)
#Inclusion of additional random effects did not improve model fit - rc_2 is still the best model

#summary(rc_2)

#group_means <- emmeans(rc_2, ~ cray_group * collect_point)
#group_means

#emmeans(rc_2, pairwise ~ cray_group | collect_point)

#plot(residuals(rc_2))

#sim_res <- simulateResiduals(rc_2)
#plot(sim_res)

#Anova(rc_2, type = 2)


#4. Evaluating whether crayfish feeding across 2 time points (before vs after) has an impact on Raup-Crick beta dispersion-----

#Renaming betadisper object from raup NMDS above
#Using output from raupcrick function - can remove 2 to make vegdist matrix
#raup_beta_diet_2tM <- betadisper(raup_dist_M2_bva, meta_cray_2tM_bva$cray_group)

#Extracting the raup-Curtis distances generated by betadisper
#meta_2tM_bva_raup <- cbind(meta_cray_2tM_bva, raup_dis = raup_beta_diet_2tM$distances)

#Comparing an ordered beta distribution and a gamma distribution
#rc_1 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name), data = meta_2tM_bva_raup, family = Gamma())

#rc_2 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name), data = meta_2tM_bva_raup, family = ordbeta())
#Choosing to test an ordered beta distribution since its values are bound between 0 and 1 (like proportions or probabilities)
#like the distances produced by betadisper whereas gamma can model continuous positive values with larger numbers

#AIC(rc_1, rc_2)
#AIC(rc_2)
#AIC for gamma is NA because gamma distribution cannot accept 0 values which are generated with Raup-Crick 
#but ordered beta seems to work

#Determining whether to include SMI as a fixed effect or a random effect
#rc_3 <- glmmTMB(raup_dis ~ cray_group * collect_point + SMI + (1 | sample_name), data = meta_2tM_bva_raup, family = ordbeta())

#rc_4 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | SMI), data = meta_2tM_bva_raup, family = ordbeta())
#rc_4 produced singular fit warning - DO NOT USE

#AIC(rc_3, rc_4, rc_2)
#Including SMI as a fixed effect did not improve model fit - drop from the model and use rc_2

#Determining whether to include richness as a fixed effect
#rc_5 <- glmmTMB(raup_dis ~ cray_group * collect_point + h0 + (1 | sample_name), data = meta_2tM_bva_raup, family = ordbeta())

#AIC(rc_5, rc_2)
#Including richness did not significantly improve fit - proceed with simpler model (rc_2)

#Testing the inclusion of additional random effects
#rc_6 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date), 
#                data = meta_2tM_bva_raup, family = ordbeta())

#rc_7 <- glmmTMB(raup_dis ~ cray_group * collect_point + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_2tM_bva_raup, family = ordbeta())

#AIC(rc_6, rc_7, rc_2)
#Inclusion of additional random effects did not improve model fit - rc_2 is still the best model

#summary(rc_2)

#group_means <- emmeans(rc_2, ~ cray_group * collect_point)
#group_means

#emmeans(rc_2, pairwise ~ cray_group | collect_point)

#plot(residuals(rc_2))

#sim_res <- simulateResiduals(rc_2)
#plot(sim_res)

#Anova(rc_2, type = 2)



##### Raup-Crick NMDS for Male Feeding Groups - 3 Time Points --> Filtering potential outliers for raupcrick function#####

#Filtering for potential outliers to see if it improves raupcrick plot
#meta_raup <- filter(meta_cray_3tM, index != "UHM810_34458" & index != "UHM818_34455" &
#                      index != "UHM832_34494")

#Making sure sample count is the same after filtering out samples
#meta_otu_inter <- intersect(meta_raup$index, otu_count$index)
#otu_table_3tM <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
#otu_table_3tM <- otu_table_3tM[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_3tM <- otu_table_3tM %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected <- intersect(tax_table_2$otu, colnames(otu_table_3tM))
#tax_table_3tM <- tax_table_2 %>% filter(otu %in% intersected)
#otu_table_3tM <- otu_table_3tM %>% select(all_of(intersected))

#raup_dist_M2 <- raupcrick(otu_table_3tM, null = "r1", nsim = 999)
#saveRDS(raup_dist_M2, file = "250109_raupcrick_dist_3tM.Rds")
#raup_dist_M2 <- readRDS("250109_raupcrick_dist_3tM.Rds")

#Checking for homogeneity of dispersion
#raup_beta_diet_M <- betadisper(raup_dist_M2, meta_raup$cray_group)
#permutest(raup_beta_diet_M)
#TukeyHSD(raup_beta_diet_M)
#significant difference in dispersion between the crayfish feeding groups, p = 0.024

#boxplot(raup_beta_diet_M)
#group 2 has slightly higher dispersion - will have to check the n = 

#group_broad_summary = meta_raup%>% 
#  group_by(cray_group) %>% 
#  tally()
#group_broad_summary
#group 1: n = 28, group 2: n = 26
#group 1 and group 2 only differ by two samples - I tried subsampling to 26 samples 
#but it did not make a difference in the dispersion between the two groups

#Checking for homogeneity of dispersion
#raup_beta_time_M <- betadisper(raup_dist_M2, meta_raup$collection_date)
#permutest(raup_beta_time_M)
#TukeyHSD(raup_beta_time_M)
#no significant difference in dispersion between sampling dates, p = 0.355 :)


#Performing the PERMANOVA
#raup_adonis_M <- adonis2(raup_dist_M2 ~ h0 * cray_group * collection_date + SMI, data = meta_raup,
#                         permutations = 999, by = "terms")
#raup_adonis_M
#significant difference in cray_group and collection_date with slight sig. of h0:collection_date

#raup_nmds_M <- metaMDS(raup_dist_M2, k = 3) #conduct non-metric multidimensional scaling

#raup_nmds_M$stress
#received following warning:
#Warning message:
#In metaMDS(raup_dist_M2, k = 3) :
#  stress is (nearly) zero: you may have insufficient data
#stress value is 7.47236e-05, which is odd
#stressplot(raup_nmds_M)


#raup_nmds_M <- metaMDS(raup_dist_M, k = 3, trymax = 150, previous.best) #conduct non-metric multidimensional scaling

#raup_nmds_M$stress
#stress = 0.08605129 when using vegdist distance matrix


#raup_nmds_df_M <- as.data.frame(raup_nmds_M[["points"]]) %>% #generate dataframe with NMDS points 
#  rownames_to_column("index") %>% 
#  left_join(meta_raup, by = join_by("index"))

#Filtering for potential outliers to see if it improves raupcrick plot
#raup_nmds_df_M <- filter(raup_nmds_df_M, index != "UHM810_34458" & index != "UHM818_34455" &
#                           index != "UHM832_34494")

#Plotting the Raup-Crick NMDS
#raup_nmds_M_plot <- ggplot(raup_nmds_df_M, aes(x = MDS1, y = MDS2, color = cray_group, shape = diet_crayfish, fill = cray_group))+
#  stat_ellipse(aes(group = cray_group, fill = cray_group), type = "t", geom = "polygon", 
#               alpha = 0.25, linetype = "twodash") +
#  scale_shape_manual(values = c(21, 24), labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding")) +
#  geom_point(size = 2.5, alpha = 0.85, color = "black") +
#  scale_fill_brewer(palette = "Pastel1") +
#  scale_color_brewer(palette = "Pastel1") + 
#  labs(title = "16S Raup-Crick - raupcrick() matrix", color = "Crayfish Feeding Group", shape = "Diet in Captivity", fill = "Crayfish Feeding Group") +
#  xlab("NMDS1") +
#  ylab("NMDS2") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
#  facet_wrap(~ collect_date_graph)
#print(raup_nmds_M_plot)




##### Playing with removing indv. that seem to have high richness #####

#Filter out individuals that have high richness
#high_rich <- c(11973, 12120, 12187, 12105, 11959, 11953, 12055, 12144, 12079, 22964)

#meta_cray_3tM_filter <- filter(meta_cray_3tM, !animal_biosample %in% high_rich)

#Using the model identified by stepwise selection that includes SMI as a fixed effect
#h0_7 <- glmmTMB(h0 ~ collect_point + cray_group + collect_point * cray_group + 
#                  SMI + (1 | sample_name) + (1 | dna_experiment_date), data = meta_cray_3tM, family = poisson())

#h0_7_filter <- glmmTMB(h0 ~ collect_point + cray_group + collect_point * cray_group + 
#                         SMI + (1 | sample_name) + (1 | dna_experiment_date), data = meta_cray_3tM_filter, family = poisson())

#AIC(h0_7, h0_7_filter)

#summary(h0_7_filter)

#plot(residuals(h0_7_filter))

#sim_res <- simulateResiduals(h0_7)
#plot(sim_res)
#testDispersion(h0_7_filter) 

#overdisp_ratio = deviance(h0_7_filter)/df.residual(h0_7_filter)
#overdisp_ratio

#group_means <- emmeans(h0_7_filter, ~ cray_group * collect_point + SMI)
#group_means

#r.squaredGLMM(h0_7_filter)




##### Bray-Curtis NMDS for Male Feeding Groups - 3 Time Points #####
#Making sure sample count is the same after filtering out samples
#meta_otu_inter <- intersect(meta_cray_3tM_filter$index, otu_count$index)
#otu_table_3tM <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
#otu_table_3tM <- otu_table_3tM[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_3tM <- otu_table_3tM %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected <- intersect(tax_table_2$otu, colnames(otu_table_3tM))
#tax_table_3tM <- tax_table_2 %>% filter(otu %in% intersected)
#otu_table_3tM <- otu_table_3tM %>% select(all_of(intersected))

#Preparing the NMDS
#set.seed(012219) #set seed for reproducibility
#bray_dist_M <- vegdist(otu_table_3tM, method = "bray", na.rm = T)

#Checking for homogeneity of dispersion
#bray_beta_diet_M <- betadisper(bray_dist_M, meta_cray_3tM_filter$cray_group)
#permutest(bray_beta_diet_M)
#TukeyHSD(bray_beta_diet_M)
#no significant difference in dispersion between the crayfish feeding groups, p = 0.744 :)

#Checking for homogeneity of dispersion
#bray_beta_time_M <- betadisper(bray_dist_M, meta_cray_3tM_filter$collection_date)
#permutest(bray_beta_time_M)
#TukeyHSD(bray_beta_time_M)
#no significant difference in dispersion between sampling dates, p = 0.665 :)

#Performing the PERMANOVA
#bray_adonis_M <- adonis2(bray_dist_M ~ h0 * cray_group + collection_date, data = meta_cray_3tM_filter,
#                         permutations = 999, by = "terms")
#bray_adonis_M



#bray_nmds_M <- metaMDS(bray_dist_M, k = 3) #conduct non-metric multidimensional scaling

#bray_nmds_M$stress

#bray_nmds_df_M <- as.data.frame(bray_nmds_M[["points"]]) %>% #generate dataframe with NMDS points 
#  rownames_to_column("index") %>% 
#  left_join(meta_cray_3tM, by = join_by("index"))

#Plotting the Bray-Curtis NMDS
#bray_nmds_M_plot <- ggplot(bray_nmds_df_M, aes(x = MDS1, y = MDS2, color = cray_group, shape = diet_crayfish, fill = cray_group))+
#  stat_ellipse(aes(group = cray_group, fill = cray_group), type = "t", geom = "polygon", 
#               alpha = 0.2, linetype = "twodash") +
#  scale_shape_manual(values = c(21, 24), labels = c("BEFORE Crayfish Feeding", "AFTER Crayfish Feeding")) +
#  geom_point(size = 2.5, alpha = 0.85, color = "black") +
#  scale_fill_brewer(palette = "Pastel1") +
#  scale_color_brewer(palette = "Pastel1") + 
#  labs(title = "16S Bray-Curtis", color = "Crayfish Feeding Group", shape = "Diet in Captivity", fill = "Crayfish Feeding Group") +
#  xlab("NMDS1") +
#  ylab("NMDS2") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
#  facet_wrap(~ collect_date_graph)
#print(bray_nmds_M_plot)



#Preparing the NMDS
#set.seed(012219) #set seed for reproducibility

#Creating random sub-sample based on betadisper results
#Random sub-sample test
#meta_sub_raup <- meta_cray_3tM %>% 
#  group_by(cray_group) %>% 
#  slice_sample(n = 26, replace = F) %>%
#  ungroup()
#subsampling to 26 samples for the smallest group size (Group 2) based on betadisper results

#Making sure sample count is the same after filtering out samples
#meta_otu_inter <- intersect(meta_sub_raup$index, otu_count$index)
#otu_table_3tM_sub <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
#otu_table_3tM_sub <- otu_table_3tM_sub[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
#otu_table_3tM_sub <- otu_table_3tM_sub %>% 
#  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
#intersected <- intersect(tax_table_2$otu, colnames(otu_table_3tM_sub))
#tax_table_2_sub <- tax_table_2 %>% filter(otu %in% intersected)
#otu_table_3tM_sub <- otu_table_3tM_sub %>% select(all_of(intersected))

#Creating the Raup-Crick distance matrix
#raup_dist_M <- vegdist(otu_table_3tM, method = "raup", na.rm = T)

#raup_dist_M2 <- raupcrick(otu_table_3tM, null = "r1", nsim = 999)
#saveRDS(raup_dist_M2, file = "250109_raupcrick_dist_3tM.Rds")
#raup_dist_M2 <- readRDS("250109_raupcrick_dist_3tM.Rds")

#Checking for homogeneity of dispersion
#raup_beta_diet_M <- betadisper(raup_dist_M, meta_cray_3tM_filter$cray_group)
#permutest(raup_beta_diet_M)
#TukeyHSD(raup_beta_diet_M)


#boxplot(raup_beta_diet_M)
#group 2 has higher dispersion

##IGNORE THIS PART##
#Original test for subsampling
#group_broad_summary = meta_cray_3tM_filter%>% 
#  group_by(cray_group) %>% 
#  tally()
#group_broad_summary
#group 1: n = 28, group 2: n = 26
#group 1 and group 2 only differ by two samples - I tried subsampling to 26 samples 
#but it did not make a difference in the dispersion between the two groups


#Checking for homogeneity of dispersion
#raup_beta_time_M <- betadisper(raup_dist_M, meta_cray_3tM_filter$collection_date)
#permutest(raup_beta_time_M)
#TukeyHSD(raup_beta_time_M)

#boxplot(raup_beta_time_M)

#Performing the PERMANOVA
#raup_adonis_M <- adonis2(raup_dist_M ~ h0 + cray_group + collection_date, data = meta_cray_3tM_filter,
#                         permutations = 999, by = "terms")
#raup_adonis_M
