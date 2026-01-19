#Researcher: Chloe Cummins
#Start Date: July 3, 2025
#Creating additional supplemental figures to include in the manuscript

#Set working directory
setwd("G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Hellbender Gut Microbiome Manuscript Stuff/Manuscript_Draft_Code_and_Figs")

#Previous code originally stored in"G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Hellbender Gut Microbiome Manuscript Stuff/Manuscript_Code_Check"
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
library(car)
library(performance)
library(lmodel2)  #for SMI calculation
library(multcompView)
library(rstatix)
library(ggsignif)
library(multcomp)
library(ggh4x)


#Importing the decontam data
otu_count <- fread("240820_HB_16S_OTU_count_table_decontam.csv")
tax_table <- fread("240820_HB_16S_taxonomy_decontam.csv")
metadata <- fread("240820_HB_16S_metadata_decontam.csv")



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
otu_table_2 <- otu_table_2 %>% dplyr::select(all_of(intersected))

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

#All of the available Hill numbers are calculated here in case they were needed for anything
#but all of my analyses mainly use richness (h0)

meta_Nash_cray$h0 <- hill_taxa(otu_table_2, q = 0, MARGIN = 1)
#q = 0 is essentially a richness metric, placing more emphasis on rare taxa
#Margin is set to the default 1 in this case since the function expects sites, or in this case samples, to be rows



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


#Creating a just male data frame for modeling 
meta_cray_feedingM <- rbind(meta_Nash_cray_g1, meta_Nash_cray_g2)
meta_cray_feedingM$cray_group <- as.factor(meta_cray_feedingM$cray_group)

#all males
meta_cray_3tM <- filter(meta_cray_feedingM, collection_date == "5/1/2023" |
                          collection_date == "11/7/2023" | collection_date == "1/23/2024")



#####Creating a boxplot comparing bacterial richness over time for crayfish feeding groups #####

#Performing the ANOVA in order to hopefully add p-value brackets to the plot
cray_anova <- aov(h0 ~ collect_date_graph * cray_group, meta_cray_3tM)
summary(cray_anova)

#Performing a TukeyHSD test since ANOVA had sig. results
cray_tukey <- TukeyHSD(x = cray_anova, conf.level = 0.95)
cray_tukey

#Plotting bacterial richness according to sample collection date and feeding group
cray_rich_date_group_plot <- ggplot(meta_cray_3tM, aes(x = collect_date_graph, y = h0, fill = cray_group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FDAE6B", "#6BAED6")) +
  labs(fill = "Crayfish Feeding Group") +
  xlab("Sample Collection Date") + 
  ylab("Effective Number of OTUs (q = 0)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
print(cray_rich_date_group_plot)


#Adding p-value brackets to the plot based on Tukey
#stat_test <- cray_tukey %>% 
#  add_xy_position(fun = "mean_sd", x = "cray_group", dodge = 0.8)

#Adding significance to the plot
#cray_rich_date_group_plot + geom_pwc(aes(group = cray_group), 
#                                     method = "tukey_hsd", label = "p.adj.signif")

#cray_rich_date_group_plot + geom_pwc(aes(group = collect_date_graph), 
#                                     method = "tukey_hsd", label = "p.adj.signif")

#test for brackets
#cray_rich_date_group_plot + geom_signif(test = "TukeyHSD", map_signif_level = T,
#                                        comparisons = list(c("Group 1", "Group 2"), c("May 2023", "Nov. 2023")))



##### Creating a plot to show trend of SMI and bacterial richness for crayfish feeding groups #####

#Creating linear model to use for regression line
cray_SMI_h0_lm <- lm(h0 ~ SMI, meta_cray_3tM)

#Plotting SMI and richness for Nashville Zoo animals
ggplot(meta_cray_3tM, aes(x = SMI, y = h0, color = cray_group)) +
  geom_point(alpha = 0.65, size = 2.5) +
  geom_abline(slope = coef(cray_SMI_h0_lm)[["SMI"]], 
              intercept = coef(cray_SMI_h0_lm)[["(Intercept)"]], color = "gray50", linewidth = 1.2) +
  #scale_shape_manual(values = c(21, 24), labels = c("Female", "Male")) +
  scale_color_manual(values = c("#FDAE6B", "#6BAED6")) +
  labs(color = "Crayfish Feeding Group") +
  xlab("Scaled Mass Index (SMI)") +
  ylab("Effective Number of OTUs (q = 0)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



##### Creating a boxplot comparing the fungal richness of the captive vs. wild mycobiome #####

#Importing the ITS data
physeq <- readRDS("physeq_ITSOTUNoDek.rds")

#Breaking the phyloseq into its separate components because I need to add a column to the metadata for plotting
#and I only need the metadata file to plot to figure
otu_mat <- as.matrix(physeq@otu_table)
tax_mat <- as.matrix(physeq@tax_table)

meta_mat <- as.matrix(physeq@sam_data)
meta_df <- as.data.frame(meta_mat)

#Adding additional information to the metadata to create groups for better plotting
meta_df <- meta_df %>% 
  mutate(group_broad = case_when(env_broad_scale == "wild" & ecoregion_III == "Interior Plateau" & bh_id_prefix == "BSC" ~ "Middle TN Wild", 
                                 env_broad_scale == "wild" & ecoregion_III == "Interior Plateau" & bh_id_prefix == "SB" ~ "Middle TN Wild",
                                 site == "Nashville Zoo" ~ "Nashville Zoo", 
                                 TRUE ~ "Middle TN Recapture"))

#Setting the group_broad variable as an ordered factor for easier plotting later on
meta_df$group_broad <- factor(meta_df$group_broad, ordered = T, levels = c("Nashville Zoo", "Middle TN Recapture", "Middle TN Wild"))

#Filtering the metadata file to remove negative controls prior to calculating alpha diversity
meta_df <- filter(meta_df, sample_name != "NTC_swab", sample_name != "Pool_neg_ctrl")

#Making sure sample count is the same after filtering out samples
otu_count <- as.data.frame(t(otu_mat))
otu_table_2 <- otu_count[-c(1,2),]

#Formatting the taxonomy table and removing any unassigned OTUs
tax_table_2 <- as.data.frame(tax_mat)
tax_table_2 <- tax_table_2 %>% 
  rownames_to_column("otu")
tax_table_2 <- filter(tax_table_2, Kingdom != "Unassigned")


#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_2))
tax_table_2 <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_2 <- otu_table_2 %>% dplyr::select(all_of(intersected))


meta_df$h0 <- hill_taxa(otu_table_2, q = 0, MARGIN = 1)
#q = 0 is essentially a richness metric, placing more emphasis on rare taxa

#Adding another column for better looking labels on the boxplot since nested labels are not easily manipulated
meta_df <- meta_df %>% 
  mutate(animal_env = case_when(env_broad_scale == "zoo" ~ "Captive", 
                                TRUE ~ "Wild–Caught"))

meta_df$animal_env <- factor(meta_df$animal_env, ordered = T, levels = c("Captive", "Wild–Caught"))

#Plotting bacterial richness according to environmental setting (zoo vs. wild) - WITH nested x-axis labels
cap_wild_rich_plot <- ggplot(meta_df, aes(interaction(group_broad, animal_env), y = h0, fill = group_broad)) +
  geom_boxplot(alpha = 0.7) +
  scale_x_discrete(guide = "axis_nested") + #to add a nested x-axis
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Site") +
  xlab("Site & Environmental Setting") + 
  ylab("Effective Number of OTUs (q = 0)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
print(cap_wild_rich_plot)


#Plotting bacterial richness according to environmental setting (zoo vs. wild) - WITHOUT nested x-axis labels
cap_wild_rich_plot_no_nest <- ggplot(meta_df, aes(x = group_broad, y = h0, fill = group_broad)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Site") +
  xlab("Site & Environmental Setting") + 
  ylab("Effective Number of OTUs (q = 0)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
print(cap_wild_rich_plot_no_nest)

#Creating a multipanel figure with the richness boxplot and fungal OTU venn diagram (must create using ITS analyses code first)
#fungi_box_venn_multi <- ggarrange(cap_wild_rich_plot_no_nest, fungi_venn_sub_plot, ncol = 1,
#                                  nrow = 2, labels = c("a)", "b)"), widths = c(0.5, 1.5), 
#                                  heights = c(0.5, 1))
#print(fungi_box_venn_multi)
