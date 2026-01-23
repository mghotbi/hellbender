#Researcher: Chloe Cummins
#Start Date: December 10, 2024
#Beginning captive vs. wild hellbender gut microbiome analyses



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
library(ggh4x)
library(DspikeIn) #this is a package created by one of our postdocs and is available on GitHub - used for core microbiome analysis
library(MuMIn)
library(MiscMetabar)
library(AICcmodavg)
library(pairwiseAdonis)

#Importing the decontam and rarefied data
otu_count <- fread("240820_HB_16S_OTU_count_table_decontam.csv")
tax_table <- fread("240820_HB_16S_taxonomy_decontam.csv")
metadata <- fread("240820_HB_16S_metadata_decontam.csv")


#####Preparing the data and creating the phyloseq object#####

#Convert the first column with sample names into row names and convert to matrices since phyloseq prefers this format
tax_mat <- tax_table %>% 
  column_to_rownames("otu") %>% 
  as.matrix()

otu_mat <- otu_count %>% 
  column_to_rownames("index") %>% 
  as.matrix()

#Metadata file needs to be kept as a dataframe to add to phyloseq
meta_df <- metadata %>% 
  column_to_rownames("index")

#Adding additional information to the metadata to create groups for better plotting
meta_df <- meta_df %>% 
  mutate(group_broad = case_when(env_broad_scale == "wild" & site == "Hiwassee River Picnic Area" ~ "East TN Wild",
                                 env_broad_scale == "environmental" & site == "Hiwassee River Picnic Area" ~ "East TN Environmental",
                                 env_broad_scale == "wild" & ecoregion_III == "Interior Plateau" & bh_id_prefix == "BSC" ~ "Middle TN Wild", 
                                 env_broad_scale == "wild" & ecoregion_III == "Interior Plateau" & bh_id_prefix == "SB" ~ "Middle TN Wild",
                                 env_broad_scale == "environmental" & ecoregion_III == "Interior Plateau" ~ "Middle TN Environmental", 
                                 site == "Chattanooga Zoo" ~ "Chattanooga Zoo",
                                 site == "Nashville Zoo" ~ "Nashville Zoo", 
                                 TRUE ~ "Middle TN Recapture"))

#Creating the phyloseq object
ASV = otu_table(otu_mat, taxa_are_rows = F)
TAX = tax_table(tax_mat)
metadata.phyloseq = sample_data(meta_df)

physeq = phyloseq(ASV, TAX, metadata.phyloseq)
physeq

#Checking the sample variables in the phyloseq object
sample_variables(physeq)

#Filtering out E.Shigella (OTU0001) before making the plots
physeq <- subset_taxa(physeq, Size != "1836290")
#as a note, I had to filter by # of reads rather than genus because there was another E. Shigella OTU other than OTU0001
physeq
physeq@tax_table@.Data

#Filtering out fecal samples and environmental samples
physeq <- subset_samples(physeq, env_medium == "cloaca swab" & env_broad_scale != "environmental")



##### Calculating the number of phyla, class, order, etc. in all HB samples #####

tax_table(physeq) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  gather("Rank", "Name", rank_names(physeq)) %>%
  na.omit() %>% # remove rows with NA value
  group_by(Rank) %>%
  summarize(ntaxa = length(unique(Name))) %>% # compute number of unique taxa
  mutate(Rank = factor(Rank, rank_names(physeq))) %>%
  arrange(Rank)



#####Creating the relative abundance plot based on Phylum - captive vs. wild all data#####

#Note: The code to create the following relative abundance plots was modified/adapted from https://rpubs.com/lgschaerer/1006964


#Transform into relative abundance of the reads
physeq <- transform_sample_counts(physeq, function(x) x/sum(x))

#Converting the phyloseq object into a data frame based on phylum rank
phy_abund <- physeq %>% tax_glom(taxrank = "phylum") %>% ##setting to the phylum level
  psmelt()  #melting data to long format for plotting
head(phy_abund)

#Filtering and modifying the data for plotting
all_data <- phy_abund %>% 
  select(phylum, Abundance, Sample, group_broad, env_broad_scale) %>% #selecting variables of interest
  filter(Abundance != 0) %>% #filtering so that abundance of taxa is at least greater than 0
  mutate(phylum = as.character(phylum))
head(all_data)

#Preparing to create relative abundance plot based on phylum
phy_plot <- all_data %>% 
  select(Sample, group_broad, env_broad_scale, phylum, Abundance) %>% #selecting all variables to be used
  group_by(Sample, group_broad, env_broad_scale) %>% 
  mutate(totalSum = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Sample, group_broad, env_broad_scale, phylum) %>% 
  summarize(Abundance = sum(Abundance), 
            totalSum, RelAb = Abundance/totalSum) %>% #calculating relative abundance based on order abundance per host
  unique() #leaving only unique observations
head(phy_plot)

#Grouping low abundance/rare taxa into Other category
phy_plot <- phy_plot %>% group_by(Sample, group_broad, env_broad_scale, phylum, totalSum) %>%
  summarise(
    Abundance = sum(Abundance),
    phylum = ifelse(RelAb < 0.03, "Other (< 3%)", phylum)) %>%
  group_by(Sample, group_broad, env_broad_scale, phylum, totalSum) %>% 
  summarize(Abundance = sum(Abundance), 
            RelAb = Abundance/totalSum) %>% #calculating relative abundance based on phylum abundance per host
  unique()
head(phy_plot)
#specifying low abundance/rare taxa (less than 3% relative abundance) as Other

#Checking to make sure that the relative abundances are as we would expect
max(phy_plot$RelAb)
mean(phy_plot$RelAb)
min(phy_plot$RelAb)
#As a note, the minimum is low because a couple of samples did not have a lot of 
#rare taxa for the Other category

length(unique(phy_plot$phylum))

#Setting up the interpolated palette for plotting
palette_length <- length(unique(phy_plot$phylum))
extend_pal <- colorRampPalette(brewer.pal(18, "Set3"))
#There is a warning here but it still creates the requested palette

#Setting up the variables as factors to allow for better plotting
phy_plot$group_broad <- factor(phy_plot$group_broad, order = T, levels = c("Nashville Zoo", "Chattanooga Zoo", 
                                                                           "East TN Wild", "Middle TN Wild", "Middle TN Recapture"))

phy_plot$env_broad_scale <- factor(phy_plot$env_broad_scale, order = T, levels = c("zoo", "wild"))

#Creating a new column based on cap/wild for better plotting
phy_plot <- phy_plot %>% 
  mutate(animal_env = case_when(env_broad_scale == "zoo" ~ "Zoo", 
                                TRUE ~ "Wild–Caught"))

phy_plot$animal_env <- factor(phy_plot$animal_env, order = T, levels = c("Zoo", "Wild–Caught"))
#converting the variable to a factor for better plotting

#for Captive labels instead of Zoo for Nash and Chatt animals
#ord_plot <- ord_plot %>% 
#  mutate(animal_env = case_when(env_broad_scale == "zoo" ~ "Captive", 
#                                TRUE ~ "Wild–Caught"))


#Plotting the relative abundance bar chart based on phylum
phy_rel_plot_all <- ggplot(phy_plot) +
  geom_col(mapping = aes(x = Sample, y = RelAb, fill = fct_relevel(phylum, c("Other (< 3%)"), after = 18)), 
           position = "stack", show.legend = TRUE, width = 1) +
  ylab("Relative Abundance") +
  xlab(NULL) +
  scale_fill_manual(values = extend_pal(palette_length)) + 
  theme_bw() +
  theme(legend.text = element_text(size = 9.25),
        legend.position = "bottom",
        legend.title = element_text(face="bold", size = 14),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(fill = "Phylum") +
  #  facet_wrap(~ group_broad) +
  #facet_grid(cols = vars(fct_relevel(group_broad, c("Nashville Zoo", "Chattanooga Zoo", "East TN Wild", "Middle TN Wild"))), 
  #           scales = "free_x", space = "free_x") +
  facet_nested(~ animal_env + group_broad, 
               scales = "free_x", space = "free_x") + #creating nested facets with the ggh4x package
  force_panelsizes(cols = c(1.32, 0.27, 0.24, 0.25, 0.31)) + 
  guides(fill = guide_legend(ncol = 8))
print(phy_rel_plot_all)

#Checking which taxa are unique amongst the groups
#ord_recap_unique <- ord_plot %>% 
#  group_by(group_broad) %>% 
#  distinct(order)
#View(ord_recap_unique)

#Calculating the SD for the phyla for each group
#phy_plot <- phy_plot %>% 
#  group_by(group_broad, phylum) %>% 
#  mutate(group_sd = sd(RelAb, na.rm = TRUE))

#Calculating the SD for overall phyla
phy_plot <- phy_plot %>% 
  group_by(phylum) %>% 
  mutate(group_sd = sd(RelAb, na.rm = TRUE))



#Changing the order of the mid TN recaps to go at the end instead of the middle
#ord_rel_plot_all +
#  facet_grid(cols = vars(fct_relevel(group_broad, c("Middle TN Recapture"), after = 4)), 
#             scales = "free_x", space = "free_x")
#print(ord_rel_plot_all)


#Note: to get the relative abundance percentages of phyla for each group, take out the "Sample" variable
#in the script above (for the all_data and phy_plot objects) so that it calculates the rel. abund. of 
#bacterial phyla based on the group rather than on a sample by sample basis (or you can remove both the sample and group variable to 
#calculate rel. abund of phyla for all hellbenders)



#####Creating the relative abundance plot based on Order - captive vs. wild all data#####

#Note: The code to create the following relative abundance plots was modified/adapted from https://rpubs.com/lgschaerer/1006964


#Transform into relative abundance of the reads
physeq <- transform_sample_counts(physeq, function(x) x/sum(x))

#Converting the phyloseq object into a data frame based on order rank
ord_abund <- physeq %>% tax_glom(taxrank = "order") %>% ##setting to the order level
  psmelt()  #melting data to long format for plotting
head(ord_abund)

#Filtering and modifying the data for plotting
all_data <- ord_abund %>% 
  select(order, Abundance, Sample, group_broad, env_broad_scale) %>% #selecting variables of interest
  filter(Abundance != 0) %>% #filtering so that abundance of taxa is at least greater than 0
  mutate(order = as.character(order))
head(all_data)

#Preparing to create relative abundance plot based on order
ord_plot <- all_data %>% 
  select(Sample, group_broad, env_broad_scale, order, Abundance) %>% #selecting all variables to be used
  group_by(Sample, group_broad, env_broad_scale) %>% 
  mutate(totalSum = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Sample, group_broad, env_broad_scale, order) %>% 
  summarize(Abundance = sum(Abundance), 
            totalSum, RelAb = Abundance/totalSum) %>% #calculating relative abundance based on order abundance per host
  unique() #leaving only unique observations
head(ord_plot)

#Grouping low abundance/rare taxa into Other category
ord_plot <- ord_plot %>% group_by(Sample, group_broad, env_broad_scale, order, totalSum) %>%
  summarise(
    Abundance = sum(Abundance),
    order = ifelse(RelAb < 0.05, "Other (< 5%)", order)) %>%
  group_by(Sample, group_broad, env_broad_scale, order, totalSum) %>% 
  summarize(Abundance = sum(Abundance), 
            RelAb = Abundance/totalSum) %>% #calculating relative abundance based on order abundance per host
  unique()
head(ord_plot)
#specifying low abundance/rare taxa (less than 5% relative abundance) as Other

#Checking to make sure that the relative abundances are as we would expect
max(ord_plot$RelAb)
mean(ord_plot$RelAb)
min(ord_plot$RelAb)
#As a note, the minimum is low because a couple of samples did not have a lot of 
#rare taxa for the Other category

length(unique(ord_plot$order))

#Setting up the interpolated palette for plotting
palette_length <- length(unique(ord_plot$order))
extend_pal <- colorRampPalette(brewer.pal(43, "Set3"))
#There is a warning here but it still creates the requested palette

#Setting up the variables as factors to allow for better plotting
ord_plot$group_broad <- factor(ord_plot$group_broad, order = T, levels = c("Nashville Zoo", "Chattanooga Zoo", 
                                                                           "East TN Wild", "Middle TN Wild", "Middle TN Recapture"))

ord_plot$env_broad_scale <- factor(ord_plot$env_broad_scale, order = T, levels = c("zoo", "wild"))

#Creating a new column based on cap/wild for better plotting
ord_plot <- ord_plot %>% 
  mutate(animal_env = case_when(env_broad_scale == "zoo" ~ "Zoo", 
                                TRUE ~ "Wild–Caught"))

ord_plot$animal_env <- factor(ord_plot$animal_env, order = T, levels = c("Zoo", "Wild–Caught"))
#converting the variable to a factor for better plotting

  #for Captive labels instead of Zoo for Nash and Chatt animals
#ord_plot <- ord_plot %>% 
#  mutate(animal_env = case_when(env_broad_scale == "zoo" ~ "Captive", 
#                                TRUE ~ "Wild–Caught"))


#Plotting the relative abundance bar chart based on order
ord_rel_plot_all <- ggplot(ord_plot) +
  geom_col(mapping = aes(x = Sample, y = RelAb, fill = fct_relevel(order, c("Other (< 5%)"), after = 43)), 
           position = "stack", show.legend = TRUE, width = 1) +
  ylab("Relative Abundance") +
  xlab(NULL) +
  scale_fill_manual(values = extend_pal(palette_length)) + 
  theme_bw() +
  theme(legend.text = element_text(size = 9.25),
        legend.position = "bottom",
        legend.title = element_text(face="bold", size = 14),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(fill = "Order") +
  #  facet_wrap(~ group_broad) +
  #facet_grid(cols = vars(fct_relevel(group_broad, c("Nashville Zoo", "Chattanooga Zoo", "East TN Wild", "Middle TN Wild"))), 
  #           scales = "free_x", space = "free_x") +
  facet_nested(~ animal_env + group_broad, 
               scales = "free_x", space = "free_x") + #creating nested facets with the ggh4x package
  force_panelsizes(cols = c(1.32, 0.27, 0.24, 0.25, 0.31)) + 
  guides(fill = guide_legend(ncol = 8))
print(ord_rel_plot_all)

#Checking which taxa are unique amongst the groups
#ord_recap_unique <- ord_plot %>% 
#  group_by(group_broad) %>% 
#  distinct(order)
#View(ord_recap_unique)


#Calculating SD for orders in each group
ord_plot <- ord_plot %>% 
  group_by(group_broad, order) %>% 
  mutate(group_sd = sd(RelAb, na.rm = TRUE))



#Changing the order of the mid TN recaps to go at the end instead of the middle
#ord_rel_plot_all +
#  facet_grid(cols = vars(fct_relevel(group_broad, c("Middle TN Recapture"), after = 4)), 
#             scales = "free_x", space = "free_x")
#print(ord_rel_plot_all)


#Note: to get the relative abundance percentages of orders for each group, take out the "Sample" variable
#in the script above (for the all_data and ord_plot objects) so that it calculates the rel. abund. of 
#bacterial orders based on the group rather than on a sample by sample basis



#####Creating the relative abundance plot based on Family - captive vs. wild all data#####

#Transform into relative abundance of the reads
#physeq <- transform_sample_counts(physeq, function(x) x/sum(x))

#Converting the phyloseq object into a data frame based on family rank
#fam_abund <- physeq %>% tax_glom(taxrank = "family") %>% ##setting to the family level
#  psmelt()  #melting data to long format for plotting
#head(fam_abund)

#Filtering and modifying the data for plotting
#all_data <- fam_abund %>% 
#  select(family, Abundance, Sample, group_broad) %>% #selecting variable of interest
#  filter(Abundance != 0) %>% #filtering so that abundance of taxa is at least greater than 0
#  mutate(family = as.character(family))
#head(all_data)

#Preparing to create relative abundance plot based on family
#fam_plot <- all_data %>% 
#  select(Sample, group_broad, family, Abundance) %>% #selecting all variables to be used
#  group_by(Sample, group_broad) %>% 
#  mutate(totalSum = sum(Abundance)) %>% 
#  ungroup() %>% 
#  group_by(Sample, group_broad, family) %>% 
#  summarize(Abundance = sum(Abundance), 
#            totalSum, RelAb = Abundance/totalSum) %>% #calculating relative abundance based on family abundance per host
#  unique() #leaving only unique observations
#head(fam_plot)

#fam_plot <- fam_plot %>% group_by(Sample, group_broad, family, totalSum) %>%
#  summarise(
#    Abundance = sum(Abundance),
#    family = ifelse(RelAb < 0.05, "Other (< 5%)", family)) %>%
#  group_by(Sample, group_broad, family, totalSum) %>% 
#  summarize(Abundance = sum(Abundance), 
#            RelAb = Abundance/totalSum) %>% #calculating relative abundance based on family abundance per host
#  unique()
#head(fam_plot)
#specifying low abundance/rare taxa (less than 5%) as Other

#Checking to make sure that the relative abundances are as we would expect
#max(fam_plot$RelAb)
#mean(fam_plot$RelAb)
#min(fam_plot$RelAb)
#As a note, the minimum is low because a couple of samples did not have a lot of 
#rare taxa for the Other category

#length(unique(fam_plot$family))

#Setting up the interpolated palette
#palette_length <- length(unique(fam_plot$family))
#extend_pal <- colorRampPalette(brewer.pal(72, "Set3"))

#Plotting the relative abundance bar chart based on family
#fam_rel_plot_all <- ggplot(fam_plot) +
#  geom_col(mapping = aes(x = Sample, y = RelAb, fill = fct_relevel(family, c("Other (< 5%)"), after = 72)), 
#           position = "stack", show.legend = TRUE, width = 1) +
#  ylab("Relative Abundance") +
#  xlab(NULL) +
#  scale_fill_manual(values = extend_pal(palette_length)) + 
#  theme_bw() +
#  theme(legend.text = element_text(size = 8),
#        legend.position = "right",
#        legend.title = element_text(face="bold", size = 14),
#        legend.key.width = unit(0.5, "cm"),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#  labs(fill = "Family") +
#  facet_wrap(~ group_broad) +
#  facet_grid(cols = vars(fct_relevel(group_broad, "Nashville Zoo")), 
#             scales = "free_x", space = "free_x") +
#  force_panelsizes(cols = c(1.32, 0.27, 0.2, 0.33, 0.25)) + 
#  guides(fill = guide_legend(ncol = 3))
#print(fam_rel_plot_all)

#Checking which taxa are unique amongst the groups
#fam_recap_unique <- fam_plot %>% 
#  group_by(group_broad) %>% 
#  distinct(family)
#View(fam_recap_unique)


#####Creating the relative abundance plot based on Genus - wild-caught only#####

#Transform into relative abundance of the reads
#physeq <- transform_sample_counts(physeq, function(x) x/sum(x))

#Filtering the metadata to only include wild-caught samples (Middle and East TN)
#physeq_wild <- subset_samples(physeq, group_broad == "Middle TN Recapture" & env_medium == "cloaca swab"| 
#                           group_broad == "Middle TN Wild" & env_medium == "cloaca swab" |
#                           group_broad == "East TN Wild")
#selecting for the pre-release samples for the 2 recaps by DNA biosample
#physeq_wild

#Converting the phyloseq object into a data frame based on genus rank
#gen_abund_wild <- physeq_wild %>% tax_glom(taxrank = "genus") %>% ##setting to the genus level
#  psmelt()  #melting data to long format for plotting
#head(gen_abund_wild)

#Filtering and modifying the data for plotting
#all_data_wild <- gen_abund_wild %>% 
#  select(genus, Abundance, Sample, group_broad) %>% #selecting variable of interest
#  filter(Abundance != 0) %>% #filtering so that abundance of taxa is at least greater than 0
#  mutate(genus = as.character(genus))
#head(all_data_wild)

#Preparing to create relative abundance plot based on genus
#gen_plot_wild <- all_data_wild %>% 
#  select(Sample, group_broad, genus, Abundance) %>% #selecting all variables to be used
#  group_by(Sample, group_broad) %>% 
#  mutate(totalSum = sum(Abundance)) %>% 
#  ungroup() %>% 
#  group_by(Sample, group_broad, genus) %>% 
#  summarize(Abundance = sum(Abundance), 
#            totalSum, RelAb = Abundance/totalSum) %>% #calculating relative abundance based on genus abundance per host
#  unique() #leaving only unique observations
#head(gen_plot_wild)

#gen_plot_wild <- gen_plot_wild %>% group_by(Sample, group_broad, genus, totalSum) %>%
#  summarise(
#    Abundance = sum(Abundance),
#    genus = ifelse(RelAb < 0.01, "Other (< 1%)", genus)) %>%
#  group_by(Sample, group_broad, genus, totalSum) %>% 
#  summarize(Abundance = sum(Abundance), 
#            RelAb = Abundance/totalSum) %>% #calculating relative abundance based on genus abundance per host
#  unique()
#head(gen_plot_wild)
#specifying low abundance/rare taxa (less than 1%) as Other

#Checking to make sure that the relative abundances are as we would expect
#max(gen_plot_wild$RelAb)
#mean(gen_plot_wild$RelAb)
#min(gen_plot_wild$RelAb)
#As a note, the minimum is low because a couple of samples did not have a lot of 
#rare taxa for the Other category

#length(unique(gen_plot_wild$genus))

#Setting up the interpolated palette
#palette_length <- length(unique(gen_plot_wild$genus))
#extend_pal <- colorRampPalette(brewer.pal(78, "Set3"))

#Plotting the relative abundance bar chart based on genus
#gen_rel_plot_wild <- ggplot(gen_plot_wild) +
#  geom_col(mapping = aes(x = Sample, y = RelAb, fill = fct_relevel(genus, c("Other (< 1%)"), after = 78)), color = "black", position = "stack", show.legend = TRUE) +
#  ylab("Relative Abundance") +
#  xlab(NULL) +
#  scale_fill_manual(values = extend_pal(palette_length)) + 
#  theme_bw() +
#  theme(legend.text = element_text(size = 8),
#        legend.position = "right",
#        legend.title = element_text(face="bold", size = 14),
#        legend.key.width = unit(0.5, "cm"),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#  labs(fill = "Genus") +
#  facet_grid(cols = vars(group_broad), scales = "free_x", space = "free_x") +
#  force_panelsizes(cols = c(1.32, 0.27, 0.2, 0.33, 0.25)) + 
#  guides(fill = guide_legend(ncol = 3))
#print(gen_rel_plot_wild)

#Checking which taxa are unique amongst the groups
#gen_recap_unique <- gen_plot %>% 
#  group_by(group_broad) %>% 
#  distinct(genus)
#View(gen_recap_unique)



#####Creating the relative abundance plot based on Genus - rewilding#####

#Transform into relative abundance of the reads
physeq <- transform_sample_counts(physeq, function(x) x/sum(x))


#Filtering the metadata to only include wild, recapture, and the pre-release samples for recaptures
physeq_rewild <- subset_samples(physeq, group_broad == "Middle TN Recapture" & env_medium == "cloaca swab"| 
                           group_broad == "Middle TN Wild" & env_medium == "cloaca swab" |
                           group_broad == "East TN Wild" | dna_biosample == "36106" |
                           dna_biosample == "36107")
#selecting for the pre-release samples for the 2 recaps by DNA biosample
physeq_rewild


#Converting the phyloseq object into a data frame based on genus rank
gen_abund_rewild <- physeq_rewild %>% tax_glom(taxrank = "genus") %>% ##setting to the genus level
  psmelt()  #melting data to long format for plotting
head(gen_abund_rewild)

#Filtering and modifying the data for plotting
all_data_rewild <- gen_abund_rewild %>% 
  select(genus, Abundance, Sample, group_broad) %>% #selecting variables of interest
  filter(Abundance != 0) %>% #filtering so that abundance of taxa is at least greater than 0
  mutate(genus = as.character(genus))
head(all_data_rewild)

#Preparing to create relative abundance plot based on genus
gen_plot_rewild <- all_data_rewild %>% 
  select(Sample, group_broad, genus, Abundance) %>% #selecting all variables to be used
  group_by(Sample, group_broad) %>% 
  mutate(totalSum = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Sample, group_broad, genus) %>% 
  summarize(Abundance = sum(Abundance), 
            totalSum, RelAb = Abundance/totalSum) %>% #calculating relative abundance based on genus abundance per host
  unique() #leaving only unique observations
head(gen_plot_rewild)

#Grouping low abundance/rare taxa into Other category
gen_plot_rewild <- gen_plot_rewild %>% group_by(Sample, group_broad, genus, totalSum) %>%
  summarise(
    Abundance = sum(Abundance),
    genus = ifelse(RelAb < 0.03, "Other (< 3%)", genus)) %>%
  group_by(Sample, group_broad, genus, totalSum) %>% 
  summarize(Abundance = sum(Abundance), 
            RelAb = Abundance/totalSum) %>% #calculating relative abundance based on genus abundance per host
  unique()
head(gen_plot_rewild)
#specifying low abundance/rare taxa (less than 3%) as Other

#Checking to make sure that the relative abundances are as we would expect
max(gen_plot_rewild$RelAb)
mean(gen_plot_rewild$RelAb)
min(gen_plot_rewild$RelAb)
#As a note, the minimum is low because a couple of samples did not have a lot of 
#rare taxa for the Other category

length(unique(gen_plot_rewild$genus))

#Setting up the interpolated palette
palette_length <- length(unique(gen_plot_rewild$genus))
extend_pal <- colorRampPalette(brewer.pal(43, "Set3"))

#Setting up the facet labels to be used
facet_labels <- as_labeller(c('Nashville Zoo' = "Pre-Release", 'Middle TN Recapture' = "Post-Release",
                              'Middle TN Wild' = "Middle TN Wild", 'East TN Wild' = "East TN Wild"))

#Plotting the relative abundance bar chart based on genus
gen_rel_plot_rewild <- ggplot(gen_plot_rewild) +
  geom_col(mapping = aes(x = Sample, y = RelAb, fill = fct_relevel(genus, c("Other (< 3%)"), after = 68)), color = "black", position = "stack", show.legend = TRUE) +
  ylab("Relative Abundance") +
  xlab(NULL) +
  scale_fill_manual(values = extend_pal(palette_length)) + 
  theme_bw() +
  theme(legend.text = element_text(size = 8),
        legend.position = "right",
        legend.title = element_text(face="bold", size = 14),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(fill = "Genus") +
  facet_grid(cols = vars(fct_relevel(group_broad, c("Nashville Zoo", "Middle TN Recapture", "Middle TN Wild", "East TN Wild"))), 
             scales = "free_x", space = "free_x", labeller = facet_labels) +
  #  force_panelsizes(cols = c(1.32, 0.27, 0.2, 0.33, 0.25)) +
  guides(fill = guide_legend(ncol = 2))
print(gen_rel_plot_rewild)


#Removing biosamples from sample name for better plotting 
gen_plot_rewild2 <- separate_wider_delim(gen_plot_rewild, cols = "Sample", delim = "_", names = c("Sample", "biosample"))

#Adding x-axis labels for samples
gen_rel_plot_rewild2 <- ggplot(gen_plot_rewild2) +
  geom_col(mapping = aes(x = Sample, y = RelAb, fill = fct_relevel(genus, c("Other (< 3%)"), after = 68)), color = "black", position = "stack", show.legend = TRUE) +
  ylab("Relative Abundance") +
  xlab(NULL) +
  scale_fill_manual(values = extend_pal(palette_length)) + 
  theme_bw() +
  theme(legend.text = element_text(size = 9.25),
        legend.position = "bottom",
        legend.title = element_text(face="bold", size = 14),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 15, vjust = 0.6, size = 7)) +
  labs(fill = "Genus") +
  facet_grid(cols = vars(fct_relevel(group_broad, c("Nashville Zoo", "Middle TN Recapture", "Middle TN Wild", "East TN Wild"))), 
             scales = "free_x", space = "free_x", labeller = facet_labels) +
  #  force_panelsizes(cols = c(1.32, 0.27, 0.2, 0.33, 0.25)) +
  guides(fill = guide_legend(ncol = 7))
print(gen_rel_plot_rewild2)

#Checking which taxa are unique amongst the groups
#gen_recap_unique <- gen_rel_plot_rewild %>% 
#  group_by(group_broad) %>% 
#  distinct(genus)
#View(gen_recap_unique)


#Note: to get the relative abundance percentages of genera for each group, take out the "Sample" variable
#in the script above (for the all_data and ord_plot objects) so that it calculates the rel. abund. of 
#bacterial genera based on the group rather than on a sample by sample basis



##### Determining how many OTUs are shared between recaps and wild residents for rewilding #####

#Subsetting the original phyloseq by broad group
ps_pre <- prune_samples(physeq_rewild@sam_data$group_broad == "Nashville Zoo", physeq_rewild) #pre-release samples
ps_recap <- prune_samples(physeq_rewild@sam_data$group_broad == "Middle TN Recapture", physeq_rewild) #recap samples
ps_mtn <- prune_samples(physeq_rewild@sam_data$group_broad == "Middle TN Wild", physeq_rewild) #Middle TN wild samples
ps_etn <- prune_samples(physeq_rewild@sam_data$group_broad == "East TN Wild", physeq_rewild) #East TN wild samples

#Filtering for OTUs that are present in the samples for each group (relative abundance greater than 0)
ps_pre <- filter_taxa(ps_pre, function(x) sum(x) > 0, prune = TRUE) #pre-release
ps_recap <- filter_taxa(ps_recap, function(x) sum(x) > 0, prune = TRUE) #recaps
ps_mtn <- filter_taxa(ps_mtn, function(x) sum(x) > 0, prune = TRUE) #Mid TN wild
ps_etn <- filter_taxa(ps_etn, function(x) sum(x) > 0, prune = TRUE) #East TN wild

#Creating the vectors of the OTUs in each group
pre_otus <- colnames(otu_table(ps_pre)) #pre-release
recap_otus <- colnames(otu_table(ps_recap)) #recaps
mtn_otus <- colnames(otu_table(ps_mtn)) #Mid TN wild
etn_otus <- colnames(otu_table(ps_etn)) #East TN wild

#Intersecting the OTU vectors to determine the shared OTUs
pre_recap <- intersect(pre_otus, recap_otus)  #between zoo pre-release samples and recaps
recap_mtn <- intersect(recap_otus, mtn_otus)  #between recaps and Mid TN wild
recap_etn <- intersect(recap_otus, etn_otus)  #between recaps and East TN wild

#Counting the number of shared OTUs between groups
length(pre_recap)
length(recap_mtn)
length(recap_etn)

#Creating tables to see which particular taxa are shared
#tax_pre_recap <- tax_table %>% filter(otu %in% pre_recap)
#tax_recap_mtn <- tax_table %>% filter(otu %in% recap_mtn)
#tax_recap_etn <- tax_table %>% filter(otu %in% recap_etn)

#Saving the shared taxa files
#write.csv(tax_pre_recap, "250202_shared_taxa_pre-release_recap.csv")
#write.csv(tax_recap_mtn, "250202_shared_taxa_recap_midTNwild.csv")
#write.csv(tax_recap_etn, "250202_shared_taxa_recap_eastTNwild.csv")

#Recoding the grouping variable for better plotting
meta_rewild <- sample_data(physeq_rewild)

meta_rewild <- data.frame(meta_rewild)

meta_rewild <- meta_rewild %>% 
  mutate(rewild_group = case_when(group_broad == "Nashville Zoo" ~ "Pre-Release", 
                                  group_broad == "Middle TN Recapture" ~ "Post-Release", 
                                  group_broad == "Middle TN Wild" ~ "Middle TN Wild", 
                                  TRUE ~ "East TN Wild"))
#Coding the groups as a factor so that they are listed in the same order as the relative abundace barplot
meta_rewild$rewild_group <- factor(meta_rewild$rewild_group, levels = c("Pre-Release", "Post-Release", 
                                                                           "Middle TN Wild", "East TN Wild"))
#Adding the recoded data back into the phyloseq object
sample_data(physeq_rewild) <- meta_rewild


#Attempting to create venn diagram for pre-release, recaps, and wild HBs
rewild_venn <- ggvenn_pq(physeq = physeq_rewild, fact = "rewild_group", set_size = 3.1, 
                         label_size = 3.1)

#Editing the venn diagram to look a little prettier
rewild_venn_plot <- rewild_venn +
  #ggplot2::scale_color_manual(c("pink3", "lightsteelblue1", "rosybrown", "blue")) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  labs(fill = "OTU count")
  #guides(fill = "none")
print(rewild_venn_plot)


#Creating a multipanel figure with the relative abundance plot and venn diagram
rewild_rel_abund_venn_multi <- ggarrange(gen_rel_plot_rewild2, rewild_venn_plot, ncol = 1, 
                                         nrow = 2, widths = c(1.1, 2), heights = c(1, 0.9),
                                         labels = c("a)", "b)"), font.label = list(size = 18))
print(rewild_rel_abund_venn_multi)

  

##### Preparing the data for NMDS plotting #####

#Converting all of the data into data frames
otu_table_2 <- as.data.frame(otu_count)

tax_table_2 <- as.data.frame(tax_table)

metadata_2 <- as.data.frame(metadata)

#Filtering the metadata to only include cloacal swab samples (no fecal or env.) for captive and wild groups
meta_all <- filter(metadata_2, env_medium == "cloaca swab")

#Adding additional information to the metadata to create groups for better plotting
meta_all <- meta_all %>% 
  mutate(group_broad = case_when(env_broad_scale == "wild" & site == "Hiwassee River Picnic Area" ~ "East TN Wild",
                                 #env_broad_scale == "environmental" & site == "Hiwassee River Picnic Area" ~ "East TN Environmental",
                                 env_broad_scale == "wild" & ecoregion_III == "Interior Plateau" & bh_id_prefix == "BSC" ~ "Middle TN Wild", 
                                 env_broad_scale == "wild" & ecoregion_III == "Interior Plateau" & bh_id_prefix == "SB" ~ "Middle TN Wild",
                                 #env_broad_scale == "environmental" & ecoregion_III == "Interior Plateau" ~ "Middle TN Environmental", 
                                 site == "Chattanooga Zoo" ~ "Chattanooga Zoo",
                                 site == "Nashville Zoo" ~ "Nashville Zoo", 
                                 TRUE ~ "Middle TN Recapture"))


#Making sure sample count is the same after filtering out samples
meta_otu_inter <- intersect(meta_all$index, otu_count$index)
otu_table_2 <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
otu_table_2 <- otu_table_2[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
otu_table_2 <- otu_table_2 %>% 
  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_2))
tax_table_2 <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_2 <- otu_table_2 %>% select(all_of(intersected))


##### Calculating Hill numbers for alpha diversity and adding to the metadata file #####

meta_all$h0 <- hill_taxa(otu_table_2, q = 0, MARGIN = 1)
#q = 0 is essentially a richness metric, placing more emphasis on rare taxa

#meta_all$h1 <- hill_taxa(otu_table_2, q = 1, MARGIN = 1)
#q = 1 is Hill-Shannon diversity, where abundant and rare taxa are equally weighted
#and evenness is accounted for

#meta_all$h2 <- hill_taxa(otu_table_2, q = 2, MARGIN = 1)
#q = 2 is inverse Simpson, which places more emphasis on common/abundant species

#meta_all$h3 <- hill_taxa(otu_table_2, q = 3, MARGIN = 1)
#q = 3 is where dominant, abundant species are given the most weight


##### Creating random sub-sample for NMDS plotting (for Raup-Crick) #####

#Note: This is an old analysis. We decided to subsample the data differently to better visualize the Nashville Zoo animals.
#Look near line 1530 for revisualization with subsampled Nashville Zoo data.

#Randomly subsampling to have equal # of samples for Raup NMDS
set.seed(050521) #set seed for reproducibility of subsampling

#Filtering for only zoo samples to subsample and add back to dataframe later 
#since the wild samples are the smallest group and all of the wild samples will need to be included
meta_sub_zoo <- filter(meta_all, group_broad == "Nashville Zoo" | group_broad == "Chattanooga Zoo")

#Subsampling only zoo data
meta_sub_zoo <- meta_sub_zoo %>% 
  group_by(group_broad) %>% 
  slice_sample(n = 12, replace = F) %>%
  ungroup()

#Filtering for wild samples to create complete subsampled data frame
meta_sub_wild <- filter(meta_all, env_broad_scale == "wild")

#Joining the zoo and wild subsets to create complete subsampled data frame for NMDS plotting
meta_sub <- rbind(meta_sub_zoo, meta_sub_wild)

#Making sure sample count is the same after filtering out samples
meta_otu_inter <- intersect(meta_sub$index, otu_count$index)
otu_table_2_sub <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
otu_table_2_sub <- otu_table_2_sub[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
otu_table_2_sub <- otu_table_2_sub %>% 
  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_2_sub))
tax_table_2_sub <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_2_sub <- otu_table_2_sub %>% select(all_of(intersected))


##### Bray-Curtis NMDS - captive vs. wild all data #####


#Calculating Bray-Curtis distances
set.seed(050521) #set seed for reproducibility
bray_dist_all <- vegdist(otu_table_2, method = "bray", na.rm = T)

#Checking for homogeneity of dispersion
#bray_beta_group <- betadisper(bray_dist_all, meta_all$group_broad)
#permutest(bray_beta_group)
#TukeyHSD(bray_beta_group)
#significant difference in dispersion based on broad group/site but ellipses are only being drawn
#based on captive vs wild so this is fine

#Checking for homogeneity of dispersion between captive and wild-caught samples
bray_beta_cw <- betadisper(bray_dist_all, meta_all$env_broad_scale)
permutest(bray_beta_cw)
TukeyHSD(bray_beta_cw)
#no significant difference in dispersion between captive vs. wild, p = 0.54

#Performing the PERMANOVA
bray_adonis_all <- adonis2(bray_dist_all ~ h0 * env_broad_scale, data = meta_all,
                         permutations = 999, by = "terms")
bray_adonis_all
#significant effect of h0, env_broad_scale, and h0:env_broad_scale without ecoregion strata

#Stratifying according to ecoregion III to account for spatial autocorrelation
bray_adonis_all_strata <- adonis2(bray_dist_all ~ h0 * env_broad_scale, data = meta_all,
                           permutations = 999, by = "terms", strata = meta_all$ecoregion_III)
bray_adonis_all_strata
#significant effect of h0, env_broad_scale, and h0:env_broad_scale with ecoregion_III strata

#As a note, I tried to stratify according to site for both the Bray-Curtis and Raup-Crick PERMANOVAs
#but the site stratification would not work for Raup-Crick (see note below), which is why I proceeded with
#ecoregion since it worked for both and to maintain consistency for the PERMANOVAs

bray_nmds_all <- metaMDS(bray_dist_all, k = 3) #conduct non-metric multidimensional scaling

bray_nmds_all$stress
#stress value is 0.1257451, which is decent

#Creating a data frame with the Bray-Curtis distances and metadata
bray_nmds_df_all <- as.data.frame(bray_nmds_all[["points"]]) %>% #generate dataframe with NMDS points 
  rownames_to_column("index") %>% 
  left_join(meta_all, by = join_by("index"))

#Plotting the Bray-Curtis NMDS
bray_nmds_all_plot <- ggplot(bray_nmds_df_all, aes(x = MDS1, y = MDS2, color = env_broad_scale, shape = group_broad, fill = env_broad_scale))+
  stat_ellipse(aes(group = env_broad_scale, fill = env_broad_scale), type = "t", geom = "polygon", 
               alpha = 0.15, linetype = "twodash") +
  scale_shape_manual(values = c(0, 2, 8, 9, 1)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_fill_manual(values = c("hotpink4","mediumpurple3")) +
  scale_color_manual(values = c("hotpink4", "mediumpurple3")) + 
  labs(title = "16S Bray-Curtis", color = "Environmental Setting", shape = "Site", fill = "Environmental Setting") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
print(bray_nmds_all_plot)

#To get hash codes for color brewer palettes
#brewer.pal(n = 8, "Dark2")

#Labelling points in NMDS to figure out which animals correspond to which points
#label <- bray_nmds_df_all$sample_name

#bray_nmds_all_plot +
#  geom_text(aes(label = label), check_overlap = TRUE)

#Adding centroid positions to the NMDS
bray_nmds_all_plot_cent <- bray_nmds_all_plot +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = env_broad_scale), size = 5, shape = 21)
print(bray_nmds_all_plot_cent)

#3D Bray-Curtis NMDS
#xaxis <- list(title = "NMDS1")
#yaxis <- list(title = "NMDS2")
#zaxis <- list(title = "NMDS3")

#bray_3d_nmds <- plot_ly(bray_nmds_df_all, x = ~MDS1, y = ~MDS2, z = ~MDS3,
#                        type = "scatter3d", mode = "markers",
#                        color = ~ env_broad_scale,
#                        colors = c("hotpink4","mediumpurple3"),
#                        #symbol = ~ group_broad, 
#                        #symbols = c("diamond", "triangle", )
#                        size = I(80)) %>% 
#  #add_markers(symbol = ~ group_broad, symbols = c("2", "8", "9", "0", "1")) %>% 
#  layout(title = "16S Bray-Curtis NMDS", scene = list(xaxis = xaxis, yaxis = yaxis, zaxis = zaxis))
#bray_3d_nmds



##### Raup-Crick NMDS - captive vs. wild subsampled data #####


#Calculating Raup-Crick distances
set.seed(050521) #set seed for reproducibility

  ###With all of the data###
#raup_dist_all <- vegdist(otu_table_2, method = "raup", na.rm = T)

#raup_dist_all2 <- raupcrick(otu_table_2, null = "r1", nsim = 999)
#saveRDS(raup_dist_all2, file = "241211_raupcrick_dist.Rds")
#raup_dist_all2 <- readRDS("241211_raupcrick_dist.Rds")

#With a subsample of the data
raup_dist_sub <- vegdist(otu_table_2_sub, method = "raup", na.rm = T)

#raup_dist_sub2 <- raupcrick(otu_table_2_sub, null = "r1", nsim = 999)
#saveRDS(raup_dist_sub2, file = "241211_raupcrick_sub_dist.Rds")
#raup_dist_sub2 <- readRDS("241211_raupcrick_sub_dist.Rds")

#Checking for homogeneity of dispersion
#raup_beta_group <- betadisper(raup_dist_sub, meta_sub$group_broad)
#permutest(raup_beta_group)
#TukeyHSD(raup_beta_group)
#no significant difference in dispersion after subsampling based on broad group/site but ellipses are only being drawn
#based on captive vs wild so this can probably be ignored for now

#Checking for homogeneity of dispersion
raup_beta_cw <- betadisper(raup_dist_sub, meta_sub$env_broad_scale)
permutest(raup_beta_cw)
TukeyHSD(raup_beta_cw)
#significant difference in dispersion between captive vs. wild using all of the data (before subsampling), p < 0.001
#after subsampling, no sig. difference between captive vs. wild, p = 0.384 

#boxplot(raup_beta_cw)
#Wild group has higher dispersion when using all of the data

#cw_summary <- meta_all %>% 
#  group_by(env_broad_scale) %>% 
#  tally()
#cw_summary
#Wild: n = 23, captive: n = 149
#Since the wild group (smaller n = ) has higher dispersion, the data should be subsampled before proceeding with PERMANOVA

#Performing the PERMANOVA
raup_adonis_sub <- adonis2(raup_dist_sub ~ env_broad_scale, data = meta_sub,
                           permutations = 999, by = "terms")
raup_adonis_sub
#env_broad_scale not significant, p = 0.058

raup_adonis_sub_strata <- adonis2(raup_dist_sub ~ env_broad_scale, data = meta_sub,
                           permutations = 999, by = "terms", strata = meta_sub$ecoregion_III)
raup_adonis_sub_strata
#I tried to stratify according to site but kept receiving a p-value of 1
#but ecoregion_III seems to work for strata - proceed with this method for now (env_broad_scale not sig., p = 0.054)

#richness (h0) was dropped from these PERMANOVAs compared to Bray-Curtis due to non-significance and potential poor fit 
#of the variable itself, where the main effect of richness resulted in a negative R2 value (-0.0447)

raup_nmds_sub <- metaMDS(raup_dist_sub, k = 3) #conduct non-metric multidimensional scaling

raup_nmds_sub$stress
#stress value is 0.1246758

raup_nmds_df_sub <- as.data.frame(raup_nmds_sub[["points"]]) %>% #generate dataframe with NMDS points 
  rownames_to_column("index") %>% 
  left_join(meta_sub, by = join_by("index"))

#Plotting the Raup-Crick NMDS
raup_nmds_sub_plot <- ggplot(raup_nmds_df_sub, aes(x = MDS1, y = MDS2, color = env_broad_scale, shape = group_broad, fill = env_broad_scale))+
  stat_ellipse(aes(group = env_broad_scale, fill = env_broad_scale), type = "t", geom = "polygon", 
               alpha = 0.15, linetype = "twodash") +
  scale_shape_manual(values = c(0, 2, 8, 9, 1)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_fill_manual(values = c("hotpink4","mediumpurple3")) +
  scale_color_manual(values = c("hotpink4", "mediumpurple3")) + 
  labs(title = "16S Raup-Crick", color = "Environmental Setting", shape = "Site", fill = "Environmental Setting") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
print(raup_nmds_sub_plot)


#Labelling points in NMDS to figure out which animals correspond to which points
#label <- raup_nmds_df_sub$sample_name

#raup_nmds_sub_plot +
#  geom_text(aes(label = label), check_overlap = TRUE)

#Adding centroid positions to the NMDS
raup_nmds_sub_plot_cent <-raup_nmds_sub_plot +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = env_broad_scale), size = 5, shape = 21)
print(raup_nmds_sub_plot_cent)

#Creating a multipanel figure of Bray and Raup with the centroids for the thesis document
bray_raup_multi <- ggarrange(bray_nmds_all_plot_cent, raup_nmds_sub_plot_cent, ncol = 1, 
                             nrow = 2, common.legend = TRUE, legend = "right")
print(bray_raup_multi)

#Creating a 3D NMDS for Raup-Crick as a supplemental figure for the thesis defense talk
#raup_3d_nmds <- plot_ly(raup_nmds_df_sub, x = ~MDS1, y = ~MDS2, z = ~MDS3,
#                        type = "scatter3d", mode = "markers",
#                        color = ~ env_broad_scale,
#                        colors = c("hotpink4","mediumpurple3"),
#                        #symbol = ~ group_broad, 
#                        #symbols = c("diamond", "triangle", )
#                        size = I(80)) %>% 
#  #add_markers(symbol = ~ group_broad, symbols = c("2", "8", "9", "0", "1")) %>% 
#  layout(title = "16S Raup-Crick NMDS", scene = list(xaxis = xaxis, yaxis = yaxis, zaxis = zaxis))
#raup_3d_nmds
  
  

##### Attempting to make a core microbiome heatmap - captive vs. wild all data #####

##### Preparing the data to create the heatmaps #####

#Convert the first column with sample names into row names and convert to matrices since phyloseq prefers this format
tax_mat_core <- tax_table %>% 
  column_to_rownames("otu") %>%
  mutate(Species = genus) %>% 
  as.matrix()

otu_mat_core <- otu_count %>% 
  column_to_rownames("index") %>% 
  as.matrix()

#Metadata file needs to be kept as a dataframe to add to phyloseq
meta_df_core <- metadata %>% 
  column_to_rownames("index")

#Adding additional information to the metadata to create groups for better plotting
meta_df_core <- meta_df_core %>% 
  mutate(group_broad = case_when(env_broad_scale == "wild" & site == "Hiwassee River Picnic Area" ~ "East TN Wild",
                                 #env_broad_scale == "environmental" & site == "Hiwassee River Picnic Area" ~ "East TN Environmental",
                                 env_broad_scale == "wild" & ecoregion_III == "Interior Plateau" & bh_id_prefix == "BSC" ~ "Middle TN Wild", 
                                 env_broad_scale == "wild" & ecoregion_III == "Interior Plateau" & bh_id_prefix == "SB" ~ "Middle TN Wild",
                                 #env_broad_scale == "environmental" & ecoregion_III == "Interior Plateau" ~ "Middle TN Environmental", 
                                 site == "Chattanooga Zoo" ~ "Chattanooga Zoo",
                                 site == "Nashville Zoo" ~ "Nashville Zoo", 
                                 TRUE ~ "Middle TN Recapture"))

#Creating the phyloseq object
ASV = otu_table(otu_mat_core, taxa_are_rows = F)
TAX = tax_table(tax_mat_core)
metadata.phyloseq = sample_data(meta_df_core)

physeq_core = phyloseq(ASV, TAX, metadata.phyloseq)
physeq_core

#Filtering out E.Shigella (OTU0001) before making the plots
physeq_core <- subset_taxa(physeq_core, Size != "1836290")
#as a note, I had to filter by # of reads rather than genus because there was another E. Shigella OTU other than OTU0001
physeq_core
physeq_core@tax_table@.Data

#Converting to relative abundance
physeq_core <- transform_sample_counts(physeq_core, function(x) x/sum(x))


##### Creating a core microbiome heatmap for wild individuals (both East and Middle TN - including recaps) #####

#Filtering out captive samples, fecal samples, and environmental samples
physeq_core_wild <- subset_samples(physeq_core, env_medium == "cloaca swab" & env_broad_scale == "wild")

#Plotting the core microbiome for wild HBs as a heatmap via Mitra's DspikeIn package
core_heatmap_wild <- plot_core_microbiome_custom(physeq_core_wild, taxrank = "family", 
#                                                output_core_csv = "241212_16S_wild_core_micro_relab.csv"
                                                )
print(core_heatmap_wild)
#The initial plot that is created is kind of busy and tries to split things into separate ASVs, but we just want to 
#see the core microbiome grouped by the family level here for interpretability, so I will try to clean up the plot a little

#Taking out the data frame produced by the function in order to get rid of the ASV prefixes on the taxa labels
core_wild_data <- core_heatmap_wild[["data"]]

#Removing the ASV prefixes (especially since my taxa represent OTUs)
core_wild_data <- separate_wider_delim(core_wild_data, cols = "Taxa", delim = "__", names = c("prefix", "Taxa"))
core_wild_data <- select(core_wild_data, -prefix)

#Replacing with the new data
core_heatmap_wild[["data"]] <- core_wild_data
print(core_heatmap_wild)
#This fixed the issue where it split families into separate ASVs and collapsed according to family
#But now, it plots according to alphabetical order instead of descending order based on prevalence

#Restructuring the plot so that taxa are arranged based on decreasing prevalence
core_heatmap_wild +
  aes(x = DetectionThreshold, y = reorder(Taxa, Prevalence), fill = Prevalence) +
  labs(x = "Detection Threshold (Relative Abundance)")


##### Creating a core microbiome heatmap for captive individuals (both Nash and Chatt Zoo) #####

#Filtering out fecal samples and environmental samples
physeq_core_captive <- subset_samples(physeq_core, env_medium == "cloaca swab" & env_broad_scale == "zoo")

#Plotting the core microbiome for captive HBs as a heatmap via Mitra's DspikeIn package
core_heatmap_captive <- plot_core_microbiome_custom(physeq_core_captive, taxrank = "family", 
                                                 #output_core_csv = "241212_16S_captive_core_micro_relab.csv"
                                                 )
print(core_heatmap_captive)

#Taking out the data frame produced by the function in order to get rid of the ASV prefixes on the taxa labels
core_captive_data <- core_heatmap_captive[["data"]]

#Removing the ASV prefixes (especially since my taxa represent OTUs)
core_captive_data <- separate_wider_delim(core_captive_data, cols = "Taxa", delim = "__", names = c("prefix", "Taxa"))
core_captive_data <- select(core_captive_data, -prefix)

#Replacing with the new data
core_heatmap_captive[["data"]] <- core_captive_data
print(core_heatmap_captive)
#This fixed the issue where it split families into separate ASVs and collapsed according to family
#But now, it plots according to alphabetical order instead of descending order based on prevalence

#Restructuring the plot so that taxa are arranged based on decreasing prevalence
core_heatmap_captive +
  aes(x = DetectionThreshold, y = reorder(Taxa, Prevalence), fill = Prevalence) +
  labs(x = "Detection Threshold (Relative Abundance)")



##### Statistical Modeling #####

#1. Evaluating whether richness is affected by captive/wild setting----

#Assessing normality of richness (response variable) among the captive and wild hellbenders
shapiro.test(meta_all$h0) #p-value < 2.2e-16, not normally distributed
#Data is not normally distributed, indicating that a GLMM may produce a better fit than an LMM

#Looking at the distribution of the response variable to double check that the stats. match
densityplot(meta_all$h0, main = "Distribution of Richness")
#right skew

#Histogram of richness
hist(meta_all$h0)

#Q-Q plot of richness
qqnorm(meta_all$h0)
qqline(meta_all$h0, col = "red")
#there seem to be some points near the top of the graph

#Q-Q plot to see which group is causing the skew in the previous plots
h0_qq <- ggplot(meta_all, aes(sample = h0, color = env_broad_scale)) + 
  stat_qq() +
  stat_qq_line() +
  theme_classic()
print(h0_qq)
#high points seem to be in both groups


## GLMM of diet's effect on richness ##

#Evaluating whether I should use a tweedie distribution or a Poisson distribution
#h0_1 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name), data = meta_all, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
h0_2 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name), data = meta_all, family = poisson())

#AIC(h0_1, h0_2)
#AICc(h0_1)
AICc(h0_2)
#summary(h0_1)
#the tweedie distribution has a lower AIC and AICc

#plot(residuals(h0_1))

#qqnorm(residuals(h0_1))
#qqline(residuals(h0_1))

#sim_res <- simulateResiduals(h0_1)
#plot(sim_res)

#I originally tested this tweedie distribution based on a suggestion by one of our postdocs
#but based on what I understand tweedie is better suited to continuous data that may have lots of zeros
#and may not be the most appropriate in this case since richness consists of discrete counts and is not zero-inflated for my data

#So I will proceed with Poisson for now and just keep tweedie in the code to avoid interfering with the numbering of the other models

#Evaluating the addition of more random effects to Poisson dist. according to potential importance
h0_3 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date),
                data = meta_all, family = poisson())

h0_4 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date) +
                  (1 | plate_ID), data = meta_all, family = poisson())

h0_5 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name) + 
                  (1 + dna_experiment_date | plate_ID), data = meta_all, family = poisson())
#h0_5 gave a model convergence error - do not use

h0_6 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date) +
                  (1 | plate_ID)  + (1 | total_length), data = meta_all, family = poisson())

AIC(h0_3, h0_4, h0_6)
#There is a warning message here because there are a couple of zoo animals that do not have complete weight/length data
#so I think that adding total length as a random effect may have dropped a couple of observations
AICc(h0_3)
AICc(h0_4)
AICc(h0_6)
#h0_6 has the lowest AIC and AICc

#Testing overdispersion
overdisp_ratio = deviance(h0_6)/df.residual(h0_6)
overdisp_ratio

summary(h0_6)

group_means <- emmeans(h0_6, ~ env_broad_scale)
group_means

emmeans(h0_6, pairwise ~ env_broad_scale)

plot(residuals(h0_6))

qqnorm(residuals(h0_6))
qqline(residuals(h0_6))
#residuals do not look great

sim_res <- simulateResiduals(h0_6)
plot(sim_res)
#KS test and dispersion test are significant while outlier test is non-sig. 

Anova(h0_6, type = 2)
#This model does not seem to fit well and appears to potentially have some issues with dispersion 
#so I will try a neg. binomial distribution


h0_7 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date), data = meta_all, family = nbinom1())

h0_8 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date), data = meta_all, family = nbinom2())

AIC(h0_7, h0_8)
AICc(h0_7)
AICc(h0_8)
#AIC and AICc lower for nbinom2 for h0_8
#Both of the negative binomial models have a lower AIC and AICc than the Poisson models

#Adding additional random and fixed effects to the model
h0_9 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date) + (1 | total_length), data = meta_all, family = nbinom2())

h0_10 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date) + (1 | total_length) + 
                  (1 | collection_date), data = meta_all, family = nbinom2())
#adding collection date as another potential random effect

#Testing total length as a fixed effect rather than random just to see how it might affect fit
h0_11 <- glmmTMB(h0 ~ env_broad_scale + total_length + (1 | sample_name) + (1 | dna_experiment_date) + 
                   (1 | collection_date), data = meta_all, family = nbinom2())

AIC(h0_9, h0_10, h0_11)
AICc(h0_9)
AICc(h0_10)
AICc(h0_11)

#These next few lines of code can be ignored because I had originally tested the model that included total length as a 
#random effect but I realized that it is probably more appropriate to include it as a fixed effect since it is continuous

summary(h0_10)

group_means <- emmeans(h0_10, ~ env_broad_scale)
group_means

emmeans(h0_10, pairwise ~ env_broad_scale)

plot(residuals(h0_10))

sim_res <- simulateResiduals(h0_10)
plot(sim_res)
#Improved model fit compared to Poisson but still problems with simulated residuals

Anova(h0_10, type = 2)

#Trying again with the tweedie to see if it helps, even though I still don't think it is the best option
#h0_12 <- glmmTMB(h0 ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date) + 
#                  (1 | collection_date), data = meta_all, family = tweedie(), 
#                control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

#AIC(h0_12)

#summary(h0_12)

#plot(residuals(h0_12))

#sim_res <- simulateResiduals(h0_12)
#plot(sim_res)
#Higher AIC and worse residuals than negative binomial


#Trying stepwise selection to see if that provides any clarity
h0_1_step <- buildglmmTMB(h0 ~ env_broad_scale + total_length + (1 | sample_name) + 
                            (1 + dna_experiment_date | plate_ID) + (1 | dna_experiment_date) + (1 | plate_ID) +
                            (1 | collection_date) + (1 | total_length), 
                          data = meta_all, family = nbinom2(), buildmerControl = 
                            buildmerControl(crit = "AIC"))
(f <- formula(h0_1_step@model))

#Testing model picked by stepwise selection
h0_13 <- glmmTMB(h0 ~ env_broad_scale + (1 | collection_date) + (1 | sample_name) + 
  (1 | total_length), data = meta_all, family = nbinom2())
#Again, I later realized that total length is more likely appropriate as a fixed effect rather than random
#which is why the model with total length as a fixed effect is tacked on at the end of this code

AIC(h0_10, h0_13)
AICc(h0_10)
AICc(h0_13)

overdisp_ratio = deviance(h0_13)/df.residual(h0_13)
overdisp_ratio

plot(residuals(h0_13))

qqnorm(residuals(h0_13))
qqline(residuals(h0_13))

sim_res <- simulateResiduals(h0_13)
plot(sim_res)
#fit is not perfect but this is as good as I have been able to get it
#As a note, I tried all of the different link functions available as well and none of them improved the model

summary(h0_13)

Anova(h0_13, type = 2)

group_means <- emmeans(h0_13, ~ env_broad_scale)
group_means

emmeans(h0_13, pairwise ~ env_broad_scale)

r.squaredGLMM(h0_13)

#testOutliers(h0_13)
#testOutliers(sim_res)
#check_outliers(sim_res, type = "default")
#check_outliers(sim_res, type = "binomial")
#check_outliers(sim_res, type = "bootstrap")

#Adding total length as a predictor/covariate based on recommendation by committee
#which makes sense since it is a continuous variable and likely more appropriate as a fixed effect
h0_13_test <- glmmTMB(h0 ~ env_broad_scale + total_length + (1 | collection_date) + 
                        (1 | sample_name), data = meta_all, family = nbinom2())

AIC(h0_13, h0_13_test)
AICc(h0_13)
AICc(h0_13_test)
#both models have similar AIC and AICc values

overdisp_ratio = deviance(h0_13_test)/df.residual(h0_13_test)
overdisp_ratio

summary(h0_13_test)

plot(residuals(h0_13_test))

qqnorm(residuals(h0_13_test))
qqline(residuals(h0_13_test))

sim_res <- simulateResiduals(h0_13_test)
plot(sim_res)

Anova(h0_13_test, type = 2)

group_means <- emmeans(h0_13_test, ~ env_broad_scale)
group_means

emmeans(h0_13_test, pairwise ~ env_broad_scale)

r.squaredGLMM(h0_13_test)


#Original Conclusion: use h0_13 based on AIC and residuals
#CHANGE IN CONCLUSION: Use h0_13_test based on comment/suggestion for including total length as a fixed effect


## Testing additional random effects in the final model - to account for potential batch effects and sample location effects ##

#Original model
h0_13_test <- glmmTMB(h0 ~ env_broad_scale + total_length + (1 | collection_date) + 
                        (1 | sample_name), data = meta_all, family = nbinom2())

#Adding plate ID as a random effect
h0_13_test2 <- glmmTMB(h0 ~ env_broad_scale + total_length + (1 | collection_date) + 
                         (1 | sample_name) + (1 | plate_ID), data = meta_all, family = nbinom2())


#Adding sequencing run number as a column in the metadata 
meta_all <- meta_all %>% 
  mutate(seq_run = case_when(plate_ID == "UHM_Hellbender_P1_CEC" ~ 1, 
                             plate_ID == "UHM_Hellbender_P2_CEC" ~ 1,
                             TRUE ~ 2))
meta_all$seq_run <- as.factor(meta_all$seq_run)


#Adding sequencing run number as a random effect
h0_13_test3 <- glmmTMB(h0 ~ env_broad_scale + total_length + (1 | collection_date) + 
                         (1 | sample_name) + (1 | seq_run), data = meta_all, family = nbinom2())

#Adding sequencing run number as a fixed effect - just for reference
h0_13_test4 <- glmmTMB(h0 ~ env_broad_scale + total_length + seq_run + (1 | collection_date) + 
                         (1 | sample_name), data = meta_all, family = nbinom2())

#Adding broad location/group as a random effect
h0_13_test5 <- glmmTMB(h0 ~ env_broad_scale + total_length + (1 | collection_date) + 
                         (1 | sample_name) + (1 | group_broad), data = meta_all, family = nbinom2())
#Here group_broad = Mid TN Wild, Mid TN Recap, Nash Zoo, Chatt Zoo, and East TN Wild

#Adding broad locale as a column in the metadata 
meta_all <- meta_all %>% 
  mutate(broad_locale = case_when(env_broad_scale == "wild" & site == "Hiwassee River Picnic Area" ~ "East TN",
                                  env_broad_scale == "wild" & ecoregion_III == "Interior Plateau" ~ "Middle TN",
                                  site == "Chattanooga Zoo" ~ "Chattanooga Zoo",
                                  site == "Nashville Zoo" ~ "Nashville Zoo", 
                                  TRUE ~ "Middle TN"))
meta_all$broad_locale <- as.factor(meta_all$broad_locale)

#Adding broad locale as a random effect
h0_13_test6 <- glmmTMB(h0 ~ env_broad_scale + total_length + (1 | collection_date) + 
                         (1 | sample_name) + (1 | broad_locale), data = meta_all, family = nbinom2())
#Here broad_locale = East TN, Chatt Zoo, Nash Zoo, and Mid TN
#Similar to broad_group but Mid TN wild and recaps are lumped together

AIC(h0_13_test, h0_13_test2, h0_13_test3, h0_13_test4, h0_13_test5, h0_13_test6)
AICc(h0_13_test)
AICc(h0_13_test2)
AICc(h0_13_test3)
AICc(h0_13_test4)
AICc(h0_13_test5)
AICc(h0_13_test6)


summary(h0_13_test)
summary(h0_13_test2)
summary(h0_13_test3)
summary(h0_13_test4)
summary(h0_13_test5)
summary(h0_13_test6)



#2. Evaluating whether Bray-Curtis beta dispersion is affected by captive/wild setting (without Nashville Zoo subsampled data)----

#Extracting the Bray-Curtis distances generated by betadisper
meta_bray_all <- cbind(meta_all, bray_dis = bray_beta_cw$distances)

#Testing for normality of the Bray-Curtis distances
shapiro.test(meta_bray_all$bray_dis) #not normally distributed, p = 4.026e-06

#Comparing an ordered beta distribution and a gamma distribution
bc_1 <- glmmTMB(bray_dis ~ env_broad_scale + (1 | sample_name), data = meta_bray_all, family = Gamma())

bc_2 <- glmmTMB(bray_dis ~ env_broad_scale + (1 | sample_name), data = meta_bray_all, family = ordbeta())
#Choosing to test an ordered beta distribution since its values are bound between 0 and 1 (like proportions or probabilities)
#like the distances produced by betadisper whereas gamma can model continuous positive values with larger numbers

AIC(bc_1, bc_2)
AICc(bc_1)
AICc(bc_2)
#AIC and AICc for ordered beta is lower 

#Determining whether to include richness as a fixed effect
bc_3 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name), data = meta_bray_all, family = ordbeta())

AIC(bc_3, bc_2)
AICc(bc_3)
AICc(bc_2)
#Including richness as a fixed effect improved model fit

#Testing the inclusion of additional random effects
bc_4 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name) + (1 | dna_experiment_date), 
                data = meta_bray_all, family = ordbeta())

bc_5 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name) + (1 | dna_experiment_date) +
                  (1 | plate_ID), data = meta_bray_all, family = ordbeta())
#bc_5 produced multiple warnings - be cautious of it

AIC(bc_4, bc_5, bc_3)
AICc(bc_4)
AICc(bc_5)
AICc(bc_3)
#Inclusion of dna exp. date with bc_4 improved fit based on AIC and AICc

#Adding more random effects
bc_6 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name) + (1 | dna_experiment_date) + 
                  (1 | collection_date), data = meta_bray_all, family = ordbeta())

#bc_7 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name) + (1 | dna_experiment_date) + 
#                  (1 | collection_date) + (1 | total_length), data = meta_bray_all, family = ordbeta())

AIC(bc_4, bc_6)
AICc(bc_4)
AICc(bc_6)

summary(bc_6)

group_means <- emmeans(bc_6, ~ env_broad_scale + h0)
group_means

emmeans(bc_6, pairwise ~ env_broad_scale + h0)

plot(residuals(bc_6))

qqnorm(residuals(bc_6))
qqline(residuals(bc_6))
#residuals look fairly decent

sim_res <- simulateResiduals(bc_6)
plot(sim_res)
#significant KS test but non-significant outlier and dispersion tests

Anova(bc_6, type = 2)

r.squaredGLMM(bc_6)

#testOutliers(bc_6)
#testOutliers(sim_res)
#check_outliers(sim_res, type = "default")
#check_outliers(sim_res, type = "binomial")
#check_outliers(sim_res, type = "bootstrap")



#3. Evaluating whether Raup-Crick beta dispersion is affected by captive/wild setting (without Nashville Zoo subsampled data)----

#Using output from vegdist function
raup_beta_cw <- betadisper(raup_dist_sub, meta_sub$env_broad_scale)

#Extracting the Raup-Crick distances generated by betadisper
meta_raup_sub <- cbind(meta_sub, raup_dis = raup_beta_cw$distances)

#Testing for normality of the Raup-Crick distances
shapiro.test(meta_raup_sub$raup_dis) #not normally distributed, p = 0.0001371

#Comparing an ordered beta distribution and a gamma distribution
#rc_1 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name), data = meta_raup_sub, family = Gamma())

rc_2 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name), data = meta_raup_sub, family = ordbeta())
#Choosing to test an ordered beta distribution since its values are bound between 0 and 1 (like proportions or probabilities)
#like the distances produced by betadisper whereas gamma can model continuous positive values with larger numbers

#AIC(rc_1, rc_2)
#Gamma distribution does not work since Raup-Crick matrix can produce 0 values as distances
AIC(rc_2)
AICc(rc_2)

#Determining whether to include richness as a fixed effect
rc_3 <- glmmTMB(raup_dis ~ h0 + env_broad_scale + (1 | sample_name), data = meta_raup_sub, family = ordbeta())

AIC(rc_3, rc_2)

summary(rc_3)
#AIC is lower when richness is included but seems to cause the model to fail somehow 
#because no values are produced in the summary output, so proceed with rc_2

#Testing the inclusion of additional random effects
rc_4 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date), 
                data = meta_raup_sub, family = ordbeta())

rc_5 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date) +
                  (1 | plate_ID), data = meta_raup_sub, family = ordbeta())
#rc_5 produced a model convergence error - do not use

AIC(rc_4, rc_5, rc_2)
AICc(rc_4)
AICc(rc_2)
#Inclusion of random effects did not significantly improve fit based on AICc


#Adding more random effects
rc_6 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name) + (1 | collection_date), 
                data = meta_raup_sub, family = ordbeta())

#rc_7 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name) + (1 | total_length), 
#                data = meta_raup_sub, family = ordbeta())

AIC(rc_2, rc_6)
AICc(rc_2)
AICc(rc_6)
#Inclusion of random effects did not significantly improve fit based on AICc
#However, since the values are somewhat similar for rc_2, rc_4, rc_6, perhaps model averaging would be appropriate here

summary(rc_2)

group_means <- emmeans(rc_2, ~ env_broad_scale)
group_means

emmeans(rc_2, pairwise ~ env_broad_scale)

plot(residuals(rc_2))

qqnorm(residuals(rc_2))
qqline(residuals(rc_2))
#residuals look a little wonky

sim_res <- simulateResiduals(rc_2)
plot(sim_res)
#but simulated residuals are fine

Anova(rc_2, type = 2)

r.squaredGLMM(rc_2)
#Got a warning message here


#testOutliers(rc_2)
#testOutliers(sim_res)
#check_outliers(sim_res, type = "default")
#check_outliers(sim_res, type = "binomial")
#check_outliers(sim_res, type = "bootstrap")



##### Double checking to make sure that glycerol stored swabs (n = 6) from Nashville Zoo do not differ from dry swabs #####

#Creating a list of individual animals with both dry swabs and glycerol swabs
gly_swabs <- c("UHM783", "UHM789", "UHM790", "UHM791", "UHM806", "UHM807")

#Filtering the metadata to contain these glycerol and dry swabs from the same set of individuals
meta_gly <- filter(meta_all, sample_name %in% gly_swabs)
#n = 6 glycerol swabs, n = 11 dry stored swabs

#Conducting a t-test for richness/alpha diversity with all of the dry and glycerol swabs
t.test(h0 ~ storage_medium, data = meta_gly)
#no significant difference in richness, p = 0.6426

#Setting up two data frames to create an equal number of swabs for a paired t-test
meta_gly2 <- filter(meta_all, sample_name %in% gly_swabs, storage_medium == "glycerol")

meta_dry2 <- filter(meta_all, sample_name %in% gly_swabs, storage_medium == "dry")

set.seed(050521) #set seed for reproducibility

#Randomly selecting 6 dry stored swabs
meta_sub_dry <- meta_dry2 %>% 
  group_by(storage_medium) %>% 
  slice_sample(n = 6, replace = F) %>%
  ungroup()

#Conducting a paired t-test
t.test(meta_gly2$h0, meta_sub_dry$h0, paired = T)
#no significant difference in richness, p = 0.6996

#Randomly selecting 1 dry stored swab for each individual animal 
#(so that every individual is represented and has a pair of dry and glycerol swabs)
set.seed(050521) #set seed for reproducibility

meta_sub_dry_samp <- meta_dry2 %>% 
  group_by(sample_name) %>% 
  slice_sample(n = 1, replace = F) %>%
  ungroup()

#conducting paired t-test
t.test(meta_gly2$h0, meta_sub_dry_samp$h0, paired = T)
#no significant difference in richness, p = 0.7383


#PERMANOVA for Bray and Raup for beta diversity

  ##Bray-Curtis##

#Making sure sample count is the same after filtering out samples
meta_otu_inter <- intersect(meta_gly$index, otu_count$index)
otu_table_2_gly <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
otu_table_2_gly <- otu_table_2_gly[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
otu_table_2_gly <- otu_table_2_gly %>% 
  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_2_gly))
tax_table_2_gly <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_2_gly <- otu_table_2_gly %>% select(all_of(intersected))

#Calculating Bray-Curtis distances
set.seed(050521) #set seed for reproducibility
bray_dist_gly <- vegdist(otu_table_2_gly, method = "bray", na.rm = T)

#Checking for homogeneity of dispersion between captive and wild-caught samples
bray_beta_gly <- betadisper(bray_dist_gly, meta_gly$storage_medium)
permutest(bray_beta_gly)
TukeyHSD(bray_beta_gly)
#significant difference in dispersion between glycerol and dry, p = 0.013

boxplot(bray_beta_gly)
#dry swabs have higher dispersion and higher n = so PERMANOVA can proceed as is

#Performing the PERMANOVA
bray_adonis_gly <- adonis2(bray_dist_gly ~ storage_medium, data = meta_gly,
                           permutations = 999, by = "terms")
bray_adonis_gly
#no significant difference between dry and glycerol, p = 0.064


  ## Raup-Crick ##

#Calculating Raup-Crick distances
set.seed(050521) #set seed for reproducibility
raup_dist_gly <- vegdist(otu_table_2_gly, method = "raup", na.rm = T)

#Checking for homogeneity of dispersion between captive and wild-caught samples
raup_beta_gly <- betadisper(raup_dist_gly, meta_gly$storage_medium)
permutest(raup_beta_gly)
TukeyHSD(raup_beta_gly)
#no significant difference in dispersion between glycerol and dry, p = 0.272

#Performing the PERMANOVA
raup_adonis_gly <- adonis2(raup_dist_gly ~ storage_medium, data = meta_gly,
                           permutations = 999, by = "terms")
raup_adonis_gly
#no significant difference between dry and glycerol, p = 0.106

#Conclusion: there seem to be no significant differences in alpha or beta diversity for 
#dry versus glycerol stored cloacal swabs



##### Subsampling Nashville Zoo data to account for potential pseudoreplication #####

#Filtering for Chatt Zoo and wild-caught animals (East TN, Middle TN, and recaps) to create a new data frame for NMDS plotting later
meta_Chatt_wild <- filter(meta_all, site != "Nashville Zoo")

#Filtering for just Nashville Zoo individuals to add back to the data frame later
meta_Nash <- filter(meta_all, site == "Nashville Zoo")

#Creating a random subsample of Nashville Zoo animals where every individual is 
#sampled only once at across random time points (to avoid pseudoreplication)
set.seed(050521) #set seed for reproducibility
sampled_data_Nash <- meta_Nash %>% 
  group_by(sample_name) %>% 
  slice_sample(n = 1, replace = F)

#Joining the randomly subsampled Nashville Zoo data frame with the Chatt Zoo and wild-caught df for NMDS plotting
meta_sub_NMDS <- rbind(sampled_data_Nash, meta_Chatt_wild)

#Making sure sample count is the same after filtering out samples
meta_otu_inter <- intersect(meta_sub_NMDS$index, otu_count$index)
otu_table_2_sub_NMDS <- otu_count %>% filter(index %in% meta_otu_inter)

#Filtering Escherichia Shigella from the OTU table
otu_table_2_sub_NMDS <- otu_table_2_sub_NMDS[,-2] #filtering out the OTU00001 column for E. Shigella

#Format the OTU table for the better plotting, joining, etc.  
otu_table_2_sub_NMDS <- otu_table_2_sub_NMDS %>% 
  column_to_rownames("index")

#Making sure that the OTUs in the OTU table and taxonomy table are the same after removing E. Shigella
intersected <- intersect(tax_table_2$otu, colnames(otu_table_2_sub_NMDS))
tax_table_2_sub_NMDS <- tax_table_2 %>% filter(otu %in% intersected)
otu_table_2_sub_NMDS <- otu_table_2_sub_NMDS %>% select(all_of(intersected))

#Coding group_broad as a factor for the pairwise PERMANOVA
meta_sub_NMDS$group_broad <- as.factor(meta_sub_NMDS$group_broad)
meta_sub_NMDS$site <- as.factor(meta_sub_NMDS$site)



##### Pairwise PERMANOVA function code #####

#Creating the function for the pairwise PERMANOVA
#based on code from https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/complex-models/
#renamed to pair.adonis2 to avoid conflicting with pairwise.adonis2 in pairwiseAdonis package
pair.adonis2 <- function(resp, fact, p.method = "none", nperm = 999) {
  require(vegan)
  resp <- as.matrix(resp)
  fact <- factor(fact)
  fun.p <- function(i, j) {
    fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, j)])
    index <- which(fact %in% levels(fact2))
    resp2 <- as.dist(resp[index, index])
    result <- adonis2(resp2 ~ fact2, permutations = nperm)
    result$`Pr(>F)`[1]
  }
  multcomp <- pairwise.table(fun.p, levels(fact), p.adjust.method = p.method)
  return(list(fact = levels(fact), p.value = multcomp, p.adjust.method = p.method))
}



##### Bray-Curtis NMDS - subsampled Nashville Zoo data #####

#Calculating Bray-Curtis distances
set.seed(050521) #set seed for reproducibility
bray_dist_sub_NMDS <- vegdist(otu_table_2_sub_NMDS, method = "bray", na.rm = T)

#Checking for homogeneity of dispersion between groups
bray_beta_group <- betadisper(bray_dist_sub_NMDS, meta_sub_NMDS$group_broad)
permutest(bray_beta_group)
TukeyHSD(bray_beta_group)
#significant difference in dispersion based on broad group/site, p = 0.001
#significant difference between East TN-Chatt Zoo and Nash Zoo-Chatt Zoo

boxplot(bray_beta_group)
#both Nash Zoo and East TN have a larger n = and higher dispersion than Chatt Zoo
#so PERMANOVA should be fine to proceed as is

#Checking for homogeneity of dispersion between captive and wild-caught samples
bray_beta_cw <- betadisper(bray_dist_sub_NMDS, meta_sub_NMDS$env_broad_scale)
permutest(bray_beta_cw)
TukeyHSD(bray_beta_cw)
#no significant difference in dispersion between captive vs. wild, p = 0.34

#Performing the PERMANOVA WITHOUT strata
bray_adonis <- adonis2(bray_dist_sub_NMDS ~ h0 * group_broad, data = meta_sub_NMDS,
                       permutations = 999, by = "terms")
bray_adonis
#significant effect of group_broad and h0:group_broad

#Performing a pairwise PERMANOVA WITHOUT strata
#bray_adonis_pair <- pairwise.adonis2(bray_dist_sub_NMDS ~ h0 * group_broad, data = meta_sub_NMDS,
#                       permutations = 999, by = "terms")
#bray_adonis_pair

#Stratifying according to site to account for spatial autocorrelation
bray_adonis_strata <- adonis2(bray_dist_sub_NMDS ~ h0 * group_broad, data = meta_sub_NMDS,
                              permutations = 999, by = "terms", strata = meta_sub_NMDS$site)
bray_adonis_strata
#significant effect of h0:group_broad with ecoregion_III strata and with site strata

#Testing PERMANOVA with cap/wild variable included for comparison - since that is ultimately the question being asked
bray_adonis_strata2 <- adonis2(bray_dist_sub_NMDS ~  env_broad_scale + h0 * group_broad, data = meta_sub_NMDS,
                               permutations = 999, by = "terms", strata = meta_sub_NMDS$site)
bray_adonis_strata2
#significant effect of env_broad_scale and h0:group_broad with site strata

#Performing a pairwise PERMANOVA with site strata
#bray_adonis_pair_strata <- pairwise.adonis2(bray_dist_sub_NMDS ~ h0 * group_broad, data = meta_sub_NMDS,
#                                  permutations = 999, strata = 'site')
#bray_adonis_pair_strata


#Pairwise PERMANOVA with hard coded function without site strata
#bray_pair <- pair.adonis2(resp = bray_dist_sub_NMDS, fact = meta_sub_NMDS$group_broad )
#bray_pair

bray_nmds <- metaMDS(bray_dist_sub_NMDS, k = 3) #conduct non-metric multidimensional scaling

bray_nmds$stress
#stress value is 0.1165133, which is decent


#Creating a data frame with the Bray-Curtis distances and metadata
bray_nmds_df <- as.data.frame(bray_nmds[["points"]]) %>% #generate dataframe with NMDS points 
  rownames_to_column("index") %>% 
  left_join(meta_all, by = join_by("index"))

#Plotting the Bray-Curtis NMDS - WITH Mid TN Recap ellipse
bray_nmds_plot <- ggplot(bray_nmds_df, aes(x = MDS1, y = MDS2, color = group_broad, shape = env_broad_scale, fill = group_broad))+
  stat_ellipse(aes(group = group_broad, fill = group_broad), type = "t", geom = "polygon", 
               alpha = 0.25, linetype = "twodash") +
  scale_shape_manual(values = c(22, 25), labels = c("Wild–Caught", "Zoo")) +
  geom_point(size = 2.5, alpha = 0.75, color = "black") +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") + 
  labs(color = "Site", shape = "Environmental Setting", fill = "Site") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate("text", x = -Inf, y = -Inf, label = "Stress = 0.12", hjust = -6.32, vjust = -1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
print(bray_nmds_plot)

#Adding centroid positions to the NMDS
bray_nmds_plot_cent <- bray_nmds_plot +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = group_broad), size = 5, shape = 21)
print(bray_nmds_plot_cent)

bray_nmds_plot_cent + guides(shape = guide_legend(order = 1),
                         size = guide_legend(order = 1),
                         color = guide_legend(order = 2),
                         fill = guide_legend(order = 2))


#Getting Brewer palette has codes for plotting
brewer.pal(n = 8, "Accent")

#Plotting the Bray-Curtis NMDS - WITHOUT Mid TN Recap ellipse 
bray_nmds_plot_MTR_no_ellipse <- ggplot(bray_nmds_df, aes(x = MDS1, y = MDS2, color = group_broad, shape = env_broad_scale, fill = group_broad))+
  stat_ellipse(data = subset(bray_nmds_df, group_broad != "Middle TN Recapture"), aes(group = group_broad, fill = group_broad), type = "t", geom = "polygon", 
               alpha = 0.25, linetype = "twodash") +
  scale_shape_manual(values = c(22, 25), labels = c("Wild–Caught", "Zoo")) +
  geom_point(size = 2.5, alpha = 0.75, color = "black") +
  scale_fill_brewer(palette = "Accent") +
  scale_color_manual(values = c("Chattanooga Zoo" = "#7FC97F", "East TN Wild" = "#BEAED4",
                                "Middle TN Recapture" = NA, "Middle TN Wild" = "#FFFF99","Nashville Zoo" = "#386CB0")) + 
  labs(color = "Site", shape = "Environmental Setting", fill = "Site") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate("text", x = -Inf, y = -Inf, label = "Stress = 0.12", hjust = -0.15, vjust = -1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
print(bray_nmds_plot_MTR_no_ellipse)

#Adding centroid positions to the NMDS
bray_nmds_plot_cent_MTR <- bray_nmds_plot_MTR_no_ellipse +
  stat_ellipse(level = 1e-10, geom = "point", data = subset(bray_nmds_df, group_broad != "Middle TN Recapture"), aes(fill = group_broad), size = 5, shape = 21)
print(bray_nmds_plot_cent_MTR)



##### Raup-Crick NMDS - subsampled Nashville Zoo data #####

#Calculating Raup-Crick distances
set.seed(050521) #set seed for reproducibility

###With all of the data###
#raup_dist_all <- vegdist(otu_table_2, method = "raup", na.rm = T)

#raup_dist_all2 <- raupcrick(otu_table_2, null = "r1", nsim = 999)
#saveRDS(raup_dist_all2, file = "241211_raupcrick_dist.Rds")
#raup_dist_all2 <- readRDS("241211_raupcrick_dist.Rds")

#With a subsample of the data
raup_dist_sub_NMDS <- vegdist(otu_table_2_sub_NMDS, method = "raup", na.rm = T)

#raup_dist_sub2 <- raupcrick(otu_table_2_sub, null = "r1", nsim = 999)
#saveRDS(raup_dist_sub2, file = "241211_raupcrick_sub_dist.Rds")
#raup_dist_sub2 <- readRDS("241211_raupcrick_sub_dist.Rds")

#Checking for homogeneity of dispersion
raup_beta_group <- betadisper(raup_dist_sub_NMDS, meta_sub_NMDS$group_broad)
permutest(raup_beta_group)
TukeyHSD(raup_beta_group)
#no significant difference in dispersion between broad group/site, p = 0.101

#Checking for homogeneity of dispersion
raup_beta_cw <- betadisper(raup_dist_sub_NMDS, meta_sub_NMDS$env_broad_scale)
permutest(raup_beta_cw)
TukeyHSD(raup_beta_cw)
#no sig. difference between captive vs. wild, p = 0.982 

#Performing the PERMANOVA
raup_adonis <- adonis2(raup_dist_sub_NMDS ~ h0 * group_broad, data = meta_sub_NMDS,
                       permutations = 999, by = "terms")
raup_adonis
#sig. effect of group_broad and h0:group_broad, while main effect of h0 is non-sig and has neg. R2

#Performing PERMANOVA with site strata
raup_adonis_strata <- adonis2(raup_dist_sub_NMDS ~ h0 * group_broad, data = meta_sub_NMDS,
                              permutations = 999, by = "terms", strata = meta_sub_NMDS$site)
raup_adonis_strata
#sig. effect of group_broad and h0:group_broad with ecoregion III strata and with site strata

#Testing PERMANOVA with cap/wild variable included for comparison - since that is ultimately the question being asked
raup_adonis_strata2 <- adonis2(raup_dist_sub_NMDS ~  env_broad_scale + h0 * group_broad, data = meta_sub_NMDS,
                               permutations = 999, by = "terms", strata = meta_sub_NMDS$site)
raup_adonis_strata2
#significant effect of env_broad_scale, group_broad, and h0:group_broad with site strata

#Pairwise PERMANOVA with hard coded function without site strata
#raup_pair <- pair.adonis2(resp = raup_dist_sub_NMDS, fact = meta_sub_NMDS$group_broad )
#raup_pair

raup_nmds <- metaMDS(raup_dist_sub_NMDS, k = 3) #conduct non-metric multidimensional scaling

raup_nmds$stress
#stress value is 0.126917

raup_nmds_df <- as.data.frame(raup_nmds[["points"]]) %>% #generate dataframe with NMDS points 
  rownames_to_column("index") %>% 
  left_join(meta_sub_NMDS, by = join_by("index"))

#Plotting the Raup-Crick NMDS - WITH Mid TN Recap ellipse
raup_nmds_plot <- ggplot(raup_nmds_df, aes(x = MDS1, y = MDS2, color = group_broad, shape = env_broad_scale, fill = group_broad))+
  stat_ellipse(aes(group = group_broad, fill = group_broad), type = "t", geom = "polygon", 
               alpha = 0.25, linetype = "twodash") +
  scale_shape_manual(values = c(22, 25), labels = c("Wild-Caught", "Zoo")) +
  geom_point(size = 2.5, alpha = 0.75, color = "black") +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") + 
  labs(color = "Site", shape = "Environmental Setting", fill = "Site") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate("text", x = -Inf, y = -Inf, label = "Stress = 0.13", hjust = -6.32, vjust = -1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
print(raup_nmds_plot)

#Adding centroid positions to the NMDS
raup_nmds_plot_cent <- raup_nmds_plot +
  stat_ellipse(level = 1e-10, geom = "point", aes(fill = group_broad), size = 5, shape = 21)
print(raup_nmds_plot_cent)

raup_nmds_plot_cent + guides(shape = guide_legend(order = 1),
                             size = guide_legend(order = 1),
                             color = guide_legend(order = 2),
                             fill = guide_legend(order = 2))

#Creating a multipanel figure of Bray and Raup with the centroids - WITH Mid TN Recap ellipse
bray_raup_multi <- ggarrange(bray_nmds_plot_cent, raup_nmds_plot_cent, ncol = 1, 
                             nrow = 2, common.legend = TRUE, legend = "right", labels = c("a)", "b)"))
print(bray_raup_multi)

#moving the shape portion of the legend prior to the color portion
bray_raup_multi + guides(shape = guide_legend(order = 1),
                         size = guide_legend(order = 1),
                         color = guide_legend(order = 2),
                         fill = guide_legend(order = 2))
print(bray_raup_multi)


#Plotting the Raup-Crick NMDS - WITHOUT Mid TN Recap ellipse
raup_nmds_plot_MTR_no_ellipse <- ggplot(raup_nmds_df, aes(x = MDS1, y = MDS2, color = group_broad, shape = env_broad_scale, fill = group_broad))+
  stat_ellipse(data = subset(raup_nmds_df, group_broad != "Middle TN Recapture"), aes(group = group_broad, fill = group_broad), type = "t", geom = "polygon", 
               alpha = 0.25, linetype = "twodash") +
  scale_shape_manual(values = c(22, 25), labels = c("Wild-Caught", "Zoo")) +
  geom_point(size = 2.5, alpha = 0.75, color = "black") +
  #scale_fill_brewer(palette = "Accent") +
  scale_fill_manual(values = c("Chattanooga Zoo" = "#7FC97F", "East TN Wild" = "#BEAED4","Middle TN Recapture" = "#FDC086", "Middle TN Wild" = "#FFFF99","Nashville Zoo" = "#386CB0")) + 
  scale_color_manual(values = c("Chattanooga Zoo" = "#7FC97F", "East TN Wild" = "#BEAED4","Middle TN Recapture" = NA, "Middle TN Wild" = "#FFFF99","Nashville Zoo" = "#386CB0")) + 
  labs(color = "Site", shape = "Environmental Setting", fill = "Site") +
  xlab("NMDS1") +
  ylab("NMDS2") +
  annotate("text", x = -Inf, y = -Inf, label = "Stress = 0.13", hjust = -0.15, vjust = -1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")
print(raup_nmds_plot_MTR_no_ellipse)

#Adding centroid positions to the NMDS
raup_nmds_plot_cent_MTR <- raup_nmds_plot_MTR_no_ellipse +
  stat_ellipse(level = 1e-10, geom = "point", data = subset(raup_nmds_df, group_broad != "Middle TN Recapture"), aes(fill = group_broad), size = 5, shape = 21)
print(raup_nmds_plot_cent_MTR)

#Creating a multipanel figure of Bray and Raup with the centroids - WITHOUT Mid TN Recap ellipse
bray_raup_multi_MTR <- ggarrange(bray_nmds_plot_cent_MTR, raup_nmds_plot_cent_MTR, ncol = 1, 
                             nrow = 2, common.legend = TRUE, legend = "right", labels = c("a)", "b)"))
print(bray_raup_multi_MTR)



##### Redoing statistical modeling with subsampled Nashville Zoo data #####

#4. Evaluating whether Bray-Curtis beta dispersion is affected by captive/wild setting (with subsampled Nashville Zoo data)----

#Extracting the Bray-Curtis distances generated by betadisper
#meta_bray_sub_nmds <- cbind(meta_sub_NMDS, bray_dis = bray_beta_cw$distances)

#Testing for normality of the Bray-Curtis distances
#shapiro.test(meta_bray_sub_nmds$bray_dis) #not normally distributed, p = 0.009266

#Comparing an ordered beta distribution and a gamma distribution
#bc_1 <- glmmTMB(bray_dis ~ env_broad_scale + (1 | sample_name), data = meta_bray_sub_nmds, family = Gamma())

#bc_2 <- glmmTMB(bray_dis ~ env_broad_scale + (1 | sample_name), data = meta_bray_sub_nmds, family = ordbeta())
#Choosing to test an ordered beta distribution since its values are bound between 0 and 1 (like proportions or probabilities)
#like the distances produced by betadisper whereas gamma can model continuous positive values with larger numbers

#AIC(bc_1, bc_2)
#AICc(bc_1)
#AICc(bc_2)
#AIC and AICc for ordered beta is slightly lower 

#Determining whether to include richness as a fixed effect
#bc_3 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name), data = meta_bray_sub_nmds, family = ordbeta())

#AIC(bc_3, bc_2)
#AICc(bc_3)
#AICc(bc_2)
#Including richness as a fixed effect improved model fit

#Testing the inclusion of additional random effects
#bc_4 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name) + (1 | dna_experiment_date), 
#                data = meta_bray_sub_nmds, family = ordbeta())

#bc_5 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_bray_sub_nmds, family = ordbeta())
#bc_4 and bc_5 produced multiple warnings - be cautious of it

#AIC(bc_4, bc_5, bc_3)
#AICc(bc_4)
#AICc(bc_5)
#AICc(bc_3)
#bc_4 and bc_5 did not improve model fit - proceed with bc_3

#Adding more random effects
#bc_6 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name) + (1 | dna_experiment_date) + 
#                  (1 | collection_date), data = meta_bray_sub_nmds, family = ordbeta())

#bc_7 <- glmmTMB(bray_dis ~ env_broad_scale + h0 + (1 | sample_name) + (1 | dna_experiment_date) + 
#                  (1 | collection_date) + (1 | total_length), data = meta_bray_sub_nmds, family = ordbeta())

#bc_6 also produced multiple warnings

#AIC(bc_3, bc_6)
#AICc(bc_3)
#AICc(bc_6)

#summary(bc_3)
#the model keeps producing NaN values for everything but the estimates
#the summary for bc_2 produces NaN values as well even though it only has one fixed and one random effect

#group_means <- emmeans(bc_3, ~ env_broad_scale + h0)
#group_means

#emmeans(bc_3, pairwise ~ env_broad_scale + h0)

#plot(residuals(bc_3))

#qqnorm(residuals(bc_3))
#qqline(residuals(bc_3))
#residuals look ok

#sim_res <- simulateResiduals(bc_3)
#plot(sim_res)
#non-significant KS test, outlier test, and dispersion tests

#Anova(bc_3, type = 2)

#r.squaredGLMM(bc_3)

#testOutliers(bc_6)
#testOutliers(sim_res)
#check_outliers(sim_res, type = "default")
#check_outliers(sim_res, type = "binomial")
#check_outliers(sim_res, type = "bootstrap")



#5. Evaluating whether Raup-Crick beta dispersion is affected by captive/wild setting (with subsampled Nashville Zoo data)----

#Extracting the Raup-Crick distances generated by betadisper
#meta_raup_sub_nmds <- cbind(meta_sub_NMDS, raup_dis = raup_beta_cw$distances)

#Testing for normality of the Raup-Crick distances
#shapiro.test(meta_raup_sub_nmds$raup_dis) #not normally distributed, p = 4.502e-10

#Comparing an ordered beta distribution and a gamma distribution
#rc_1 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name), data = meta_raup_sub_nmds, family = Gamma())

#rc_2 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name), data = meta_raup_sub_nmds, family = ordbeta())
#Choosing to test an ordered beta distribution since its values are bound between 0 and 1 (like proportions or probabilities)
#like the distances produced by betadisper whereas gamma can model continuous positive values with larger numbers

#AIC(rc_1, rc_2)
#Gamma distribution does not work since Raup-Crick matrix can produce 0 values as distances
#AIC(rc_2)
#AICc(rc_2)

#Determining whether to include richness as a fixed effect
#rc_3 <- glmmTMB(raup_dis ~ h0 + env_broad_scale + (1 | sample_name), data = meta_raup_sub_nmds, family = ordbeta())

#AIC(rc_3, rc_2)
#AICc(rc_3)
#AICc(rc_2)

#summary(rc_3)


#Testing the inclusion of additional random effects
#rc_4 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date), 
#                data = meta_raup_sub_nmds, family = ordbeta())

#rc_5 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name) + (1 | dna_experiment_date) +
#                  (1 | plate_ID), data = meta_raup_sub_nmds, family = ordbeta())

#AIC(rc_4, rc_5, rc_2)
#AICc(rc_4)
#AICc(rc_5)
#AICc(rc_2)
#Inclusion of random effects did not significantly improve fit based on AICc


#Adding more random effects
#rc_6 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name) + (1 | collection_date), 
#                data = meta_raup_sub_nmds, family = ordbeta())

#rc_7 <- glmmTMB(raup_dis ~ env_broad_scale + (1 | sample_name) + (1 | total_length), 
#                data = meta_raup_sub_nmds, family = ordbeta())

#AIC(rc_2, rc_6)
#AICc(rc_2)
#AICc(rc_6)
#Inclusion of random effects did not significantly improve fit based on AICc

#summary(rc_2)
#env_broad_scale not significant

#group_means <- emmeans(rc_2, ~ env_broad_scale)
#group_means

#emmeans(rc_2, pairwise ~ env_broad_scale)

#plot(residuals(rc_2))

#qqnorm(residuals(rc_2))
#qqline(residuals(rc_2))
#residuals look a little wonky

#sim_res <- simulateResiduals(rc_2)
#plot(sim_res)
#but simulated residuals are fine

#Anova(rc_2, type = 2)

#r.squaredGLMM(rc_2)
#Got a warning message here and the R2 value is pretty low overall

#testOutliers(rc_2)
#testOutliers(sim_res)
#check_outliers(sim_res, type = "default")
#check_outliers(sim_res, type = "binomial")
#check_outliers(sim_res, type = "bootstrap")
