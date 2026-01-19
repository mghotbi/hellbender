#Researcher: Chloe Cummins
#Start Date: January 31, 2025
#Analyzing the fecal mycobiome of hellbenders by creating a rel. abund. bar chart



#Set working directory to manuscript code check folder
setwd("G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Hellbender Gut Microbiome Manuscript Stuff/Manuscript_Draft_Code_and_Figs")

#Code originally stored in "G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Hellbender Gut Microbiome Manuscript Stuff/Manuscript_Code_Check"
#and in "G:/Shared drives/Microbiome Ecology Lab Shared Drive/Chloe Cummins PC backup/Database/Data/Thesis_ITS_Analyses"

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
library(MiscMetabar)

#Importing the ITS data with no contamination and spiked species removed
physeq <- readRDS("physeq_ITSOTUNoDek.rds")

#Breaking the phyloseq into its separate components because I need to add a column to the metadata for plotting
#Once I add that column to the metadata, I will join everything back together again
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

#Joining everything back together to create the updated phyloseq object
ASV = otu_table(otu_mat, taxa_are_rows = T)
TAX = tax_table(tax_mat)
metadata.phyloseq = sample_data(meta_df)

physeq_updated = phyloseq(ASV, TAX, metadata.phyloseq)
physeq_updated


#Transform into relative abundance of the reads
physeq_fungi <- transform_sample_counts(physeq_updated, function(x) x/sum(x))


#Filtering the metadata in the phyloseq to make sure negative controls are removed
physeq_fungi <- subset_samples(physeq_updated, sample_name != "NTC_swab" & 
                                 sample_name != "Pool_neg_ctrl")
physeq_fungi


#Converting the phyloseq object into a data frame based on genus rank
gen_abund_fungi <- physeq_fungi %>% tax_glom(taxrank = "Genus") %>% ##setting to the genus level
  psmelt()  #melting data to long format for plotting
head(gen_abund_fungi)
#As a note, the unassigned OTUs did not have to be removed prior to this because they are listed as blank for genus
#and automatically get dropped during the tax_glom step

#Filtering and modifying the data for plotting
all_data_fungi <- gen_abund_fungi %>% 
  select(Genus, Abundance, Sample, group_broad) %>% #selecting variables of interest
  filter(Abundance != 0) %>% #filtering so that abundance of taxa is at least greater than 0
  mutate(Genus = as.character(Genus))
head(all_data_fungi)

#Preparing to create relative abundance plot based on genus
gen_plot_fungi <- all_data_fungi %>% 
  select(Sample, group_broad, Genus, Abundance) %>% #selecting all variables to be used
  group_by(Sample, group_broad) %>% 
  mutate(totalSum = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Sample, group_broad, Genus) %>% 
  summarize(Abundance = sum(Abundance), 
            totalSum, RelAb = Abundance/totalSum) %>% #calculating relative abundance based on genus abundance per host
  unique() #leaving only unique observations
head(gen_plot_fungi)

#Grouping low abundance/rare taxa into Other category
gen_plot_fungi <- gen_plot_fungi %>% group_by(Sample, group_broad, Genus, totalSum) %>%
  summarise(
    Abundance = sum(Abundance),
    Genus = ifelse(RelAb < 0.03, "Other (< 3%)", Genus)) %>%
  group_by(Sample, group_broad, Genus, totalSum) %>% 
  summarize(Abundance = sum(Abundance), 
            RelAb = Abundance/totalSum) %>% #calculating relative abundance based on genus abundance per host
  unique()
head(gen_plot_fungi)
#specifying low abundance/rare taxa (less than 3%) as Other

#Checking to make sure that the relative abundances are as we would expect
max(gen_plot_fungi$RelAb)
mean(gen_plot_fungi$RelAb)
min(gen_plot_fungi$RelAb)

length(unique(gen_plot_fungi$Genus))

#Setting up the interpolated palette
palette_length <- length(unique(gen_plot_fungi$Genus))
extend_pal <- colorRampPalette(brewer.pal(29, "Set1"))


#Plotting the relative abundance bar chart based on genus
gen_rel_plot_fungi <- ggplot(gen_plot_fungi) +
  geom_col(mapping = aes(x = Sample, y = RelAb, fill = fct_relevel(Genus, c("Other (< 3%)"), after = 29)), 
           color = "black", position = "stack", show.legend = TRUE, alpha = 0.6) +
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
  facet_grid(cols = vars(fct_relevel(group_broad, c("Nashville Zoo", "Middle TN Recapture", "Middle TN Wild"))), 
             scales = "free_x", space = "free_x") +
  #  force_panelsizes(cols = c(1.32, 0.27, 0.2, 0.33, 0.25)) +
  guides(fill = guide_legend(ncol = 2))
print(gen_rel_plot_fungi)


#Adding animal ID labels to relative abundance plot

#Removing biosamples from sample name for better plotting 
gen_plot_fungi2 <- separate_wider_delim(gen_plot_fungi, cols = "Sample", delim = "_", names = c("Sample", "biosample"))

#Adding x-axis labels for samples
gen_rel_plot_fungi2 <- ggplot(gen_plot_fungi2) +
  geom_col(mapping = aes(x = Sample, y = RelAb, fill = fct_relevel(Genus, c("Other (< 3%)"), after = 29)), 
           color = "black", position = "stack", show.legend = TRUE, alpha = 0.6) +
  ylab("Relative Abundance") +
  xlab(NULL) +
  scale_fill_manual(values = extend_pal(palette_length)) + 
  theme_bw() +
  theme(legend.text = element_text(size = 9.25),
        legend.position = "right",
        legend.title = element_text(face="bold", size = 14),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 15, vjust = 0.6, size = 7)) +
  labs(fill = "Genus") +
  facet_grid(cols = vars(fct_relevel(group_broad, c("Nashville Zoo", "Middle TN Recapture", "Middle TN Wild"))), 
             scales = "free_x", space = "free_x") +
  #  force_panelsizes(cols = c(1.32, 0.27, 0.2, 0.33, 0.25)) +
  guides(fill = guide_legend(ncol = 2))
print(gen_rel_plot_fungi2)

#Note: to get the relative abundance percentages of genera for each group, take out the "Sample" variable
#in the script above (for the all_data and ord_plot objects) so that it calculates the rel. abund. of 
#bacterial genera based on the group rather than on a sample by sample basis



##### Determining how many OTUs were shared between the different groups #####

#Removing "unassigned" OTUs from the phyloseq object
physeq_fungi_sub <- subset_taxa(physeq_fungi, Kingdom != "Unassigned")

#Subsetting the original phyloseq according to group/env. setting
ps_zoo <- prune_samples(physeq_fungi_sub@sam_data$group_broad == "Nashville Zoo", physeq_fungi_sub)
ps_recap <- prune_samples(physeq_fungi_sub@sam_data$group_broad == "Middle TN Recapture", physeq_fungi_sub)
ps_mtn <- prune_samples(physeq_fungi_sub@sam_data$group_broad == "Middle TN Wild", physeq_fungi_sub)

#Filtering for OTUs that are present in the samples for each group (relative abundance greater than 0)
ps_zoo <- filter_taxa(ps_zoo, function(x) sum(x) > 0, prune = TRUE)
ps_recap <- filter_taxa(ps_recap, function(x) sum(x) > 0, prune = TRUE)
ps_mtn <- filter_taxa(ps_mtn, function(x) sum(x) > 0, prune = TRUE)

#Creating the vectors of the OTUs in each group
zoo_otus <- row.names(otu_table(ps_zoo))
recap_otus <- row.names(otu_table(ps_recap))
mtn_otus <- row.names(otu_table(ps_mtn))

#Intersecting the OTU vectors to determine the shared OTUs
zoo_recap <- intersect(zoo_otus, recap_otus)
recap_mtn <- intersect(recap_otus, mtn_otus)
zoo_mtn <- intersect(zoo_otus, mtn_otus)

#Counting the number of shared OTUs between groups
length(zoo_recap)
length(recap_mtn)
length(zoo_mtn)

#Creating tables to see which particular taxa are shared
tax_table <- tax_mat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "otu")

tax_zoo_recap <- tax_table %>% filter(otu %in% zoo_recap)
tax_recap_mtn <- tax_table %>% filter(otu %in% recap_mtn)
tax_zoo_mtn <- tax_table %>% filter(otu %in% zoo_mtn)

#Saving the shared taxa files
#write.csv(tax_zoo_recap, "250203_shared_fungi_zoo_recap.csv")
#write.csv(tax_recap_mtn, "250203_shared_fungi_recap_midTNWild.csv")
#write.csv(tax_zoo_mtn, "250203_shared_fungi_zoo_midTNWild.csv")



##### Shared OTUs venn diagrams #####


## Without Unassigned OTUs ##

#Attempting to create venn diagram for pre-release, recaps, and wild HBS
fungi_venn_sub <- ggvenn_pq(physeq = physeq_fungi_sub, fact = "group_broad")

#Editing the venn diagram to look a little prettier
fungi_venn_sub_plot <- fungi_venn_sub +
  #ggplot2::scale_color_manual(c("pink3", "lightsteelblue1", "rosybrown", "blue")) +
  scale_fill_distiller(palette = "RdPu", direction = 1) +
  labs(fill = "OTU count")
  #guides(fill = "none")
print(fungi_venn_sub_plot)


## With Unassigned OTUs ##

#Naming new phyloseq object
#physeq_fungi_all <- physeq_fungi
#physeq_fungi_all

#Attempting to create venn diagram for pre-release, recaps, and wild HBS
#fungi_venn_all <- ggvenn_pq(physeq = physeq_fungi_all, fact = "group_broad")

#Editing the venn diagram to look a little prettier
#fungi_venn_all_plot <- fungi_venn_all +
#  #ggplot2::scale_color_manual(c("pink3", "lightsteelblue1", "rosybrown", "blue")) +
#  scale_fill_distiller(palette = "RdPu", direction = 1) +
#  labs(fill = "OTU count")
#guides(fill = "none")
#print(fungi_venn_all_plot)



##### Calculating Hill numbers for alpha diversity and adding to the metadata file #####

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
otu_table_2 <- otu_table_2 %>% select(all_of(intersected))



meta_df$h0 <- hill_taxa(otu_table_2, q = 0, MARGIN = 1)
#q = 0 is essentially a richness metric, placing more emphasis on rare taxa

#meta_df$h1 <- hill_taxa(otu_table_2, q = 1, MARGIN = 1)
#q = 1 is Hill-Shannon diversity, where abundant and rare taxa are equally weighted
#and evenness is accounted for

#meta_df$h2 <- hill_taxa(otu_table_2, q = 2, MARGIN = 1)
#q = 2 is inverse Simpson, which places more emphasis on common/abundant species

#meta_df$h3 <- hill_taxa(otu_table_2, q = 3, MARGIN = 1)
#q = 3 is where dominant, abundant species are given the most weight

#Calculating average richness and SD for the groups and by captive vs. wild
meta_df %>% 
  group_by(group_broad) %>% 
  summarize(mean = mean(h0, na.rm = T), 
            sd = sd(h0, na.rm = T))

meta_df %>% 
  group_by(env_broad_scale) %>% 
  summarize(mean = mean(h0, na.rm = T), 
            sd = sd(h0, na.rm = T))



##### Calculating the number of phyla, class, order, etc. in all HB samples #####

tax_table(physeq_fungi_sub) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  gather("Rank", "Name", rank_names(physeq_fungi_sub)) %>%
  na.omit() %>% # remove rows with NA value
  group_by(Rank) %>%
  summarize(ntaxa = length(unique(Name))) %>% # compute number of unique taxa
  mutate(Rank = factor(Rank, rank_names(physeq_fungi_sub))) %>%
  arrange(Rank)
