

##### SCRIPT TO GENERATE CYTOSCAPE'S INPUT FILES #####

### Load packages ###

library(dplyr)
library(tidyr)
library(readr)

### Set a working directory
setwd("~/Documents/MTSU/Hellbender_project_ChloeCummins")

### 1.  Load OTU table
otu_table <- read.csv("250508_16S_OTU_count_table_decontam.csv", check.names = FALSE)
otu_table <- otu_table %>% rename(sample_id = index)

### 2.  Load metadata table
metadata <- read.csv("250508_16S_metadata_decontam_Lluvia_network_MOD.csv", stringsAsFactors = FALSE)
colnames(metadata)[colnames(metadata) == "index"] <- "sample_id"

### 2.  Load taxonomy table
taxonomy <- read.csv("250508_16S_taxonomy_table_decontam_Lluvia_network.csv", stringsAsFactors = FALSE)

### 3.  Create EDGES input file
edges <- otu_table %>%
  pivot_longer(-sample_id, names_to = "otu_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  left_join(taxonomy %>% rename(otu_id = OTU), by = "otu_id") %>%  
  group_by(sample_id, family) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  rename(from = sample_id, to = family)

write.csv(edges, "cytoscape_edges_family.csv", row.names = FALSE)

### 4.  Create NODES input file
otu_long <- otu_table %>%
  pivot_longer(-sample_id, names_to = "otu_id", values_to = "count") %>%
  left_join(taxonomy %>% rename(otu_id = OTU), by = "otu_id")

### 4.1.  Aggregate counts by family and sample
family_sample_counts <- otu_long %>%
  group_by(family, sample_id) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = sample_id, values_from = count, values_fill = 0)

### 4.2. Compute ZOO/WILD totals per genus.
### I'm inlcuding log10 calculations to manage the size of the nodes.
env_totals <- otu_long %>%
  group_by(family, sample_id) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  left_join(metadata %>% select(sample_id, env_broad_scale), by = "sample_id") %>%
  group_by(family, env_broad_scale) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = env_broad_scale, values_from = total, values_fill = 0) %>%
  mutate(
    absolute_abundance = zoo + wild,
    log10_abundance = log10(absolute_abundance +0.1),
    log10_wild = log10(wild + 0.1),  # +0.1 avoids log(0)
    log10_zoo  = log10(zoo + 0.1)
  )

### 4.3. Combine counts plus env totals
family_nodes <- family_sample_counts %>%
  left_join(env_totals, by = "family") %>%
  mutate(Type = "Family") %>%
  rename(name = family)

### 4.4. Sample nodes
sample_nodes <- metadata %>%
  rename(name = sample_id) %>%
  mutate(Type = "Sample") %>%
  select(name, Type, env_broad_scale, group_broad)

### 4.5. Combine family plus nodes.
### NOTE: since we want to visualize in the network family nodes as pies (proportion of wild counts and proportion of zoo counts),
### we will keep those columns that we create in the above step 4.2.
nodes <- bind_rows(
  family_nodes %>% select(name, Type, zoo, wild, everything()),  # keep sample columns for pies
  sample_nodes %>% mutate(zoo = NA, wild = NA)
)

write.csv(nodes, "cytoscape_nodes_pie_family.csv", row.names = FALSE)
