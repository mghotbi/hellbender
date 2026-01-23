# Prepare coding environment & load data #####

#load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(data.table)
  library(decontam)
  library(vegan)
  library(ggtext)
})

#prioritize dplyr functions
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

#get date for file saving 
date <- format(Sys.Date(), "%y%m%d")


#set working directory and random seed
set.seed(11301996)

setwd("H:/Shared drives/Microbiome Ecology Lab Shared Drive/Jason Dallas PC backup/mentoring/chloe_cummins/240819_16S_mothur_output")

tax.file = fread("240819_hellbender_16S_taxonomy.tsv")
count.file = fread('240819_hellbender_16S_count_table.csv')
metadata.file = fread('240819_hellbender_16S_metadata.csv', stringsAsFactors=T)

###  find the OTUs that are Tetragenoccocus

tetra.tax <- tax.file %>% 
  filter(genus == 'Tetragenococcus') #196 OTUs assigned to Tetragenococcus, only 1 (Otu00003) is actually abundant represents 99.7% of all Tetra reads

tetra_intersect <- intersect(tetra.tax$OTU, colnames(count.file)) #create a list of the Tetra OTUs

###  format count table

otu.mat <- count.file %>% 
  select(-c(label, numOtus)) %>% #remove bioinformatics metadata 
  select(-all_of(tetra_intersect)) %>% 
  rename('index' = 'Group') %>%
  mutate(index = case_when(row_number() == 1 ~'NTC.swab_35815',
                           row_number() == 2 ~'NTC.swab_36113',
                           row_number() == 3 ~'Pool.neg.ctrl_34523',
                           row_number() == 4 ~'Pos.control_34522',
                           row_number() == 214 ~'poolNTC.swab_34624',
                           row_number() == 215 ~'spiked.swab_34623',
                           row_number() == 216 ~'spiked1ul.swab_35814',
                           TRUE ~ index)) %>% 
  column_to_rownames('index') #ensure dataframe contains only numeric columns 


#format taxonomy file 
tax <- tax.file %>% 
  filter(OTU %in% names(otu.mat)) %>%   #retain only otus with sequencing data (in taxonomy file) 
  filter(!genus == 'Tetragenococcus') %>%  #remove the spike in genus
  rename('otu' = 'OTU')

###   update and format metadata file

metadata <- metadata.file %>% 
  mutate(sample_name.x = case_when(row_number() == 233 ~'NTC.swab',
                           row_number() == 238 ~'NTC.swab',
                           row_number() == 114 ~'Pool.neg.ctrl',
                           row_number() == 214 ~'poolNTC.swab',
                           row_number() == 232 ~'spiked1ul.swab',
                           TRUE ~ sample_name.x)) %>% 
  unite(col = 'index', c('sample_name.x', 'sequence_biosample'), sep = '_') %>% 
  filter(index %in% rownames(otu.mat)) %>% #retain samples that have undergone sequencing
  mutate(across(all_of(c('diet_crayfish', 'env_broad_scale', 'site')), as.factor))

###  filter the otu matrix to only include samples in metadata

otu.mat <- otu.mat %>% 
  rownames_to_column('index') %>% 
  filter(index %in% metadata$index) %>% 
  column_to_rownames('index')

#clean up working environment
keep.objs <- c("otu.mat",
               "tax",
               "metadata",
               "select",
               "filter",
               "date")
remove.objs <- subset(ls(), !(ls() %in% keep.objs)) %>% 
  subset(. != "keep.objs")
rm(list = remove.objs) 


# Remove rare taxa #####

#filter rare taxa using abundance cutoff 
otu.mat %>% 
  select(where( ~ is.numeric(.x) && sum(.x) > 10)) -> otu.mat.abund

#clean up working environment
keep.objs <- subset(keep.objs, !keep.objs %in% c("otu.mat")) %>% 
  c(., "otu.mat.abund")
remove.objs <- subset(ls(), !(ls() %in% keep.objs)) %>% 
  subset(. != "keep.objs")
rm(list = remove.objs) 

#     Perform decontamination ####

#Variable that determines if computationally expensive calculation should be repeated
calculate <- "Y"
#file naming parameters
environ.desig <- "decontam.out" #file name in R environment
file.suffix <- "decontam.out.rds" #file suffix in windows environment 


#either load previous iteration or perform new computation based on user input 
if(calculate == "N"){
  #first check if they file already exists in R environment
  if(environ.desig %in% ls()){
    print("file already exists in Environment")
    #next check if there is a file created today in the working directory
  } else if (paste(date, file.suffix, sep = "_") %in% list.files()){
    assign(environ.desig, readRDS(paste(date, file.suffix, sep = "_")))
    print("loaded file (generated today) from working directory")
    #next load the last file that was generated
  } else if (any(grepl(file.suffix, list.files()))) {
    #get the date the last file was generated
    subset(list.files(), grepl(file.suffix, list.files())) %>%
      str_replace(., file.suffix, "") %>%
      str_replace(., "_", "") %>%
      subset(., . != "") %>%
      as.Date(., format = "%y%m%d") %>%
      sort() %>%
      tail(n = 1) %>%
      format("%y%m%d") -> lastfile.date
    #use it to get the file name
    subset(list.files(), grepl(lastfile.date, list.files())) %>%
      subset(., grepl(file.suffix, .)) -> lastfile.name
    #load file from working directory
    assign(environ.desig, readRDS(lastfile.name))
    print(paste("loaded file (generated ",as.Date(lastfile.date, format = "%y%m%d"),") from working directory", sep = ""))
    #finally indicate there is no file to be loaded if that is the case 
  } else {
    print("saved objection representation does not exist in working directory")    
  }
} else if(calculate == "Y"){
  
  ####conduct decontamination via decontam
  print("Performing decontamination")
  
  #Reorder sequence metadata to match OTU matrix 
  otu.mat.abund %>%
    rownames_to_column("index") %>%
    select(c(1)) %>% 
    left_join(.,
              metadata) -> metadata 
  
  #join OTU feature table with sequencing metadata to perform decontamination
  otu.mat.abund %>% 
    rownames_to_column("index") %>% 
    right_join(metadata, ., by = "index") %>% 
    #remove any samples without a known DNA concentration
    filter(!is.na(ampliconlibrary_quantification_ng.ul)) %>%
    #annotate samples are NTC or actual samples with boolean
    mutate(sample_type = (grepl("pool", index, ignore.case = T)|
                            grepl("swab", index, ignore.case = T)),
           .after = index) -> otus.seqdat
  
  #remove rownames - artefact from OTU matrix
  row.names(otus.seqdat) <- NULL
  
  #specify sequence of threshold probability values to iterate through 
  threshold <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)

  
  #generate empty lists for loop output items 
  sapply(c("decontam.results.list",
           'decontam.summary',
           "plots.sampletype",
           "plots.persample"),
         function(k) assign(k, list(), envir = .GlobalEnv))
  
  #format #sequence data for decontam analysis
  otus.seqdat %>% 
    column_to_rownames("index") %>% 
    select(starts_with("Otu")) %>% 
    as.matrix() -> decontam.seqs 
  
  # conduct decontamination across range of threshold values 
  for(j in 1:length(threshold)){
    
    #report threshold value to track loop 
    print(paste("Threshold:", threshold[j], sep = " "))
    
    #create var for sys.time() and track time
    t.start <- Sys.time()
    
    #identify contaminating taxa via decontam 
    decontam.results <- isContaminant(
      seqtab = decontam.seqs,
      neg = otus.seqdat$sample_type,
      conc = otus.seqdat$ampliconlibrary_quantification_ng.ul,
      method = "combined",
      threshold = threshold[j],
      normalize = TRUE,
      detailed = T) %>%
      #retain only contaminating taxa 
      filter(contaminant == T) %>% 
      #add taxonomy data to results 
      rownames_to_column("otu") %>%
      left_join(., tax, by = "otu")
    
    #save results into list
    decontam.results.list[[j]] <- decontam.results
    
    #name each dataframe according to the data set and threshold value 
    names(decontam.results.list)[j] <- paste("Threshold", threshold[j], sep=": ")
    
    #extract list of contaminant OTUs 
    contaminants <- decontam.results.list[[j]]$otu
    
    #generate long dataframe and classify taxa according to decontam results 
    df.otu.long <- otus.seqdat %>%
      pivot_longer(starts_with("Otu"), names_to = "otu", values_to = "reads") %>% 
      mutate(otu.type = ifelse(otu %in% contaminants, "contamination", "sample taxa")) %>% 
      group_by(sample_type) %>% 
      mutate(reads = as.numeric(reads),
             total.reads = sum(reads))  %>% 
      ungroup() %>% 
      group_by(index) %>% 
      mutate(sample.reads = sum(reads)) %>%
      ungroup() %>% 
      mutate(otu.type = as.factor(str_to_title(otu.type)))
    
    #Ensure plots always have a full legend 
    if(!("contamination" %in% levels(df.otu.long$otu.type))){
      levels(df.otu.long$otu.type)[length(levels(df.otu.long$otu.type))+1] <- "Contamination"
    }
    
    #plot decontam results according to sample type 
    #generate dataframe for plotting
    df.plot.sampletype <- df.otu.long %>%  
      group_by(index, otu.type) %>%
      mutate(otu.id.reads = sum(reads),
             otu.id.proportion = otu.id.reads/sample.reads) %>%
      sample_n(1) %>%   
      group_by(sample_type, otu.type) %>% 
      mutate(mean.otu.id.pro = mean(otu.id.proportion),
             se.otu.id.pro = sd(otu.id.proportion)/sqrt(length(otu.id.proportion))) %>%
      mutate(sample_type = ifelse(sample_type, "No Template Control", "Sample")) %>% 
      ungroup()
    #render plot
    plot.sampletype <- df.plot.sampletype %>% 
      ggplot(aes(x = sample_type, y = -mean.otu.id.pro, fill = otu.type)) +
      geom_errorbar(aes(ymax = -(mean.otu.id.pro + se.otu.id.pro),
                        ymin = -(mean.otu.id.pro - se.otu.id.pro)),
                    position = position_dodge(width = 0.9), width = 0.25) +
      geom_col(color = "black", position = position_dodge(width = 0.9)) + 
      scale_y_reverse(limits = c(0,-1), labels = scales::percent) +
      scale_fill_manual(drop = FALSE,
                        values = c("Contamination" = "#F8766D", "Sample Taxa" = "#00BFC4")) +
      labs(title = paste("Threshold =", threshold[j], sep = " "),
           x = "Sample Type",
           y = "Sample Composition") + 
      theme_classic() +
      theme(plot.title = element_text(size = 12),
            legend.title = element_blank())
    
    #extract summary values for decontamination
    decontam.summary[[j]] <- df.plot.sampletype %>%
      mutate(threshold = threshold[j]) %>%
      as.data.frame()
    
    #save plots in a named list
    plots.sampletype[[j]] <- plot.sampletype
    names(plots.sampletype)[j] <- paste("Threshold", threshold[j], sep=": ")
    
    #plot decontam results on a per-sample basis 
    df.persample <- df.otu.long %>%
      group_by(index) %>% 
      mutate(sample.reads = sum(reads)) %>% 
      ungroup() %>%
      mutate(index = as.factor(index),
             sample.reads = as.numeric(sample.reads),
             index = fct_reorder(index, sample.reads, .desc = T)) %>%
      group_by(index, otu.type) %>%
      mutate(read.counts = sum(reads)) %>% 
      sample_n(1) %>% 
      select(-otu) 
    
    plot.persample <- df.persample %>%
      filter(read.counts > 0) %>% 
      ggplot(aes(x = index, y = read.counts, fill = otu.type)) + 
      geom_col(aes(y = read.counts-100000, x = index, fill = otu.type), color = "black", alpha = 1, linewidth = 0.5) + 
      geom_col(width = 1, show.legend = F) +
      scale_y_continuous(expand = c(0,0)) + 
      coord_cartesian(ylim = c(0,as.numeric(quantile(df.persample$read.counts, 0.999)))) + 
      labs(title = paste("Threshold =", threshold[j], sep = " ")) + 
      scale_fill_manual(values = c("Contamination" = "#F8766D", "Sample Taxa" = "#00BFC4")) +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 12))
    
    #save plots in a named list
    plots.persample[[j]] <- plot.persample
    names(plots.persample)[j] <- paste("Threshold", threshold[j], sep=": ")
    
    #report time to track loop 
    t.deco <- abs(t.start-Sys.time())
    print(paste("Loop", " ", j, "/", length(threshold), "    ", "Elapsed time: ", round(t.deco, 4), sep = ""))
    
  }
  
  # generate plot to select optimal threshold value by examining
  # how many taxa are removed at each step change in threshold value
  summary.df <- decontam.summary %>% 
    bind_rows() %>%
    select(-c(otu, reads)) %>% 
    filter(otu.type == "Contamination") %>% 
    group_by(threshold, sample_type) %>% 
    mutate(sample.type.reads = total.reads, 
           thresh.contam.reads = sum(otu.id.reads)) %>% 
    sample_n(1) %>%
    ungroup() %>% 
    mutate(thresh.contam.pro = thresh.contam.reads/sample.type.reads,
           delta.thresh.contam = -(thresh.contam.pro - shift(thresh.contam.pro, n = 2))) %>% 
    mutate(delta.thresh.contam = ifelse(is.na(delta.thresh.contam), thresh.contam.pro, delta.thresh.contam))
  summary.df %>% 
    ggplot(aes(x = as.factor(threshold), y = delta.thresh.contam, fill = sample_type)) +
    geom_col(position = position_dodge(), color = "black") +
    scale_y_reverse(labels = scales::percent) +
    labs(y = "&Delta; Total Reads",
         x = "Threshold Value") +
    scale_fill_brewer(palette = "Paired") +
    theme_classic() +
    theme(legend.title = element_blank(), 
          legend.position = "bottom", 
          text = element_text(size = rel(4.5)),
          legend.text = element_text(size = 14),
          legend.spacing.x = unit(0.5, 'cm'),
          axis.title.x = element_text(margin = margin(10,0,0,0)),
          axis.title.y = ggtext::element_markdown(margin = margin(0,15,0,10))
    ) -> decontam.summary.plot
  
  #store decontam output in tidy list
  decontam.out <- list(
    'decontam.results.list' = decontam.results.list,
    'plots.persample' = plots.persample, 
    'plots.sampletype' = plots.sampletype,
    'decontam.summary' = decontam.summary,
    'summary.df' = summary.df,
    'decontam.summary.plot' = decontam.summary.plot
  )
  
  #save decontam output as an .rds object to prevent recursive computational costs
  saveRDS(decontam.out, file = paste(date, "decontam.out.rds", sep = "_"))
  
}

#examine results of decontamination
(decontam.summary.plot<-decontam.out$decontam.summary.plot)

plots.sampletype
plots.persample

#analytically select threshold for final decontamination 
#
# i.e., first threshold where proportionally more seqs are removed from NTC than sample libraries 
# after removing 10% of ntc library seqs (background noise) 
decontam.out$summary.df %>% 
  group_by(sample_type) %>% 
  mutate(reads.remove.cumulative = cumsum(thresh.contam.pro)) %>% 
  select(threshold, sample_type, thresh.contam.pro, reads.remove.cumulative) %>% 
  mutate(ntc.reads.remove.cumulative = ifelse(sample_type == "No Template Control", reads.remove.cumulative, NA)) %>% 
  group_by(threshold) %>% 
  fill(ntc.reads.remove.cumulative) %>% 
  ungroup() %>% 
  filter(ntc.reads.remove.cumulative > .1) %>% 
  select(threshold, sample_type, thresh.contam.pro) %>% 
  mutate(sample_type = as.factor(sample_type),
         sample_type = fct_recode(sample_type, `ntc` = "No Template Control", `sample` = "Sample")) %>% 
  pivot_wider(names_from = "sample_type", values_from = 3) %>% 
  filter(ntc > sample) %>% 
  arrange(threshold) %>% 
  slice_head(n = 1) %>% 
  pull(threshold) -> threshold.decontam
threshold.decontam.index <- paste("Threshold:", threshold.decontam); threshold.decontam.index #threshold set to 0.45
#As a note, the code above said that analytically the best threshold cutoff was 0.45, which is what was used
#However, this still contained the E. coli genus which produced a very diff. plot for the 0.55 and 0.65 cutoffs
#Look back at Jason's message in Slack on August 20th, 2024 for more info

#extract decontam results for selected threshold value
decontam.final <- decontam.out[["decontam.results.list"]][[threshold.decontam.index]]

#remove contaminant taxa from OTU abundance matrix
contaminants <- decontam.final$otu
otu.mat.abund.decontam <- select(otu.mat.abund, -any_of(contaminants))


#generate figures for supplementary materials 
decontam.out$plots.persample[[threshold.decontam.index]] +
  labs(x = "Samples<br><span style = 'font-size:11pt'>(arranged by sequencing depth)</span>",
       y = "16s Reads") +
  theme(legend.title = element_blank(),
        plot.title = element_blank(), 
        legend.position = c(0.8,0.9), 
        text = element_text(size = rel(4.5)),
        legend.text = element_text(size = 14),
        legend.spacing.x = unit(0.5, 'cm'),
        axis.title.x = element_markdown(margin = margin(10,0,0,0)),
        axis.title.y = element_text(margin = margin(0,15,0,10))
  ) -> persample.plot 

decontam.out$plots.sampletype[[threshold.decontam.index]] +
  labs(y = "% Reads<br><span style = 'font-size:12pt'>(Mean &plusmn; SE)</span>") +
  theme(legend.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "none",
        text = element_text(size = rel(4.5)),
        legend.text = element_text(size = 14),
        legend.spacing.x = unit(0.5, 'cm'),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(margin = margin(0,15,0,10))
  ) -> sampletype.plot

decontam.summary.plot <- decontam.summary.plot +
  theme(plot.margin = margin(0,20,20,35),
        legend.position = c(0.1, 0.9))


#multipanel figure
ggpubr::ggarrange(
  decontam.summary.plot,
  ggpubr::ggarrange(sampletype.plot,
                    persample.plot,
                    ncol = 2,
                    align = "hv",
                    labels = c("b", 'c'),
                    font.label = list(size = 24, color = "black", face = "bold", family = "sans")
  ),
  labels = c("a", ''),
  font.label = list(size = 24, color = "black", face = "bold", family = "sans"),
  nrow = 2
) 


#clean up working environment
keep.objs <- subset(keep.objs, !keep.objs %in% c("otu.mat.abund", "decontam.out")) %>%
  c(., "otu.mat.abund.decontam") %>%
  unique()
remove.objs <- subset(ls(), !(ls() %in% keep.objs)) %>%
  subset(. != "keep.objs")
rm(list = remove.objs)



######     Perform rarefaction & remove no template controls ####

#check the number of samples retained across a range of sequencing depths 
range.rarefy <- seq(0, 50)*1000   #sequencing depth range
rarefy.summary <- list()
for(i in 1:length(range.rarefy)){
  rarefy.depth <- range.rarefy[i]
  otu.mat.abund.decontam %>% 
    as.data.frame() %>% 
    filter(rowSums(.) >= rarefy.depth) %>% 
    nrow() %>%
    as.numeric() -> samplesretained
  data.frame(
    samples.retained = samplesretained,
    depth = rarefy.depth,
    per.samples.retained = samplesretained/nrow(otu.mat.abund.decontam)
  ) -> rarefy.summary[[i]]
}
do.call(rbind, rarefy.summary) -> rarefy.summary

view(rarefy.summary)

#rarefy to 6k results in 88% of samples retained

#generate plot to select the most appropriate rarefaction depth 
rarefy.summary %>% 
  ggplot(aes(x = depth, y = per.samples.retained)) + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) + 
  #drawing horizontal line at 0.87850467% of samples retained
  geom_segment(aes(x=min(depth),xend=max(depth),y=0.87850467,yend=0.87850467), lty = "dashed", color = "red", linewidth = 1) +
  #drawing vertical line at 6k rarefacation depth
  geom_segment(aes(y=0,yend=1,x=6000,xend=6000), lty = "dashed" , color = "blue", linewidth = 1) + 
  geom_line(linewidth = 1.5) +
  labs(y = "Samples Retained",
       x = "Rarefaction Depth") + 
  theme_classic() +
  theme(text = element_text(size = rel(4.5)),
        axis.title.y = element_text(margin = margin(0,15,0,5)),
        axis.title.x = element_text(margin = margin(10,0,10,0)))


#proceed with rarefaction
rare.depth <- 6000 #explicitly state
suppressWarnings({otu.mat.abund.decontam %>% 
    rrarefy(rare.depth) %>%  #rarefy to 6k
    as.data.frame() %>% 
    filter(rare.depth == rowSums(.)) %>% 
    #remove any no template controls from the OTU table
    rownames_to_column("swab_label") %>% 
    mutate(sample_type = !(grepl("neg", swab_label, ignore.case = T)|
                             grepl("nc", swab_label, ignore.case = T))) %>% 
    filter(sample_type) %>% 
    select(-sample_type) %>% 
    mutate(swab_label = fct_recode(swab_label, `CMFP39.1` = "CMFP39")) %>% 
    column_to_rownames("swab_label") %>% 
    #store as matrix 
    as.matrix()}) -> otu.mat.rare


#collect fecal sample names for those samples which passed bioinformatics pipeline
final.sample.names <- rownames(otu.mat.rare)

#collect otu labels for those otus which passed bioinformatics pipeline
final.otus <- colnames(otu.mat.rare)

#subset taxonomy file a final time
tax.final <- filter(tax, otu %in% final.otus)

#clean up working environment
keep.objs <- keep.objs %>% 
  subset(!. %in% c('otu.mat.abund.decontam')) %>% 
  c("otu.mat.rare",
    "final.sample.names",
    "final.otus",
    "tax.final")
remove.objs <- subset(ls(), !(ls() %in% keep.objs)) %>%
  subset(. != "keep.objs")
rm(list = remove.objs)

#####   load and format metadata ####

#all samples that passed bioinformatics are in master metadata files 
length(setdiff(rownames(otu.mat.rare), metadata$index)) #all samples are present in both datasets


#filter taxonomy data final time 
tax.final <- tax.final %>% 
  filter(otu %in% colnames(otu.mat.rare))

#join metadata with otu data to produce a dataframe with metadata and count data from rarefied dataset

otu.mat.rare %>% 
  as.data.frame() %>% 
  rownames_to_column("index") %>%
  right_join(metadata, ., by = "index") %>% 
  filter(sample_or_blank == 'sample') -> final.df



#export data as .csv objects  
write.csv(tax.final, file = paste(date, "hellbender_taxonomy.csv", sep = "_"), row.names = F) #taxonomy data 
write.csv(final.df, file = paste(date, "hellbender_decontam_rarefied_dataset.csv", sep = "_"), row.names = F) #metadata and OTU data 

#   attempt other methods of decontamination ####

#create phyloseq

library(phyloseq)

otu_mat <- as.matrix(otu.mat.abund)

taxonomy.mat = tax %>% 
  column_to_rownames('otu') 

tax_mat <- as.matrix(taxonomy.mat)

metadata.2  <- metadata %>% 
  column_to_rownames('index')

OTU = otu_table(otu_mat, taxa_are_rows = F)
TAX = tax_table(tax_mat)
metadata.phyloseq = sample_data(metadata.2)

hellbender.ps = phyloseq(OTU, TAX, metadata.phyloseq)

df <- as.data.frame(sample_data(hellbender.ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(hellbender.ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_or_blank)) + geom_point() + scale_y_continuous(trans='log10')


