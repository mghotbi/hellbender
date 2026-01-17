# hellbender
***Effects of environmental setting and diet on the gut microbial ecology of eastern hellbenders (Cryptobranchus alleganiensis alleganiensis)***

Why this script exists

Many microbiome studies summarize results at higher taxonomic ranks (e.g. family or order). However, ecological patterns observed at these levels may be driven by different underlying OTUs, which can mask biologically meaningful differences.

This script was developed to address the following questions:

Are similar family- or order-level abundances across groups driven by the same OTUs, or by different OTUs within the same family?
Which OTUs form the core microbiome within and across host or environmental groups?
Do taxa traditionally associated with a given function (e.g. Chitinophagaceae) share OTUs across groups, or only share taxonomic labels?


---
What the script does

This script performs OTU-level core microbiome analysis and family-specific OTU visualization using a phyloseq object as input.

Specifically, it:
Builds a phyloseq object from:

OTU count table
Taxonomy table
Sample metadata

Defines core OTUs per group based on prevalence (DspikeIn::plot_core_microbiome_custom(); https://github.com/mghotbi/DspikeIn)

Core = OTU present in ≥ X proportion of samples within a group
Threshold is user-controlled (prev_core_threshold)
Identifies OTUs belonging to a target family/order
Example: Chitinophagales / Chitinophagaceae

Can be extended to all core families
Calculates OTU prevalence by group
Outputs a table suitable for Supplementary Materials
Generates OTU-level heatmaps
Values = log10(relative abundance + pseudocount)
Rows = OTUs
Columns = samples
Annotated by group
Hierarchical clustering with safety checks (no errors if only one OTU)
Saves outputs as PDF + PNG heatmaps

Prevalence tables

Short explanation text for methods

